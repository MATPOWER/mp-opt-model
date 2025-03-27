function [x, f, eflag, output, lambda] = qcqps_gurobi(H, b, Q, C, k, l1, u1, A, l2, u2, xmin, xmax, x0, opt)
% qcqps_gurobi - Quadratically Constrained Quadratic Program Solver based on GUROBI.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QCQPS_GUROBI(H, B, Q, C, K, L1, U1, A, L2, U2, XMIN, XMAX, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = QCQPS_GUROBI(PROBLEM)
%   A wrapper function providing a standardized interface for using
%   GUROBI to solve the following (possibly non-convex) QCQP (quadratically 
%   constrained quadratic programming) problem:
%
%       min X'*H*X + B'*X
%        X
%
%   subject to
%
%    L1(i) <= X'*Q{i}*X + C(i,:)*X + K(i) <= U1(i),  i = 1,2,...,NQ   (quadratic constraints)
%                 L2 <= A*X <= U2                                     (linear constraints)
%                XMIN <= X <= XMAX                                    (variable bounds)
%
%   Inputs (all optional except H, B, Q, C, K, L1 and U1):
%       H : matrix (possibly sparse) of quadratic cost coefficients
%       B : vector of linear cost coefficients
%       Q : NQ x 1 cell array of sparse quadratic matrices for quadratic constraints
%       C : matrix (posibly sparse) of linear parameters of quadratic constraints
%       K : vector of constant parameters of quadratic constraints
%       L1, U1: define the lower an upper bounds on the quadratic constraints
%       A, L2, U2 : define the optional linear constraints. Default
%           values for the elements of L and U are -Inf and Inf,
%           respectively.
%       XMIN, XMAX : optional lower and upper bounds on the
%           X variables, defaults are -Inf and Inf, respectively.
%       X0 : optional starting value of optimization vector X
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%               3 = even more verbose progress output
%           grb_opt - options struct for GUROBI, value in verbose
%                   overrides these options
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, B, Q, C, k, l1, u1, A, l2, u2, xmin, xmax, x0, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : exit flag
%           1 = converged
%           0 or negative values = solver specific failure codes
%       OUTPUT : output struct with the following fields:
%           alg - algorithm code of solver used
%           (others) - algorithm specific fields
%       LAMBDA : struct containing the Langrange and Kuhn-Tucker
%           multipliers on the constraints, with fields:
%           mu_l - lower (left-hand) limit on linear constraints
%           mu_u - upper (right-hand) limit on linear constraints
%           mu_l_quad - (QCQP only) lower (left-hand) limit on quadratic 
%                       constraints
%           mu_u_quad - (QCQP only) upper (right-hand) limit on quadratic 
%                       constraints
%           lower - lower bound on optimization variables
%           upper - upper bound on optimization variables
%
% See also qcqps_master, gurobi_options, gurobi.

%   MP-Opt-Model
%   Copyright (c) 2019-2023, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia
%
%   This file is part of MP-Opt-Model..
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%%----- input argument handling  -----
%% gather inputs
if nargin == 1 && isstruct(H)       %% problem struct
    p = H;
    if isfield(p, 'opt'),   opt = p.opt;    else,   opt = [];   end
    if isfield(p, 'x0'),    x0 = p.x0;      else,   x0 = [];    end
    if isfield(p, 'xmax'),  xmax = p.xmax;  else,   xmax = [];  end
    if isfield(p, 'xmin'),  xmin = p.xmin;  else,   xmin = [];  end
    if isfield(p, 'u2'),    u2 = p.u2;      else,   u2 = [];    end
    if isfield(p, 'l2'),    l2 = p.l2;      else,   l2 = [];    end
    if isfield(p, 'A'),     A = p.A;        else,   A = [];     end
    if isfield(p, 'u1'),    u1 = p.u1;      else,   u1 = [];    end
    if isfield(p, 'l1'),    l1 = p.l1;      else,   l1 = [];    end
    if isfield(p, 'k'),     k = p.k;        else,   k = [];     end
    if isfield(p, 'C'),     C = p.C;        else,   C = [];     end
    if isfield(p, 'Q'),     Q = p.Q;        else,   Q = {};     end
    if isfield(p, 'b'),     b = p.b;        else,   b = [];     end
    if isfield(p, 'H'),     H = p.H;        else,   H = [];     end
else                                %% individual args
    if nargin < 14
        opt = [];
        if nargin < 13
            x0 = [];
            if nargin < 12
                xmax = [];
                if nargin < 11
                    xmin = [];
                    if nargin < 8
                        A = [];
                        l2 = [];
                        u2 = [];
                    end
                end
            end
        end
    end
end

%% define nx, set default values for missing optional inputs
if isempty(Q)
    if isempty(H) || ~any(any(H))
        if isempty(C) && isempty(k) && isempty(A) && isempty(xmin) && isempty(xmax)
            error('qcqps_gurobi: LP problem must include constraints or variable bounds');
        else
            if ~isempty(A) && ~isempty(C)
                if size(A,2) == size(C,2)
                    nx = size(A, 2);
                else
                    error('qcqp_gubori: number of columns of A and C must agree')
                end
            elseif ~isempty(xmin)
                nx = length(xmin);
            else    % if ~isempty(xmax)
                nx = length(xmax);
            end
        end
        H = sparse(nx,nx);
    else
        nx = size(H, 1);
    end
    if isempty(C) || (~isempty(C) && (isempty(l1) || all(l1 == -Inf)) && ...
                                     (isempty(u1) || all(u1 == Inf)))
        C = sparse(0,nx);           %% no l1 & u1 limits => no quadratic constraints
        k = [];
    end
    nrowC = size(C, 1);             %% number of original quadratic constraints
    if isempty(u1)                  %% By default, quadratic inequalities are ...
        u1 = Inf(nrowC, 1);         %% ... unbounded above and ...
    end
    if isempty(l1)
        l1 = -Inf(nrowC, 1);        %% ... unbounded below.
    end
else
    [nrowQ, ncolQ] = size(Q);
    if isa(Q,'cell') && ncolQ == 1
        size_Q  = cell2mat(cellfun(@(x)(size(x)), Q, 'UniformOutput', false));
        if ~isempty(H)
            if abs(sum((1/nrowQ)*size_Q(:)) - sum(size(H)))  > 1e-10
                error('qcqp_gurobi: Dimensions of matrices H and Q{i}, i=1,2,...,%d must agree.', nrowQ)
            end
            nx = size(H, 1);
        else
            if abs(sum((1/nrowQ)*size_Q(:,2)) - length(b))  > 1e-10
                error('qcqp_gurobi: Dimensions of matrices Q{i}, i=1,2,...,%d and vector b must agree.', nrowQ)
            end
            nx = length(b);
        end
    else
        error('qcqp_gurobi: Input argument Q must be an N x 1 cell array')
    end
end
if isempty(b)
    b = zeros(nx, 1);
end
nrowA = size(A, 1);                %% number of original linear constraints
if isempty(u2)                   %% By default, linear inequalities are ...
    u2 = Inf(nrowA, 1);             %% ... unbounded above and ...
end
if isempty(l2)
    l2 = -Inf(nrowA, 1);            %% ... unbounded below.
end
if isempty(xmin)                %% By default, optimization variables are ...
    xmin = -Inf(nx, 1);             %% ... unbounded below and ...
end
if isempty(xmax)
    xmax = Inf(nx, 1);              %% ... unbounded above.
end
if isempty(x0)
    x0 = zeros(nx, 1);
end

%% default options
if ~isempty(opt) && isfield(opt, 'verbose') && ~isempty(opt.verbose)
    verbose = opt.verbose;
else
    verbose = 0;
end

%% set up options struct for Gurobi
if ~isempty(opt) && isfield(opt, 'grb_opt') && ~isempty(opt.grb_opt)
    g_opt = gurobi_options(opt.grb_opt);
else
    g_opt = gurobi_options;
end
if verbose > 1
    g_opt.LogToConsole = 1;
    g_opt.OutputFlag = 1;
    if verbose > 2
        g_opt.DisplayInterval = 1;
    else
        g_opt.DisplayInterval = 100;
    end
else
    g_opt.LogToConsole = 0;
    g_opt.OutputFlag = 0;
end
g_opt.NonConvex = 2;
g_opt.QCPDual = 1;      %% Turn on the multipliers for quadratic constraints

if ~issparse(A)
    A = sparse(A);
end
if isempty(A)
   A = spalloc(1, nx, 0);
end
if issparse(b)
    b = full(b);
end

%% split up quadratic constraints
[ieq_quad, igt_quad, ilt_quad, Q_quad, C_quad, b_quad] = convert_quad_constraint(Q, C, k, l1, u1);

%% split up linear constraints
[ieq_lin, igt_lin, ilt_lin, A_lin, b_lin] = convert_lin_constraint(A, l2, u2);

%% grab some dimensions
neq_quad = length(ieq_quad);                       %% number of quadratic equalities
niq_quad = length(ilt_quad) + length(igt_quad);    %% number of quadratic inequalities
neq_lin = length(ieq_lin);                         %% number of linear equalities
niq_lin = length(ilt_lin) + length(igt_lin);       %% number of linear inequalities

%% set up model
if ~isempty(Q_quad)
    m.quadcon   = cell2struct([ cellfun(@(x)(0.5*x), Q_quad, 'UniformOutput', false),                               ...
                                mat2cell(reshape(C_quad', prod(size(C_quad)) ,[]), nx*ones(neq_quad+niq_quad,1)),   ...
                                num2cell(b_quad,2),                                                                 ...
                                cellstr(char([double('=')*ones(neq_quad,1); double('<')*ones(niq_quad,1)]))], ...
                              [  "Qc"  , ...
                                 "q"   , ...
                                 "rhs" , ...
                                 "sense"  ],                                                                  ...
                              2);
end
m.A         = A_lin;
m.rhs       = b_lin;
m.sense     = char([ double('=')*ones(1,neq_lin) double('<')*ones(1,niq_lin) ]);
m.lb        = xmin;
m.ub        = xmax;
m.obj       = b;

%% Call the solver
isemptyQ = cell2mat(cellfun(@(x)(isempty(x) || ~any(any(x))), Q_quad, 'UniformOutput', false));
if sum(isemptyQ) == (neq_quad + niq_quad)   %% No quadratic terms in quadratic constraints (linear constraints)
    if isempty(H) || ~any(any(H))
        lpqcqp = 'LP';
    else
        lpqcqp = 'QP';
        if ~issparse(H)
            H = sparse(H);
        end
        m.Q = 0.5 * H;
    end
else
    lpqcqp = 'QCQP';
    if ~isempty(H) || any(any(H))
        if ~issparse(H)
            H = sparse(H);
        end
        m.Q = 0.5 * H;
    end
end
if verbose
    alg_names = {
        'automatic',
        'primal simplex',
        'dual simplex',
        'interior point',
        'concurrent',
        'deterministic concurrent',
        'deterministic concurrent simplex'
    };
    vn = gurobiver;
    fprintf('Gurobi Version %s -- %s %s solver\n', ...
        vn, alg_names{g_opt.Method+2}, lpqcqp);
end
results = gurobi(m, g_opt);

%% Check for status of the optimization run and prepare output
switch results.status
    case 'LOADED'           %% 1
        eflag = -1;
    case 'OPTIMAL'          %% 2, optimal solution found
        eflag = 1;
    case 'INFEASIBLE'       %% 3
        eflag = -3;
    case 'INF_OR_UNBD'      %% 4
        eflag = -4;
    case 'UNBOUNDED'        %% 5
        eflag = -5;
    case 'CUTOFF'           %% 6
        eflag = -6;
    case 'ITERATION_LIMIT'  %% 7
        eflag = -7;
    case 'NODE_LIMIT'       %% 8
        eflag = -8;
    case 'TIME_LIMIT'       %% 9
        eflag = -9;
    case 'SOLUTION_LIMIT'   %% 10
        eflag = -10;
    case 'INTERRUPTED'      %% 11
        eflag = -11;
    case 'NUMERIC'          %% 12
        eflag = -12;
    case 'SUBOPTIMAL'       %% 13
        eflag = -13;
    case 'INPROGRESS'       %% 14
        eflag = -14;
    case 'USER_OBJ_LIMIT'   %% 15
        eflag = -15;
    case 'WORK_LIMIT'       %% 15
        eflag = -16;
    case 'MEM_LIMIT'        %% 16
        eflag = -17;
    otherwise
        eflag = 0;
end
if nargout > 3
    output = results;
end

%% check for empty results (in case optimization failed)
if ~isfield(results, 'x') || isempty(results.x)
    x = NaN(nx, 1);
    lam.lower   = NaN(nx, 1);
    lam.upper   = NaN(nx, 1);
else
    x = results.x;
    lam.lower   = zeros(nx, 1);
    lam.upper   = zeros(nx, 1);
end
if ~isfield(results, 'objval') || isempty(results.objval)
    f = NaN;
else
    f = results.objval;
end
if ~isfield(results, 'pi') || isempty(results.pi)
    pi  = NaN(length(m.rhs), 1);
else
    pi  = results.pi;
end
if ~isfield(results, 'rc') || isempty(results.rc)
    rc  = NaN(nx, 1);
else
    rc  = results.rc;
end

kl = find(rc > 0);   %% lower bound binding
ku = find(rc < 0);   %% upper bound binding
lam.lower(kl)   =  rc(kl);
lam.upper(ku)   = -rc(ku);

[mu_l, mu_u] = convert_lin_constraint_multipliers(-pi(1:neq_lin), -pi(neq_lin+(1:niq_lin)), ieq_lin, igt_lin, ilt_lin);

if ~isempty(Q_quad) 
    if ~isfield(results, 'qcpi') || isempty(results.qcpi)
        qcpi = NaN(length(m.quadcon), 1); % Current version of Gurobi (11.0.3) does not return multipliers for non-convex qcqp. See https://docs.gurobi.com/projects/optimizer/en/current/reference/attributes/constraintquadratic.html#qcpi
    else
        qcpi = results.qcpi;
    end
    [mu_l_quad, mu_u_quad] = convert_lin_constraint_multipliers(-qcpi(1:neq_quad), -qcpi(neq_quad+(1:niq_quad)), ieq_quad, igt_quad, ilt_quad);  % WGV: we reuse the function for linear constraints, but a 'convert_quad_constraint_multipliers' is intended here

    lambda = struct( ...
    'mu_l'      , mu_l, ...
    'mu_u'      , mu_u, ...
    'mu_l_quad' , mu_l_quad, ...
    'mu_u_quad' , mu_u_quad, ...
    'lower'     , lam.lower, ...
    'upper'     , lam.upper ...
     );
else
    lambda = struct( ...
    'mu_l'      , mu_l, ...
    'mu_u'      , mu_u, ...
    'lower'     , lam.lower, ...
    'upper'     , lam.upper ...
     );
end