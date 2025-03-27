function [x, f, eflag, output, lambda] = qcqps_knitro(H, b, Q, C, l1, u1, A, l2, u2, xmin, xmax, x0, opt)
% qcqps_knitro - Quadratically Constrained Quadratic Program Solver based on 
% Artelys KNITRO.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QCQPS_KNITRO(H, B, Q, C, K, L1, U1, A, L2, U2, XMIN, XMAX, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = QCQPS_KNITRO(PROBLEM)
%   A wrapper function providing a standardized interface for using
%   KNITRO to solve the following (possibly non-convex) QCQP (quadratically 
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
%       C : matrix (posibly sparse) of linear term of quadratic constraints
%       K : vector of constant terms of quadratic constraints
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
% See also qcqps_master, artelys_knitro_options, qcqp_knitro

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
    if isfield(p, 'C'),     C = p.C;        else,   C = [];     end
    if isfield(p, 'Q'),     Q = p.Q;        else,   Q = {};     end
    if isfield(p, 'b'),     b = p.b;        else,   b = [];     end
    if isfield(p, 'H'),     H = p.H;        else,   H = [];     end
else                                %% individual args
    if nargin < 13
        opt = [];
        if nargin < 12
            x0 = [];
            if nargin < 11
                xmax = [];
                if nargin < 10
                    xmin = [];
                    if nargin < 7
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
        if isempty(C) && isempty(A) && isempty(xmin) && isempty(xmax)
            error('qcqps_knitro: LP problem must include constraints or variable bounds');
        else
            if ~isempty(A) && ~isempty(C)
                if size(A,2) == size(C,2)
                    nx = size(A, 2);
                else
                    error('qcqp_knitro: number of columns of A and C must agree')
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
                error('qcqp_knitro: Dimensions of matrices H and Q{i}, i=1,2,...,%d must agree.', nrowQ)
            end
            nx = size(H, 1);
        else
            if abs(sum((1/nrowQ)*size_Q(:,2)) - length(b))  > 1e-10
                error('qcqp_knitro: Dimensions of matrices Q{i}, i=1,2,...,%d and vector b must agree.', nrowQ)
            end
            nx = length(b);
        end
    else
        error('qcqp_knitro: Input argument Q must be an N x 1 cell array')
    end
end
if ~isempty(H) || any(any(H))
    if ~issparse(H)
        H = sparse(H);
    end
end
if isempty(b)
    b = zeros(nx, 1);
end
if isempty(A)
   A = spalloc(1, nx, 0);
end
nrowA = size(A, 1);                %% number of original linear constraints
if isempty(u2)                  %% By default, linear inequalities are ...
    u2 = Inf(nrowA, 1);             %% ... unbounded above and ...
end
if isempty(l2)
    l2 = -Inf(nrowA, 1);            %% ... unbounded below.
end
if isempty(xmin)                %% By default, optimization variables are ...
    xmin = -Inf(nx, 1);         %% ... unbounded below and ...
end
if isempty(xmax)
    xmax = Inf(nx, 1);          %% ... unbounded above.
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

%% set up options struct for Knitro
if ~isempty(opt) && isfield(opt, 'knitro_opt') && ~isempty(opt.knitro_opt)
    kn_opt = artelys_knitro_options(opt.knitro_opt);
else
    kn_opt = knitro_options;
end
if verbose > 1
    opt.outlev = 3;
    if verbose > 2
        opt.outlev = 4;
    end
else
    opt.outlev = 0;
end
if verbose
    alg_names = {
        'automatic',
        'interior point direct',
        'active set',
        'sequential QP'
    };
    vn = knitrover;
    fprintf('Artelys Knitro Version %s -- %s %s solver\n', ...
        vn, alg_names{kn_opt.Method+1}, lpqcqp);
end

%% split up quadratic constraints
[ieq_quad, igt_quad, ilt_quad, Qe, Ce, sc, Qi, Ci, ubi] = ...
    convert_quad_constraint(Q, C, l1, u1);

%% split up linear constraints
if ~issparse(A)
    A = sparse(A);
end
if issparse(b)
    b = full(b);
end

[ieq_lin, igt_lin, ilt_lin, Ae, be, Ai, bi] = ...
    convert_lin_constraint(A, l2, u2);

%% grab some dimensions and adjust constraints
neq_quad = length(ieq_quad);                       %% number of quadratic equalities
niq_quad = length(ilt_quad) + length(igt_quad);    %% number of quadratic inequalities
neq_lin = length(ieq_lin);                         %% number of linear equalities
niq_lin = length(ilt_lin) + length(igt_lin);       %% number of linear inequalities

Qi_quad  = vertcat(cell(niq_lin,1), Qi);    % both linear/quadratic constraints must be passed as quadratic constraints of the form:
Qeq_quad = vertcat(cell(neq_lin,1), Qe);    %
A_quad   = [Ai; Ci];                        %          1/2 * x' * Qi_quad * x + A_quad * x <= b_quad      (inequality constraints)
Aeq_quad = [Ae, Ce];                        %          1/2 * x' * Qeq_quad * x + Aeq_quad * x = be_quad   (equality constraints)
b_quad   = [bi; ubi];                       %           
beq_quad = [be; sc];                        %

%% Call the solver
[x, f, exitflag, output, Lambda] = ...
            knitro_qcqp(H, b, Qi_quad, A_quad, b_quad, Qeq_quad, Aeq_quad, beq_quad, xmin, xmax, x0, [], kn_opt);

%% Extract multipliers
[mu_l, mu_u] = convert_lin_constraint_multipliers(Lambda.eqlin(1:neq_lin), ...
                    Lambda.ineqlin(1:niq_lin), ieq_lin, igt_lin, ilt_lin);

if (neq_quad + niq_quad) > 0 
    [mu_l_quad, mu_u_quad] = convert_lin_constraint_multipliers(Lambda.eqlin(neq_lin+1:end), ...
                    Lambda.ineqlin(niq_lin+1:end), ieq_quad, igt_quad, ilt_quad);

    lambda = struct( ...
    'mu_l'      , mu_l, ...
    'mu_u'      , mu_u, ...
    'mu_l_quad' , mu_l_quad, ...
    'mu_u_quad' , mu_u_quad, ...
    'lower'     , -1 * Lambda.lower, ...
    'upper'     , Lambda.upper ...
     );
else
    lambda = struct( ...
    'mu_l'      , mu_l, ...
    'mu_u'      , mu_u, ...
    'lower'     , -1 * lam.lower, ...
    'upper'     , lam.upper ...
     );
end

if exitflag == 0
    eflag = exitflag + 1;
else
    eflag = exitflag;
end