function [x, f, eflag, output, lambda] = qcqps_knitro(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt)
% qcqps_knitro - Quadratically Constrained Quadratic Program Solver based on
% Artelys Knitro.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QCQPS_KNITRO(H, C, Q, B, LQ, UQ, A, L, U, XMIN, XMAX, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = QCQPS_KNITRO(PROBLEM)
%   A wrapper function providing a standardized interface for using
%   KNITRO to solve the following (possibly non-convex) QCQP (quadratically
%   constrained quadratic programming) problem:
%
%       min 1/2 X'*H*X + C'*X
%        X
%
%   subject to
%
%       LQ(i) <= 1/2 X'*Q{i}*X + B(i,:)*X <= UQ(i), i = 1,2,...,nq
%                           (quadratic constraints)
%       L <= A*X <= U       (linear constraints)
%       XMIN <= X <= XMAX   (variable bounds)
%
%   Inputs (all optional except H, C, Q, B, LQ, and UQ):
%       H : matrix (possibly sparse) of quadratic cost coefficients
%       C : vector of linear cost coefficients
%       Q : nq x 1 cell array of sparse quadratic matrices for quadratic constraints
%       B : matrix (possibly sparse) of linear term of quadratic constraints
%       LQ, UQ: define the lower an upper bounds on the quadratic constraints
%       A, L, U : define the optional linear constraints. Default
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
%           knitro_opt - options struct for Artelys Knitro, value in verbose
%                        overrides these options
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt
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
%           mu_lq - lower (left-hand) limit on quadratic constraints
%           mu_uq - upper (right-hand) limit on quadratic constraints
%           lower - lower bound on optimization variables
%           upper - upper bound on optimization variables
%
%   Calling syntax options:
%       [x, f, exitflag, output, lambda] = ...
%           qcqps_knitro(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt)
%
%       x = qcqps_knitro(H, c, Q, B, lq, uq)
%       x = qcqps_knitro(H, c, Q, B, lq, uq, A, l, u)
%       x = qcqps_knitro(H, c, Q, B, lq, uq, A, l, u, xmin, xmax)
%       x = qcqps_knitro(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0)
%       x = qcqps_knitro(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt)
%       x = qcqps_knitro(problem), where problem is a struct with fields:
%                       H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'c', are optional, and problem with
%                       linear costs must include constraints
%       x = qcqps_knitro(...)
%       [x, f] = qcqps_knitro(...)
%       [x, f, exitflag] = qcqps_knitro(...)
%       [x, f, exitflag, output] = qcqps_knitro(...)
%       [x, f, exitflag, output, lambda] = qcqps_knitro(...)
%
%   Example: (problem from https://docs.gurobi.com/projects/examples/en/current/examples/matlab/qcp.html)
%       H = [];
%       c = [-1;0;0];
%       Q = { sparse([2 0 0; 0 2 0; 0 0 -2]), ...
%             sparse([2 0 0; 0 0 -2; 0 -2 0]) };
%       B = zeros(2,3);
%       lq = [-Inf;-Inf];
%       uq = [0; 0];
%       A = [1 1 1];
%       l = 1;
%       u = 1;
%       xmin = zeros(3,1);
%       xmax = Inf(3,1);
%       x0 = zeros(3,1);
%       opt = struct('verbose', 2);
%       [x, f, s, out, lambda] = ...
%           qcqps_knitro(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt);
%
% See also qcqps_master, artelys_knitro_options, qcqp_knitro

%   MP-Opt-Model
%   Copyright (c) 2019-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%   and Ray Zimmerman, PSERC Cornell
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
    if isfield(p, 'u'),     u = p.u;        else,   u = [];     end
    if isfield(p, 'l'),     l = p.l;        else,   l = [];     end
    if isfield(p, 'A'),     A = p.A;        else,   A = [];     end
    if isfield(p, 'uq'),    uq = p.uq;      else,   uq = [];    end
    if isfield(p, 'lq'),    lq = p.lq;      else,   lq = [];    end
    if isfield(p, 'B'),     B = p.B;        else,   B = [];     end
    if isfield(p, 'Q'),     Q = p.Q;        else,   Q = {};     end
    if isfield(p, 'c'),     c = p.c;        else,   c = [];     end
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
                        l = [];
                        u = [];
                    end
                end
            end
        end
    end
end

%% define nx, set default values for missing optional inputs
if isempty(Q)
    if isempty(H) || ~any(any(H))
        if isempty(c)
            error('qcqps_knitro: H or c could be empty, but not both.')
        end
        if isempty(B) && isempty(A) && isempty(xmin) && isempty(xmax)
            error('qcqps_knitro: LP problem must include constraints or variable bounds.');
        else
            if ~isempty(A) && ~isempty(B)
                if size(A,2) == size(B,2)
                    nx = size(A, 2);
                else
                    error('qcqps_knitro: number of columns of A and B must agree.')
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
    if isempty(B) || (~isempty(B) && (isempty(lq) || all(lq == -Inf)) && ...
                                     (isempty(uq) || all(uq == Inf)))
        B = sparse(0,nx);       %% no lq & uq limits => no quadratic constraints
    end
    nrowB = size(B, 1);         %% number of original quadratic constraints
    if isempty(uq)              %% By default, quadratic inequalities are ...
        uq = Inf(nrowB, 1);     %% ... unbounded above and ...
    end
    if isempty(lq)
        lq = -Inf(nrowB, 1);    %% ... unbounded below.
    end
else
    [nrowQ, ncolQ] = size(Q);
    if isa(Q,'cell') && ncolQ == 1
        size_Q  = cell2mat(cellfun(@(x)(size(x)), Q, 'UniformOutput', false));
        if sum(size_Q(:)) ~= 2*size_Q(1,1)*nrowQ
            error('qcqps_knitro: All matrices Q{i}, i=1,...,%d must be square of the same size.', nrowQ)
        end
        if ~isempty(B)
            [nrowB, ncolB] = size(B);
            if nrowB ~= nrowQ
                error('qcqps_knitro: Dimension mismatch between rows of Q (%d) and B (%d).', nrowQ, nrowB)
            end
            if size_Q(1,1) ~= ncolB
                error('qcqps_knitro: Dimension mismatch between columns of Q{i} (%d) and B (%d).', size_Q(1,1), ncolB)
            end
        else
            B = sparse(nrowQ, size_Q(1,1));
        end
        if isempty(H) && isempty(c)
            error('qcqps_knitro: H or c could be empty, but not both.')
        else
            if ~isempty(H)
                if abs(sum((1/nrowQ)*size_Q(:)) - sum(size(H)))  > 1e-10
                    error('qcqps_knitro: Dimensions of matrices Q{i}, i=1,...,%d and H must agree.', nrowQ)
                end
            else
                H = sparse(size_Q(1,1), size_Q(1,2));
            end
            if ~isempty(c)
                if size_Q(1,1) ~= numel(c)
                    error('qcqps_knitro: Dimensions of matrices Q{i}, i=1,...,%d and vector c must agree.', nrowQ)
                end
            else
                c = sparse(1, size_Q(1,1));
            end            
        end
        if ~isempty(A)
            if size_Q(1,1) ~= size(A,2)
                error('qcqps_knitro: Dimensions of matrices Q{i}, i=1,...,%d and matrix A must agree.', nrowQ)
            end
        else
            A = sparse(0, size_Q(1,1));
        end
        if ~isempty(xmin)
            if size_Q(1,1) ~= numel(xmin)
                error('qcqps_knitro: Dimensions of matrices Q{i}, i=1,...,%d and lower bound must agree.', nrowQ)
            end
        end
        if ~isempty(xmax)
            if size_Q(1,1) ~= numel(xmax)
                error('qcqps_knitro: Dimensions of matrices Q{i}, i=1,...,%d and upper bound must agree.', nrowQ)
            end
        end
        nx = size_Q(1,1);
    else
        error('qcqps_knitro: Input argument Q must be column vector cell array.')
    end
end
if isempty(c)
    c = zeros(nx, 1);
end
nrowA = size(A, 1);         %% number of original linear constraints
if isempty(u)               %% By default, linear inequalities are ...
    u = Inf(nrowA, 1);      %% ... unbounded above and ...
end
if isempty(l)
    l = -Inf(nrowA, 1);     %% ... unbounded below.
end
if isempty(xmin)            %% By default, optimization variables are ...
    xmin = -Inf(nx, 1);     %% ... unbounded below and ...
end
if isempty(xmax)
    xmax = Inf(nx, 1);      %% ... unbounded above.
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
    kn_opt.outlev = 3;
    if verbose > 2
        kn_opt.outlev = 4;
    end
else
    kn_opt.outlev = 0;
end

if ~issparse(A)
    A = sparse(A);
end
if issparse(c)
    c = full(c);
end

%% split up quadratic constraints
if ~isempty(Q)
    [ieq_quad, igt_quad, ilt_quad, Qe, Be, de, Qi, Bi, di] = ...
        convert_quad_constraint(Q, B, lq, uq);
    %% grab some dimensions
    neq_quad = length(ieq_quad);                       %% number of quadratic equalities
    niq_quad = length(ilt_quad) + length(igt_quad);    %% number of quadratic inequalities
else
    Qi = {}; Qe = {};
    neq_quad = 0; niq_quad = 0;
end

%% split up linear constraints
if isempty(Q) && ~isempty(B)
    [ieq_lin, igt_lin, ilt_lin, Ae, be, Ai, bi] = convert_lin_constraint([A; B], [l; lq], [u; uq]);
else
    [ieq_lin, igt_lin, ilt_lin, Ae, be, Ai, bi] = convert_lin_constraint(A, l, u);
end
%% grab some dimensions
neq_lin = length(ieq_lin);                         %% number of linear equalities
niq_lin = length(ilt_lin) + length(igt_lin);       %% number of linear inequalities


%% Call the solver
isemptyQ = cell2mat(cellfun(@(x)(isempty(x) || ~any(any(x))), [Qe; Qi], 'UniformOutput', false));
if sum(isemptyQ) == (neq_quad + niq_quad)   %% No quadratic terms in quadratic constraints (linear constraints)
    if isempty(H) || ~any(any(H))
        lpqcqp = 'LP';
    else
        lpqcqp = 'QP';
        if ~issparse(H)
            H = sparse(H);
        end
    end
else
    lpqcqp = 'QCQP';
    if ~isempty(H) || any(any(H))
        if ~issparse(H)
            H = sparse(H);
        end
    end
end
if verbose
    alg_names = {
        'automatic',
        'interior point direct',
        'interior point conjugate gradient',
        'active set',
        'sequential QP'
        };
    vn = knitrover;
    fprintf('Artelys Knitro Version %s -- %s %s solver\n', ...
        vn, alg_names{kn_opt.algorithm+1}, lpqcqp);
end

switch lpqcqp
    case {'QCQP'}
        Qi_quad  = vertcat(cell(niq_lin,1), Qi);    % both linear/quadratic constraints must be passed as quadratic constraints of the form:
        Qeq_quad = vertcat(cell(neq_lin,1), Qe);    %
        A_quad   = [Ai; Bi];                        %          1/2 * x' * Qi_quad * x + A_quad * x <= b_quad      (inequality constraints)
        Aeq_quad = [Ae, Be];                        %          1/2 * x' * Qeq_quad * x + Aeq_quad * x = be_quad   (equality constraints)
        b_quad   = [bi; di];                        %
        beq_quad = [be; de];                        %

        [x, f, exitflag, output, Lambda] = ...
            knitro_qcqp(H, c, Qi_quad, A_quad, b_quad, Qeq_quad, Aeq_quad, beq_quad, ...
            xmin, xmax, x0, [], kn_opt);
    case {'QP'}
        [x, f, exitflag, output, Lambda] = ...
            knitro_qp(H, c, Ai, bi, Ae, be, xmin, xmax, x0, [], kn_opt);
    case {'LP'}
        [x, f, exitflag, output, Lambda] = ...
            knitro_lp(c, Ai, bi, Ae, be, xmin, xmax, x0, [], kn_opt);
end

%% Extract multipliers
[mu_l, mu_u] = convert_constraint_multipliers(Lambda.eqlin(1:neq_lin), ...
    Lambda.ineqlin(1:niq_lin), ieq_lin, igt_lin, ilt_lin);

if (neq_quad + niq_quad) > 0
    [mu_lq, mu_uq] = convert_constraint_multipliers( ...
        Lambda.eqlin(neq_lin+1:end), Lambda.ineqlin(niq_lin+1:end), ...
        ieq_quad, igt_quad, ilt_quad);

    lambda = struct( ...
    'mu_l'      , mu_l, ...
    'mu_u'      , mu_u, ...
    'mu_lq'     , mu_lq, ...
    'mu_uq'     , mu_uq, ...
    'lower'     , -1 * Lambda.lower, ...
    'upper'     , Lambda.upper ...
     );
else
    lambda = struct( ...
    'mu_l'      , mu_l, ...
    'mu_u'      , mu_u, ...
    'lower'     , -1 * Lambda.lower, ...
    'upper'     , Lambda.upper ...
     );
end

if exitflag == 0
    eflag = exitflag + 1;
else
    eflag = exitflag;
end