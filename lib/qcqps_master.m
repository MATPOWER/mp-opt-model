function [x, f, eflag, output, lambda] = qcqps_master(H, b, Q, C, k, l1, u1, A, l2, u2, xmin, xmax, x0, opt)
% qcqps_master - Quadratically Constrained Quadratic Program Solver wrapper function.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QCQPS_MASTER(H, B, Q, C, K, L1, U1, A, L2, U2, XMIN, XMAX, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = QCQPS_MASTER(PROBLEM)
%   A common wrapper function for various QCQP solvers.
%   Solves the following QCQP (quadratically constrained quadratic programming) problem:
%
%       min 1/2 X'*H*X + B'*X
%        X
%
%   subject to
%
%    L1(i) <= 1/2 X'*Q{i}*X + C(i,:)*X + K(i) <= U1(i),  i = 1,2,...,NQ   (quadratic constraints)
%                 L2 <= A*X <= U2                                         (linear constraints)
%                XMIN <= X <= XMAX                                        (variable bounds)
%
%   Inputs (all optional except H, B, Q, C, K, L1, and U1):
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
%           alg ('DEFAULT') : determines which solver to use
%               'DEFAULT' :  automatic, first available of Gurobi,
%                       KNITRO, IPOPT, MIPS
%               'MIPS'    : MIPS, MATPOWER Interior Point Solver
%                        pure MATLAB implementation of a primal-dual
%                        interior point method, if mips_opt.step_control = 1
%                        (or alg=250) it uses MIPS-sc, a step controlled
%                        variant of MIPS
%               'GUROBI'  : Gurobi
%               'KNITRO'  : Artelys Knitro, requires Artelys Knitro solver
%                           https://www.artelys.com/solvers/knitro/
%               'IPOPT'   : IPOPT, requires MEX interface to IPOPT
%                           solver, https://github.com/coin-or/Ipopt
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           grb_opt     - options struct for GUROBI
%           knitro_opt  - options struct for Artelys Knitro
%           ipopt_opt   - options struct for IPOPT
%           mips_opt    - options struct for QCQPS_MIPS
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
%           mu_l_quad - lower (left-hand) limit on quadratic constraints
%           mu_u_quad - upper (right-hand) limit on quadratic constraints
%           lower - lower bound on optimization variables
%           upper - upper bound on optimization variables
%
%   Note the calling syntax is almost identical to that of QUADPROG
%   from MathWorks' Optimization Toolbox. The main difference is that
%   the linear constraints are specified with A, L, U instead of
%   A, B, Aeq, Beq.
%
%   Calling syntax options:
%       [x, f, exitflag, output, lambda] = ...
%           qps_master(H, c, A, l, u, xmin, xmax, x0, opt)
%
%       x = qps_master(H, c, A, l, u)
%       x = qps_master(H, c, A, l, u, xmin, xmax)
%       x = qps_master(H, c, A, l, u, xmin, xmax, x0)
%       x = qps_master(H, c, A, l, u, xmin, xmax, x0, opt)
%       x = qps_master(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = qps_master(...)
%       [x, f] = qps_master(...)
%       [x, f, exitflag] = qps_master(...)
%       [x, f, exitflag, output] = qps_master(...)
%       [x, f, exitflag, output, lambda] = qps_master(...)
%
%   Example: (problem from from https://v8doc.sas.com/sashtml/iml/chap8/sect12.htm)
%       H = [   1003.1  4.3     6.3     5.9;
%               4.3     2.2     2.1     3.9;
%               6.3     2.1     3.5     4.8;
%               5.9     3.9     4.8     10  ];
%       c = zeros(4,1);
%       A = [   1       1       1       1;
%               0.17    0.11    0.10    0.18    ];
%       l = [1; 0.10];
%       u = [1; Inf];
%       xmin = zeros(4,1);
%       x0 = [1; 0; 0; 1];
%       opt = struct('verbose', 2);
%       [x, f, s, out, lambda] = qps_master(H, c, A, l, u, xmin, [], x0, opt);

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

%% check for quadratic constraints
is_qcqp = ~isempty(Q);

%% default options
if ~isempty(opt) && isfield(opt, 'alg') && ~isempty(opt.alg)
    alg = opt.alg;
    %% convert integer codes to string values
    if ~ischar(alg)
        switch alg
            case 0
                alg = 'DEFAULT';
            case 100
                alg = 'BPMPD';
            case 200
                alg = 'MIPS';
                opt.mips_opt.step_control = 0;
            case 250
                alg = 'MIPS';
                opt.mips_opt.step_control = 1;
            case 300
                alg = 'OT';
            case 400
                alg = 'IPOPT';
            case 500
                alg = 'CPLEX';
            case 600
                alg = 'MOSEK';
            case 700
                alg = 'GUROBI';
            otherwise
                error('qps_master: %d is not a valid algorithm code', alg);
        end
    end
else
    alg = 'DEFAULT';
end
if strcmp(alg, 'DEFAULT')
    if have_feature('ipopt')
        alg = 'IPOPT';
    elseif have_feature('gurobi')       %% use Gurobi by default, if available
        alg = 'GUROBI';
    elseif have_feature('knitromatlab')    %% if not, then Artelys Knitro, if available
        alg = 'KNITRO';
    elseif have_feature('fmincon') && have_feature('matlab')   %% if not, then Opt Tbx, if available in MATLAB
        alg = 'OT';
    else                            %% otherwise MIPS
        alg = 'MIPS';
    end
end

%%----- call the appropriate solver  -----
if is_qcqp
    switch alg
        case 'GUROBI'
            [x, f, eflag, output, lambda] = ...
                qcqps_gurobi(H, b, Q, C, k, l1, u1, A, l2, u2, xmin, xmax, x0, opt);
        case 'KNITRO'
            [x, f, eflag, output, lambda] = ...
                qcqps_knitro(H, b, Q, C, k, l1, u1, A, l2, u2, xmin, xmax, x0, opt);
        case {'MOSEK'}
            [x, f, eflag, output, lambda] = ...
                qcqps_mosek(H, b, Q, C, k, l1, u1, A, l2, u2, xmin, xmax, x0, opt);
        otherwise
            id_ = find(alg == '_');
            if ~isempty(id_)         %% QCQP treated as a general nonlinear program
                opt.alg = alg(1 : id_ - 1);  % name of nonlinear solver
            else
                opt.alg = alg;       %% Default solver
            end

            % compute parameters for constraints function evaluation
            [ieq_quad, igt_quad, ilt_quad, Qe, Ce, lbe, Qi, Ci, lbi] = convert_quad_constraint(Q, C, k, l1, u1);
            QQ = struct('blkQe', blkdiag(Qe{:}), 'blkQi', blkdiag(Qi{:}));
            CC = struct('Ce', Ce, 'Ci', Ci);
            bb = struct('be', lbe, 'bi', lbi);

            % compute parameters for Hessian evaluation
            matQi = cell2mat(Qi);
            matQe = cell2mat(Qe);

            %% run solver
            f_fcn = @(x)qcqp_nlp_costfcn(x, H, b);
            gh_fcn = @(x)qcqp_nlp_consfcn(x, QQ, CC, bb);
            hess_fcn = @(x, lambda, cost_mult)qcqp_nlp_hessfcn(x, lambda, H, matQi, matQe, cost_mult);
            [x, f, eflag, output, Lambda] = ...
                nlps_master(f_fcn, x0, A, l2, u2, xmin, xmax, gh_fcn, hess_fcn, opt);
            
            if ~isfield(Lambda, 'eqnonlin')
                Lambda.eqnonlin =  zeros(length(lbe), 1);
            end
            if ~isfield(Lambda, 'ineqnonlin')
                Lambda.ineqnonlin = zeros(length(lbi), 1);
            end
            
            % gather multipliers for quadratic constraints
            [mu_l_quad, mu_u_quad] = convert_lin_constraint_multipliers(Lambda.eqnonlin, Lambda.ineqnonlin, ieq_quad, igt_quad, ilt_quad);  % WGV: we reuse the function for linear constraints, but a 'convert_quad_constraint_multipliers' is intended here
            
            lambda = struct( ...
                'mu_l'      , Lambda.mu_l, ...
                'mu_u'      , Lambda.mu_u, ...
                'mu_l_quad' , mu_l_quad, ...
                'mu_u_quad' , mu_u_quad, ...
                'lower'     , Lambda.lower, ...
                'upper'     , Lambda.upper ...
                 );
    end
else
    [x, f, eflag, output, lambda] = qps_master(H, b, A, l2, u2, xmin, xmax, x0, opt);
end
if ~isfield(output, 'alg') || isempty(output.alg)
    output.alg = alg;
end