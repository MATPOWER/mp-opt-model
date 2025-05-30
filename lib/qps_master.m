function [x, f, eflag, output, lambda] = qps_master(H, c, A, l, u, xmin, xmax, x0, opt)
% qps_master - Quadratic Program Solver wrapper function.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QPS_MASTER(H, C, A, L, U, XMIN, XMAX, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = QPS_MASTER(PROBLEM)
%   A common wrapper function for various QP solvers.
%   Solves the following QP (quadratic programming) problem:
%
%       min 1/2 X'*H*X + C'*X
%        X
%
%   subject to
%
%       L <= A*X <= U       (linear constraints)
%       XMIN <= X <= XMAX   (variable bounds)
%
%   Inputs (all optional except H, C, A and L):
%       H : matrix (possibly sparse) of quadratic cost coefficients
%       C : vector of linear cost coefficients
%       A, L, U : define the optional linear constraints. Default
%           values for the elements of L and U are -Inf and Inf,
%           respectively.
%       XMIN, XMAX : optional lower and upper bounds on the
%           X variables, defaults are -Inf and Inf, respectively.
%       X0 : optional starting value of optimization vector X
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           alg ('DEFAULT') : determines which solver to use, can be either
%                   a string (new-style) or a numerical alg code (old-style)
%               'DEFAULT' : (or 0) automatic, first available of Gurobi,
%                       CPLEX, MOSEK, Opt Tbx (if MATLAB), HIGHS, GLPK (LPs
%                       only), BPMPD, MIPS
%               'MIPS'    : (or 200) MIPS, MATPOWER Interior Point Solver
%                        pure MATLAB implementation of a primal-dual
%                        interior point method, if mips_opt.step_control = 1
%                        (or alg=250) it uses MIPS-sc, a step controlled
%                        variant of MIPS
%               'BPMPD'   : (or 100) BPMPD_MEX
%               'CLP'     : CLP
%               'CPLEX'   : (or 500) CPLEX
%               'GLPK'    : GLPK, (LP problems only, i.e. empty H matrix)
%               'GUROBI'  : (or 700) Gurobi
%               'HIGHS'   : HiGHS, https://highs.dev
%               'IPOPT'   : (or 400) IPOPT, requires MEX interface to IPOPT
%                           solver, https://github.com/coin-or/Ipopt
%               'MOSEK'   : (or 600) MOSEK
%               'OSQP'    : OSQP, https://osqp.org
%               'OT'      : (or 300) Optimization Toolbox, QUADPROG or LINPROG
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           bp_opt      - options vector for BP
%           clp_opt     - options vector for CLP
%           cplex_opt   - options struct for CPLEX
%           glpk_opt    - options struct for GLPK
%           grb_opt     - options struct for GUROBI
%           highs_opt   - options struct for HIGHS
%           ipopt_opt   - options struct for IPOPT
%           knitro_opt  - options struct for KNITRO
%           linprog_opt - options struct for LINPROG
%           mips_opt    - options struct for QPS_MIPS
%           mosek_opt   - options struct for MOSEK
%           osqp_opt    - options struct for OSQP
%           quadprog_opt - options struct for QUADPROG
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, opt
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
%   Example: (problem from https://v8doc.sas.com/sashtml/iml/chap8/sect12.htm)
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
%   Copyright (c) 2010-2025, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MP-Opt-Model.
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
    if isfield(p, 'c'),     c = p.c;        else,   c = [];     end
    if isfield(p, 'H'),     H = p.H;        else,   H = [];     end
else                                %% individual args
    if nargin < 9
        opt = [];
        if nargin < 8
            x0 = [];
            if nargin < 7
                xmax = [];
                if nargin < 6
                    xmin = [];
                end
            end
        end
    end
end

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
            case 800
                alg = 'KNITRO';
            otherwise
                error('qps_master: %d is not a valid algorithm code', alg);
        end
    end
else
    alg = 'DEFAULT';
end
if ~isempty(opt) && isfield(opt, 'verbose') && ~isempty(opt.verbose)
    verbose = opt.verbose;
else
    verbose = 0;
end
if strcmp(alg, 'DEFAULT')
    if have_feature('gurobi')       %% use Gurobi by default, if available
        alg = 'GUROBI';
    elseif have_feature('cplex')    %% if not, then CPLEX, if available
        alg = 'CPLEX';
    elseif have_feature('mosek')    %% if not, then MOSEK, if available
        alg = 'MOSEK';
    elseif have_feature('quadprog') && have_feature('matlab')   %% if not, then Opt Tbx, if available in MATLAB
        alg = 'OT';
    elseif have_feature('highs')    %% if not, then HiGHS, if available
        alg = 'HIGHS';
    elseif have_feature('knitro')   %% if not, then Artelys Knitro, if available
        alg = 'KNITRO';
    elseif ~nnz(H) && have_feature('glpk')  %% if not, and
        alg = 'GLPK';               %% prob is LP (not QP), then GLPK, if available
    elseif have_feature('bpmpd')    %% if not, then BPMPD_MEX, if available
        alg = 'BPMPD';
    else                            %% otherwise MIPS
        alg = 'MIPS';
    end
end

%%----- call the appropriate solver  -----
switch alg
    case 'BPMPD'                    %% use BPMPD_MEX
        [x, f, eflag, output, lambda] = ...
            qps_bpmpd(H, c, A, l, u, xmin, xmax, x0, opt);

        if eflag == -99
            if verbose
                fprintf('         Retrying with QPS_MIPS solver ...\n\n');
            end
            %% save (incorrect) solution from BPMPD
            bpmpd = struct('x', x, 'f', f, 'eflag', eflag, ...
                            'output', output, 'lambda', lambda);
            opt.alg = 'MIPS';
            [x, f, eflag, output, lambda] = ...
                qps_master(H, c, A, l, u, xmin, xmax, x0, opt);
            output.bpmpd = bpmpd;
        end
    case 'CLP'
        [x, f, eflag, output, lambda] = ...
            qps_clp(H, c, A, l, u, xmin, xmax, x0, opt);
    case 'CPLEX'
        [x, f, eflag, output, lambda] = ...
            qps_cplex(H, c, A, l, u, xmin, xmax, x0, opt);
    case 'GLPK'
        [x, f, eflag, output, lambda] = ...
            qps_glpk(H, c, A, l, u, xmin, xmax, x0, opt);
    case 'GUROBI'
        [x, f, eflag, output, lambda] = ...
            qps_gurobi(H, c, A, l, u, xmin, xmax, x0, opt);
    case 'HIGHS'
        [x, f, eflag, output, lambda] = ...
            qps_highs(H, c, A, l, u, xmin, xmax, x0, opt);
    case 'IPOPT'
        [x, f, eflag, output, lambda] = ...
            qps_ipopt(H, c, A, l, u, xmin, xmax, x0, opt);
    case 'MIPS'
        %% set up options
        if ~isempty(opt) && isfield(opt, 'mips_opt') && ~isempty(opt.mips_opt)
            mips_opt = opt.mips_opt;
        else
            mips_opt = [];
        end
        mips_opt.verbose = verbose;

        %% call solver
        [x, f, eflag, output, lambda] = ...
            qps_mips(H, c, A, l, u, xmin, xmax, x0, mips_opt);
    case 'MOSEK'
        [x, f, eflag, output, lambda] = ...
            qps_mosek(H, c, A, l, u, xmin, xmax, x0, opt);
    case 'OSQP'
        [x, f, eflag, output, lambda] = ...
            qps_osqp(H, c, A, l, u, xmin, xmax, x0, opt);
    case 'OT'                   %% use QUADPROG or LINPROG from Opt Tbx ver 2.x+
        [x, f, eflag, output, lambda] = ...
            qps_ot(H, c, A, l, u, xmin, xmax, x0, opt);
    case 'KNITRO'
        [x, f, eflag, output, lambda] = ...
            qps_knitro(H, c, A, l, u, xmin, xmax, x0, opt);
    otherwise
        fcn = ['qps_' lower(alg)];
        if exist([fcn '.m']) == 2
            [x, f, eflag, output, lambda] = ...
                feval(fcn, H, c, A, l, u, xmin, xmax, x0, opt);
        else
            error('qps_master: ''%s'' is not a valid algorithm code', alg);
        end
end
if ~isfield(output, 'alg') || isempty(output.alg)
    output.alg = alg;
end
