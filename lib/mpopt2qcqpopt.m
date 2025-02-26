function qcqpopt = mpopt2qcqpopt(mpopt, model, alg)
% mpopt2qcqpopt - Create/modify qcqps_master options struct from ``mpopt``.
% ::
%
%   QCQPOPT = MPOPT2QCQPOPT(MPOPT, MODEL)
%   QCQPOPT = MPOPT2QCQPOPT(MPOPT, MODEL, ALG)
%
%   Uses a MATPOWER options struct, MPOPT, to create or modify a
%   QCQPS_MASTER options struct.
%
%   Inputs (default values in parentheses):
%       MPOPT : MATPOWER options struct
%       MODEL ('QCQP') : (optional) one of the following model types, required
%               for selection of solver in case ALG is 'DEFAULT' (solver
%               precedence for each model type list in parentheses):
%           'LP'   - linear program with all continuous variables
%                   (GUROBI, CPLEX, MOSEK, OT (if MATLAB), GLPK, BPMPD, MIPS)
%           'QP'   - quadratic program with all continuous variables
%                   (GUROBI, CPLEX, MOSEK, OT (if large-scale alg available),
%                    BPMPD, MIPS)
%           'MILP' - LP with mixed integer/continuous variables
%                   (GUROBI, CPLEX, MOSEK, OT, GLPK)
%           'QCQP' - (default) QCQP with all continuous variables
%                   (GUROBI, IPOPT)
%       ALG ('opf.ac') : (optional) 'opf.ac', 'most', or any valid value of
%               OPT.alg for QCQPS_MASTER. The first option indicates 
%               that it should be taken from MPOPT.opf.ac.solver.
%
%   Output:
%       QCQPOPT : an options struct for use by QCQPS_MASTER and friends
%
% See also qcqps_master, mpoption.

%   MATPOWER
%   Copyright (c) 2019-2024, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% set default args
if nargin < 3
    alg = '';
    if nargin < 2
        model = '';
    end
end
if isempty(model)
    model = 'QCQP';
else
    model = upper(model);
end

%% get ALG from mpopt, if necessary
switch alg
    case {'opf.ac', ''}
        alg = upper(mpopt.opf.ac.solver);
    otherwise
        alg = upper(alg);
end

%% default solver
switch alg
    case {'DEFAULT', 0}
        if have_feature('knitro')
            alg = 'KNITRO';     %% use Artelys Knitro by default, if available
        elseif have_feature('gurobi')
            alg = 'GUROBI';     %% if not, then IPOPT, if available
        elseif have_feature('ipopt')
            alg = 'IPOPT';      %% if not, then IPOPT, if available
        else
            if model(1) ~= 'M'  %% not a mixed-integer problem
                alg = 'MIPS';   %% otherwise MIPS, if applicable
            else
                error('mpopt2qcqpopt: Sorry, no solver available for %s models', model);
            end
        end
end

%% create QCQPS_MASTER options struct
qcqpopt = struct('alg', alg, 'verbose', mpopt.verbose);
switch alg
    case {'MIPS', 200, 250}
        %% set up options
        qcqpopt.mips_opt = mpopt.mips;
        if qcqpopt.mips_opt.feastol == 0      %% = MPOPT.opf.violation by default
            qcqpopt.mips_opt.feastol = mpopt.opf.violation;
        end
    case {'IPOPT', 400}
        qcqpopt.ipopt_opt = ipopt_options([], mpopt);
    case {'GUROBI', 700}
        qcqpopt.grb_opt = gurobi_options([], mpopt);
    case {'KNITRO', 800}
        qcqpopt.knitro_opt = artelys_knitro_options([], mpopt);
    case {'OT', 300}
        if isfield(mpopt, 'linprog') && ~isempty(mpopt.linprog)
            qcqpopt.linprog_opt = mpopt.linprog;
        end
        if isfield(mpopt, 'quadprog') && ~isempty(mpopt.quadprog)
            qcqpopt.quadprog_opt = mpopt.quadprog;
        end
        if isfield(mpopt, 'intlinprog') && ~isempty(mpopt.intlinprog)
            qcqpopt.intlinprog_opt = mpopt.intlinprog;
        end
end