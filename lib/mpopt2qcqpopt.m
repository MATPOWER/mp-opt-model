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
%           'QCQP' - (default) QCQP with all continuous variables
%                   (GUROBI, IPOPT)
%       ALG ('opf.ac') : (optional) 'opf.ac', or any valid value of
%               OPT.alg for QCQPS_MASTER. The first option indicates
%               that it should be taken from MPOPT.opf.ac.solver.
%
%   Output:
%       QCQPOPT : an options struct for use by QCQPS_MASTER and friends
%
% See also qcqps_master, mpoption.

%   MP-Opt-Model
%   Copyright (c) 2015-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

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
            alg = 'MIPS';   %% otherwise MIPS, if applicable
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
end
