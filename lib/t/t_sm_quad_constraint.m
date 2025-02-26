function om = t_sm_quad_constraint(solver)
% This is a simple example to test the new functionality of Mp-Opt-Model 
% that handles Quadratically Constrained Quadratic Problems (QCQP) via the
% qcqp master solver (qcqps_master)
%
% The problem to solve is as follows:
%
%  minimize     -x1^2 - 2 x2^2 - x3^2 - 0.5 x1 x2 - 0.5 x1 x3
%  subject to   8 x1 + 14 x2 + 7 x3 = 56
%               -x1^2 - x2^2 - x3^2 <= -25
%               x1 >= 0, x2 >= 0, x3 >= 0
%
%  The global minimum is at (0, 0, 8), with final objective = -64.0.
%  Depending on the starting point, solvers such as Knitro may converge to 
%  an alternate local solution at (7, 0, 0), with objective = -49.0.

%   MP-Opt-Model
%   Copyright (c) 2019-2023, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia
%
%   This file is part of MP-Opt-Model..
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

% Gather input
if nargin < 1
    solver = 'GUROBI';     % Default. Current options are: GUROBI and IPOPT
end


% Initialize optimization model object
om = opt_model();
om.init_set_types();
om.prob_type = 'QCQP';

% Prepare parameters for building objective function and constraints
n = 3;                                      % number of variables
H = sparse([-2 -1 -1; -1 -4 0; -1 0 -2]);   % quadratic term of objective
B = zeros(n,1);                             % linear term of objective
Q = sparse(-2*eye(n));                      % quadratic term of quadratic constraint
C = spalloc(1,n,0);                         % linear term of quadratic constraint
k = 0;                                      % constant term of quadratic constraint
l1 = - inf;                                 % lower bound of quadratic constraint
u1 = -25;                                   % upper bound of quadratic constraint
A = [8 14 7];                               % Coefficient matrix of linear constraint
l2 = 56;                                    % lower bound of linear constraint
u2 = 56;                                    % upper bound of linear constraint
xmin = zeros(n,1);                          % lower bound of decision variables
xmax = inf(n,1);                            % upper bound of decision variables

% Add optimization variables and their bounds
x0 = [2; 2; 2];                             % Initial point
om.var.add('X', n, x0, xmin, xmax);

% Add objective function
om.qdc.add(om.var, 'cost', H, B);

% Add quadratic constraint (inequality)
om.qcn.add(om.var, 'quad_cons', l1, u1, {Q}, C, k);

% Add linear constraint
om.lin.add(om.var, 'lin_cons', A, l2, u2);

% Define optimization options

opt = struct('alg', solver, 'verbose', 2);
switch solver
    case {'GUROBI'}
        
    case {'IPOPT'}
        opt.ipopt_opt = ipopt_options();
end

% Solve the problem (internal call to qcqps_master)
om.solve(opt);

% Print results
fprintf('============================================ \n')

if isfield(om.soln.output, 'status')
    fprintf('Problem status: %s \n', om.soln.output.status);
end

if om.soln.eflag
    fprintf('Optimal solution: \n'); disp(om.soln.x);
    fprintf('With objective: %2.3f \n', om.soln.f);
end
fprintf('============================================ \n')
