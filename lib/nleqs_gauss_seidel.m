function [x, f, eflag, output, J] = nleqs_gauss_seidel(fcn, x0, opt)
%NLEQS_GAUSS_SEIDEL  Nonlinear Equation Solver based on Gauss-Seidel method.
%   [X, F, EXITFLAG, OUTPUT] = NLEQS_GAUSS_SEIDEL(FCN, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT] = NLEQS_GAUSS_SEIDEL(PROBLEM)
%   A function providing a standardized interface for using Gauss-Seidel
%   method to solve the nonlinear equation f(x) = 0, beginning from a
%   starting point x0.
%
%   Inputs:
%       FCN : handle to function that evaluates the function f(x) to
%           be solved. Calling syntax for this function is:
%               f = FCN(x)
%       X0 : starting value, x0, of vector x
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           max_it (1000) - maximum number of iterations for Gauss-Seidel method
%           tol (1e-8) - tolerance on Inf-norm of f(x)
%           gs_opt - options struct for Gauss-Seidel method, with field:
%               x_update_fcn (required) - handle to function that performs
%                   the Gauss-Seidel update step, with the following calling
%                   syntax:  x = x_update_fcn(x, f);
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: fcn, x0, opt
%
%   Outputs (all optional, except X):
%       X : solution vector x
%       F : final function value, f(x)
%       EXITFLAG : exit flag
%           1 = converged
%           0 or negative values = solver specific failure codes
%       OUTPUT : output struct with the following fields:
%           alg - algorithm code of solver used ('GS')
%           iterations - number of iterations performed
%           hist - struct array with trajectories of the following:
%                   normf
%           message - exit message
%
%   Note the calling syntax is almost identical to that of FSOLVE from
%   MathWorks' Optimization Toolbox. The function for evaluating the
%   nonlinear function is identical.
%
%   Calling syntax options:
%       [x, f, exitflag, output] = nleqs_gauss_seidel(fcn, x0);
%       [x, f, exitflag, output] = nleqs_gauss_seidel(fcn, x0, opt);
%       x = nleqs_gauss_seidel(problem);
%               where problem is a struct with fields: fcn, x0, opt
%               and all fields except 'fcn' and 'x0' are optional
%       x = nleqs_gauss_seidel(...);
%       [x, f] = nleqs_gauss_seidel(...);
%       [x, f, exitflag] = nleqs_gauss_seidel(...);
%       [x, f, exitflag, output] = nleqs_gauss_seidel(...);
%
%   Example: (problem from Christi Patton Luks, https://www.youtube.com/watch?v=pJG4yhtgerg)
%       function f = f2(x)
%       f = [  x(1)^2 +   x(1)*x(2)   - 10;
%              x(2)   + 3*x(1)*x(2)^2 - 57  ];
%
%       function x = x_update_fcn2(x, f)
%       x(1) = sqrt(10 - x(1)*x(2));
%       x(2) = sqrt((57-x(2))/3/x(1));
%
%       problem = struct( ...
%           'fcn', @(x)f2(x), ...
%           'x0',  [0; 0], ...
%           'opt', struct( ...
%               'verbose', 2, ...
%               'gs_opt', struct( ...
%                   'x_update_fcn', @x_update_fcn2 )));
%       [x, f, exitflag, output] = nleqs_gauss_seidel(problem);
%
%   See also NLEQS_MASTER.

%   MP-Opt-Model
%   Copyright (c) 1996-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%%----- input argument handling  -----
%% gather inputs
if nargin == 1 && isstruct(fcn) %% problem struct
    p = fcn;
    fcn = p.fcn;
    x0 = p.x0;
    if isfield(p, 'opt'),   opt = p.opt;    else,   opt = [];   end
else                            %% individual args
    if nargin < 3
        opt = [];
    end
end
nx = size(x0, 1);           %% number of variables

%% set default options
opt0 = struct(  'verbose', 0, ...
                'max_it', 1000, ...
                'tol', 1e-8 );
if isempty(opt)
    opt = opt0;
end
if isfield(opt, 'verbose') && ~isempty(opt.verbose)
    verbose = opt.verbose;
else
    verbose = opt0.verbose;
end
if isfield(opt, 'max_it') && opt.max_it     %% not empty or zero
    max_it = opt.max_it;
else
    max_it = opt0.max_it;
end
if isfield(opt, 'tol') && opt.tol           %% not empty or zero
    tol = opt.tol;
else
    tol = opt0.tol;
end
if isfield(opt, 'gs_opt') && isfield(opt.gs_opt, 'x_update_fcn')
    x_update_fcn = opt.gs_opt.x_update_fcn;
else
    error('nleqs_gauss_seidel: required ''gs_opt.x_update_fcn'' option missing');
end

%% initialize
eflag = 0;
i = 0;
x = x0;
hist(max_it+1) = struct('normf', 0);

%% evaluate f(x0)
f = fcn(x);

%% check tolerance
normf = norm(f, inf);
if verbose > 1
    fprintf('\n it     max residual');
    fprintf('\n----  ----------------');
    fprintf('\n%3d     %10.3e', i, normf);
end
if normf < tol
    eflag = 1;
    msg = sprintf('Gauss-Seidel method converged in %d iterations.', i);
    if verbose > 1
        fprintf('\nConverged!\n');
    end
end

%% save history
hist(i+1).normf = normf;

%% do Gauss-Seidel iterations
while (~eflag && i < max_it)
    %% update iteration counter
    i = i + 1;

    %% update x
    x = x_update_fcn(x, f);

    %% evalute f(x)
    f = fcn(x);

    %% check for convergence
    normf = norm(f, inf);
    if verbose > 1
        fprintf('\n%3d     %10.3e', i, normf);
    end

    %% save history
    hist(i+1).normf = normf;

    if normf < tol
        eflag = 1;
        msg = sprintf('Gauss-Seidel method converged in %d iterations.', i);
    end
end
if eflag ~= 1
    msg = sprintf('Gauss-Seidel method did not converge in %d iterations.', i);
end
if verbose
    fprintf('\n%s\n', msg);
end
if nargout > 3
    output = struct('alg', 'GS', ...
                    'iterations', i, ...
                    'hist', hist(1:i+1), ...
                    'message', msg  );
    if nargout > 4
        J = [];
    end
end
