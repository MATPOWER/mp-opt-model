function [f, df, d2f] = qcqp_nlp_costfcn(x, H, b)
% qcqp_nlp_costfcn - Evaluates objective function, gradient and Hessian
% of a quadratically constrained quadratic program (QCQP) solved using a 
% solver with nonlinear programming features
% ::
%
%   [F, DF, D2F] = QCQP_NLP_COSTFCN(X, H, b)
%
%   Objective function evaluation routine, suitable for use with MIPS,
%   FMINCON, etc. Computes objective function value, gradient and Hessian
%   of the quadratic function:
%
%                           1/2 X'*H*X + B'*X
%
%   Inputs:
%     X : optimization vector
%     H : matrix (possibly sparse) of quadratic cost coefficients
%     B : vector of linear cost coefficients
%
%   Outputs:
%     F   : value of objective function
%     DF  : (optional) gradient of objective function (column vector)
%     D2F : (optional) Hessian of objective function (sparse matrix)
%
%   Examples:
%       f = qcqp_nlp_costfcn(x, H, b);
%       [f, df] = qcqp_nlp_costfcn(x, H, b);
%       [f, df, d2f] = qcqp_nlp_costfcn(x, H, b);
%
% See also qcqp_nlp_consfcn, qcqp_nlp_hessfcn, qcqps_master.

%   MP-Opt-Model
%   Copyright (c) 2019-2023, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%%  ----- check inputs -----
if isempty(H)
    nx = length(b);
    H = sparse(nx,nx);
end

%%  ----- evaluate objective function -----
if nargout == 3
    f = 1/2 * x' * H * x + b' * x;
    df = H * x + b;
    d2f = H;
elseif nargout == 2
    f = 1/2 * x' * H * x + b' * x;
    df = H * x + b;
else
    f = 1/2 * x' * H * x + b' * x;
end
