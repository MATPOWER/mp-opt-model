function [p, dp] = pne_pfcn_arc_len(x, xp, step, zp)
%PNE_PFCN_ARC_LEN  Parameterized Nonlinear EQ natural parameterization fcn
%   [P, DP] = PNE_PFCN_ARC_LEN(X, XP, STEP)
%   [P, DP] = PNE_PFCN_ARC_LEN(X, XP, STEP, ZP)
%
%   Inputs:
%       X    : solution vector x (last element is continuation parameter lambda)
%       XP   : previous solution vector
%       STEP : continuation parameter step size
%       ZP   : normalized tangent vector at XP
%
%   Outputs:
%       P : value of parameterization function
%       DP : Jacobian of paramerization function (transpose of gradient)
%
%   See also PNE_PFCN_NATURAL, PNE_PFCN_PSEUDO_ARC_LEN.

%   MP-Opt-Model
%   Copyright (c) 2013-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% function
p = sum((x - xp).^2) - step^2;

%% derivative
if nargout > 1
    dp = 2 * [x - xp]';
    if x(end) == xp(end)    %% first step
        dp(end) = 1.0;      %% avoid singular Jacobian from dp = 0
    end
end
