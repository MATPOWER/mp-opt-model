function ef = pne_nose_event(cx, opt)
%PNE_NOSE_EVENT  Event function to detect the nose point
%   EF = PNE_NOSE_EVENT(CX, OPT)
%
%   PNE_MASTER event function to detect the nose point of the continuation
%   curve, based on the sign of the lambda component of the tangent vector.
%   
%   Inputs:
%       CX : struct containing info about current point (continuation soln)
%       OPT - PNES_MASTER options struct
%
%   Outputs:
%       EF : event function value

%   MP-Opt-Model
%   Copyright (c) 2016-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% event function value is dlam, the last element of the
%% normalized tangent vector at the current soln
ef = cx.z(end);
