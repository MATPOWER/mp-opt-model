function ef = pne_target_lam_event(cb_data, cx)
%PNE_TARGET_LAM_EVENT  Event function to detect a target lambda value
%   EF = PNE_TARGET_LAM_EVENT(CB_DATA, CX)
%
%   PNES_MASTER event function to detect the completion of the continuation
%   curve or another target value of lambda.
%
%   Inputs:
%       CB_DATA : struct of data for callback functions
%       CX : struct containing info about current point (continuation soln)
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

%% event function value is a scalar equal to: current lambda - target lambda
target = cb_data.opt.stop_at;
if ischar(target)       %% 'FULL' trace requested
    if cx.z(end) >= 0   %% before the nose point ...
        target = -1;    %% prevent early termination (e.g. itr 1 rollback to 0)
    else                %% after the nose point ...
        target = 0;     %% terminate at lam = 0
    end
end
ef = cx.x(end) - target;
