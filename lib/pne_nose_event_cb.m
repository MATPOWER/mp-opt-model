function [nx, cx, done, rollback, evnts, opt, results] = pne_nose_event_cb(...
        k, nx, cx, px, done, rollback, evnts, opt, results)
%PNE_NOSE_EVENT_CB  Callback to handle NOSE events
%   [NX, CX, DONE, ROLLBACK, EVNTS, OPT, RESULTS] = 
%       PNE_NOSE_EVENT_CB(K, NX, CX, PX, DONE, ROLLBACK, EVNTS, ...
%                               OPT, CB_ARGS, RESULTS)
%
%   Callback to handle NOSE events, triggered by event function
%   PNE_NOSE_EVENT to indicate the nose point of the continuation curve.
%
%   This function sets the msg field of the event when the nose point has
%   been found, raises the DONE.flag and sets the DONE.msg.
%
%   See PNE_DEFAULT_CALLBACK for details of the input and output arguments.

%   MP-Opt-Model
%   Copyright (c) 2016-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% skip if initialize, finalize or done
if k <= 0 || done.flag
    return;
end

%% handle event
if ~rollback || nx.step == 0
    for i = 1:length(evnts)
        if strcmp(evnts(i).name, 'NOSE') && evnts(i).zero
            if nx.step == 0
                evnts(i).msg = ...
                    sprintf('Nose point eliminated by limit induced bifurcation at %d continuation steps, lambda = %.4g.', k, nx.x(end));
            else
                evnts(i).msg = ...
                    sprintf('Reached limit in %d continuation steps, lambda = %.4g.', k, nx.x(end));
            end

            %% the following conditional is only necessary if we also allow
            %% finding the location of the nose-point without terminating
            if ischar(opt.stop_at) && strcmp(opt.stop_at, 'NOSE');
                done.flag = 1;
                done.msg = evnts(i).msg;
            end
            break;
        end
    end
end
