function [nx, cx, s] = pne_callback_nose(k, nx, cx, px, s, opt)
%PNE_CALLBACK_NOSE  Callback to handle NOSE events
%   [NX, CX, S] = PNE_CALLBACK_NOSE(K, NX, CX, PX, S, OPT)
%
%   Callback to handle NOSE events, triggered by event function
%   PNE_EVENT_NOSE to indicate the nose point of the continuation curve.
%
%   This function sets the msg field of the event when the nose point has
%   been found, raises the S.done flag and sets S.done_msg.
%
%   See PNE_CALLBACK_DEFAULT for details of the input and output arguments.

%   MP-Opt-Model
%   Copyright (c) 2016-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% skip if initialize, finalize or done
if k <= 0 || s.done
    return;
end

%% handle event
if ~s.rollback || nx.step == 0
    for i = 1:length(s.events)
        if strcmp(s.events(i).name, 'NOSE') && s.events(i).zero
            if nx.step == 0
                msg = ...
                    sprintf('Nose point eliminated by limit induced bifurcation at %d continuation steps, lambda = %.4g.', k, nx.x(end));
            else
                msg = ...
                    sprintf('Reached limit in %d continuation steps, lambda = %.4g.', k, nx.x(end));
            end

            %% the following conditional is only necessary if we also allow
            %% finding the location of the nose-point without terminating
            if ischar(opt.stop_at) && strcmp(opt.stop_at, 'NOSE') || ...
                ~ischar(opt.stop_at) && nx.x(end) < opt.stop_at
                s.done = 1;
                s.done_msg = msg;
            end
            break;
        end
    end
end
