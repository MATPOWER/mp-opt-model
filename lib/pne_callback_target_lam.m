function [nx, cx, s] = pne_callback_target_lam(k, nx, cx, px, s, opt)
%PNE_CALLBACK_TARGET_LAM  Callback to handle TARGET_LAM events
%   [NX, CX, S] = PNE_CALLBACK_TARGET_LAM(K, NX, CX, PX, S, OPT)
%
%   Callback to handle TARGET_LAM events, triggered by event function
%   PNE_EVENT_TARGET_LAM to indicate that a target lambda value has been
%   reached or that the full continuation curve has been traced.
%
%   This function sets the msg field of the event when the target lambda has
%   been found, raises the S.done flag and sets S.done_msg. If the current
%   or predicted next step overshoot the target lambda, it adjusts the step
%   size to be exactly what is needed to reach the target, and sets the
%   parameterization for that step to be the natural parameterization.
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

%% make stop_at numerical, 0 = FULL, -Inf = NOSE
target = opt.stop_at;
if ischar(target)
    if strcmp(target, 'FULL')
        target = 0;
    else                %% NOSE
        target = -Inf;
    end
end

%% handle event
event_detected = 0;
for i = 1:length(s.evnts)
    if strcmp(s.evnts(i).name, 'TARGET_LAM')
        event_detected = 1;
        if s.evnts(i).zero  %% prepare to terminate
            s.done = 1;
            if target == 0      %% FULL
                s.done_msg = sprintf('Traced full continuation curve in %d continuation steps', k);
            else                %% target lambda value
                s.done_msg = sprintf('Reached desired lambda %g in %d continuation steps', ...
                    target, k);
            end
        else                    %% set step-size & parameterization to terminate next time
            cx.this_parm = @pne_pfcn_natural;   %% change to natural parameterization
            if target == 0      %% FULL
                cx.this_step = cx.x(end);
                s.evnts(i).msg = sprintf('%s\n  step %d to overshoot full trace, reduce step size and set natural param', s.evnts(i).msg, k);
            else                %% target lambda value
                cx.this_step = target - cx.x(end);
                s.evnts(i).msg = sprintf('%s\n  step %d to overshoot target lambda, reduce step size and set natural param', s.evnts(i).msg, k);
            end
        end
        break;
    end
end

%% otherwise, check if predicted lambda of next step will overshoot
%% (by more than remaining distance, to play it safe)
if ~event_detected && ~s.rollback
    if isempty(nx.this_step)
        step = nx.default_step;
    else
        step = nx.this_step;
    end

    %% predictor step
    x_hat = nx.x + step * nx.z;

    if target == 0          %% FULL
        if x_hat(end) < -nx.x(end)
            nx.this_step = nx.x(end);
            nx.this_parm = @pne_pfcn_natural;   %% change to natural parameterization
            if opt.verbose > 3
                fprintf('  step %d prediction to overshoot full trace, set next step to natural param w/reduced size\n', k+1);
            end
        end
    elseif target > 0       %% target lambda value
        if x_hat(end) > target + (target - nx.x(end))
            nx.this_step = target - nx.x(end);
            nx.this_parm = @pne_pfcn_natural;   %% change to natural parameterization
            if opt.verbose > 3
                fprintf('  step %d prediction to overshoot target lambda, set next step to natural param w/reduced size\n', k+1);
            end
        end
    end
end
