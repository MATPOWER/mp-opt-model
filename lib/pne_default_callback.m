function [nx, cx, done, rollback, evnts, cb_data, results] = ...
    pne_default_callback(k, nx, cx, px, done, rollback, evnts, ...
                            cb_data, cb_args, results)
%PNE_DEFAULT_CALLBACK   Default callback function for PNES_MASTER
%   [NX, CX, DONE, ROLLBACK, EVNTS, CB_DATA, RESULTS] = 
%       PNE_DEFAULT_CALLBACK(K, NX, CX, PX, DONE, ROLLBACK, EVNTS, ...
%                               CB_DATA, CB_ARGS, RESULTS)
%
%   Default callback function used by PNES_MASTER that collects the resulst and
%   optionally, plots the nose curve. Inputs and outputs are defined below,
%   with the RESULTS argument being optional, used only for the final call
%   when K is negative.
%
%   Inputs:
%       K - continuation step iteration count
%       NX - next state (corresponding to proposed next step), struct with
%            the following fields:
%           x_hat - solution vector from predictor
%           x - solution vector from corrector
%           default_step - default step size
%           default_parm - default parameterization
%           this_step - step size for this step only
%           this_parm - paramterization for this step only
%           step - current step size
%           parm - current parameterization
%           events - struct array, event log
%           cb - user state, for callbacks (replaces CB_STATE), the user may
%               add fields containing any information the callback function
%               would like to pass from one invokation to the next, taking
%               care not to step on fields being used by other callbacks,
%               such as the 'default' field used by this default callback
%           ef - cell array of event function values
%       CX - current state (corresponding to most recent successful step)
%            (same structure as NX)
%       PX - previous state (corresponding to last successful step prior to CX)
%       DONE - struct, with flag to indicate termination of numerical
%           continuation and reason, with fields:
%           flag - termination flag, 1 => terminate, 0 => continue
%           msg - string containing reason for termination
%       ROLLBACK - scalar flag to indicate that the current step should be
%           rolled back and retried with a different step size, etc.
%       EVNTS - struct array listing any events detected for this step,
%           see CPF_DETECT_EVENTS for details
%       CB_DATA - struct containing potentially useful "static" data,
%           with the following fields (all based on internal indexing):
%       CB_ARGS - arbitrary data structure containing callback arguments
%       RESULTS - initial value of output struct to be assigned to
%           CONT field of OUTPUT struct returned by PNES_MASTER
%
%   Outputs:
%       (all are updated versions of the corresponding input arguments)
%       NX - user state ('cb' field ) should be updated here if ROLLBACK
%           is false
%       CX - may contain updated 'this_step' or 'this_parm' values to be used
%           if ROLLBACK is true
%       DONE - callback may have requested termination and set the msg field
%       ROLLBACK - callback can request a rollback step, even if it was not
%           indicated by an event function
%       EVNTS - msg field for a given event may be updated
%       CB_DATA - this data should only be modified if the underlying problem
%           has been changed (e.g. fcn has been altered) and should always
%           be followed by a step of zero length, i.e. set NX.this_step to 0
%           It is the job of any callback modifying CB_DATA to ensure that
%           all data in CB_DATA is kept consistent.
%       RESULTS - updated version of RESULTS input arg
%
%   This function is called in three different contexts, distinguished by
%   the value of K, as follows:
%   (1) initial - called with K = 0, without RESULTS input/output args,
%           after base power flow, before 1st continuation step.
%   (2) iterations - called with K > 0, without RESULTS input/output args,
%           at each iteration, after predictor-corrector step
%   (3) final - called with K < 0, with RESULTS input/output args, after
%           exiting predictor-corrector loop, inputs identical to last
%           iteration call, except K which is negated
%
%   User Defined PNE Callback Functions:
%       The user can define their own callback functions which take
%       the same form and are called in the same contexts as
%       PNE_DEFAULT_CALLBACK. These are specified via the MATPOWER
%       option 'pne.user_callback'. This option can be a string containing
%       the name of the callback function, or a struct with the following
%       fields, where all but the first are optional:
%           'fcn'       - string with name of callback function
%           'priority'  - numerical value specifying callback priority
%                (default = 20, see CPF_REGISTER_CALLBACK for details)
%           'args'      - arbitrary value (any type) passed to the callback
%                         as CB_ARGS each time it is invoked
%       Multiple user callbacks can be registered by assigning a cell array
%       of such strings and/or structs to the 'pne.user_callback' option.
%
%   See also PNES_MASTER, CPF_REGISTER_CALLBACK.

%   MP-Opt-Model
%   Copyright (c) 2013-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% skip if rollback, except if it is a FINAL call
if rollback && k > 0
    return;
end

%% initialize variables
step    = nx.step;
x       = nx.x;
x_hat   = nx.x_hat;

%%-----  initialize/update state/results  -----
if k == 0       %% INITIAL call
    %% initialize state
    cxx = struct(   'x_hat',    x_hat, ...
                    'x',        x, ...
                    'steps',    step, ...
                    'iterations', 0     );
    nxx = cxx;
    cx.cb.default = cxx;    %% update current callback state
    nx.cb.default = nxx;    %% update next callback state
else
    nxx = nx.cb.default;    %% get next callback state
    if k > 0    %% ITERATION call
        %% update state
        nxx.x_hat   = [nxx.x_hat    x_hat];
        nxx.x       = [nxx.x        x];
        nxx.steps   = [nxx.steps    step];
        nxx.iterations = k;
        nx.cb.default = nxx;    %% update next callback state
    else            %% FINAL call
        %% assemble results struct
        results.x_hat       = nxx.x_hat;
        results.x           = nxx.x;
        results.steps       = nxx.steps;
        results.iterations  = -k;
        results.max_lam     = max(results.x(end, :));
    end
end

%%-----  plot continuation curve  -----
%% initialize plotting options
plt = cb_data.opt.plot;
plot_idx_default = 0;
if plt.level
    xf = plt.xfcn;
    yf = plt.yfcn;

    if isempty(plt.idx) && ~isfield(nxx, 'plot_idx_default')   %% no index specified
        idx = length(x) - 1;    %% last one before lambda
        idx_e = idx;            %% an external index
        
        %% save it to keep it from changing in subsequent calls
        plot_idx_default = idx_e;
    else
        if isempty(plt.idx)
            idx_e = nxx.plot_idx_default;   %% external index, saved
        else
            idx_e = plt.idx;                %% external index, provided
        end
        if 0%idx_e are all valid
            %kk = first invalid entry in idx_e
            error('pne_default_callback: %d is not a valid index for OPT.plot.idx', idx_e(kk(1)));
        end
        idx = idx_e;    %%convert from idx_e;
        if any(idx == 0)
            kk = find(idx == 0);
            error('pne_default_callback: %d is not a valid index for OPT.plot.idx', idx_e(kk(1)));
        end
    end
    nplots = length(idx_e);

    %% set bounds for plot axes
    xmin = 0;
    xmax = max([max(xf(nxx.x_hat)); max(xf(nxx.x))]);
    ymin = min([min(min(yf(nxx.x_hat, idx))); min(min(yf(nxx.x, idx)))]);
    ymax = max([max(max(yf(nxx.x_hat, idx))); max(max(yf(nxx.x, idx)))]);
    if xmax < xmin + cb_data.opt.step / 100;
        xmax = xmin + cb_data.opt.step / 100;
    end
    if ymax - ymin < 2e-5;
        ymax = ymax + 1e-5;
        ymin = ymin - 1e-5;
    end
    xmax = xmax * 1.05;
    ymax = ymax + 0.05 * (ymax-ymin);
    ymin = ymin - 0.05 * (ymax-ymin);

    %%-----  INITIAL call  -----
    if k == 0
        %% save default plot idx in the state so we don't have to detect it
        %% each time, since we don't want it to change in the middle of the run
        if plot_idx_default
            cx.cb.default.plot_idx_default = plot_idx_default;
        end
        
        %% initialize continuation curve plot
        axis([xmin xmax ymin ymax]);
        plot(xf(cxx.x_hat(:,1)), yf(cxx.x_hat(:,1), idx), '-', 'Color', [0.25 0.25 1]);
        if nplots > 1
            plot_title = plt.title2;
        else
            plot_title = plt.title;
        end
        title(sprintf(plot_title, idx_e));
        xlabel(plt.xlabel);
        ylabel(plt.ylabel);
        hold on;
    %%-----  ITERATION call  -----
    elseif k > 0
        %% plot single step of the continuation curve
        if plt.level > 1
            axis([xmin xmax ymin ymax]);
            for kk = 1:nplots
                %% prediction line
                plot([xf(nxx.x(:,k)); xf(nxx.x_hat(:,k+1))], ...
                    [yf(nxx.x(:, k), idx(kk)); yf(nxx.x_hat(:,k+1), idx(kk))], '-', ...
                    'Color', 0.85*[1 0.75 0.75]);
                %% correction line
                plot([xf(nxx.x_hat(:,k+1)); xf(nxx.x(:,k+1))], ...
                    [yf(nxx.x_hat(:,k+1),idx(kk)); yf(nxx.x(:,k+1),idx(kk))], '-', ...
                    'Color', 0.85*[0.75 1 0.75]);
                %% prediciton point
                plot(xf(nxx.x_hat(:,k+1)), yf(nxx.x_hat(:,k+1),idx(kk)), 'x', ...
                    'Color', 0.85*[1 0.75 0.75]);
                %% correction point
                plot(xf(nxx.x(:,k+1))', yf(nxx.x(:,k+1),idx(kk))', '-o', ...
                    'Color', [0.25 0.25 1]);
                drawnow;
            end
            if plt.level > 2
                pause;
            end
        end
    %%-----  FINAL call  -----
    else    % k < 0
        %% finish final lambda-V nose curve plot
        axis([xmin xmax ymin ymax]);
        %% curve of corrected points
        if isprop(gca, 'ColorOrderIndex')
            set(gca, 'ColorOrderIndex', 1); %% start over with color 1
        end
        hp = plot(xf(nxx.x)', yf(nxx.x, idx)',  '-');
        if nplots > 1
            leg = cell(nplots, 1);
            for kk = 1:nplots
                leg{kk} = sprintf(plt.legend, idx_e(kk));
            end
            legend(hp, leg);
        end
        hold off;
    end
end
