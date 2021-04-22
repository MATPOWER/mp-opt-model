function [nx, cx, s] = pne_callback_default(k, nx, cx, px, s, opt)
%PNE_CALLBACK_DEFAULT   Default callback function for PNES_MASTER
%   [NX, CX, S] = PNE_CALLBACK_DEFAULT(K, NX, CX, PX, S, OPT)
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
%           cbs - callback state, the user may add fields containing
%               any information the callback function would like to pass from
%               one invokation to the next, taking care not to step on fields
%               being used by other callbacks, such as the 'default' field
%               used by this default callback
%           ef - cell array of event function values
%       CX - current state (corresponding to most recent successful step)
%            (same structure as NX)
%       PX - previous state (corresponding to last successful step prior to CX)
%       S - struct for various flags, etc., with fields:
%           done - termination flag, 1 => terminate, 0 => continue
%           done_msg - char array with reason for termination
%           rollback - flag to indicate that the current step should be
%               rolled back and retried with a different step size, etc.
%           events - struct array listing events detected for this step,
%               see PNE_DETECT_EVENTS for details
%           results - output struct to be assigned to 'cont' field of
%               OUTPUT struct returned by PNES_MASTER
%       OPT - PNES_MASTER options struct
%
%   Outputs:
%       (all are updated versions of the corresponding input arguments)
%       NX - update this state if S.rollback is false,
%           e.g. user callback state ('cbs' field ), etc.
%       CX - update this state if S.rollback is true,
%           e.g. 'this_step' or 'this_parm'
%       S - struct for various flags, etc.
%           done - can request termination by setting to 1
%           done_msg - can set termination reason here
%           rollback - can request a rollback step, even if it was not
%               indicated by an event function
%           events - msg field for a given event may be updated
%           results - updated output struct to be assigned to 'cont' field
%               of OUTPUT struct returned by PNES_MASTER
%
%   This function is called in three different contexts, distinguished by
%   the value of K, as follows:
%   (1) initial - called with K = 0, after initial corrector step,
%           before 1st continuation step.
%   (2) iterations - called with K > 0, at each iteration, after
%           predictor-corrector step
%   (3) final - called with K < 0, after exiting predictor-corrector loop,
%           same as last iteration call, with negated K and updated PX, CX
%
%   User Defined PNE Callback Functions:
%       The user can define their own callback functions which take
%       the same form and are called in the same contexts as
%       PNE_CALLBACK_DEFAULT. These are specified via OPT.callbacks.
%       OPT.callbacks takes the same form as the MY_CBACKS input to
%       PNE_REGISTER_CALLBACKS.
%
%   See also PNES_MASTER, PNE_REGISTER_CALLBACKS.

%   MP-Opt-Model
%   Copyright (c) 2013-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% skip if rollback, except if it is a FINAL call
if s.rollback && k > 0
    return;
end

%% initialize variables
step    = nx.step;
x       = nx.x;
x_hat   = nx.x_hat;

%% get output function
if isempty(opt.output_fcn)
    output_fcn = @pne_output_fcn_default;
else
    output_fcn = opt.output_fcn;
end

%%-----  initialize/update state/results  -----
if k == 0       %% INITIAL call
    %% initialize state
    [names, vals] = output_fcn(x, x_hat);
    cbs = struct(   'iterations',   k, ...
                    'steps',        step, ...
                    'lam_hat',      x_hat(end), ...
                    'lam',          x(end)  );
    for j = 1:length(names)
        cbs.(names{j}) = vals{j};
    end
    cx.cbs.default = cbs;   %% update current callback state
    nx.cbs.default = cbs;   %% update next callback state
else
    cbs = nx.cbs.default;   %% get next callback state
    if k > 0    %% ITERATION call
        %% update state
        [names, vals] = output_fcn(x, x_hat);
        cbs.iterations = k;
        cbs.steps   = [ cbs.steps    step       ];
        cbs.lam_hat = [ cbs.lam_hat x_hat(end)  ];
        cbs.lam     = [ cbs.lam     x(end)      ];
        for j = 1:length(names)
            cbs.(names{j}) = [ cbs.(names{j})   vals{j} ];
        end
        nx.cbs.default = cbs;   %% update next callback state
    else            %% FINAL call
        %% assemble results struct
        names = output_fcn();
        r = s.results;
        r.steps         = cbs.steps;
        r.iterations    = -k;
        r.lam_hat       = cbs.lam_hat;
        r.lam           = cbs.lam;
        r.max_lam       = max(cbs.lam);
        for j = 1:length(names)
            r.(names{j}) = cbs.(names{j});
        end
        s.results = r;
    end
end

%%-----  plot continuation curve  -----
%% initialize plotting options
plt = opt.plot;
plot_idx_default = 0;
if plt.level && (k >= 0 || isempty(s.warmstart))
    %% set functions to use for horizontal and vertical coordinates
    if isempty(plt.xfcn)        %% default horizontal coord is last
        xf = @(x)x;             %% element of x (i.e. parameter lambda)
    else
        xf = plt.xfcn;
    end
    if isempty(plt.yfcn)        %% default vertical coord simply uses
        yf = @(x,idx)x(idx, :); %% element idx of x
    else
        yf = plt.yfcn;
        if isempty(plt.idx) && isempty(plt.idx_default)
            error('pne_callback_default: PNE plotting options require specifying either ''idx'' or ''idx_default'' when supplying custom ''yfcn''');
        end
    end

    if isempty(plt.idx) && ~isfield(cbs, 'plot_idx_default')   %% no index specified
        if isempty(plt.idx_default)
            idx = length(x) - 1;        %% last one before parameter lambda
        else
            idx = plt.idx_default();    %% get default from provided function
        end

        %% save it to keep it from changing in subsequent calls
        plot_idx_default = idx;
    else
        if isempty(plt.idx)
            idx = cbs.plot_idx_default; %% index, saved
        else
            idx = plt.idx;              %% index, provided
        end
    end
    nplots = length(idx);

    %% get plot data, initialize bounds for plot axes
    xx  = xf(cbs.(plt.xname));
    yy  = yf(cbs.(plt.yname), idx);
    xmin = 0;
    xmax = max(xx);
    ymin = min(min(yy));
    ymax = max(max(yy));
    if plt.level > 1
        xxh = xf(cbs.([plt.xname '_hat']));
        yyh = yf(cbs.([plt.yname '_hat']), idx);
        xmax = max(xmax, max(xxh));
        ymin = min(ymin, min(min(yyh)));
        ymax = max(ymax, max(max(yyh)));
    else
        xxh = [];
        yyh = [];
    end

    %% adjust bounds for plot axes
    if xmax < xmin + opt.step / 100;
        xmax = xmin + opt.step / 100;
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
            cx.cbs.default.plot_idx_default = plot_idx_default;
        end
        
        %% initialize continuation curve plot
        axis([xmin xmax ymin ymax]);
        plot(xx(:,1), yy(:,1), '-', 'Color', [0.25 0.25 1]);
        if nplots > 1
            plot_title = plt.title2;
        else
            plot_title = plt.title;
        end
        title(sprintf(plot_title, idx));
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
                plot([xx(:,k); xxh(:,k+1)], ...
                    [yy(kk, k); yyh(kk, k+1)], '-', ...
                    'Color', 0.85*[1 0.75 0.75]);
                %% correction line
                plot([xxh(:,k+1); xx(:,k+1)], ...
                    [yyh(kk, k+1); yy(kk, k+1)], '-', ...
                    'Color', 0.85*[0.75 1 0.75]);
                %% prediciton point
                plot(xxh(:,k+1), yyh(kk, k+1), 'x', ...
                    'Color', 0.85*[1 0.75 0.75]);
                %% correction point
                plot(xx(:,k+1), yy(kk, k+1), '-o', ...
                    'Color', [0.25 0.25 1]);
                drawnow;
            end
            if plt.level > 2
                pause;
            end
        end
    %%-----  FINAL call  -----
    else    % k < 0
        %% finish final continuation curve plot
        axis([xmin xmax ymin ymax]);
        %% curve of corrected points
        if isprop(gca, 'ColorOrderIndex')
            set(gca, 'ColorOrderIndex', 1); %% start over with color 1
        end
        hp = plot(xx', yy',  '-');
        if nplots > 1
            leg = cell(nplots, 1);
            for kk = 1:nplots
                leg{kk} = sprintf(plt.legend, idx(kk));
            end
            legend(hp, leg);
        end
        hold off;
    end
end


function [names, vals] = pne_output_fcn_default(x, x_hat)
%% [names, vals] = pne_output_fcn_default(x, x_hat)
%% names = pne_output_fcn_default()
names = {'x_hat', 'x'};
if nargin
    vals = {x_hat, x};
end
