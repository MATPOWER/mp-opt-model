function varargout = pnes_master(fcn, x0, opt)
%PNES_MASTER  Parameterized Nonlinear Equation Solver wrapper function.
%
%   UNFINISHED - need to update all of the help text
%
%   Inputs:
%   [X, F, EXITFLAG, OUTPUT, JAC] = PNES_MASTER(FCN, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT, JAC] = PNES_MASTER(PROBLEM)
%   A common wrapper function for numerical continuation methods for
%   solving parameterized nonlinear equations. Traces the solutions of
%   a parameterized nonlinear equation f(x) = 0, beginning from a starting
%   point x0, where f(x) has dimension n and x has dimension n+1.
%
%   Inputs:
%       FCN : handle to function that evaluates the function f(x) to
%           be solved and its (optionally, depending on the selected
%           solver) Jacobian, J(x). Calling syntax for this function is:
%               f = FCN(x)
%               [f, J] = FCN(x)
%           For a parameterized function, f is n x 1, x is (n+1) x 1,
%           and J is the n x (n+1) matrix of partial derivatives of
%           f (rows) w.r.t. x (cols).
%       X0 : starting value, x0, of vector x ((n+1) x 1)
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           solve_base (1) : determines whether or not to run an
%               initial corrector stage on the base solution
%
%           should I put all of the nleqs_master opts in a
%               nleqs_opts field?
%
%           something_alg ('DEFAULT') : overall algorithm
%               'DEFAULT' : automatic, current default is ...
%
%           alg ('DEFAULT') : determines which solver to use for corrector
%               'DEFAULT' : automatic, current default is NEWTON
%               'NEWTON'  : standard, full-Jacobian Newton's method
%               'CORE'    : core algorithm, with arbitrary update function
%               'FD'      : fast-decoupled Newton's method
%               'FSOLVE'  : FSOLVE, MATLAB Optimization Toolbox
%               'GS'      : Gauss-Seidel method
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           max_it (0) - maximum number of iterations
%                       (0 means use solver's own default)
%           tol (0) - termination tolerance on f(x)
%                       (0 means use solver's own default)
%           core_sp - solver parameters struct for NLEQS_CORE, required
%               when alg = 'CORE' (see NLEQS_CORE for details)
%           fd_opt - options struct for fast-decoupled Newton, NLEQS_FD_NEWTON
%           fsolve_opt - options struct for FSOLVE
%           gs_opt - options struct for Gauss-Seidel method, NLEQS_GAUSS_SEIDEL
%           newton_opt - options struct for Newton's method, NLEQS_NEWTON
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
%           alg - algorithm code of solver used
%           (others) - algorithm specific fields
%       JAC : final Jacobian matrix, J(x)
%
%   Note the calling syntax is almost identical to that of FSOLVE from
%   MathWorks' Optimization Toolbox. The function for evaluating the
%   nonlinear function and Jacobian is identical.
%
%   Calling syntax options:
%       [x, f, exitflag, output, jac] = nleqs_master(fcn, x0);
%       [x, f, exitflag, output, jac] = nleqs_master(fcn, x0, opt);
%       x = nleqs_master(problem);
%               where problem is a struct with fields: fcn, x0, opt
%               and all fields except 'fcn' and 'x0' are optional
%       x = nleqs_master(...);
%       [x, f] = nleqs_master(...);
%       [x, f, exitflag] = nleqs_master(...);
%       [x, f, exitflag, output] = nleqs_master(...);
%       [x, f, exitflag, output, jac] = nleqs_master(...);
%
%   Example: (problem from https://www.chilimath.com/lessons/advanced-algebra/systems-non-linear-equations/)
%       function [f, J] = f1(x)
%       f = [  x(1)   + x(2) - 1;
%             -x(1)^2 + x(2) + 5    ];
%       if nargout > 1
%           J = [1 1; -2*x(1) 1];
%       end
%
%       problem = struct( ...
%           'fcn',    @(x)f1(x), ...
%           'x0',       [0; 0], ...
%           'opt',      struct('verbose', 2) ...
%       );
%       [x, f, exitflag, output, jac] = nleqs_master(problem);
%
%   See also PNES_...

%   MP-Opt-Model
%   Copyright (c) 2013-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell,
%   Shrirang Abhyankar, Argonne National Laboratory,
%   and Alexander Flueck, IIT
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
    if isfield(p, 'opt'),   opt = p.opt;    else,   opt = struct(); end
else                            %% individual args
    if nargin < 3
        opt = struct();
    end
end

%% default options
if isfield(opt, 'verbose') && opt.verbose > 4
    nleqs_opt_verbose = 2;
else
    nleqs_opt_verbose = 0;
end
dopts = struct( ...
    'alg',              'DEFAULT', ...  %% algorithm
    'verbose',          0, ...
    'nleqs_opt',        struct('verbose', nleqs_opt_verbose), ...
    'solve_base',       1, ...          %% UNFINISHED
    'parameterization', 3, ...          %% 1 - natural, 2 - arc len, 3 - pseudo arc len
    'stop_at',          'NOSE', ...     %% 'NOSE', 'FULL', <lam_stop>
    'step',             0.05, ...       %% UNFINISHED
    'step_min',         1e-4, ...       %% UNFINISHED
    'step_max',         0.2, ...        %% UNFINISHED
    'adapt_step',       0, ...          %% UNFINISHED
    'adapt_step_ws',    1, ...          %% UNFINISHED
    'adapt_step_damping', 0.7, ...      %% UNFINISHED
    'adapt_step_tol',   1e-3, ...       %% UNFINISHED
    'default_event_tol',1e-3, ...       %% UNFINISHED
    'target_lam_tol',   0, ...          %% UNFINISHED
    'nose_tol',         0, ...          %% UNFINISHED
    'events',           {{}}, ...       %% UNFINISHED
    'callbacks',        {{}}, ...       %% UNFINISHED
    'output_fcn',       [], ...         %% custom output fcn, default callback
    'warmstart',        [], ...         %% default warm start state
    'plot',             struct( ...     %% used by pne_callback_default() for plotting
        'level',        0, ...          %% 0 - no plot, 1 - final, 2 - steps, 3 - steps w/pause
        'idx',          [], ...         %% index of quantity to plot, passed to yfcn()
        'idx_default',  [], ...         %% fcn to provide default value for idx, if none provided
        'xfcn',         [], ...         %% UNFINISHED
        'yfcn',         [], ...         %% UNFINISHED
        'title',        'Value of Variable %d', ... %% UNFINISHED
        'title2',       'Value of Multiple Variables', ...  %% UNFINISHED
        'xname',        'lam', ...      %% name of output field holding x vals
        'yname',        'x', ...        %% name of output field holding y vals
        'xlabel',       '\lambda', ...              %% UNFINISHED
        'ylabel',       'Variable Value', ...       %% UNFINISHED
        'legend',       'Variable %d' ...           %% UNFINISHED
    ) ...
);
opt = nested_struct_copy(dopts, opt);
%% use opt.default_event_tol for NOSE and TARGET_LAM events, unless specified
if opt.target_lam_tol == 0
    opt.target_lam_tol = opt.default_event_tol;
end
if opt.nose_tol == 0
    opt.nose_tol = opt.default_event_tol;
end

%% initialize
warmstarted = ~isempty(opt.warmstart);
s = struct( ...         %% container struct for various variables, flags
    'done',     0, ...      %% flag indicating continuation has terminated
    'warmstart',[], ...     %% warm start state to return when done (to pass
                            ...%% to subsequent warm-started call to PNES_MASTER)
    'done_msg', '', ...     %% termination message
    'rollback', 0, ...      %% flag to indicate a step must be rolled back
    'evnts',    [], ...     %% struct array for detected events
    'results',  []  );      %% results struct

%% register event and callback functions
switch opt.stop_at
    case 'NOSE'
        my_events = {{'NOSE', @pne_event_nose, opt.nose_tol}, opt.events{:}};
        my_cbacks = {{@pne_callback_nose, 51}, opt.callbacks{:}};
    case 'FULL'
        my_events = {{'TARGET_LAM', @pne_event_target_lam, opt.target_lam_tol}, opt.events{:}};
        my_cbacks = {{@pne_callback_target_lam, 50}, opt.callbacks{:}};
    otherwise   %% numeric, stop at target lam or nose point, whichever is 1st
        my_events = {{'TARGET_LAM', @pne_event_target_lam, opt.target_lam_tol}, ...
                     {'NOSE', @pne_event_nose, opt.nose_tol}, ...
                        opt.events{:}};
        my_cbacks = {{@pne_callback_target_lam, 50}, ...
                     {@pne_callback_nose, 51}, ...
                        opt.callbacks{:}};
end
my_cbacks{end+1} = {@pne_callback_default, 0};
reg_ev = pne_register_events(my_events, opt);   %% registered event functions
reg_cb = pne_register_callbacks(my_cbacks);     %% registered callback functions
nef = length(reg_ev);   %% number of registered event functions
ncb = length(reg_cb);   %% number of registered callback functions

t0 = tic;                       %% start timing

%% initialize continuation step counter
if warmstarted
    cont_steps = opt.warmstart.cont_steps + 1;
    if opt.verbose
        fprintf('... CONTINUATION RESUMED\n');
    end
else
    cont_steps = 0;
    if opt.verbose
        v = mpomver('all');
        fprintf('\nMP-Opt-Model Version %s, %s', v.Version, v.Date);
        fprintf(' -- Predictor/Corrector Continuation Method\n');
    end
end

%% solve corrector step for base point
if opt.solve_base && ~warmstarted
    pfcn = @(xx)pne_pfcn_natural(xx, x0, 0);
    [x, f, exitflag, out] = nleqs_master(@(xx)pne_corrector_fcn(xx, fcn, pfcn), x0, opt.nleqs_opt);
    if exitflag
        if opt.verbose > 1
            fprintf('step %3d  :                          lambda = %6.3f, %2d corrector steps\n', cont_steps, x0(end), out.iterations);
        end
    else
        s.done = 1;
        s.done_msg = sprintf('base solution did not converge in %d iterations', out.iterations);
        if opt.verbose
            fprintf('%s\n', s.done_msg);
        end
    end
else
    x = x0;     %% ignored for warmstart, overwritten by warmstart.cx
end

%% initialize numerical continuation
if ~s.done
    locating = 0;   %% flag to indicate that an event interval was detected,
                    %% but the event has not yet been located
    rb_cnt_ef = 0;  %% counter for rollback steps triggered by event function intervals
    rb_cnt_cb = 0;  %% counter for rollback steps triggered directly by callbacks

    if warmstarted
        manual_direction_switch = 0;
        ws = opt.warmstart;
        step = 0;
        parm = ws.parm;
        direction = ws.direction;
        default_parm = ws.default_parm;
        default_step = ws.default_step;
        cbx = ws.cbx;
        evnts = ws.events;
        x = ws.x;
        z = ws.z;
        if manual_direction_switch
            %% decide whether to switch directions
            reply = input('Switch directions? Y/N [N]:','s');
            if strcmp(upper(reply), 'Y')
                direction = -direction;
            end
        elseif isfield(ws, 'dir_from_jac_eigs') && ws.dir_from_jac_eigs
            [~, J] = fcn(x);
            eigs_opt.tol = 1e-3;
            eigs_opt.maxit = 2*length(x);
            direction = sign(z(end) * ...
                        min(real(eigs(J(:,1:end-1), 1, 'SR', eigs_opt))));
        end

        if opt.adapt_step   %% hey, maybe slow down, things may have changed
            default_step = default_step * opt.adapt_step_ws;
        end
    else
        %% initialize parameterization function
        switch opt.parameterization
            case 1
                parm = @pne_pfcn_natural;           %% NAT
            case 2
                parm = @pne_pfcn_arc_len;           %% ARC
            case 3
                parm = @pne_pfcn_pseudo_arc_len;    %% PAL
            otherwise
                error('pnes_master: OPT.parameterization (= %d) must be 1, 2, or 3', parm);
        end

        %% finish initializing tangent vector
        direction = 1;  %% increasing lambda
        z0 = zeros(length(x0), 1); z0(end) = direction; %% direction of positive lambda
        z = pne_tangent(x, x, z0, fcn, parm, direction);

        step = opt.step;
        default_step = step;
        default_parm = parm;
        cbx = [];
        evnts = [];
    end

    %% initialize state for current continuation step
    cx = struct( ...        %% current state
        'x_hat',        x, ...      %% predicted solution value
        'x',            x, ...      %% corrected solution value
        'z',            z, ...      %% normalized tangent vector
        'default_step', default_step, ...   %% default step size
        'default_parm', default_parm, ...   %% default parameterization
        'this_step', [], ...        %% step size for this step only
        'this_parm', [], ...        %% parameterization for this step only
        'step', step, ...           %% current step size
        'parm', parm, ...           %% current parameterization
        'events', evnts, ...        %% event log
        'cb', cbx, ...              %% user state, for callbacks
        'ef', [] ...                %% event function values
    );

    %% initialize event function values
    cx.ef = cell(nef, 1);
    for k = 1:nef
        cx.ef{k} = reg_ev(k).fcn(cx, opt);
    end

    if warmstarted
        %% initialize state for previous continuation step
        px = cx;
        px.x = ws.xp;   %% use warm start value for solution value
        px.z = ws.zp;   %% use warm start value for tangent
    else
        %% invoke callbacks - "initialize" context
        for k = 1:ncb
            [nx, cx, s] = reg_cb(k).fcn(cont_steps, cx, cx, cx, s, opt);
        end
        cont_steps = cont_steps + 1;

        %% check for case with base and target the same
        if opt.solve_base
            fb = f(1:end-1);
        else
            fb = fcn(x0);
            exitflag = 1;
        end
        xt = x0;
        xt(end) = 1;
        ft = fcn(xt);
        if norm(fb - ft, Inf) < 1e-12
            s.done = 1;
            s.done_msg = 'base and target functions are identical';
        end

        %% initialize state for previous continuation step
        px = cx;
    end
end

%%-----  run numerical continuation  -----
while ~s.done
    %% initialize next candidate with current state
    nx = cx;

    %% predictor step
    nx.x_hat = cx.x + cx.step * cx.z;

    %% corrector step
    pfcn = @(xx)cx.parm(xx, cx.x, cx.step, cx.z);
    [nx.x, f, exitflag, out] = nleqs_master(@(xx)pne_corrector_fcn(xx, fcn, pfcn), nx.x_hat, opt.nleqs_opt);
    if ~exitflag    %% corrector failed
        s.done = 1;
        s.done_msg = sprintf('Corrector did not converge in %d iterations.', out.iterations);
        if opt.verbose
            fprintf('step %3d  : %s stepsize = %-9.3g lambda = %6.3f  corrector did not converge in %d iterations\n', cont_steps, pne_ptag(cx.parm), cx.step, nx.x(end), out.iterations);
        end
        cont_steps = max(cont_steps - 1, 1);    %% go back to last step, but not to 0
        break;
    end

    %% compute new tangent direction, based on current or prev state: tx
    if nx.step == 0     %% if this is a re-do step, cx and nx are the same
        tx = px;            %% so use px as the previous state
    else                %% otherwise
        tx = cx;            %% use cx as the previous state
    end
    nx.z = pne_tangent(nx.x, tx.x, tx.z, fcn, nx.parm, direction);

    %% detect events
    for k = 1:nef
        nx.ef{k} = reg_ev(k).fcn(nx, opt);  %% update event functions
    end
    [s.rollback, s.evnts, nx.ef] = pne_detect_events(reg_ev, nx.ef, cx.ef, nx.step);

    %% adjust step-size to locate event function zero, if necessary
    if s.rollback               %% current step overshot
        %% rollback and initialize next step size based on rollback and previous
        rx = nx;                    %% save state we're rolling back from
        rx_evnts = s.evnts;         %% and critical event info
        cx.this_step = s.evnts.step_scale * rx.step;
        cx.this_parm = rx.parm;     %% keep same parameterization as last step
        locating = 1;               %% enter "locating" mode (or stay in it)
        rb_cnt_ef = rb_cnt_ef + 1;  %% increment rollback counter for ef intervals
        if rb_cnt_ef > 26
            s.done = 1;
            s.done_msg = sprintf('Could not locate %s event!', s.evnts.name);
        end
        if opt.verbose > 3
            loc_msg = sprintf('OVERSHOOT  : f = [%g, <<%g>>], step <-- %.4g', ...
                        cx.ef{s.evnts.eidx}(s.evnts.idx(1)), ...
                        rx.ef{s.evnts.eidx}(s.evnts.idx(1)), cx.this_step);
        end
    elseif locating
        if s.evnts(1).zero      %% found the zero!
            %% step size will be reset to previously used default step size
            locating = 0;           %% exit "locating" mode
            rb_cnt_ef = 0;          %% reset rollback counter for ef intervals
            if opt.verbose > 3
                loc_msg = sprintf('ZERO!      : f = %g, step <-- %.4g', ...
                    nx.ef{rx_evnts.eidx}(rx_evnts.idx(1)), nx.default_step);
            end
        else                    %% prev rollback undershot
            %% initialize next step size based on critical event function
            %% values from prev rollback step and current step
            rx_ef = rx.ef{rx_evnts.eidx}(rx_evnts.idx(1));
            cx_ef = nx.ef{rx_evnts.eidx}(rx_evnts.idx(1));
            step_scale = cx_ef / (cx_ef - rx_ef);
            nx.this_step = step_scale * (rx.step - nx.step);
            rb_cnt_ef = 0;          %% reset rollback counter for ef intervals
            if opt.verbose > 3
                loc_msg = sprintf('UNDERSHOOT : f [<<%g>>, %g], step <-- %.4g', ...
                    cx_ef, rx_ef, nx.this_step);
            end
        end
    else
        loc_msg = '';
        direction = 1;
    end

    %% invoke callbacks - "iterations" context
    rb = s.rollback;
    for k = 1:ncb
        [nx, cx, s] = reg_cb(k).fcn(cont_steps, nx, cx, px, s, opt);
    end
    if ~rb && s.rollback    %% rollback triggered by callback (vs event function interval)
        rb_cnt_cb = rb_cnt_cb + 1;  %% increment rollback counter for callbacks
        if rb_cnt_cb > 26
            s.done = 1;
            s.done_msg = 'Too many rollback steps triggered by callbacks!';
        end
    else
        rb_cnt_cb = 0;              %% reset rollback counter for callbacks
    end

    %% print iteration information
    if opt.verbose > 1
        %% set label for rollback step counter
        if rb_cnt_ef
            sub_step = char('a' + rb_cnt_ef - 1);
        elseif rb_cnt_cb
            sub_step = char('A' + rb_cnt_cb - 1);
        else
            sub_step = ' ';
        end

        fprintf('step %3d%s : %s stepsize = %-9.3g lambda = %6.3f', cont_steps, sub_step, pne_ptag(cx.parm), cx.step, nx.x(end));
        if opt.verbose < 5
            fprintf('  %2d corrector steps', out.iterations);
        end
        if s.rollback
            fprintf(' ^ ROLLBACK\n');
        else
            fprintf('\n');
        end
        if opt.verbose > 3 && ~isempty(loc_msg)
            fprintf('    LOCATING -- %s\n', loc_msg);
        end
    end

    %% log events
    for k = 1:length(s.evnts)
        if s.evnts(k).log
            e = struct( 'k', cont_steps, ...
                        'name', s.evnts(k).name, ...
                        'idx', s.evnts(k).idx, ...
                        'msg', s.evnts(k).msg   );
            if isempty(nx.events)
                nx.events = e;
            else
                nx.events(end+1) = e;
            end
        end
        if (opt.verbose > 2 && s.evnts(k).log) || ...
                (opt.verbose > 3 && s.evnts(k).eidx)
            fprintf('    %s\n', s.evnts(k).msg);
        end
    end

    %% adapt stepsize if requested and not terminating, locating a zero
    %% or re-doing a step after changing the problem data
    if opt.adapt_step && ~s.done && ~locating && ~s.evnts(1).zero && nx.step ~= 0
        pred_error = norm(nx.x - nx.x_hat, Inf);

        %% new nominal step size is current size * tol/err, but we reduce
        %% the change from the current size by a damping factor and limit
        %% increases to a factor of 2
        step_scale = min(2, 1 + opt.adapt_step_damping * ...
                        (opt.adapt_step_tol/pred_error - 1));
        nx.default_step = nx.step * step_scale;

        %% limit step-size
        if nx.default_step > opt.step_max
            nx.default_step = opt.step_max;
        end
        if nx.default_step < opt.step_min
            nx.default_step = opt.step_min;
        end
    end

    %% if this is a normal step
    if ~s.rollback
        px = cx;    %% save current state before update
        cx = nx;    %% update current state to next candidate
%         if cont_steps > 1000
%             s.done = 1;
%             s.done_msg = 'RUNAWAY!';
%         end
        if ~s.done
            cont_steps = cont_steps + 1;
        end
    end

    %% set step size and parameterization, from one-time or defaults
    if isempty(cx.this_step)
        cx.step = cx.default_step;
    else
        cx.step = cx.this_step;
        cx.this_step = [];      %% disable for next time
    end
    if isempty(cx.this_parm)
        cx.parm = cx.default_parm;
    else
        cx.parm = cx.this_parm;
        cx.this_parm = [];      %% disable for next time
    end
end     %% while ~s.done

%% prepare to exit
if isempty(s.warmstart)
    %% invoke callbacks - "final" context
    s.results = struct();   %% initialize results struct
    for k = 1:ncb
        [nx, cx, s] = reg_cb(k).fcn(-cont_steps, nx, cx, px, s, opt);
    end
    output = s.results;
    output.done_msg = s.done_msg;
    output.events = cx.events;  %% copy eventlog to results
    output.corrector = out;     %% output from last corrector run

    if opt.verbose
        fprintf('CONTINUATION TERMINATION: %s\n', s.done_msg);
    end
else
    %% save warmstart values
    ws = s.warmstart;
    ws.cont_steps = cont_steps;
    ws.direction = direction;

    %% from state at current step
    ws.x = cx.x;            %% state from current step
    ws.z = cx.z;            %% tangent vector from current step
    ws.parm = cx.parm;
    ws.default_parm = cx.default_parm;
    ws.default_step = cx.default_step;
    ws.events = cx.events;
    ws.cbx = cx.cb;

    %% from state at previous step
    ws.xp = px.x;           %% state from previous step
    ws.zp = px.z;           %% tangent vector from previous step

    output.warmstart = ws;
    output.iterations = max(ws.cbx.default.iterations);
    output.max_lam = max(ws.cbx.default.lam);
    output.events = ws.events;
    if opt.verbose
        fprintf('%s : CONTINUATION SUSPENDED ...\n', s.done_msg);
    end
end

%% output arguments
if nargout > 4
    [f, J] = fcn(cx.x);
elseif nargout > 1
    f = fcn(cx.x);
end
varargout{1} = cx.x;
if nargout > 1
    varargout{2} = f;
    if nargout > 2
        varargout{3} = exitflag;
        if nargout > 3
            varargout{4} = output;
            if nargout > 4
                varargout{5} = J;
            end
        end
    end
end

function [fp, dfp] = pne_corrector_fcn(x, fcn, pfcn)
if nargout < 2
    fp = [ fcn(x); pfcn(x) ];
else
    [f, df] = fcn(x);
    [p, dp] = pfcn(x);
    fp = [f; p];
    dfp = [df; dp];
end

function z = pne_tangent(x, xp, zp, fcn, parm, direction)
pfcn = @(xx)parm(xx, xp, 0, zp);
[f, df] = fcn(x);
[p, dp] = pfcn(x);
rhs = [ zeros(length(f), 1); direction ];
z = [df; dp] \ rhs;
z = z / norm(z);    %% normalize it

function ptag = pne_ptag(parm)
switch func2str(parm)
    case 'pne_pfcn_natural'
        ptag = 'NAT';
    case 'pne_pfcn_arc_len'
        ptag = 'ARC';
    case 'pne_pfcn_pseudo_arc_len'
        ptag = 'PAL';
end
