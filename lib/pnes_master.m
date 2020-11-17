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
dopts = struct( ...
    'alg',              'DEFAULT', ...  %% algorithm
    'verbose',          0, ...
    'nleqs_opt',        struct('verbose', 0), ...
    'solve_base',       1, ...          %% UNFINISHED
    'parameterization', 3, ...          %% 1 - natural, 2 - arc len, 3 - pseudo arc len
    'stop_at',          'NOSE', ...     %% 'NOSE', 'FULL', <lam_stop>
    'step',             0.05, ...       %% UNFINISHED
    'step_min',         1e-4, ...       %% UNFINISHED
    'step_max',         0.2, ...        %% UNFINISHED
    'adapt_step',       0, ...          %% UNFINISHED
    'adapt_step_damping', 0.7, ...      %% UNFINISHED
    'adapt_step_tol',   1e-3, ...       %% UNFINISHED
    'default_event_tol',1e-3, ...       %% UNFINISHED
    'target_lam_tol',   0, ...          %% UNFINISHED
    'nose_tol',         0, ...          %% UNFINISHED
    'events',           {{}}, ...       %% UNFINISHED
    'callbacks',        {{}}, ...       %% UNFINISHED
    'plot',             struct( ...     %% used by pne_default_callback() for plotting
        'level',        0, ...          %% UNFINISHED
        'idx',          [], ...         %% UNFINISHED
        'xfcn',         @pne_plot_xfcn, ...         %% UNFINISHED
        'yfcn',         @pne_plot_yfcn, ...         %% UNFINISHED
        'title',        'Value of Variable %d', ... %% UNFINISHED
        'title2',       'Value of Multiple Variables', ...  %% UNFINISHED
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
done = struct('flag', 0, 'msg', '');
z0 = zeros(length(x0), 1);      %% (n+1) x 1, zeros

%% register event and callback functions
if ischar(opt.stop_at) && strcmp(opt.stop_at, 'NOSE');
    my_events = {{'NOSE', @pne_nose_event, opt.nose_tol}, opt.events{:}};
    my_cbacks = {{@pne_nose_event_cb, 51}, opt.callbacks{:}};
else        %% FULL or target lambda
    my_events = {{'TARGET_LAM', @pne_target_lam_event, opt.target_lam_tol}, opt.events{:}};
    my_cbacks = {{@pne_target_lam_event_cb, 50}, opt.callbacks{:}};
end
my_cbacks{end+1} = {@pne_default_callback, 0};
reg_ev = pne_register_events(my_events, opt);   %% registered event functions
reg_cb = pne_register_callbacks(my_cbacks);     %% registered callback functions
nef = length(reg_ev);   %% number of registered event functions
ncb = length(reg_cb);   %% number of registered callback functions

t0 = tic;                       %% start timing

if opt.verbose
    v = mpomver('all');
    fprintf('\nMP-Opt-Model Version %s, %s', v.Version, v.Date);
    fprintf(' -- Predictor/Corrector Continuation Method\n');
end

%% solve corrector step for base point
if opt.solve_base
    pfcn = @(xx)pne_pfcn_natural(xx, x0, 0);
    [x, f, exitflag, out] = nleqs_master(@(xx)pne_corrector_fcn(xx, fcn, pfcn), x0, opt.nleqs_opt);
    if exitflag
        if opt.verbose > 1
            fprintf('step %3d  :                      lambda = %6.3f, %2d corrector steps\n', 0, 0, out.iterations);
        end
    else
        done.flag = 1;
        done.msg = sprintf('base solution did not converge in %d iterations', out.iterations);
        if opt.verbose
            fprintf('%s\n', done.msg);
        end
    end
else
    x = x0;
end

%% initialize numerical continuation
if ~done.flag
    step = opt.step;
    cont_steps = 0; %% continuation step counter
    rollback = 0;   %% flag to indicate that a step must be rolled back
    locating = 0;   %% flag to indicate that an event interval was detected,
                    %% but the event has not yet been located
    rb_cnt_ef = 0;  %% counter for rollback steps triggered by event function intervals
    rb_cnt_cb = 0;  %% counter for rollback steps triggered directly by callbacks

    %% initialize parameterization function
    switch opt.parameterization
        case 1
            parm = @pne_pfcn_natural;
        case 2
            parm = @pne_pfcn_arc_len;
        case 3
            parm = @pne_pfcn_pseudo_arc_len;
        otherwise
            error('pnes_master: OPT.parameterization (= %d) must be 1, 2, or 3', parm);
    end

    %% initialize tangent: z = dx
    direction = 1;
    z = zeros(length(x0), 1); z(end) = direction;   %% direction of positive lambda
    z = pne_tangent(x, x, z, fcn, parm, direction);

    %% initialize state for current continuation step
    cx = struct( ...        %% current state
        'x_hat',        x, ...      %% predicted solution value
        'x',            x, ...      %% corrected solution value
        'z',            z, ...      %% normalized tangent vector
        'default_step', step, ...   %% default step size
        'default_parm', parm, ...   %% default parameterization
        'this_step', [], ...        %% step size for this step only
        'this_parm', [], ...        %% parameterization for this step only
        'step', step, ...           %% current step size
        'parm', parm, ...           %% current parameterization
        'events', [], ...           %% event log
        'cb', struct(), ...         %% user state, for callbacks
        'ef', {cell(nef, 1)} ...    %% event function values
    );

    %% initialize event function values
    for k = 1:nef
        cx.ef{k} = reg_ev(k).fcn(cx, opt);
    end

    %% invoke callbacks - "initialize" context
    for k = 1:ncb
        [nx, cx, done, rollback, evnts, opt] = reg_cb(k).fcn( ...
            cont_steps, cx, cx, cx, done, 0, [], opt);
    end
    
    %% check for case with base and target the same
    if 0
        done.flag = 1;
        done.msg = 'no difference between base and target';
    end
    
    cont_steps = cont_steps + 1;
    px = cx;    %% initialize state for previous continuation step
end

%%-----  run numerical continuation  -----
while ~done.flag
    %% initialize next candidate with current state
    nx = cx;

    %% predictor step
    nx.x_hat = cx.x + cx.step * cx.z;
%x_hat = nx.x_hat

    %% corrector step
    pfcn = @(xx)cx.parm(xx, cx.x, cx.step, cx.z);
    [nx.x, f, exitflag, out] = nleqs_master(@(xx)pne_corrector_fcn(xx, fcn, pfcn), nx.x_hat, opt.nleqs_opt);
    if ~exitflag    %% corrector failed
        done.flag = 1;
        done.msg = sprintf('Corrector did not converge in %d iterations.', out.iterations);
        if opt.verbose
            fprintf('step %3d  : stepsize = %-9.3g lambda = %6.3f  corrector did not converge in %d iterations\n', cont_steps, cx.step, nx.x(end), out.iterations);
        end
        cont_steps = max(cont_steps - 1, 1);    %% go back to last step, but not to 0
        break;
    end

%x = nx.x

    %% compute new tangent direction, based on a previous state: tx
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
    [rollback, evnts, nx.ef] = pne_detect_events(reg_ev, nx.ef, cx.ef, nx.step);

    %% adjust step-size to locate event function zero, if necessary
    if rollback                 %% current step overshot
        %% rollback and initialize next step size based on rollback and previous
        rx = nx;                    %% save state we're rolling back from
        rx.evnts = evnts;           %% and critical event info
        cx.this_step = evnts.step_scale * rx.step;
        cx.this_parm = rx.parm;     %% keep same parameterization as last step
        locating = 1;               %% enter "locating" mode (or stay in it)
        rb_cnt_ef = rb_cnt_ef + 1;  %% increment rollback counter for ef intervals
        if rb_cnt_ef > 26
            done.flag = 1;
            done.msg = sprintf('Could not locate %s event!', evnts.name);
        end
        if opt.verbose > 3
            loc_msg = sprintf('OVERSHOOT  : f = [%g, <<%g>>], step <-- %.4g', ...
                        cx.ef{evnts.eidx}(evnts.idx(1)), ...
                        rx.ef{evnts.eidx}(evnts.idx(1)), cx.this_step);
        end
    elseif locating
        if evnts(1).zero        %% found the zero!
            %% step size will be reset to previously used default step size
            locating = 0;           %% exit "locating" mode
            rb_cnt_ef = 0;          %% reset rollback counter for ef intervals
            if opt.verbose > 3
                loc_msg = sprintf('ZERO!      : f = %g, step <-- %.4g', ...
                    nx.ef{rx.evnts.eidx}(rx.evnts.idx(1)), nx.default_step);
            end
        else                    %% prev rollback undershot
            %% initialize next step size based on critical event function
            %% values from prev rollback step and current step
            rx_ef = rx.ef{rx.evnts.eidx}(rx.evnts.idx(1));
            cx_ef = nx.ef{rx.evnts.eidx}(rx.evnts.idx(1));
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
    rb = rollback;
    for k = 1:ncb
        [nx, cx, done, rollback, evnts, opt] = reg_cb(k).fcn( ...
            cont_steps, nx, cx, px, done, rollback, evnts, opt);
    end
    if ~rb && rollback      %% rollback triggered by callback (vs event function interval)
        rb_cnt_cb = rb_cnt_cb + 1;  %% increment rollback counter for callbacks
        if rb_cnt_cb > 26
            done.flag = 1;
            done.msg = 'Too many rollback steps triggered by callbacks!';
        end
    else
        if ~done.flag && evnts(1).zero
            %% decide whether to switch directions
            reply = input('Switch directions? Y/N [N]:','s');
            if strcmp(upper(reply), 'Y')
                direction = -direction;
            end
        end
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
        if opt.verbose > 4
            fprintf('step %3d%s : stepsize = %-9.3g lambda = %6.3f', cont_steps, sub_step, cx.step, nx.x(end));
        else
            fprintf('step %3d%s : stepsize = %-9.3g lambda = %6.3f  %2d corrector steps', cont_steps, sub_step, cx.step, nx.x(end), out.iterations);
        end
        if rollback
            fprintf(' ^ ROLLBACK\n');
        else
            fprintf('\n');
        end
        if opt.verbose > 3 && ~isempty(loc_msg)
            fprintf('    LOCATING -- %s\n', loc_msg);
        end
    end

    %% log events
    for k = 1:length(evnts)
        if evnts(k).log
            e = struct( 'k', cont_steps, ...
                        'name', evnts(k).name, ...
                        'idx', evnts(k).idx, ...
                        'msg', evnts(k).msg   );
            if isempty(nx.events)
                nx.events = e;
            else
                nx.events(end+1) = e;
            end
        end
        if (opt.verbose > 2 && evnts(k).log) || ...
                (opt.verbose > 3 && evnts(k).eidx)
            fprintf('    %s\n', evnts(k).msg);
        end
    end

    %% adapt stepsize if requested and not terminating, locating a zero
    %% or re-doing a step after changing the problem data
    if opt.adapt_step && ~done.flag && ~locating && ~evnts(1).zero && nx.step ~= 0
        pred_error = norm(nx.x - nx.x_hat, inf);

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
    if ~rollback
        px = cx;    %% save current state before update
        cx = nx;    %% update current state to next candidate
%         if cont_steps > 1000
%             done.flag = 1;
%             done.msg = 'RUNAWAY!';
%         end
        if ~done.flag
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
end     %% while ~done.flag


%% invoke callbacks - "final" context
results = struct();     %% initialize results struct
for k = 1:ncb
    [nx, cx, done, rollback, evnts, opt, results] = ...
        reg_cb(k).fcn(-cont_steps, nx, cx, px, ...
            done, rollback, evnts, opt, results);
end
results.events = cx.events;     %% copy eventlog to results
results.done_msg = done.msg;
out.cont = results;
if opt.verbose
    fprintf('CONTINATION TERMINATION: %s\n', done.msg);
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
            varargout{4} = out;
            if nargout > 4
                varargout{5} = J;
            end
        end
    end
end

% function [f, J] = non_param_fcn(x, fcn)
% if nargout < 2
%     f = fcn([x; 0]);
% else
%     [f, J] = fcn([x; 0]);
%     J(:, end) = [];     %% delete last col
% end

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

function xx = pne_plot_xfcn(x)
xx = x(end, :);

function yy = pne_plot_yfcn(x, idx)
yy = x(idx, :);
