function t_om_solve_pne(quiet)
%T_OM_SOLVE_PNE  Tests of PNE solvers via OM.SOLVE().

%   MP-Opt-Model
%   Copyright (c) 2010-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin < 1
    quiet = 0;
end

%%  alg         name        check       opts
cfg = {
    {'DEFAULT', 'default',  []          []  },
};

n = 154;

t_begin(n*length(cfg), quiet);

for k = 1:length(cfg)
    alg   = cfg{k}{1};
    name  = cfg{k}{2};
    check = cfg{k}{3};
    opts  = cfg{k}{4};
    if ~isempty(check) && ~have_feature(check)
        t_skip(n, sprintf('%s not installed', name));
    else
        %% default start point
        y0 = [-1;0];
        opt = struct( ...
            'verbose', 0, ...
            'alg', alg, ...
            'stop_at', 0.7, ...
            'parameterization', 1, ...
            'nleqs_opt', struct('tol', 1e-9), ...
            'adapt_step', 0, ...
            'step_max', 10, ...
            'target_lam_tol', 1e-6, ...
            'nose_tol', 1e-6, ...
            'adapt_step_tol', 1e-2 );
%         opt.plot = struct('level', 2, 'idx', 2);
        om = opt_model;
        om.add_var('y', 2, y0);
        om.add_var('lam', 1, 0);
        om.add_nln_constraint('f', 2, 1, @f1p, [], {'y', 'lam'});

        t = sprintf('%s - TARGET_LAM = 0.7 (natural) : ', name);
        [x, f, e, out, jac] = om.solve(opt);
        it = 14;
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [-1.931782106; -1.268217894; 0.7], 6, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_is(out.cont.x(:,1), [-3;4;0], 8, [t 'out.cont.x(:,1)']);
        t_is(out.cont.max_lam, 0.7, 10, [t 'out.cont.max_lam']);
        t_is(out.cont.iterations, it, 12, [t 'out.cont.iterations']);
        t_is(size(out.cont.lam), [1,it+1], 12, [t 'size(out.cont.lam)']);
        t_is(size(out.cont.lam_hat), [1,it+1], 12, [t 'size(out.cont.lam_hat)']);
        t_is(size(out.cont.x), [3,it+1], 12, [t 'size(out.cont.x)']);
        t_is(size(out.cont.x_hat), [3,it+1], 12, [t 'size(out.cont.x_hat)']);
        t_is(size(out.cont.steps), [1,it+1], 12, [t 'size(out.cont.steps)']);
        t_is(length(out.cont.events), 1, 12, [t 'length(out.cont.events)']);
        t_is(out.cont.events.k, it, 12, [t 'out.cont.events.k']);
        t_is(out.cont.events.idx, 1, 12, [t 'out.cont.events.idx']);
        t_ok(strcmp(out.cont.events.name, 'TARGET_LAM'), [t 'out.cont.events.name']);
        t_ok(strcmp(out.cont.done_msg, sprintf('Reached desired lambda 0.7 in %d continuation steps', it)), [t 'out.cont.done_msg']);

        t = sprintf('%s - TARGET_LAM = 0.7 (arc len) : ', name);
        opt.adapt_step = 1;
        opt.parameterization = 2;
        [x, f, e, out, jac] = om.solve(opt);
        it = 9;
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [-1.931782106; -1.268217894; 0.7], 6, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_is(out.cont.x(:,1), [-3;4;0], 8, [t 'out.cont.x(:,1)']);
        t_is(out.cont.max_lam, 0.7, 10, [t 'out.cont.max_lam']);
        t_is(out.cont.iterations, it, 12, [t 'out.cont.iterations']);
        t_is(size(out.cont.lam), [1,it+1], 12, [t 'size(out.cont.lam)']);
        t_is(size(out.cont.lam_hat), [1,it+1], 12, [t 'size(out.cont.lam_hat)']);
        t_is(size(out.cont.x), [3,it+1], 12, [t 'size(out.cont.x)']);
        t_is(size(out.cont.x_hat), [3,it+1], 12, [t 'size(out.cont.x_hat)']);
        t_is(size(out.cont.steps), [1,it+1], 12, [t 'size(out.cont.steps)']);
        t_is(length(out.cont.events), 1, 12, [t 'length(out.cont.events)']);
        t_is(out.cont.events.k, it, 12, [t 'out.cont.events.k']);
        t_is(out.cont.events.idx, 1, 12, [t 'out.cont.events.idx']);
        t_ok(strcmp(out.cont.events.name, 'TARGET_LAM'), [t 'out.cont.events.name']);
        t_ok(strcmp(out.cont.done_msg, sprintf('Reached desired lambda 0.7 in %d continuation steps', it)), [t 'out.cont.done_msg']);

        t = sprintf('%s - TARGET_LAM = 0.7 (pseudo arc len) : ', name);
        opt.parameterization = 3;
        [x, f, e, out, jac] = om.solve(opt);
        it = 9;
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [-1.931782106; -1.268217894; 0.7], 6, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_is(out.cont.x(:,1), [-3;4;0], 8, [t 'out.cont.x(:,1)']);
        t_is(out.cont.max_lam, 0.7, 10, [t 'out.cont.max_lam']);
        t_is(out.cont.iterations, it, 12, [t 'out.cont.iterations']);
        t_is(size(out.cont.lam), [1,it+1], 12, [t 'size(out.cont.lam)']);
        t_is(size(out.cont.lam_hat), [1,it+1], 12, [t 'size(out.cont.lam_hat)']);
        t_is(size(out.cont.x), [3,it+1], 12, [t 'size(out.cont.x)']);
        t_is(size(out.cont.x_hat), [3,it+1], 12, [t 'size(out.cont.x_hat)']);
        t_is(size(out.cont.steps), [1,it+1], 12, [t 'size(out.cont.steps)']);
        t_is(length(out.cont.events), 1, 12, [t 'length(out.cont.events)']);
        t_is(out.cont.events.k, it, 12, [t 'out.cont.events.k']);
        t_is(out.cont.events.idx, 1, 12, [t 'out.cont.events.idx']);
        t_ok(strcmp(out.cont.events.name, 'TARGET_LAM'), [t 'out.cont.events.name']);
        t_ok(strcmp(out.cont.done_msg, sprintf('Reached desired lambda 0.7 in %d continuation steps', it)), [t 'out.cont.done_msg']);

        t = sprintf('%s - FULL : ', name);
        opt.stop_at = 'FULL';
        [x, f, e, out, jac] = om.solve(opt);
        it = 34;
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [2;-1;0], 10, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_is(out.cont.x(:,1), [-3;4;0], 8, [t 'out.cont.x(:,1)']);
        t_is(out.cont.max_lam, 1.04127275, 8, [t 'out.cont.max_lam']);
        t_is(out.cont.iterations, it, 12, [t 'out.cont.iterations']);
        t_is(size(out.cont.lam), [1,it+1], 12, [t 'size(out.cont.lam)']);
        t_is(size(out.cont.lam_hat), [1,it+1], 12, [t 'size(out.cont.lam_hat)']);
        t_is(size(out.cont.x), [3,it+1], 12, [t 'size(out.cont.x)']);
        t_is(size(out.cont.x_hat), [3,it+1], 12, [t 'size(out.cont.x_hat)']);
        t_is(size(out.cont.steps), [1,it+1], 12, [t 'size(out.cont.steps)']);
        t_is(length(out.cont.events), 1, 12, [t 'length(out.cont.events)']);
        t_is(out.cont.events.k, it, 12, [t 'out.cont.events.k']);
        t_is(out.cont.events.idx, 1, 12, [t 'out.cont.events.idx']);
        t_ok(strcmp(out.cont.events.name, 'TARGET_LAM'), [t 'out.cont.events.name']);
        t_ok(strcmp(out.cont.done_msg, sprintf('Traced full continuation curve in %d continuation steps', it)), [t 'out.cont.done_msg']);

        t = sprintf('%s - NOSE (arc len): ', name);
        opt.stop_at = 'NOSE';
        opt.parameterization = 2;
        [x, f, e, out, jac] = om.solve(opt);
        it = 18;
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [-0.5; -4.75; 1.04166667], 6, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_is(out.cont.x(:,1), [-3;4;0], 8, [t 'out.cont.x(:,1)']);
        t_is(out.cont.max_lam, 1.04166666667, 10, [t 'out.cont.max_lam']);
        t_is(out.cont.iterations, it, 12, [t 'out.cont.iterations']);
        t_is(size(out.cont.lam), [1,it+1], 12, [t 'size(out.cont.lam)']);
        t_is(size(out.cont.lam_hat), [1,it+1], 12, [t 'size(out.cont.lam_hat)']);
        t_is(size(out.cont.x), [3,it+1], 12, [t 'size(out.cont.x)']);
        t_is(size(out.cont.x_hat), [3,it+1], 12, [t 'size(out.cont.x_hat)']);
        t_is(size(out.cont.steps), [1,it+1], 12, [t 'size(out.cont.steps)']);
        t_is(length(out.cont.events), 1, 12, [t 'length(out.cont.events)']);
        t_is(out.cont.events.k, it, 12, [t 'out.cont.events.k']);
        t_is(out.cont.events.idx, 1, 12, [t 'out.cont.events.idx']);
        t_ok(strcmp(out.cont.events.name, 'NOSE'), [t 'out.cont.events.name']);
        t_ok(strcmp(out.cont.done_msg, sprintf('Reached limit in %d continuation steps, lambda = 1.042.', it)), [t 'out.cont.done_msg']);

        t = sprintf('%s - NOSE (pseudo arc len): ', name);
        opt.parameterization = 3;
        [x, f, e, out, jac] = om.solve(opt);
        it = 18;
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [-0.5; -4.75; 1.04166667], 6, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_is(out.cont.x(:,1), [-3;4;0], 8, [t 'out.cont.x(:,1)']);
        t_is(out.cont.max_lam, 1.04166666667, 10, [t 'out.cont.max_lam']);
        t_is(out.cont.iterations, it, 12, [t 'out.cont.iterations']);
        t_is(size(out.cont.lam), [1,it+1], 12, [t 'size(out.cont.lam)']);
        t_is(size(out.cont.lam_hat), [1,it+1], 12, [t 'size(out.cont.lam_hat)']);
        t_is(size(out.cont.x), [3,it+1], 12, [t 'size(out.cont.x)']);
        t_is(size(out.cont.x_hat), [3,it+1], 12, [t 'size(out.cont.x_hat)']);
        t_is(size(out.cont.steps), [1,it+1], 12, [t 'size(out.cont.steps)']);
        t_is(length(out.cont.events), 1, 12, [t 'length(out.cont.events)']);
        t_is(out.cont.events.k, it, 12, [t 'out.cont.events.k']);
        t_is(out.cont.events.idx, 1, 12, [t 'out.cont.events.idx']);
        t_ok(strcmp(out.cont.events.name, 'NOSE'), [t 'out.cont.events.name']);
        t_ok(strcmp(out.cont.done_msg, sprintf('Reached limit in %d continuation steps, lambda = 1.042.', it)), [t 'out.cont.done_msg']);

        t = sprintf('%s - NOSE (opp dir) : ', name);
        y0 = [1;-1];
        om.set_params('var', 'y', 'v0', y0);
        [x, f, e, out, jac] = om.solve(opt);
        it = 20;
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [-0.5; -4.75; 1.04166667], 5, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_is(out.cont.x(:,1), [2;-1;0], 8, [t 'out.cont.x(:,1)']);
        t_is(out.cont.max_lam, 1.04166666667, 10, [t 'out.cont.max_lam']);
        t_is(out.cont.iterations, it, 12, [t 'out.cont.iterations']);
        t_is(size(out.cont.lam), [1,it+1], 12, [t 'size(out.cont.lam)']);
        t_is(size(out.cont.lam_hat), [1,it+1], 12, [t 'size(out.cont.lam_hat)']);
        t_is(size(out.cont.x), [3,it+1], 12, [t 'size(out.cont.x)']);
        t_is(size(out.cont.x_hat), [3,it+1], 12, [t 'size(out.cont.x_hat)']);
        t_is(size(out.cont.steps), [1,it+1], 12, [t 'size(out.cont.steps)']);
        t_is(length(out.cont.events), 1, 12, [t 'length(out.cont.events)']);
        t_is(out.cont.events.k, it, 12, [t 'out.cont.events.k']);
        t_is(out.cont.events.idx, 1, 12, [t 'out.cont.events.idx']);
        t_ok(strcmp(out.cont.events.name, 'NOSE'), [t 'out.cont.events.name']);
        t_ok(strcmp(out.cont.done_msg, sprintf('Reached limit in %d continuation steps, lambda = 1.042.', it)), [t 'out.cont.done_msg']);

        t = sprintf('%s - FULL warmstart, before SWITCH : ', name);
        opt.parameterization = 3;
        opt.adapt_step = 1;
        opt.adapt_step_ws = 0.25;
        opt.stop_at = 'FULL';
        opt.adapt_step_tol = 5e-3;
        opt.callbacks = {@pne_callback_test1};
        opt.events = {{'SWITCH!', @pne_event_test1, 1e-6}};
        opt.output_fcn = @pne_output_fcn_test1;
        % opt.verbose = 4;
        % opt.plot = struct('level', 2, 'yname', 'y', 'idx', [1;2]);
        y0 = [-1;0];
        om.set_params('var', 'y', 'v0', y0);
        [x, f, e, out, jac] = om.solve(opt);
        it = 10;
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [-2; -1; 2/3], 6, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_ok(isfield(out, 'warmstart') && ~isempty(out.warmstart), [t 'out.warmstart exists']);
        t_is(out.warmstart.cont_steps, it, 12, [t 'out.warmstart.cont_steps']);
        t_is(out.warmstart.default_step, 0.64432407, 8, [t 'out.warmstart.default_step']);
        t_ok(isa(out.warmstart.parm, 'function_handle'), [t 'out.warmstart.parm is function']);
        t_ok(isa(out.warmstart.default_parm, 'function_handle'), [t 'out.warmstart.default_parm is function']);
        t_ok(isstruct(out.warmstart.cbx), [t 'out.warmstart.cbx is struct']);
        t_ok(isstruct(out.warmstart.cbx.default), [t 'out.warmstart.cbx.default is struct']);
        t_is(out.warmstart.cbx.default.iterations, it, 12, [t 'out.warmstart.cbx.default.iterations']);
        t_ok(isstruct(out.warmstart.events), [t 'out.warmstart.events is struct']);
        t_is(out.warmstart.events.k, 10, 12, [t 'out.warmstart.events.k']);
        t_is(out.warmstart.events.idx, 1, 12, [t 'out.warmstart.events.idx']);
        t_ok(strcmp(out.warmstart.events.name, 'SWITCH!'), [t 'out.warmstart.events.name']);
        t_ok(strcmp(out.warmstart.events.msg, 'ZERO detected for SWITCH! event'), [t 'out.warmstart.events.msg']);
        t_ok(~isfield(out, 'cont'), [t 'out.cont does not exist']);

        t = sprintf('%s - FULL warmstart, after SWITCH : ', name);
        y0 = [-1;0];
        om = opt_model;
        om.add_var('z', 1, 0);
        om.add_var('y', 2, y0);
        om.add_var('lam', 1, 0);
        om.add_nln_constraint('f', 2, 1, @f2p, [], {'y', 'lam'});
        om.add_nln_constraint('g', 1, 1, @g2p, [], {'z', 'lam'});
        opt.x0 = [x(end)/2; x];
        opt.warmstart = out.warmstart;
        opt.events = {{'SNOUT!', @pne_event_test2, 1e-6}};
        opt.callbacks = {};
        [x, f, e, out, jac] = om.solve(opt);
        it = 49;
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [0; 6; -5; 0], 6, [t 'x - final']);
        t_is(f, [0;0;0], 10, [t 'f']);
        t_is(out.cont.y(:,1), [-3;4], 8, [t 'out.cont.y(:,1)']);
        t_is(out.cont.max_lam, 1.04166666667, 10, [t 'out.cont.max_lam']);
        t_is(out.cont.iterations, it, 12, [t 'out.cont.iterations']);
        t_is(size(out.cont.lam), [1,it+1], 12, [t 'size(out.cont.lam)']);
        t_is(size(out.cont.lam_hat), [1,it+1], 12, [t 'size(out.cont.lam_hat)']);
        t_is(size(out.cont.y), [2,it+1], 12, [t 'size(out.cont.y)']);
        t_is(size(out.cont.y_hat), [2,it+1], 12, [t 'size(out.cont.y_hat)']);
        t_is(size(out.cont.steps), [1,it+1], 12, [t 'size(out.cont.steps)']);
        t_is(length(out.cont.events), 3, 12, [t 'length(out.cont.events)']);
        t_is(out.cont.events(1).k, 10, 12, [t 'out.cont.events(1).k']);
        t_is(out.cont.events(1).idx, 1, 12, [t 'out.cont.events(1).idx']);
        t_ok(strcmp(out.cont.events(1).name, 'SWITCH!'), [t 'out.cont.events(1).name']);
        t_ok(strcmp(out.cont.events(1).msg, 'ZERO detected for SWITCH! event'), [t 'out.cont.events(1).msg']);
        t_is(out.cont.events(2).k, 28, 12, [t 'out.cont.events(2).k']);
        t_is(out.cont.events(2).idx, 1, 12, [t 'out.cont.events(2).idx']);
        t_ok(strcmp(out.cont.events(2).name, 'SNOUT!'), [t 'out.cont.events(2).name']);
        t_ok(strcmp(out.cont.events(2).msg, 'ZERO detected for SNOUT! event'), [t 'out.cont.events(2).msg']);
        t_is(out.cont.events(3).k, it, 12, [t 'out.cont.events(3).k']);
        t_is(out.cont.events(3).idx, 1, 12, [t 'out.cont.events(3).idx']);
        t_ok(strcmp(out.cont.events(3).name, 'TARGET_LAM'), [t 'out.cont.events(3).name']);
        t_ok(strcmp(out.cont.events(3).msg, 'ZERO detected for TARGET_LAM event'), [t 'out.cont.events(3).msg']);
        t_ok(strcmp(out.cont.done_msg, sprintf('Traced full continuation curve in %d continuation steps', it)), [t 'out.cont.done_msg']);
    end
end

t_end;


%% 2-d problem with 2 solutions
%% from https://www.chilimath.com/lessons/advanced-algebra/systems-non-linear-equations/
function [f, J] = f1(x)
f = [  x(1)   + x(2) - 1;
      -x(1)^2 + x(2) + 5    ];
if nargout > 1
    J = sparse([1 1; -2*x(1) 1]);
end

%% parameterized 2-d problem
%% based on https://www.chilimath.com/lessons/advanced-algebra/systems-non-linear-equations/
function [f, J] = f1p(xx)
[y, lam] = deal(xx{:});
if nargout < 2
    f = f1(y);
else
    [f, J] = f1(y);
    J = [J [6;0]];
end
f = f + [6*lam; 0];

%% another parameterized 2-d problem, based on f1p
function [f, J] = f2p(xx)
[y, lam] = deal(xx{:});
y = [y(2)+2; y(1)-2];
if nargout < 2
    f = f1p({y,lam});
else
    [f, JJ] = f1p({y,lam});
    J = [JJ(:, 2) JJ(:, 1) JJ(:, 3)];
end

%% another parameterized 2-d equation
function [g, J] = g2p(xx)
[z, lam] = deal(xx{:});
g = 2*z-lam;
if nargout >= 2
    J = [2 -1];
end

%% example custom event function 1 (target lambda == 2/3)
function ef = pne_event_test1(cx, opt)
tlam = 2/3;
ef = cx.x(end) - tlam;

%% example custom event function 2 (nose point)
function ef = pne_event_test2(cx, opt)
ef = cx.z(end);

%% example custom callback function (exit for warmstart on SWITCH! event)
function [nx, cx, s] = pne_callback_test1(k, nx, cx, px, s, opt)
if k <= 0 || s.done, return; end    %% skip if initialize, finalize or done
tlam = 2/3;
for i = 1:length(s.evnts)
    if strcmp(s.evnts(i).name, 'SWITCH!')
        if s.evnts(i).zero  %% prepare to terminate
            s.done = 1;
            s.done_msg = sprintf('Reached switching point in %d continuation steps', k);
            s.warmstart = struct(); %% signal that we want to exit, then resume
        else                    %% set step-size & parameterization to terminate next time
            cx.this_parm = @pne_pfcn_natural;   %% change to natural parameterization
            cx.this_step = tlam - cx.x(end);
            s.evnts(i).msg = sprintf('%s\n  step %d to overshoot lambda = %g, reduce step size and set natural param', s.evnts(i).msg, k, tlam);
        end
        break;
    end
end

function rv = pne_output_fcn_test1(cbx, x, x_hat)
%% cbx     = pne_output_fcn_test1(cbx, x, x_hat)
%% results = pne_output_fcn_test1(cbx, results)
if nargin == 3      %% store values in callback state
    if length(x) == 3
        k = [1;2];
    else %% == 3
        k = [2;3];
    end
    rv = cbx;
    if isfield(cbx, 'y')    %% append values (ITERATION)
        rv.lam_hat = [ rv.lam_hat   x_hat(end)  ];
        rv.lam     = [ rv.lam       x(end)      ];
        rv.y_hat   = [ rv.y_hat     x_hat(k)    ];
        rv.y       = [ rv.y         x(k)        ];
    else                    %% initialize values (INITIAL)
        rv.lam_hat = x_hat(end);
        rv.lam     = x(end);
        rv.y_hat   = x_hat(k);
        rv.y       = x(k);
    end
else                        %% copy fields to results (FINAL)
    rv = x;
    rv.lam_hat = cbx.lam_hat;
    rv.lam     = cbx.lam;
    rv.max_lam = max(cbx.lam);
    rv.y_hat   = cbx.y_hat;
    rv.y       = cbx.y;
end


% %% another 2-d problem
% %% from Christi Patton Luks, https://www.youtube.com/watch?v=pJG4yhtgerg
% function [f, J] = f2(x)
% m = 10;
% m1 = 0; %2;
% m2 = 0; %5;
% f = [  x(1)^2 +   x(1)*x(2)   - 10 + m1;
%        (x(2)   + 3*x(1)*x(2)^2 - 57) / m  ] + m2;
% if nargout > 1
%     J = sparse([    2*x(1)+x(2)    x(1);
%                     (3*x(2)^2) / m (6*x(1)*x(2)+1) / m    ]);
% end
