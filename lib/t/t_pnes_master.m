function t_pnes_master(quiet)
%T_PNES_MASTER  Tests of PNE solvers via PNES_MASTER().

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

n = 44;

t_begin(n*length(cfg), quiet);

for k = 1:length(cfg)
    alg   = cfg{k}{1};
    name  = cfg{k}{2};
    check = cfg{k}{3};
    opts  = cfg{k}{4};
    if ~isempty(check) && ~have_feature(check)
        t_skip(n, sprintf('%s not installed', name));
    else
        t = sprintf('%s - FULL : ', name);
        x0 = [-1;0;0];
        opt = struct( ...
            'verbose', 0, ...
            'alg', alg, ...
            'stop_at', 'FULL', ...
            'nleqs_opt', struct('tol', 1e-8), ...
            'adapt_step', 1, ...
            'step_max', 10, ...
            'target_lam_tol', 1e-6, ...
            'nose_tol', 1e-6, ...
            'adapt_step_tol', 1e-2 );
%         opt.plot = struct('level', 2, 'idx', 2);
        [x, f, e, out, jac] = pnes_master(@f1p, x0, opt);
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [2;-1;0], 10, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_is(out.cont.x(:,1), [-3;4;0], 8, [t 'out.cont.x(:,1) : x - initial']);
        t_is(out.cont.max_lam, 1.04127275, 8, [t 'out.cont.max_lam']);
        t_is(out.cont.iterations, 34, 12, [t 'out.cont.iterations']);
        t_is(length(out.cont.events), 1, 12, [t 'length(out.cont.events)']);
        t_is(out.cont.events.k, 34, 12, [t 'out.cont.events.k']);
        t_is(out.cont.events.idx, 1, 12, [t 'out.cont.events.idx']);
        t_ok(strcmp(out.cont.events.name, 'TARGET_LAM'), [t 'out.cont.events.name']);
        t_ok(strcmp(out.cont.done_msg, 'Traced full continuation curve in 34 continuation steps'), [t 'out.cont.done_msg']);

        t = sprintf('%s - TARGET_LAM : ', name);
        opt.stop_at = 0.7;
        p = struct('fcn', @f1p, 'x0', x0, 'opt', opt);
        [x, f, e, out, jac] = pnes_master(p);
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [-1.931782106; -1.268217894; 0.7], 6, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_is(out.cont.x(:,1), [-3;4;0], 8, [t 'out.cont.x(:,1) : x - initial']);
        t_is(out.cont.max_lam, 0.7, 10, [t 'out.cont.max_lam']);
        t_is(out.cont.iterations, 10, 12, [t 'out.cont.iterations']);
        t_is(length(out.cont.events), 1, 12, [t 'length(out.cont.events)']);
        t_is(out.cont.events.k, 10, 12, [t 'out.cont.events.k']);
        t_is(out.cont.events.idx, 1, 12, [t 'out.cont.events.idx']);
        t_ok(strcmp(out.cont.events.name, 'TARGET_LAM'), [t 'out.cont.events.name']);
        t_ok(strcmp(out.cont.done_msg, 'Reached desired lambda 0.7 in 10 continuation steps'), [t 'out.cont.done_msg']);

        t = sprintf('%s - NOSE : ', name);
        p.opt.stop_at = 'NOSE';
        [x, f, e, out, jac] = pnes_master(p);
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [-0.5; -4.75; 1.04166667], 6, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_is(out.cont.x(:,1), [-3;4;0], 8, [t 'out.cont.x(:,1) : x - initial']);
        t_is(out.cont.max_lam, 1.04166666667, 10, [t 'out.cont.max_lam']);
        t_is(out.cont.iterations, 18, 12, [t 'out.cont.iterations']);
        t_is(length(out.cont.events), 1, 12, [t 'length(out.cont.events)']);
        t_is(out.cont.events.k, 18, 12, [t 'out.cont.events.k']);
        t_is(out.cont.events.idx, 1, 12, [t 'out.cont.events.idx']);
        t_ok(strcmp(out.cont.events.name, 'NOSE'), [t 'out.cont.events.name']);
        t_ok(strcmp(out.cont.done_msg, 'Reached limit in 18 continuation steps, lambda = 1.042.'), [t 'out.cont.done_msg']);

        t = sprintf('%s - NOSE (opp dir) : ', name);
        p.x0 = [1;-1;0];
        [x, f, e, out, jac] = pnes_master(p);
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [-0.5; -4.75; 1.04166667], 5, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_is(out.cont.x(:,1), [2;-1;0], 8, [t 'out.cont.x(:,1) : x - initial']);
        t_is(out.cont.max_lam, 1.04166666667, 10, [t 'out.cont.max_lam']);
        t_is(out.cont.iterations, 20, 12, [t 'out.cont.iterations']);
        t_is(length(out.cont.events), 1, 12, [t 'length(out.cont.events)']);
        t_is(out.cont.events.k, 20, 12, [t 'out.cont.events.k']);
        t_is(out.cont.events.idx, 1, 12, [t 'out.cont.events.idx']);
        t_ok(strcmp(out.cont.events.name, 'NOSE'), [t 'out.cont.events.name']);
        t_ok(strcmp(out.cont.done_msg, 'Reached limit in 20 continuation steps, lambda = 1.042.'), [t 'out.cont.done_msg']);
    end
end

t_end;

% lam = 1;
% m = 6;
% xs = [-m:0.2:m];
% ys = xs;
% [xx, yy] = meshgrid(xs, ys);
% zz0 = zeros([size(xx), 2]);
% zz = zz0;
% for i = 1:length(xs)
%     for j = 1:length(ys)
%         zz(i, j, :) = f1p([xx(i, j); yy(i, j); lam]);
% %         zz(i, j, :) = f2([xx(i, j); yy(i, j)]);
%     end
% end
% 
% figure
% ax = gca;
% surf(ax, xx, yy, squeeze(zz(:, :, 1)))
% hold on
% surf(ax, xx, yy, squeeze(zz(:, :, 2)))
% surf(ax, m*[-1 -1; 1 1], m*[-1 1; -1 1], [0 0; 0 0])
% hold off


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
function [f, J] = f1p(x)
if nargout < 2
    f = f1(x(1:2));
else
    [f, J] = f1(x(1:2));
    J = [J [6;0]];
end
f = f + [6*x(3); 0];

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
