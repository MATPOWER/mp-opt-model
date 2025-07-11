function t_mm_solve_leqs(quiet)
% t_mm_solve_leqs - Tests of LEQ solvers via mp.opt_model.solve.

%   MP-Opt-Model
%   Copyright (c) 2010-2025, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin < 1
    quiet = 0;
end

%%  alg     name        check   opts
cfg = {
    {'',    'default',  []      []  },
    {'\',   '\',        []      []  },
    {'LU',  'LU',       []      []  },
%     {'PARDISO',  'PARDISO',       []      []  },
};

A1 = [2 1; -3 5];
b1 = [8; 1];
x1 = [3;2];
n1 = length(x1);

A2 = sparse([2 -1 0; -3 1 -2; 0 5 -4]);
b2 = [-5; 1; -7];
x2 = [-2; 1; 3];
n2 = length(x2);

ijs3 = [
	1	1	16;
	10	1	-16;
	2	2	17.064846416382252;
	8	2	-17.064846416382252;
	3	3	17.361111111111111;
	12	3	-17.361111111111111;
	4	4	16;
	16	4	-16;
	5	5	17.064846416382252;
	14	5	-17.064846416382252;
	6	6	39.995382210855354;
	7	6	-10.869565217391305;
	11	6	-11.764705882352940;
	6	7	-10.869565217391305;
	7	7	16.751918158567776;
	8	7	-5.882352941176470;
	2	8	-17.064846416382252;
	7	8	-5.882352941176470;
	8	8	32.867834278193641;
	9	8	-9.920634920634921;
	8	9	-9.920634920634921;
	9	9	23.809523809523810;
	10	9	-13.888888888888889;
	1	10	-16;
	9	10	-13.888888888888889;
	10	10	36.100069013112488;
	11	10	-6.211180124223603;
	6	11	-11.764705882352940;
	10	11	-6.211180124223603;
	11	11	17.975886006576545;
	3	12	-17.361111111111111;
	12	12	39.995382210855354;
	13	12	-10.869565217391305;
	17	12	-11.764705882352940;
	12	13	-10.869565217391305;
	13	13	16.751918158567776;
	14	13	-5.882352941176470;
	5	14	-17.064846416382252;
	13	14	-5.882352941176470;
	14	14	32.867834278193641;
	15	14	-9.920634920634921;
	14	15	-9.920634920634921;
	15	15	23.809523809523810;
	16	15	-13.888888888888889;
	4	16	-16;
	15	16	-13.888888888888889;
	16	16	36.100069013112488;
	17	16	-6.211180124223603;
	12	17	-11.764705882352940;
	16	17	-6.211180124223603;
	17	17	17.975886006576545;
];
A3 = sparse(ijs3(:, 1), ijs3(:, 2), ijs3(:, 3));
b3 = [
	1.63;
	0.85;
	0.6;
	1.63;
	0.85;
	0;
	-0.9;
	0;
	-1;
	0;
	-1.25;
	0;
	-0.9;
	0;
	-1;
	0;
	-1.25;
];
n3 = size(A3, 2);

n = 10;
diff_tool = 'bbdiff';
show_diff_on_fail = false;

reps = {};

t_begin(12+n*length(cfg), quiet);

for k = 1:length(cfg)
    alg   = cfg{k}{1};
    name  = cfg{k}{2};

    t = sprintf('%s - full A matrix : ', name);
    mm = mp.opt_model;
    mm.var.add('x', n1, zeros(n1, 1));
    mm.lin.add(mm.var, 'A', A1, b1, b1);
    if isempty(alg)
        opt = struct();
    else
        opt = struct('leq_opt', struct('solver', alg));
    end
    x = mm.solve(opt);
    t_is(x, x1, 14, [t 'x']);
    t_ok(~mm.has_parsed_soln(), [t 'has_parsed_soln() is false']);

    t = sprintf('%s - sparse (nearly-singular) A matrix : ', name);
    mm = mp.opt_model;
    mm.var.add('x', n3, zeros(n3, 1));
    mm.lin.add(mm.var, 'A', A3, b3, b3);
    if isempty(alg)
        opt = struct('leq_opt', struct('thresh', 1e9));
    else
        opt = struct('leq_opt', struct('solver', alg, 'thresh', 1e5));
    end
    warn_state = warning;
    warning('off', 'all');  %% turn of (near-)singular matrix warnings
    [x, f, e] = mm.solve(opt);
    warning(warn_state);
    t_is(f, A3*x-b3, 14, [t 'f']);
    t_is(e, 0, 14, [t 'exitflag = 0']);
    t_ok(~mm.has_parsed_soln(), [t 'has_parsed_soln() is false']);

    t = sprintf('%s - sparse A matrix : ', name);
    mm = mp.opt_model;
    mm.var.add('x', n2, zeros(n2, 1));
    mm.lin.add(mm.var, 'A12', A2(1:2, :), b2(1:2), b2(1:2));
    mm.lin.add(mm.var, 'A3', A2(3, :), b2(3), b2(3));
    opt = struct('leq_opt', struct('solver', alg), 'parse_soln', 1);
    [x, f, e, out, jac] = mm.solve(opt);
    t_is(x, x2, 14, [t 'x']);
    t_is(f, A2*x-b2, 14, [t 'f']);
    t_is(e, 1, 14, [t 'exitflag = 1']);
    t_str_match(out.alg, alg, [t 'output']);
    t_is(jac, A2, 14, [t 'jac']);
    opt.parse_soln = 0;
end

t = 'mm.soln.';
t_is(mm.soln.x, x, 14, [t 'x']);
t_is(mm.soln.f, f, 14, [t 'f']);
t_is(mm.soln.eflag, e, 14, [t 'eflag']);
t_str_match(mm.soln.output.alg, out.alg, [t 'output.alg']);
t_is(mm.soln.jac, jac, 14, [t 'jac']);

t = 'mm.var.get_soln(mm.soln, ''x'') : ';
t_is(mm.var.get_soln(mm.soln, 'x'), x, 14, [t 'x']);

t = 'mm.lin.get_soln(mm.var, mm.soln, ''A12'') : ';
f12 = mm.lin.get_soln(mm.var, mm.soln, 'A12');
t_is(f12{1}, f(1:2), 14, [t 'A12 * x - u12']);

t = 'mm.lin.get_soln(mm.var, mm.soln, ''f'', ''A3'') : ';
g = mm.lin.get_soln(mm.var, mm.soln, 'f', 'A3');
t_is(g, f(3), 14, [t 'A3 * x - u3']);

t = 'mm.lin.get_soln(mm.var, mm.soln, ''l_Ax'', ''A3'') : ';
g = mm.lin.get_soln(mm.var, mm.soln, 'l_Ax', 'A3');
t_is(g, f(3), 14, [t 'l3 - A3 * x']);

t = 'parse_soln : ';
t_ok(mm.has_parsed_soln(), [t 'has_parsed_soln() is true']);
t_is(mm.var.soln.val.x, mm.var.get_soln(mm.soln, 'x'), 14, [t 'var.val.x']);

t = 'disp_soln';
rn = fix(1e9*rand);
[pathstr, name, ext] = fileparts(which('t_opt_model'));
fname = 't_mm_solve_leqs_display_soln';
fname_e = fullfile(pathstr, 'display_soln', sprintf('%s.txt', fname));
fname_g = sprintf('%s_%d.txt', fname, rn);
[fd, msg] = fopen(fname_g, 'wt');   %% open solution file
if fd == -1
    error('t_mm_solve_leqs: could not create %d : %s', fname, msg);
end
mm.display_soln(fd);    %% write out solution
fclose(fd);
if ~t_file_match(fname_g, fname_e, t, reps, 1);
    fprintf('  compare these 2 files:\n    %s\n    %s\n', fname_g, fname_e);
    if show_diff_on_fail
        cmd = sprintf('%s %s %s', diff_tool, fname_g, fname_e);
        [status, result] = system(cmd);
        keyboard
    end
end

t_end;
