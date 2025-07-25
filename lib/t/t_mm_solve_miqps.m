function t_mm_solve_miqps(quiet)
% t_mm_solve_miqps - Tests of MIQP solvers via mp.opt_model.solve.

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

algs = {'DEFAULT', 'CPLEX', 'MOSEK', 'GUROBI', 'HIGHS', 'GLPK', 'OT'};
names = {'DEFAULT', 'CPLEX', 'MOSEK', 'Gurobi', 'HiGHS', 'glpk', 'intlin/lin/quadprog'};
check = {@have_milp_solver, 'cplex', 'mosek', 'gurobi', 'highs', 'glpk', 'intlinprog'};
does_qp   = [0 1 1 1 1 0 1];
does_miqp = [0 1 1 1 0 0 0];
if have_feature('gurobi') || have_feature('cplex') || have_feature('mosek') ...
        || have_feature('highs') || have_feature('quadprog')
    does_qp(1) = 1;
end
if have_feature('gurobi') || have_feature('cplex') || have_feature('mosek')
    does_miqp(1) = 1;
end

n = 33;
nmiqp = 20;
nmiqp_soln = 30;
diff_tool = 'bbdiff';
show_diff_on_fail = false;

reps = {};

t_begin(nmiqp_soln+n*length(algs), quiet);

diff_alg_warn_id = 'optim:linprog:WillRunDiffAlg';
if have_feature('quadprog') && have_feature('quadprog', 'vnum') == 7.005
    s1 = warning('query', diff_alg_warn_id);
    warning('off', diff_alg_warn_id);
end

for k = 1:length(algs)
    if ~isempty(check{k}) && ...
            (ischar(check{k}) && ~have_feature(check{k}) || ...
             isa(check{k}, 'function_handle') && ~check{k}())
        t_skip(n, sprintf('%s not installed', names{k}));
    else
        opt = struct('verbose', 0, 'alg', algs{k});
        mpopt = struct( ...
            'verbose', 0, ...
            'opf', struct( ...
                'violation', 5e-6 ), ...
            'cplex', struct( ...
                'lpmethod', 0, ...
                'qpmethod', 0, ...
                'opt', 0 ), ...
            'mosek', struct( ...
                'lp_alg', 0, ...
                'max_it', 0, ...
                'gap_tol', 0, ...
                'max_time', 0, ...
                'num_threads', 0, ...
                'opt', 0 ) ...
        );
%         if have_feature('highs')
%             opt.highs_opt.primal_feasibility_tolerance = 1e-10;
%             opt.highs_opt.dual_feasibility_tolerance = 1e-10;
%             opt.highs_opt.ipm_optimality_tolerance = 1e-12;
%             opt.highs_opt.primal_residual_tolerance = 1e-10;
%             opt.highs_opt.dual_residual_tolerance = 1e-10;
%         end
        if have_feature('cplex')
            % alg = 0;        %% default uses barrier method with NaN bug in lower lim multipliers
            alg = 2;        %% use dual simplex
            mpopt.cplex.lpmethod = alg;
            mpopt.cplex.qpmethod = min([4 alg]);
            opt.cplex_opt = cplex_options([], mpopt);
        end
        if have_feature('mosek')
%             sc = mosek_symbcon;
%             alg = sc.MSK_OPTIMIZER_DUAL_SIMPLEX;    %% use dual simplex
%             alg = sc.MSK_OPTIMIZER_INTPNT;          %% use interior point
%             mpopt.mosek.lp_alg = alg;
            mpopt.mosek.gap_tol = 1e-10;
%             mpopt.mosek.opts.MSK_DPAR_INTPNT_TOL_PFEAS = 1e-10;
%             mpopt.mosek.opts.MSK_DPAR_INTPNT_TOL_DFEAS = 1e-10;
%             mpopt.mosek.opts.MSK_DPAR_INTPNT_TOL_INFEAS = 1e-10;
%             mpopt.mosek.opts.MSK_DPAR_INTPNT_TOL_REL_GAP = 1e-10;
            vnum = have_feature('mosek', 'vnum');
            if vnum >= 8
%                 mpopt.mosek.opts.MSK_DPAR_INTPNT_QO_TOL_PFEAS = 1e-10;
%                 mpopt.mosek.opts.MSK_DPAR_INTPNT_QO_TOL_DFEAS = 1e-10;
%                 mpopt.mosek.opts.MSK_DPAR_INTPNT_QO_TOL_INFEAS = 1e-10;
%                 mpopt.mosek.opts.MSK_DPAR_INTPNT_QO_TOL_MU_RED = 1e-10;
                mpopt.mosek.opts.MSK_DPAR_INTPNT_QO_TOL_REL_GAP = 1e-10;
            end
            opt.mosek_opt = mosek_options([], mpopt);
        end
        opt_r = opt;
        opt_r.relax_integer = 1;

% opt.verbose = 3;
        t = sprintf('%s - 2-d ILP : ', names{k});
        %% from MOSEK 6.0 Guided Tour, section  7.13.1, https://docs.mosek.com/6.0/toolbox/node009.html
        c = [-2; -3];
        A = sparse([195 273; 4 40]);
        u = [1365; 140];
        xmax = [4; Inf];
        vtype = 'I';
        mm = mp.opt_model;
        mm.var.add('x', 2, [], [], xmax, vtype);
        mm.qdc.add(mm.var, 'c', [], c);
        mm.lin.add(mm.var, 'Ax', A, [], u);
        [x, f, s, out, lam] = mm.solve(opt);
        t_is(s, 1, 12, [t 'success']);
        t_is(x, [4; 2], 12, [t 'x']);
        t_is(f, -14, 12, [t 'f']);
        t_ok(~mm.has_parsed_soln(), [t 'has_parsed_soln() is false']);
        t = sprintf('%s - 2-d ILP (integer relaxed) : ', names{k});
        [x, f, s, out, lam] = mm.solve(opt_r);
        t_is(s, 1, 12, [t 'success']);
        t_is(x, [2.441860465; 3.255813953], 8, [t 'x']);
        t_is(f, -14.651162791, 8, [t 'f']);

        t = sprintf('%s - 6-d ILP : ', names{k});
        %% from https://doi.org/10.1109/TASE.2020.2998048
        c = [1; 2; 3; 1; 2; 3];
        A = [1 3 5 1 3 5;
             2 1.5 5 2 0.5 1];
        l = [26; 16];
        xmin = zeros(6, 1);
        xmax = 3 * ones(6, 1);
        vtype = 'I';
        mm = mp.opt_model;
        mm.var.add('x', 6, [], xmin, xmax, vtype);
        mm.qdc.add(mm.var, 'c', [], c);
        mm.lin.add(mm.var, 'Ax', A, l, u);
        [x, f, s, out, lam] = mm.solve(opt);
        t_is(s, 1, 12, [t 'success']);
        t_ok(norm(x - [1; 0; 3; 0; 0; 2], Inf) < 1e-12 || ...
             norm(x - [0; 0; 3; 1; 0; 2], Inf) < 1e-12 || ...
             norm(x - [0; 0; 3; 0; 2; 1], Inf) < 1e-12, [t 'x']);
        t_is(f, 16, 12, [t 'f']);
        t = sprintf('%s - 6-d ILP (integer relaxed) : ', names{k});
        [x, f, s, out, lam] = mm.solve(opt_r);
        t_is(s, 1, 12, [t 'success']);
        t_is([x([1;2;4;5]); x(3)+x(6)], [0; 0; 0; 0; 5.2], 7, [t 'x']);
        t_is(f, 15.6, 7, [t 'f']);

        if does_miqp(k)
            t = sprintf('%s - 4-d MIQP : ', names{k});
            %% from cplexmiqpex.m, CPLEX_Studio_Academic124/cplex/examples/src/matlab/cplexmiqpex.m
            %% Note: This is a lame example; the integer relaxed problem already
            %%       has an integer feasible solution, so this is actually just
            %%       a simple QP. -RDZ 10/29/24
            H = sparse([ 33   6    0    0;
                          6  22   11.5  0;
                          0  11.5 11    0;
                          0   0    0    0]);
            c = [-1 -2 -3 -1]';
            Aineq = [-1  1  1 10;
               1 -3  1  0];
            bineq = [20  30]';
            Aeq   = [0  1  0 -3.5];
            beq   =  0;
            xmin    = [ 0;   0;   0; 2];
            xmax    = [40; Inf; Inf; 3];
            A = sparse([Aeq; Aineq]);
            l = [beq; -Inf; -Inf];
            u = [beq; bineq];
            vtype = 'CCCI';
            mm = mp.opt_model;
            mm.var.add('x', 4, [], xmin, xmax, vtype);
            mm.qdc.add(mm.var, 'c', H, c);
            mm.lin.add(mm.var, 'Ax', A, l, u);
            [x, f, s, out, lam] = mm.solve(opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [7; 7; 0; 2], 7, [t 'x']);
            t_is(f, 1618.5, 4, [t 'f']);
            t_is(lam.mu_l, [466; 0; 0], 6, [t 'lam.mu_l']);
            t_is(lam.mu_u, [0; 272; 0], 6, [t 'lam.mu_u']);
            t_is(lam.lower, [0; 0; 349.5; 4350], 5, [t 'lam.lower']);
            t_is(lam.upper, [0; 0; 0; 0], 7, [t 'lam.upper']);
            t = sprintf('%s - 4-d MIQP (integer relaxed) : ', names{k});
            [x, f, s, out, lam] = mm.solve(opt_r);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [7; 7; 0; 2], 7, [t 'x']);
            t_is(f, 1618.5, 4, [t 'f']);
            t_is(lam.mu_l, [466; 0; 0], 6, [t 'lam.mu_l']);
            t_is(lam.mu_u, [0; 272; 0], 6, [t 'lam.mu_u']);
            t_is(lam.lower, [0; 0; 349.5; 4350], 5, [t 'lam.lower']);
            t_is(lam.upper, [0; 0; 0; 0], 7, [t 'lam.upper']);

            t = sprintf('%s - 6-d IQP : ', names{k});
            %% from Bragin, et. al. https://doi.org/10.1007/s10957-014-0561-3
            H = sparse(1:6, 1:6, [1 0.2 1 0.2 1 0.2], 6, 6);
            a = [-5 1 -5 1 -5 1];
            A = [a/5; a];
            u = [-48; -250];
            xmin = zeros(6, 1);
            vtype = 'I';
            mm = mp.opt_model;
            mm.var.add('x', 6, [], xmin, [], vtype);
            mm.qdc.add(mm.var, 'c', H, []);
            mm.lin.add(mm.var, 'Ax', A, [], u);
            [x, f, s, out, lam] = mm.solve(opt);
            t_is(s, 1, 12, [t 'success']);
            t_ok(norm(x - [16; 0; 17; 0; 17; 0], Inf) < 1e-7 || ...
                 norm(x - [17; 0; 16; 0; 17; 0], Inf) < 1e-7 || ...
                 norm(x - [17; 0; 17; 0; 16; 0], Inf) < 1e-7, [t 'x']);
            t_is(f, 417, 6, [t 'f']);
            t = sprintf('%s - 6-d IQP (integer relaxed) : ', names{k});
            [x, f, s, out, lam] = mm.solve(opt_r);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [50;0;50;0;50;0]/3, 8, [t 'x']);
            t_is(f, 1250/3, 6, [t 'f']);
        else
            t_skip(nmiqp, sprintf('%s does not handle MIQP problems', names{k}));
        end
% opt.verbose = 0;
    end
end

if have_milp_solver()
    t = 'mm.soln.';
    c = [-2; -3];
    A = sparse([195 273; 4 40]);
    u = [1365; 140];
    xmax = [4; Inf];
    vtype = 'I';
    mm = mp.opt_model;
    mm.var.add('x', 2, [], [], xmax, vtype);
    mm.qdc.add(mm.var, 'c', [], c);
    mm.lin.add(mm.var, 'Ax', A, [], u);
    opt.parse_soln = 1;
    [x, f, s, out, lam] = mm.solve(opt);
    t_is(mm.soln.x, x, 14, [t 'x']);
    t_is(mm.soln.f, f, 14, [t 'f']);
    t_is(mm.soln.eflag, s, 14, [t 'eflag']);
    t_str_match(mm.soln.output.alg, out.alg, [t 'output.alg']);
    t_is(mm.soln.lambda.lower, lam.lower, 14, [t 'mm.soln.lambda.lower']);
    t_is(mm.soln.lambda.upper, lam.upper, 14, [t 'mm.soln.lambda.upper']);
    t_is(mm.soln.lambda.mu_l, lam.mu_l, 14, [t 'mm.soln.lambda.mu_l']);
    t_is(mm.soln.lambda.mu_u, lam.mu_u, 14, [t 'mm.soln.lambda.mu_u']);

    t = 'mm.var.get_soln(mm.soln, ''x'') : ';
    [x1, mu_l, mu_u] = mm.var.get_soln(mm.soln, 'x');
    t_is(x1, x, 14, [t 'x']);
    t_is(mu_l, lam.lower, 14, [t 'mu_l']);
    t_is(mu_u, lam.upper, 14, [t 'mu_u']);

    t = 'mm.var.get_soln(mm.soln, ''mu_u'', ''x'') : ';
    t_is(mm.var.get_soln(mm.soln, 'mu_u', 'x'), lam.upper, 14, [t 'mu_l']);

    t = 'mm.lin.get_soln(mm.var, mm.soln, ''Ax'') : ';
    [g, mu_u] = mm.lin.get_soln(mm.var, mm.soln, 'Ax');
    t_is(g{1}, A*x-u, 14, [t 'A * x - u']);
    t_ok(all(isinf(g{2}) & g{2} < 0), [t 'l - A * x']);
    t_is(mu_u, lam.mu_u, 14, [t 'mu_u']);

    t = 'mm.lin.get_soln(mm.var, mm.soln, {''mu_u'', ''mu_l'', ''g''}, ''Ax'') : ';
    [mu_u, mu_l, g] = mm.lin.get_soln(mm.var, mm.soln, {'mu_u', 'mu_l', 'g'}, 'Ax');
    t_is(g{1}, A*x-u, 14, [t 'A * x - u']);
    t_ok(all(isinf(g{2}) & g{2} < 0), [t 'l - A * x']);
    t_is(mu_l, lam.mu_l, 14, [t 'mu_l']);
    t_is(mu_u, lam.mu_u, 14, [t 'mu_u']);

    t = 'mm.qdc.get_soln(mm.var, mm.soln, ''c'') : ';
    [f1, df, d2f] = mm.qdc.get_soln(mm.var, mm.soln, 'c');
    t_is(sum(f1), f, 14, [t 'f']);
    t_is(df, c, 14, [t 'df']);
    t_is(d2f, zeros(2,1), 14, [t 'd2f']);

    t = 'mm.qdc.get_soln(mm.var, mm.soln, ''df'', ''c'') : ';
    df = mm.qdc.get_soln(mm.var, mm.soln, 'df', 'c');
    t_is(df, c, 14, [t 'df']);

    t = 'parse_soln : ';
    t_ok(mm.has_parsed_soln(), [t 'has_parsed_soln() is true']);
    t_is(mm.var.soln.val.x, mm.var.get_soln(mm.soln, 'x'), 14, [t 'var.val.x']);
    t_is(mm.var.soln.mu_l.x, mm.var.get_soln(mm.soln, 'mu_l', 'x'), 14, [t 'var.mu_l.x']);
    t_is(mm.var.soln.mu_u.x, mm.var.get_soln(mm.soln, 'mu_u', 'x'), 14, [t 'var.mu_u.x']);
    t_is(mm.lin.soln.mu_l.Ax, mm.lin.get_soln(mm.var, mm.soln, 'mu_l', 'Ax'), 14, [t 'lin.mu_l.Ax']);
    t_is(mm.lin.soln.mu_u.Ax, mm.lin.get_soln(mm.var, mm.soln, 'mu_u', 'Ax'), 14, [t 'lin.mu_u.Ax']);

    t = 'disp_soln';
    rn = fix(1e9*rand);
    [pathstr, name, ext] = fileparts(which('t_opt_model'));
    fname = 't_mm_solve_miqps_display_soln';
    fname_e = fullfile(pathstr, 'display_soln', sprintf('%s.txt', fname));
    fname_g = sprintf('%s_%d.txt', fname, rn);
    [fd, msg] = fopen(fname_g, 'wt');   %% open solution file
    if fd == -1
        error('t_mm_solve_miqps: could not create %d : %s', fname, msg);
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
else
    t_skip(nmiqp_soln, 'no MILP/MIQP solver installed');
end

if have_feature('quadprog') && have_feature('quadprog', 'vnum') == 7.005
    warning(s1.state, diff_alg_warn_id);
end

t_end;

function TorF = have_milp_solver()
TorF = have_feature('cplex') || have_feature('highs') || ...
    have_feature('glpk') || have_feature('gurobi') || ...
    have_feature('intlinprog') || have_feature('mosek');

