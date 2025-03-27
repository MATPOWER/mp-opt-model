function t_om_solve_qcqps(quiet)
% t_om_solve_qcqps - Tests of QCQP solvers via opt_model.solve.

%   MP-Opt-Model
%   Copyright (c) 2010-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin < 1
    quiet = 0;
end

algs = {'DEFAULT', 'KNITRO', 'GUROBI', 'IPOPT', 'KNITRO_NLP', 'FMINCON', 'MIPS'};
names = {'DEFAULT', 'KNITRO', 'GUROBI', 'IPOPT', 'KNITRO_NLP', 'FMINCON', 'MIPS'};
check = {[], 'knitro', 'gurobi', 'knitro', 'ipopt', 'knitro', 'fmincon', 'MIPS'};
does_nonconvex = [1 1 1 1 1 1 1];

nqcqpconv = 2;
nqcqpnonconv = 1;
nproblems = nqcqpconv+nqcqpnonconv;
ntests_om_out = 9;
ntests_om_mthds = 25;
nskipsolver = ntests_om_out*nproblems;
nskipqcqpnonconv = ntests_om_out*nqcqpnonconv;
show_diff_on_fail = false;

reps = {};

t_begin(nskipsolver*length(algs)+ntests_om_mthds, quiet);

diff_alg_warn_id = 'optim:linprog:WillRunDiffAlg';
if have_feature('quadprog') && have_feature('quadprog', 'vnum') == 7.005
    s1 = warning('query', diff_alg_warn_id);
    warning('off', diff_alg_warn_id);
end

for k = 1:length(algs)
    if ~isempty(check{k}) && ~have_feature(check{k})
        t_skip(n, sprintf('%s not installed', names{k}));
    else
        opt = struct('verbose', 0, 'alg', algs{k});
        mpopt = struct( ...
            'verbose', 0, ...
            'opf', struct( ...
                'violation', 1e-8 ), ...
            'knitro', struct( ...
                 'tol_x', 1e-10, ...
                 'tol_f', 1e-10, ...
                 'maxit', 0, ...
                 'opt', 0), ...
            'gurobi', struct( ...
                 'method', -1, ...
                 'timelimit', Inf, ...
                 'threads', 0, ...
                 'opts', [], ...
                 'opt_fname', [], ...
                 'opt', 0), ...
            'mosek', struct( ...
                'lp_alg', 0, ...
                'max_it', 0, ...
                'gap_tol', 0, ...
                'max_time', 0, ...
                'num_threads', 0, ...
                'opt', 0 ) ...            
        );
        if have_feature('knitro')
            opt.knitro_opt = artelys_knitro_options([],  mpopt);
            opt.knitro_opt.ncvx_qcqp_init = 1;
        end
        if have_feature('gurobi')
            opt.grb_opt = gurobi_options([], mpopt);
            opt.grb_opt.BarQCPConvTol = 1e-10;
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
    
        %% 1) From https://docs.gurobi.com/projects/examples/en/current/examples/matlab/qcp.html
        t = sprintf('%s - convex 3-d QCQP with linear objective: ', names{k});
        H = [];
        B = [-1;0;0];
        Q = cell(2,1);
        Q{1} = sparse([2 0 0; 0 2 0; 0 0 -2]);
        Q{2} = sparse([2 0 0; 0 0 -2; 0 -2 0]);
        C = zeros(2,3);
        K = [0;0];
        l1 = [-inf;-inf];
        u1 = [0; 0];
        A = [1 1 1];
        l2 = 1;
        u2 = 1;
        xmin = zeros(3,1);
        xmax = inf(3,1);
        x0 = zeros(3,1);
        om = opt_model();
        om.init_set_types();
        om.var.add('x', 3, x0, xmin, xmax);
        om.qdc.add(om.var, 'cost', H, B);
        om.qcn.add(om.var, 'QFx', l1, u1, Q, C, K);
        om.lin.add(om.var, 'Ax', A, l2, u2);
        [x, f, s, out, lam] = om.solve(opt);
        t_is(s, 1, 12, [t 'success']);
        t_is(x, [3.91577; 1.78203; 4.3022]*1e-1, 6, [t 'x']);
        t_is(f, -0.391577, 6, [t 'f']);
        t_is(lam.lower, [0.0482;0.0092;0.1658]*1e-8, 6, [t 'lam.lower']);
        t_is(lam.upper, [0;0;0], 6, [t 'lam.upper']);
        t_is(lam.mu_l, 0, 6, [t 'lam.mu_l']);
        t_is(lam.mu_u, 0.391577, 6, [t 'lam.mu_u']);
        t_is(lam.mu_l_quad, [0; 0], 5, [t 'lam.mu_l_quad']);
        t_is(lam.mu_u_quad, [0.227544; 0.549342], 5, [t 'lam.mu_u_quad']);
        
        %% 2) From https://docs.mosek.com/latest/toolbox/examples-list.html#doc-example-file-qcqo1-m
        t = sprintf('%s - convex 3-d QCQP with quad objective: ', names{k});
        H = sparse([2 0 -1; 0 0.2 0; -1 0 2]);
        B = [0;-1;0];
        Q = {sparse([-2 0 0.2; 0 -2 0; 0.2 0 -0.2])};
        C = [1 1 1];
        K = 0;
        l1 = 1;
        u1 = Inf;
        xmin = zeros(3,1);
        xmax = inf(3,1);
        x0 = zeros(3,1);
        om = opt_model();
        om.init_set_types();
        om.var.add('x', 3, x0, xmin, xmax);
        om.qdc.add(om.var, 'cost', H, B);
        om.qcn.add(om.var, 'QFx', l1, u1, Q, C, K);
        [x, f, s, out, lam] = om.solve(opt);
        t_is(s, 1, 12, [t 'success']);
        t_is(x, [4.4880; 9.3192; 6.7411]*1e-1, 4, [t 'x']);
        t_is(f, -0.4918, 4, [t 'f']);
        t_is(lam.lower, [0.2228;0.1073;0.1483]*1e-10, 6, [t 'lam.lower']);
        t_is(lam.upper, [0;0;0], 6, [t 'lam.upper']);
        t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
        t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
        t_is(lam.mu_l_quad, 0.9419, 4, [t 'lam.mu_l_quad']);
        t_is(lam.mu_u_quad, 0, 4, [t 'lam.mu_l_quad']);

        if does_nonconvex(k)
            %% 3) From "examples" folder of Knitro (exampleQCQP1)
            t = sprintf('%s - nonconvex 3-d QCQP : ', names{k});
            H = sparse([-2 -1 -1; -1 -4 0; -1 0 -2]);
            B = zeros(3,1);
            Q = {sparse(-2*eye(3))};
            C = spalloc(1,3,0);
            K = 0;
            l1 = -inf;
            u1 = -25;
            A = [8 14 7];
            l2 = 56;
            u2 = 56;
            xmin = zeros(3,1);
            xmax = inf(3,1);
            x0 = [0;0;20];
            om = opt_model();
            om.init_set_types();
            om.var.add('x', 3, x0, xmin, xmax);
            om.qdc.add(om.var, 'cost', H, B);
            om.qcn.add(om.var, 'QFx', l1, u1, Q, C, K);
            om.lin.add(om.var, 'Ax', A, l2, u2);
            [x, f, s, out, lam] = om.solve(opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [0; 0; 8], 5, [t 'x']);
            t_is(f, -64, 4, [t 'f']);
            if strcmp(algs{k}, 'GUROBI')
                t_skip(6, [t 'lam : GUROBI version 12.0 and lower does not return multipliers for nonconvex qcqp']); %% See: https://docs.gurobi.com/projects/optimizer/en/current/reference/attributes/constraintquadratic.html#qcpi
            else
                t_is(lam.lower, [10.285714;32;0], 5, [t 'lam.lower']);
                t_is(lam.upper, [0;0;0], 6, [t 'lam.upper']);
                t_is(lam.mu_l, 0, 5, [t 'lam.mu_l']);
                t_is(lam.mu_u, 2.28571, 5, [t 'lam.mu_u']);
                t_is(lam.mu_l_quad, 0, 4, [t 'lam.mu_l_quad']);
                t_is(lam.mu_u_quad, 0, 4, [t 'lam.mu_u_quad']);
            end
        else
            t_skip(nskipqcqpnonconv, sprintf('%s does not handle nonconvex QCQP problems', names{k}));
        end

    end
end

%% From https://docs.mosek.com/latest/toolbox/examples-list.html#doc-example-file-qcqo1-m
opt.alg = 'MIPS';

t = sprintf('%s - om.sonl: ', opt.alg);
H = sparse([2 0 -1; 0 0.2 0; -1 0 2]);
B = [0;-1;0];
Q = {sparse([-2 0 0.2; 0 -2 0; 0.2 0 -0.2])};
C = [1 1 1];
K = 0;
l1 = 1;
u1 = Inf;
xmin = zeros(3,1);
xmax = inf(3,1);
x0 = zeros(3,1);
om = opt_model();
om.init_set_types();
om.var.add('x', 3, x0, xmin, xmax);
om.qdc.add(om.var, 'cost', H, B);
om.qcn.add(om.var, 'QFx', l1, u1, Q, C, K);
opt.parse_soln = 1;
[x, f, s, out, lam] = om.solve(opt);
t_is(om.soln.x, x, 4, [t 'x']);
t_is(om.soln.f, f, 4, [t 'f']);
t_is(om.soln.eflag, s, 12, [t 'success']);
t_str_match(om.soln.output.alg, out.alg, [t 'output.alg']);
t_is(om.soln.lambda.lower, lam.lower, 6, [t 'lam.lower'])
t_is(om.soln.lambda.upper, lam.upper, 6, [t 'lam.upper']);
t_ok(isempty(om.soln.lambda.mu_l), [t 'lam.mu_l']);
t_ok(isempty(om.soln.lambda.mu_u), [t 'lam.mu_u']);
t_is(om.soln.lambda.mu_l_quad, lam.mu_l_quad, 4, [t 'lam.mu_l_quad']);
t_is(om.soln.lambda.mu_u_quad, lam.mu_u_quad, 4, [t 'lam.mu_l_quad']);

t = sprintf('%s - om.var.get_soln(soln, ''x'') : ', opt.alg);
[x1, mu_l, mu_u] = om.var.get_soln(om.soln, 'x');
t_is(x1, x, 4, [t 'x']);
t_is(mu_l, lam.lower, 6, [t 'mu_l']);
t_is(mu_u, lam.upper, 6, [t 'mu_u']);

t = sprintf('%s - om.var.get_soln(soln, ''mu_l'', ''x'') : ', opt.alg);
t_is(om.var.get_soln(om.soln, 'mu_l', 'x'), lam.lower, 6, [t 'mu_l']);

t = sprintf('%s - om.qcn.get_soln(var, soln, ''QFx'') : ', opt.alg);
[g, mu_l] = om.qcn.get_soln(om.var, om.soln, 'QFx');
%t_is(g{1}, 1/2*x'*Q{:}*x+C*x+K-u1, 8, [t '1/2*x''*Q*x + C*x + K - u']);
t_ok(isequal(g{1}, 1/2*x'*Q{:}*x+C*x+K-u1), [t '1/2*x''*Q*x + C*x + K - u']);
t_is(g{2}, l1-(1/2*x'*Q{:}*x+C*x+K), 8, [t 'l - (1/2*x''*Q*x + C*x + K)']);
t_is(mu_l, lam.mu_l_quad, 8, [t 'mu_l_quad']);

t = sprintf('%s - om.qcn.get_soln(var, soln, {''mu_u_quad'', ''mu_l_quad'', ''QFx_u''}, ''QFx'') : ', opt.alg);
[mu_u, mu_l, g] = om.qcn.get_soln(om.var, om.soln, {'mu_u', 'mu_l', 'QFx_u'}, 'QFx');
%t_is(g, 1/2*x'*Q{:}*x+C*x+K-u1, 8, [t '1/2*x''*Q*x + C*x + K - u']);
t_ok(isequal(g, 1/2*x'*Q{:}*x+C*x+K-u1), [t '1/2*x''*Q*x + C*x + K - u']);
t_is(mu_l, lam.mu_l_quad, 8, [t 'mu_l_quad']);
t_is(mu_u, lam.mu_u_quad, 8, [t 'mu_u_quad']);

t = sprintf('%s - parse_soln : ', opt.alg);
t_ok(om.has_parsed_soln(), [t 'has_parsed_soln() is true']);
t_is(om.var.soln.val.x, om.get_soln('var', 'x'), 14, [t 'var.val.x']);
if om.has_parsed_soln()
    t_is(om.qcn.soln.mu_l.QFx, mu_l, 14, [t 'mu_l_quad']);
    t_is(om.qcn.soln.mu_u.QFx, mu_u, 14, [t 'mu_u_quad']);
else
    t_skip(2, [t 'has_parsed_soln() is false'])
end

t = 'disp_soln';
rn = fix(1e9*rand);
[pathstr, name, ext] = fileparts(which('t_opt_model'));
fname = 't_om_solve_qcqps_display_soln';
fname_e = fullfile(pathstr, 'display_soln', sprintf('%s.txt', fname));
fname_g = sprintf('%s_%d.txt', fname, rn);
[fd, msg] = fopen(fname_g, 'wt');   %% open solution file
if fd == -1
    error('t_om_solve_qps: could not create %d : %s', fname, msg);
end
om.display_soln(fd);    %% write out solution
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