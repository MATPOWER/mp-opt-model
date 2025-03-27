function t_qcqps_master(quiet)
% t_qcqps_master - Tests of LP/QP/QCQP solvers via qcqps_master.

%   MP-Opt-Model
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin < 1
    quiet = 0;
end

%           1         2       3      4      5       6      7          8        9       10         11      12      13     14         15 
algs = {'DEFAULT', 'BPMPD', 'MIPS', 250, 'IPOPT', 'OT', 'FMINCON', 'CPLEX', 'MOSEK', 'GUROBI', 'KNITRO', 'CLP', 'GLPK', 'OSQP', 'KNITRO_NLP'};
names = {'DEFAULT', 'BPMPD_MEX', 'MIPS', 'sc-MIPS', 'IPOPT', 'linprog/quadprog', 'fmincon', 'CPLEX', 'MOSEK', 'GUROBI', 'KNITRO', 'CLP', 'glpk', 'OSQP', 'KNITRO_NLP'};
check = {[], 'bpmpd', [], [], 'ipopt', 'quadprog', 'fmincon', 'cplex', 'mosek', 'gurobi', 'knitro', 'clp', 'glpk', 'osqp', 'knitro'};
%                1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
does_lp       = [1 1 1 1 1 1 0 1 1 1  1  1  1  1  0];
does_qp       = [1 1 1 1 1 1 0 1 1 1  1  1  0  1  0];
does_qcqp     = [1 0 1 1 1 0 1 1 0 1  1  0  0  0  1];
does_nonconv  = [1 1 1 1 1 0 1 0 0 1  1  0  0  0  1];

nlp_feas = 1;
nlp_nonfeas = 1;
nqpconv = 4;
nqpnonconv = 0;
nqcqpconv = 2;
nqcqpnonconv = 1;
nproblems_feas = nqpconv+nqpnonconv+nqcqpconv+nqcqpnonconv;
nproblems_nonfeas = nlp_nonfeas;
ntests_feas = 9;
ntests_nonfeas = 1;
nskiplp_feas = (ntests_feas-2)*nlp_feas;
nskiplp_nonfeas = ntests_nonfeas*nlp_nonfeas;
nskipqpnonconv = ntests_feas*nqpnonconv;
nskipqp = nskipqpnonconv+ntests_feas*nqpconv;
nskipqcqpnonconv = ntests_feas*nqcqpnonconv;
nskipqcqp = nskipqcqpnonconv+ntests_feas*nqcqpconv;
nskipsolver = ntests_feas*nproblems_feas + ...
              (ntests_feas-2)*nlp_feas + ntests_nonfeas*nproblems_nonfeas;

t_begin(nskipsolver*length(algs), quiet);

diff_alg_warn_id = 'optim:linprog:WillRunDiffAlg';
if have_feature('quadprog') && have_feature('quadprog', 'vnum') == 7.005
    s1 = warning('query', diff_alg_warn_id);
    warning('off', diff_alg_warn_id);
end

for k = 1:length(algs)
    if ~isempty(check{k}) && ~have_feature(check{k})
        t_skip(nskipsolver, sprintf('%s not installed', names{k}));
    else
        opt = struct('verbose', 0, 'alg', algs{k});
        mpopt = struct( ...
            'verbose', 0, ...
            'opf', struct( ...
                'violation', 1e-6 ), ...
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
                'opt', 0 ), ...
            'gurobi', struct( ...
                 'method', -1, ...
                 'timelimit', Inf, ...
                 'threads', 0, ...
                 'opts', [], ...
                 'opt_fname', [], ...
                 'opt', 0), ...
             'knitro', struct( ...
                 'tol_x', 1e-10, ...
                 'tol_f', 1e-10, ...
                 'maxit', 0, ...
                 'opt', 0) ...
        );
        opt.mips_opt.comptol = 1e-8;
%         if have_feature('linprog')
%             opt.linprog_opt.Algorithm = 'interior-point';
%             opt.linprog_opt.Algorithm = 'active-set';
%             opt.linprog_opt.Algorithm = 'simplex';
%             opt.linprog_opt.Algorithm = 'dual-simplex';
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
        if have_feature('osqp')
            opt.osqp_opt.polish = 1;
%             opt.osqp_opt.alpha = 1;
%             opt.osqp_opt.eps_abs = 1e-8;
%             opt.osqp_opt.eps_rel = 1e-10;
%             opt.osqp_opt.eps_prim_inf = 1e-8;
%             opt.osqp_opt.eps_dual_inf = 1e-8;
        end
        if have_feature('gurobi')
            opt.grb_opt = gurobi_options([], mpopt);
            opt.grb_opt.BarQCPConvTol = 1e-8;
        end
        if have_feature('knitro')
            opt.knitro_opt = artelys_knitro_options([],  mpopt);
            opt.knitro_opt.ncvx_qcqp_init = 1;
        end

        if does_lp(k)
            t = sprintf('%s - 3-d LP : ', names{k});
            %% 1) based on example from 'doc linprog'
            b = [-5; -4; -6];
            A = [1 -1  1;
                -3  -2  -4;
                3  2  0];
            l2 = [-Inf; -42; -Inf];
            u2 = [20; Inf; 30];
            xmin = [0; 0; 0];
            xmax = [Inf; Inf; Inf];
            x0 = [];
            [x, f, s, out, lam] = qcqps_master([], b, [], [], [], [], A, l2, u2, xmin, xmax, x0, opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [0; 15; 3], 6, [t 'x']);
            t_is(f, -78, 6, [t 'f']);
            t_is(lam.mu_l, [0;1.5;0], 7, [t 'lam.mu_l']);
            t_is(lam.mu_u, [0;0;0.5], 7, [t 'lam.mu_u']);
            if strcmp(algs{k}, 'CLP') && ~have_feature('opti_clp')
                t_skip(2, [t 'lam.lower/upper : MEXCLP does not return multipliers on var bounds']);
            else
                t_is(lam.lower, [1;0;0], 7, [t 'lam.lower']);
                t_is(lam.upper, zeros(size(x)), 7, [t 'lam.upper']);
            end
            t_skip(2, [t 'lam.mu_l/u_quad: no multipliers of quad constraints for LPs'])

            %% 2) Infeasible LP problem
            t = sprintf('%s - infeasible LP : ', names{k});
            p = struct('A', sparse([1 1]), 'b', [1;1], 'u2', -1, 'xmin', [0;0], 'opt', opt);
            [x, f, s, out, lam] = qcqps_master(p);
            t_ok(s <= 0, [t 'no success']);
        else
            t_skip(nskiplp_feas+nskiplp_nonfeas, sprintf('%s does not handle LP problems', names{k}));
        end

        if does_qp(k)
            t = sprintf('%s - unconstrained 3-d convex QP : ', names{k});
            %% 3) From http://www.akiti.ca/QuadProgEx0Constr.html
            H = [5 -2 -1; -2 4 3; -1 3 5];
            b = [2; -35; -47];
            x0 = [0; 0; 0];
            [x, f, s, out, lam] = qcqps_master(H, b, [], [], [], [], [], [], [], [], [], [], opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [3; 5; 7], 7, [t 'x']);
            t_is(f, -249, 11, [t 'f']);
            t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
            t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
            t_is(lam.lower, zeros(size(x)), 13, [t 'lam.lower']);
            t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);
            t_skip(2, [t 'lam.mu_l/u_quad: no multipliers of quad constraints for QPs'])
        
            t = sprintf('%s - constrained 2-d convex QP : ', names{k});
            %% 4) Example from 'doc quadprog'
            H = [   1   -1;
                    -1  2   ];
            b = [-2; -6];
            A = [   1   1;
                    -1  2;
                    2   1   ];
            l2 = [];
            u2 = [2; 2; 3];
            xmin = [0; 0];
            xmax = [Inf; Inf];
            x0 = [];
            [x, f, s, out, lam] = qcqps_master(H, b, [], [], [], [], A, l2, u2, xmin, xmax, x0, opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [2; 4]/3, 7, [t 'x']);
            t_is(f, -74/9, 6, [t 'f']);
            t_is(lam.mu_l, [0;0;0], 13, [t 'lam.mu_l']);
            t_is(lam.mu_u, [28;4;0]/9, 4, [t 'lam.mu_u']);
            if strcmp(algs{k}, 'CLP') && ~have_feature('opti_clp')
                t_skip(2, [t 'lam.lower/upper : MEXCLP does not return multipliers on var bounds']);
            else
                t_is(lam.lower, zeros(size(x)), 7, [t 'lam.lower']);
                t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);
            end
            t_skip(2, [t 'lam.mu_l/u_quad: no multipliers of quad constraints for QPs'])

            t = sprintf('%s - constrained 4-d convex QP : ', names{k});
            %% 5) From https://v8doc.sas.com/sashtml/iml/chap8/sect12.htm
            H = [   1003.1  4.3     6.3     5.9;
                    4.3     2.2     2.1     3.9;
                    6.3     2.1     3.5     4.8;
                    5.9     3.9     4.8     10  ];
            b = zeros(4,1);
            A = [   1       1       1       1;
                    0.17    0.11    0.10    0.18    ];
            l2 = [1; 0.10];
            u2 = [1; Inf];
            xmin = zeros(4,1);
            xmax = inf(4,1);
            x0 = [1; 0; 0; 1];
            [x, f, s, out, lam] = qcqps_master(H, b, [], [], [], [], A, l2, u2, xmin, xmax, x0, opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [0; 2.8; 0.2; 0]/3, 5, [t 'x']);
            t_is(f, 3.29/3, 6, [t 'f']);
            t_is(lam.mu_l, [6.58;0]/3, 6, [t 'lam.mu_l']);
            t_is(lam.mu_u, [0;0], 13, [t 'lam.mu_u']);
            if strcmp(algs{k}, 'CLP') && ~have_feature('opti_clp')
                t_skip(2, [t 'lam.lower/upper : MEXCLP does not return multipliers on var bounds']);
            else
                t_is(lam.lower, [2.24;0;0;1.7667], 4, [t 'lam.lower']);
                t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);
            end
            t_skip(2, [t 'lam.mu_l/u_quad: no multipliers of quad constraints for this QP'])
            
            %% 6) Same previous passing a struct
            t = sprintf('%s - (struct) constrained 4-d convex QP : ', names{k});
            p = struct('H', H, 'A', A, 'l2', l2, 'u2', u2, 'xmin', xmin, 'x0', x0, 'opt', opt);
            [x, f, s, out, lam] = qcqps_master(p);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [0; 2.8; 0.2; 0]/3, 5, [t 'x']);
            t_is(f, 3.29/3, 6, [t 'f']);
            t_is(lam.mu_l, [6.58;0]/3, 6, [t 'lam.mu_l']);
            t_is(lam.mu_u, [0;0], 13, [t 'lam.mu_u']);
            if strcmp(algs{k}, 'CLP') && ~have_feature('opti_clp')
                t_skip(2, [t 'lam.lower/upper : MEXCLP does not return multipliers on var bounds']);
            else
                t_is(lam.lower, [2.24;0;0;1.7667], 4, [t 'lam.lower']);
                t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);
            end
            t_skip(2, [t 'lam.mu_l/u_quad: no multipliers of quad constraints for this QP'])
        else
            t_skip(nskipqp, sprintf('%s does not handle QP problems', names{k}));
        end

        if does_qcqp(k)
            %% 7) From https://docs.gurobi.com/projects/examples/en/current/examples/matlab/qcp.html
            t = sprintf('%s - convex 3-d QCQP with lin objective: ', names{k});
            H = [];
            B = [-1;0;0];
            Q = cell(2,1);
            Q{1} = sparse([2 0 0; 0 2 0; 0 0 -2]);
            Q{2} = sparse([2 0 0; 0 0 -2; 0 -2 0]);
            C = zeros(2,3);
            K = [0;0];
            l1 = [-inf;-inf] - K;
            u1 = [0; 0] - K;
            A = [1 1 1];
            l2 = 1;
            u2 = 1;
            xmin = zeros(3,1);
            xmax = inf(3,1);
            x0 = zeros(3,1);
            [x, f, s, out, lam] = qcqps_master(H, B, Q, C, l1, u1, A, l2, u2, xmin, xmax, x0, opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [3.91577; 1.78203; 4.3022]*1e-1, 6, [t 'x']);
            t_is(f, -0.391577, 6, [t 'f']);
            t_is(lam.lower, [0.0482;0.0092;0.1658]*1e-8, 6, [t 'lam.lower']);
            t_is(lam.upper, [0;0;0], 6, [t 'lam.upper']);
            t_is(lam.mu_l, 0, 6, [t 'lam.mu_l']);
            t_is(lam.mu_u, 0.391577, 6, [t 'lam.mu_u']);
            t_is(lam.mu_l_quad, [0; 0], 5, [t 'lam.mu_l_quad']);
            t_is(lam.mu_u_quad, [0.227544; 0.549342], 5, [t 'lam.mu_u_quad']);
            
            %% 8) From https://docs.mosek.com/latest/toolbox/examples-list.html#doc-example-file-qcqo1-m
            t = sprintf('%s - convex 3-d QCQP with quad objective: ', names{k});
            H = sparse([2 0 -1; 0 0.2 0; -1 0 2]);
            B = [0;-1;0];
            Q = {sparse([-2 0 0.2; 0 -2 0; 0.2 0 -0.2])};
            C = [1 1 1];
            K = 0;
            l1 = 1 - K;
            u1 = Inf - K;
            xmin = zeros(3,1);
            xmax = inf(3,1);
            x0 = zeros(3,1);
            [x, f, s, out, lam] = qcqps_master(H, B, Q, C, l1, u1, [], [], [], xmin, xmax, x0, opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [4.4880; 9.3192; 6.7411]*1e-1, 4, [t 'x']);
            t_is(f, -0.4918, 4, [t 'f']);
            t_is(lam.lower, [0.2228;0.1073;0.1483]*1e-10, 6, [t 'lam.lower']);
            t_is(lam.upper, [0;0;0], 6, [t 'lam.upper']);
            t_is(lam.mu_l_quad, 0.9419, 4, [t 'lam.mu_l_quad']);
            t_is(lam.mu_u_quad, 0, 4, [t 'lam.mu_l_quad']);

            %% 9) From "examples" folder of Knitro (exampleQCQP1)
            t = sprintf('%s - nonconvex 3-d QCQP : ', names{k});
            if does_nonconv(k)
                H = sparse([-2 -1 -1; -1 -4 0; -1 0 -2]);
                B = zeros(3,1);
                Q = {sparse(-2*eye(3))};
                C = spalloc(1,3,0);
                K = 0;
                l1 = -inf - K;
                u1 = -25 - K;
                A = [8 14 7];
                l2 = 56;
                u2 = 56;
                xmin = zeros(3,1);
                xmax = inf(3,1);
                x0 = [0;0;20];
                [x, f, s, out, lam] = qcqps_master(H, B, Q, C, l1, u1, A, l2, u2, xmin, xmax, x0, opt);
                t_is(s, 1, 12, [t 'success']);
                t_is(x, [0; 0; 8], 5, [t 'x']);
                t_is(f, -64, 4, [t 'f']);
                if strcmp(algs{k}, 'GUROBI')
                    t_skip(2, [t 'lam : QCQP-GUROBI version 11.0.3  does not return multipliers for nonconvex qcqp']);
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
        else
            t_skip(nskipqcqp, sprintf('%s does not handle QCQP problems', names{k}));
        end
    end
end

if have_feature('quadprog') && have_feature('quadprog', 'vnum') == 7.005
    warning(s1.state, diff_alg_warn_id);
end

t_end;