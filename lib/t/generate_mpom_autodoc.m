function generate_mpom_autodoc(install_dir)
% generate_mpom_autodoc - Generate the stubs and symlinks for Ref Manual.
% ::
%
%   generate_mpom_autodoc(install_dir)
%
% Inputs:
%   install_dir (char array) : path to the install directory for the package
%
% Creates the .rst stubs and symlinks to the source files for all functions
% and classes to be included in the Reference Manual. Creates all of the
% inputs (lists of functions and classes) to pass to generate_autodoc_stubs
% and generate_source_symlinks.

%   MP-Opt-Model
%   Copyright (c) 2023-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin < 1
    matpower_dir = '~/dev/projects/mp-opt-model/';
end

sphinx_src_dir = [matpower_dir 'docs/sphinx/source/'];

lib_classes = { ...
    '@mp_idx_manager', ...
    '@opt_model', ...
};
lib_mp_classes = { ...
    'set_manager', ...
    'sm_variable', ...
};
lib_fcns = {
    'convert_lin_constraint_multipliers', ...
    'convert_lin_constraint', ...
    'cplex_options', ...
    'glpk_options', ...
    'gurobi_options', ...
    'gurobiver', ...
    'have_feature_bpmpd', ...
    'have_feature_catchme', ...
    'have_feature_clp', ...
    'have_feature_cplex', ...
    'have_feature_evalc', ...
    'have_feature_fmincon_ipm', ...
    'have_feature_fmincon', ...
    'have_feature_fsolve', ...
    'have_feature_glpk', ...
    'have_feature_gurobi', ...
    'have_feature_intlinprog', ...
    'have_feature_ipopt_auxdata', ...
    'have_feature_ipopt', ...
    'have_feature_isequaln', ...
    'have_feature_knitro', ...
    'have_feature_knitromatlab', ...
    'have_feature_ktrlink', ...
    'have_feature_linprog_ds', ...
    'have_feature_linprog', ...
    'have_feature_mosek', ...
    'have_feature_opti_clp', ...
    'have_feature_optim', ...
    'have_feature_optimoptions', ...
    'have_feature_osqp', ...
    'have_feature_quadprog_ls', ...
    'have_feature_quadprog', ...
    'have_feature_sdpt3', ...
    'have_feature_sedumi', ...
    'have_feature_yalmip', ...
    'ipopt_options', ...
    'miqps_cplex', ...
    'miqps_glpk', ...
    'miqps_gurobi', ...
    'miqps_master', ...
    'miqps_mosek', ...
    'miqps_ot', ...
    'mosek_options', ...
    'mosek_symbcon', ...
    'mpomver', ...
    'mpopt2nleqopt', ...
    'mpopt2nlpopt', ...
    'mpopt2pneopt', ...
    'mpopt2qpopt', ...
    'nested_struct_copy', ...
    'nleqs_core', ...
    'nleqs_fd_newton', ...
    'nleqs_fsolve', ...
    'nleqs_gauss_seidel', ...
    'nleqs_master', ...
    'nleqs_newton', ...
    'nlp_consfcn', ...
    'nlp_costfcn', ...
    'nlp_hessfcn', ...
    'nlps_fmincon', ...
    'nlps_ipopt', ...
    'nlps_knitro', ...
    'nlps_master', ...
    'osqp_options', ...
    'osqpver', ...
    'pne_callback_default', ...
    'pne_callback_nose', ...
    'pne_callback_target_lam', ...
    'pne_detect_events', ...
    'pne_detected_event', ...
    'pne_event_nose', ...
    'pne_event_target_lam', ...
    'pne_pfcn_arc_len', ...
    'pne_pfcn_natural', ...
    'pne_pfcn_pseudo_arc_len', ...
    'pne_register_callbacks', ...
    'pne_register_events', ...
    'pnes_master', ...
    'qps_bpmpd', ...
    'qps_clp', ...
    'qps_cplex', ...
    'qps_glpk', ...
    'qps_gurobi', ...
    'qps_ipopt', ...
    'qps_master', ...
    'qps_mosek', ...
    'qps_osqp', ...
    'qps_ot', ...
};
lib_t_fcns = {
    'test_mp_opt_model', ...
    'nleqs_master_ex1', ...
    'nleqs_master_ex2', ...
    'nlps_master_ex1', ...
    'nlps_master_ex2', ...
    'pne_ex1', ...
    'qp_ex1', ...
    't_miqps_master', ...
    't_nested_struct_copy', ...
    't_nleqs_master', ...
    't_nlps_master', ...
    't_om_solve_leqs', ...
    't_om_solve_miqps', ...
    't_om_solve_nleqs', ...
    't_om_solve_nlps', ...
    't_om_solve_pne', ...
    't_om_solve_qps', ...
    't_opt_model', ...
    't_pnes_master', ...
    't_qps_master', ...
};
%     'generate_mpom_autodoc', ...

in = struct(...
    'class', struct(...
        'destdir', 'classes', ...
        'gh_base_url', 'https://github.com/MATPOWER/mp-opt-model/blob/master', ...
        'list', struct(...
            'mod', {'mp_opt_model', 'mp_opt_model.+mp'}, ...
            'src_path', {'lib', 'lib'}, ...
            'names', {lib_classes, lib_mp_classes} ...
        ) ...
    ), ...
    'function', struct(...
        'destdir', 'functions', ...
        'gh_base_url', 'https://github.com/MATPOWER/mp-opt-model/blob/master', ...
        'list', struct(...
            'mod', {'mp_opt_model', 'mp_opt_model'}, ...
            'src_path', {'lib', 'lib/t'}, ...
            'names', {lib_fcns, lib_t_fcns} ...
        ) ...
    ) ...
);

%% create stubs and symlinks for reference manual
generate_autodoc_stubs(in, sphinx_src_dir);
generate_source_symlinks(in, [sphinx_src_dir 'matlab-source/'], '../../../../../');
