function rv = have_feature(tag, rtype)
%HAVE_FEATURE  Test for optional functionality / version info.
%   TORF = HAVE_FEATURE(TAG)
%   TORF = HAVE_FEATURE(TAG, TOGGLE)
%   VER_STR = HAVE_FEATURE(TAG, 'vstr')
%   VER_NUM = HAVE_FEATURE(TAG, 'vnum')
%   DATE    = HAVE_FEATURE(TAG, 'date')
%   INFO    = HAVE_FEATURE(TAG, 'all')
%   HAVE_FEATURE(TAG, 'clear_cache')
%   HAVE_FEATURE('all', 'clear_cache')
%
%   Returns availability, version and release information for optional
%   functionality. All information is cached, and the cached values
%   returned on subsequent calls. If the functionality exists, an attempt is
%   made to determine the release date and version number. The second
%   argument defines which value is returned, as follows:
%       <none>      1 = optional functionality is available, 0 = not available
%       'vstr'      version number as a string (e.g. '3.11.4')
%       'vnum'      version number as numeric value (e.g. 3.011004)
%       'date'      release date as a string (e.g. '20-Jan-2015')
%       'all'       struct with fields named 'av' (for 'availability'), 'vstr',
%                   'vnum' and 'date', and values corresponding to the above,
%                   respectively.
%
%   For functionality that is not available, all calls with a string-valued
%   second argument will return an empty value.
%
%   Alternatively, the optional functionality specified by TAG can be toggled
%   OFF or ON by calling HAVE_FEATURE with a numeric second argument TOGGLE
%   with one of the following values:
%       0 - turn OFF the optional functionality
%       1 - turn ON the optional functionality (if available)
%      -1 - toggle the ON/OFF state of the optional functionality
%
%   Finally, passing 'clear_cache' as the second argument will cause the
%   cached information to be cleared for the specified TAG or, if the first
%   argument is 'all', for all optional functionality. When calling with
%   'clear_cache' no return value is defined.
%
%   Possible values for input TAG and their meanings:
%       bpmpd       - BP, BPMPD interior point LP/QP solver
%       clp         - CLP, LP/QP solver(https://github.com/coin-or/Clp)
%         opti_clp  - version of CLP distributed with OPTI Toolbox
%                       (https://www.inverseproblem.co.nz/OPTI/)
%       cplex       - CPLEX, IBM ILOG CPLEX Optimizer
%       fmincon     - FMINCON, solver from Optimization Toolbox
%         fmincon_ipm - FMINCON with Interior Point solver, from Opt Tbx 4.x+
%       fsolve      - FSOLVE, nonlinear equation solver from Opt Toolbox
%       glpk        - GLPK, GNU Linear Programming Kit
%       gurobi      - GUROBI, Gurobi solver (https://www.gurobi.com/)
%       intlinprog  - INTLINPROG, MILP solver from Optimization
%                     Toolbox 7.0 (R2014a)+
%       ipopt       - IPOPT, NLP solver
%                       (https://github.com/coin-or/Ipopt)
%       knitro      - Artelys Knitro, NLP solver
%                     (https://www.artelys.com/solvers/knitro/)
%         knitromatlab - Artelys Knitro, version 9.0.0+
%         ktrlink      - KNITRO, version < 9.0.0 (requires Opt Tbx)
%       linprog     - LINPROG, LP solver from Optimization Toolbox
%         linprog_ds - LINPROG with dual-simplex solver
%                       from Optimization Toolbox 7.1 (R2014b) +
%       matlab      - code is running under MATLAB, as opposed to Octave
%       mosek       - MOSEK, LP/QP solver (https://www.mosek.com/)
%       octave      - code is running under GNU Octave, as opposed to MATLAB
%       optim       - Optimization Toolbox
%       optimoptions - OPTIMOPTIONS, option setting funciton for Optim Tbx 6.3+
%       osqp        - OSQP (Operator Splitting QP) solver (https://osqp.org)
%       pardiso     - PARDISO, Parallel Sparse Direct & Iterative Linear Solver
%                       (https://pardiso-project.org)
%       quadprog    - QUADPROG, QP solver from Optimization Toolbox 2.x +
%         quadprog_ls - QUADPROG with large-scale interior point convex solver
%                       from Optimization Toolbox 6.x +
%       sdpt3       - SDPT3 SDP solver (https://github.com/sqlp/sdpt3)
%       sedumi      - SeDuMi SDP solver (http://sedumi.ie.lehigh.edu)
%       yalmip      - YALMIP SDP modeling platform (https://yalmip.github.io)
%
%     Functionality related to MATPOWER
%       e4st        - E4ST (https://e4st.com/)
%       minopf      - MINOPF, MINOPF, MINOS-based OPF solver
%       most        - MOST, MATPOWER Optimal Scheduling Tool
%       pdipmopf    - PDIPMOPF, primal-dual interior point method OPF solver
%       scpdipmopf  - SCPDIPMOPF, step-controlled PDIPM OPF solver
%       sdp_pf      - SDP_PF applications of semi-definite programming
%                     relaxation of power flow equations
%       smartmarket - RUNMARKET and friends, for running an energy auction
%       syngrid     - SynGrid, Synthetic Grid Creation for MATPOWER
%       tralmopf    - TRALMOPF, trust region based augmented Langrangian
%                     OPF solver
%
%   Examples:
%       if have_feature('minopf')
%           results = runopf(mpc, mpoption('opf.ac.solver', 'MINOPF'));
%       end

%   Private tags for internal use only:
%       catchme         - support for 'catch me' syntax in try/catch constructs
%       evalc           - support for evalc() function
%       ipopt_auxdata   - support for ipopt_auxdata(), required by 3.11 and later
%       lu_vec          - support for lu(..., 'vector') syntax
%       pardiso_legacy  - PARDISO v5, individual MEX files for factor, solve, etc
%       pardiso_object  - PARDISO v6 and later, object interface
%       regexp_split    - support for 'split' argument to regexp()
%       rithmaticker    - used for testing HAVE_FEATURE
%
%   The following calling syntaxes are also implemented to set and get the
%   entire cache struct and are used during testing only.
%       CACHE = HAVE_FEATURE('all', 'get_cache')
%       HAVE_FEATURE(CACHE, 'set_cache')

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

persistent h_f_cache;

action = 'D';                   %% detecting functionality (default)
if nargin > 1
    if isnumeric(rtype)
        action = 'T';           %% toggling functionality
        on_off = rtype;
        if on_off < 0                   %% flip the toggle
            TorF = have_feature(tag);
            on_off = ~TorF;
        end
    else
        switch lower(rtype)
            case 'get_cache'
                action = 'C';   %% getting cache
                rv = h_f_cache;
            case 'set_cache'
                action = 'C';   %% setting cache
                h_f_cache = tag;
            case 'clear_cache'
                action = 'C';   %% clearing cache
                if strcmpi(tag, 'all')  %% delete all fields
                    h_f_cache = struct();
                else                    %% delete field to force single re-check
                    if isfield(h_f_cache, tag)
                        h_f_cache = rmfield(h_f_cache, tag);
                    end
                end
        end
    end
end

if action == 'T'            %% change availability
    if on_off                   %% turn on if available
        if isfield(h_f_cache, tag)  %% delete field to force single re-check
            h_f_cache = rmfield(h_f_cache, tag);
        end
    else                        %% turn off
        if ~isfield(h_f_cache, tag)     %% not yet been checked
            TorF = have_feature(tag);   %% cache result first
        end
        h_f_cache.(tag).av = 0;         %% then turn off
    end
    TorF = have_feature(tag);           %% return availability
                                        %% (recheck if ON, cached 0 if OFF)
elseif action == 'D'        %% detect availability
    %% info not yet cached?
    if ~isfield(h_f_cache, tag)
        %%-----  determine installation status, version number, etc.  -----
        %% initialize default values
        TorF = 0;
        vstr = '';
        rdate = '';

        %% check for feature
% switch tag
%     case {'bpmpd', 'cplex', 'clp', 'opti_clp', 'e4st', 'fmincon', 'fmincon_ipm', ...
%         'fsolve', 'glpk', 'gurobi', 'intlinprog', 'ipopt', 'ipopt_auxdata', ...
%         'knitro', 'knitromatlab', 'ktrlink', 'linprog', ...
%         'linprog_ds', 'matlab', 'minopf', 'mosek', 'most', 'octave', 'optim', ...
%         'optimoptions', 'osqp', 'pardiso', 'pardiso_object', 'pardiso_legacy', ...
%         'quadprog', 'quadprog_ls', 'rithmaticker', ...
%         'smartmarket', 'syngrid', 'pdipmopf', 'scpdipmopf', 'tralmopf', ...
%         'sdp_pf', 'sdpt3', 'sedumi', 'yalmip', ...
%         'catchme', 'evalc', 'lu_vec', 'regexp_split'}
% fprintf('+++++ %s\n', tag);
        fcn = ['have_feature_' tag];
        if isempty(which(fcn))
            warning('have_feature: unknown functionality ''%s''', tag);
            vstr = 'unknown';
        else
            [TorF, vstr, rdate] = feval(fcn);
        end
%     otherwise
% fprintf('----- %s\n', tag);
%         have_fcn_legacy('all', 'clear_cache');
%         v = have_fcn_legacy(tag, 'all');
%         [TorF, vstr, rdate] = deal(v.av, v.vstr, v.date);
% end

        %% assign values to cache
        h_f_cache.(tag).av   = TorF;
        h_f_cache.(tag).vstr = vstr;
        if isempty(vstr)
            h_f_cache.(tag).vnum = [];
        else
            h_f_cache.(tag).vnum = vstr2num(vstr);
        end
        h_f_cache.(tag).date = rdate;
    end
end

%% extract desired values from cache
if action ~= 'C' || nargout
    if nargin < 2 || action == 'T'
        rv = h_f_cache.(tag).av;
    else
        switch lower(rtype)
            case 'vstr'
                rv = h_f_cache.(tag).vstr;
            case 'vnum'
                rv = h_f_cache.(tag).vnum;
            case 'date'
                rv = h_f_cache.(tag).date;
            case 'all'
                rv = h_f_cache.(tag);
        end
    end
end

function num = vstr2num(vstr)
% Converts version string to numerical value suitable for < or > comparisons
% E.g. '3.11.4' -->  3.011004
pat = '\.?(\d+)';
[s,e,tE,m,t] = regexp(vstr, pat);
b = 1;
num = 0;
for k = 1:length(t)
    num = num + b * str2num(t{k}{1});
    b = b / 1000;
end
