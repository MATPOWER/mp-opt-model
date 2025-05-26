classdef opt_model < handle
% mp.opt_model - MP-Opt-Model optimization model class.
% ::
%
%   OM = MP.OPT_MODEL
%   OM = MP.OPT_MODEL(S)
%
%   This class implements the optimization model object used to encapsulate
%   a given optimization problem formulation. It allows for access to
%   optimization variables, constraints and costs in named blocks, keeping
%   track of the ordering and indexing of the blocks as variables,
%   constraints and costs are added to the problem.
%
%   Below are the list of available methods for use with the Opt Model class.
%   Please see the help on each individual method for more details:
%
%   Return index structure for variables, linear and nonlinear constraints
%   and costs:
%       get_idx
%
%   Determine model type:
%       problem_type
%
%   Solve the model, access, and display the solution:
%       solve
%       display_soln
%       parse_soln
%       is_solved
%       has_parsed_soln
%
%   Retreive user data in the model object:
%       get_userdata
%
%   Display the object (called automatically when you omit the semicolon
%   at the command-line):
%       display
%
%   Return the value of an individual field:
%       get
%
%   Make a shallow copy of the object:
%       copy
%
%   The following is the structure of the data in the Opt-Model object.
%   Each field of .idx or .data is a struct whose field names are the names
%   of the corresponding blocks of vars, constraints or costs (found in
%   order in the corresponding .order field). The description next to these
%   fields gives the meaning of the value for each named sub-field.
%   E.g. om.var.data.v0.Pg contains a vector of initial values for the 'Pg'
%   block of variables.
%
%   om
%       .var        - data for optimization variable sets that make up
%                     the full optimization variable x
%           .idx
%               .i1 - starting index within x
%               .iN - ending index within x
%               .N  - number of elements in this variable set
%           .N      - total number of elements in x
%           .NS     - number of variable sets or named blocks
%           .data   - bounds and initial value data
%               .v0 - vector of initial values
%               .vl - vector of lower bounds
%               .vu - vector of upper bounds
%               .vt - scalar or vector of variable types
%                       'C' = continuous
%                       'I' = integer
%                       'B' = binary
%           .order  - struct array of names/indices for variable
%                     blocks in the order they appear in x
%               .name   - name of the block, e.g. Pg
%               .idx    - indices for name, {2,3} => Pg(2,3)
%           .cache - cache for previously assembled aggregate parameters
%               .v0  - aggregate vector of variable initial values
%               .vl  - aggregate vector of variable lower bounds
%               .vu  - aggregate vector of variable upper bounds
%               .vt  - aggregate vector of variable types
%       .nle        - data for nonlinear equality constraints that make up the
%                     full set of nonlinear constraints ghne(x)
%           .idx
%               .i1 - starting index within ghne(x)
%               .iN - ending index within ghne(x)
%               .N  - number of elements in this constraint set
%           .N      - total number of elements in ghne(x)
%           .NS     - number of nonlinear constraint sets or named blocks
%           .data   - data for nonlinear constraints
%               .fcn - function handle for constraints/gradient evaluation
%               .hess - function handle for Hessian evaluation
%               .vs - cell array of variable sets that define the xx for
%                     this constraint block
%           .order  - struct array of names/indices for nonlinear constraint
%                     blocks in the order they appear in ghne(x)
%               .name   - name of the block, e.g. Pmis
%               .idx    - indices for name, {2,3} => Pmis(2,3)
%       .nli        - data for nonlinear inequality constraints that make up the
%                     full set of nonlinear constraints ghni(x)
%           .idx
%               .i1 - starting index within ghni(x)
%               .iN - ending index within ghni(x)
%               .N  - number of elements in this constraint set
%           .N      - total number of elements in ghni(x)
%           .NS     - number of nonlinear constraint sets or named blocks
%           .data   - data for nonlinear constraints
%               .fcn - function handle for constraints/gradient evaluation
%               .hess - function handle for Hessian evaluation
%               .vs - cell array of variable sets that define the xx for
%                     this constraint block
%           .order  - struct array of names/indices for nonlinear constraint
%                     blocks in the order they appear in ghni(x)
%               .name   - name of the block, e.g. Pmis
%               .idx    - indices for name, {2,3} => Pmis(2,3)
%       .lin        - data for linear constraints that make up the
%                     full set of linear constraints ghl(x)
%           .idx
%               .i1 - starting index within ghl(x)
%               .iN - ending index within ghl(x)
%               .N  - number of elements in this constraint set
%           .N      - total number of elements in ghl(x)
%           .NS     - number of linear constraint sets or named blocks
%           .data   - data for l <= A*xx <= u linear constraints
%               .A  - sparse linear constraint matrix
%               .l  - left hand side vector, bounding A*x below
%               .u  - right hand side vector, bounding A*x above
%               .vs - cell array of variable sets that define the xx for
%                     this constraint block
%           .order  - struct array of names/indices for linear constraint
%                     blocks in the order they appear in ghl(x)
%               .name   - name of the block, e.g. Pmis
%               .idx    - indices for name, {2,3} => Pmis(2,3)
%           .cache - cache for previously assembled aggregate parameters
%               .A  - aggregate sparse linear constraint matrix
%               .l  - aggregate left hand side vector, bounding A*x below
%               .u  - aggregate right hand side vector, bounding A*x above
%       .qdc       - quadratic costs
%           .idx
%               .i1 - starting row index within quadratic costs
%               .iN - ending row index within quadratic costs
%               .N  - number of rows in this quadratic cost block
%           .N      - total number of rows in quadratic costs
%           .NS     - number of quadratic cost blocks
%           .data   - data for each quadratic cost block
%               .Q  - sparse matrix (or vector) of quadratic cost coefficients
%               .c  - column vector of linear cost coefficients
%               .k  - constant term
%               .vs - cell array of variable sets that define xx for this
%                     quadratic cost block, where sizes of Q, c, k conform to xx
%           .order  - struct array of names/indices for quadratic cost blocks
%                     in the order the were added
%               .name   - name of the block, e.g. R
%               .idx    - indices for name, {2,3} => R(2,3)
%           .cache - cache for previously assembled aggregate parameters
%               .Q  - aggregate sparse matrix of quadratic cost coefficients
%               .c  - aggregate column vector of linear cost coefficients
%               .k  - aggregate constant term
%       .nlc       - general nonlinear costs
%           .idx
%               .i1 - starting row index within nonlinear costs
%               .iN - ending row index within nonlinear costs
%               .N  - number of rows in this nonlinear cost block
%                     (always equal to 1 for nonlinear cost blocks)
%           .N      - total number of rows in nonlinear costs
%           .NS     - number of nonlinear cost blocks
%           .data   - data for each nonlinear cost block
%               .fcn - function handle for cost, gradient, Hessian evaluation
%               .vs - cell array of variable sets that define xx for this
%                     nonlinear cost block, where xx is the input to the
%                     evaluation function
%           .order  - struct array of names/indices for nonlinear cost blocks
%                     in the order they were added
%               .name   - name of the block, e.g. R
%               .idx    - indices for name, {2,3} => R(2,3)
%       .prob_type  - used to cache the return value of problem_type()
%       .soln       - struct containing the output of the solve() method
%                     with the following fields
%           .eflag  - exit flag
%           .output - output struct with the following fields
%               .alg     - algorithm code of solver used
%               (others) - solver specific fields
%           .x      - solution vector
%           .f      - final (objective) function value
%           .jac    - final Jacobian matrix (if available, for LEQ/NLEQ probs)
%           .lambda - Lagrange and Kuhn-Tucker multipliers on constraints
%       .userdata   - any user defined data
%           .(user defined fields)
%
% See also mp.set_manager.

%   MP-Opt-Model
%   Copyright (c) 2008-2025, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

    properties
        var             %% variables
        lin             %% linear constraints
        qcn             %% quadratic constraints
        nle             %% nonlinear equality constraints
        nli             %% nonlinear inequality constraints
        qdc             %% quadratic costs
        nlc             %% general nonlinear costs
        prob_type = ''; %% problem type
        userdata = struct();

        %% results of solve()
        soln = struct( ...
            'eflag', [], ...    %% exit flag
            'output', [], ...   %% algorithm code & solver-specific fields
            'x', [], ...        %% solution vector
            'f', [], ...        %% final (objective) function value
            'jac', [], ...      %% Jacobian (if available) for LEQ/NLEQ
            'lambda', [] );     %% constraint shadow prices
    end     %% properties

    methods
        function om = opt_model(varargin)
            % Constructor.
            % ::
            %
            %   om = mp.opt_model()
            %   om = mp.opt_model(a_struct)
            %   om = mp.opt_model(an_om)

            if nargin > 0
                s = varargin{1};
                if isa(s, 'mp.opt_model') || isa(s, 'opt_model')
                    if have_feature('octave')
                        s1 = warning('query', 'Octave:classdef-to-struct');
                        warning('off', 'Octave:classdef-to-struct');
                    end
                    props = fieldnames(s);
                    if have_feature('octave')
                        warning(s1.state, 'Octave:classdef-to-struct');
                    end
                    [~, k] = ismember('set_types', props);
                    if k
                        props(k) = [];  %% remove 'set_types'
                    end
                    for k = 1:length(props)
                        om = copy_prop(s, om, props{k});
                    end
                elseif isstruct(s)
                    props = fieldnames(om);
                    for k = 1:length(props)
                        if isfield(s, props{k})
                            om = copy_prop(s, om, props{k});
                        end
                    end
                else
                    error('mp.opt_model: input must be an ''mp.opt_model'' object or a struct');
                end
            end

            if isempty(om.var)      %% skip for copy constructor
                om.var = mp.sm_variable('VARIABLES');
                om.lin = mp.sm_lin_constraint('LINEAR CONSTRAINTS');
                om.qcn = mp.sm_quad_constraint('QUADRATIC CONSTRAINTS');
                om.nle = mp.sm_nln_constraint('NONLIN EQ CONSTRAINTS');
                om.nli = mp.sm_nln_constraint('NONLIN INEQ CONSTRAINTS');
                om.qdc = mp.sm_quad_cost('QUADRATIC COSTS');
                om.nlc = mp.sm_nln_cost('GEN NONLIN COSTS');
            end
        end

        function st = get_set_types(om)
            % List of names of properties of set types managed by this class.
            % ::
            %
            %   st = om.get_set_types();
            %
            % Output:
            %   st (cell array) : list of set types, namely:
            %       ``{'var', 'lin', 'qcn', 'nle', 'nli', 'qdc', 'nlc'}``

            st = {'var', 'lin', 'qcn', 'nle', 'nli', 'qdc', 'nlc'};
        end

        function new_om = copy(om)
            % Duplicate the object.

            %% delete old 'params' (cached parameters) fields
            %% to avoid clash with newer params() method
            fn = om.get_set_types();
            for f = 1:length(fn)
                if isfield(om.(fn{f}), 'params')
                    om.(fn{f}) = rmfield(om.(fn{f}), 'params');
                end
            end

            %% initialize copy
            new_om = eval(class(om));   %% create new object

            %% copy properties/fields
            if have_feature('octave')
                s1 = warning('query', 'Octave:classdef-to-struct');
                warning('off', 'Octave:classdef-to-struct');
            end
            props = fieldnames(om);
            if have_feature('octave')
                warning(s1.state, 'Octave:classdef-to-struct');
            end
            for k = 1:length(props)
                new_om = copy_prop(om, new_om, props{k});
            end
        end

        function varargout = get_idx(om, varargin)
            % get_idx - Returns the idx struct for vars, lin/nonlin constraints, costs.
            % ::
            %
            %   VV = OM.GET_IDX()
            %   [VV, LL] = OM.GET_IDX()
            %   [VV, LL, NNE] = OM.GET_IDX()
            %   [VV, LL, NNE, NNI] = OM.GET_IDX()
            %   [VV, LL, NNE, NNI, QQ] = OM.GET_IDX()
            %   [VV, LL, NNE, NNI, QQ, NNC] = OM.GET_IDX()
            %   [VV, LL, NNE, NNI, QQ, NNC, QQCN] = OM.GET_IDX()
            %
            %   Returns a structure for each with the beginning and ending
            %   index value and the number of elements for each named block.
            %   The 'i1' field (that's a one) is a struct with all of the
            %   starting indices, 'iN' contains all the ending indices and
            %   'N' contains all the sizes. Each is a struct whose fields are
            %   the named blocks.
            %
            %   Alternatively, you can specify the type of named set(s) directly
            %   as inputs ...
            %
            %   [IDX1, IDX2, ...] = OM.GET_IDX(SET_TYPE1, SET_TYPE2, ...);
            %   VV = OM.GET_IDX('var');
            %   [LL, NNE, NNI] = OM.GET_IDX('lin', 'nle', 'nli');
            %
            %   The specific type of named set being referenced is
            %   given by the SET_TYPE inputs, with the following valid options:
            %       SET_TYPE = 'var'   => variable set
            %       SET_TYPE = 'lin'   => linear constraint set
            %       SET_TYPE = 'qcn'   => quadratic constraint set
            %       SET_TYPE = 'nle'   => nonlinear equality constraint set
            %       SET_TYPE = 'nli'   => nonlinear inequality constraint set
            %       SET_TYPE = 'qdc'   => quadratic cost set
            %       SET_TYPE = 'nnc'   => nonlinear cost set
            %
            %   Examples:
            %       [vv, ll, nne] = om.get_idx();
            %       [vv, ll, qq] = om.get_idx('var', 'lin', 'qdc');
            %
            %       For a variable block named 'z' we have ...
            %           vv.i1.z - starting index for 'z' in optimization vector x
            %           vv.iN.z - ending index for 'z' in optimization vector x
            %           vv.N.z  - number of elements in 'z'
            %
            %       To extract a 'z' variable from x:
            %           z = x(vv.i1.z:vv.iN.z);
            %
            %       To extract the multipliers on a linear constraint set
            %       named 'foo', where mu_l and mu_u are the full set of
            %       linear constraint multipliers:
            %           mu_l_foo = mu_l(ll.i1.foo:ll.iN.foo);
            %           mu_u_foo = mu_u(ll.i1.foo:ll.iN.foo);
            %
            %       The number of nonlinear equality constraints in a set named 'bar':
            %           nbar = nne.N.bar;
            %         (note: the following is preferable ...
            %           nbar = om.nle.get_N('bar');
            %         ... if you haven't already called get_idx to get nne.)
            %
            %       If 'z', 'foo' and 'bar' are indexed sets, then you can
            %       replace them with something like 'z(i,j)', 'foo(i,j,k)'
            %       or 'bar(i)' in the examples above.

            if nargin == 1
                varargout{1} = om.var.idx;
                if nargout > 1
                    varargout{2} = om.lin.idx;
                    if nargout > 2
                        varargout{3} = om.nle.idx;
                        if nargout > 3
                            varargout{4} = om.nli.idx;
                            if nargout > 4
                                varargout{5} = om.qdc.idx;
                                if nargout > 5
                                    varargout{6} = om.nlc.idx;
                                    if nargout > 6
                                        varargout{7} = om.qcn.idx;
                                    end
                                end
                            end
                        end
                    end
                end
            else
                for k = nargout:-1:1
                    varargout{k} = om.(varargin{k}).idx;
                end
            end
        end

        function rv = get_userdata(om, name)
            % get_userdata - Used to retrieve values of user data.
            % ::
            %
            %   VAL = OM.GET_USERDATA(NAME) returns the value specified by the given name
            %   or an empty matrix if userdata with NAME does not exist.
            %
            %   This function allows the user to retrieve any arbitrary data that was
            %   saved in the object for later use. Data for a given NAME is saved by
            %   assigning it to OM.userdata.(NAME).
            %
            %   This can be useful, for example, when using a user function to add
            %   variables or constraints, etc. Suppose some special indexing is
            %   constructed when adding some variables or constraints. This indexing data
            %   can be stored and used later to "unpack" the results of the solved case.
            %
            % See also mp.opt_model.

            if isfield(om.userdata, name)
                rv = om.userdata.(name);
            else
                rv = [];
            end
        end

        function prob = problem_type(om, recheck)
            % problem_type - Return a string identifying the type of mathematical program
            % ::
            %
            %   PROB_TYPE = OM.PROBLEM_TYPE()
            %   PROB_TYPE = OM.PROBLEM_TYPE(RECHECK)
            %
            %   Returns a string identifying the type of mathematical program
            %   represented by the current model, based on the variables, costs,
            %   and constraints that have been added to the model. Used to
            %   automatically select an appropriate solver.
            %
            %   Linear and nonlinear equations are models with no costs, no inequality
            %   constraints, and an equal number of continuous variables and equality
            %   constraints. If the number of variables in a nonlinear equation model
            %   is one more than the number of constraints, it is a parameterized
            %   nonlinear equation.
            %
            %   Outputs:
            %       PROB_TYPE : problem type, one of the following strings:
            %           'LEQ'   - linear equations
            %           'NLEQ'  - nonlinear equations
            %           'PNE'   - parameterized nonlinear equations
            %           'LP'    - linear program
            %           'QP'    - quadratic program
            %           'NLP'   - nonlinear program
            %           'MILP'  - mixed-integer linear program
            %           'MIQP'  - mixed-integer quadratic program
            %           'MINLP' - mixed-integer nonlinear program
            %
            %   The output value is cached for future calls, but calling with a true
            %   value for the optional RECHECK argument will force it to recheck in
            %   case the problem type has changed due to modifying the variables,
            %   constraints or costs in the model.
            %
            % See also mp.opt_model.
            
            if isempty(om.prob_type) || nargin > 1 && recheck
                nleN = om.nle.get_N();      %% nonlinear equalities
                nliN = om.nli.get_N();      %% nonlinear inequalities
                nlcN = om.nlc.get_N();      %% general nonlinear costs
                qdcN = om.qdc.get_N();      %% quadratic costs
                qcnN = om.qcn.get_N();      %% quadratic constraints
                linN = om.lin.get_N();      %% linear constraints
                varN = om.var.get_N();      %% variables
                if varN == 0
                    prob = '';
                elseif nlcN || qdcN         %% problem has costs
                    if nliN || nleN || nlcN %% nonlinear
                        prob = 'NLP';           %% nonlinear program
                    elseif qcnN                 %% quadratically constrained quadratic program
                        prob = 'QCQP';
                    else                    %% linear constraints, no general nonlinear costs
                        %% get quadratic cost coefficients
                        H = om.qdc.params(om.var);
                        if isempty(H) || ~any(any(H))
                            prob = 'LP';        %% linear program
                        else
                            prob = 'QP';        %% quadratic program
                        end
                    end
                else                    %% problem has no costs
                    if nliN
                        error('mp.opt_model.problem_type: invalid problem - nonlinear inequality constraints with no costs');
                    end
                    if nleN + qcnN + linN == varN || nleN + qcnN + linN == varN - 1   %% square (or almost) system
                        if linN > 0
                            %% get lower & upper bounds
                            [A, l, u] = om.lin.params(om.var);
                            if any(l ~= u)
                                error('mp.opt_model.problem_type: invalid problem - linear inequality constraints with no costs');
                            end
                        end
                        if nleN + qcnN + linN == varN  %% square system
                            if (nleN + qcnN) ~= 0
                                prob = 'NLEQ';      %% square nonlinear set of equations
                            else
                                prob = 'LEQ';       %% square linear set of equations
                            end
                        elseif nleN + linN + 1 == varN  %% square + 1 extra (parameterization) variable
                            if nleN
                                prob = 'PNE';       %% parameterized nonlinear set of equations
                            else
                                prob = 'PLEQ';      %% parameterized linear set of equations
                                error('mp.opt_model.problem_type: invalid problem - PNE not implemented for for linear constraints only');
                            end
                        else
                            error('mp.opt_model.problem_type: invalid problem - PNE must have num of vars = num of constraints + 1');
                        end
                    else
                        error('mp.opt_model.problem_type: invalid problem - non-square system with no costs');
                    end
                end
                if om.is_mixed_integer() && ~strcmp(prob, 'NLEQ')
                    prob = ['MI' prob];
                end
                om.prob_type = prob;    %% cache it
            else
                prob = om.prob_type;    %% return cached type
            end
        end

        function TorF = is_mixed_integer(om)
            % is_mixed_integer - Return true if model is mixed integer, false otherwise.
            % ::
            %
            %   TorF = OM.IS_MIXED_INTEGER()
            %
            %   Outputs:
            %       TorF : 1 or 0, indicating whether any of the variables are
            %              binary or integer
            %
            % See also mp.opt_model.

            TorF = 0;
            if om.var.get_N()
                for k = 1:length(om.var.order)
                    t = om.var.data.vt.(om.var.order(k).name);
                    if iscell(t)
                        for j = 1:length(t(:))
                            if any(t{j} ~= 'C')
                                TorF = 1;
                                break;
                            end
                        end
                    else
                        if any(t ~= 'C')
                            TorF = 1;
                            break;
                        end
                    end
                end
            end
        end

        function TorF = is_solved(om)
            % Return true if model has been solved.

            TorF = ~isempty(om.soln.eflag);
        end

        function [x, f, eflag, output, lambda] = solve(om, opt)
            % solve - Solve the optimization model.
            % ::
            %
            %   X = OM.SOLVE()
            %   [X, F] = OM.SOLVE()
            %   [X, F, EXITFLAG] = OM.SOLVE()
            %   [X, F, EXITFLAG, OUTPUT] = OM.SOLVE()
            %   [X, F, EXITFLAG, OUTPUT, JAC] = OM.SOLVE()      (LEQ/NLEQ problems)
            %   [X, F, EXITFLAG, OUTPUT, LAMBDA] = OM.SOLVE()   (other problem types)
            %   [X ...] = OM.SOLVE(OPT)
            %
            %   Solves the optimization model using one of the following, depending
            %   on the problem type: QPS_MASTER, MIQPS_MASTER, NLPS_MASTER, NLEQS_MASTER.
            %
            %   Inputs:
            %       OPT : optional options structure with the following fields,
            %           all of which are also optional (default values shown in
            %           parentheses)
            %           alg ('DEFAULT') : determines which solver to use, list of relevant
            %                   problem types are listed in parens next to each
            %               'DEFAULT' : automatic, depending on problem type, uses the
            %                       the first available of:
            %                   LP - Gurobi, CPLEX, MOSEK, linprog (if MATLAB), HIGHS, GLPK,
            %                           BPMPD, MIPS
            %                   QP - Gurobi, CPLEX, MOSEK, quadprog (if MATLAB), HIGHS,
            %                           BPMPD, MIPS
            %                   MILP - Gurobi, CPLEX, MOSEK, Opt Tbx (intlingprog), HIGHS,
            %                           GLPK
            %                   MIQP - Gurobi, CPLEX, MOSEK
            %                   NLP - MIPS
            %                   MINLP - Artelys Knitro (not yet implemented)
            %                   LEQ - built-in backslash operator
            %                   NLEQ - Newton's method
            %               'BPMPD'   : (LP, QP) BPMPD_MEX
            %               'CLP'     : (LP, QP) CLP
            %               'CPLEX'   : (LP, QP, MILP, MIQP) CPLEX
            %               'FD'      : (NLEQ) fast-decoupled Newon's method
            %               'FMINCON' : (NLP) FMINCON, MATLAB Optimization Toolbox
            %               'FSOLVE'  : (NLEQ) FSOLVE, MATLAB Optimization Toolbox
            %               'GLPK'    : (LP, MILP) GLPK
            %               'GS'      : (NLEQ) Gauss-Seidel
            %               'GUROBI'  : (LP, QP, MILP, MIQP) Gurobi
            %               'HIGHS'   : (LP, QP, MILP) HiGHS, https://highs.dev
            %               'IPOPT'   : (LP, QP, NLP) IPOPT, requires MEX interface to IPOPT solver
            %                           https://github.com/coin-or/Ipopt
            %               'KNITRO'  : (NLP, MINLP) Artelys Knitro, requires Artelys Knitro solver
            %                           https://www.artelys.com/solvers/knitro/
            %               'MIPS'    : (LP, QP, NLP) MIPS, MATPOWER Interior Point Solver
            %                        pure MATLAB implementation of a primal-dual
            %                        interior point method, if mips_opt.step_control = 1
            %                        it uses MIPS-sc, a step controlled variant of MIPS
            %               'MOSEK'   : (LP, QP, MILP, MIQP) MOSEK
            %               'NEWTON'  : (NLEQ) Newton's method
            %               'OSQP'    : (LP, QP) OSQP, https://osqp.org
            %               'OT'      : (LP, QP, MILP) MATLAB Optimization Toolbox,
            %                           LINPROG, QUADPROG or INTLINPROG
            %           verbose (0) - controls level of progress output displayed
            %               0 = no progress output
            %               1 = some progress output
            %               2 = verbose progress output
            %           bp_opt      - options vector for BP (BPMPD)
            %           clp_opt     - options vector for CLP
            %           cplex_opt   - options struct for CPLEX
            %           fd_opt      - options struct for fast-decoupled Newton's method,
            %                           nleqs_fd_newton()
            %           fmincon_opt - options struct for FMINCON
            %           fsolve_opt  - options struct for FSOLVE
            %           glpk_opt    - options struct for GLPK
            %           grb_opt     - options struct for GUROBI
            %           gs_opt      - options struct for Gauss-Seidel method,
            %                           nleqs_gauss_seidel()
            %           highs_opt   - options struct for HIGHS
            %           intlinprog_opt - options struct for INTLINPROG
            %           ipopt_opt   - options struct for IPOPT
            %           knitro_opt  - options struct for Artelys Knitro
            %           leq_opt     - options struct for MPLINSOLVE, with optional fields
            %               'solver' and 'opt' corresponding to respective MPLINSOLVE args,
            %               and 'thresh' specifying a threshold on the absolute value of
            %               any element X, above which EXITFLAG will be set to 0
            %           linprog_opt - options struct for LINPROG
            %           mips_opt    - options struct for MIPS
            %           mosek_opt   - options struct for MOSEK
            %           newton_opt  - options struct for Newton method, NLEQS_NEWTON
            %           osqp_opt    - options struct for OSQP
            %           quadprog_opt - options struct for QUADPROG
            %           parse_soln (0) - flag that specifies whether or not to call
            %               the PARSE_SOLN method and place the return values in OM.soln.
            %           price_stage_warn_tol (1e-7) - tolerance on the objective fcn
            %               value and primal variable relative match required to avoid
            %               mis-match warning message if mixed integer price computation
            %               stage is not skipped
            %           relax_integer (0) - relax integer constraints, if true
            %           skip_prices (0) - flag that specifies whether or not to skip the
            %               price computation stage for mixed integer problems, in which
            %               the problem is re-solved for only the continuous variables,
            %               with all others being constrained to their solved values
            %           x0 (empty)  - optional initial value of x, overrides value
            %               stored in model (ignored by some solvers)
            %
            %   Outputs:
            %       X : solution vector
            %       F : final (objective) function value
            %       EXITFLAG : exit flag
            %           1 = converged
            %           0 or negative values = solver specific failure codes
            %       OUTPUT : output struct with the following fields:
            %           alg - algorithm code of solver used
            %           et  - elapsed time (sec)
            %           (others) - solver specific fields
            %       JAC : final Jacobian matrix (if available, for LEQ/NLEQ problems)
            %       LAMBDA : (for all non-NLEQ problem types) struct containing the
            %           Langrange and Kuhn-Tucker multipliers on the constraints, with
            %           fields:
            %           eqnonlin - nonlinear equality constraints
            %           ineqnonlin - nonlinear inequality constraints
            %           mu_l - lower (left-hand) limit on linear constraints
            %           mu_u - upper (right-hand) limit on linear constraints
            %           lower - lower bound on optimization variables
            %           upper - upper bound on optimization variables
            %
            % See also mp.opt_model, qps_master, miqps_master, nlps_master, nleqs_master,
            % pnes_master, mp_linsolve.
            
            t0 = tic;       %% start timer
            if nargin < 2
                opt = struct();
            end
            % opt.parse_soln = 1;
            
            %% call appropriate solver
            pt = om.problem_type();
            switch pt
                case 'LEQ'          %% LEQ   - linear equations
                    if isfield(opt, 'leq_opt')
                        if isfield(opt.leq_opt, 'solver')
                            leq_solver = opt.leq_opt.solver;
                        else
                            leq_solver = '';
                        end
                        if isfield(opt.leq_opt, 'opt')
                            leq_opt = opt.leq_opt.opt;
                        else
                            leq_opt = struct();
                        end
                        if isfield(opt.leq_opt, 'thresh')
                            leq_thresh = opt.leq_opt.thresh;
                        else
                            leq_thresh = 0;
                        end
                    else
                        leq_solver = '';
                        leq_opt = struct();
                        leq_thresh = 0;
                    end
            
                    [A, b] = om.lin.params(om.var);
                    if leq_thresh           %% check for failure
                        %% set up to trap non-singular matrix warnings
                        [lastmsg, lastid] = lastwarn;
                        lastwarn('');
            
                        x = mplinsolve(A, b, leq_solver, leq_opt);
            
                        [msg, id] = lastwarn;
                        %% Octave is not consistent in assigning proper warning id,
                        %% so we just check for presence of *any* warning
                        if ~isempty(msg) || max(abs(x)) > leq_thresh
                            eflag = 0;
                        else
                            eflag = 1;
                        end
                    else                    %% no failure check
                        x = mplinsolve(A, b, leq_solver, leq_opt);
                        eflag = 1;
                    end
                    f = A*x - b;
                    output = struct('alg', leq_solver);
                    lambda = A;     %% jac
                case {'NLEQ', 'PNE'}    %% NLEQ, PNE - nonlinear equations
                    if isfield(opt, 'x0')
                        x0 = opt.x0;
                    else
                        x0 = om.var.params();
                    end
            
                    fcn = @(x)nleq_fcn_(om, x);
                    switch pt
                        case 'NLEQ' %% NLEQ - nonlinear equation
                            [x, f, eflag, output, lambda] = nleqs_master(fcn, x0, opt);
                        case 'PNE'  %% PNE - parameterized nonlinear equation
                            [x, f, eflag, output, lambda] = pnes_master(fcn, x0, opt);
                    end
                case {'MINLP', 'NLP'}
                    mixed_integer = strcmp(pt(1:2), 'MI') && ...
                        (~isfield(opt, 'relax_integer') || ~opt.relax_integer);
                    if mixed_integer    %% MINLP - mixed integer non-linear program
                        error('mp.opt_model.solve: not yet implemented for ''MINLP'' problems.')
                    else                %% NLP   - nonlinear program
                        %% optimization vars, bounds, types
                        [x0, xmin, xmax] = om.var.params();
                        if isfield(opt, 'x0')
                            x0 = opt.x0;
                        end
            
                        %% run solver
                        [A, l, u] = om.lin.params(om.var);
                        f_fcn = @(x)nlp_costfcn(om, x);
                        gh_fcn = @(x)nlp_consfcn(om, x);
                        hess_fcn = @(x, lambda, cost_mult)nlp_hessfcn(om, x, lambda, cost_mult);
                        [x, f, eflag, output, lambda] = ...
                            nlps_master(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
                    end
                otherwise
                    %% get parameters
                    [HH, CC, C0] = om.qdc.params(om.var);
                    [Q, B, ll, uu] = om.qcn.params(om.var);
                    [A, l, u] = om.lin.params(om.var);
                    mixed_integer = strcmp(pt(1:2), 'MI') && ...
                        (~isfield(opt, 'relax_integer') || ~opt.relax_integer);
            
                    if mixed_integer
                        %% optimization vars, bounds, types
                        [x0, xmin, xmax, vtype] = om.var.params();
                        if isfield(opt, 'x0')
                            x0 = opt.x0;
                        end
            
                        %% run solver
                        if isempty(Q)          %% MILP, MIQP - mixed integer linear/quadratic program
                            [x, f, eflag, output, lambda] = ...
                                miqps_master(HH, CC, A, l, u, xmin, xmax, x0, vtype, opt);
                        else                   %% MIQCQP - mixed integer quadratically constrained quadratic program
                            % To be implemented ...
                            % [x, f, eflag, output, lambda] = ...
                            %    miqcqps_master(HH, CC, Q, B, k, ll, uu, A, l, u, xmin, xmax, x0, vtype, opt);
                        end
                    else                %% LP, QP - linear/quadratic program
                        %% optimization vars, bounds, types
                        [x0, xmin, xmax] = om.var.params();
                        if isfield(opt, 'x0')
                            x0 = opt.x0;
                        end
            
                        %% run solver
                        if isempty(Q)          %% LP, QP - linear/quadratic program
                            [x, f, eflag, output, lambda] = ...
                                qps_master(HH, CC, A, l, u, xmin, xmax, x0, opt);
                        else                   %% QCQP - quadratically constrained quadratic program
                            [x, f, eflag, output, lambda] = ...
                                qcqps_master(HH, CC, Q, B, ll, uu, A, l, u, xmin, xmax, x0, opt);
                        end
                    end
                    f = f + C0;
            end
            
            %% store solution
            om.soln.eflag = eflag;
            om.soln.x = x;
            om.soln.f = f;
            om.soln.output = output;
            if isstruct(lambda)
                om.soln.lambda = lambda;
            else
                om.soln.jac = lambda;
            end
            
            %% parse solution
            if isfield(opt, 'parse_soln') && opt.parse_soln
                om.parse_soln(true);
            end
            om.soln.output.et = toc(t0);    %% stop timer
        end
        function TorF = has_parsed_soln(om)
            % Return true if model has a parsed solution.

            TorF = om.var.has_parsed_soln();
        end

        function ps = parse_soln(om, stash)
            % parse_soln - Parse solution vector and shadow prices by all named sets.
            % ::
            %
            %   PS = OM.PARSE_SOLN()
            %   OM.PARSE_SOLN(STASH)
            %
            %   For a solved model, OM.PARSE_SOLN() returns a struct of parsed
            %   solution vector and shadow price values for each named set of
            %   variables and constraints. The returned PS (parsed solution) struct
            %   has the following format, where each of the terminal elements is a
            %   struct with fields corresponding to the respective named sets:
            %
            %   Output:
            %       PS
            %           .var
            %               .val
            %               .mu_l
            %               .mu_u
            %           .lin
            %               .mu_l
            %               .mu_u
            %           .qcn
            %               .mu_l
            %               .mu_u
            %           .nle
            %               .lam
            %           .nli
            %               .mu
            %
            %   The value of each element in the returned struct can be obtained
            %   via the GET_SOLN method as well, but using PARSE_SOLN is generally
            %   more efficient if a complete set of values is needed.
            %
            %   If the optional STASH input argument is present and true, the fields
            %   of the return struct are copied to OM.SOLN.
            %
            % See also get_soln.
            
            if nargin < 2
                stash = false;
            end
            
            if ~om.is_solved()
                error('mp.opt_model.parse_soln: model not solved');
            end
            
            %% var
            ps = struct('var', om.var.parse_soln(om.soln, stash));
            
            %% lin
            ps_lin = om.lin.parse_soln(om.soln, stash);
            if ~isempty(ps_lin)
                ps.lin = ps_lin;
            end
            
            %% qcn
            ps_qcn = om.qcn.parse_soln(om.soln, stash);
            if ~isempty(ps_qcn)
                ps.qcn = ps_qcn;
            end
            
            %% nle
            ps_nle = om.nle.parse_soln(om.soln, true, stash);
            if ~isempty(ps_nle)
                ps.nle = ps_nle;
            end
            
            %% nli
            ps_nli = om.nli.parse_soln(om.soln, false, stash);
            if ~isempty(ps_nli)
                ps.nli = ps_nli;
            end
            
            %%-----  DEPRECATED  -----
            %% if requested, stash the result directly in om.soln
            %% (they are already stashed in the soln property of each set type)
            % if stash
            %     om.soln = nested_struct_copy(om.soln, ps, struct('copy_mode', '='));
            % end
        end

        function display(om, varargin)
            % display - Displays the object.
            %
            % Called when semicolon is omitted at the command-line. Displays the details
            % of the variables, constraints, costs included in the model.
            %
            % See also mp.opt_model.

            if nargin < 2
                more_set_types = {};
            else
                more_set_types = varargin{1};
            end
            
            %% display details of each set type
            set_types = om.get_set_types();
            set_types = horzcat(set_types, more_set_types);
            % set_types = {'var', 'nle', 'nli', 'lin', 'qcn', 'qdc', 'nlc', more_set_types{:}};
            fprintf('\n');
            for k = 1:length(set_types)
                om.(set_types{k}).display(set_types{k});
            end
            
            %% user data
            fields = fieldnames(om.userdata);
            if ~isempty(fields)
                fprintf('\nUSER DATA\n')
                fprintf('=========\n')
                fprintf('  name                               size       class\n');
                fprintf(' ------------------------------   -----------  --------------------\n');
                for k = 1:length(fields)
                    f = om.userdata.(fields{k});
                    [m, n] = size(f);
                    fprintf('  %-31s %5dx%-5d   %s\n', fields{k}, m, n, class(f));
                end
            else
                fprintf('USER DATA                   :  <none>\n');
            end
        end

        function om = display_soln(om, varargin)
            % display_soln - Display solution values.
            % ::
            %
            %   OM.DISPLAY_SOLN()
            %   OM.DISPLAY_SOLN(SET_TYPE)
            %   OM.DISPLAY_SOLN(SET_TYPE, NAME)
            %   OM.DISPLAY_SOLN(SET_TYPE, NAME, IDX)
            %   OM.DISPLAY_SOLN(FID)
            %   OM.DISPLAY_SOLN(FID, SET_TYPE)
            %   OM.DISPLAY_SOLN(FID, SET_TYPE, NAME)
            %   OM.DISPLAY_SOLN(FID, SET_TYPE, NAME, IDX)
            %
            %   Displays the model solution, including values, bounds and shadow
            %   prices for variables and linear constraints, values and shadow
            %   prices for nonlinear constraints, and individual cost components.
            %
            %   Results are displayed for each SET_TYPE or specified SET_TYPE and
            %   for each named/indexed set or a specified NAME/IDX.
            %
            %   Inputs:
            %       SET_TYPE - one of the following, specifying the type of set:
            %           'var' - variables
            %           'lin' - linear constraints
            %           'nle' - nonlinear equality constraints
            %           'nli' - nonlinear inequality constraints
            %           'nlc' - nonlinear costs
            %           'qdc' - quadratic costs
            %         or
            %           a cell array of one or more of the above
            %         or
            %           '' or 'all' - indicating to display all
            %       NAME - (optional) char array specifying the name of the set
            %       IDX  - (optional) cell array specifying the indices of the set
            %
            %   Examples:
            %       om.display_soln('var');
            %       om.display_soln({'nle', 'nli'});
            %       om.display_soln('var', 'P');
            %       om.display_soln('lin', 'lin_con_1');
            %       om.display_soln('nle', 'nle_con_b', {2,3});
            %
            % See also get_soln, parse_soln.

            %% input arg handling
            if nargin < 2 || ischar(varargin{1})
                fid = 1;
                args = varargin;
            else
                fid = varargin{1};
                args = varargin(2:end);
            end
            nargs = length(args);
            
            set_type = 'all';
            name = [];
            idx = [];
            if nargs >= 1
                set_type = args{1};
                if nargs >= 2
                    name = args{2};
                    if nargs >= 3
                        idx = args{3};
                    end
                end
            end
            
            %% print header
            if om.is_solved()
                if strcmp(set_type, 'all')
                    set_types = om.get_set_types(); %% all set types
                elseif ~iscell(set_type)
                    set_types = {set_type}; %% make set_type a cell array of char arrays
                else
                    set_types = set_type;
                end
            
                for ss = 1:length(set_types)
                    st = set_types{ss};
                    om_st = om.(st);
            
                    switch st
                    case 'var'
                        om_st.display_soln(om.soln, fid, args{2:end});
                    case 'nle'
                        om_st.display_soln(om.var, om.soln, 1, fid, args{2:end});
                    case 'nli'
                        om_st.display_soln(om.var, om.soln, 0, fid, args{2:end});
                    otherwise
                        om_st.display_soln(om.var, om.soln, fid, args{2:end});
                    end
                end             %% loop over set types
            else
                fprintf(fid, 'Not a solved model.\n');
            end
        end
    end     %% methods
end         %% classdef

function d = copy_prop(s, d, prop)
    if isa(s.(prop), 'mp.sm_quad_cost_legacy')
        d.(prop) = s.(prop).copy('mp.sm_quad_cost');
    elseif isa(s.(prop), 'mp.set_manager')
        d.(prop) = s.(prop).copy();
    elseif isa(d.(prop), 'mp.set_manager')
        d.(prop) = nested_struct_copy( ...
            d.(prop), s.(prop));
    else
        d.(prop) = s.(prop);
    end
end

%% system of nonlinear and linear equations
function [f, J] = nleq_fcn_(om, x)
    flin = []; Jlin = [];
    fqcn = []; Jqcn = [];
    fnln = []; Jnln = [];
    if om.lin.get_N()
        [flin, ~, Jlin] = om.lin.eval(om.var, x);
    end
    if nargout > 1
        if om.qcn.get_N()
            [fqcn, Jqcn] = om.qcn.eval(om.var, x);
        end
        if om.nle.get_N()
            [fnln, Jnln] = om.nle.eval(om.var, x);
        end
        J = [Jnln; Jqcn; Jlin];
    else
        if om.qcn.get_N()
            fqcn = om.qcn.eval(om.var, x);
        end
        if om.nle.get_N()
            fnln = om.nle.eval(om.var, x);
        end
    end
    f = [fnln; fqcn; flin];
end
