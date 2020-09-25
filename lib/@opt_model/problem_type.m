function prob = problem_type(om, recheck)
%PROBLEM_TYPE  Return a string identifying the type of mathematical program
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
%   constraints.
%
%   Outputs:
%       PROB_TYPE : problem type, one of the following strings:
%           'LEQ'   - linear equations
%           'NLEQ'  - nonlinear equations
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
%   See also OPT_MODEL

%   MP-Opt-Model
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if isempty(om.prob_type) || nargin > 1 && recheck
    nleN = om.getN('nle');      %% nonlinear equalities
    nliN = om.getN('nli');      %% nonlinear inequalities
    nlcN = om.getN('nlc');      %% general nonlinear costs
    qdcN = om.getN('qdc');      %% quadratic costs
    linN = om.getN('lin');      %% linear constraints
    varN = om.getN('var');      %% variables
    if nlcN || qdcN         %% problem has costs
        if nliN || nleN || nlcN %% nonlinear
            prob = 'NLP';           %% nonlinear program
        else                    %% linear constraints, no general nonlinear costs
            %% get quadratic cost coefficients
            H = om.params_quad_cost();
            if isempty(H) || ~any(any(H))
                prob = 'LP';        %% linear program
            else
                prob = 'QP';        %% quadratic program
            end
        end
    else                    %% problem has no costs
        if nliN
            error('@opt_model/problem_type: invalid problem - nonlinear inequality constraints with no costs');
        end
        if nleN + linN == varN  %% square system
            if linN > 0
                %% get lower & upper bounds
                [A, l, u] = om.params_lin_constraint();
                if any(l ~= u)
                    error('@opt_model/problem_type: invalid problem - linear inequality constraints with no costs');
                end
            end
            if nleN
                prob = 'NLEQ';      %% square nonlinear set of equations
            else
                prob = 'LEQ';       %% square linear set of equations
            end
        else
            error('@opt_model/problem_type: invalid problem - non-square system with no costs');
        end
    end
    if om.is_mixed_integer() && ~strcmp(prob, 'NLEQ')
        prob = ['MI' prob];
    end
    om.prob_type = prob;    %% cache it
else
    prob = om.prob_type;    %% return cached type
end
