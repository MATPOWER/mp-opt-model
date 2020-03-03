function obj = init_indexed_name(obj, set_type, name, dim_list)
%INIT_INDEXED_NAME  Initializes the dimensions for an indexed named set.
%
%   OBJ.INIT_INDEXED_NAME(SET_TYPE, NAME, DIM_LIST)
%
%   Initializes the dimensions for an indexed named variable, constraint
%   or cost set.
%
%   Variables, constraints and costs are referenced in OPT_MODEL in terms
%   of named sets. The specific type of named set being referenced is
%   given by SET_TYPE, with the following valid options:
%       SET_TYPE = 'var'   => variable set
%       SET_TYPE = 'lin'   => linear constraint set
%       SET_TYPE = 'nle'   => nonlinear equality constraint set
%       SET_TYPE = 'nli'   => nonlinear inequality constraint set
%       SET_TYPE = 'cost'  => cost set
%
%   Indexed Named Sets
%
%   A variable, constraint or cost set can be identified by a single NAME,
%   such as 'Pmismatch', or by a name that is indexed by one or more indices,
%   such as 'Pmismatch(3,4)'. For an indexed named set, before adding the
%   indexed variable, constraint or cost sets themselves, the dimensions of
%   the indexed set must be set by calling INIT_INDEXED_NAME, where
%   DIM_LIST is a cell array of the dimensions.
%
%   Examples:
%       %% linear constraints with indexed named set 'R(i,j)'
%       obj.init_indexed_name('lin', 'R', {2, 3});
%       for i = 1:2
%         for j = 1:3
%           obj.add_lin_constraint('R', {i, j}, A{i,j}, ...);
%         end
%       end
%
%   See also OPT_MODEL, ADD_VAR, ADD_LIN_CONSTRAINT, ADD_NLN_CONSTRAINT,
%            ADD_QUAD_COST, ADD_NLN_COST and ADD_LEGACY_COST.

%   MATPOWER
%   Copyright (c) 2008-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% check for valid type for named set
st_label = obj.valid_named_set_type(set_type);
if st_label
    ff = set_type;
else
    error('@opt_model/init_indexed_name: ''%s'' is not a valid SET_TYPE, must be one of ''var'', ''lin'', ''nle'', ''nli'', ''cost''', set_type);
end

%% prevent duplicate name in set of specified type
if isfield(obj.(ff).idx.N, name)
    error('@opt_model/init_indexed_name: %s set named ''%s'' already exists', ...
        st_label, name);
end

%% use column vector if single dimension
if length(dim_list) == 1
    dim_list = {dim_list{:}, 1};
end

%% add general info about this named set
zero_vector = zeros(dim_list{:});
obj.(ff).idx.i1.(name)  = zero_vector;  %% starting index
obj.(ff).idx.iN.(name)  = zero_vector;  %% ending index
obj.(ff).idx.N.(name)   = zero_vector;  %% number of vars/constraints/costs
