function om = add_named_set(om, set_type, name, idx, N, varargin)
%ADD_NAMED_SET  Adds a named set of variables/constraints/costs to the model.
%
%   -----  PRIVATE METHOD  -----
%
%   This method is intended to be a private method, used internally by
%   the public methods ADD_VAR, ADD_LIN_CONSTRAINT, ADD_NLN_CONSTRAINT
%   ADD_QUAD_COST, ADD_NLN_COST and ADD_LEGACY_COST.
%
%   Variable Set
%       OM.ADD_NAMED_SET('var', NAME, IDX_LIST, N, V0, VL, VU, VT);
%
%   Linear Constraint Set
%       OM.ADD_NAMED_SET('lin', NAME, IDX_LIST, N, A, L, U, VARSETS);
%
%   Nonlinear Inequality Constraint Set
%       OM.ADD_NAMED_SET('nle', NAME, IDX_LIST, N, FCN, HESS, COMPUTED_BY, VARSETS);
%
%   Nonlinear Inequality Constraint Set
%       OM.ADD_NAMED_SET('nli', NAME, IDX_LIST, N, FCN, HESS, COMPUTED_BY, VARSETS);
%
%   Quadratic Cost Set
%       OM.ADD_NAMED_SET('qdc', NAME, IDX_LIST, N, CP, VARSETS);
%
%   General Nonlinear Cost Set
%       OM.ADD_NAMED_SET('nlc', NAME, IDX_LIST, N, FCN, VARSETS);
%
%   Legacy Cost Set
%       OM.ADD_NAMED_SET('cost', NAME, IDX_LIST, N, CP, VARSETS);
%
%   See also OPT_MODEL, ADD_VAR, ADD_LIN_CONSTRAINT, ADD_NLN_CONSTRAINT
%            ADD_QUAD_COST and ADD_NLN_COST.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% check for valid type for named set
st_label = om.valid_named_set_type(set_type);
if st_label
    ff = set_type;
    om_ff = om.(ff);
    om.(ff) = [];
else
    error('@opt_model/add_named_set: ''%s'' is not a valid SET_TYPE, must be one of ''var'', ''lin'', ''nle'', ''nli'', ''qdc'', ''nlc'', ''cost''', set_type);
end

%% add general indexing info about this named set
if isempty(idx)     %% simple named set
    %% prevent duplicate name in set of specified type
    if isfield(om_ff.idx.N, name)
        error('@opt_model/add_named_set: %s set named ''%s'' already exists', st_label, name);
    end

    %% add indexing info about this set
    om_ff.idx.i1.(name)  = om_ff.N + 1; %% starting index
    om_ff.idx.iN.(name)  = om_ff.N + N; %% ending index
    om_ff.idx.N.(name)   = N;           %% number in set
    om_ff.N  = om_ff.idx.iN.(name);     %% number of elements of this type
    om_ff.NS = om_ff.NS + 1;            %% number of sets of this type
    om_ff.order(om_ff.NS).name = name;  %% add name to ordered list of sets
    om_ff.order(om_ff.NS).idx  = {};
else                %% indexed named set
    %% calls to substruct() are relatively expensive, so we pre-build the
    %% struct for addressing numeric array fields
    %% sn = substruct('.', name, '()', idx);
    sn = struct('type', {'.', '()'}, 'subs', {name, idx});  %% num array field

    %% prevent duplicate name in set of specified type
    if subsref(om_ff.idx.i1, sn) ~= 0
        str = '%d'; for m = 2:length(idx), str = [str ',%d']; end
        nname = sprintf(['%s(' str, ')'], name, idx{:});
        error('@opt_model/add_named_set: %s set named ''%s'' already exists', st_label, nname);
    end

    %% add indexing info about this set
    om_ff.idx.i1  = subsasgn(om_ff.idx.i1, sn, om_ff.N + 1);    %% starting index
    om_ff.idx.iN  = subsasgn(om_ff.idx.iN, sn, om_ff.N + N);    %% ending index
    om_ff.idx.N   = subsasgn(om_ff.idx.N,  sn, N);              %% number in set
    om_ff.N  = subsref(om_ff.idx.iN, sn);   %% number of elements of this type
    om_ff.NS = om_ff.NS + 1;                %% number of sets of this type
    om_ff.order(om_ff.NS).name = name;      %% add name to ordered list of sets
    om_ff.order(om_ff.NS).idx  = idx;       %% add indices to ordered list of sets
end

om.(ff) = om_ff;
