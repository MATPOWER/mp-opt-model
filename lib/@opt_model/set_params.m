function om = set_params(om, st, name, idx, params, vals)
% set_params - Modifies parameters for variable, cost or constraint in model
% ::
%
%   This method can be used to modify parameters for an existing variable,
%   constraint or cost in the model.
%
%   OM.SET_PARAMS(SET_TYPE, NAME, PARAMS, VALS)
%   OM.SET_PARAMS(SET_TYPE, NAME, IDX, PARAMS, VALS)
%
%   Inputs:
%       SET_TYPE : one of 'var', 'lin', 'nle', 'nli', 'nlc', 'qdc' for
%           variables, linear constraints, nonlinear equality constraints,
%           nonlinear inequality constraints, general nonlinear costs,
%           and quadratic costs, respectively
%       NAME : name of set
%       IDX : index of named set (for an indexed set)
%       PARAMS : can be one of three options:
%           1 - 'all', indicating that VALS is a cell array whose elements
%               correspond to the input parameters of the respective
%               add_*() method
%           2 - the name of a PARAM, VAL is the value of that parameter
%           3 - a cell array of PARAM names, VALS is a cell array of
%               corresponding values
%           Note: Changing the dimension of a 'var' is not allowed
%               and changing the #1 ('all') is the only option for 'nle', 'nli', and 'nlc'
%       VALS : new value or cell array of new values for PARAMS
%
%   Valid PARAM names:
%       var - N, v0, vl, vu, vt
%       lin - A, l, u, vs
%       nle - N, fcn, hess, include, vs
%       nli - N, fcn, hess, include, vs
%       nlc - N, fcn, vs
%       qdc - Q, c, k, vs
%
%   Examples:
%       om.set_params('var', 'Pg', 'v0', Pg0);
%       om.set_params('lin', 'y', {2,3}, {'l', 'u'}, {l, u});
%       om.set_params('nle', 'Pmis', 'all', {N, @fcn, @hess, vs});
%
% See also opt_model, add_var, add_lin_constraint, add_nln_constraint,
% add_quad_cost, add_nln_cost.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin < 6
    vals = params;
    params = idx;
    idx = {};
end

%% create default list of parameters to update based on set type & inputs
switch st
    case 'var'
        om.var.set_params(name, idx, params, vals);
        return;
    case 'lin'
        om.lin.set_params(om.var, name, idx, params, vals);
        return;
    case {'nle', 'nli'}
        default_params = {'N', 'fcn', 'hess', 'vs'};
    case 'nlc'
        default_params = {'N', 'fcn', 'vs'};
    case 'qdc'
        om.qdc.set_params(om.var, name, idx, params, vals);
        return;
    otherwise
        error('opt_model.set_params: ''%s'' is not a valid SET_TYPE', st);
end

%% standardize provided arguments in cell arrays params, vals
is_all = 0;     %% flag to indicate all params for set are being replaced
if ischar(params)
    if strcmp(params, 'all')
        is_all = 1;
        np = length(vals);      %% number of parameter values provided
        params = default_params(1:np);
    else
        np = 1;                 %% number of parameter values provided
        params = {params};
        vals = {vals};
    end
else
    np = length(vals);          %% number of parameter values provided
end

if ~isempty(idx)
    %% calls to substruct() are relatively expensive, so we pre-build the
    %% struct for addressing cell array fields
    %% sc = substruct('.', name, '{}', idx);
    sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
end

switch st
    case 'var'
    case 'lin'
    case 'qdc'
    case {'nle', 'nli'}
        %% get current parameters
        if isempty(idx)
            [N0, fcn, hess, vs, include] = om.params_nln_constraint(st(3) == 'e', name, idx);
        else
            [N0, fcn, hess, vs] = om.params_nln_constraint(st(3) == 'e', name, idx);
            include = '';
        end
        if isempty(vs), vs = {vs}; end
        p = struct('N', N0, 'fcn', fcn, 'hess', hess, 'vs', vs);    %% current parameters
        u = struct('N',  0, 'fcn',   0, 'hess',    0, 'vs',  0);    %% which ones to update

        %% replace with new parameters
        for k = 1:np
            p.(params{k}) = vals{k};
            u.(params{k}) = 1;
        end
        N = p.N;

        %% set missing default params for 'all'
        if is_all
            u.N    = 1;         %% always update N
            u.fcn  = 1;         %% alwaus update fcn
            u.hess = 1;         %% alwaus update hess
            if np < 4
                p.vs = {};
                u.vs = 1;       %% update vs
            end
        end

        %% check consistency of parameters
        %% no dimension change unless 'all'
        if N ~= N0 && ~is_all
            error('opt_model.set_params: dimension change for ''%s'' ''%s'' not allowed except for ''all''', st, nameidxstr(name, idx));
        end

        %% included constraints not yet implemented
        if ~isempty(include)
            error('opt_model.set_params: modifications for ''%s'' ''%s'' not (yet) supported since it includes evaluation of other constraints', st, nameidxstr(name, idx));
        end

        %% convert vs to struct
        if u.vs
            p.vs = om.varsets_cell2struct(p.vs);
        end

        %% assign new parameters
        if isempty(idx)     %% simple named set
            for k = 2:length(default_params)
                pn = default_params{k};     %% param name
                if u.(pn)   %% assign new val for this parameter
                    om.(st).data.(pn).(name) = p.(pn);
                end
            end
        else                %% indexed named set
            for k = 2:length(default_params)
                pn = default_params{k};     %% param name
                if u.(pn)   %% assign new val for this parameter
                    om.(st).data.(pn) = subsasgn(om.(st).data.(pn), sc, p.(pn));
                end
            end
        end
    case 'nlc'
        %% get current parameters
        [N0, fcn, vs] = om.params_nln_cost(name, idx);
        if isempty(vs), vs = {vs}; end
        p = struct('N', N0, 'fcn', fcn, 'vs', vs);  %% current parameters
        u = struct('N',  0, 'fcn',   0, 'vs',  0);  %% which ones to update

        %% replace with new parameters
        for k = 1:np
            p.(params{k}) = vals{k};
            u.(params{k}) = 1;
        end
        N = p.N;

        %% set missing default params for 'all'
        if is_all
            u.N   = 1;          %% always update N
            u.fcn = 1;          %% alwaus update fcn
            if np < 3
                p.vs = {};
                u.vs = 1;       %% update vs
            end
        end

        %% check consistency of parameters
        %% no dimension change unless 'all'
        if N ~= N0 && ~is_all
            error('opt_model.set_params: dimension change for ''%s'' ''%s'' not allowed except for ''all''', st, nameidxstr(name, idx));
        end

        %% vector valued costs not yet implemented
        if N ~= 1
            error('opt_model.set_params: vector value for ''%s'' ''%s'' not yet implemented', st, nameidxstr(name, idx));
        end

        %% assign new parameters
        if isempty(idx)     %% simple named set
            for k = 2:length(default_params)
                pn = default_params{k};     %% param name
                if u.(pn)   %% assign new val for this parameter
                    om.nlc.data.(pn).(name) = p.(pn);
                end
            end
        else                %% indexed named set
            for k = 2:length(default_params)
                pn = default_params{k};     %% param name
                if u.(pn)   %% assign new val for this parameter
                    om.nlc.data.(pn) = subsasgn(om.nlc.data.(pn), sc, p.(pn));
                end
            end
        end
    otherwise
        error('opt_model.set_params: ''%s'' is not a valid SET_TYPE', st);
end

%% update dimensions and indexing, if necessary
dN = N - N0;
if is_all && dN
    %% calls to substruct() are relatively expensive, so we pre-build the
    %% struct for addressing num array fields
    %% sn = substruct('.', name, '()', idx);
    sn = struct('type', {'.', '()'}, 'subs', {'', 1});  %% num array field
    om_ff = om.(st);
    update = 0;             %% not yet reached set being updated
    update_i1 = 0;          %% flag to indicate whether to update i1
    for k = 1:om_ff.NS
        o = om_ff.order(k);
        if ~update && strcmp(o.name, name) && isequal(o.idx, idx)
            update = 1;     %% arrived at set being updated
        end
        if update
            if isempty(o.idx)   %% simple named set
                if update_i1
                    om_ff.idx.i1.(o.name) = om_ff.idx.i1.(o.name) + dN;
                else
                    om_ff.idx.N.(o.name) = om_ff.idx.N.(o.name) + dN;
                end
                om_ff.idx.iN.(o.name) = om_ff.idx.iN.(o.name) + dN;
            else                %% indexed named set
                sn(1).subs = o.name;
                sn(2).subs = o.idx;
                if update_i1
                    v = subsref(om_ff.idx.i1, sn);
                    om_ff.idx.i1 = subsasgn(om_ff.idx.i1, sn, v + dN);
                else
                    v = subsref(om_ff.idx.N, sn);
                    om_ff.idx.N = subsasgn(om_ff.idx.N, sn, v + dN);
                end
                v = subsref(om_ff.idx.iN, sn);
                om_ff.idx.iN = subsasgn(om_ff.idx.iN, sn, v + dN);
            end
            update_i1 = 1;  %% update i1 from here on out
        end
    end
    om_ff.N = om_ff.N + dN;
    om.(st) = om_ff;
end


function str = nameidxstr(name, idx)
str = sprintf('%s%s', name, idxstr(idx));

function str = idxstr(idx)
if isempty(idx)
    str = '';
elseif length(idx) == 1
    str = sprintf('(%d)', idx{1});
else
    str = ['(' sprintf('%d', idx{1}) sprintf(',%d', idx{2:end}) ')'];
end
