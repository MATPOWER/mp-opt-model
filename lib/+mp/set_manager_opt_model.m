classdef set_manager_opt_model < mp.set_manager
% mp.set_manager_opt_model -  MP Set Manager base class for opt_model fields.
% ::
%
%   sm = mp.set_manager_opt_model(label)
%
% Implements functionality to handle parameters and solution data for
% set types used to implement properties of the opt_model class.
%
% mp.set_manager_opt_model Methods:
%   * params - *(abstract)* return set-type-specific parameter data
%   * set_params - *(abstract)* modify set-type-specific parameter data
%
% By convention, ``sm`` is the variable name used for mp.set_manager_opt_model
% objects.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

    methods
        function obj = set_manager_opt_model(varargin)
            % Constructor.
            % ::
            %
            %   sm = mp.set_manager_opt_model(label)

            obj@mp.set_manager(varargin{:});
        end

        function rv = params(obj, name, idx)
            % Return set-type-specific parameters.
            % ::
            %
            %   [...] = sm.params()
            %   [...] = sm.params(name)
            %   [...] = sm.params(name, idx_list)
            %
            % .. note:: This abstract method must be implemented by a
            %   subclass.
            %
            % Returns set-type-specific parameters for the full set, if called
            % without input arguments, or for a specific named or named and
            % indexed subset.
            %
            % Inputs:
            %   name (char array) : *(optional)* name of subset
            %   idx_list (cell array) : *(optional)* index list for subset
            %
            % Outputs are determined by the implementing subclass.
            %
            % See also mp.set_manager.add, set_params.
        end

        function obj = set_params(obj, name, idx, params, vals)
            % Modify parameter data.
            % ::
            %
            %   sm.set_params(name, params, vals)
            %   sm.set_params(name, idx, params, vals)
            %
            % .. note:: This abstract method must be implemented by a
            %   subclass.
            %
            % This method can be used to modify set-type-specific parameters
            % for an existing subset.
            %
            % Inputs:
            %   name (char array) : name of subset/block of entities to modify
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of entities modify (for an indexed subset)
            %   params : can be one of three options:
            %
            %       - ``'all'`` - indicates that ``vals`` is a cell array
            %         whose elements correspond to the input parameters of
            %         the :meth:`add() <mp.set_manager.add>` method
            %       - name of a parameter - ``val`` is the value of that
            %         parameter
            %       - cell array of parameter names - ``vals`` is a cell array
            %         of corresponding values
            %   vals : new value or cell array of new values corresponding to
            %       ``params``
            %
            % Valid parameter names are defined by the implementing subclass.
            %
            % See also mp.set_manager.add, params.
        end
    end     %% methods

    methods (Access=protected)
        function [is_all, np, params, vals] = set_params_std_args(obj, default_params, params, vals)
            % Standardize input args for use in subclass set_params method.
            % ::
            %
            %   [is_all, np, params, vals] = sm.set_params_std_args(default_params, params, vals)
            %
            % See also set_params.

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
        end

        function obj = set_params_update_dims(obj, dN, name, idx)
            % Update parameter dimensions for use in subclass set_params method.
            % ::
            %
            %   sm.set_params_update_dims(dN, name, idx)
            %
            % See also set_params.

            %% calls to substruct() are relatively expensive, so we pre-build the
            %% struct for addressing num array fields
            %% sn = substruct('.', name, '()', idx);
            sn = struct('type', {'.', '()'}, 'subs', {'', 1});  %% num array field
            update = 0;             %% not yet reached set being updated
            update_i1 = 0;          %% flag to indicate whether to update i1
            for k = 1:obj.NS
                o = obj.order(k);
                if ~update && strcmp(o.name, name) && isequal(o.idx, idx)
                    update = 1;     %% arrived at set being updated
                end
                if update
                    if isempty(o.idx)   %% simple named set
                        if update_i1
                            obj.idx.i1.(o.name) = obj.idx.i1.(o.name) + dN;
                        else
                            obj.idx.N.(o.name) = obj.idx.N.(o.name) + dN;
                        end
                        obj.idx.iN.(o.name) = obj.idx.iN.(o.name) + dN;
                    else                %% indexed named set
                        sn(1).subs = o.name;
                        sn(2).subs = o.idx;
                        if update_i1
                            v = subsref(obj.idx.i1, sn);
                            obj.idx.i1 = subsasgn(obj.idx.i1, sn, v + dN);
                        else
                            v = subsref(obj.idx.N, sn);
                            obj.idx.N = subsasgn(obj.idx.N, sn, v + dN);
                        end
                        v = subsref(obj.idx.iN, sn);
                        obj.idx.iN = subsasgn(obj.idx.iN, sn, v + dN);
                    end
                    update_i1 = 1;  %% update i1 from here on out
                end
            end
            obj.N = obj.N + dN;
        end

        function default_tags = get_soln_default_tags(obj)
            % Return default tags for use in get_soln method.
            % ::
            %
            %   default_tags = sm.get_soln_default_tags()
            %
            % .. note:: This protected abstract method must be implemented by
            %    a subclass.
            %
            % Output:
            %   default_tags (cell array) : tags defining the default outputs
            %       of get_soln()
            %
            % See also get_soln.

            default_tags = {};
        end

        function [tags, name, idx, N, i1, iN] = get_soln_std_args(obj, tags, name, idx)
            % Standardize input args for use in subclass get_soln method.
            % ::
            %
            %   [tags, name, idx, N, i1, iN] = sm.get_soln_std_args(tags, name, idx)
            %
            % See also get_soln.

            %% input arg handling
            if nargin == 2              %% obj.get_soln(soln, name)
                idx = [];
                name = tags;
                tags = {};
            elseif nargin == 3
                if ischar(name)         %% obj.get_soln(soln, tags, name)
                    idx = [];
                else                    %% obj.get_soln(soln, name, idx)
                    idx = name;
                    name = tags;
                    tags = {};
                end
            end

            %% set up tags for default outputs
            if isempty(tags)
                tags = obj.get_soln_default_tags();
            elseif ~iscell(tags)
                tags = { tags };
            end

            %% set up indexing
            if isempty(idx)         %% simple named set
                N = obj.idx.N.(name);
                i1 = obj.idx.i1.(name);         %% starting row index
                iN = obj.idx.iN.(name);         %% ending row index
            else                    %% indexed named set
                %% calls to substruct() are relatively expensive, so we pre-build the
                %% structs for addressing cell and numeric array fields, updating only
                %% the subscripts before use
                sn = struct('type', {'.', '()'}, 'subs', {name, idx});  %% num array field
                N = subsref(obj.idx.N, sn);
                i1 = subsref(obj.idx.i1, sn);   %% starting row index
                iN = subsref(obj.idx.iN, sn);   %% ending row index
            end
        end

        function ps = parse_soln_fields(obj, params)
            % Parse solution fields for subclass parse_soln method.
            % ::
            %
            %   ps = sm.parse_soln_fields(params)
            %
            % Parse solution fields.
            %
            % Input:
            %   params (struct) : struct array with fields:
            %
            %           - ``src`` - values from full solution struct
            %           - ``dst`` - name of destination field in parsed
            %             solution struct
            %
            % Output:
            %   ps (struct) : parsed solution struct, fields names match
            %       those provided in ``params.dst``, values are structs
            %       whose field names correspond to the named subsets in the
            %       set

            %% calls to substruct() are relatively expensive, so we pre-build the
            %% structs for addressing cell and numeric array fields, updating only
            %% the subscripts before use
            persistent sn;
            persistent sc;
            if isempty(sc)
                sc = struct('type', {'.', '{}'}, 'subs', {'', 1});  %% cell array field
            end
            if isempty(sn)
                sn = struct('type', {'.', '()'}, 'subs', {'', 1});  %% num array field
            end

            ps = struct();      %% parsed solution

            np = length(params);
            have_param = zeros(np, 1);
            for j = 1:np
                have_param(j) = ~isempty(params(j).src);
            end
            for k = 1:obj.NS
                name = obj.order(k).name;
                idx = obj.order(k).idx;
                if isempty(idx)
                    N = obj.idx.N.(name);
                else
                    sn(1).subs = name;
                    sn(2).subs = idx;
                    N = subsref(obj.idx.N, sn);
                    need_init = all([idx{:}] == 1);
                end
                if N
                    for j = 1:np
                        if have_param(j)    %% parameter is available
                            dname = params(j).dst;  %% destination field name
                            if isempty(idx)
                                i1 = obj.idx.i1.(name);
                                iN = obj.idx.iN.(name);
                                ps.(dname).(name)  = params(j).src(i1:iN);
                            else
                                if need_init
                                    param_names = fieldnames(obj.data);
                                    ps.(dname).(name) = cell(size(obj.data.(param_names{1}).(name)));
                                end
                                i1 = subsref(obj.idx.i1, sn);    %% starting row index
                                iN = subsref(obj.idx.iN, sn);    %% ending row index
                                sc(1).subs = name;
                                sc(2).subs = idx;
                                ps.(dname) = ...
                                    subsasgn(ps.(dname), sc, params(j).src(i1:iN));
                            end
                        end
                    end     %% for j
                end
            end     %% for k
        end
    end     %% methods (Access=protected)
end         %% classdef
