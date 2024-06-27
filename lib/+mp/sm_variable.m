classdef sm_variable < mp.set_manager
% mp.sm_variable -  MP Set Manager class for variables.
% ::
%
%   var = mp.sm_variable()
%   var = mp.sm_variable(label)
%
% MP Set Manager class for variables. Manages variable initial values,
% lower and upper bounds, and variable type, along with indexing.
%
% By convention, ``var`` is the variable name used for mp.sm_variable objects.
%
% mp.sm_variable Properties:
%   * cache - struct for caching aggregated parameters for the set
%
% mp.sm_variable Methods:
%   * sm_variable - constructor
%   * add - add a subset of variables, with initial value, bounds, and var type
%   * params - return initial values, lower bounds, upper bounds, and var type
%   * set_params - modify parameter data
%   * varsets_len - return the total number of variables in ``varsets``
%   * varsets_cell2struct - convert ``varsets`` from cell array to struct array
%
% See also mp.set_manager.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

    properties
        % struct for caching aggregated parameters for variables
%         cache = [];
    end     %% properties

    methods
        function obj = sm_variable(varargin)
            % Constructor.
            % ::
            %
            %   var = mp.sm_variable(label)

            es = struct();  %% empty struct
            obj@mp.set_manager(varargin{:});
            obj.data = struct( ...
                'v0', es, ...
                'vl', es, ...
                'vu', es, ...
                'vt', es );
        end

        function obj = add(obj, name, idx, varargin)
            % Add a subset of variables with initial value, bounds, type.
            % ::
            %
            %   var.add(name, N, v0, vl, vu, vt)
            %   var.add(name, N, v0, vl, vu)
            %   var.add(name, N, v0, vl)
            %   var.add(name, N, v0)
            %   var.add(name, N)
            %
            %   var.add(name, idx_list, N, v0, vl, vu, vt)
            %   var.add(name, idx_list, N, v0, vl, vu)
            %   var.add(name, idx_list, N, v0, vl)
            %   var.add(name, idx_list, N, v0)
            %   var.add(name, idx_list, N)
            %
            % Add a named, and possibly indexed, subset of variables to
            % the set.
            %
            % Inputs:
            %   name (char array) : name of subset/block of variables to add
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of variables to add (for an indexed subset)
            %   N (integer) : number of variables in the subset
            %   v0 (double) : *(optional, default = 0)* scalar or
            %       :math:`N \times 1` vector of variable initial values
            %   vl (double) : *(optional, default = -Inf)* scalar or
            %       :math:`N \times 1` vector of variable lower bounds
            %   vu (double) : *(optional, default = Inf)* scalar or
            %       :math:`N \times 1` vector of variable upper bounds
            %   vt (char array) : *(optional, default =* ``'C'`` *)* scalar or
            %       :math:`1 \times N` char array of variable types, where
            %       accepted types are:
            %
            %       - ``'C'`` - continuous
            %       - ``'I'`` - integer
            %       - ``'B'`` - binary
            %
            % Examples::
            %
            %   var.add('V', nb, V0, Vmin, Vmax, 'C');
            %
            %   var.init_indexed_name('x', {2, 3});
            %   for i = 1:2
            %       for j = 1:3
            %           var.add('x', {i, j}, nx(i,j), ...);
            %       end
            %   end
            %
            % See also params, set_params.

            %% call parent to handle standard indexing
            add@mp.set_manager(obj, name, idx, varargin{:});

            %% set up default args
            if iscell(idx)          %% indexed named set
                N = varargin{1};
                args = varargin(2:end);
            else                    %% simple named set
                N = idx;
                idx = {};
                args = varargin;
            end
            nargs = length(args);

            %% prepare data
            v0 = []; vl = []; vu = []; vt = [];
            if nargs >= 1
                v0 = args{1};
                if nargs >= 2
                    vl = args{2};
                    if nargs >= 3
                        vu = args{3};
                        if nargs >= 4
                            vt = args{4};
                        end
                    end
                end
            end
            if isempty(v0)
                v0 = zeros(N, 1);   %% init to zero by default
            elseif N > 1 && length(v0) == 1     %% expand from scalar as needed
                v0 = v0 * ones(N, 1);
            end
            if isempty(vl)
                vl = -Inf(N, 1);    %% unbounded below by default
            elseif N > 1 && length(vl) == 1     %% expand from scalar as needed
                vl = vl * ones(N, 1);
            end
            if isempty(vu)
                vu = Inf(N, 1);     %% unbounded above by default
            elseif N > 1 && length(vu) == 1     %% expand from scalar as needed
                vu = vu * ones(N, 1);
            end
            if isempty(vt) && N > 0
                vt = 'C';           %% all continuous by default
            end

            %% assign data
            if isempty(idx)
                obj.data.v0.(name) = v0;        %% initial value
                obj.data.vl.(name) = vl;        %% lower bound
                obj.data.vu.(name) = vu;        %% upper bound
                obj.data.vt.(name) = vt;        %% variable type
            else
                %% calls to substruct() are relatively expensive, so we
                %% pre-build the struct for addressing cell array fields
                %% sc = substruct('.', name, '{}', idx);
                sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
                obj.data.v0 = subsasgn(obj.data.v0, sc, v0);    %% initial value
                obj.data.vl = subsasgn(obj.data.vl, sc, vl);    %% lower bound
                obj.data.vu = subsasgn(obj.data.vu, sc, vu);    %% upper bound
                obj.data.vt = subsasgn(obj.data.vt, sc, vt);    %% variable type
            end
            if ~isempty(obj.cache)  %% clear cache of aggregated params
                obj.cache = [];
            end
        end

        function [v0, vl, vu, vt] = params(obj, name, idx)
            % Return initial values, lower bounds, upper bounds and var type.
            % ::
            %
            %   [v0, vl, vu] = var.params()
            %   [v0, vl, vu] = var.params(name)
            %   [v0, vl, vu] = var.params(name, idx_list)
            %   [v0, vl, vu, vt] = var.params(...)
            %
            % Returns the initial value, lower bound, upper bound, and
            % variable type for the full set of variables, if called without
            % input arguments, or for a specific named or named and indexed
            % subset. Values for the full set are cached for subsequent calls.
            %
            % Inputs:
            %   name (char array) : *(optional)* name of subset
            %   idx_list (cell array) : *(optional)* index list for subset
            %
            % Outputs:
            %   v0 (double) : column vector of variable initial values
            %   vl (double) : column vector of variable lower bounds
            %   vu (double) : column vector of variable upper bounds
            %   vt (char array) : char array of variable types:
            %
            %       - ``'C'`` - continuous
            %       - ``'I'`` - integer
            %       - ``'B'`` - binary
            %
            % Examples::
            %
            %   [x0, xmin, xmax] = var.params();
            %   [Pg0, Pmin, Pmax] = var.params('Pg');
            %   [zij0, zijmin, zijmax, ztype] = var.params('z', {i, j});
            %
            % See also add, set_params.

            if nargout > 3
                have_vt = 1;
            else
                have_vt = 0;
            end
            if nargin < 2       %% aggregate
                if isempty(obj.cache)       %% build the aggregate
                    v0 = []; vl = []; vu = []; vt = char([]);
                    %% calls to substruct() are relatively expensive, so we pre-build the
                    %% structs for addressing cell and numeric array fields, updating only
                    %% the subscripts before use
                    sc = struct('type', {'.', '{}'}, 'subs', {'', 1});  %% cell array field
                    for k = 1:obj.NS
                        name = obj.order(k).name;
                        idx = obj.order(k).idx;
                        if isempty(idx)
                            v0 = [ v0; obj.data.v0.(name) ];
                            vl = [ vl; obj.data.vl.(name) ];
                            vu = [ vu; obj.data.vu.(name) ];
                            if have_vt
                                N = obj.idx.N.(name);
                                vt0 = obj.data.vt.(name);
                                if isscalar(vt0) && N > 1 
                                    vt = [ vt char(vt0 * ones(1, N)) ];
                                else
                                    vt = [ vt vt0 ];
                                end
                            end
                        else
                            % (calls to substruct() are relatively expensive ...
                            % sc = substruct('.', name, '{}', idx);
                            % ... so replace it with these more efficient lines)
                            sc(1).subs = name;
                            sc(2).subs = idx;
                            v0 = [ v0; subsref(obj.data.v0, sc) ];
                            vl = [ vl; subsref(obj.data.vl, sc) ];
                            vu = [ vu; subsref(obj.data.vu, sc) ];
                            if have_vt
                                % (calls to substruct() are relatively expensive ...
                                % sn = substruct('.', name, '()', idx);
                                % ... so replace it with these more efficient lines)
                                sn = sc; sn(2).type = '()';
                                N = subsref(obj.idx.N, sn);
                                vt0 = subsref(obj.data.vt, sc);
                                if isscalar(vt0) && N > 1 
                                    vt = [ vt char(vt0 * ones(1, N)) ];
                                else
                                    if ~isempty(vt0)
                                        vt = [ vt vt0 ];
                                    end
                                end
                            end
                        end
                    end

                    %% cache aggregated parameters
                    obj.cache = struct('v0', v0, 'vl', vl, 'vu', vu, 'vt', vt);
                else                    %% return cached values
                    v0 = obj.cache.v0;
                    vl = obj.cache.vl;
                    vu = obj.cache.vu;
                    vt = obj.cache.vt;
                end
            else                %% individual set
                if isfield(obj.idx.N, name)
                    if nargin < 3 || isempty(idx)
                        v0 = obj.data.v0.(name);
                        vl = obj.data.vl.(name);
                        vu = obj.data.vu.(name);
                        if have_vt
                            N = obj.idx.N.(name);
                            vt0 = obj.data.vt.(name);
                            if isscalar(vt0) && N > 1 
                                vt = char(vt0 * ones(1, N));
                            else
                                vt = vt0;
                            end
                        end
                    else
                        % (calls to substruct() are relatively expensive ...
                        % sc = substruct('.', name, '{}', idx);
                        % ... so replace it with these more efficient lines)
                        sc = struct('type', {'.', '{}'}, 'subs', {name, idx});
                        v0 = subsref(obj.data.v0, sc);
                        vl = subsref(obj.data.vl, sc);
                        vu = subsref(obj.data.vu, sc);
                        if have_vt
                            % (calls to substruct() are relatively expensive ...
                            % sn = substruct('.', name, '()', idx);
                            % ... so replace it with these more efficient lines)
                            sn = sc; sn(2).type = '()';
                            N = subsref(obj.idx.N, sn);
                            vt0 = subsref(obj.data.vt, sc);
                            if isscalar(vt0) && N > 1 
                                vt = char(vt0 * ones(1, N));
                            else
                                vt = vt0;
                            end
                        end
                    end
                else
                    v0 = [];
                    vl = [];
                    vu = [];
                    if have_vt
                        vt = [];
                    end
                end
            end
        end

        function obj = set_params(obj, name, idx, params, vals)
            % set_params - Modify parameter data.
            % ::
            %
            %   var.set_params(name, params, vals)
            %   var.set_params(name, idx, params, vals)
            %
            % This method can be used to modify parameters for an existing
            % subset of variables.
            %
            % Inputs:
            %   name (char array) : name of subset/block of variables to modify
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of variables modify (for an indexed subset)
            %   params : can be one of three options:
            %
            %       - ``'all'`` - indicates that ``vals`` is a cell array
            %         whose elements correspond to the input parameters of
            %         the add() method
            %       - name of a parameter - ``val`` is the value of that
            %         parameter
            %       - cell array of parameter names - ``vals`` is a cell array
            %         of corresponding values
            %   vals : new value or cell array of new values corresponding to
            %       ``params``
            %
            % Valid parameter names are ``N``, ``v0``, ``vl``, ``vu``, ``vt``.
            %
            % .. note::
            %    Changing the dimension of a variable subset is not allowed.
            %
            % Examples::
            %
            %   var.set_params('Pg', 'v0', Pg0);
            %   var.set_params('y', {2,3}, {'vl', 'vu'}, {yl, yu});
            %   var.set_params('Va', 'all', {N, va0, val, vau, vat});
            %
            % See also add, params.

            if nargin < 5
                vals = params;
                params = idx;
                idx = {};
            end

            %% create default list of parameters to update based on set type & inputs
            default_params = {'N', 'v0', 'vl', 'vu', 'vt'};

            %% standardize provided arguments in cell arrays params, vals
            [is_all, np, params, vals] = ...
                obj.set_params_std_args(default_params, params, vals);

            %% get current parameters
            [v0, vl, vu, vt] = obj.params(name, idx);
            N0 = obj.get_N(name, idx);
            p = struct('N', N0, 'v0', v0, 'vl', vl, 'vu', vu, 'vt', vt);    %% current parameters
            u = struct('N',  0, 'v0',  0, 'vl',  0, 'vu',  0, 'vt',  0);    %% which ones to update

            %% replace with new parameters
            for k = 1:np
                p.(params{k}) = vals{k};
                u.(params{k}) = 1;
            end
            N = p.N;

            %% set missing default params for 'all'
            if is_all
                if np < 5
                    p.vt = 'C';
                    u.vt = 1;               %% update vt
                    if np < 4
                        p.vu = Inf(N, 1);
                        u.vu = 1;           %% update vu
                        if np < 3
                            p.vl = -Inf(N, 1);
                            u.vl = 1;       %% update vl
                            if np < 2
                                p.v0 = zeros(N, 1);
                                u.v0 = 1;   %% update v0
                            end
                        end
                    end
                end
            end

            %% check consistency of parameters
            %% no dimension change
            if N ~= N0
                error('sm_variable.set_params: dimension change for ''%s'' not allowed', obj.nameidxstr(name, idx));
            end

            %% check sizes of new values of v0, vl, vu, vt
            for pn = {'v0', 'vl', 'vu', 'vt'}
                if u.(pn{1})
                    nn = length(p.(pn{1}));
                    if nn ~= N
                        if nn == 0
                            switch pn{1}
                                case 'v0'
                                    p.(pn{1}) = zeros(N, 0);
                                case 'vl'
                                    p.(pn{1}) = -Inf(N, 0);
                                case 'vu'
                                    p.(pn{1}) =  Inf(N, 0);
                                case 'vt'
                                    p.(pn{1}) = 'C';
                            end
                        elseif nn == 1
                            if pn{1} ~= 'vt'
                                p.(pn{1}) = p.(pn{1}) * ones(N, 1);   %% expand from scalar
                            end
                        else
                            error('sm_variable.set_params: parameter ''%s'' ''%s'' should have length %d (or 1)', obj.nameidxstr(name, idx), pn{1}, N);
                        end
                    end
                end
            end

            %% assign new parameters
            if isempty(idx)     %% simple named set
                for k = 2:length(default_params)
                    pn = default_params{k};     %% param name
                    if u.(pn)   %% assign new val for this parameter
                        obj.data.(pn).(name) = p.(pn);
                    end
                end
            else                %% indexed named set
                %% calls to substruct() are relatively expensive, so we pre-build the
                %% struct for addressing cell array fields
                %% sc = substruct('.', name, '{}', idx);
                sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
                for k = 2:length(default_params)
                    pn = default_params{k};     %% param name
                    if u.(pn)   %% assign new val for this parameter
                        obj.data.(pn) = subsasgn(obj.data.(pn), sc, p.(pn));
                    end
                end
            end

            %% clear cached parameters
            obj.cache = [];

            %% update dimensions and indexing, if necessary
            dN = N - N0;
            if is_all && dN
                obj.set_params_update_dims(dN, name, idx);
            end
        end

        function nv = varsets_len(obj, vs)
            % varsets_len - Return the total number of variables in ``varsets``.
            % ::
            %
            %   nv = var.varsets_len(varsets)
            %
            % Return the total number of elements in the variable sub-vector
            % specified by ``varsets``.
            %
            % Input:
            %   varsets (cell or struct array) : cell array of names of
            %       variable subsets, or a struct array of ``name``, ``idx``
            %       pairs of indexed named subsets of variables
            %
            % Output:
            %   nv (integer) : total number of elements in the variable
            %       sub-vector specified by ``varsets``.
            %
            % See also varsets_cell2struct.

            persistent sn;
            if isempty(vs)
                nv = obj.N;
            else
                nv = 0;
            
                %% calls to substruct() are relatively expensive, so we pre-build the
                %% struct for addressing numeric array fields, updating only
                %% the subscripts before use
                if isempty(sn)
                    sn = struct('type', {'.', '()'}, 'subs', {'', 1});
                end
            
                for v = 1:length(vs)
                    idx = vs(v).idx;
                    if isempty(idx)
                        N = obj.idx.N.(vs(v).name);
                    else
                        % (calls to substruct() are relatively expensive ...
                        % sn = substruct('.', vs(v).name, '()', vs(v).idx);
                        % ... so replace it with these more efficient lines)
                        sn(1).subs = vs(v).name;
                        sn(2).subs = idx;
                        N = subsref(obj.idx.N, sn);
                    end
                    nv = nv + sum(N(:));
                end
            end
        end
    end     %% methods

    methods (Static)
        function vs = varsets_cell2struct(vs)
            % varsets_cell2struct - Convert ``varsets`` from cell array to struct array.
            % ::
            %
            %   varsets = mp.sm_variable.varsets_cell2struct(varsets)
            %
            % Converts ``varsets`` from a cell array to a struct array,
            % if necessary.
            %
            % Input:
            %   varsets (cell or struct array) : cell array of names of
            %       variable subsets, or a struct array of ``name``, ``idx``
            %       pairs of indexed named subsets of variables
            %
            % Output:
            %   varsets (struct array) : struct array of ``name``, ``idx``
            %       pairs of indexed named subsets of variables
            %
            % See also varsets_len.

            %% convert varsets from cell to struct array if necessary
            if ~isempty(vs) && iscell(vs)
                empty_cells = cell(1, length(vs));
                [empty_cells{:}] = deal({});    %% empty cell arrays
                vs = struct('name', vs, 'idx', empty_cells);
            end
        end
    end     %% methods (Static)
end         %% classdef
