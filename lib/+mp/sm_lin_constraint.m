classdef sm_lin_constraint < mp.set_manager
% mp.sm_lin_constraint -  MP Set Manager class for linear constraints.
% ::
%
%   lin = mp.sm_lin_constraint()
%   lin = mp.sm_lin_constraint(label)
%
% MP Set Manager class for linear constraints of the form
%
% .. math:: \l \le \AA \x \le \u
%   :label: eq_lin_form
%
% Manages constraint parameters :math:`\AA, \l, \u`, along with indexing.
%
% By convention, ``lin`` is the variable name used for mp.sm_lin_constraint
% objects.
%
% mp.sm_lin_constraint Properties:
%   * cache - struct for caching aggregated parameters for the set
%
% mp.sm_lin_constraint Methods:
%   * sm_lin_constraint - constructor
%   * add - add a subset of linear constraints, with parameters :math:`\AA, \l, \u`
%   * params - build and return linear constraint parameters :math:`\AA, \l, \u`
%   * set_params - modify linear constraint parameter data
%   * eval - evaluate individual or full set of linear constraints
%   * parse_soln - parse solution for linear constraints
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
        % struct for caching aggregated parameters for linear constraints
        cache = [];
    end     %% properties

    methods
        function obj = sm_lin_constraint(varargin)
            % Constructor.
            % ::
            %
            %   sm = mp.sm_lin_constraint(label)

            es = struct();  %% empty struct
            obj@mp.set_manager(varargin{:});
            obj.data = struct( ...
                'A', es, ...
                'l', es, ...
                'u', es, ...
                'tr', es, ...
                'vs', es );
        end

        function obj = add(obj, var, name, idx, varargin)
            % Add a subset of linear constraints, with parameters :math:`\AA, \l, \u`.
            % ::
            %
            %   lin.add(var, name, A, l, u);
            %   lin.add(var, name, A, l, u, vs);
            %   lin.add(var, name, A, l, u, vs, tr);
            %
            %   lin.add(var, name, idx_list, A, l, u);
            %   lin.add(var, name, idx_list, A, l, u, vs);
            %   lin.add(var, name, idx_list, A, l, u, vs, tr);
            %
            % Add a named, and possibly indexed, subset of linear constraints
            % to the set, of the form :math:`\l \le \AA \x \le \u`, where
            % :math:`\x` is a vector made up of the variables specified in
            % the optional ``vs`` *(in the order given)*. This allows the
            % :math:`\AA` matrix to be defined in terms of only the relevant
            % variables without the need to manually create a lot of properly
            % located zero columns.
            %
            % There is also the option to take and store the transpose of the
            % constraint matrix :math:`\AA`. For constraints with many more
            % columns than rows, this can sometimes save significant memory.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   name (char array) : name of subset/block of constraints to add
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of constraints to add (for an indexed subset)
            %   A (double) : constraint coefficient matrix :math:`\AA`
            %   l (double) : *(optional, default =* ``-Inf`` *)* constraint
            %       left-hand side vector :math:`\l`, or scalar which is
            %       expanded to a vector
            %   u (double) : *(optional, default =* ``Inf`` *)* constraint
            %       right-hand side vector :math:`\u`, or scalar which is
            %       expanded to a vector
            %   vs (cell or struct array) : *(optional, default* ``{}`` *)*
            %       variable set defining vector :math:`\x` for this
            %       constraint subset; can be either a cell array of names of
            %       variable subsets, or a struct array of ``name``, ``idx``
            %       pairs of indexed named subsets of variables; order of
            %       ``vs`` determines order of blocks in :math:`\x`; if
            %       empty, :math:`\x` is assumed to be the full variable vector
            %   tr (boolean) : *(optional, default* ``false`` *)* if true, it
            %       means that :math:`\trans{\AA}` is supplied/stored rather
            %       than :math:`\AA`
            %
            % Examples::
            %
            %   lin.add(var, 'vl', Avl, lvl, uvl, {'Pg', 'Qg'});
            %
            %   lin.init_indexed_name('R', {2, 3});
            %   for i = 1:2
            %       for j = 1:3
            %           lin.add(var, 'R', {i, j}, A{i,j}, ...);
            %       end
            %   end
            %
            % See also params, set_params, eval.

            %% set up default args
            if iscell(idx)          %% indexed named set
                A = varargin{1};
                args = varargin(2:end);
            else                    %% simple named set
                A = idx;
                idx = {};
                args = varargin;
            end
            nargs = length(args);

            %% prepare data
            l = []; u = []; vs = {}; tr = 0;
            if nargs >= 1
                l = args{1};
                if nargs >= 2
                    u = args{2};
                    if nargs >= 3
                        vs = args{3};
                        if nargs >= 4
                            tr = args{4};
                        end
                    end
                end
            end

            if tr
                [M, N] = size(A);
            else
                [N, M] = size(A);
            end

            %% call parent to handle standard indexing
            if isempty(idx)
                add@mp.set_manager(obj, name, N, args{:});
            else
                add@mp.set_manager(obj, name, idx, N, args{:});
            end

            if isempty(l)                   %% default l is -Inf
                l = -Inf(N, 1);
            elseif N > 1 && length(l) == 1  %% expand from scalar as needed
                l = l * ones(N, 1);
            end
            if isempty(u)                   %% default u is Inf
                u = Inf(N, 1);
            elseif N > 1 && length(u) == 1  %% expand from scalar as needed
                u = u * ones(N, 1);
            end

            %% check sizes
            if size(l, 1) ~= N || size(u, 1) ~= N
                error('mp.sm_lin_constraint.add: sizes of A, l and u must match');
            end

            %% convert varsets from cell to struct array if necessary
            vs = mp.sm_variable.varsets_cell2struct(vs);
            nv = var.varsets_len(vs);   %% number of variables

            %% check consistency of varsets with size of A
            if M ~= nv
                error('mp.sm_lin_constraint.add: number of columns of A does not match\nnumber of variables, A is %d x %d, nv = %d\n', N, M, nv);
            end

            %% assign data
            if isempty(idx)
                obj.data.A.(name)  = A;
                obj.data.l.(name)  = l;
                obj.data.u.(name)  = u;
                obj.data.tr.(name)  = tr;
                obj.data.vs.(name) = vs;
            else
                %% calls to substruct() are relatively expensive, so we pre-build the
                %% struct for addressing cell array fields
                %% sc = substruct('.', name, '{}', idx);
                sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
                obj.data.A  = subsasgn(obj.data.A, sc, A);
                obj.data.l  = subsasgn(obj.data.l, sc, l);
                obj.data.u  = subsasgn(obj.data.u, sc, u);
                obj.data.tr  = subsasgn(obj.data.tr, sc, tr);
                obj.data.vs = subsasgn(obj.data.vs, sc, vs);
            end
            if ~isempty(obj.cache)  %% clear cache of aggregated params
                obj.cache = [];
            end
        end

        function [A, l, u, vs, i1, iN, tr] = params(obj, var, name, idx)
            % Build and return linear constraint parameters :math:`\AA, \l, \u`.
            % ::
            %
            %   [A, l, u] = lin.params(var)
            %   [A, l, u] = lin.params(var, name)
            %   [A, l, u] = lin.params(var, name, idx_list)
            %   [A, l, u, vs] = lin.params(...)
            %   [A, l, u, vs, i1, in] = lin.params(...)
            %   [A, l, u, vs, i1, in, tr] = lin.params(name ...)
            %
            % With no input parameters, it assembles and returns the parameters
            % for the aggregate linear constraints from all linear constraint
            % sets added using add(). The values of these parameters are
            % cached for subsequent calls. The parameters are
            % :math:`\AA, \l, \u`, where the linear constraint is of the form
            % in :eq:`eq_lin_form`.
            %
            % If a name or name and index list are provided then it simply
            % returns the parameters for the corresponding set. It can also
            % optionally return the variable sets used by this constraint set
            % (the size of :math:`\AA` will be consistent with this variable
            % set), the starting and ending row indices of the subset within
            % the full aggregate constraint matrix, and flag indicating whether
            % the transpose of the constraint matrix was  supplied/stored for
            % this subset.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   name (char array) : *(optional)* name of subset
            %   idx_list (cell array) : *(optional)* index list for subset
            %
            % Outputs:
            %   A (double) : constraint coefficient matrix :math:`\AA`
            %   l (double) : constraint left-hand side vector :math:`\l`
            %   u (double) : constraint right-hand side vector :math:`\u`
            %   vs (struct array) : variable set, ``name``, ``idx`` pairs
            %       specifying the set of variables defining vector :math:`\x`
            %       for this constraint subset; order of ``vs`` determines
            %       order of blocks in :math:`\x`
            %   i1 (integer) : index of 1st row of specified subset in full set
            %   iN (integer) : index of last row of specified subset in full set
            %   tr (boolean) : if true, it means that :math:`\trans{\AA}` was
            %       supplied/stored rather than :math:`\AA`
            %
            % Examples::
            %
            %   [A, l, u] = lin.params(var);
            %   [A, l, u, vs, i1, iN] = lin.params(var, 'Pmis');
            %
            % See also add.

            if nargin > 2       %% individual set
                if nargin < 4
                    idx = {};
                end
                if isempty(idx)
                    if numel(obj.idx.i1.(name)) == 1     %% simple named set
                        A = obj.data.A.(name);
                        l = obj.data.l.(name);
                        u = obj.data.u.(name);
                        if nargout > 3
                            vs = obj.data.vs.(name);
                            if nargout > 5
                                i1 = obj.idx.i1.(name);      %% starting row index
                                iN = obj.idx.iN.(name);      %% ending row index
                                if nargout > 6
                                    if isfield(obj.data, 'tr')
                                        tr = obj.data.tr.(name);
                                    else
                                        tr = 0;
                                    end
                                end
                            end
                        end
                    else                                    %% indexing required
                        error('mp.sm_lin_constraint.params: linear constraint set ''%s'' requires an IDX_LIST arg', name);
                    end
                else
                    % (calls to substruct() are relatively expensive ...
                    % s = substruct('.', name, '{}', idx);
                    % ... so replace it with these more efficient lines)
                    sc = struct('type', {'.', '{}'}, 'subs', {name, idx});
                    A = subsref(obj.data.A, sc);
                    l = subsref(obj.data.l, sc);
                    u = subsref(obj.data.u, sc);
                    if nargout > 3
                        vs = subsref(obj.data.vs, sc);
                        if nargout > 5
                            sn = sc; sn(2).type = '()';     %% num array field
                            i1 = subsref(obj.idx.i1, sn);   %% starting row index
                            iN = subsref(obj.idx.iN, sn);   %% ending row index
                            if nargout > 6
                                if isfield(obj.data, 'tr')
                                    tr = subsref(obj.data.tr, sc);
                                else
                                    tr = 0;
                                end
                            end
                        end
                    end
                end
            else                %% aggregate
                cache = obj.cache;
                if isempty(cache)       %% build the aggregate
                    nx = var.N;         %% number of variables
                    nlin = obj.N;       %% number of linear constraints
                    if obj.NS < 25 || obj.NS < 100 && nx < 300
                        %% METHOD 1: Add sparse matrices (original method)
                        At = sparse(nx, nlin);  %% transpose of constraint matrix
                        u = Inf(nlin, 1);       %% upper bound
                        l = -u;                 %% lower bound

                        %% fill in each piece
                        for k = 1:obj.NS
                            name = obj.order(k).name;
                            idx  = obj.order(k).idx;
                            [Ak, lk, uk, vs, i1, iN, tr] = obj.params(var, name, idx);
                            if tr
                                [nk, mk] = size(Ak);        %% size of Ak
                            else
                                [mk, nk] = size(Ak);        %% size of Ak
                            end
                            if mk
                                Akt_full = sparse(nx, nlin);
                                if isempty(vs)
                                    if nk == nx     %% full size
                                        if tr
                                            Akt_full(:, i1:iN) = Ak;
                                        else
                                            Akt_full(:, i1:iN) = Ak';
                                        end
                                    else            %% vars added since adding this constraint set
                                        if tr
                                            Akt_full(1:nk, i1:iN) = Ak;
                                        else
                                            Ak_all_cols = sparse(mk, nx);
                                            Ak_all_cols(:, 1:nk) = Ak;
                                            Akt_full(:, i1:iN) = Ak_all_cols';
                                        end
                                    end
                                else
                                    jj = var.varsets_idx(vs);   %% indices for var set
                                    if tr
                                        Akt_full(jj, i1:iN) = Ak;
                                    else
                                        Ak_all_cols = sparse(mk, nx);
                                        Ak_all_cols(:, jj) = Ak;
                                        Akt_full(:, i1:iN) = Ak_all_cols';
                                    end
                                end
                                At = At + Akt_full;
                                l(i1:iN) = lk;
                                u(i1:iN) = uk;
                            end
                        end
                        A = At';
                    else
                        %% METHOD 2: construct using single call to sparse()
                        A_ijv = cell(obj.NS, 3); %% indices/values to construct A
                        u = Inf(nlin, 1);       %% upper bound
                        l = -u;                 %% lower bound

                        %% fill in each piece
                        for k = 1:obj.NS
                            name = obj.order(k).name;
                            idx  = obj.order(k).idx;
                            [Ak, lk, uk, vs, i1, iN, tr] = obj.params(var, name, idx);
                            if tr
                                [nk, mk] = size(Ak);        %% size of Ak
                            else
                                [mk, nk] = size(Ak);        %% size of Ak
                            end
                            if mk
                                % find nonzero sub indices and values
                                if tr
                                    [j, i, v] = find(Ak);
                                    if nk == 1  %% force col vectors for single col Ak
                                        i = i'; j = j'; v = v';
                                    end
                                else
                                    [i, j, v] = find(Ak);
                                    if mk == 1  %% force col vectors for single row Ak
                                        i = i'; j = j'; v = v';
                                    end
                                end

                                if isempty(vs)
                                    A_ijv(k,:) = {i+(i1-1), j, v};
                                else
                                    jj = var.varsets_idx(vs)';  %% indices for var set
                                    A_ijv(k,:) = {i+(i1-1), jj(j), v};
                                end
                                l(i1:iN) = lk;
                                u(i1:iN) = uk;
                            end
                        end
                        A = sparse( vertcat(A_ijv{:,1}), ...
                                    vertcat(A_ijv{:,2}), ...
                                    vertcat(A_ijv{:,3}), nlin, nx);
                    end

                    %% cache aggregated parameters
                    obj.cache = struct('A', A, 'l', l, 'u', u);
                else                    %% return cached values
                    A = cache.A;
                    l = cache.l;
                    u = cache.u;
                end
                if nargout > 3
                    vs = {};
                end
            end
        end

        function obj = set_params(obj, var, name, idx, params, vals)
            % Modify linear constraint parameter data.
            % ::
            %
            %   lin.set_params(var, name, params, vals)
            %   lin.set_params(var, name, idx, params, vals)
            %
            % This method can be used to modify parameters for an existing
            % subset of linear constraints.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   name (char array) : name of subset/block of linear constraints
            %       to modify
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of linear constraints to modify (for an indexed subset)
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
            % Valid parameter names are ``A``, ``l``, ``u``, ``vs``, ``tr``.
            %
            % Examples::
            %
            %   lin.set_params(var, 'y', {2,3}, {'l', 'u'}, {l, u});
            %   lin.set_params(var, 'Pmis', 'all', {A, l, u, vs});
            %
            % See also add, params.

            if nargin < 6
                vals = params;
                params = idx;
                idx = {};
            end

            %% create default list of parameters to update based on set type & inputs
            default_params = {'A', 'l', 'u', 'vs', 'tr'};

            %% standardize provided arguments in cell arrays params, vals
            [is_all, np, params, vals] = ...
                obj.set_params_std_args(default_params, params, vals);

            %% get current parameters
            [A, l, u, vs, ~, ~, tr] = obj.params(var, name, idx);
            if tr
                [M0, N0] = size(A);
            else
                [N0, M0] = size(A);
            end
            if isempty(vs), vs = {vs}; end
            p = struct('A', A, 'l', l, 'u', u, 'vs', vs, 'tr', tr); %% current parameters
            u = struct('A', 0, 'l', 0, 'u', 0, 'vs',  0, 'tr', 0);  %% which ones to update

            %% replace with new parameters
            for k = 1:np
                p.(params{k}) = vals{k};
                u.(params{k}) = 1;
            end

            %% set missing default params for 'all'
            if p.tr
                [M, N] = size(p.A);
            else
                [N, M] = size(p.A);
            end
            if is_all
                u.A = 1;            %% always update A
                u.l = 1;            %% alwaus update l
                if np < 5
                    if p.tr
                        p.tr = 0;
                        u.tr = 1;       %% update tr
                        [N, M] = deal(M, N);
                    end
                    if np < 4
                        p.vs = {};
                        u.vs = 1;       %% update vs
                        if np < 3
                            p.u = Inf(N, 1);
                            u.u = 1;    %% update u
                        end
                    end
                end
            end

            %% check consistency of parameters
            %% no dimension change unless 'all'
            if N ~= N0 && ~is_all
                error('mp.sm_lin_constraint.set_params: dimension change for ''%s'' not allowed except for ''all''', obj.nameidxstr(name, idx));
            end
            %% no transpose change unless providing A
            if u.tr && ~u.A
                error('mp.sm_lin_constraint.set_params: update to ''tr'' for ''%s'' requires update to ''A''', obj.nameidxstr(name, idx));
            end

            %% check sizes of new values of l, u
            for pn = {'l', 'u'}
                if u.(pn{1})
                    nn = length(p.(pn{1}));
                    if nn ~= N
                        if nn == 0
                            switch pn{1}
                                case 'l'
                                    p.(pn{1}) = -Inf(N, 0);
                                case 'u'
                                    p.(pn{1}) =  Inf(N, 0);
                            end
                        elseif nn == 1
                            p.(pn{1}) = p.(pn{1}) * ones(N, 1);   %% expand from scalar
                        else
                            error('mp.sm_lin_constraint.set_params: parameter ''%s'' ''%s'' should have length %d (or 1)', obj.nameidxstr(name, idx), pn{1}, N);
                        end
                    end
                end
            end

            %% check consistency of A and vs
            if u.A || u.vs
                p.vs = mp.sm_variable.varsets_cell2struct(p.vs);
                nv = var.varsets_len(p.vs);     %% number of variables
                if M ~= nv
                    error('mp.sm_lin_constraint.set_params: for ''%s'' number of columns of ''A'' (%d) must be consistent with ''vs'' (%d)', obj.nameidxstr(name, idx), M, nv);
                end
            end

            %% assign new parameters
            if isempty(idx)     %% simple named set
                for k = 1:length(default_params)
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
                for k = 1:length(default_params)
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

        function [Ax_u, l_Ax, A] = eval(obj, var, x, name, idx)
            % Evaluate individual or full set of linear constraints.
            % ::
            %
            %   Ax_u = lin.eval(var, x)
            %   Ax_u = lin.eval(var, x, name)
            %   Ax_u = lin.eval(var, x, name, idx_list)
            %   [Ax_u, l_Ax] = lin.eval(...)
            %   [Ax_u, l_Ax, A] = lin.eval(...)
            %
            % For a given value of the variable vector :math:`\x`, this method
            % evaluates the linear constraints for an individual subset, if
            % name or name and index list are provided, otherise, for the full
            % set of constraints.
            %
            % Returns :math:`\AA \x - \u`, and optionally :math:`\l - \AA \x`
            % and :math:`\AA`.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   x (double) : full :math:`n_x \times 1` variable vector :math:`\x`
            %   name (char array) : name of subset/block of linear constraints
            %       to evaluate
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of linear constraints to evaluate (for an indexed subset)
            %
            % Outputs:
            %   Ax_u (double) : value of :math:`\AA \x - \u`
            %   l_Ax (double) : *(optional)* value of :math:`\l - \AA \x`
            %   A (double) : *(optional)* constraint matrix :math:`\AA`
            %
            % See also add, params.

            if obj.N
                %% collect cost parameters
                if nargin < 4                       %% full set
                    [A, l, u, vs] = obj.params(var);
                    tr = 0;
                    N = 1;
                elseif nargin < 5 || isempty(idx)   %% name, no idx provided
                    dims = size(obj.idx.i1.(name));
                    if prod(dims) == 1              %% simple named set
                        [A, l, u, vs, ~, ~, tr] = obj.params(var, name);
                        N = obj.get_N(name);
                    else
                        error('mp.sm_lin_constraint.eval: linear constraint set ''%s'' requires an IDX_LIST arg', name)
                    end
                else                                %% indexed named set
                    [A, l, u, vs, ~, ~, tr] = obj.params(var, name, idx);
                    N = obj.get_N(name, idx);
                end

                %% assemble appropriately-sized x vector
                xx = var.varsets_x(x, vs, 'vector');

                %% compute constraints
                if tr
                    Ax = (xx' * A)';
                    if nargout > 2
                        A = A';
                    end
                else
                    Ax = A * xx;
                end
                Ax_u = Ax - u;
                if nargout > 1
                    l_Ax = l - Ax;
                end
            else
                Ax_u = [];
                if nargout > 1
                    l_Ax = [];
                    if nargout > 2
                        A = [];
                    end
                end
            end
        end

        function ps = parse_soln(obj, soln)
            % Parse solution for linear constraints.
            % ::
            %
            %   ps = lin.parse_soln(soln)
            %
            % Parse a full solution struct into parts corresponding to
            % individual linear constraint subsets.
            %
            % Input:
            %   soln (struct) : full solution struct with these fields
            %       (among others):
            %
            %           - ``x`` - variable values
            %           - ``lambda`` - constraint shadow prices, struct with
            %             fields:
            %
            %               - ``eqnonlin`` - nonlinear equality constraints
            %               - ``ineqnonlin`` - nonlinear inequality constraints
            %               - ``mu_l`` - linear constraint lower bounds
            %               - ``mu_u`` - linear constraint upper bounds
            %               - ``lower`` - variable lower bounds
            %               - ``upper`` - variable upper bounds
            %
            % Output:
            %   ps (struct) : parsed solution, struct where each field listed
            %       below is a struct whos names are the names of the relevant
            %       linear constraint subsets and values are scalars for named
            %       sets, arrays for named/indexed sets:
            %
            %           - ``mu_l`` - constraint lower bound shadow prices
            %           - ``mu_u`` - constraint upper bound shadow prices

            ps = [];
            if obj.get_N()
                if isfield(soln.lambda, 'mu_l')
                    if isfield(soln.lambda, 'mu_u')
                        params = struct('src', {soln.lambda.mu_l, soln.lambda.mu_u}, ...
                                        'dst', {'mu_l', 'mu_u'});
                    else
                        params = struct('src', soln.lambda.mu_l, 'dst', 'mu_l');
                    end
                else
                    if isfield(soln.lambda, 'mu_u')
                        params = struct('src', soln.lambda.mu_u, 'dst', 'mu_u');
                    else
                        params = [];
                    end
                end
                if ~isempty(params)
                    ps = obj.parse_soln_fields(params);
                end
            end
        end
    end     %% methods
end         %% classdef
