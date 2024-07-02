classdef sm_nln_cost < mp.set_manager
% mp.sm_nln_cost -  MP Set Manager class for general nonlinear costs.
% ::
%
%   nlc = mp.sm_nln_cost()
%   nlc = mp.sm_nln_cost(label)
%
% MP Set Manager class for general nonlinear costs :math:`f(\x)` where
% :math:`\x` is an :math:`n_x \times 1` vector.
%
% Manages general nonlinear cost sets and their indexing.
%
% By convention, ``nlc`` is the variable name used for mp.sm_nln_cost objects.
%
% mp.sm_nln_cost Methods:
%   * sm_nln_cost - constructor
%   * add - add a subset of general nonlinear costs
%   * params - return general nonlinear cost parameters
%   * set_params - modify general nonlinear cost parameter data
%   * eval - evaluate individual or full set of general nonlinear costs
%
% See also mp.set_manager.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

    methods
        function obj = sm_nln_cost(varargin)
            % Constructor.
            % ::
            %
            %   sm = mp.set_manager(label)

            es = struct();  %% empty struct
            obj@mp.set_manager(varargin{:});
            obj.data = struct( ...
                'fcn', es, ...
                'vs', es );
        end

        function obj = add(obj, var, name, idx, varargin)
            % Add a subset of nonlinear costs.
            % ::
            %
            %   nlc.add(var, name, N, fcn);
            %   nlc.add(var, name, N, fcn, vs);
            %
            %   nlc.add(var, name, idx_list, N, fcn);
            %   nlc.add(var, name, idx_list, N, fcn, vs);
            %
            % Add a named, and possibly indexed, subset of general nonlinear
            % costs :math:`f(\x)` to the set, where :math:`\x` is an
            % :math:`n_x \times 1` vector made up of the variables specified
            % in the optional ``vs`` *(in the order given)*. This allows the
            % cost function to be defined in terms of only the relevant
            % variables without the need to manually account for the locations
            % of other variable sets.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   name (char array) : name of subset/block of costs to add
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of costs to add (for an indexed subset)
            %   N (integer) : dimension of cost function; currently this must
            %       equal 1; vector-valued cost functions are not yet
            %       implemented
            %   fcn (function handle) : handle to function that evaluates the
            %       cost :math:`f(\x)`, its gradient :math:`\der{f}{\x}`, and
            %       Hessian :math:`\der{^2 f}{\x^2}` as described below
            %   vs (cell or struct array) : *(optional, default* ``{}`` *)*
            %       variable set defining vector :math:`\x` for this
            %       cost subset; can be either a cell array of names of
            %       variable subsets, or a struct array of ``name``, ``idx``
            %       pairs of indexed named subsets of variables; order of
            %       ``vs`` determines order of blocks in :math:`\x`; if
            %       empty, :math:`\x` is assumed to be the full variable vector
            %
            % **Cost Function Implmentation** : ``fcn``
            %
            % For a cost function :math:`f(\x)`, the ``fcn`` input should point
            % to a function with the following interface::
            %
            %   f = fcn(x)
            %   [f, df] = fcn(x)
            %   [f, df, d2f] = fcn(x)
            %
            % Input:
            %   x (double or cell array of doubles) : variable :math:`\x`
            %       corresponding to the full variable vector, or to a set of
            %       sub-vectors defined by the variable set provided in ``vs``
            %
            % Outputs:
            %   f (double) : scalar value of cost function :math:`f(\x)`
            %   df (double) : *(optional)* :math:`n_x \times 1` cost gradient,
            %       transpose of :math:`f_\x = \der{f}{\x}`
            %   d2f (double) : *(optional)* :math:`n_x \times n_x` Hessian
            %       matrix :math:`f_{\x\x} = \der{}{\x}(\trans{f_\x})`
            %
            % The input argument ``x`` takes one of two forms. If the
            % cost set is added with varset input ``vs`` empty or missing, then
            % ``x`` will be the full variable vector. Otherwise it will be a
            % cell array of vectors corresponding to the variable sets
            % specified in ``vs``.
            %
            % Examples::
            %
            %   fcn1 = @(x)my_cost_function1(x, other_args)
            %   fcn2 = @(x)my_cost_function2(x, other_args)
            %   nlc.add(var, 'mycost1', 1, fcn1);
            %   nlc.add(var, 'mycost2', 1, fcn2, {'Vm', 'Pg', 'Qg', 'z'});
            %
            %   nlc.init_indexed_name('c', {2, 3});
            %   for i = 1:2
            %     for j = 1:3
            %       nlc.add(var, 'c', {i, j}, 1, fcn(i,j), ...);
            %     end
            %   end
            %
            % See also params, set_params, eval.

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
            vs = {};
            if nargs >= 1
                fcn = args{1};
                if nargs >= 2
                    vs = args{2};
                end
            end

            if N ~= 1
                error('mp.sm_nln_cost.add: not yet implemented for vector valued functions (i.e. N currently must equal 1)');
            end

            %% convert varsets from cell to struct array if necessary
            vs = mp.sm_variable.varsets_cell2struct(vs);

            %% call parent to handle standard indexing
            if isempty(idx)
                add@mp.set_manager(obj, name, N, args{:});
            else
                add@mp.set_manager(obj, name, idx, N, args{:});
            end

            %% assign data
            if isempty(idx)
                obj.data.fcn.(name)  = fcn;
                obj.data.vs.(name) = vs;
            else
                %% calls to substruct() are relatively expensive, so we pre-build the
                %% struct for addressing cell array fields
                %% sc = substruct('.', name, '{}', idx);
                sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
                obj.data.fcn  = subsasgn(obj.data.fcn, sc, fcn);
                obj.data.vs   = subsasgn(obj.data.vs, sc, vs);
            end
        end

        function [N, fcn, vs] = params(obj, var, name, idx)
            % Return general nonlinear cost parameters.
            % ::
            %
            %   [N, fcn] = nlc.params(var, name)
            %   [N, fcn] = nlc.params(var, name, idx_list)
            %   [N, fcn, vs] = nlc.params(...)
            %
            % Returns the parameters for the general nonlinear cost subset
            % corresponding to the name or name and index list provided. The
            % parameters are those supplied via the add method when the
            % subset was added.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   name (char array) : name of subset
            %   idx_list (cell array) : *(optional)* index list for subset
            %
            % Outputs:
            %   N (integer) : dimension of cost function; currently equals 1;
            %       vector-valued cost functions are not yet implemented
            %   fcn (function handle) : handle to function that evaluates the
            %       cost :math:`f(\x)`, its gradient :math:`\trans{f_\x}`, and
            %       Hessian :math:`f_{\x\x}`; see add() for details
            %   vs (struct array) : variable set, ``name``, ``idx`` pairs
            %       specifying the set of variables defining vector :math:`\x`
            %       for this cost subset; order of ``vs`` determines
            %       order of blocks in :math:`\x`
            %
            % See also add, eval.

            if nargin < 4
                idx = {};
            end

            if isempty(idx)
                dims = size(obj.idx.i1.(name));
                if prod(dims) ~= 1
                    error('mp.sm_nln_cost.params: general nonlinear cost set ''%s'' requires an IDX_LIST arg', name)
                end
                N = obj.idx.N.(name);
                fcn = obj.data.fcn.(name);
                vs = obj.data.vs.(name);
            else
                %% calls to substruct() are relatively expensive, so we pre-build the
                %% structs for addressing cell and numeric array fields
                %% sn = substruct('.', name, '()', idx);
                %% sc = substruct('.', name, '{}', idx);
                sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
                sn = sc; sn(2).type = '()';                             %% num array field

                N = subsref(obj.idx.N, sn);
                fcn = subsref(obj.data.fcn, sc);
                vs = subsref(obj.data.vs, sc);
            end
        end

        function obj = set_params(obj, var, name, idx, params, vals)
            % Modify general nonlinear cost parameter data.
            % ::
            %
            %   nlc.set_params(name, params, vals)
            %   nlc.set_params(name, idx, params, vals)
            %
            % This method can be used to modify parameters for an existing
            % subset of general nonlinear costs.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   name (char array) : name of subset/block of general nonlinear
            %       costs to modify
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of general nonlinear costs to modify (for an indexed subset)
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
            % Valid parameter names are ``N``, ``fcn``, ``vs``.
            %
            % Examples::
            %
            %   nlc.set_params(var, 'y', {2,3}, {'fcn'}, {fcn});
            %   nlc.set_params(var, 'Pg', 'all', {N, fcn, vs});
            %
            % See also add, params.

            if nargin < 6
                vals = params;
                params = idx;
                idx = {};
            end

            %% create default list of parameters to update based on set type & inputs
            default_params = {'N', 'fcn', 'vs'};

            %% standardize provided arguments in cell arrays params, vals
            [is_all, np, params, vals] = ...
                obj.set_params_std_args(default_params, params, vals);

            %% get current parameters
            [N0, fcn, vs] = obj.params(var, name, idx);
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
                error('mp.sm_nln_cost.set_params: dimension change for ''%s'' not allowed except for ''all''', obj.nameidxstr(name, idx));
            end

            %% vector valued costs not yet implemented
            if N ~= 1
                error('mp.sm_nln_cost.set_params: vector value for ''%s'' not yet implemented', obj.nameidxstr(name, idx));
            end

            %% convert vs to struct
            if u.vs
                p.vs = mp.sm_variable.varsets_cell2struct(p.vs);
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

            %% update dimensions and indexing, if necessary
            dN = N - N0;
            if is_all && dN
                obj.set_params_update_dims(dN, name, idx);
            end
        end

        function [f, df, d2f] = eval(obj, var, x, name, idx)
            % Evaluate individual or full set of general nonlinear costs.
            % ::
            %
            %   f = nlc.eval(var, x ...)
            %   [f, df] = nlc.eval(var, x ...)
            %   [f, df, d2f] = nlc.eval(var, x ...)
            %   [f, df, d2f] = nlc.eval(var, x, name)
            %   [f, df, d2f] = nlc.eval(var, x, name, idx_list)
            %
            % For a given value of the variable vector :math:`\x`, this method
            % evaluates the general nonlinear cost function and optionally
            % its derivatives for an individual subset, if name or name and
            % index list are provided, otherise, for the full set of costs.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   x (double) : full :math:`n_x \times 1` variable vector :math:`\x`
            %   name (char array) : *(optional)* name of subset/block of
            %       general nonlinear costs to evaluate
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of general nonlinear costs to evaluate (for an indexed
            %       subset)
            %
            % Outputs:
            %   f (double) : scalar value of cost function :math:`f(\x)`
            %   df (double) : *(optional)* :math:`n_x \times 1` cost gradient,
            %       transpose of :math:`f_\x = \der{f}{\x}`
            %   d2f (double) : *(optional)* :math:`n_x \times n_x` Hessian
            %       matrix :math:`f_{\x\x} = \der{}{\x}(\trans{f_\x})`
            %
            % See also add, params.

            %% initialize
            if nargin < 5
                idx = {};
            end

            if nargin < 4                       %% full set
                f = 0;
                nx = var.N;             %% number of variables
                if nargout > 1
                    df = zeros(nx, 1);  %% column vector
                    if nargout > 2
                        d2f = sparse(nx, nx);   %% symmetric, so transpose (used for speed)
                    end                         %% is equal
                end

                for k = 1:obj.NS
                    name = obj.order(k).name;
                    idx  = obj.order(k).idx;
                    [N, fcn, vs] = obj.params(var, name, idx);
                    if N ~= 1
                        error('mp.sm_nln_cost.eval: not yet implemented for vector valued functions');
                    end
                    xx = var.varsets_x(x, vs);
                    if nargout == 3
                        [fk, dfk, d2fk] = fcn(xx);  %% evaluate kth cost, gradient, Hessian
                    elseif nargout == 2
                        [fk, dfk] = fcn(xx);        %% evaluate kth cost and gradient
                    else
                        fk = fcn(xx);               %% evaluate kth cost
                    end

                    f = f + fk;

                    if nargout > 1              %% assemble gradient
                        nk = length(dfk);
                        if isempty(vs)          %% all rows of x
                            if nk == nx
                                df = df + dfk;
                                if nargout > 2  %% assemble Hessian
                                    d2f = d2f + d2fk;
                                end
                            else                %% must have added vars since adding
                                                %% this cost set
                                df(1:nk) = df(1:nk) + dfk;
                                if nargout > 2      %% assemble Hessian
                                    d2fk_all_cols = sparse(nk, nx);
                                    d2fk_all_cols(:, 1:nk) = d2fk;
                                    d2fk_full = sparse(nx, nx);
                                    d2f(:, 1:nk) = d2f(:, 1:nk) + d2fk_all_cols';
                                end
                            end
                        else                    %% selected rows of x
                            jj = var.varsets_idx(vs);   %% indices for var set
                            df(jj) = df(jj) + dfk;
                            if nargout > 2      %% assemble Hessian
                                d2fk_all_cols = sparse(nk, nx);
                                d2fk_all_cols(:, jj) = d2fk;
                                d2fk_full = sparse(nx, nx);
                                d2f(:, jj) = d2f(:, jj) + d2fk_all_cols';
                            end
                        end
                    end
                end
            else                                %% individual named set
                dims = size(obj.idx.i1.(name));
                if ~isempty(idx) || prod(dims) == 1 %% indexed, or simple named set
                    [N, fcn, vs] = obj.params(var, name, idx);
                    if N ~= 1
                        error('mp.sm_nln_cost.eval: not yet implemented for vector valued functions');
                    end
                    xx = var.varsets_x(x, vs);
                    if nargout == 3
                        [f, df, d2f] = fcn(xx);     %% evaluate kth cost, gradient, Hessian
                    elseif nargout == 2
                        [f, df] = fcn(xx);          %% evaluate kth cost and gradient
                    else
                        f = fcn(xx);                %% evaluate kth cost
                    end
                elseif nargout == 1             %% sum over all indices for name
                    done = 0;
                    f = 0;          %% initialize cumulative cost
                    idx = num2cell(ones(size(dims))); %% initialize idx
                    while ~done     %% call eval_nln_cost() recursively
                        f = f + sum(obj.eval(var, x, name, idx));

                        %% increment idx
                        D = length(dims);
                        idx{D} = idx{D} + 1;    %% increment last dimension
                        for d = D:-1:2          %% increment next dimension, if necessary
                            if idx{d} > dims(d)
                                idx{d} = 1;
                                idx{d-1} = idx{d-1} + 1;
                            end
                        end
                        if idx{1} > dims(1)     %% check if done
                            done = 1;
                        end
                    end
                else
                    error('mp.sm_nln_cost.eval: general nonlinear cost set ''%s'' requires an IDX_LIST arg when requesting DF output', name)
                end
            end
        end
    end     %% methods
end         %% classdef