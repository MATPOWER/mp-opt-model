classdef sm_quad_constraint < mp.set_manager_opt_model
% mp.sm_quad_constraint -  MP Set Manager class for quadratic constraints.
% ::
%
%   qcn = mp.sm_quad_constraint()
%   qcn = mp.sm_quad_constraint(label)
%
% MP Set Manager class for quadratic constraints of the form
%
% .. math:: l(i) <= 1/2 X'*Q{i}*X + C(i,:)*X <= u(i),  i = 1,2,...,NQ 
%   :label: 
%
% Manages quadratic constraint sets and their indexing.
%
% By convention, ``qcn`` is the variable name used for mp.sm_quad_constraint 
% objects.
%
% mp.sm_quad_constraint Properties:
%   * cache - struct for caching aggregated parameters for the set
%
% mp.sm_quad_constraint Methods:
%   * sm_quad_constraint - constructor
%   * add - add a subset of quadratic constraints
%   * params - build and return quadratic constraint parameters
%   * set_params - modify quadratic constraint parameter data
%   * eval - evaluate individual or full set of quadratic constraints
%   * display_soln - display solution values for quadratic constraints
%   * get_soln - fetch solution values for specific named/indexed subsets
%   * parse_soln - parse solution for quadratic constraints
%
% See also mp.set_manager, mp.set_manager_opt_model.

%   MP-Opt-Model
%   Copyright (c) 2019-2023, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia
%
%   This file is part of MP-Opt-Model..
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

    properties
        % struct for caching aggregated parameters for quadratic constraints
        cache = [];
    end

    methods
        function obj = sm_quad_constraint(varargin)
            % Constructor.
            % ::
            %
            %   qcn = mp.sm_quad_constraint(label)
            
            es = struct();  %% empty struct
            obj@mp.set_manager_opt_model(varargin{:});
            obj.data = struct( ...
                'Q', es, ...
                'C', es, ...
                'l', es, ...
                'u', es, ...
                'vs', es);
        end

        function obj = add(obj, var, name, idx, varargin)
            % Add a subset of NQ quadratic constraints.
            % ::
            %
            %   qcn.add(var, name, Q, C, l, u);
            %   qcn.add(var, name, Q, C, l, u, vs);
            %
            %   Indexed Named Sets:
            %   qcn.add(var, name, idx_list, Q, C, l, u);
            %   qcn.add(var, name, idx_list, Q, C, l, u, vs);            
            %
            % Add a named, and possibly indexed, subset of quadratic constraints
            % to the set, of the form :math:``, 
            % where :math:`\x` is a vector made up of the variables specified
            % in the optional ``vs`` *(in the order given)*. This allows the
            % :math:`\AA` matrix an the set of {Qi} matrices to be defined in 
            % terms of only the relevant variables without the need to manually 
            % create a lot of properly located zero columns.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   name (char array or cell array) : name(s) of subset/block of
            %       constraints to add
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of constraints to add (for an indexed subset)
            %   Q (cell vector): NQ x 1 cell array of sparse quadratic matrices
            %       for quadratic constraints. Each element of the array must
            %       have the same size based on whether the full vector of 
            %       variables or an (optional) provided varset.
            %   C (double) : matrix (posibly sparse) of linear term of quadratic
            %        constraints. Each row of the matrix is the linear term of 
            %        each quadratic onstraint.
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
            %
            % Examples::
            %
            %   qcn.add(var, 'my_quad', Q1, C1, l1, u1, {'my_var1', 'my_var2'});
            %
            %   qcn.init_indexed_name('my_set', {3, 2})
            %   for i = 1:3
            %       for j = 1:2
            %           qcn.add(var, 'my_set', {i, j}, Q{i,j}, C{i,j}, l{i,j}, ...);
            %       end
            %   end
            %
            % See also params, set_params, eval.

            %% set up default args
            if isfield(obj.idx.N, name) && ~isscalar(obj.idx.N.(name))
                %% indexed named set
                Q = varargin{1};
                args = varargin(2:end);
            else
                %% simple named set
                Q = idx;
                idx = {};
                args = varargin;
            end
            nargs = length(args);

            %% prepare data
            C = []; l = []; u = []; vs = {};
            if nargs >= 1
                C = args{1};
                if nargs >= 2
                    l = args{2};
                    if nargs >= 3
                        u = args{3};
                        if nargs >= 4
                            vs = args{4};
                        end
                    end
                end
            end

            %% check bounds
            [MQ, NQ] = size(Q);
            [MQi, NQi] = size(Q{1});
            [MC, NC] = size(C);

            if isempty(l)
                l = -inf(MC, 1);
            elseif numel(l) ~= MC
                error('mp.sm_quad_constraint.add: l (%d x 1) must be a column vector with %d elements \n', length(l), MC);
            elseif MC > 1 && numel(l) == 1
                l = l*ones(MC, 1);
            end

            if isempty(u)
                u = inf(MC, 1);
            elseif numel(u) ~= MC
                error('mp.sm_quad_constraint.add: l (%d x 1) must be a column vector with %d elements \n', length(u), MC);
            elseif MC > 1 && numel(u) == 1
                u = u*ones(MC, 1);
            end

            %% convert varsets from cell to struct array if necessary
            vs = mp.sm_variable.varsets_cell2struct(vs);
            nv = var.varsets_len(vs);   %% number of variables

            %% Check parameters
            if MQi ~= NQi
                error('mp.sm_quad_constraint.add: Q_i (%d x %d) must be a square matrix', MQi, NQi);
            end
            if MQ
                if NQ ~= 1
                    error('mp.sm_quad_constraint.add: Q (%d x %d) must be a column cell array (or empty)', MQ, NQ);
                end
            end
            if MC && NC ~= nv
                error('mp.sm_quad_constraint.add: C (%d x %d) must be a matrix with %d columns', MC, NC, nv);
            end
            if MQ
                if NC && NC ~= MQi
                    error('mp.sm_quad_constraint.add: dimensions of Q (%d x %d) and C (%d x %d) are not compatible', MQ, NQ, MC, NC);
                end
                nx = MQi;
            else
                if nv && ~MC
                    error('mp.sm_quad_constraint.add: Q and C cannot both be empty');
                end
                nx = MC;
            end
            N = MC;

            if nx ~= nv
                error('mp.sm_quad_constraint.add: dimensions of Q (%d x %d), C (%d x %d), and K (%d x %d) do not match\nnumber of variables (%d)\n', MQ, NQ, MC, NC, MK, NK, nv);
            end

            %% call parent to handle standard indexing
            if isempty(idx)
                add@mp.set_manager_opt_model(obj, name, N, args{:});
            else
                add@mp.set_manager_opt_model(obj, name, idx, N, args{:});
            end

            %% assign data (wgv: this replaces the old add_named_set method)
            if isempty(idx)
                obj.data.Q.(name)  = Q;
                obj.data.C.(name)  = C;
                obj.data.l.(name)  = l;
                obj.data.u.(name)  = u;
                obj.data.vs.(name) = vs;
            else
                %% calls to substruct() are relatively expensive, so we pre-build the
                %% struct for addressing cell array fields
                %% sc = substruct('.', name, '{}', idx);
                sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
                obj.data.Q  = subsasgn(obj.data.Q, sc, Q);
                obj.data.C  = subsasgn(obj.data.C, sc, C);
                obj.data.l  = subsasgn(obj.data.l, sc, l);
                obj.data.u  = subsasgn(obj.data.u, sc, u);
                obj.data.vs = subsasgn(obj.data.vs, sc, vs);
            end
            if ~isempty(obj.cache)  %% clear cache of aggregated params
                obj.cache = [];
            end
        end

        function [Qblk, C, l, u, vs, i1, iN] = params(obj, var, name, idx, isblk)
            % Returns the constraint parameters for a set of quadratic constraints
            % ::
            %
            %   [Qblk, C, l, u] = qcn.params(var)
            %   [Qblk, C, l, u] = qcn.params(var, name)
            %   [Qblk, C, l, u] = qcn.params(var, name, idx)
            %   [Qblk, C, l, u, vs] = qcn.params(...)
            %   [Qblk, C, l, u, vs, i1, iN] = qcn.params(name, ...)
            %
            % With no input parameters, it assembles and returns the parameters
            % for the aggregate parameters from all quadratic constraint sets added
            % using the method add(). The values of these parameters are cached
            % for subsequent calls. The parameters are Qblk, C, k, l, u, and K. If
            % input isblk is set to 1, it calculated a set of  assembled quadratic 
            % constraints of the form:
            %
            %    F(X) = 1/2 * DIAG( BLKDIAG(X)' * QBLK * BLKDIAG(X) )  +  C * X
            %
            % Here Qblk is the block diagonal matrix formed from the individual
            % quadratic matrices of the set of constraints. When isblk is set to 0,
            % then Qblk is simply the vertical stack of all cell arrays of quadratic 
            % matrices of all sets. If a name or name and index list are provided 
            % then it simply returns the parameters for the corresponding named set, 
            % whether as a block diagonal matrix or not depending on parameter isblk.
            % It can also optionally return the variable sets used by this constraint
            % set (the size of Qblk and C will be consistent with this variable
            % set), and the starting and ending row indices of the subset within
            % the full aggregate constraint matrix.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   name (char array) : *(optional)* name of subset
            %   idx (cell array) : *(optional)* index list for subset
            %   isblk *(optional default: 0)* : indicator for building Qblk parameter 
            %       set to 1 for building a block diagonal matrix with the
            %       matrices stored in the indicated quadratic constraint set
            %       Set to 0 for biulding a cell array formed by stacking vertically
            %       the matrices stored in the indicated contrained set
            %
            % Outputs:
            %   Qblk (double or cell array) : constraint coefficient matrix
            %       or array depending on the value of input isblk
            %   C (double) : contraint coefficient matrix for linear terms
            %   l (double) : constraint left-hand side vector
            %   u (double) : constraint right-hand side vector
            %   vs (struct array) : variable set, ``name``, ``idx`` pairs
            %       specifying the set of variables defining vector :math:`\x`
            %       for this constraint subset; order of ``vs`` determines
            %       order of blocks in :math:`\x`
            %   i1 (integer) : index of 1st row of specified subset in full set
            %   iN (integer) : index of last row of specified subset in full set
            %
            % Examples::
            % 
            %   [Qblk, C, l, u] = qcn.params(var)
            %   [Qblk, C, l, u] = qcn.params(var, 'my_set')
            %
            % See also add.

            if nargin > 2       %% individual set
                if nargin < 5
                    isblk = 0;
                    if nargin < 4
                        idx = {};
                    end
                end
                if isempty(idx)                 %% name, no index provided
                    if numel(obj.idx.i1.(name)) == 1     %% simple named set
                        Qblk = obj.data.Q.(name);  % Cell array of quadratic matrices
                        if isblk               % Single diagonal matrix
                            if issparse(Qblk{1})  % Sparse block diagonal matrix
                                Qblk = blkdiag(Qblk{:});
                            else               % Dense [row col var] matrix
                                % Pending ...
                            end
                        end
                        C = obj.data.C.(name);
                        l = obj.data.l.(name);
                        u = obj.data.u.(name);
                        if nargout > 4
                            vs = obj.data.vs.(name);
                            if nargout > 6
                                i1 = obj.idx.i1.(name);      %% starting row index
                                iN = obj.idx.iN.(name);      %% ending row index
                            end
                        end
                    else                                    %% indexing required
                        error('mp.sm_lin_constraint.params: set of quadratic constraints ''%s'' requires an IDX_LIST arg', name);
                    end
                else                            %% indexed named set
                    % (calls to substruct() are relatively expensive ...
                    % s = substruct('.', name, '{}', idx);
                    % ... so replace it with these more efficient lines)
                    sc = struct('type', {'.', '{}'}, 'subs', {name, idx});
                    Qblk  = subsref(obj.data.Q, sc);
                    C     = subsref(obj.data.C, sc);
                    l     = subsref(obj.data.l, sc);
                    u     = subsref(obj.data.u, sc);
                    if isblk
                        Qblk = blkdiag(Qblk{:});
                    end
                    if nargout > 3
                        vs = subsref(obj.data.vs, sc);
                        if nargout > 5
                            sn = sc; sn(2).type = '()';     %% num array field
                            i1 = subsref(obj.idx.i1, sn);   %% starting row index
                            iN = subsref(obj.idx.iN, sn);   %% ending row index
                        end
                    end
                end
            else                %% aggregate
                cache = obj.cache;
                if isempty(cache)
                    % Initialize parameters for aggregate
                    nx = var.N;              %% number of full set of variables
                    nquad = obj.N;           %% number of quadratic constraints
                    Qblk = cell(nquad, 1);   %% cell array of quadratic constraints
                    C = sparse(nquad, nx);   %% matrix of linear components of quadratic constraints
                    u = Inf(nquad, 1);       %% upper bound
                    l = -u;                  %% lower bound
                    
                    for j = 1:obj.NS   %% For each set of quadratic constraints
                        name = obj.order(j).name;
                        idx  = obj.order(j).idx;
                        [Qj, Cj, lj, uj, vsj, i1, iN] = obj.params(var, name, idx);
                        
                        if isempty(vsj)   % full nx vars
                            Qblk(i1:iN) = Qj;
                            C(i1:iN,:) = Cj;
                        else              % subset of vars
                            n = obj.get_N(name, idx); % number of constraints in current set
                            jj = var.varsets_idx(vsj); % indices for var set
                            rowid = repmat(jj', length(jj), 1);
                            colid = repmat(jj, length(jj), 1); colid = colid(:);
                            rowid = repmat(rowid, n, 1);
                            colid = repmat(colid, n, 1);
                            vals = cellfun(@(x)(x(:)), Qj, 'UniformOutput', false);
                            Qaux = mat2cell([rowid colid cell2mat(vals)], length(jj)^2*ones(n,1));
                            Qblk(i1:iN) = cellfun(@(x)(sparse(x(:,1), x(:,2), x(:,3), nx, nx)), Qaux, 'UniformOutput', false);
                            Caux = sparse(n, length(jj));
                            Caux(:, 1:length(jj)) = Cj;
                            C(i1:iN, jj) = Caux;
                        end
                        l(i1:iN,:) = lj;
                        u(i1:iN,:) = uj;
                    end
                    %% cache aggregated parameters
                    obj.cache = struct('C', C, 'l', l, 'u', u);
                    obj.cache.Qblk = Qblk;
                else
                    Qblk = cache.Qblk;
                    C = cache.C;
                    l = cache.l;
                    u = cache.u;
                end
                if nargout > 3
                    vs = {};
                end
            end
        end

        function [QFx_u, JQF, QFx, l_QFx] = eval(obj, var, x, name, idx)
            % Evaluate individual or full set of quadratic constraints.
            % ::
            %
            %   QFx_u = qcn.eval(var, x)
            %   QFx_u = qcn.eval(var, x, name)
            %   QFx_u = qcn.eval(var, x, name, idx)
            %   [QFx_u, JQF] = qcn.eval(...)
            %   [QFx_u, JQF, QFx] = qcn.eval(...)
            %   [QFx_u, JQF, QFx, l_QFx] = qcn.eval(...)
            %
            % For a given value of the variable vector x, this method evaluates
            % the quadratic constraints for an individual subset, if name or name 
            % and index list are provided, otherwise, for the full set of constraints.
            %
            % The constraints are of the form:
            %
            %      l_i <= (1/2)*x'*Q_i*x + c_i*x  <= u_i , for all i = 1,2,...,NQ
            %
            % or in compact form for a subset or full set of constraints:
            % 
            %                   l <= QFX <= u                       
            %
            % Returns QFX - u, and optionally l - QFX and the jacobian JQF.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   x (double) : full n x 1 variable vector x
            %   name (char array) : name of subset/block of quadratic constraints
            %       to evaluate
            %   idx (cell array) : *(optional)* index list for subset/block
            %       of quadratic constraints to evaluate (for an indexed subset)
            %
            % Outputs:
            %   QFx_u (double) : value of QFX - u 
            %   l_QFx (double) : *(optional)* value of l - QFX
            %   QFx (double) : *(optional)* value of QFX
            %   JQF (double) : *(optional)* Jacobian of quadratic constraints            
            % 
            % Examples::
            %
            %   QFx_u = qcn.eval(var, x)
            %   [QFx_u, l_QFx] = qcn.eval(var, x)
            %   [QFx_u, l_QFx, JQF] = qcn.eval(var, x)
            %   [QFx_u, l_QFx, JQF] = qcn.eval(var, x, 'my_set')
            %   [QFx_u, l_QFx, JQF] = qcn.eval(var, x, 'my_set', {3,2})
            %
            % See also add, params.

            if obj.N
                %% collect constraint parameters
                if nargin < 4                       %% full set
                    [Qblk_aux, C, l, u, vs] = obj.params(var);
                    Qblk = blkdiag(Qblk_aux{:});
                    Nq = obj.N; 
                elseif nargin < 5 || isempty(idx)   %% name, no idx provided
                    dims = size(obj.idx.i1.(name));
                    if prod(dims) == 1              %% simple named set
                        [Qblk, C, l, u, vs] = obj.params(var, name, {}, 1);
                        Nq = obj.get_N(name);
                    else
                        error('mp.sm_quad_constraint.eval: quadratic constraint set ''%s'' requires an IDX_LIST arg', name)
                    end
                else                                %% indexed named set
                    [Qblk, C, l, u, vs] = obj.params(var, name, idx, 1);
                    Nq = obj.get_N(name, idx);
                end

                %% assemble block diagonal matrix from x vector
                x = var.varsets_x(x, vs, 'vector');
                xx = mat2cell(repmat(sparse(x'), Nq, 1), ones(Nq,1));
                blkx = blkdiag(xx{:});

                %% Compute quadratic constraints
                if isempty(C)
                    QFX = 1/2 * diag(blkx * Qblk * blkx');
                    QFx_u = QFX - u;
                else
                    QFX = 1/2 * diag(blkx * Qblk * blkx') + C * x;
                    QFx_u = QFX - u;
                end

                if nargout > 1 %% Jacobian is requested
                    Qx = obj.blkprod2vertcat(blkx, Qblk, length(x));
                    JQF = Qx + C;
                    if nargout > 2
                        QFx = QFX;
                        if nargout > 3
                            l_QFx = l - QFX;
                        end
                    end
                end
            else
                QFx_u = [];
                if nargout > 1
                    JQF = [];
                    if nargout > 2
                        QFx = [];
                        if nargout > 3
                            l_QFx = [];
                        end
                    end
                end
            end
        end

        function obj = set_params(obj, var, name, idx, params, vals)
            % Modify quadratic constraint parameter data.
            % ::
            %
            %   qcn.set_params(var, name, params, vals)
            %   qcn.set_params(var, name, idx, params, vals)
            %
            % This method can be used to modify parameters for an existing
            % subset of quadratic constraints.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   name (char array) : name of subset/block of quadratic 
            %       constraints to modify
            %   idx (cell array) : *(optional)* index list for subset/block
            %       of quadratic constraints to modify (for an indexed subset)
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
            % Valid parameter names are ``Q``, ``C``, ``l``, ``u``, ``vs``.
            %
            % Examples::
            %
            %   qcn.set_params(var, 'y', {2,3}, {'l', 'u'}, {l, u});
            %   qcn.set_params(var, 'Pmis', 'all', {Q, C, l, u, vs});
            %
            % See also add, params.

            if nargin < 6
                vals = params;
                params = idx;
                idx = {};
            end

            %% create default list of parameters to update based on set type & inputs
            default_params = {'Q', 'C', 'l', 'u', 'vs'};

            %% standardize provided arguments in cell arrays params, vals
            [is_all, np, params, vals] = ...
                obj.set_params_std_args(default_params, params, vals);

            %% get current parameters
            [Q, C, l , u, vs] = obj.params(var, name, idx);
            MC0 = size(C, 1);
            MQ0 = size(Q, 1);
            if isempty(vs), vs = {vs}; end
            p = struct('C', C, 'l', l, 'u', u, 'vs', vs); p.Q = Q; %% current parameters
            u = struct('Q', 0, 'C', 0, 'l', 0, 'u', 0, 'vs',  0);  %% which ones to update

            %% replace with new parameters
            for j = 1:np
                p.(params{j}) = vals{j};
                u.(params{j}) = 1;
            end

            %% set missing default params for 'all'
            [MC, NC] = size(p.C);
            [MQ, NQ] = size(p.Q);
            if is_all
                u.Q = 1;            %% always update Q
                u.C = 1;            %% always update C
                u.l = 1;            %% always update l
                if np < 6
                    p.vs = {};
                    u.vs = 1;       %% update vs
                    if np < 5
                        p.u = Inf(MC, 1);
                        u.u = 1;    %% update u
                    end
                end
            end

            %% check consistency of parameters            
            %% Q must be a MQ x 1 cell array
            if MQ ~= MC || NQ ~= 1
                error('mp.sm_quad_constraint.set_params : dimension mismatch between rows/cols of ''C'' and ''Q''');
            end

            %% no dimension change unless 'all'
            if (MC ~= MC0 && ~is_all) || (MQ ~= MQ0 && ~is_all)
                error('mp.sm_quad_constraint.set_params: dimension change for ''%s'' not allowed except for ''all''', obj.nameidxstr(name, idx));
            end         

            %% check sizes of new values of l, u, and k
            for pn = {'l', 'u'}
                if u.(pn{1})
                    nn = length(p.(pn{1}));
                    if nn ~= MC
                        if nn == 0
                            switch pn{1}
                                case 'l'
                                    p.(pn{1}) = -Inf(MC, 0);
                                case 'u'
                                    p.(pn{1}) =  Inf(MC, 0);
                            end
                        elseif nn == 1
                            p.(pn{1}) = p.(pn{1}) * ones(MC, 1);   %% expand from scalar
                        else
                            error('mp.sm_quad_constraint.set_params: parameter ''%s'' ''%s'' should have length %d (or 1)', obj.nameidxstr(name, idx), pn{1}, NC);
                        end
                    end
                end
            end           
            
            %% check consistency of C and vs
            p.vs = mp.sm_variable.varsets_cell2struct(p.vs);
            nv = var.varsets_len(p.vs);     %% number of variables
            if u.C || u.vs
                if NC ~= nv
                    error('mp.sm_quad_constraint.set_params: for ''%s'' number of columns of ''C'' (%d) must be consistent with ''vs'' (%d)', obj.nameidxstr(name, idx), MC, nv);
                end
            end

            %% check consistency of Q
            if u.Q
                if iscell(p.Q)
                    sz = cellfun(@(x)(size(x)), p.Q, 'UniformOutput', false);
                else
                    error('mp.sm_quad_constraint.set_pamars: paramer ''Q'' must be a %d x 1 cell array', MQ);
                end
                if sum(sum(cell2mat(sz)))/(2*MQ) ~= nv
                    error('mp.sm_quad_constraint.set_params: for ''%s'' number of columns of ''Q'' (%d) must be consistent with ''vs'' (%d)', obj.nameidxstr(name, idx), MQ, nv);
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
            dMC = MC - MC0;            
            if is_all && dMC
                obj.set_params_update_dims(dMC, name, idx);
            end
        end

        function obj = display_soln(obj, var, soln, varargin)
            % Display solution values for quadratic constraints.
            % ::
            %
            %   qcn.display_soln(var, soln)
            %   qcn.display_soln(var, soln, name)
            %   qcn.display_soln(var, soln, name, idx)
            %   qcn.display_soln(var, soln, fid)
            %   qcn.display_soln(var, soln, fid, name)
            %   qcn.display_soln(var, soln, fid, name, idx)
            %
            % Displays the solution values for all quadratic constraints (default)
            % or an individual named or named/indexed subset.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   soln (struct) : full solution struct with these fields
            %       (among others):
            %
            %           - ``eflag`` - exit flag, 1 = success, 0 or negative =
            %             solver-specific failure code
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
            %   fid (fileID) : fileID of open file to write to (default is
            %       1 for standard output)
            %   name (char array) : *(optional)* name of individual subset
            %   idx (cell array) : *(optional)* indices of individual subset

            [fid, name, idx, idxs, hdr1] = obj.display_soln_std_args(varargin{:});

            if obj.N
                [Q, C, vl, vu] = obj.params(var);
                Nq = length(vl);
                Qblk = blkdiag(Q{:});
                xx = mat2cell(repmat(sparse(soln.x'), Nq, 1), ones(Nq,1));
                blkx = blkdiag(xx{:});
                v = diag(1/2 * blkx * Qblk * blkx' + C * soln.x);
                if isempty(soln.lambda)
                    mu_l_quad = NaN(size(v));
                    mu_u_quad = mu_l_quad;
                else
                    mu_l_quad = soln.lambda.mu_l_quad;
                    mu_u_quad = soln.lambda.mu_u_quad;
                end

                %% print header rows
                hdr2 = {'   mu_lb     lb       val      ub      mu_ub', ...
                        ' -------- -------- -------- -------- --------' };
                obj.display_soln_print_headers(fid, hdr1, hdr2);

                %% print data
                none = '- ';
                for k = 1:length(idxs)
                    obj.display_soln_print_row(fid, idxs(k));

                    if isnan(mu_l_quad(idxs(k)))
                        mu_lb = sprintf( ' ');
                    elseif abs(mu_l_quad(idxs(k))) < obj.mu_thresh()
                        mu_lb = sprintf(none);
                    else
                        mu_lb = obj.sprintf_num(8, mu_l_quad(idxs(k)));
                    end
                    if isnan(mu_u_quad(idxs(k)))
                        mu_ub = sprintf( ' ');
                    elseif abs(mu_u_quad(idxs(k))) < obj.mu_thresh()
                        mu_ub = sprintf(none);
                    else
                        mu_ub = obj.sprintf_num(8, mu_u_quad(idxs(k)));
                    end
                    if vl(idxs(k)) < -obj.num_inf()
                        lb = sprintf(none);
                    else
                        lb = obj.sprintf_num(8, vl(idxs(k)));
                    end
                    if vu(idxs(k)) > obj.num_inf()
                        ub = sprintf(none);
                    else
                        ub = obj.sprintf_num(8, vu(idxs(k)));
                    end
                    fprintf(fid, '%9s%9s%9s%9s%9s\n', ...
                        mu_lb, lb, obj.sprintf_num(8, v(idxs(k))), ub, mu_ub);
                end

                %% print footer rows
                fprintf(fid, '%s\n', [hdr1{2} hdr2{2}]);
                fprintf(fid, '%7s %-28s%9s%9s%9s%9s%9s\n', '', 'Min', ...
                    obj.sprintf_num(8, min(mu_l_quad)), obj.sprintf_num(8, min(vl)), ...
                    obj.sprintf_num(8, min(v)), ...
                    obj.sprintf_num(8, min(vu)), obj.sprintf_num(8, min(mu_u_quad)));
                fprintf(fid, '%7s %-28s%9s%9s%9s%9s%9s\n', '', 'Max', ...
                    obj.sprintf_num(8, max(mu_l_quad)), obj.sprintf_num(8, max(vl)), ...
                    obj.sprintf_num(8, max(v)), ...
                    obj.sprintf_num(8, max(vu)), obj.sprintf_num(8, max(mu_u_quad)));
                fprintf(fid, '\n');
            end
        end

        function varargout = get_soln(obj, var, soln, varargin)
            % Fetch solution values for specific named/indexed subsets.
            % ::
            %
            %   vals = qcn.get_soln(var, soln, name)
            %   vals = qcn.get_soln(var, soln, name, idx)
            %   vals = qcn.get_soln(var, soln, tags, name)
            %   vals = qcn.get_soln(var, soln, tags, name, idx)
            %
            % Returns named/indexed quadratic constraint results for a solved
            % model, evaluated at the solution found.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   soln (struct) : full solution struct with these fields
            %       (among others):
            %
            %           - ``eflag`` - exit flag, 1 = success, 0 or negative =
            %             solver-specific failure code
            %           - ``x`` - variable values
            %           - ``lambda`` - constraint shadow prices, struct with
            %             fields:
            %
            %               - ``eqnonlin`` - nonlinear equality constraints
            %               - ``ineqnonlin`` - nonlinear inequality constraints
            %               - ``mu_l`` - linear constraint lower bounds
            %               - ``mu_u`` - linear constraint upper bounds
            %               - ``mu_l_quad`` - quadratic constraint lower bounds
            %               - ``mu_u_quad`` - quadratic constraint upper bounds
            %               - ``lower`` - variable lower bounds
            %               - ``upper`` - variable upper bounds
            %   tags (char array or cell array of char arrays) : names of
            %       desired outputs, default is ``{'g', 'mu_l', 'mu_u'}``
            %       with valid values:
            %
            %           - ``'g'`` - 2 element cell array with constraint values
            %             QFx - u and l - QFx, respectively, where QFx is
            %             the quadratic form related to the quadratic
            %             constraints being fetched
            %           - ``'QFx_u'`` or ``'f'`` - constraint values QFx - u
            %           - ``'l_QFx'`` - constraint values l - QFx
            %           - ``'mu_l'`` - shadow price on l - QFx
            %           - ``'mu_u'`` - shadow price on QFx - u
            %   name (char array) : name of the subset
            %   idx (cell array) : *(optional)* indices of the subset
            %
            % Outputs:
            %     : Variable number of outputs corresponding to ``tags`` input.
            %       If ``tags`` is empty or not specified, the calling context
            %       will define the number of outputs, returned in order of
            %       default tags.
            %
            % Example::
            %
            %     [g, mu_l, mu_u] = qcn.get_soln(var, soln, 'flow');
            %     mu_l_quad_Pmis_5_3 = qcn.get_soln(var, soln, 'mu_l', 'Pmis', {5,3});
            %
            % For a complete set of solution values, using the parse_soln()
            % method may be more efficient.
            %
            % See also parse_soln.

            %% input arg handling
            [tags, name, idx, N, i1, iN] = obj.get_soln_std_args(varargin{:});

            %% get outputs
            varargout = cell(1, nargout);
            if N && ~isempty(soln.eflag)
                if any(ismember({'g', 'QFx_u', 'l_QFx'}, tags(1:nargout)))
                    g = cell(1,4);
                    [g{:}] = obj.eval(var, soln.x, name, idx);
                end
                for k = 1:nargout
                    switch tags{k}
                        case 'g'
                            varargout{k} = g;
                        case 'QFx_u'
                            varargout{k} = g{1};
                        case 'l_QFx'
                            varargout{k} = g{4};
                        case 'f'
                            varargout{k} = soln.f(i1:iN);
                        case 'mu_l'
                            varargout{k} = soln.lambda.mu_l_quad(i1:iN);
                        case 'mu_u'
                            varargout{k} = soln.lambda.mu_u_quad(i1:iN);
                        otherwise
                            error('mp.sm_quad_constraint.get_soln: unknown tag ''%s''', tags{k});
                    end
                end
            end
        end

        function ps = parse_soln(obj, soln, stash)
            % Parse solution for quadratic constraints.
            % ::
            %
            %   ps = qcn.parse_soln(soln)
            %
            % Parse a full solution struct into parts corresponding to
            % individual quadratic constraint subsets.
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
            %               - ``mu_l_quad`` - quadratic constraint lower bounds
            %               - ``mu_u_quad`` - quadratic constraint upper bounds
            %               - ``lower`` - variable lower bounds
            %               - ``upper`` - variable upper bounds
            %   stash (boolean) : if true, store return value in :attr:`soln`
            %       property
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
                if isfield(soln.lambda, 'mu_l_quad')
                    if isfield(soln.lambda, 'mu_u_quad')
                        params = struct('src', {soln.lambda.mu_l_quad, soln.lambda.mu_u_quad}, ...
                                        'dst', {'mu_l', 'mu_u'});
                    else
                        params = struct('src', soln.lambda.mu_l_quad, 'dst', 'mu_l');
                    end
                else
                    if isfield(soln.lambda, 'mu_u_quad')
                        params = struct('src', soln.lambda.mu_u_quad, 'dst', 'mu_u');
                    else
                        params = [];
                    end
                end
                if ~isempty(params)
                    ps = obj.parse_soln_fields(params);
                end
            end

            if nargin > 2 && stash
                obj.soln = ps;
            end
        end
    end     %% methods

    methods (Static)
        function M = blkprod2vertcat(blk1, blk2, n)
            % Compute all the products of two block diagonal matrices
            % ::
            %
            %   M = qcn.blkprod2vertcat(blk1, blk2, n)
            %
            % Take two block diagonal matrices and returns a matrix formed 
            % by stacking vertically the result of all products of the two.
            % Here the matrices are assumed to have compatible sizes, that
            % is, the number of columns of the set of matrices used to form
            % the first block diagonal matrix must be the same as the nuber
            % of rows of the matrices used to for the second block diagonal
            % matrix.
            %
            % Inputs:
            %   blk1 : block diagonal matrix formed from a set of matrices, 
            %          each of size m1 x n columns
            %   blk2 : block diagonal matrix formed from a set of matrices,
            %          each of size n x m2
            %      n : compatible dimension of the two block diagonal matrices
            %          used to find the number of blocks to be multiplied
            %
            % Outputs:
            %   M : (m1*n)x(m2) matrix holding a vertical stack of the
            %       resulting products between the two block diagonal
            %       matrices
            %
            % Examples:
            % 
            %   xx = blkdiag(repmat(x, m, 1));  %% x \in R^n, xx \in R^(mxn)
            %   QQ = blkdiag(Q{:});  %% Q is a cell array of m matrices \in R^(nxn)   
            %   M = qcn.blkprod2vertcat(xx, QQ, n) %% M is a n x n matrix
            %
            % See also eval.

            [row1, col1] = size(blk1);
            [row2, col2] = size(blk2);
            
            if (col1) ~= (row2)
                error('sm_quad_constraint.blkprod2vertcat: number of columns of elements in blk1 (%d) do not match the number of rows of elements in blk2 (%d) \n', col1/n, row2/n);
            end
            
            blkprod = blk1 * blk2;    % Product of blocks
            N = col1/n;               % Number of blocks in blkprod
            
            % id_block = sparse(logical(ones(row1/N, col2/N)));
            % id_blkprod = mat2cell(repmat(id_block, N, 1), (row1/N)*ones(N,1));
            % id_blkprod = blkdiag(id_blkprod{:});
            % 
            % M = blkprod(id_blkprod);
            % M = reshape(M, N, [])';

            %                 rows                  cols
            id_start = [(0:N-1)*(row1/N)+1 ; (0:N-1)*(col2/N)+1 ];
            id_end   = [   (1:N)*(row1/N)  ;    (1:N)*(col2/N)  ];

            id_start = mat2cell(id_start(:), 2*ones(N,1));
            id_end   = mat2cell(id_end(:), 2*ones(N,1));

            M = cellfun(@(x,y)(blkprod(x(1):y(1), x(2):y(2))), id_start, id_end, 'UniformOutput', false);

            M = cell2mat(M);
        end
    end

    methods (Access=protected)
        function default_tags = get_soln_default_tags(obj)
            % Return default tags for get_soln().
            % ::
            %
            %   default_tags = nln.get_soln_default_tags()
            %
            % Output:
            %   default_tags (cell array) : tags defining the default outputs
            %       of get_soln(), namely ``{'g', 'lam', 'dg'}``
            %
            % See also get_soln.

            default_tags = {'g', 'mu_l', 'mu_u'};
        end
    end     %% methods (Access=protected)
end         %% classdef