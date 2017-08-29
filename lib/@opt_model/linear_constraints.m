function [A, l, u] = linear_constraints(om)
%LINEAR_CONSTRAINTS  Builds and returns the full set of linear constraints.
%   [A, L, U] = OM.LINEAR_CONSTRAINTS()
%   Builds the full set of linear constraints based on those added by
%   ADD_LIN_CONSTRAINTS.
%
%       L <= A * x <= U
%
%   Example:
%       [A, l, u] = om.linear_constraints();
%
%   See also OPT_MODEL, ADD_LIN_CONSTRAINTS.

%   MATPOWER
%   Copyright (c) 2008-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.


%% initialize A, l and u
nnzA = 0;
s = struct('type', {'.', '{}'}, 'subs', {'', 1});
for k = 1:om.lin.NS
    name = om.lin.order(k).name;
    idx  = om.lin.order(k).idx;
    if isempty(idx)
        nnzA = nnzA + nnz(om.lin.data.A.(name));
    else
        % (calls to substruct() are relatively expensive ...
        % s = substruct('.', name, '{}', idx);
        % ... so replace it with these more efficient lines)
        s(1).subs = name;
        s(2).subs = idx;
        nnzA = nnzA + nnz(subsref(om.lin.data.A, s));
    end
end
At = sparse([], [], [], om.var.N, om.lin.N, nnzA);  %% use A transpose for speed
u = Inf(om.lin.N, 1);
l = -u;

%% fill in each piece
s2 = s;
s(2).type = '()';
s1 = s;
for k = 1:om.lin.NS
    name = om.lin.order(k).name;
    idx  = om.lin.order(k).idx;
    if isempty(idx)
        N = om.lin.idx.N.(name);
    else
        % (calls to substruct() are relatively expensive ...
        % s1 = substruct('.', name, '()', idx);
        % s2 = substruct('.', name, '{}', idx);
        % ... so replace them with these more efficient lines)
        s1(1).subs = name;
        s1(2).subs = idx;
        s2(1).subs = name;
        s2(2).subs = idx;
        N = subsref(om.lin.idx.N, s1);
    end
    if N                                %% non-zero number of rows to add
        if isempty(idx)
            Ak = om.lin.data.A.(name);          %% A for kth linear constrain set
            i1 = om.lin.idx.i1.(name);          %% starting row index
            iN = om.lin.idx.iN.(name);          %% ending row index
            vs = om.lin.data.vs.(name);         %% var sets
        else
            Ak = subsref(om.lin.data.A, s2);    %% A for kth linear constrain set
            i1 = subsref(om.lin.idx.i1, s1);    %% starting row index
            iN = subsref(om.lin.idx.iN, s1);    %% ending row index
            vs = subsref(om.lin.data.vs, s2);   %% var sets
        end
        if isempty(vs)          %% full rows
            if size(Ak,2) == om.var.N
                At(:, i1:iN) = Ak';     %% assign as columns in transpose for speed
            else                %% must have added vars since adding
                                %% this constraint set
                At(1:size(Ak,2), i1:iN) = Ak';  %% assign as columns in transpose for speed
            end
        else                    %% selected columns
            jj = om.varsets_idx(vs);    %% column indices for var set
            Ai = sparse(N, om.var.N);
            Ai(:, jj) = Ak;
            At(:, i1:iN) = Ai';     %% assign as columns in transpose for speed
        end

        if isempty(idx)
            l(i1:iN) = om.lin.data.l.(name);
            u(i1:iN) = om.lin.data.u.(name);
        else
            l(i1:iN) = subsref(om.lin.data.l, s2);
            u(i1:iN) = subsref(om.lin.data.u, s2);
        end
    end
end
A = At';
