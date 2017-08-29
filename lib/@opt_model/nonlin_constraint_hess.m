function d2G = nonlin_constraint_hess(om, x, lam, iseq)
%NONLIN_CONSTRAINT_HESS  Builds and returns Hessian of nonlinear constraints.
%   D2G = OM.NONLIN_CONSTRAINT_HESS(X, LAM, ISEQ)
%   Builds the Hessian of the full set of nonlinear equality or inequality
%   constraints for given values of the optimization vector X and dual
%   variables LAM, based on constraints added by ADD_NLN_CONSTRAINTS.
%
%       g(X) = 0
%       h(X) <= 0
%
%   Example:
%       d2G = om.nonlin_constraint_hess(x, lam, 1)
%       d2H = om.nonlin_constraint_hess(x, lam, 0)
%
%   See also OPT_MODEL, ADD_NLN_CONSTRAINTS, NONLIN_CONSTRAINTS.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% get constraint type
if iseq         %% equality constraints
    om_nlx = om.nle;
else            %% inequality constraints
    om_nlx = om.nli;
end

%% initialize d2G (use transpose for speed on older versions of Matlab)
d2Gt = sparse(om.var.N, om.var.N);

%% fill in each piece
s = struct('type', {'.', '()'}, 'subs', {'', 1});
s1 = s;
s2 = s;
s2(2).type = '{}';
for k = 1:om_nlx.NS
    name = om_nlx.order(k).name;
    idx  = om_nlx.order(k).idx;
    if isempty(idx)
        if ~isfield(om_nlx.data.hess, name)
            continue;   %% skip, there is no function handle stored here,
                        %% the function value for this named set was included
                        %% in the value computed by a previous named set
        end
        N = om_nlx.idx.N.(name);    %% number of constraint functions
                                    %% evaluated for this named set
        if isfield(om_nlx.data.include, name)
            N = N + sum(om_nlx.data.include.(name).N);
        end
    else
        % (calls to substruct() are relatively expensive ...
        % s1 = substruct('.', name, '()', idx);
        % s2 = substruct('.', name, '{}', idx);
        % ... so replace them with these more efficient lines)
        s1(1).subs = name;
        s1(2).subs = idx;
        s2(1).subs = name;
        s2(2).subs = idx;
        N = subsref(om_nlx.idx.N, s1);
    end
    if N                                %% non-zero number of rows
        if isempty(idx)
            d2G_fcn = om_nlx.data.hess.(name);  %% Hessian fcn for kth constraint set
            i1 = om_nlx.idx.i1.(name);          %% starting row index
            iN = i1 + N - 1;                    %% ending row index
            vs = om_nlx.data.vs.(name);         %% var sets
        else
            d2G_fcn = subsref(om_nlx.data.hess, s2);  %% Hessian fcn for kth constraint set
            i1 = subsref(om_nlx.idx.i1, s1);    %% starting row index
            iN = subsref(om_nlx.idx.iN, s1);    %% ending row index
            vs = subsref(om_nlx.data.vs, s2);   %% var sets
        end
        xx = om.varsets_x(x, vs);
        d2Gk = d2G_fcn(xx, lam(i1:iN));     %% evaluate kth Hessian
        
        nk = size(d2Gk, 2);
        if isempty(vs)          %% all rows of x
            if nk == om.var.N
                d2Gkt_full = d2Gk';
            else                %% must have added vars since adding
                                %% this constraint set
                d2Gk_all_cols = sparse(nk, om.var.N);
                d2Gk_all_cols(:, 1:nk) = d2Gk;
                d2Gkt_full = sparse(om.var.N, om.var.N);
                d2Gkt_full(:, 1:nk) = d2Gk_all_cols';
            end
        else                    %% selected rows of x
            jj = om.varsets_idx(vs);    %% indices for var set
            d2Gk_all_cols = sparse(nk, om.var.N);
            d2Gk_all_cols(:, jj) = d2Gk;
            d2Gkt_full = sparse(om.var.N, om.var.N);
            d2Gkt_full(:, jj) = d2Gk_all_cols';
        end
        d2Gt = d2Gt + d2Gkt_full;
    end
end
d2G = d2Gt';
