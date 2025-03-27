function [ieq, igt, ilt, Qe, Ce, be, Qi, Ci, bi] = convert_quad_constraint(Q, C, l, u)
% convert_quad_constraint - Convert a set of bounded quadratic constraints to 
% equality/inequality pair
% ::
%
%   [ieq, igt, ilt, Qe, Ce, be, Qi, Ci, bi] = convert_quad_constraint(Q, C, l, u)
%   [ieq, igt, ilt, Q, C, b] = convert_quad_constraint(Q, C, l, u)
%
% Convert a set of constraints of the form:
%
%       l(j) <= x'*Q{j}*x + C(j,:)*x <= u(j),  j = 1,2,...,NQ
%
% to:
%
%       x'*Qe{j}*x + Ce(j,:)*x  =  be(j),  j = 1,2,...,NQe
%       x'*Qi{j}*x + Ci(j,:)*x  <= bi(j),  j = 1,2,...,NQi   , NQ = NQe+NQi
%
% where:
%
%       Qe = Q(ieq)
%       Ce = C(ieq,:)
%       be = u(ieq)
%       Qi = [Q(ilt); -Q(igt)]
%       Ci = [C(ilt,:); -C(igt,:)]
%       bi = [u(ilt);  -l(igt)]
%
% Alternatively, the returned cells, matrices, and RHS vectors can be stacked
% into a single set with the equalities first, then the inequalities.
%
%       Q = [Qe; Qi]
%       C = [Ce; Ci]
%       b = [be; bi]
% 
% Inputs:
%   Q (double) : NQ x 1, cell array with sparse n x n symmetric matrices 
%                holding the quadratic parameters of the quadratic constraints
%   C (double) : NQ x n, matrix whose rows represent the linear parameters of
%                the quadratic constraints
%   l (double) : NQ x 1, quadratic constraint lower bound vector
%   u (double) : NQ x 1, quadratic constraint upper bound vector
%
% Outputs:
%   ieq (integer) : vector of indices of equality constraints
%   igt (integer) : vector of indices of lower bounded inequality constraints
%   ilt (integer) : vector of indices of upper bounded inequality constraints
%   Qe (double)   : cell array of quadratic parameters for equality constraints
%   Ae (double)   : matrix of linear parameters for equality constraints
%   be (double)   : equality constraint RHS
%   Qi (double)   : cell array of quadratic parameters for inequality constraints
%   Ai (double)   : matrix of linear terms for inequality constraints
%   bi (double)   : inequality constraint RHS
%   Q (double)    : stacked cell of quadratic parameters
%   C (double)    : stacked constraint matrix of linear parameters
%   b (double)    : stacked RHS
%
% See also convert_quad_constraint_multipliers.

%   MP-Opt-Model
%   Copyright (c) 2019-2023, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia
%
%   This file is part of MP-Opt-Model..
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

eq = abs(u-l) <= eps;
ieq = find( eq );               %% equality
igt = find( l > -1e10 & ~eq );  %% inequality w/finite lower bound
ilt = find( u <  1e10 & ~eq );  %% inequality w/finite upper bound

if nargout < 9      %% stacked case
    Qe = [Q([ieq; ilt]); cellfun(@(x)(-1*x), Q(igt), 'UniformOutput', false)];
    Ce = [C([ieq; ilt],:); -C(igt,:)];
    be = [u([ieq; ilt]);  -l(igt)];
else                %% equality/inequality pairs
    Qe = Q(ieq);
    Ce = C(ieq,:);
    be = u(ieq);
    Qi = [Q(ilt); cellfun(@(x)(-1*x), Q(igt), 'UniformOutput', false)];
    Ci = [C(ilt,:); -C(igt,:)];
    bi = [u(ilt);  -l(igt)];
end