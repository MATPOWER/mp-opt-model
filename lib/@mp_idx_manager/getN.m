function N = getN(obj, set_type, name, idx)
%GETN  Returns the number of variables, constraints or cost rows.
%   N = OBJ.GETN(SET_TYPE)
%   N = OBJ.GETN(SET_TYPE, NAME)
%   N = OBJ.GETN(SET_TYPE, NAME, IDX)
%
%   Returns either the total number of variables/constraints/cost rows
%   or the number corresponding to a specified named block.
%
%   Examples:
%       N = obj.getN('var')     : total number of variables
%       N = obj.getN('lin')     : total number of linear constraints
%       N = obj.getN('nle')     : total number of nonlin equality constraints
%       N = obj.getN('nli')     : total number of nonlin inequality constraints
%       N = obj.getN('qdc')     : total number of quadratic cost rows
%       N = obj.getN('nlc')     : total number of general nonlinear cost rows
%       N = obj.getN('var', name)   : # of variables in named set
%       N = obj.getN('lin', name)   : # of linear constraints in named set
%       N = obj.getN('nle', name)   : # of nonlin eq cons. in named set
%       N = obj.getN('nli', name)   : # of nonlin ineq cons. in named set
%       N = obj.getN('qdc', name)   : # of quadratic cost rows in named set
%       N = obj.getN('nlc', name)   : # of gen nonlin cost rows in named set
%       N = obj.getN('var', name, idx)  : # of variables in indexed named set
%
%   See also OPT_MODEL.

%   MATPOWER
%   Copyright (c) 2008-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 3
    N = obj.(set_type).N;
else
    if isfield(obj.(set_type).idx.N, name)
        if nargin < 4
            idx = {};
        end
        s1 = substruct('.', name, '()', idx);
        N = subsref(obj.(set_type).idx.N, s1);
    else
        N = 0;
    end
end
