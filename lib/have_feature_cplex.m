function [TorF, vstr, rdate] = have_feature_cplex()
%HAVE_FEATURE_CPLEX  Detect availability/version info for CPLEX
%
%   Feature detection function implementing 'cplex' tag for HAVE_FEATURE
%   to detect availability/version of CPLEX (IBM ILOG CPLEX Optimizer).
%
%   See also HAVE_FEATURE, QPS_MASTER, MIQPS_MASTER, CPLEX, CPLEXQP, CPLEXLP.

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

TorF = 0;
vstr = '';
rdate = '';
if exist('cplexqp', 'file')
    %% it's installed, but we need to check for MEX for this arch
    p = which('cplexqp');   %% get the path
    len = length(p) - length('cplexqp.p');
    w = what(p(1:len));             %% look for mex files on the path
    for k = 1:length(w.mex)
        if regexp(w.mex{k}, 'cplexlink[^\.]*');
            TorF = 1;
            break;
        end
    end
end
if TorF
    try
        cplex = Cplex('null');
        vstr = cplex.getVersion;
    catch
        TorF = 0;
    end
end
