function [TorF, vstr, rdate] = have_feature_clp()
%HAVE_FEATURE_CLP  Detect availability/version info for CLP
%
%   Used by HAVE_FEATURE.

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

tmp = have_feature('opti_clp', 'all');
vstr = '';
rdate = '';
if tmp.av   %% have opti_clp
    TorF = tmp.av;
    vstr = tmp.vstr;
    rdate = tmp.date;
elseif exist('clp','file') == 2 && exist('mexclp','file') == 3
    TorF = 1;
else
    TorF = 0;
end
