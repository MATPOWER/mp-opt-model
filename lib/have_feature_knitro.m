function [TorF, vstr, rdate] = have_feature_knitro()
%HAVE_FEATURE_KNITRO  Detect availability/version info for Artelys Knitro
%
%   Used by HAVE_FEATURE.

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

tmp = have_feature('knitromatlab', 'all');
if tmp.av
    TorF = tmp.av;
    vstr = tmp.vstr;
    rdate = tmp.date;
else
    tmp = have_feature('ktrlink', 'all');
    if tmp.av
        TorF = tmp.av;
        vstr = tmp.vstr;
        rdate = tmp.date;
    else
        TorF = 0;
        vstr = '';
        rdate = '';
    end
end
