function [TorF, vstr, rdate] = have_feature_evalc()
%HAVE_FEATURE_EVALC  Detect availability/version info for EVALC
%
%   Used by HAVE_FEATURE.

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% ignore evalc in Octave (I don't remember why)
v = have_feature('matlab', 'all');
if v.av
    TorF = 1;
    vstr = v.vstr;
    rdate = v.date;
else
    TorF = 0;
    vstr = '';
    rdate = '';
end
