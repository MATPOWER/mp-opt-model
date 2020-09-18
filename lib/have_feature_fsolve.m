function [TorF, vstr, rdate] = have_feature_fsolve()
%HAVE_FEATURE_FSOLVE  Detect availability/version info for FSOLVE
%
%   Used by HAVE_FEATURE.

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

TorF = exist('fsolve', 'file') == 2;
if TorF
    if have_feature('matlab')
        v = have_feature('optim', 'all');
        vstr = v.vstr;
        rdate = v.date;
    else
        v = have_feature('octave', 'all');
        vstr = v.vstr;
        rdate = v.date;
    end
else
    vstr = '';
    rdate = '';
end
