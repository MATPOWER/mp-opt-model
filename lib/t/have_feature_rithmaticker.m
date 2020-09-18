function [TorF, vstr, rdate] = have_feature_rithmaticker()
%HAVE_FEATURE_RITHMATICKER  Detect availability/version info for rithmaticker
%
%   Used by HAVE_FEATURE.

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

TorF = exist('rithmaticker', 'file') == 2;
if TorF
    vstr = '3.1.4';
    rdate = datestr([2019 5 30 0 0 0], 'dd-mmm-yyyy');
else
    vstr = '';
    rdate = '';
end
