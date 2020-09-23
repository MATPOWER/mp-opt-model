function [TorF, vstr, rdate] = have_feature_ipopt()
%HAVE_FEATURE_IPOPT  Detect availability/version info for IPOPT
%
%   Feature detection function implementing 'ipopt' tag for HAVE_FEATURE
%   to detect availability/version of IPOPT, a nonlinear programming
%   solver from COIN-OR (https://github.com/coin-or/Ipopt).
%
%   See also HAVE_FEATURE, NLPS_MASTER, QPS_MASTER, IPOPT.

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

TorF = exist('ipopt', 'file') == 3;
vstr = '';
rdate = '';
if TorF
    if have_feature('evalc')
        str = evalc('qps_ipopt([],[1; 1],[1 1],[2],[2],[1; 1],[1; 1],[1; 1],struct(''verbose'', 2))');
        pat = 'Ipopt version ([^\s,]+)';
        [s,e,tE,m,t] = regexp(str, pat);
        if ~isempty(t)
            vstr = t{1}{1};
            if vstr2num(vstr) >= 3.011 && ~exist('ipopt_auxdata', 'file')
                TorF = 0;
                warning('Improper installation of IPOPT. Version %s detected, but IPOPT_AUXDATA.M is missing.', vstr);
            end
        end
    else
        try
            x = feval('qps_ipopt', [],[1; 1],[1 1],[2],[2],[1; 1],[1; 1],[1; 1],struct('verbose', 0));
            if ~isequal(x, [1;1])
                TorF = 0;
            end
        catch
            TorF = 0;
        end
    end
end

function num = vstr2num(vstr)
% Converts version string to numerical value suitable for < or > comparisons
% E.g. '3.11.4' -->  3.011004
pat = '\.?(\d+)';
[s,e,tE,m,t] = regexp(vstr, pat);
b = 1;
num = 0;
for k = 1:length(t)
    num = num + b * str2num(t{k}{1});
    b = b / 1000;
end
