function mosek_log_callback(h, str)
% mosek_log_callback  Used to redirect console output via mp.logger.
% ::
%
%   mosek_log_callback(handle, str)
%
% Used as a log callback in miqps_mosek and qps_mosek to redirect MOSEK
% console output via mp.logger.
%
% See also mp.logger, miqps_mosek, qps_mosek.

%   MP-Opt-Model
%   Copyright (c) 2026, Ray Zimmerman
%   by Ray Zimmerman
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

mp_printf(str);
