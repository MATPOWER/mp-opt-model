function t_have_fcn_legacy(quiet)
%T_HAVE_FCN_LEGACY  Tests for HAVE_FCN.

%   MP-Opt-Model
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin < 1
    quiet = 0;
end

n_tests = 20;

t_begin(n_tests, quiet);

%% save cache
saved = have_fcn_legacy('all', 'get_cache');

%% test set_cache/get_cache
t = 'set_cache/get_cache';
e = struct('foo', 2, 'bar', 3, 'baz', 'boz');
have_fcn_legacy(e, 'set_cache');
t_ok(isequal(have_fcn_legacy('all', 'get_cache'), e), t);
have_fcn_legacy(64, 'set_cache');
t_ok(isequal(have_fcn_legacy('all', 'get_cache'), 64), t);

%% clear cache
have_fcn_legacy(struct(), 'set_cache');

%% Matlab vs. Octave
if exist('OCTAVE_VERSION', 'builtin') == 5
    t_ok(have_fcn_legacy('octave'), 'Octave');
    t_ok(~have_fcn_legacy('matlab'), 'not Matlab');
else
    t_ok(~have_fcn_legacy('octave'), 'not Octave');
    t_ok(have_fcn_legacy('matlab'), 'Matlab');
end

t = '<MPOM>/lib/t/t_have_fcn must not be in path';
t_ok(exist('rithmaticker') ~= 2, t);

%% find path to this test file
cwd = pwd;
[p, n, e] = fileparts(which('t_have_fcn_legacy'));

%% initially not available
t = 'have_fcn_legacy(''rithmaticker'')';
t_ok(have_fcn_legacy('rithmaticker') == 0, [t ' : not available']);

%% switch dir so it is available
cd(fullfile(p, 't_have_fcn'));
cwd2 = pwd;

t = '<MPOM>/lib/t/t_have_fcn must be in path';
t_ok(exist('rithmaticker') == 2, t);

%% still not available (cached)
t = 'have_fcn_legacy(''rithmaticker'')';
t_ok(have_fcn_legacy('rithmaticker') == 0, [t ' : still not available (cached)']);

%% clear cache, check again
have_fcn_legacy('rithmaticker', 'clear_cache');
t_ok(have_fcn_legacy('rithmaticker') == 1, [t ' : available (cache cleared)']);

cd(cwd);
t = 'successful switch to original working directory';
t_ok(strcmp(cwd, pwd), t);

t = 'have_fcn_legacy(''rithmaticker'')';
t_ok(have_fcn_legacy('rithmaticker') == 1, [t ' : still available (cached)']);

t = 'have_fcn_legacy(''rithmaticker'', ''vstr'')';
t_ok(strcmp(have_fcn_legacy('rithmaticker', 'vstr'), '3.1.4'), t);

t = 'have_fcn_legacy(''rithmaticker'', ''vnum'')';
t_is(have_fcn_legacy('rithmaticker', 'vnum'), 3.001004, 12, t);

t = 'have_fcn_legacy(''rithmaticker'', ''date'')';
t_ok(strcmp(have_fcn_legacy('rithmaticker', 'date'), '30-May-2019'), t);

t = 'have_fcn_legacy(''rithmaticker'', ''all'')';
rv = have_fcn_legacy('rithmaticker', 'all');
t_ok(isstruct(rv), [t ' : isstruct']);
t_is(rv.av, 1, 12, [t ' : av']);
t_ok(strcmp(rv.vstr, '3.1.4'), [t ' : vstr']);
t_is(rv.vnum, 3.001004, 12, [t ' : vnum']);
t_ok(strcmp(rv.date, '30-May-2019'), [t ' : date']);

%% clear cache, check again
t = 'have_fcn_legacy(''rithmaticker'')';
have_fcn_legacy('all', 'clear_cache');
t_ok(have_fcn_legacy('rithmaticker') == 0, [t ' : not available (cache cleared)']);

%% restore cache
have_fcn_legacy(saved, 'set_cache');

t_end;
