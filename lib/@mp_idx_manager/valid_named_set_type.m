function str = valid_named_set_type(obj, set_type)
%VALID_NAMED_SET_TYPE  Returns a label for the given named set type or empty.
%
%   -----  PRIVATE METHOD  -----
%
%   STR = OBJ.VALID_NAMED_SET_TYPE(SET_TYPE)
%
%   Returns a string corresponding to the type of the named set for
%   valid values of SET_TYPE, otherwise an empty string. Can be used to
%   check whether SET_TYPE is a valid type.
%
%   Valid types and their corresponding return values are:
%       SET_TYPE    STR
%       'var'     'variable'
%       'lin'     'linear constraint'
%       'nle'     'nonlinear equality constraint'
%       'nli'     'nonlinear inequality constraint'

%   MATPOWER
%   Copyright (c) 2017-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if isfield(obj.set_types, set_type)
    str = obj.set_types.(set_type);
else
    str = '';
end
