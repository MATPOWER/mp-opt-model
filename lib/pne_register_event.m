function event_list = pne_register_event(event_list, name, fcn, tol, locate)
%PNE_REGISTER_EVENT  Register event functions=
%   EVENT_LIST = PNE_REGISTER_EVENT(EVENT_LIST, NAME, FCN, TOL, LOCATE)
%
%   Registers an event function to be called by PNES_MASTER.
%
%   Inputs:
%       EVENT_LIST : struct array containing info about registered event fcns
%           with fields name, fcn, tol, and locate corresponding to the
%           respective inputs below
%       NAME : string containing event name
%       FCN : string containing name of event function, returning numerical
%             scalar or vector value that changes sign at location of the event
%       TOL : scalar or vector of same dimension as event function return value
%             of tolerance for detecting the event, i.e. abs(val) <= tol
%       LOCATE : flag indicating whether the event requests a rollback step
%                to locate the event function zero
%
%   Outputs:
%       EVENT_LIST : updated struct containing info about registered event fcns
%                    with new entry appended at end

%   MP-Opt-Model
%   Copyright (c) 2016-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% the event function data to be registered
e = struct( ...
        'name', name, ...
        'fcn', fcn, ...
        'tol', tol, ...
        'locate', locate ...
    );

%% convert function names to function handles, as necessary
if ~isa(e.fcn, 'function_handle')
    e.fcn = str2func(e.fcn);
end

%% register to list of event functions
if isempty(event_list)
    event_list = e;
else
    nef = length(event_list);
    for k = 1:nef
        if strcmp(event_list(k).name, name)
            error('cpf_register_event: duplicate event name: ''%s''', name);
        end
    end
    event_list(nef+1) = e;
end
