function cb_list = pne_register_callback(cb_list, fcn, priority, args)
%PNE_REGISTER_CALLBACK  Register callback functions for PNES_MASTER
%   CB_LIST = PNE_REGISTER_CALLBACK(CB_LIST, FCN, PRIORITY)
%
%   Registers a callback function to be called by PNES_MASTER.
%
%   Inputs:
%       CB_LIST : struct containing info about registered callback fcns
%       FCN : string containing name of callback function
%       PRIORITY : number that determines order of execution for multiple
%                  callback functions, where higher numbers run first,
%                  default priority is 20, where the standard callbacks
%                  are called with the following priority:
%                       pne_nose_event_cb       51
%                       pne_target_lam_event_cb 50
%                       pne_default_callback    0
%^      ARGS : arguments to be passed to the callback each time it is invoked
%
%   Outputs:
%       CB_LIST : updated struct containing info about registered callback fcns
%
%   UNFINISHED - need to update user-defined callback documentation
%   User Defined PNES_MASTER Callback Functions:
%       The user can define their own callback functions which take
%       the same form and are called in the same contexts as
%       PNE_DEFAULT_CALLBACK. These are specified via the MATPOWER
%       option 'cpf.user_callback'. This option can be a string containing
%       the name of the callback function, or a struct with the following
%       fields, where all but the first are optional:
%           'fcn'       - string with name of callback function
%           'priority'  - numerical value specifying callback priority
%                (default = 20, see PNE_REGISTER_CALLBACK for details)
%           'args'      - arbitrary value (any type) passed to the callback
%                         as CB_ARGS each time it is invoked
%       Multiple user callbacks can be registered by assigning a cell array
%       of such strings and/or structs to the 'cpf.user_callback' option.
%
%   See also PNES_MASTER, PNE_DEFAULT_CALLBACK.

%   MP-Opt-Model
%   Copyright (c) 2016-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% default priority
if nargin < 4
    args = [];
    if nargin < 3
        priority = [];
    end
end
if isempty(priority)
    priority = 20;
end

cb = struct( ...
        'fcn', fcn, ...
        'priority', priority, ...
        'args', args ...
    );
if ~isa(cb.fcn, 'function_handle')
    cb.fcn = str2func(cb.fcn);
end

%% add it to the list
if isempty(cb_list)
    cb_list = cb;           %% first one
else
    ncb = length(cb_list) + 1;
    cb_list(ncb) = cb;      %% append
    
    %% sort by descending priority
    p = cell(ncb, 1);
    [p{:}] = deal(cb_list.priority);
    [junk, i] = sort(cell2mat(p), 'descend');
    cb_list = cb_list(i);
end
