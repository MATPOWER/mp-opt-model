function display(om, varargin)
% display - Displays the object.
%
% Called when semicolon is omitted at the command-line. Displays the details
% of the variables, constraints, costs included in the model.
%
% See also mp.opt_model.

%   MP-Opt-Model
%   Copyright (c) 2008-2025, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin < 2
    more_set_types = {};
else
    more_set_types = varargin{1};
end

%% display details of each set type
set_types = om.get_set_types();
set_types = horzcat(set_types, more_set_types);
fprintf('\n');
for k = 1:length(set_types)
    om.(set_types{k}).display(set_types{k});
end

%% user data
fields = fieldnames(om.userdata);
if ~isempty(fields)
    fprintf('\nUSER DATA\n')
    fprintf('=========\n')
    fprintf('  name                               size       class\n');
    fprintf(' ------------------------------   -----------  --------------------\n');
    for k = 1:length(fields)
        f = om.userdata.(fields{k});
        [m, n] = size(f);
        fprintf('  %-31s %5dx%-5d   %s\n', fields{k}, m, n, class(f));
    end
else
    fprintf('USER DATA                   :  <none>\n');
end
