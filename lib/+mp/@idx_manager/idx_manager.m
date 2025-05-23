classdef idx_manager < handle
% mp.idx_manager - |MATPOWER| Index Manager abstract class
% ::
%
%   A MATPOWER Index Manager object can be used to manage the indexing of
%   various named and indexed blocks of various set types, such as variables,
%   constraints, etc. This class helps keep track of the ordering and
%   indexing of the various blocks as they are added to the object.
%
%   Properties
%       userdata    - a struct containing arbitrary data added by the user
%
%   Public Methods
%       copy - Makes a shallow copy of the object.
%
%       display - (must be implemented in the subclass) Displays the
%           object (called automatically when you omit the semicolon at
%           the command-line).
%
%       display_set - Prints to screen the indexing details for the
%           specified set type. Intended to be called by DISPLAY method.
%
%       get_idx - Returns index structure(s) for specified set type(s),
%           with starting/ending indices and number of elements for
%           each named (and optionally indexed) block.
%
%       get_userdata - Retreives values of user data stored in the object.
%
%       get - Return the value of any individual field.
%
%       get_set_types - (must be implemented in the subclass)
%           Returns a cell array of set types, that is, the names of the
%           properties containing the mp.set_manager objects.
%
%   The following is the structure of the data in the object, using a set
%   type named 'var' for illustration. Each field of .idx or .data is a
%   struct whose field names are the names of the corresponding blocks of
%   elements of that type (e.g. variables, constraints, etc.). They are
%   found in order in the corresponding .order field. The description next
%   to these fields gives the meaning of the value for each named sub-field.
%   E.g. obj.var.data.v0.Pg contains a vector of initial values for the 'Pg'
%   block of variables.
%
%   obj
%       .var        - data for 'var' set type, e.g. variable sets that
%                     make up the full optimization variable x
%           .idx
%               .i1 - starting index within x
%               .iN - ending index within x
%               .N  - number of elements in this variable set
%           .N      - total number of elements in x
%           .NS     - number of variable sets or named blocks
%           .data   - additional set-type-specific data for each block
%           .order  - struct array of names/indices for variable
%                     blocks in the order they appear in x
%               .name   - name of the block, e.g. Pg
%               .idx    - indices for name, {2,3} => Pg(2,3)
%       .userdata   - any user defined data
%           .(user defined fields)

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%    es = struct();

    properties
        userdata = struct();
    end     %% properties

    methods
        function obj = idx_manager(s)
            % Constructor.
            % ::
            %
            %   obj = mp.idx_manager()
            %   obj = mp.idx_manager(a_struct)
            %   obj = mp.idx_manager(an_obj)

            if nargin > 0
                if isa(s, 'mp.idx_manager')
                    %% this copy constructor will not be inheritable under
                    %% Octave until the fix has been included for:
                    %%      https://savannah.gnu.org/bugs/?52614
                    if have_feature('octave')
                        s1 = warning('query', 'Octave:classdef-to-struct');
                        warning('off', 'Octave:classdef-to-struct');
                    end
                    props = fieldnames(s);
                    if have_feature('octave')
                        warning(s1.state, 'Octave:classdef-to-struct');
                    end
                    for k = 1:length(props)
                        obj = copy_prop(s, obj, props{k});
                    end
                elseif isstruct(s)
                    props = fieldnames(obj);
                    for k = 1:length(props)
                        if isfield(s, props{k})
                            obj = copy_prop(s, obj, props{k});
                        end
                    end
                else
                    error('mp.idx_manager.idx_manager: input must be an ''mp.idx_manager'' object or a struct');
                end
            end
        end

        function new_obj = copy(obj)
            % Duplicate the object.

            %% initialize copy
            new_obj = eval(class(obj));  %% create new object

            %% copy properties/fields
            if have_feature('octave')
                s1 = warning('query', 'Octave:classdef-to-struct');
                warning('off', 'Octave:classdef-to-struct');
            end
            props = fieldnames(obj);
            if have_feature('octave')
                warning(s1.state, 'Octave:classdef-to-struct');
            end
            for k = 1:length(props)
                new_obj = copy_prop(obj, new_obj, props{k});
            end
        end

        function display_set(obj, stype, sname)
            % Display indexing information for a given set type.

            obj.(stype).display(stype);
        end

        varargout = get_idx(obj, varargin)

        rv = get_userdata(obj, name)

        val = get(obj, varargin)
    end     %% methods
end         %% classdef

function d = copy_prop(s, d, prop)
    if isa(s.(prop), 'mp.set_manager')
        d.(prop) = s.(prop).copy();
    elseif isa(d.(prop), 'mp.set_manager')
        d.(prop) = nested_struct_copy( ...
            d.(prop), s.(prop));
    else
        d.(prop) = s.(prop);
    end
end
