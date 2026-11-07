function s = updateStruct(s, paths, values, allowNew)
%% AI GENERATED
%UPDATESTRUCT Set values deep inside a nested structure from a list of paths.
%
%   s = updateStruct(s, paths, values)
%   s = updateStruct(s, paths, values, allowNew)
%
%   Inputs
%     s        - struct to update (e.g. const)
%     paths    - cell array of paths. Each path is either
%                  * a dot-separated char/string, e.g. 'a.b.c' (relative to s), or
%                  * a cell array of field names, e.g. {'a','b','c'}
%                A single path may also be given on its own.
%     values   - cell array of new values, same length as paths.
%                (A single value may be given on its own if there is one path.)
%     allowNew - (optional, default true) if false, error when a path does not
%                already exist instead of creating new fields.
%
%   Example
%     const.a.b.c = 'test';
%     const.a.d.f = 1;
%     const = updateStruct(const, {'a.b.c', 'a.d.f', 'a.g.h'}, {'new', 42, [1 2 3]});
%
%     % Cell-array paths also work:
%     const = updateStruct(const, {{'a','b','c'}}, {'again'});

    if nargin < 4 || isempty(allowNew)
        allowNew = true;
    end

    % --- normalise inputs -------------------------------------------------
    if ischar(paths) || isstring(paths)
        paths = {paths};
    end

    if ~iscell(values)
        values = {values};
    end

    if numel(paths) ~= numel(values)
        error('updateStruct:sizeMismatch', ...
            'paths (%d) and values (%d) must have the same number of elements.', ...
            numel(paths), numel(values));
    end

    % --- apply each update ------------------------------------------------
    for k = 1:numel(paths)
        p = paths{k};

        if ischar(p) || isstring(p)
            parts = strsplit(char(p), '.');
        elseif iscell(p)
            parts = cellfun(@char, p, 'UniformOutput', false);
        else
            error('updateStruct:badPath', 'Path %d must be a string or cell array of field names.', k);
        end

        if isempty(parts) || any(cellfun(@isempty, parts))
            error('updateStruct:badPath', 'Path %d is empty or has an empty field name.', k);
        end

        if ~allowNew && ~fieldPathExists(s, parts)
            error('updateStruct:missingField', ...
                'Path "%s" does not exist in the structure.', strjoin(parts, '.'));
        end

        % setfield builds intermediate structs automatically
        s = setfield(s, parts{:}, values{k}); 
    end
end

function tf = fieldPathExists(s, parts)
    tf = true;
    for i = 1:numel(parts)
        if isstruct(s) && isscalar(s) && isfield(s, parts{i})
            s = s.(parts{i});
        else
            tf = false;
            return;
        end
    end
end