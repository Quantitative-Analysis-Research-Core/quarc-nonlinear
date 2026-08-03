function varargout = RQA(varargin) %#ok<STOUT,VANUS>
%RQA (removed) Use rqa, crqa, jrqa or mdrqa.
%   The combined RQA entry point is not part of this library. It took a TYPE
%   argument and dispatched internally; the four analyses now have their own
%   functions with validated arguments.
%
%   RQA(x,'RQA',...)    -> rqa(x, delay, dim, param, threshold, ...)
%   RQA(x,'cRQA',...)   -> crqa(...)
%   RQA(x,'jRQA',...)   -> jrqa(...)
%   RQA(x,'mdRQA',...)  -> mdrqa(...)
%
%   The argument lists differ, so this is not a mechanical substitution. See
%   the help for each function.
%
%   The previous implementation remains in the archived NONAN Library until
%   31 December 2026 if results need reproducing exactly.
%
%   See also RQA, CRQA, JRQA, MDRQA.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

error('RQA:removed', ...
    ['RQA has been replaced by rqa, crqa, jrqa and mdrqa, which take ' ...
     'validated arguments instead of a TYPE string. The argument lists ' ...
     'differ, so this is not a mechanical substitution; see the help for ' ...
     'each. The previous implementation remains in the archived NONAN ' ...
     'Library until 31 December 2026.']);
end
