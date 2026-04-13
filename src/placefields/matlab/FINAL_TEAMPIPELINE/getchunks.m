function [d, id] = getchunks(a, opt)
%GETCHUNKS Get the number of repetitions that occur in consecutive chunks.
%   C = GETCHUNKS(A) returns an array of n elements, where n is the number
%   of consecutive chunks (2 or more repetitions) in A, and each element is
%   the number of repetitions in each chunk. A can be LOGICAL, any
%   numeric vector, or CELL array of strings. It can also be a character
%   array (see below, for its special treatment).
%
%   [C, I] = GETCHUNKS(A) also returns the indices of the beginnings of the
%   chunks.
%
%   If A is a character array, then it finds words (consecutive
%   non-spaces), returning the number of chararcters in each word and the
%   indices to the beginnings of the words.
%
%   GETCHUNKS(A, OPT) accepts an optional argument OPT, which can be any of
%   the following three:
%
%       '-reps'  : return repeating chunks only. (default)
%       '-full'  : return chunks including single-element chunks.
%       '-alpha' : (for CHAR arrays) only consider alphabets and numbers as
%                  part of words. Punctuations and symbols are regarded as
%                  spaces.
%
%   Examples:
%     A = [1 2 2 3 4 4 4 5 6 7 8 8 8 8 9];
%     getchunks(A)
%       ans =
%           2   3   4
%
%
%     B = 'This is a generic (simple) sentence';
%     [C, I] = getchunks(B)
%       C =
%            4     2     1     7     8     8
%       I =
%            1     6     9    11    19    28
%
%
%     [C, I] = getchunks(B, '-alpha')
%       C =
%            4     2     1     7     6     8
%       I =
%            1     6     9    11    20    28
%
%   See also HIST, HISTC.
%
%   VERSIONS:
%     v1.0 - first version
%     v1.0.1 - added option '-alpha'
%
% Copyright 2009 The MathWorks, Inc.
%--------------------------------------------------------------------------
% Error checking
%--------------------------------------------------------------------------
error(nargchk(1, 2, nargin));
if ndims(a) > 2 || min(size(a)) > 1
  error('Input must be a 2-D vector');
end
alphanumeric = false;
fullList     = false;
%--------------------------------------------------------------------------
% Process options
%--------------------------------------------------------------------------
if nargin == 2
  if ~ischar(opt)
    error('Additional argument must be a string array');
  end
  
  % Allow for partial arguments
  possibleOptions = ['-full '; '-reps '; '-alpha'];
  iOpt = strmatch(lower(opt), possibleOptions);
  
  if isempty(iOpt) || length(iOpt) > 1
    error('Invalid option. Allowed option: ''-full'', ''-reps'', ''-alpha''');
  else
    switch iOpt
      
      case 1  % '-full'
        % Include single-element chunks
        fullList = true;
        if ischar(a)
          fprintf('''-full'' option not applicable to CHAR arrays.\n');
        end
        
      case 2  % '-reps'
        % Only find 2 or more repeating blocks
        fullList = false;
        
      case 3  % '-alpha'
        % For char arrays, only consider alphabets and numbers as part of
        % words. Punctuations and symbols are regarded as space.
        alphanumeric = true;
        if ~ischar(a)
          fprintf('''-alpha'' option only applicable to CHAR arrays.\n');
        end
        
    end
  end
end
%--------------------------------------------------------------------------
% Convert to a row vector for STRFIND
%--------------------------------------------------------------------------
a = a(:)';
%--------------------------------------------------------------------------
% Deal with differet classes
%--------------------------------------------------------------------------
switch class(a)
  
  case 'double'
    % Leave as is
    
  case {'logical', 'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', 'single'}
    % Convert to DOUBLE
    a = double(a);
    
  case 'char'
    if alphanumeric % Get alphabet and number locations
      try % call C-helper function directly (slightly faster)
        a = isletter(a) | ismembc(a, 48:57);
      catch %#ok<CTCH>
        a = isletter(a) | ismember(a, 48:57);
      end
      
    else  % Get non-space locations
      a = ~isspace(a);  
    end
  
  case 'cell'
    % Convert cell array of strings into unique numbers
    if all(cellfun('isclass', a, 'char'))
      [tmp, tmp, a] = unique(a); %#ok<ASGLU>
    else
      error('Cell arrays must be array of strings.');
    end
    
  otherwise
    error('Invalid type. Allowed type: CHAR, LOGICAL, NUMERIC, and CELL arrays of strings.');
end
%--------------------------------------------------------------------------
% Character arrays (now LOGICAL) are dealt differently
%--------------------------------------------------------------------------
if islogical(a)
  % Pad the array
  a  = [false, a, false];
  % Here's a very convoluted engine
  b  = diff(a);
  id = strfind(b, 1);
  d  = strfind(b, -1) - id;
%--------------------------------------------------------------------------
% Everything else (numeric arrays) are processed here
else
  % Pad the array
  a                 = [NaN, a, NaN];
  % Here's more convoluted code
  b                 = diff(a);
  b1                = b;  % to be used in fullList (below)
  ii                = true(size(b));
  ii(strfind(b, 0)) = false;
  b(ii)             = 1;
  c                 = diff(b);
  id                = strfind(c, -1);
  
  % Get single-element chunks also
  if fullList
  
    % And more convoluted code
    b1(id)          = 0;
    ii2             = find(b1(1:end-1));
    d               = [strfind(c, 1) - id + 1, ones(1, length(ii2))];
    id              = [id,ii2];
    [id,tmp]        = sort(id);
    d               = d(tmp);
    
  else
    
    d               = strfind(c, 1) - id + 1;
    
  end
end