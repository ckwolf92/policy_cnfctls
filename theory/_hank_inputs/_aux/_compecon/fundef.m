% FUNDEF Creates a function family definition structure and/or performs checks
% USAGE:
%   g=FUNDEF({bastype1,p11,p12,...},{bastype2,p21,p22,...},...);
% or
%   [g,errorstr]=FUNDEF({bastype1,p11,p12,...},{bastype2,p21,p22,...},...);
%
% The first syntax will generate an error message and halt execution if
%   an error is encountered. The second syntax sets g=[] and returns an
%   error message if an error in encountered, otherwise it returns a 
%   function definition structure and an empty string (test using if isempty(errorstr)).
%
% FUNDEF accepts any number of input arguments, each of which must be a cell array
% containing the following information:
%   bastype     : the string prefix referencing the function family name 
%                 (e.g., 'cheb', 'spli')
%   p1, p2, ... : the parameters required by the specified family
%
%  The family definition parameters required by a specific family are
%  documented in that family's DEF file. For example the CHEB family requires n, the
%  degree of the polynomial (1 plus its order), and a and b, the left and right endpoints.
%
% Examples:
%       cdef=funbas({'cheb',10,0,1});
%   defines a 1-D polynomial of order 9 on the interval [0,1].
%       cdef=funbas({'cheb',10,0,1},{'cheb',5,-1,1});
%   defines a 2-D polynomial:
%     the first dimension is order 9 on the interval [0,1],
%     the second dimension is order 4 on the interval [-1 1].
%       cdef=funbas({'cheb',10,0,1},{'spli',[-1;1],5});
%   defines a 2-D function:
%     the first dimension is an order 9 polynomialon the interval [0,1],
%     the second dimension is a cubic spline with 5 evenly spaced breakpoints 
%       on the interval [-1 1].
%
% To test for a valid coefficient structure use
%   [g,errorstr]=FUNDEF(G);
% On return g=G if no errors are found, [] if an error is found.
%   ERRORSTR contains a string describing the error and is empty 
%   when no errors are detected (use if isempty(errorstr) to check).
%
% See also: FUNEVAL, FUNBAS, CHEBDEF, SPLIDEF, LINDEF, FOURDEF

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [g,errorstr]=FUNDEF(varargin);

errorstr=lasterr;    % clears the last error
errorstr='';
g=[];
% !!!!!!!!!! perform consistency checks !!!!!!!!!!!!!!
if isstruct(varargin{1})
  if nargin>1               
    errorstr='Pass a single function structure for consistency check'; 
    if nargout<2, error(errorstr); else, return; end
  end
  c=varargin{1};
  f=fieldnames(c);
  if ~strcmp(f,{'d';'n';'bastype';'parms'})
    errorstr='Structure has improper fields'; 
    if nargout<2, error(errorstr); else, return; end
  end
  [m,d]=size(c.n)
  if m~=1 | d~=c.d
    errorstr='''n'' field of improper size'; 
    if nargout<2, error(errorstr); else, return; end
  end
  [m,d]=size(c.bastype);
  if m~=1 | d~=c.d
    errorstr='''bastype'' field of improper size'; 
    if nargout<2, error(errorstr); else, return; end
  end
  [m,d]=size(c.parms);
  if m~=1 | d~=c.d
    errorstr='parameter field of improper size'; 
    if nargout<2, error(errorstr); else, return; end
  end
  if size(c.vals,1)~=prod(n)
    errorstr='Rows in ''vals'' field does not equal product of ''n'' field';
    if nargout<2, error(errorstr); else, return; end
  end
  if ~iscell('c.bastype')
    errorstr='''bastype'' field should be 1xd cell array (cell string)';
    if nargout<2, error(errorstr); else, return; end
  end
  if ~iscell('c.parms')
    errorstr='''parms'' field should be 1xd cell array'; 
    if nargout<2, error(errorstr); else, return; end
  end
  g=c;  % no errors found

% !!!!!!!! Create a new coefficient structure !!!!!!!!

else 
  errorstr='';
  g=[];
  if isa(varargin{1},'cell') d=length(varargin); else, d=1; varargin(1)=varargin{:}; end
  n=zeros(1,d);
  b=zeros(1,d);
  a=zeros(1,d);
  parms=cell(1,d);

  for j=1:d
    bastype{j}=varargin{j}{1};
    varargin{j}(1)=[];
  
    if ~exist([bastype{j} 'bas'])
      errorstr=['Cannot find M-file ' lower(bastype{j}) 'bas.m']; 
      if nargout<2, error(errorstr); else, return; end
    end
  end 
  
  warning=[];
  for j=1:d
    auxname=[bastype{j} 'def'];
    if exist(auxname,'file');  % check if DEF file exists
      m=nargin(lower(auxname));
      if m>0 & m<length(varargin{j})
        errorstr=['Too many parameters for ' lower(bastype{j}) 'def'];
        if nargout<2, error(errorstr); else, return; end
      end
      eval(['[n(j),a(j),b(j),parms{j}]=' bastype{j} 'def(varargin{j}{:});']);%,'errorstr=lasterr;');
      if ~isempty(errorstr)
        if nargout<2, error(errorstr); else, return; end
      end  
    else
      warning={warning;upper(auxname)};
    end
  end
  if ~isempty(warning)
    disp('Error:')
    disp('The following parameter definition files could not be found: ')
    disp(warning{2:end});
    error(' ');
  end  
  
  % Assign values to fields
  g=struct('d',[],'n',[],'a',[],'b',[],'bastype',[],'parms',[]);
  g.d=d;
  g.n=n;
  g.a=a;
  g.b=b;
  g.bastype=bastype;
  g.parms=parms;
end
