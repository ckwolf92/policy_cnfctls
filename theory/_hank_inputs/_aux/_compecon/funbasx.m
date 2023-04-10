% FUNBASX Creates basis structures for function evaluation
% Note: Use FUNBAS to obtain single basis matrices
% USAGE
%   B=funbasx(fspace,x,order,bformat);
% INPUTS
%   fspace  : a structure defining a family of functions (use FUNDEF to create this)
%   x       : an mxd matrix or a 1xd cell array of columns vectors
%                   (default: created by calling FUNNODE)
%   order   : a d-column matrix (a single column is expanded)
%                   (default: zeros(1,d)))
%   bformat : a string: 'tensor', 'direct' or 'expanded'
%                   defaults: 'expanded' if d=1, otherwise
%                               'tensor' if x is a cell array
%                               'direct' if x is a matrix 
% OUTPUTS
%   B : a basis structure (defined below)
%   x : the computed evaluation points if none are passed
%
%   A basis structure has the following fields:
%     B.vals        - a cell array containing the basis data 
%                         (see exception below)
%     B.format      - 'tensor', 'direct', 'expanded'
%     B.order       - the orders of differentiation ('expanded' format) 
%                   - smallest orders of differentiation
%                       and number of bases ('tensor' and 'direct' formats)
%
%     ORDER Determines the # of basis matrices created 
%           (i.e., the # of elements in B.bases)
%       for 'tensor' and 'direct' formats ORDER should be
%         1xd if only a single basis matrix is needed in each dimension
%         2xd specifying the minimum and maximum orders in each dimension
%         kxd listing of all desired basis matrices 
%              (in this form only the min and max of ORDER are used)
%       for 'expanded' format
%         kxd listing of all desired basis matrices
%
%   If d=1 the format will be 'expanded' regardless of bformat
%
% USES: gridmake, funnode, funbconv
%
% See also: FUNDEF, FUNNODE, FUNEVAL, CHEBBAS, SPLIBAS, LINBAS, FOURBAS

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [B,x]=funbasx(fspace,x,order,bformat)
  if nargin<1 | ~isstruct(fspace) 
     error('A coefficient structure must be specified')
  end
  if nargin<2 x=[]; end
  if nargin<3 | isempty(order), order=0; end 
  if nargin<4, bformat=[]; end

  % Determine the problem dimension
  d=length(fspace.n);
  
  % Expand ORDER if it has a single columns
  if d>1 & size(order,2)==1, order=order*ones(1,d); end 
  
% Initialize basis structure 
  m=size(order,1);
  if m>1 
     minorder=min(order)+zeros(1,d); 
     numbases=(max(order)-minorder)+1;
  else 
    minorder=order+zeros(1,d);
    numbases=ones(1,d); 
  end
 
  B=struct('vals',[],'order',minorder,'format',bformat);
  B.vals=cell(max(numbases),d);
  
  if isempty(x) x=funnode(fspace); end

  if isempty(bformat)
    if isa(x,'cell') bformat='tensor';
    else             bformat='direct';
    end
  end
  
  if d>1
    if ~isa(x,'cell') & strcmp(bformat,'tensor')
      error('Must pass a cell array to form a tensor format basis structure')
    end
    if isa(x,'cell') & strcmp(bformat,'direct')   
      % it would be more efficient in this case to
      % use the tensor form to compute the bases and then
      % to use indexing to expand to the direct form
%      x=gridmake(x); % convert to grid for direct form
    end
  end

  if isa(x,'cell') B.format='tensor';
  else             B.format='direct';
  end

  % Compute basis matrices
  switch B.format
  case 'tensor'
    for j=1:d
      if (m>1) orderj=unique(order(:,j));
      else, orderj=order(1,j);
      end
      if length(orderj)==1
       B.vals{1,j}=feval([fspace.bastype{j} 'bas'],fspace.parms{j}{:},x{j},orderj);
      else
       B.vals(orderj-minorder(j)+1,j)=feval([fspace.bastype{j} 'bas'],fspace.parms{j}{:},x{j},orderj);
      end
    end
  case 'direct'
    for j=1:d      
      if (m>1) orderj=unique(order(:,j));
      else, orderj=order(1,j);
      end
      if length(orderj)==1
        B.vals{1,j}=feval([fspace.bastype{j} 'bas'],fspace.parms{j}{:},x(:,j),orderj);
      else
        B.vals(orderj-minorder(j)+1,j)=feval([fspace.bastype{j} 'bas'],fspace.parms{j}{:},x(:,j),orderj);
      end
    end
  end
  
  % d=1, switch to expanded format
  if size(B.vals,2)==1      
    B.format='expanded';
    B.order=order;
    B.vals=B.vals(order+(1-min(order)));
    return
  end
  
% Create expanded format
  switch bformat
  case 'expanded'
    B=funbconv(B,order,'expanded');
  case 'direct'
    if isa(x,'cell'), B=funbconv(B,order,'direct'); end
  end
  