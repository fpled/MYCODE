function varargout = gmshFCBAbedsiderail(Q,L,clQ,clL,filename,indim,varargin)
% function varargout = gmshFCBAbedsiderail(Q,L,clQ,clL,filename,indim,varargin)
% Q : QUADRANGLE
% L : cell of LIGNE
% clQ, clL : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getindim(Q) by default)

if nargin<4
    clL = clQ;
end
if nargin<6
    indim = getindim(Q);
end

if ~iscell(L)
    L = {L};
end
if isscalar(clL)
    clL = repmat(clL,1,length(L));
end

G = GMSHFILE();
if nargin>=5 && ischar(filename)
    G = setfile(G,filename);
end

PQ = getvertices(Q);

G = createpoints(G,PQ,clQ,1:4);

numpoints = 4+(1:2*length(L));
numlines = 5+(1:length(L));
for j=1:length(L)
    GI = gmshfile(L{j},clL(j),numpoints([2*j-1,2*j]),numlines(j));
    G = G+GI;
end
G = createcontour(G,[1:3,numpoints(end),4],1:5,1);
G = createplanesurface(G,1,1);
G = embedcurvesinsurface(G,numlines,1);

if ischarin('recombine',varargin)
    G = recombinesurface(G,1);
end

varargin = delonlycharin('recombine',varargin);

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(Q):-1:getdim(Q)-n+1,varargin{:});
