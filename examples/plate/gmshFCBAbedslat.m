function varargout = gmshFCBAbedslat(Q,P,clQ,clP,filename,indim,varargin)
% function varargout = gmshFCBAbedslat(Q,P,clQ,clP,filename,indim,varargin)
% Q : QUADRANGLE
% P : POINT
% clQ, clP : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getindim(Q) by default)

if nargin<4
    clP = clQ;
end
if nargin<6
    indim = getindim(Q);
end

PQ = getvertices(Q);

G = gmshfile(P,clP,1);
G = createpoints(G,PQ,clQ,[5,2:4]);
G = createcontour(G,1:5,1:5,1);
G = createplanesurface(G,1,1);
if ischarin('recombine',varargin)
    G = recombinesurface(G,1);
end

varargin = delonlycharin('recombine',varargin);

if nargin>=5 && ischar(filename)
    G = setfile(G,filename);
end

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(Q):-1:getdim(Q)-n+1,varargin{:});
