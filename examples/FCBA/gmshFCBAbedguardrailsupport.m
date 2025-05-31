function varargout = gmshFCBAbedguardrailsupport(Q,L,clQ,clL,filename,indim,varargin)
% function varargout = gmshFCBAbedguardrailsupport(Q,L,clQ,clL,filename,indim,varargin)
% Q : QUADRANGLE
% L : LIGNE
% clQ, clL : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getindim(Q) by default)

if nargin<4
    clL = clQ;
end
if nargin<6
    indim = getindim(Q);
end

PQ = getvertices(Q);

G = gmshfile(L,clL,[2 1],1);
G = createpoints(G,PQ,clQ,[6,3:5]);
G = createcontour(G,2:6,2:6,1);
G = createplanesurface(G,1,1);
if ischarin('recombine',varargin)
    G = recombinesurface(G,1);
end
G = embedcurveinsurface(G,1,1);

varargin = delonlycharin('recombine',varargin);

if nargin>=5 && ischar(filename)
    G = setfile(G,filename);
end

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(Q):-1:getdim(Q)-n+1,varargin{:});
