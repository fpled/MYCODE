function varargout = gmshFCBAbedguardrailsupport(Q,L_slat,clQ,clL_slat,filename,indim,varargin)
% function varargout = gmshFCBAbedguardrailsupport(Q,L_slat,clQ,clL_slat,filename,indim,varargin)
% Q : QUADRANGLE
% L_slat : LIGNE
% clQ, clslat : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getindim(Q) by default)

if nargin<4
    clL_slat = clQ;
end
if nargin<8
    indim = getindim(Q);
end

G = GMSHFILE();
if nargin>=5 && ischar(filename)
    G = setfile(G,filename);
end

PQ = getvertices(Q);
G = createpoints(G,PQ,clQ,1:4);
G = createcontour(G,1:7,1:7,1);
G = createplanesurface(G,1,1);

if ischarin('recombine',varargin)
    G = recombinesurface(G,1);
end

varargin = delonlycharin('recombine',varargin);

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(Q):-1:getdim(Q)-n+1,varargin{:});
