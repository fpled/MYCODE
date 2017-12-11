function varargout = gmshFCBAdesk12(Q,La,Lb,clQ,clLa,clLb,filename,indim,varargin)
% function varargout = gmshFCBAdesk12(Q,La,Lb,clQ,clLa,clLb,filename,indim,varargin)
% Q : QUADRANGLE
% La, Lb : LIGNE
% clQ, clLa, clLb : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getdim(Q) by default)

if nargin<5
    clLa = clQ;
end
if nargin<6
    clLb = clQ;
end
if nargin<8
    indim = getdim(Q);
end

G = GMSHFILE();
if nargin>=7 && ischar(filename)
    G = setfile(G,filename);
end

PQ = getvertices(Q);
PLa = getvertices(La);
PLb = getvertices(Lb);

G = createpoints(G,PQ,clQ,1:4);
G = createpoints(G,PLa,clLa,7:8);
G = createpoints(G,PLb,clLb,5:6);

G = createcontour(G,[1 5 2:4],1:5,1);
G = createplanesurface(G,1,1);

seg = [5 6;7 8];
G = createlines(G,seg,6:7);
G = embedlinesinsurface(G,6:7,1);

if ischarin('recombine',varargin)
    G = recombinesurface(G,1);
end

varargin = delonlycharin('recombine',varargin);

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(Q):-1:getdim(Q)-n+1,varargin{:});
