function varargout = gmshFCBAdesk3(Q,I,L1,L2,LbQ,CiQeI,CiI,clQ,clI,clL1,clL2,clLbQ,clCiQeI,clCiI,filename,indim,varargin)
% function varargout = gmshFCBAdesk3(Q,I,L1,L2,PbQ,CiQeI,CiI,clQ,clI,clL1,clL2,clPbQ,clCiQeI,clCiI,filename,indim,varargin)
% Q : QUADRANGLE
% I : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE
% L1, L2 : LIGNE
% LbQ : LIGNE
% CiQeI, CiI : CIRCLE
% clQ, clI, clL1, clL2, clLbQ, clCiQeI, clCiI : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getdim(Q) by default)

if nargin<9
    clI = clQ;
end
if nargin<10
    clL1 = clQ;
end
if nargin<11
    clL2 = clQ;
end
if nargin<12
    clLbQ = clQ;
end
if nargin<13
    clCiQeI = clQ;
end
if nargin<14
    clCiI = clI;
end
if nargin<16
    indim = getdim(Q);
end

if ~iscell(LbQ)
    LbQ = {LbQ};
end
if length(clLbQ)==1
    clLbQ = repmat(clLbQ,1,length(LbQ));
end

G = GMSHFILE();
if nargin>=15 && ischar(filename)
    G = setfile(G,filename);
end

PQ = getvertices(Q);
numpoints = 1:(length(PQ)+2*length(LbQ));
numlines = 1:(length(PQ)+2*length(LbQ));
G = createpoints(G,PQ,clQ,numpoints(1:5:end));
npoints = setdiff(numpoints,numpoints(1:5:end));
for k=1:length(LbQ)
    PbQ = getvertices(LbQ{k});
    G = createpoints(G,PbQ,clLbQ(k),npoints([2*k-1 2*k]));
end
G = createcontour(G,numpoints,numlines,1);

numpoints = numpoints(end)+(1:5);
numlines = numlines(end)+(1:4);
GCiI = gmshfile(CiI,clCiI,numpoints(1),numpoints(2:end),numlines,2);
G = G+GCiI;
G = createplanesurface(G,2,2);

numpoints = numpoints(end)+(1:5);
numlines = numlines(end)+(1:4);
if isa(I,'DOMAIN') || isa(I,'QUADRANGLE')
    GI = gmshfile(I,clI,numpoints(1:end-1),numlines,3);
elseif isa(I,'CIRCLE') || isa(I,'ELLIPSE')
    GI = gmshfile(I,clI,numpoints(1),numpoints(2:end),numlines,3);
end
G = G+GI;
G = createplanesurface(G,[2 3],3);
if ischarin('recombine',varargin)
    G = recombinesurface(G,3);
end

numpoints = numpoints(end)+(1:5);
numlines = numlines(end)+(1:4);
CiQeI = gmshfile(CiQeI,clCiQeI,numpoints(1),numpoints(2:end),numlines,4);
G = G+CiQeI;
G = createplanesurface(G,4,4);

G = createplanesurface(G,[1 3 4],1);

PL1 = getvertices(L1);
PL2 = getvertices(L2);
numpoints = numpoints(end)+(1:4);
numlines = numlines(end)+(1:2);
seg = [1:2;3:4];
seg = numpoints(seg);
G = createpoints(G,PL1,clL1,numpoints(1:2));
G = createpoints(G,PL2,clL2,numpoints(3:4));
G = createlines(G,seg,numlines);
G = embedlinesinsurface(G,numlines,1);

if ischarin('recombine',varargin)
    G = recombinesurface(G,1);
end

varargin = delonlycharin('recombine',varargin);

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(Q):-1:getdim(Q)-n+1,varargin{:});
