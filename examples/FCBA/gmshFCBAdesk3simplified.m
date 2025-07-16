function varargout = gmshFCBAdesk3simplified(Q,I,L1,L2,PbQ,PiQeI,PiI,clQ,clI,clL1,clL2,clPbQ,clPiQeI,clPiI,filename,indim,varargin)
% function varargout = gmshFCBAdesk3simplified(Q,I,L1,L2,PbQ,PiQeI,PiI,clQ,clI,clL1,clL2,clPbQ,clPiQeI,clPiI,filename,indim,varargin)
% Q : QUADRANGLE
% I : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE
% L1, L2 : LIGNE
% PbQ, PiQeI, PiI : POINT
% clQ, clI, clL1, clL2, clPbQ, clPiQeI, clPiI : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getindim(Q) by default)

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
    clPbQ = clQ;
end
if nargin<13
    clPiQeI = clQ;
end
if nargin<14
    clPiI = clI;
end
if nargin<16
    indim = getindim(Q);
end

if ~iscell(PbQ)
    PbQ = {PbQ};
end
if isscalar(clPbQ)
    clPbQ = repmat(clPbQ,1,length(PbQ));
end

if ~iscell(PiQeI)
    PiQeI = {PiQeI};
end
if isscalar(clPiQeI)
    clPiQeI = repmat(clPiQeI,1,length(PiQeI));
end

if ~iscell(PiI)
    PiI = {PiI};
end
if isscalar(clPiI)
    clPiI = repmat(clPiI,1,length(PiI));
end

G = GMSHFILE();
if nargin>=15 && ischar(filename)
    G = setfile(G,filename);
end

PQ = getvertices(Q);
numpoints = 1:(length(PQ)+length(PbQ));
numlines = 1:(length(PQ)+length(PbQ));
G = createpoints(G,PQ,clQ,numpoints(1:3:end));
G = createpoints(G,PbQ,clPbQ,setdiff(numpoints,numpoints(1:3:end)));
G = createcontour(G,numpoints,numlines,1);

numpoints = numpoints(end)+(1:5);
numlines = numlines(end)+(1:4);
if isa(I,'DOMAIN') || isa(I,'QUADRANGLE')
    GI = gmshfile(I,clI,numpoints(1:end-1),numlines,2);
elseif isa(I,'CIRCLE') || isa(I,'ELLIPSE')
    GI = gmshfile(I,clI,numpoints(1),numpoints(2:end),numlines,2);
end
G = G+GI;
G = createplanesurface(G,2,2);
if ~isempty([PiI{:}])
    numpoints = numpoints(end)+(1:length(PiI));
    G = createpoints(G,PiI,clPiI,numpoints);
    G = embedpointsinsurface(G,numpoints,2);
end
if ischarin('recombine',varargin)
    G = recombinesurface(G,2);
end

G = createplanesurface(G,[1 2],1);
if ~isempty([PiQeI{:}])
    numpoints = numpoints(end)+(1:length(PiQeI));
    G = createpoints(G,PiQeI,clPiQeI,numpoints);
    G = embedpointsinsurface(G,numpoints,1);
end

PL1 = getvertices(L1);
PL2 = getvertices(L2);
numpoints = numpoints(end)+(1:4);
numlines = numlines(end)+(1:2);
seg = [1:2;3:4];
seg = numpoints(seg);
G = createpoints(G,PL1,clL1,numpoints(1:2));
G = createpoints(G,PL2,clL2,numpoints(3:4));
G = createlines(G,seg,numlines);
G = embedcurvesinsurface(G,numlines,1);

if ischarin('recombine',varargin)
    G = recombinesurface(G,1);
end

varargin = delonlycharin('recombine',varargin);

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(Q):-1:getdim(Q)-n+1,varargin{:});
