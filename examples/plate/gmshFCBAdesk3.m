function varargout = gmshFCBAdesk3(Q,I,PL,PbQ,PiQeI,PiI,clQ,clI,clPL,clPbQ,clPiQeI,clPiI,filename,indim,varargin)
% function varargout = gmshFCBAdesk3(Q,I,PL,PbQ,PiQeI,PiI,clQ,clI,clPL,clPbQ,clPiQeI,clPiI,filename,indim,varargin)
% Q : QUADRANGLE
% I : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE
% PL, PbQ, PiQeI, PiI : POINT
% clQ, clI, clPL, clPbQ, clPiQeI, clPiI : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getdim(D) by default)

if nargin<8
    clI = clQ;
end
if nargin<9
    clPL = clQ;
end
if nargin<10
    clPbQ = clQ;
end
if nargin<11
    clPiQeI = clQ;
end
if nargin<12
    clPiI = clI;
end
if nargin<14
    indim = getdim(Q);
end

if ~iscell(PbQ)
    PbQ = {PbQ};
end
if length(clPbQ)==1
    clPbQ = repmat(clPbQ,1,length(PbQ));
end

if ~iscell(PiQeI)
    PiQeI = {PiQeI};
end
if length(clPiQeI)==1
    clPiQeI = repmat(clPiQeI,1,length(PiQeI));
end

if ~iscell(PiI)
    PiI = {PiI};
end
if length(clPiI)==1
    clPiI = repmat(clPiI,1,length(PiI));
end

G = GMSHFILE();
if nargin>=13 && ischar(filename)
    G = setfile(G,filename);
end

numpoints = 1:length(PbQ);
numlines = 1:length(PbQ);
G = createpoints(G,PbQ,clPbQ,numpoints);
G = createcontour(G,numpoints,numlines,1);

numpoints = numpoints(end)+(1:5);
numlines = numlines(end)+(1:4);
if isa(I,'DOMAIN') || isa(I,'QUADRANGLE')
    GI = gmshfile(I,clI,numpoints(1:end-1),numlines,2);
elseif isa(I,'CIRCLE') || isa(I,'ELLIPSE')
    GI = gmshfile(I,clI,numpoints(1),numpoints(2:end),numlines,2);
end
G = G+GI;

numlineloop = [1:length(PbQ),-numlines];
G = createlineloop(G,numlineloop,3);
G = createplanesurface(G,3,1);
G = createplanesurface(G,2,2);
if ~isempty([PiI{:}])
    numpoints = numpoints(end)+(1:length(PiI));
    G = createpoints(G,PiI,clPiI,numpoints);
    G = embedpointsinsurface(G,numpoints,2);
end
if ischarin('recombine',varargin)
    G = recombinesurface(G,2);
end

if ~isempty([PiQeI{:}])
    numpoints = numpoints(end)+(1:length(PiQeI));
    G = createpoints(G,PiQeI,clPiQeI,numpoints);
    G = embedpointsinsurface(G,numpoints,1);
end
if ischarin('recombine',varargin)
    G = recombinesurface(G,1);
end

numpoints = numpoints(end)+(1:4);
numlines = numlines(end)+(1:2);
seg = [1:2;3:4];
seg = numpoints(seg);
G = createpoints(G,PL,clPL,numpoints);
G = createlines(G,seg,numlines);
G = embedlinesinsurface(G,numlines,1);

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(Q):-1:getdim(Q)-n+1);