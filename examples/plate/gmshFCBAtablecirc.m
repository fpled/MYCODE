function varargout = gmshFCBAtablecirc(C,Q,I,PbC,PiCeQ,PiQeI,PiI,clC,clQ,clI,clPbC,clPiCeQ,clPiQeI,clPiI,filename,indim,varargin)
% function varargout = gmshFCBAtablecirc(C,Q,I,PbC,PiCeQ,PiQeI,PiI,clC,clI,clPbC,clPiCeQ,clPiQeI,clPiI,filename,indim,varargin)
% C : CIRCLE
% Q : QUADRANGLE
% I : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE
% PbC, PiCeQ, PiQeI, PiI : POINT
% clC, clQ, clI, clPbC, clPiCeQ, clPiQeI, clPiI : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getdim(D) by default)

if nargin<9
    clQ = clC;
end
if nargin<10
    clI = clC;
end
if nargin<11
    clPbC = clC;
end
if nargin<12
    clPiCeQ = clC;
end
if nargin<13
    clPiQeI = clQ;
end
if nargin<14
    clPiI = clI;
end
if nargin<15
    indim = getdim(C);
end

if ~iscell(PbC)
    PbC = {PbC};
end
if length(clPbC)==1
    clPbC = repmat(clPbC,1,length(PbC));
end

if ~iscell(PiCeQ)
    PiCeQ = {PiCeQ};
end
if length(clPiCeQ)==1
    clPiCeQ = repmat(clPiCeQ,1,length(PiCeQ));
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
if nargin>=15 && ischar(filename)
    G = setfile(G,filename);
end

numcenter = 1;
numpoints = 1+(1:length(PbC));
numlines = 1:length(PbC);
numlineloop = numlines;
% PC = getvertices(C);
Pc = double(getcenter(C));
G = createpoint(G,Pc,clC,numcenter);
G = createpoints(G,PbC,clPbC,numpoints);
G = createcirclecontour(G,numcenter,numpoints,numlines,1);

numpoints = numpoints(end)+(1:4);
numlines = numlines(end)+(1:4);
GQ = gmshfile(Q,clQ,numpoints,numlines,2);
G = G+GQ;
numlineloop = [numlineloop,-numlines];
G = createlineloop(G,numlineloop,3);
G = createplanesurface(G,3,1);
if ~isempty([PiCeQ{:}])
    numpoints = numpoints(end)+(1:length(PiCeQ));
    G = createpoints(G,PiCeQ,clPiCeQ,numpoints);
    G = embedpointsinsurface(G,numpoints,1);
end
if ischarin('recombine',varargin)
    G = recombinesurface(G,1);
end

numpoints = numpoints(end)+(1:5);
numlineloop = numlines;
numlines = numlines(end)+(1:4);
if isa(I,'DOMAIN') || isa(I,'QUADRANGLE')
    GI = gmshfile(I,clI,numpoints(1:end-1),numlines,4);
elseif isa(I,'CIRCLE') || isa(I,'ELLIPSE')
    GI = gmshfile(I,clI,numpoints(1),numpoints(2:end),numlines,4);
end
G = G+GI;
G = createplanesurface(G,4,3);
if ~isempty([PiI{:}])
    numpoints = numpoints(end)+(1:length(PiI));
    G = createpoints(G,PiI,clPiI,numpoints);
    G = embedpointsinsurface(G,numpoints,3);
end
if ischarin('recombine',varargin)
    G = recombinesurface(G,3);
end
numlineloop = [numlineloop,-numlines];
G = createlineloop(G,numlineloop,5);
G = createplanesurface(G,5,2);
if ~isempty([PiQeI{:}])
    numpoints = numpoints(end)+(1:length(PiQeI));
    G = createpoints(G,PiQeI,clPiQeI,numpoints);
    G = embedpointsinsurface(G,numpoints,2);
end
if ischarin('recombine',varargin)
    G = recombinesurface(G,2);
end

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(C):-1:getdim(C)-n+1);
