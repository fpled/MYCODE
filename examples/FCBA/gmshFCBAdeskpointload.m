function varargout = gmshFCBAdeskpointload(Q1,Q2,Q3,Q5a,Q5b,I,PbQ3,PiQ3eI,PiI,clQ1,clQ2,clQ3,clQ5a,clQ5b,clI,clPbQ3,clPiQ3eI,clPiI,filename,indim,varargin)
% function varargout = gmshFCBAdeskpointload(Q1,Q2,Q3,Q5a,Q5b,I,PbQ3,PiQ3eI,PiI,clQ1,clQ2,clQ3,clQ5a,clQ5b,clI,clPbQ3,clPiQ3eI,clPiI,filename,indim,varargin)
% Q1, Q2, Q3, Q5a, Q5b : QUADRANGLE
% I : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE
% PbQ3, PiQ3eI, PiI : POINT
% clQ1, clQ2, clQ3, clQ5a, clQ5b, clI, clPbQ3, clPiQ3eI, clPiI : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, max([getindim(Q1),getindim(Q2),getindim(Q3),getindim(Q5a),getindim(Q5b),getindim(I)]) by default)

if nargin<15
    clI = clQ3;
end
if nargin<16
    clPbQ3 = clQ3;
end
if nargin<17
    clPiQ3eI = clQ3;
end
if nargin<18
    clPiI = clI;
end
if nargin<20
    indim = max([getindim(Q1),getindim(Q2),getindim(Q3),getindim(Q5a),getindim(Q5b),getindim(I)]);
end

if ~iscell(PbQ3)
    PbQ3 = {PbQ3};
end
if isscalar(clPbQ3)
    clPbQ3 = repmat(clPbQ3,1,length(PbQ3));
end

if ~iscell(PiQ3eI)
    PiQ3eI = {PiQ3eI};
end
if isscalar(clPiQ3eI)
    clPiQ3eI = repmat(clPiQ3eI,1,length(PiQ3eI));
end

if ~iscell(PiI)
    PiI = {PiI};
end
if isscalar(clPiI)
    clPiI = repmat(clPiI,1,length(PiI));
end

G = GMSHFILE();
if nargin>=19 && ischar(filename)
    G = setfile(G,filename);
end

numpoints = 1:4;
numlines = 1:4;
numlineloop = 1;  
GQ5b = gmshfile(Q5b,clQ5b,numpoints,numlines,numlineloop,6,varargin{:});
G = G+GQ5b;

numpoints = numpoints(end)+(1:4);
numlines = numlines(end)+(1:4);
numlineloop = numlineloop(end)+1; 
GQ5a = gmshfile(Q5a,clQ5a,numpoints,numlines,numlineloop,5,varargin{:});
G = G+GQ5a;

PQ1 = getvertices(Q1);
numpoints = numpoints(end)+(1:4);
numlines = numlines(end)+(1:5);
numlineloop = numlineloop(end)+1; 
G = createpoints(G,PQ1,clQ1,numpoints);
G = createcontour(G,[9 2 10:12],numlines,numlineloop);
G = createplanesurface(G,numlineloop,1);
G = embedcurvesinsurface(G,[2 6],1);
if ischarin('recombine',varargin)
    G = recombinesurface(G,1);
end

PQ2 = getvertices(Q2);
numpoints = numpoints(end)+(1:4);
numlines = numlines(end)+(1:5);
numlineloop = numlineloop(end)+1;
G = createpoints(G,PQ2,clQ2,numpoints);
G = createcontour(G,[13 1 14:16],numlines,numlineloop);
G = createplanesurface(G,numlineloop,2);
G = embedcurvesinsurface(G,[4 8],2);
if ischarin('recombine',varargin)
    G = recombinesurface(G,2);
end

PQ3 = getvertices(Q3);
numpoints = numpoints(end)+(1:(length(PQ3)+length(PbQ3)));
numlines = numlines(end)+(1:(length(PQ3)+length(PbQ3)));
numlineloop = numlineloop(end)+(1:2);
G = createpoints(G,PQ3,clQ3,numpoints(1:3:end));
G = createpoints(G,PbQ3,clPbQ3,setdiff(numpoints,numpoints(1:3:end)));
G = createcontour(G,numpoints,numlines,numlineloop(1));

numpoints = numpoints(end)+(1:5);
numlines = numlines(end)+(1:4);
if isa(I,'DOMAIN') || isa(I,'QUADRANGLE')
    GI = gmshfile(I,clI,numpoints(1:end-1),numlines,numlineloop(2),4);
elseif isa(I,'CIRCLE') || isa(I,'ELLIPSE')
    GI = gmshfile(I,clI,numpoints(1),numpoints(2:end),numlines,numlineloop(2),4);
end
G = G+GI;
if ~isempty([PiI{:}])
    numpoints = numpoints(end)+(1:length(PiI));
    G = createpoints(G,PiI,clPiI,numpoints);
    G = embedpointsinsurface(G,numpoints,4);
end
if ischarin('recombine',varargin)
    G = recombinesurface(G,4);
end

G = createplanesurface(G,numlineloop(1:2),3);
G = embedcurvesinsurface(G,[12 17],3);
if ~isempty([PiQ3eI{:}])
    numpoints = numpoints(end)+(1:length(PiQ3eI));
    G = createpoints(G,PiQ3eI,clPiQ3eI,numpoints);
    G = embedpointsinsurface(G,numpoints,3);
end
if ischarin('recombine',varargin)
    G = recombinesurface(G,3);
end

varargin = delonlycharin('recombine',varargin);

n = max(nargout,1);
varargout = cell(1,n);
dim = max([getdim(Q1),getdim(Q2),getdim(Q3),getdim(Q5a),getdim(Q5b),getdim(I)]);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
