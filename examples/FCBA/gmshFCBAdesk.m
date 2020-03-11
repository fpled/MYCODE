function varargout = gmshFCBAdesk(Q1,Q2,Q3,Q5a,Q5b,I,LbQ3,GCiQ3eI,CiI,clQ1,clQ2,clQ3,clQ5a,clQ5b,clI,clLbQ3,clCiQ3eI,clCiI,filename,indim,varargin)
% function varargout = gmshFCBAdesk(Q1,Q2,Q3,Q5a,Q5b,I,LbQ3,CiQ3eI,CiI,clQ1,clQ2,clQ3,clQ5a,clQ5b,clI,clLbQ3,clCiQ3eI,clCiI,filename,indim,varargin)
% Q1, Q2, Q3, Q5a, Q5b : QUADRANGLE
% I : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE
% LbQ3 : LIGNE
% CiQ3eI, CiI : CIRCLE
% clQ1, clQ2, clQ3, clQ5a, clQ5b, clI, clLbQ3, clCiQ3eI, clCiI : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, max([getindim(Q1),getindim(Q2),getindim(Q3),getindim(Q5a),getindim(Q5b),getindim(I)]) by default)

if nargin<15
    clI = clQ3;
end
if nargin<16
    clLbQ3 = clQ3;
end
if nargin<17
    clCiQ3eI = clQ3;
end
if nargin<18
    clCiI = clI;
end
if nargin<20
    indim = max([getindim(Q1),getindim(Q2),getindim(Q3),getindim(Q5a),getindim(Q5b),getindim(I)]);
end

if ~iscell(LbQ3)
    LbQ3 = {LbQ3};
end
if length(clLbQ3)==1
    clLbQ3 = repmat(clLbQ3,1,length(LbQ3));
end

G = GMSHFILE();
if nargin>=19 && ischar(filename)
    G = setfile(G,filename);
end

numpoints = 1:4;
numlines = 1:4;
numlineloop = 1;  
GQ5b = gmshfile(Q5b,clQ5b,numpoints,numlines,numlineloop,8,varargin{:});
G = G+GQ5b;

numpoints = numpoints(end)+(1:4);
numlines = numlines(end)+(1:4);
numlineloop = numlineloop(end)+1; 
GQ5a = gmshfile(Q5a,clQ5a,numpoints,numlines,numlineloop,7,varargin{:});
G = G+GQ5a;

PQ1 = getvertices(Q1);
numpoints = numpoints(end)+(1:4);
numlines = numlines(end)+(1:5);
numlineloop = numlineloop(end)+1; 
G = createpoints(G,PQ1,clQ1,numpoints);
G = createcontour(G,[9 2 10:12],numlines,numlineloop);
G = createplanesurface(G,numlineloop,1);
G = embedlinesinsurface(G,[2 6],1);
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
G = embedlinesinsurface(G,[4 8],2);
if ischarin('recombine',varargin)
    G = recombinesurface(G,2);
end

PQ3 = getvertices(Q3);
numpoints = numpoints(end)+(1:(length(PQ3)+2*length(LbQ3)));
numlines = numlines(end)+(1:(length(PQ3)+2*length(LbQ3)));
numlineloop = numlineloop(end)+(1:4);
G = createpoints(G,PQ3,clQ3,numpoints(1:5:end));
npoints = setdiff(numpoints,numpoints(1:5:end));
for k=1:length(LbQ3)
    PbQ3 = getvertices(LbQ3{k});
    G = createpoints(G,PbQ3,clLbQ3(k),npoints([2*k-1 2*k]));
end
G = createcontour(G,numpoints,numlines,numlineloop(1));

numpoints = numpoints(end)+(1:5);
numlines = numlines(end)+(1:4);
GCiI = gmshfile(CiI,clCiI,numpoints(1),numpoints(2:end),numlines,numlineloop(2),4,varargin{:});
G = G+GCiI;

numpoints = numpoints(end)+(1:5);
numlines = numlines(end)+(1:4);
if isa(I,'DOMAIN') || isa(I,'QUADRANGLE')
    GI = gmshfile(I,clI,numpoints(1:end-1),numlines,numlineloop(3));
elseif isa(I,'CIRCLE') || isa(I,'ELLIPSE')
    GI = gmshfile(I,clI,numpoints(1),numpoints(2:end),numlines,numlineloop(3));
end
G = G+GI;
G = createplanesurface(G,numlineloop(2:3),5);
if ischarin('recombine',varargin)
    G = recombinesurface(G,5);
end

numpoints = numpoints(end)+(1:5);
numlines = numlines(end)+(1:4);
GCiQ3eI = gmshfile(GCiQ3eI,clCiQ3eI,numpoints(1),numpoints(2:end),numlines,numlineloop(4),6,varargin{:});
G = G+GCiQ3eI;

G = createplanesurface(G,numlineloop([1 3 4]),3);
G = embedlinesinsurface(G,[12 17],3);
if ischarin('recombine',varargin)
    G = recombinesurface(G,3);
end

varargin = delonlycharin('recombine',varargin);

n=max(nargout,1);
varargout = cell(1,n);
dim = max([getdim(Q1),getdim(Q2),getdim(Q3),getdim(Q5a),getdim(Q5b),getdim(I)]);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
