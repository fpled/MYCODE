function varargout = gmshFCBAtablecirc(C,I,Pi,Pe,clC,clI,clPi,clPe,filename,indim,varargin)
% function varargout = gmshFCBAtablecirc(C,I,Pi,Pe,clC,clI,clPi,clPe,filename,indim)
% C : CIRCLE
% I : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE
% Pi, Pe : POINT
% clD, clI, clPi, clPe : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getdim(D) by default)

if nargin<10
    indim = getdim(C);
end
if nargin<8
    clPe = clC;
end
if nargin<7
    clPi = clC;
end
if nargin<6
    clI = clC;
end

if ~iscell(Pe)
    Pe = {Pe};
end
if length(clPe)==1
    clPe = repmat(clPe,1,length(Pe));
end

G = GMSHFILE();
if nargin>=9 && ischar(filename)
    G = setfile(G,filename);
end

numcenter = 1;
numpoints = 1+(1:length(Pe));
numlines = 1:length(Pe);
numlineloop = 5;
% PC = getvertices(C);
Pc = double(getcenter(C));
G = createpoint(G,Pc,clC,numcenter);
G = createpoints(G,Pe,clPe,numpoints);
G = createcirclecontour(G,numcenter,numpoints,numlines,numlineloop);

if ~iscell(I)
    I = {I};
end
if length(clI)==1
    clI = repmat(clI,1,length(I));
end

if ~iscell(Pi)
    Pi = {Pi};
end
if length(clPi)==1
    clPi = repmat(clPi,1,length(Pi));
end

numpoints = numpoints(end)+(1:5);
numlines = numlines(end)+(1:5);
numlineloop = 1:6;
for j=1:length(I)
    numlineloop = [numlineloop,-numlines(1:end-1)];
    if isa(I{j},'DOMAIN') || isa(I{j},'QUADRANGLE')
        GI = gmshfile(I{j},clI(j),numpoints(1:end-1),numlines(1:end-1),numlines(end));
    elseif isa(I{j},'CIRCLE') || isa(I{j},'ELLIPSE')
        GI = gmshfile(I{j},clI(j),numpoints(1),numpoints(2:end),numlines(1:end-1),numlines(end));
    end
    G = G+GI;
    G = createplanesurface(G,numlines(end),j+1);
    if ischarin('recombine',varargin)
        G = recombinesurface(G,j+1);
    end
    numpoints = numpoints+5;
    numlines = numlines+5;
end
G = createlineloop(G,numlineloop,numlines(end));
G = createplanesurface(G,numlines(end),1);
numberembeddedpoints = numpoints(end)+(1:length(Pi));
G = createpoints(G,Pi,clPi,numberembeddedpoints);
G = embedpointsinsurface(G,numberembeddedpoints,1);
if ischarin('recombine',varargin)
    G = recombinesurface(G,1);
end

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(C):-1:getdim(C)-n+1);
