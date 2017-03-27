function varargout = gmshFCBAdesk12(Q,PL,PbQ,clQ,clPL,clPbQ,filename,indim,varargin)
% function varargout = gmshFCBAdesk12(Q,PL,PbQ,clQ,clPL,clPbQ,filename,indim,varargin)
% Q : QUADRANGLE
% PL, PbQ : POINT
% clQ, clPL : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getdim(D) by default)

if nargin<5
    clPL = clQ;
end
if nargin<6
    clPbQ = clQ;
end
if nargin<8
    indim = getdim(Q);
end

if ~iscell(PL)
    PL = {PL};
end
if length(clPL)==1
    clPL = repmat(clPL,1,length(PL));
end

G = GMSHFILE();
if nargin>=7 && ischar(filename)
    G = setfile(G,filename);
end

numpoints = 1:length(PbQ);
numlines = 1:length(PbQ);
G = createpoints(G,PbQ,clPbQ,numpoints);
G = createcontour(G,numpoints,numlines,1);

G = createplanesurface(G,1,1);
numpoints = numpoints(end)+(1:3);
numlines = numlines(end)+(1:2);
seg = [2 numpoints(3);numpoints(1) numpoints(2)];
G = createpoints(G,PL,clPL,numpoints);
G = createlines(G,seg,numlines);
G = embedlinesinsurface(G,numlines,1);

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(Q):-1:getdim(Q)-n+1);
