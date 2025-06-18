function varargout = gmshFCBAdesk3D1(D,Qa,Qb,clD,clQa,clQb,filename,indim,varargin)
% function varargout = gmshFCBAdesk3D1(D,Qa,Qb,clD,clQa,clQb,filename,indim,varargin)
% D : DOMAIN
% Qa, Qb : QUADRANGLE
% clD, clQa, clQb : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getindim(D) by default)

if nargin<5
    clQa = clD;
end
if nargin<6
    clQb = clD;
end
if nargin<8
    indim = getindim(D);
end

G = GMSHFILE();
if nargin>=7 && ischar(filename)
    G = setfile(G,filename);
end

PD = getvertices(D);
PQa = getvertices(Qa);
PQb = getvertices(Qb);

G = createpoints(G,PD,clD,1:8);
G = createpoints(G,PQa,clQa,13:16);
G = createpoints(G,PQb,clQb,9:12);

G = createcontour(G,[1 2 6 5],1:4,1);
G = createplanesurface(G,1,1);
if ischarin('recombine',varargin)
    G = recombinesurface(G,1);
end

G = createlines(G,[[2 3];[3 7];[7 6]],5:7);
G = createcurveloop(G,[5:7 -2],2);
G = createplanesurface(G,2,2);
if ischarin('recombine',varargin)
    G = recombinesurface(G,2);
end

G = createlines(G,[[3 4];[4 8];[8 7]],8:10);
G = createcurveloop(G,[8:10 -6],3);
G = createplanesurface(G,3,3);
if ischarin('recombine',varargin)
    G = recombinesurface(G,3);
end

G = createlines(G,[[4 9];[9 10];[10 11];[11 12];[12 9];[10 1];[5 8]],11:17);
G = createlines(G,[[13 14];[14 15];[15 16];[16 13]],18:21);

G = createcurveloop(G,[1 5 8 11 12 16],5);
G = createplanesurface(G,5,5);
if ischarin('recombine',varargin)
    G = recombinesurface(G,5);
end

G = createcurveloop(G,[3 17 10 7],6);
G = createplanesurface(G,6,6);
if ischarin('recombine',varargin)
    G = recombinesurface(G,6);
end

G = createcurveloop(G,12:15,7);
G = createplanesurface(G,7,7);
if ischarin('recombine',varargin)
    G = recombinesurface(G,7);
end

G = createcurveloop(G,18:21,8);
G = createplanesurface(G,8,8);
if ischarin('recombine',varargin)
    G = recombinesurface(G,8);
end

G = createcurveloop(G,[-16 13:15 -11 9 -17 4],4);
G = createplanesurface(G,[4,8],4);
if ischarin('recombine',varargin)
    G = recombinesurface(G,4);
end

G = createsurfaceloop(G,1:8,1);
G = createvolume(G,1,1);

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(D):-1:getdim(D)-n+1,varargin{:});
