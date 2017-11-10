function varargout = gmshFCBAdesk3D2(D,Qa,Qb,clD,clQa,clQb,filename,indim,varargin)
% function varargout = gmshFCBAdesk3D2(D,Qa,Qb,clD,clQa,clQb,filename,indim,varargin)
% D : DOMAIN
% Qa, Qb : QUADRANGLE
% clD, clQa, clQb : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getdim(D) by default)

if nargin<5
    clQa = clD;
end
if nargin<6
    clQb = clD;
end
if nargin<8
    indim = getdim(D);
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

G = createlines(G,[[2 9];[9 10];[10 11];[11 12];[12 9];[10 3];[3 7];[7 6]],5:12);
G = createlines(G,[[3 4];[4 8];[8 7]],13:15);

G = createlineloop(G,[13:15 -11],3);
G = createplanesurface(G,3,3);

G = createlines(G,[[1 4];[8 5]],16:17);
G = createlineloop(G,[16 14 17 4],4);
G = createplanesurface(G,4,4);

G = createlineloop(G,[1 5 6 10 13 -16],5);
G = createplanesurface(G,5,5);

G = createlineloop(G,[3 -17 15 12],6);
G = createplanesurface(G,6,6);

G = createlineloop(G,6:9,7);
G = createplanesurface(G,7,7);

G = createlines(G,[[13 14];[14 15];[15 16];[16 13]],18:21);
G = createlineloop(G,18:21,8);
G = createplanesurface(G,8,8);

G = createlineloop(G,[2 -12:-10 7:9 -5],2);
G = createplanesurface(G,[2,8],2);

G = createsurfaceloop(G,1:8,1);
G = createvolume(G,1,1);

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(D):-1:getdim(D)-n+1);
