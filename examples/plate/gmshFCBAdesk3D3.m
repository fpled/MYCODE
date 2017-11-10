function varargout = gmshFCBAdesk3D3(D,I,Q1,Q2,LbD,CiDeI,CiI,clD,clI,clQ1,clQ2,clLbD,clCiDeI,clCiI,filename,indim,varargin)
% function varargout = gmshFCBAdesk33D(D,I,Q1,Q2,LbD,CiDeI,CiI,clD,clI,clQ1,clQ2,clLbD,clCiDeI,clCiI,filename,indim,varargin)
% D : DOMAIN
% I : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE
% Q1, Q2 : QUADRANGLE
% LbD : LIGNE
% CbQ, CiQeI, CiI : CIRCLE
% clD, clI, clQ1, clQ2, clLbD, clCiDeI, clCiI : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getdim(D) by default)

if nargin<9
    clI = clD;
end
if nargin<10
    clQ1 = clD;
end
if nargin<11
    clQ2 = clD;
end
if nargin<12
    clLbD = clD;
end
if nargin<13
    clCiDeI = clD;
end
if nargin<14
    clCiI = clI;
end
if nargin<16
    indim = getdim(D);
end

if ~iscell(LbD)
    LbD = {LbD};
end
if length(clLbD)==1
    clLbD = repmat(clLbD,1,length(LbD));
end

G = GMSHFILE();
if nargin>=15 && ischar(filename)
    G = setfile(G,filename);
end

PD = getvertices(D);
numpoints = 1:8;
G = createpoints(G,PD,clD,numpoints);
for k=1:length(LbD)
    numpoints = numpoints(end)+(1:2);
    PbD = getvertices(LbD{k});
    G = createpoints(G,PbD,clLbD(k),numpoints);
end

numlines = 1:4;
G = createcontour(G,[1 9 10 5],1:4,1);
G = createplanesurface(G,1,1);

numlines = numlines(end)+(1:3);
G = createlines(G,[[9 11];[11 12];[12 10]],5:7);
G = createlineloop(G,[5:7 -2],2);
G = createplanesurface(G,2,2);

numlines = numlines(end)+(1:3);
G = createlines(G,[[11 13];[13 14];[14 12]],8:10);
G = createlineloop(G,[8:10 -6],3);
G = createplanesurface(G,3,3);

numlines = numlines(end)+(1:3);
G = createlines(G,[[13 15];[15 16];[16 14]],11:13);
G = createlineloop(G,[11:13 -9],4);
G = createplanesurface(G,4,4);

numlines = numlines(end)+(1:3);
G = createlines(G,[[15 2];[2 6];[6 16]],14:16);
G = createlineloop(G,[14:16 -12],5);
G = createplanesurface(G,5,5);

numlines = numlines(end)+(1:3);
G = createlines(G,[[2 17];[17 18];[18 6]],17:19);
G = createlineloop(G,[17:19 -15],6);
G = createplanesurface(G,6,6);

numlines = numlines(end)+(1:3);
G = createlines(G,[[17 19];[19 20];[20 18]],20:22);
G = createlineloop(G,[20:22 -18],7);
G = createplanesurface(G,7,7);

numlines = numlines(end)+(1:3);
G = createlines(G,[[19 21];[21 22];[22 20]],23:25);
G = createlineloop(G,[23:25 -21],8);
G = createplanesurface(G,8,8);

numlines = numlines(end)+(1:3);
G = createlines(G,[[21 23];[23 24];[24 22]],26:28);
G = createlineloop(G,[26:28 -24],9);
G = createplanesurface(G,9,9);

numlines = numlines(end)+(1:3);
G = createlines(G,[[23 3];[3 7];[7 24]],29:31);
G = createlineloop(G,[29:31 -27],10);
G = createplanesurface(G,10,10);

numlines = numlines(end)+(1:3);
G = createlines(G,[[3 25];[25 26];[26 7]],32:34);
G = createlineloop(G,[32:34 -30],11);
G = createplanesurface(G,11,11);

numlines = numlines(end)+(1:3);
G = createlines(G,[[25 27];[27 28];[28 26]],35:37);
G = createlineloop(G,[35:37 -33],12);
G = createplanesurface(G,12,12);

numlines = numlines(end)+(1:3);
G = createlines(G,[[27 29];[29 30];[30 28]],38:40);
G = createlineloop(G,[38:40 -36],13);
G = createplanesurface(G,13,13);

numlines = numlines(end)+(1:3);
G = createlines(G,[[29 31];[31 32];[32 30]],41:43);
G = createlineloop(G,[41:43 -39],14);
G = createplanesurface(G,14,14);

numlines = numlines(end)+(1:3);
G = createlines(G,[[31 4];[4 8];[8 32]],44:46);
G = createlineloop(G,[44:46 -42],15);
G = createplanesurface(G,15,15);

numlines = numlines(end)+(1:3);
G = createlines(G,[[4 33];[33 34];[34 8]],47:49);
G = createlineloop(G,[47:49 -45],16);
G = createplanesurface(G,16,16);

numlines = numlines(end)+(1:3);
G = createlines(G,[[33 35];[35 36];[36 34]],50:52);
G = createlineloop(G,[50:52 -48],17);
G = createplanesurface(G,17,17);

numlines = numlines(end)+(1:3);
G = createlines(G,[[35 37];[37 38];[38 36]],53:55);
G = createlineloop(G,[53:55 -51],18);
G = createplanesurface(G,18,18);

numlines = numlines(end)+(1:3);
G = createlines(G,[[37 39];[39 40];[40 38]],56:58);
G = createlineloop(G,[56:58 -54],19);
G = createplanesurface(G,19,19);

numlines = numlines(end)+(1:2);
G = createlines(G,[[39 1];[5 40]],59:60);
G = createlineloop(G,[59 -4 60 -57],20);
G = createplanesurface(G,20,20);

numpoints = numpoints(end)+(1:4);
numlines = numlines(end)+(1:4);
GQ1 = gmshfile(Q1,clQ1,numpoints,numlines,21);
G = G+GQ1;
G = createplanesurface(G,21,21);

numpoints = numpoints(end)+(1:4);
numlines = numlines(end)+(1:4);
GQ2 = gmshfile(Q2,clQ2,numpoints,numlines,22);
G = G+GQ2;
G = createplanesurface(G,22,22);

G = createlineloop(G,[1 5:3:59],23);
G = createplanesurface(G,[21,22,23],23);

numpoints = numpoints(end)+(1:5);
numlines = numlines(end)+(1:4);
GCiI = gmshfile(CiI,clCiI,numpoints(1),numpoints(2:end),numlines,24);
G = G+GCiI;
G = createplanesurface(G,24,24);

numpoints = numpoints(end)+(1:5);
numlines = numlines(end)+(1:4);
if isa(I,'DOMAIN') || isa(I,'QUADRANGLE')
    GI = gmshfile(I,clI,numpoints(1:end-1),numlines,25);
elseif isa(I,'CIRCLE') || isa(I,'ELLIPSE')
    GI = gmshfile(I,clI,numpoints(1),numpoints(2:end),numlines,25);
end
G = G+GI;
G = createplanesurface(G,[24 25],25);

numpoints = numpoints(end)+(1:5);
numlines = numlines(end)+(1:4);
GCiDeI = gmshfile(CiDeI,clCiDeI,numpoints(1),numpoints(2:end),numlines,26);
G = G+GCiDeI;
G = createplanesurface(G,26,26);

G = createlineloop(G,[60 58:-3:7 3],28);
G = createplanesurface(G,[25 26 28],27);

G = createsurfaceloop(G,1:27,1);
G = createvolume(G,1,1);

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(D):-1:getdim(D)-n+1);
