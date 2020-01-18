function varargout = gmshFCBAdesk3D(D1,D2,D3,D5a,D5b,I,LbD3,CiD3eI,CiI,clD1,clD2,clD3,clD5a,clD5b,clI,clLbD3,clCiD3eI,clCiI,filename,indim,varargin)
% function varargout = gmshFCBAdesk3D(D1,D2,D3,D5a,D5b,I,LbD3,CiD3eI,CiI,clD1,clD2,clD3,clD5a,clD5b,clI,clLbD3,clCiD3eI,clCiI,filename,indim,varargin)
% D1, D2, D3, D5a, D5b : DOMAIN
% I : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE
% LbD3 : LIGNE
% CiD3eI, CiI : CIRCLE
% clD1, clD2, clD3, clD5a, clD5b, clI, clLbD3, clCiD3eI, clCiI : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, max([getindim(D1),getindim(D2),getindim(D3),getindim(D5a),getindim(D5b),getindim(I)]) by default)

if nargin<15
    clI = clD3;
end
if nargin<16
    clLbD3 = clD3;
end
if nargin<17
    clCiD3eI = clD3;
end
if nargin<18
    clCiI = clI;
end
if nargin<20
    indim = max([getindim(D1),getindim(D2),getindim(D3),getindim(D5a),getindim(D5b),getindim(I)]);
end

if ~iscell(LbD3)
    LbD3 = {LbD3};
end
if length(clLbD3)==1
    clLbD3 = repmat(clLbD3,1,length(LbD3));
end

G = GMSHFILE();
if nargin>=19 && ischar(filename)
    G = setfile(G,filename);
end

PD5b = getvertices(D5b);
numpoints = 1:8;
numlines = 1:12;
numlineloop = 1:6;
numsurface = 1:6;
G = createpoints(G,PD5b,clD5b,numpoints);
G = createcontour(G,numpoints([1 2 6 5]),numlines(1:4),numlineloop(1));
G = createplanesurface(G,numlineloop(1),numsurface(1));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(1));
end

G = createlines(G,numpoints([[2 3];[3 7];[7 6]]),numlines(5:7));
G = createlineloop(G,[numlines(5:7) -numlines(2)],numlineloop(2));
G = createplanesurface(G,numlineloop(2),numsurface(2));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(2));
end

G = createlines(G,numpoints([[3 4];[4 8];[8 7]]),numlines(8:10));
G = createlineloop(G,[numlines(8:10) -numlines(6)],numlineloop(3));
G = createplanesurface(G,numlineloop(3),numsurface(3));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(3));
end

G = createlines(G,numpoints([[4 1];[5 8]]),numlines(11:12));
G = createlineloop(G,[numlines(4) -numlines(11) numlines(9) -numlines(12)],numlineloop(4));
G = createplanesurface(G,numlineloop(4),numsurface(4));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(4));
end

G = createlineloop(G,numlines([1 5 8 11]),numlineloop(5));
G = createplanesurface(G,numlineloop(5),numsurface(5));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(5));
end

G = createlineloop(G,numlines([3 12 10 7]),numlineloop(6));
G = createplanesurface(G,numlineloop(6),numsurface(6));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(6));
end

G = createsurfaceloop(G,numsurface(1:6),5);
G = createvolume(G,5,5);

PD5a = getvertices(D5a);
numpoints = numpoints(end)+(1:8);
numlines = numlines(end)+(1:12);
numlineloop = numlineloop(end)+(1:6);
numsurface = numsurface(end)+(1:6);
G = createpoints(G,PD5a,clD5a,numpoints);
G = createcontour(G,numpoints([1 2 6 5]),numlines(1:4),numlineloop(1));
G = createplanesurface(G,numlineloop(1),numsurface(1));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(1));
end

G = createlines(G,numpoints([[2 3];[3 7];[7 6]]),numlines(5:7));
G = createlineloop(G,[numlines(5:7) -numlines(2)],numlineloop(2));
G = createplanesurface(G,numlineloop(2),numsurface(2));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(2));
end

G = createlines(G,numpoints([[3 4];[4 8];[8 7]]),numlines(8:10));
G = createlineloop(G,[numlines(8:10) -numlines(6)],numlineloop(3));
G = createplanesurface(G,numlineloop(3),numsurface(3));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(3));
end

G = createlines(G,numpoints([[4 1];[5 8]]),numlines(11:12));
G = createlineloop(G,[numlines(4) -numlines(11) numlines(9) -numlines(12)],numlineloop(4));
G = createplanesurface(G,numlineloop(4),numsurface(4));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(4));
end

G = createlineloop(G,numlines([1 5 8 11]),numlineloop(5));
G = createplanesurface(G,numlineloop(5),numsurface(5));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(5));
end

G = createlineloop(G,numlines([3 12 10 7]),numlineloop(6));
G = createplanesurface(G,numlineloop(6),numsurface(6));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(6));
end

G = createsurfaceloop(G,numsurface(1:6),4);
G = createvolume(G,4,4);

PD1 = getvertices(D1);
numpoints = numpoints(end)+(1:8);
numlines = numlines(end)+(1:13);
numlineloop = numlineloop(end)+(1:6);
numsurface = numsurface(end)+(1:6);
G = createpoints(G,PD1,clD1,numpoints);
G = createcontour(G,numpoints([1 2 6 5]),numlines(1:4),numlineloop(1));
G = createplanesurface(G,numlineloop(1),numsurface(1));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(1));
end

G = createlines(G,numpoints([[2 3];[3 7];[7 6]]),numlines(5:7));
G = createlineloop(G,[numlines(5:7) -numlines(2)],numlineloop(2));
G = createplanesurface(G,numlineloop(2),numsurface(2));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(2));
end

G = createlines(G,numpoints([[3 4];[4 8];[8 7]]),numlines(8:10));
G = createlineloop(G,[numlines(8:10) -numlines(6)],numlineloop(3));
G = createplanesurface(G,numlineloop(3),numsurface(3));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(3));
end

G = createlines(G,[[numpoints(4) 3];[2 numpoints(1)];numpoints([5 8])],numlines(11:13));
G = createlineloop(G,[numlines([1 5 8 11]) -5 numlines(12)],numlineloop(5));
G = createplanesurface(G,numlineloop(5),numsurface(5));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(5));
end

G = createlineloop(G,numlines([3 13 10 7]),numlineloop(6));
G = createplanesurface(G,numlineloop(6),numsurface(6));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(6));
end

G = createlineloop(G,[-numlines(12) 2 -7 -6 -numlines(11) numlines(9) -numlines(13) numlines(4)],numlineloop(4));
G = createplanesurface(G,[numlineloop(4) 8],numsurface(4));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(4));
end

G = createsurfaceloop(G,[numsurface(1:6) 2 8],1);
G = createvolume(G,1,1);

PD2 = getvertices(D2);
numpoints = numpoints(end)+(1:8);
numlines = numlines(end)+(1:13);
numlineloop = numlineloop(end)+(1:6);
numsurface = numsurface(end)+(1:6);
G = createpoints(G,PD2,clD2,numpoints(1:8));
G = createcontour(G,numpoints([1 2 6 5]),numlines(1:4),numlineloop(1));
G = createplanesurface(G,numlineloop(1),numsurface(1));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(1));
end

G = createlines(G,[[numpoints(2) 1];[4 numpoints(3)];numpoints([[3 7];[7 6]])],numlines(5:8));
G = createlines(G,numpoints([[3 4];[4 8];[8 7]]),numlines(9:11));
G = createlineloop(G,[numlines(9:11) -numlines(7)],numlineloop(3));
G = createplanesurface(G,numlineloop(3),numsurface(3));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(3));
end

G = createlines(G,numpoints([[1 4];[8 5]]),numlines(12:13));
G = createlineloop(G,numlines([12 10 13 4]),numlineloop(4));
G = createplanesurface(G,numlineloop(4),numsurface(4));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(4));
end

G = createlineloop(G,[numlines([1 5]) -11 numlines(6) numlines(9) -numlines(12)],numlineloop(5));
G = createplanesurface(G,numlineloop(5),numsurface(5));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(5));
end

G = createlineloop(G,[numlines(3) -numlines(13) numlines([11 8])],numlineloop(6));
G = createplanesurface(G,numlineloop(6),numsurface(6));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(6));
end

G = createlineloop(G,[[numlines(2) -numlines(8):-numlines(6)] 9 -12 4 -numlines(5)],numlineloop(2));
G = createplanesurface(G,[numlineloop(2) 10],numsurface(2));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(2));
end

G = createsurfaceloop(G,[numsurface(1:6) 4 10],2);
G = createvolume(G,2,2);

PD3 = getvertices(D3);
numpoints = numpoints(end)+(1:(length(PD3)+2*length(LbD3)));
G = createpoints(G,PD3,clD3,numpoints(1:length(PD3)));
for k=1:length(LbD3)
    PbD3 = getvertices(LbD3{k});
    G = createpoints(G,PbD3,clLbD3(k),numpoints(length(PD3)+[2*k-1 2*k]));
end

numlines = numlines(end)+(1:60);
numlineloop = numlineloop(end)+(1:25);
numsurface = numsurface(end)+(1:25);
G = createcontour(G,numpoints([1 9 10 5]),numlines(1:4),numlineloop(1));
G = createplanesurface(G,numlineloop(1),numsurface(1));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(1));
end

G = createlines(G,numpoints([[9 11];[11 12];[12 10]]),numlines(5:7));
G = createlineloop(G,[numlines(5:7) -numlines(2)],numlineloop(2));
G = createplanesurface(G,numlineloop(2),numsurface(2));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(2));
end

G = createlines(G,numpoints([[11 13];[13 14];[14 12]]),numlines(8:10));
G = createlineloop(G,[numlines(8:10) -numlines(6)],numlineloop(3));
G = createplanesurface(G,numlineloop(3),numsurface(3));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(3));
end

G = createlines(G,numpoints([[13 15];[15 16];[16 14]]),numlines(11:13));
G = createlineloop(G,[numlines(11:13) -numlines(9)],numlineloop(4));
G = createplanesurface(G,numlineloop(4),numsurface(4));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(4));
end

G = createlines(G,numpoints([[15 2];[2 6];[6 16]]),numlines(14:16));
G = createlineloop(G,[numlines(14:16) -numlines(12)],numlineloop(5));
G = createplanesurface(G,numlineloop(5),numsurface(5));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(5));
end

G = createlines(G,numpoints([[2 17];[17 18];[18 6]]),numlines(17:19));
G = createlineloop(G,[numlines(17:19) -numlines(15)],numlineloop(6));
G = createplanesurface(G,numlineloop(6),numsurface(6));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(6));
end

G = createlines(G,numpoints([[17 19];[19 20];[20 18]]),numlines(20:22));
G = createlineloop(G,[numlines(20:22) -numlines(18)],numlineloop(7));
G = createplanesurface(G,numlineloop(7),numsurface(7));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(7));
end

G = createlines(G,numpoints([[19 21];[21 22];[22 20]]),numlines(23:25));
G = createlineloop(G,[numlines(23:25) -numlines(21)],numlineloop(8));
G = createplanesurface(G,numlineloop(8),numsurface(8));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(8));
end

G = createlines(G,numpoints([[21 23];[23 24];[24 22]]),numlines(26:28));
G = createlineloop(G,[numlines(26:28) -numlines(24)],numlineloop(9));
G = createplanesurface(G,numlineloop(9),numsurface(9));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(9));
end

G = createlines(G,numpoints([[23 3];[3 7];[7 24]]),numlines(29:31));
G = createlineloop(G,[numlines(29:31) -numlines(27)],numlineloop(10));
G = createplanesurface(G,numlineloop(10),numsurface(10));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(10));
end

G = createlines(G,numpoints([[3 25];[25 26];[26 7]]),numlines(32:34));
G = createlineloop(G,[numlines(32:34) -numlines(30)],numlineloop(11));
G = createplanesurface(G,numlineloop(11),numsurface(11));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(11));
end

G = createlines(G,numpoints([[25 27];[27 28];[28 26]]),numlines(35:37));
G = createlineloop(G,[numlines(35:37) -numlines(33)],numlineloop(12));
G = createplanesurface(G,numlineloop(12),numsurface(12));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(12));
end

G = createlines(G,numpoints([[27 29];[29 30];[30 28]]),numlines(38:40));
G = createlineloop(G,[numlines(38:40) -numlines(36)],numlineloop(13));
G = createplanesurface(G,numlineloop(13),numsurface(13));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(13));
end

G = createlines(G,numpoints([[29 31];[31 32];[32 30]]),numlines(41:43));
G = createlineloop(G,[numlines(41:43) -numlines(39)],numlineloop(14));
G = createplanesurface(G,numlineloop(14),numsurface(14));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(14));
end

G = createlines(G,numpoints([[31 4];[4 8];[8 32]]),numlines(44:46));
G = createlineloop(G,[numlines(44:46) -numlines(42)],numlineloop(15));
G = createplanesurface(G,numlineloop(15),numsurface(15));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(15));
end

G = createlines(G,numpoints([[4 33];[33 34];[34 8]]),numlines(47:49));
G = createlineloop(G,[numlines(47:49) -numlines(45)],numlineloop(16));
G = createplanesurface(G,numlineloop(16),numsurface(16));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(16));
end

G = createlines(G,numpoints([[33 35];[35 36];[36 34]]),numlines(50:52));
G = createlineloop(G,[numlines(50:52) -numlines(48)],numlineloop(17));
G = createplanesurface(G,numlineloop(17),numsurface(17));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(17));
end

G = createlines(G,numpoints([[35 37];[37 38];[38 36]]),numlines(53:55));
G = createlineloop(G,[numlines(53:55) -numlines(51)],numlineloop(18));
G = createplanesurface(G,numlineloop(18),numsurface(18));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(18));
end

G = createlines(G,numpoints([[37 39];[39 40];[40 38]]),numlines(56:58));
G = createlineloop(G,[numlines(56:58) -numlines(54)],numlineloop(19));
G = createplanesurface(G,numlineloop(19),numsurface(19));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(19));
end

G = createlines(G,numpoints([[39 1];[5 40]]),numlines(59:60));
G = createlineloop(G,[numlines(59) -numlines(4) numlines(60) -numlines(57)],numlineloop(20));
G = createplanesurface(G,numlineloop(20),numsurface(20));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(20));
end

G = createlineloop(G,numlines([1 5:3:59]),numlineloop(21));
G = createplanesurface(G,[18 24 numlineloop(21)],numsurface(21));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(21));
end

numpoints = numpoints(end)+(1:5);
numberlines = numlines(end)+(1:4);
GCiI = gmshfile(CiI,clCiI,numpoints(1),numpoints(2:end),numberlines,numlineloop(22),numsurface(22),varargin{:});
G = G+GCiI;

numpoints = numpoints(end)+(1:5);
numberlines = numberlines(end)+(1:4);
if isa(I,'DOMAIN') || isa(I,'QUADRANGLE')
    GI = gmshfile(I,clI,numpoints(1:end-1),numberlines,numlineloop(23));
elseif isa(I,'CIRCLE') || isa(I,'ELLIPSE')
    GI = gmshfile(I,clI,numpoints(1),numpoints(2:end),numberlines,numlineloop(23));
end
G = G+GI;
G = createplanesurface(G,numlineloop([22 23]),numsurface(23));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(23));
end

numpoints = numpoints(end)+(1:5);
numberlines = numberlines(end)+(1:4);
GCiDeI = gmshfile(CiD3eI,clCiD3eI,numpoints(1),numpoints(2:end),numberlines,numlineloop(24),numsurface(24),varargin{:});
G = G+GCiDeI;

G = createlineloop(G,numlines([60 58:-3:7 3]),numlineloop(25));
G = createplanesurface(G,numlineloop([23 24 25]),numsurface(25));
if ischarin('recombine',varargin)
    G = recombinesurface(G,numsurface(25));
end

G = createsurfaceloop(G,[numsurface(1:25) 18 24],3);
G = createvolume(G,3,3);

n=max(nargout,1);
varargout = cell(1,n);
dim = max([getdim(D1),getdim(D2),getdim(D3),getdim(D5a),getdim(D5b),getdim(I)]);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
