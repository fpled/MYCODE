function Apatchsep = calc_patchmatrix_sep(Spart,PC,typepatch,varargin)
%function Apatch = calc_patchmatrix_sep(Spart,PC,typepatch,varargin)
% ('beh' pour diffusion, 'D' pour Dirichlet et 'N' pour Neumann)
%calc_patchmatrix_sep(Spart,PC,typepatch,num,'beh','k',ktp)
if(nargin==2)
    error('Veillez entre les type de patch : ''behavior'' pour diffusion \n''D'' pour Dirichlet et ''N'' pour Neumann  ');
end

if(nargin>=3)
    if strcmp(typepatch,'behavior')
        ktp = getcharin('k',varargin);
        
        if(~isa(ktp,'cell'))
            fprintf('pas prevu')
        end
        %
        for k=2:size(Spart,2)+1
            L1 = setphi(PCTPMATRIX(PC,1),getphi(ktp{k}{1}));
            ak = BILINFORM(1,1,getphi0(ktp{k}{1}),0);
            Apatch{k} = ak{Spart{k}}(:,:)*L1;
            for i=2:getm(ktp{k})
                Li = setphi(PCTPMATRIX(PC,1),getphi(ktp{k}{i}));
                ak = BILINFORM(1,1,getphi0(ktp{k}{i}),0);
                Apatch{k} = Apatch{k}+ak{Spart{k}}(:,:)*Li;
            end
            Apatch{k} = calc_ximasse(Apatch{k});
            Apatchsep{k} = SEPMATRIX(Apatch{k});
        end
    end
    
    if(strcmp(typepatch,'D'))
        A = PCTPMATRIX(PC,zeros(getnbddlfree(Spart),getnbddlfree(Spart)));
        L=PCTPMATRIX(PC,1);
        %
        psi=getcharin('psi',varargin)
        VS = getfuns(psi);
        
        for i=1:getm(psi)
            
            a1 = MULTILINFORM([1,1,0,0]);
            a2 = MULTILINFORM([0,1,0,1]);
            a3 = MULTILINFORM([1,0,1,0]);
            a4 = MULTILINFORM([0,0,1,1]);
            
            %Li = getphi(VS{i},i);
            %Li = PCMATRIX(Li,[size(Li,2),1],getpcgroup(PCphi,i));
            Li = getphi(VS{i});
            Li = setphi(L,Li)
            for j=1:getm(psi)
                
                L2=PCTPMATRIX(PCphi,1);
                Lj = getphi(VS{j});
                Lj = setphi(L2,Lj);
                lilj = Li*Lj;
                
                
                
                A = A + a1{S}(:,:,getphi0(VS{i}),getphi0(VS{j}))*lilj+...
                    a2{S}(:,:,getphi0(VS{i}),getphi0(VS{j}))*lilj+...
                    a3{S}(:,:,getphi0(VS{i}),getphi0(VS{j}))*lilj+...
                    a4{S}(:,:,getphi0(VS{i}),getphi0(VS{j}))*lilj;
            end
        end
        
        Apatchsep = calc_ximasse(A);
    end
    
    
end
end
