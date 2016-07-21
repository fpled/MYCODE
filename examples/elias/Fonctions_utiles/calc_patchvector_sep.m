function bpatchsep = calc_patchvector_sep(Spart,PC,typepatch,varargin)

%function bpatch = calc_patchvector_sep(Spart,PC,typepatch,varargin)
% ('beh' pour diffusion, 'D' pour Dirichlet et 'N' pour Neumann)
%calc_patchvector_sep(Spart,PC,typepatch,num,'beh','k',ktp)
if(nargin==2)
    error('Veillez entre les type de patch : ''behavior'' pour diffusion \n''D'' pour Dirichlet et ''N'' pour Neumann  ');
end

if(nargin>=3)
    if strcmp(typepatch,'behavior')
        lk = LINFORM(0);
        
        for k=2:size(Spart,2)
            bpatch{k} = lk{Spart{k}}(:)*one(PC);
            bpatchsep{k} = SEPVECTOR(bpatch{k});
        end
        
    end
end
end
