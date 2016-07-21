function [errorModes solution]=test_karloeve(S,RV,V,rs,VTest,rsModel,varargin);
basisType=getparam(S,'basistype');
p=getparam(S,'order');
tb=getparam(S,'typebase');
mi=getparam(S,'modei');
mf=getparam(S,'modef');
regroup=getparam(S,'regroup');
if basisType==1
    if regroup==0
        [X,PC] = PCTPMODEL(RV,'order',p,'pcg','typebase',tb,'nomasse');
    else
        [X,PC]=PCTPMODEL(RV,'order',p,'pcg','typebase',tb,...
            'groups',{[1,2,3,4],[5,6,7,8]},'nomasse');
    end
else
    RV=RANDVARS(RVUNIFORM(0,1),8);
    rs=(rs+ones(size(rs)))/2;
    rsModel=(rsModel+ones(size(rsModel)))/2;
    if basisType==2
        ne=getparam(S,'nmultelem');
        [X PC]=PCTPMODEL(RV,'order',p,'pcg','fedim',1:8,'typebase',tb,'femesh',...
            {ne,ne,ne,ne,ne,ne,ne,ne},'nomasse');
    else
        w=getparam(S,'waveres');
        [X PC]=PCTPMODEL(RV,'order',p,'pcg','waveletdim',1:8,'typebase',tb,...
            'waveletlevel',{w,w,w,w,w,w,w,w},'nomasse');
    end
end
solution = cell(8,1);
for kk=mi:mf
    u=PCTPMATRIXSUM(PC);
    fs=V(:,kk);
    maxRank=10;
    fsm=fs;
    for m=1:maxRank
        [lopt]=get_optimal_l_stat(PC,fsm,rs);
        u=u+lopt;
        us=evaluate_rank_one_fn(PC,lopt,rs);
        if isempty(us)
            fprintf('Max rank is %d\n',m-1);
            break
        else
            fsm=fsm-us;
        end
    end
    solution{kk}=u;
end
errorModes=cell(8,1);
for kk=mi:mf
    error=[];
    modelapprox=zeros(size(VTest,1),1);
    fsModel=VTest(:,kk);
    for m=1:getm(solution{kk})
        modelapprox=modelapprox+evaluate_rank_one_fn(PC,solution{kk}{m},rsModel);
        error(m)=norm((fsModel-modelapprox),2)/norm(fsModel,2);
    end
    errorModes{kk}=error;
end

