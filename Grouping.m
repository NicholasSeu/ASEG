function [H_UL_re] = Grouping(SysSet,H_UL,Theta)
% Grouping 
switch SysSet.GroupMethod
    case 'RAND'
        group=RAND_Grouping(SysSet);
    case 'ASEG'
        group=ASEG_Grouping(SysSet,Theta);
    case 'GREED'
        group=GREED_Grouping(SysSet,H_UL);
    case 'SUS'
        group=SUS_Grouping(SysSet,H_UL);
    case 'SUS-G1'
        group=SUS_Grouping(SysSet,H_UL);
    case 'SEGA'
        group=SEGA_Grouping(SysSet,H_UL,RAND_Grouping(SysSet));
    otherwise
        warning('Unexpected GroupMethod type.')
end
H_UL_re={};
for gro=1:SysSet.GroupNum
    group_num=sum(group(gro,:)>0);%Each group of nonzero elements
    H_UL_re{gro}=H_UL(:,group(gro,1:group_num));%Channels are rearranged by group
end
end
 
function [group]=RAND_Grouping(SysSet)
reshape_r=SysSet.GroupNum;%The number of rows represents the number of groups
reshape_c=fix(SysSet.K/reshape_r)+(SysSet.K/reshape_r-fix(SysSet.K/reshape_r)>0);%The number of columns represents group members
group=reshape([randperm(SysSet.K),zeros(1,reshape_r*reshape_c-SysSet.K)],reshape_r,reshape_c);%Round it up and group it
end
 
function [group]=ASEG_Grouping(SysSet,Theta)
[~,order]=sort(Theta);
reshape_r=SysSet.GroupNum;
reshape_c=fix(SysSet.K/reshape_r)+(SysSet.K/reshape_r-fix(SysSet.K/reshape_r)>0);
group=reshape([order.',zeros(1,reshape_r*reshape_c-SysSet.K)],reshape_r,reshape_c);
end
 
function [group]=SEGA_Grouping(SysSet,H_UL,pre_group)
if (SysSet.GroupNum==1)
    group=pre_group;
    return;
end
Elite=30;
reshape_r=SysSet.GroupNum;
reshape_c=fix(SysSet.K/reshape_r)+(SysSet.K/reshape_r-fix(SysSet.K/reshape_r)>0);
Group=zeros(reshape_r,reshape_c,Elite*3);
for i=1:Elite %10 Elites
    Group(:,:,i)=pre_group;
end
rate_pop=zeros(1,Elite*3);
rate_rec=zeros(SysSet.ITER);
for i=1:SysSet.ITER
    for e=1:Elite
        mut=randi([1 Elite],1,1);
        for mu=1:mut%mutation
            tar1r=randi([1 reshape_r-1],1,1);
            tar1c=randi([1 reshape_c-1],1,1);
            tar2r=randi([tar1r reshape_r],1,1);
            tar2c=randi([tar1c reshape_c],1,1);
            Group(:,:,e+Elite)=Group(:,:,e);
            Group(tar1r,tar1c,e+Elite)=Group(tar2r,tar2c,e);
            Group(tar2r,tar2c,e+Elite)=Group(tar1r,tar1c,e);
        end
        mut=randi([1 Elite],1,1);
        for mu=1:mut%mutation
            tar1r=randi([1 reshape_r-1],1,1);
            tar1c=randi([1 reshape_c-1],1,1);
            tar2r=randi([tar1r reshape_r],1,1);
            tar2c=randi([tar1c reshape_c],1,1);
            Group(:,:,e+2*Elite)=Group(:,:,e);
            Group(tar1r,tar1c,e+2*Elite)=Group(tar2r,tar2c,e);
            Group(tar2r,tar2c,e+2*Elite)=Group(tar1r,tar1c,e);
        end
    end
 
    for pop=1:Elite*3
        H_UL_re={};
        rate_group=zeros(1,SysSet.GroupNum);
        for gro=1:SysSet.GroupNum
            group_num=sum(Group(gro,:,pop)>0);
            H_UL_re{gro}=H_UL(:,Group(gro,1:group_num,pop));
        end
        for gro=1:SysSet.GroupNum
            P=PreCode_Gen(SysSet,H_UL_re{gro}.');
            rate_group(gro)=Rate_Calcu(SysSet,H_UL_re{gro}.',P);
            rate_group(gro)=rate_group(gro)*size(H_UL_re{gro},2)/SysSet.K;
        end
        rate_pop(pop)=sum(rate_group);
    end
    [~,order]=sort(rate_pop,'descend');
    Group(:,:,1:Elite)=Group(:,:,order(1:Elite));
    rate_rec(i)=rate_pop(order(1));
    if(i>50)
        if(max(rate_rec(i-20:i))<1.001*max(rate_rec(1:(i-30))))
            break;
        end
    end
end
group=Group(:,:,1);
end
 
function [group]=GREED_Grouping(SysSet,H_UL)
%ZHAO2017TWC
rest_users=1:SysSet.K;
reshape_r=SysSet.GroupNum;
reshape_c=fix(SysSet.K/reshape_r)+(SysSet.K/reshape_r-fix(SysSet.K/reshape_r)>0);
Group=zeros(reshape_r,reshape_c);
for gro=1:reshape_r
    while((sum(Group(gro,:)>0)<reshape_c) && sum(rest_users>0))
        H_UL_re={};
        I_Max=0;
        R_Max=0;
        next=sum(Group(gro,:)>0)+1;
        for i=1:size(rest_users,2)
            Group(gro,next)=rest_users(i);
            H_UL_re{gro}=H_UL(:,Group(gro,1:next));
            P=PreCode_Gen(SysSet,H_UL_re{gro}.');
            rate=Rate_Calcu(SysSet,H_UL_re{gro}.',P);
            if(rate>=R_Max)
                R_Max=rate;                I_Max=i;
            end
            Group(gro,next)=0;
        end
        Group(gro,next)=rest_users(I_Max);
        rest_users(I_Max)=[];
    end
end
group=Group;
end
 
function [group]=SUS_Grouping(SysSet,H_UL)
% YOO2006JSAC & LEE2018TCOM
H=H_UL.';
cal_T=1:SysSet.K;% Eq-17
reshape_r=SysSet.GroupNum;
reshape_c=fix(SysSet.K/reshape_r)+(SysSet.K/reshape_r-fix(SysSet.K/reshape_r)>0);
G_So=zeros(reshape_r,reshape_c);% Eq-19
for gro=1:reshape_r
    gb=zeros(size(H));
    while((sum(G_So(gro,:)>0)<reshape_c) && sum(cal_T>0))
        k_max=0;%Eq-22
        g_k_Max=0;
        G_i=sum(G_So(gro,:)>0)+1;%Eq-18
        g=zeros(size(H));
        norm_g=zeros(1,size(cal_T,2));
        for k=1:size(cal_T,2)
            eq20=0;
            for j=1:(G_i-1)
                eq20=H(cal_T(k),:)*gb(j,:)'*gb(j,:)/((norm(gb(j,:),'fro'))^2);
            end
            g(k,:)=H(cal_T(k),:)-eq20;
            norm_g(k)=norm(g(k,:));
            if(norm(g(k,:))>g_k_Max)
                k_max=k;
                g_k_Max=norm(g(k,:));
            end
        end
        gb(G_i,:)=g(k_max,:);
        G_So(gro,G_i)=cal_T(k_max);
        cal_T(k_max)=[];
    end
end
group=G_So;
end
 
function [Simu_Rate] = Rate_Calcu(SysSet,H_DL,P)
reci=SysSet.Power*(abs(H_DL*P).^2);
signal=SysSet.Power*diag(abs(H_DL*P).^2);
intef=reci-diag(signal);
SINR=(signal)./(SysSet.n_var+sum(intef,2));
rate=log2(1+SINR);
Simu_Rate=sum(rate);
end

