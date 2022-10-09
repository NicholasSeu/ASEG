function [H_UL_re] = Grouping(SysSet,H_UL,Theta)
% Grouping 用户分组
% EPGA_Grouping与GREED_Grouping使用了Simulation中的Rate_Calcu
%   此处显示详细说明
switch SysSet.GroupMethod
    case 'RAND'%随机分组，直接reshape
        group=RAND_Grouping(SysSet);
    case 'ASEG'%快速分组，按角度排序，然后按顺序分给各组
        group=ASEG_Grouping(SysSet,Theta);
    case 'GREED'%使用了Simulation中的Rate_Calcu
        group=GREED_Grouping(SysSet,H_UL);
    case 'SUS'
        group=SUS_Grouping(SysSet,H_UL);
    case 'SUS-G1'
        group=SUS_Grouping(SysSet,H_UL);
    case 'SEGA'%使用了Simulation中的Rate_Calcu
        group=SEGA_Grouping(SysSet,H_UL,RAND_Grouping(SysSet));
    case 'ASEG-GA'%使用了Simulation中的Rate_Calcu
        group=SEGA_Grouping(SysSet,H_UL,ASEG_Grouping(SysSet,Theta));
    case 'CVX'%新建文件夹
        group=CVX_Grouping(SysSet,H_UL);
    otherwise
        warning('Unexpected GroupMethod type.')
end
%H_UL_re=zeros(SysSet.M,reshape_c,reshape_r);
H_UL_re={};
for gro=1:SysSet.GroupNum
    group_num=sum(group(gro,:)>0);%每组非零元素
    H_UL_re{gro}=H_UL(:,group(gro,1:group_num));%信道按组重排
end
end

function [group]=RAND_Grouping(SysSet)
reshape_r=SysSet.GroupNum;%行数代表组数
reshape_c=fix(SysSet.K/reshape_r)+(SysSet.K/reshape_r-fix(SysSet.K/reshape_r)>0);%列数代表组员
group=reshape([randperm(SysSet.K),zeros(1,reshape_r*reshape_c-SysSet.K)],reshape_r,reshape_c);%fix取整，分组
end

function [group]=ASEG_Grouping(SysSet,Theta)
[~,order]=sort(Theta);
reshape_r=SysSet.GroupNum;%行数代表组数
reshape_c=fix(SysSet.K/reshape_r)+(SysSet.K/reshape_r-fix(SysSet.K/reshape_r)>0);%列数代表组员
group=reshape([order.',zeros(1,reshape_r*reshape_c-SysSet.K)],reshape_r,reshape_c);%fix取整，分组
end

function [group]=SEGA_Grouping(SysSet,H_UL,pre_group)
% 参考鲍嘉龙毕设精英保留遗传算法
if (SysSet.GroupNum==1)
    group=pre_group;
    return;
end
Elite=30;
reshape_r=SysSet.GroupNum;%行数代表组数
reshape_c=fix(SysSet.K/reshape_r)+(SysSet.K/reshape_r-fix(SysSet.K/reshape_r)>0);%列数代表组员
Group=zeros(reshape_r,reshape_c,Elite*3);%组，城员，精英
for i=1:Elite %10个精英
    %Group(:,:,i)=reshape([order.',zeros(1,reshape_r*reshape_c-SysSet.K)],reshape_r,reshape_c);%fix取整，分组
    Group(:,:,i)=pre_group;
end
rate_pop=zeros(1,Elite*3);
rate_rec=zeros(SysSet.ITER);
for i=1:SysSet.ITER
    for e=1:Elite
        mut=randi([1 Elite],1,1);%随机变异次数
        for mu=1:mut%mutation
            tar1r=randi([1 reshape_r-1],1,1);
            tar1c=randi([1 reshape_c-1],1,1);
            tar2r=randi([tar1r reshape_r],1,1);
            tar2c=randi([tar1c reshape_c],1,1);
            Group(:,:,e+Elite)=Group(:,:,e);
            Group(tar1r,tar1c,e+Elite)=Group(tar2r,tar2c,e);
            Group(tar2r,tar2c,e+Elite)=Group(tar1r,tar1c,e);
        end
        mut=randi([1 Elite],1,1);%随机变异次数
        for mu=1:mut%mutation第二轮变异
            tar1r=randi([1 reshape_r-1],1,1);
            tar1c=randi([1 reshape_c-1],1,1);
            tar2r=randi([tar1r reshape_r],1,1);
            tar2c=randi([tar1c reshape_c],1,1);
            Group(:,:,e+2*Elite)=Group(:,:,e);
            Group(tar1r,tar1c,e+2*Elite)=Group(tar2r,tar2c,e);
            Group(tar2r,tar2c,e+2*Elite)=Group(tar1r,tar1c,e);
        end
    end
    %速率计算
    for pop=1:Elite*3
        H_UL_re={};
        rate_group=zeros(1,SysSet.GroupNum);
        for gro=1:SysSet.GroupNum
            group_num=sum(Group(gro,:,pop)>0);%每组非零元素
            H_UL_re{gro}=H_UL(:,Group(gro,1:group_num,pop));%信道按组重排
        end
        for gro=1:SysSet.GroupNum
            P=PreCode_Gen(SysSet,H_UL_re{gro}.');
            rate_group(gro)=Rate_Calcu(SysSet,H_UL_re{gro}.',P);
            rate_group(gro)=rate_group(gro)*size(H_UL_re{gro},2)/SysSet.K;%资源分配
        end
        rate_pop(pop)=sum(rate_group);
    end
    %排序择优
    [~,order]=sort(rate_pop,'descend');
    Group(:,:,1:Elite)=Group(:,:,order(1:Elite));
    rate_rec(i)=rate_pop(order(1));
    if(i>50)
        if(max(rate_rec(i-20:i))<1.001*max(rate_rec(1:(i-30))))
            break;%如果连续多次优化幅度过小
        end
    end
end
group=Group(:,:,1);
%plot(rate_rec);
%给出结果
end

function [group]=GREED_Grouping(SysSet,H_UL)
%ZHAO2017TWC
rest_users=1:SysSet.K;
reshape_r=SysSet.GroupNum;%行数代表组数
reshape_c=fix(SysSet.K/reshape_r)+(SysSet.K/reshape_r-fix(SysSet.K/reshape_r)>0);%列数代表组员
Group=zeros(reshape_r,reshape_c);
for gro=1:reshape_r%对于每一组
    while((sum(Group(gro,:)>0)<reshape_c) && sum(rest_users>0))%组人数不满且有剩余用户时
        H_UL_re={};
        I_Max=0;
        R_Max=0;
        next=sum(Group(gro,:)>0)+1;%下一位编号
        for i=1:size(rest_users,2)%遍历所有剩余用户
            Group(gro,next)=rest_users(i);
            H_UL_re{gro}=H_UL(:,Group(gro,1:next));
            P=PreCode_Gen(SysSet,H_UL_re{gro}.');
            rate=Rate_Calcu(SysSet,H_UL_re{gro}.',P);
            if(rate>=R_Max)
                R_Max=rate;%rate恒正，故第一位用户一定可以加入备选
                I_Max=i;
            end
            Group(gro,next)=0;
        end%万一加入用户反而使得和速率下降呢？那也得加，你还能不让人上网了？
        Group(gro,next)=rest_users(I_Max);%用户入组
        %选择拥有最大的{加入组后，组中和速率}的用户
        rest_users(I_Max)=[];%从剩余用户中删除
    end
end
group=Group;
end

function [group]=SUS_Grouping(SysSet,H_UL)
% YOO2006JSAC & LEE2018TCOM
H=H_UL.';
cal_T=1:SysSet.K;% Eq-17
reshape_r=SysSet.GroupNum;%行数代表组数
reshape_c=fix(SysSet.K/reshape_r)+(SysSet.K/reshape_r-fix(SysSet.K/reshape_r)>0);%列数代表组员
G_So=zeros(reshape_r,reshape_c);% Eq-19
for gro=1:reshape_r%对于每一组
    gb=zeros(size(H));%bracket永久工
    while((sum(G_So(gro,:)>0)<reshape_c) && sum(cal_T>0))%组人数不满且有剩余用户时
        k_max=0;%Eq-22
        g_k_Max=0;
        G_i=sum(G_So(gro,:)>0)+1;%Eq-18 下一位编号,空组next为1
        g=zeros(size(H));%临时工
        norm_g=zeros(1,size(cal_T,2));
        for k=1:size(cal_T,2)%遍历所有剩余用户
            eq20=0;
            for j=1:(G_i-1)
                eq20=H(cal_T(k),:)*gb(j,:)'*gb(j,:)/((norm(gb(j,:),'fro'))^2);
            end
            g(k,:)=H(cal_T(k),:)-eq20;
            norm_g(k)=norm(g(k,:));
            if(norm(g(k,:))>g_k_Max)%如果范数更大，则备选入组
                k_max=k;
                g_k_Max=norm(g(k,:));
            end
        end%万一加入用户反而使得和速率下降呢？那也得加，你还能不让人上网了？
        gb(G_i,:)=g(k_max,:);
        G_So(gro,G_i)=cal_T(k_max);%用户入组
        %选择拥有最大的{加入组后，组中和速率}的用户
        cal_T(k_max)=[];%从剩余用户中删除
    end
end
group=G_So;
end

function [Simu_Rate] = Rate_Calcu(SysSet,H_DL,P)
%Rate_Calcu 计算整组和速率
%   此处显示详细说明
reci=SysSet.Power*(abs(H_DL*P).^2);%接收信号
signal=SysSet.Power*diag(abs(H_DL*P).^2);%信号部分
intef=reci-diag(signal);%干扰部分
SINR=(signal)./(SysSet.n_var+sum(intef,2));%每个用户接收到的SINR
rate=log2(1+SINR);%(User_Posi(:,1).^2+User_Posi(:,2).^2).^0.5
Simu_Rate=sum(rate);
end
