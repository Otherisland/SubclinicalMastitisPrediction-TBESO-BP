function [bestf, bestx, curve] = SO(varargin)
%varargin包括popsize,maxgen,lb,ub,dim
% global optimizer   % 获得提供的优化算法
popsize = cell2mat(varargin(1));   % 种群大小
maxgen = cell2mat(varargin(2));   % 最大迭代次数
lb = cell2mat(varargin(3));   % 变量下界,low
ub = cell2mat(varargin(4));   % 变量上界,upper
dim = cell2mat(varargin(5));  % 变量个数,dimension
fobj = cell2mat(varargin(6));  % 适应度函数,fobj


%initial 
vec_flag=[1,-1];
Threshold=0.25;%食物质量阈值
Thresold2= 0.6;%温度阈值
C1=0.5;%可优化
C2=.05;%可优化
C3=2;%可优化

%% 初始化种群

pos=lb+rand(popsize,dim).*(ub-lb);
X=pos;

%计算种群中每条蛇的适应度
for i=1:popsize
 fitness(i)=feval(fobj,X(i,:));
 %feval('函数名',输入值)
end
[GYbest, gbest] = min(fitness);
bestx = X(gbest,:);

% lb=lb*ones(1,dim);
% ub=ub*ones(1,dim);

%% Diving the snakes into two equal groups males and females
Nm=round(popsize/2);%eq.(2&3)
Nf=popsize-Nm;
%种群赋值
Xm=X(1:Nm,:);
Xf=X(Nm+1:popsize,:);
%适应度赋值
fitness_m=fitness(1:Nm);
fitness_f=fitness(Nm+1:popsize);
%最佳雌雄
[fitnessBest_m, gbest1] = min(fitness_m);
Xbest_m = Xm(gbest1,:);
[fitnessBest_f, gbest2] = min(fitness_f);
Xbest_f = Xf(gbest2,:);
%% 
for t = 1:maxgen
    % 温度系数
    Temp=exp(-((t)/maxgen));  %eq.(4)
  %食物质量系数
  Q=C1*exp(((t-maxgen)/(maxgen)));%eq.(5)
    if Q>1        Q=1;    end
    % Exploration Phase (no Food)
if Q<Threshold %食物质量低，exploration 全局寻优
    for i=1:Nm
        for j=1:1:dim
                rand_leader_index = floor(Nm*rand()+1);
                X_randm = Xm(rand_leader_index, :);
                flag_index = floor(2*rand()+1);
                Flag=vec_flag(flag_index);
                Am=exp(-fitness_m(rand_leader_index)/(fitness_m(i)+eps));%eq.(7)
                Xnewm(i,j)=X_randm(j)+Flag*C2*Am*(ub(j)-lb(j)*rand+lb(j));%eq.(6)
        end
    end
    for i=1:Nf
        for j=1:1:dim
                rand_leader_index = floor(Nf*rand()+1);
                X_randf = Xf(rand_leader_index, :);
                flag_index = floor(2*rand()+1);
                Flag=vec_flag(flag_index);
                Af=exp(-fitness_f(rand_leader_index)/(fitness_f(i)+eps));%eq.(9)
                Xnewf(i,j)=X_randf(j)+Flag*C2*Af*(ub(j)-lb(j))*rand+lb(j);%eq.(8)
        end
    end
else %Exploitation Phase (Food Exists) 食物充足时，局部寻优、斗争或繁衍
    if Temp>Thresold2  %hot 热 exploitation 局部寻优 
        for i=1:Nm
            flag_index = floor(2*rand()+1);
            Flag=vec_flag(flag_index);
            for j=1:1:dim
                Xnewm(i,j)=bestx(j)+C3*Flag*Temp*rand*(bestx(j)-Xm(i,j));%eq.(10)
            end
        end
        for i=1:Nf
            flag_index = floor(2*rand()+1);
            Flag=vec_flag(flag_index);
            for j=1:1:dim
                Xnewf(i,j)=bestx(j)+Flag*C3*Temp*rand*(bestx(j)-Xf(i,j));%eq.(10)
            end
        end
    else %cold 寒冷时，繁衍或斗争
        if rand>0.6 %fight 斗争， rand值决定进入哪种形态
            for i=1:Nm
                for j=1:1:dim
                    FM=exp(-(fitnessBest_f)/(fitness_m(i)+eps));%eq.(13)
                    Xnewm(i,j)=Xm(i,j) +C3*FM*rand*(Q*Xbest_f(j)-Xm(i,j));%eq.(11)
                    
                end
            end
            for i=1:Nf
                for j=1:1:dim
                    FF=exp(-(fitnessBest_m)/(fitness_f(i)+eps));%eq.(14)
                    Xnewf(i,j)=Xf(i,j)+C3*FF*rand*(Q*Xbest_m(j)-Xf(i,j));%eq.(12)
                end
            end
        else %mating 繁衍
            for i=1:Nm
                for j=1:1:dim
                    Mm=exp(-fitness_f(i)/(fitness_m(i)+eps));%eq.(17)
                    Xnewm(i,j)=Xm(i,j) +C3*rand*Mm*(Q*Xf(i,j)-Xm(i,j));%eq.(15
                end
            end
            for i=1:Nf
                for j=1:1:dim
                    Mf=exp(-fitness_m(i)/(fitness_f(i)+eps));%eq.(18)
                    Xnewf(i,j)=Xf(i,j) +C3*rand*Mf*(Q*Xm(i,j)-Xf(i,j));%eq.(16)
                end
            end
            %下蛋
            flag_index = floor(2*rand()+1);
            %孵蛋
            egg=vec_flag(flag_index);
            if egg==1;
                %代替最差的
                [GYworst, gworst] = max(fitness_m);
                Xnewm(gworst,:)=lb+rand(1, dim).*(ub-lb);%eq.(19)
                [GYworst, gworst] = max(fitness_f);
                Xnewf(gworst,:)=lb+rand(1, dim).*(ub-lb);%eq.(20)
            end
        end
    end
end

    %% 更新最优雄蛇
    for j=1:Nm
        Xnewm(j,:)=min(max(Xnewm(j,:), lb), ub);
        y = feval(fobj,Xnewm(j,:));
        if y<fitness_m(j)
            fitness_m(j)=y;
            Xm(j,:)= Xnewm(j,:);
        end
    end
    
    [Ybest1,gbest1] = min(fitness_m);
    
    %% 更新最优雌蛇
    for j=1:Nf
         Xnewf(j,:)=min(max(Xnewf(j,:), lb), ub);
        y = feval(fobj,Xnewf(j,:));
        if y<fitness_f(j)
            fitness_f(j)=y;
            Xf(j,:)= Xnewf(j,:);
        end
    end

    [Ybest2,gbest2] = min(fitness_f);
    
    if Ybest1<fitnessBest_m
        Xbest_m = Xm(gbest1,:);
        fitnessBest_m=Ybest1;
    end
    if Ybest2<fitnessBest_f
        Xbest_f = Xf(gbest2,:);
        fitnessBest_f=Ybest2;
        
    end
    if Ybest1<Ybest2
        curve(t)=min(Ybest1);
    else
        curve(t)=min(Ybest2);
    
    end
    %% 返回最优解
    if fitnessBest_m<fitnessBest_f
        GYbest=fitnessBest_m;
        bestx=Xbest_m;%最优Y对应的X,pos
    else
        GYbest=fitnessBest_f;
        bestx=Xbest_f;%最优Y对应的X,pos
    end
    
end
bestf= GYbest; %Ybest1 Ybest2 里最小的那个



