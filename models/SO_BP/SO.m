function [curve, bestx, bestf] = SO(fobj, varargin)

popsize = cell2mat(varargin(1));
maxgen = cell2mat(varargin(2));
lb = cell2mat(varargin(3));
ub = cell2mat(varargin(4));
dim = cell2mat(varargin(5));

%initial 
vec_flag=[1,-1];
Threshold=0.25;
Thresold2= 0.6;
C1=0.5;
C2=.05;
C3=2;

%% initial population
X=lb+rand(popsize,dim).*(ub-lb);
parfor i=1:popsize
 fitness(i)=feval(fobj,X(i,:));
end
[GYbest, gbest] = min(fitness);
bestx = X(gbest,:);

%% Diving the snakes into two equal groups males and females
Nm=round(popsize/2);%eq.(2&3)
Nf=popsize-Nm;

Xm=X(1:Nm,:);
Xf=X(Nm+1:popsize,:);

fitness_m=fitness(1:Nm);
fitness_f=fitness(Nm+1:popsize);

[fitnessBest_m, gbest1] = min(fitness_m);
Xbest_m = Xm(gbest1,:);
[fitnessBest_f, gbest2] = min(fitness_f);
Xbest_f = Xf(gbest2,:);
%% 
for t = 1:maxgen

    Temp=exp(-((t)/maxgen));  %eq.(4)

  Q=C1*exp(((t-maxgen)/(maxgen)));%eq.(5)
    if Q>1        Q=1;    end
    % Exploration Phase (no Food)
if Q<Threshold 
    parfor i=1:Nm
        for j=1:1:dim
            rand_leader_index = floor(Nm*rand()+1);
            X_randm = Xm(rand_leader_index, :);
            flag_index = floor(2*rand()+1);
            Flag=vec_flag(flag_index);
            Am=exp(-fitness_m(rand_leader_index)/(fitness_m(i)+eps));%eq.(7)
            Xnewm(i,j)=X_randm(j)+Flag*C2*Am*((ub(j)-lb(j))*rand+lb(j));%eq.(6)
        end
    end
    parfor i=1:Nf
        for j=1:1:dim
            rand_leader_index = floor(Nf*rand()+1);
            X_randf = Xf(rand_leader_index, :);
            flag_index = floor(2*rand()+1);
            Flag=vec_flag(flag_index);
            Af=exp(-fitness_f(rand_leader_index)/(fitness_f(i)+eps));%eq.(9)
            Xnewf(i,j)=X_randf(j)+Flag*C2*Af*((ub(j)-lb(j))*rand+lb(j));%eq.(8)
        end
    end
else %Exploitation Phase (Food Exists)
    if Temp>Thresold2  %hot
        parfor i=1:Nm
            flag_index = floor(2*rand()+1);
            Flag=vec_flag(flag_index);
            for j=1:1:dim
                Xnewm(i,j)=bestx(j)+C3*Flag*Temp*rand*(bestx(j)-Xm(i,j));%eq.(10)
            end
        end
        parfor i=1:Nf
            flag_index = floor(2*rand()+1);
            Flag=vec_flag(flag_index);
            for j=1:1:dim
                Xnewf(i,j)=bestx(j)+Flag*C3*Temp*rand*(bestx(j)-Xf(i,j));%eq.(10)
            end
        end
    else %cold
        if rand>0.6 %fight
            parfor i=1:Nm
                for j=1:1:dim
                    FM=exp(-(fitnessBest_f)/(fitness_m(i)+eps));%eq.(13)
                    Xnewm(i,j)=Xm(i,j) +C3*FM*rand*(Q*Xbest_f(j)-Xm(i,j));%eq.(11)
                    
                end
            end
            parfor i=1:Nf
                for j=1:1:dim
                    FF=exp(-(fitnessBest_m)/(fitness_f(i)+eps));%eq.(14)
                    Xnewf(i,j)=Xf(i,j)+C3*FF*rand*(Q*Xbest_m(j)-Xf(i,j));%eq.(12)
                end
            end
        else %mating
            parfor i=1:Nm
                for j=1:1:dim
                    Mm=exp(-fitness_f(i)/(fitness_m(i)+eps));%eq.(17)
                    Xnewm(i,j)=Xm(i,j) +C3*rand*Mm*(Q*Xf(i,j)-Xm(i,j));%eq.(15
                end
            end
            parfor i=1:Nf
                for j=1:1:dim
                    Mf=exp(-fitness_m(i)/(fitness_f(i)+eps));%eq.(18)
                    Xnewf(i,j)=Xf(i,j) +C3*rand*Mf*(Q*Xm(i,j)-Xf(i,j));%eq.(16)
                end
            end

            flag_index = floor(2*rand()+1);

            egg=vec_flag(flag_index);
            if egg==1;

                [GYworst, gworst] = max(fitness_m);
                Xnewm(gworst,:)=lb+rand(1, dim).*(ub-lb);%eq.(19)
                [GYworst, gworst] = max(fitness_f);
                Xnewf(gworst,:)=lb+rand(1, dim).*(ub-lb);%eq.(20)
            end
        end
    end
end


    parfor j=1:Nm
        Xnewm(j,:)=min(max(Xnewm(j,:), lb), ub);
        y = feval(fobj,Xnewm(j,:));
        if y<fitness_m(j)
            fitness_m(j)=y;
            Xm(j,:)= Xnewm(j,:);
        end
    end
    
    [Ybest1,gbest1] = min(fitness_m);
    

    parfor j=1:Nf
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

    if fitnessBest_m<fitnessBest_f
        GYbest=fitnessBest_m;
        bestx=Xbest_m;
    else
        GYbest=fitnessBest_f;
        bestx=Xbest_f;
    end
    
end
bestf= GYbest;
    disp('running... ...')
    disp('End Of Run... ...')
end

%_______________________________________________________________



