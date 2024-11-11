%___________________________________________________________________%
%  Snake Optimizer (SO) source codes version 1.0                    %
%                                                                   %
%  Developed in MATLAB R2021b                                       %
%                                                                   %
%  Author and programmer:  Fatma Hashim & Abdelazim G. Hussien      %
%                                                                   %
%         e-Mail: fatma_hashim@h-eng.helwan.edu.eg                  %
%                 abdelazim.hussien@liu.se                          %
%                 aga08@fayoum.edu.eg                               %
%                                                                   %
%                                                                   %
%   Main paper: Fatma Hashim & Abdelazim G. Hussien                 %
%               Knowledge-based Systems                             %
%               in press,                                           %
%               DOI: 10.1016/j.knosys.2022.108320                   %
%                                                                   %
%___________________________________________________________________%
function [ fval,Xfood,gbest_t] = BEESO(fhd,pos,varargin)
%initial 
N = cell2mat(varargin(1));   % 种群大小
T = cell2mat(varargin(2));   % 最大迭代次数
lb = cell2mat(varargin(3));   % 变量下界,low
ub = cell2mat(varargin(4));   % 变量上界,upper
dim = cell2mat(varargin(5));  % 变量个数,dimension

vec_flag=[1,-1];
Threshold=0.25;
Thresold2= 0.6;
C1=0.5;
C2=0.05;
C3=2;
% X=lb+rand(N,dim)*(ub-lb);
X=pos;
freq=1/dim;
for i=1:N
 fitness(i)=feval(fhd,X(i,:)) ;   
end
[GYbest, gbest] = min(fitness);
Xfood = X(gbest,:);
%Diving the swarm into two equal groups males and females
Nm=round(N/2);%eq.(2&3)
Nf=N-Nm;
Xm=X(1:Nm,:);
Xf=X(Nm+1:N,:);
fitness_m=fitness(1:Nm);
fitness_f=fitness(Nm+1:N);
[fitnessBest_m, gbest1] = min(fitness_m);
Xbest_m = Xm(gbest1,:);
[fitnessBest_f, gbest2] = min(fitness_f);
Xbest_f = Xf(gbest2,:);
for t = 1:T
    Temp=exp(-((t)/T));  %eq.(4)
    Q=C1*exp(((t-T)/(T)));%eq.(5)
    if Q>1        
        Q=1;    
    end
    % Exploration Phase (no Food)
[~, gmworst] = max(fitness_m);
Xworst_m = Xm(gmworst,:);
[~, gfworst] = max(fitness_f);
Xworst_f = Xf(gfworst,:);
if Q<Threshold
     for i=1:Nm
        for j=1:1:dim
            rand_leader_index = floor(Nm*rand()+1);
            X_randm = Xm(rand_leader_index, :);
            flag_index = floor(2*rand()+1);
            Flag=vec_flag(flag_index);
            Am=exp(-fitness_m(rand_leader_index)/(fitness_m(i)+eps));%eq.(7)
            Xnewm1(i,j)=X_randm(j)+Flag*C2*Am*((ub-lb)*rand+lb);%eq.(6)
            Xnewm2(i,j)=Xm(i,j)+rand*(Xbest_m(j)-Xm(i,j))-rand*(Xworst_m(j)-Xm(i,j));%%双向搜索
        end
            fit1=feval(fhd, Xnewm1(i,:)) ;
            fit2=feval(fhd, Xnewm2(i,:)) ;
            if fit1<fit2
                Xnewm(i,:)=Xnewm1(i,:);
            else 
                Xnewm(i,:)=Xnewm2(i,:);
            end
    end
    for i=1:Nf
        for j=1:1:dim
            rand_leader_index = floor(Nf*rand()+1);
            X_randf = Xf(rand_leader_index, :);
            flag_index = floor(2*rand()+1);
            Flag=vec_flag(flag_index);
            Af=exp(-fitness_f(rand_leader_index)/(fitness_f(i)+eps));%eq.(9)
            Xnewf1(i,j)=X_randf(j)+Flag*C2*Af*((ub-lb)*rand+lb);%eq.(8)
            Xnewf2(i,j)=Xf(i,j)+rand*(Xbest_f(j)-Xf(i,j))-rand*(Xworst_f(j)-Xf(i,j));
        end
            fit3=feval(fhd, Xnewf1(i,:)) ;
            fit4=feval(fhd, Xnewf2(i,:)) ;
            if fit3<fit4
                Xnewf(i,:)=Xnewf1(i,:);
            else 
                Xnewf(i,:)=Xnewf2(i,:);
            end
    end
else %Exploitation Phase (Food Exists)
    if Temp>Thresold2  %hot
        for i=1:Nm
            flag_index = floor(2*rand()+1);
            Flag=vec_flag(flag_index);
            for j=1:1:dim
                Xnewm(i,j)=Xfood(j)+C3*Flag*Temp*rand*(Xfood(j)-Xm(i,j));%eq.(10)
            end
        end
        for i=1:Nf
            flag_index = floor(2*rand()+1);
            Flag=vec_flag(flag_index);
            for j=1:1:dim
                Xnewf(i,j)=Xfood(j)+Flag*C3*Temp*rand*(Xfood(j)-Xf(i,j));%eq.(10)
            end
        end
    else %cold
        if rand>0.6 %fight
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
        else%mating
            for i=1:Nm
                for j=1:1:dim
                    Mm=exp(-fitness_f(i)/(fitness_m(i)+eps));%eq.(17)
                    Xnewm(i,j)=Xm(i,j) +C3*rand*Mm*(Q*Xf(i,j)-Xm(i,j));%eq.(15)
                end
            end
            for i=1:Nf
                for j=1:1:dim
                    Mf=exp(-fitness_m(i)/(fitness_f(i)+eps));%eq.(18)
                    Xnewf(i,j)=Xf(i,j) +C3*rand*Mf*(Q*Xm(i,j)-Xf(i,j));%eq.(16)
                end
            end
            flag_index = floor(2*rand()+1);
            egg=vec_flag(flag_index);
            if egg == 1
              X=[Xnewm;Xnewf];
              for i=1:N
                 FIT1(i)=feval(fhd,X(i,:));
              end
                [PAIXU,I1]=sort(FIT1);
                a=ceil(N/2);
                %%变异操作
                F=1/2*(sin(2*pi*freq*t+pi)*(t/T)+1);
%                   F=0.5;
                for i=1:a
                  dx = randperm(N);  
                  q = dx(1);        k = dx(2);        p = dx(3);  
                %要保证与i不同  
                if q == i  
                    q  = dx(4);  
                else
                    if k == i
                    k = dx(4);  
                    else
                        if p == i
                           p = dx(4);  
                        end  
                    end 
               end       
                 %变异的个体来自三个随机父代 
                  for j=1:dim
                    Xnew_son(I1(i),j) = X(p,j) + F*(X(q,j) - X(k,j));     %p,j,k为三个不同的1-Np的随机数   
                  end
                   fit_mson=feval(fhd,Xnew_son(I1(i),:));
                   if fit_mson< PAIXU(i)
                      Xnewm(I1(i),:)= Xnew_son(I1(i),j);
                   end
                end
                for i=a+1:N
                    for j=1:dim
                        z1=I1(i);
                        if rand<0.5
                          Xnew1(z1,j)=Xfood(j)+sign(rand-0.50)*(lb+rand*(ub-lb));
                        else
                          Xnew1(z1,j)=X(z1,j)+sign(rand-0.50)*(lb+rand*(ub-lb));
                        end
                    end
                    FIT3=feval(fhd,Xnew1(z1,:));
                        if FIT3<FIT1(z1)
                            X(z1,:)=Xnew1(z1,:);
                        else
                            X(z1,:)=lb+rand*(ub-lb);
                        end
                end
                for i=1:Nm
                    Xnewm(i,:)=X(i,:);
                    Xnewf(i,:)=X(i+Nm,:);
                end
             end
        end
    end
end
 for i=1:Nm
    ym(i)= feval(fhd,Xnewm(i,:));
    yf(i)= feval(fhd,Xnewf(i,:));
  end
   [~,I1]=sort(ym);
   a=ceil(0.1*Nm);
   for i=1:a
      m_elite(i,:)=Xnewm(I1(i),:);
   end

   for j=1:Nm
      for i=1:dim
         obl_Xnewm(j,i)=rand*( max(m_elite(:,i))+min(m_elite(:,i)))-Xnewm(j,i);
            if  obl_Xnewm(j,i)>ub||obl_Xnewm(j,i)<lb
                obl_Xnewm(j,i)=rand*(max(m_elite(:,i))-min(m_elite(:,i)))+min(m_elite(:,i));
            end
      end
      fit1(j)=feval(fhd,obl_Xnewm(j,:));
      if fit1(j)<ym(j)   
         Xnewm(j,:)=obl_Xnewm(j,:);
      end
      Flag4ub=Xnewm(j,:)>ub;
      Flag4lb=Xnewm(j,:)<lb;
      Xnewm(j,:)=(Xnewm(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
      y=feval(fhd,Xnewm(j,:));
      if y<fitness_m(j)
         fitness_m(j)=y;
         Xm(j,:)= Xnewm(j,:);
      end
  end
 
 [Ybest1,gbest1] = min(fitness_m);
 
  [~,I2]=sort(yf);
  a=ceil(0.1*Nf);
   for i=1:a
      f_elite(i,:)=Xnewf(I2(i),:);
   end
    for j=1:Nf
         for i=1:dim
            obl_Xnewf(j,i)=rand*(max(f_elite(:,i))+min(f_elite(:,i)))-Xnewf(j,i);
            if  obl_Xnewf(j,i)>ub||obl_Xnewf(j,i)<lb
                obl_Xnewf(j,i)=rand*(max(f_elite(:,i))-min(f_elite(:,i)))+min(f_elite(:,i));
            end
         end
         fit2(j)=feval(fhd,obl_Xnewf(j,:));
         if fit2(j)<yf(j)   
            Xnewf(j,:)=obl_Xnewf(j,:);
         end
        Flag4ub=Xnewf(j,:)>ub;
        Flag4lb=Xnewf(j,:)<lb;
        Xnewf(j,:)=(Xnewf(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        y =feval(fhd,Xnewf(j,:));
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
        gbest_t(t)=min(Ybest1);
    else
        gbest_t(t)=min(Ybest2);
        
    end
    if fitnessBest_m<fitnessBest_f
        GYbest=fitnessBest_m;
        Xfood=Xbest_m;
    else
        GYbest=fitnessBest_f;
        Xfood=Xbest_f;
    end
    
end
fval = GYbest;
end





