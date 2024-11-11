function pop=tent_chaos_init(N,dim,lb,ub)
%n_pop 种群规模
%n_var 维数
%lb 下界
%ub 上界

% Tent混沌映射序列
tent=1.1; %tent混沌系数
Tent_mapping=rand(N,dim);%随机序列
for i=1:N
    for j=2:dim
        if Tent_mapping(i,j-1)<tent
            Tent_mapping(i,j)=Tent_mapping(i,j-1)/tent;
        elseif Tent_mapping(i,j-1)>=tent
            Tent_mapping(i,j)=(1-Tent_mapping(i,j-1))/(1-tent);
        end
    end
end

factor=Tent_mapping;

% figure
% plot(z(:,1),"black");
% hold on;

%% 初始化种群

% pop = lb + z*(ub - lb);  % 初始化种群
pop = lb + factor.*(ub - lb); %初始化种群

figure
scatter(pop(:,1), pop(:,2), 'red')
title("Tent映射初始化种群")
xlabel("x")
ylabel("y")
box on;
