function pop=tent_chaos_init(n_pop,n_var,lb,ub)
%n_pop 种群规模
%n_var 维数
%lb 下界
%ub 上界

% Tent混沌映射序列
z = rand(n_pop, n_var); % 随机序列
for i=1:n_pop
    for j=1:n_var
        if z(i,j)<0.5
            z(i,j) = 2*z(i,j);
        elseif z(i)>=0.5
            z(i,j) = 2*(1-z(i,j));
        end
    end
end
% 
% figure
% plot(z(:,1),"black");
% hold on;

%% 初始化种群

% pop = lb + z*(ub - lb);  % 初始化种群
pop = lb + z.*(ub - lb); %初始化种群

% figure
% scatter(pop(:,1), pop(:,2), 'red')
% title("Tent映射初始化种群")
% xlabel("x")
% ylabel("y")
% box on;
