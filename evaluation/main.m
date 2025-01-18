clear;clc;close all
addpath('Optimization Algorthims');
addpath('My Optimization Algorthims');
addpath('Apperance');
%% initial
pop_num=100; % Number of search agents 种群数量
Max_iter=500; % Maximum numbef of iterations 最大迭代次数
dim = 30; % 可选 2, 10, 30, 50, 100
%% choose the function by its name
Function_num=1; % 使用方程的名字，对应 Functions_details 文件
Function_name=strcat('F',num2str(Function_num));
[lb,ub,dim,fobj]=Get_Functions_cec2017(Function_num,dim);  %得到具体的方程即目标函数，变量的维度，变量的上下限
lb=-100*ones(1,dim);
ub=100*ones(1,dim);
%% Initialization of comparison variables for each algorithm
Time_compare=[];      %算法的运行时间比较
Fival_compare=[];       %算法的最终目标比较
curve_compare=[];     %算法的收敛曲线比较
results_compare = [];  %统计标准差，平均值，最优值等结果
rank_sum_results_compare = [];   %统计秩和检验结果
box_plot_compare= [];  %统计箱型图结果
%% names of these algorithms
name_all=[];     %算法的名称记录
%% Calling algorithm

%% the algorithms of myself
%% 蛇优化算法 SO 
t1=clock;
iter=1;
[fMin_SO,bestX_SO,SO_curve]=SO(pop_num,Max_iter,lb,ub,dim,fobj);     
t2=clock;
time_SO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
Fival_compare=[Fival_compare,fMin_SO];
Time_compare=[Time_compare,time_SO(end)];
curve_compare=[curve_compare;SO_curve];
%mean std min max
name_all{1,iter}='SO';
iter=iter+1;


%% TCM_SO
t1=clock;
[fMin_TCM_SO,bestX_TCM_SO,TCM_SO_curve]=TCM_SO(pop_num,Max_iter,lb,ub,dim,fobj);%TBESO算法
t2=clock;
time_TCM_SO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
Fival_compare=[Fival_compare,fMin_TCM_SO];
Time_compare=[Time_compare,time_TCM_SO(end)];
curve_compare=[curve_compare;TCM_SO_curve];
name_all{1,iter}='TCM-SO';
iter=iter+1;

%% BDS_SO
t1=clock;
[fMin_BDS_SO,bestX_BDS_SO,BDS_SO_curve]=BDS_SO(pop_num,Max_iter,lb,ub,dim,fobj);%TBESO算法
t2=clock;
time_BDS_SO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
Fival_compare=[Fival_compare,fMin_BDS_SO];
Time_compare=[Time_compare,time_BDS_SO(end)];
curve_compare=[curve_compare;BDS_SO_curve];
name_all{1,iter}='BDS-SO';
iter=iter+1;

%% EOBL_SO
t1=clock;
[fMin_EOBL_SO,bestX_EOBL_SO,EOBL_SO_curve]=EOBL_SO(pop_num,Max_iter,lb,ub,dim,fobj);%TBESO算法
t2=clock;
time_EOBL_SO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
Fival_compare=[Fival_compare,fMin_EOBL_SO];
Time_compare=[Time_compare,time_EOBL_SO(end)];
curve_compare=[curve_compare;EOBL_SO_curve];
name_all{1,iter}='EOBL-SO';
iter=iter+1;

%% TB_SO
t1=clock;
[fMin_TB_SO,bestX_TB_SO,TB_SO_curve]=TB_SO(pop_num,Max_iter,lb,ub,dim,fobj);%TBESO算法
t2=clock;
time_TB_SO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
Fival_compare=[Fival_compare,fMin_TB_SO];
Time_compare=[Time_compare,time_TB_SO(end)];
curve_compare=[curve_compare;TB_SO_curve];
name_all{1,iter}='TB-SO';
iter=iter+1;

%% TE_SO
t1=clock;
[fMin_TE_SO,bestX_TE_SO,TE_SO_curve]=TE_SO(pop_num,Max_iter,lb,ub,dim,fobj);%TBESO算法
t2=clock;
time_TE_SO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
Fival_compare=[Fival_compare,fMin_TE_SO];
Time_compare=[Time_compare,time_TE_SO(end)];
curve_compare=[curve_compare;TE_SO_curve];
name_all{1,iter}='TE-SO';
iter=iter+1;

%% BE_SO
t1=clock;
[fMin_BE_SO,bestX_BE_SO,BE_SO_curve]=BE_SO(pop_num,Max_iter,lb,ub,dim,fobj);%TBESO算法
t2=clock;
time_BE_SO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
Fival_compare=[Fival_compare,fMin_BE_SO];
Time_compare=[Time_compare,time_BE_SO(end)];
curve_compare=[curve_compare;BE_SO_curve];
name_all{1,iter}='BE-SO';
iter=iter+1;

%% TEBSO 算法
t1=clock;
[fMin_TBESO,bestX_TBESO,TBESO_curve]=TBESO(pop_num,Max_iter,lb,ub,dim,fobj);%TBESO算法
t2=clock;
time_TBESO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
Fival_compare=[Fival_compare,fMin_TBESO];
Time_compare=[Time_compare,time_TBESO(end)];
curve_compare=[curve_compare;TBESO_curve];
name_all{1,iter}='TBESO';
iter=iter+1;

str(mean(bestX_TBESO),4),' TBESO_curve=',num2str(mean(TBESO_curve),4)]);

% Mean
disp('Function_name Mean:');
for i=1:length(name_all)
    % display([cell2mat(name_all(i)),': ']);
    box_plot_compare = [box_plot_compare;Fival_compare(i)]; %统计箱型图结果
    result_msmm=[min(Fival_compare(i));std(Fival_compare(i));mean(Fival_compare(i));median(Fival_compare(i));max(Fival_compare(i))];
    disp([cell2mat(name_all(i)),': 最优值:',num2str(result_msmm(1)),' 标准差:',num2str(result_msmm(2)),' 平均值:',num2str(result_msmm(3)),' 中值:',num2str(result_msmm(4)),' 最差值:',num2str(result_msmm(5))]);
    results_compare = [results_compare,result_msmm];
end