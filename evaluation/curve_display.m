clear;clc;close all
addpath('Optimization Algorthims');
addpath('My Optimization Algorthims');
addpath('Apperance');
%% initial
pop_num=30; % Number of search agents 种群数量
Max_iter=300; % Maximum numbef of iterations 最大迭代次数
dim = 30; % 可选 2, 10, 30, 50, 100
func_count=30;
% 1-3 Unimodal Functions 50
% 4-10 Simple Multimodal Functions 100
% 11-20 Hybrid Functions 200
% 21-30 Composition Functions

%% 曲线外观列表
% 加载保存的数据
load('cec2017/Apperance/style_info.mat');
%随机取用一个标志
% line_markers=markers{randi(length(markers))};
% marker_list=markers(randperm(length(markers)));
marker_list=markers;

for indexi=1:func_count
    %% choose the function by its name
    Function_num=indexi; % 使用方程的名字，对应 Functions_details 文件
    Function_name=strcat('F',num2str(Function_num));
    [lb,ub,dim,fobj]=Get_Functions_cec2017(Function_num,dim);  %得到具体的方程即目标函数，变量的维度，变量的上下限
    lb=lb.*ones(1,dim);
    ub=ub.*ones(1,dim);
    %% Initialization of comparison variables for each algorithm
    Time_compare=[];      %算法的运行时间比较
    Fival_compare=[];       %算法的最终目标比较
    curve_compare=[];     %算法的收敛曲线比较
    bestX_compare=[];   %算法
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

    % disp('Loaded colors');
    % disp(colors);

    %% 结果绘图
    figure('Position',[500 500 660 290])
    % 绘制搜索空间
    subplot(1,2,1);
    % func_plot(Function_name);
    % 画出目标函数的图示,只能画到三维，目标函数dim的设置可能很多维
    % x=lb/2:0.1:ub/2;  %x轴,y轴范围
    % y=x;
    % L_num=length(x);
    % f_value=[];      %对应x,y的函数值
    % for i=1:L_num
    %     for j=1:L_num
    %         f_value(i,j)=fobj([x(i),y(j)]);
    %     end
    % end
    % surfc(x,y,f_value,'LineStyle','none');
    % title('Test function')
    % xlabel('x_1');
    % ylabel('x_2');
    % zlabel([Function_name,'( x_1 , x_2 )'])
    % grid off

    % Convergence curve
    subplot(1,2,2);
    %% 画迭代过程曲线
    step=Max_iter/15;%指定步长
    for N=1:length(Fival_compare)
        plot(1:step:Max_iter,log(curve_compare(N,1:step:Max_iter)),...
            'LineStyle', line_styles{mod(N,2)+1}, ...
            'Marker',marker_list{N}, ...
            'LineWidth', line_widths(N), ...
            'Color',color_all(N,:));
        hold on
    end
    % plot(curve_compare(length(Fival_compare)-1,:),'LineWidth',2)
    % hold on
    % 设置 x 轴刻度为原始范围
    title(Function_name);
    xlabel('Iteration');
    ylabel('log(Fitness)');
    % x_tick(0:5:15)=num2cell(0:step*5:Max_iter);
    y = cell(1,length(0:step:Max_iter));
    y(1:5:end)=num2cell(0:step*5:Max_iter);
    h=gca;  %获取句柄
    h.XTickLabelRotation=0; % y轴刻度标签要旋转改成YTickLabelRotation
    % XTickLabelRotation 我理解为刻度按标签的旋转角度属性
    % -90 是旋转度数，正负就是顺时针或逆时针旋转
    set(h,'xtick',0:step:Max_iter,'xticklabels',y);
    grid on
    box on
    legend(name_all)
    %以下可以显示
    % display(['The best solution obtained by SSA is : ', num2str(bestX)]);
    % display(['The best optimal value of the objective funciton found by SSA is : ', num2str(fMin)]);
    %%  智能优化算法和内点法比较
    % t1=clock;
    % lb= lb.*ones( 1,dim );    %     约束上限
    % ub= ub.*ones( 1,dim );    %   约束下限   
    % x0 = zeros( 1, dim );   %随机初始化n个种群
    % [xsol,fval] = fmincon(fobj,x0,[],[],[],[],lb,ub);
    % time_interpoint=(clock-t1);         %内点法
    % Fival_compare=[Fival_compare,fval];
    % Time_compare=[Time_compare,time_interpoint(end)];
    % name_all{1,iter}='inter-point';
    % iter=iter+1;
    %%  运行值和最终目标函数比较
    % figure(3)
    % color=color_list(randperm(length(color_list)),:);
    % width=0.4; %柱状图宽度
    % for  i=1:length(Fival_compare) 
    %    %set(bar(i,Fival_compare(i),width),'FaceColor',color(i,:),'EdgeColor',[0,0,0],'LineWidth',1)
    %    set(bar(i,Fival_compare(i),width),'EdgeColor',[0,0,0],'LineWidth',1)
    %    hold on
    % 
    %    %在柱状图 x,y 基础上 绘制误差 ,low为下误差，high为上误差，LineStyle 误差图样式，'Color' 误差图颜色  
    %    % 'LineWidth', 线宽,'CapSize',误差标注大小
    %    % errorbar(i, y(i), low(i), high(i), 'LineStyle', 'none', 'Color', color(i+3,:), 'LineWidth', 1.5,'CapSize',18);
    % end
    % % set(bar(length(Fival_compare)-1,Fival_compare(length(Fival_compare)-1),width),'EdgeColor',[0,0,0],'LineWidth',2)
    % % hold on
    % ylabel('obj-value')
    % ax=gca;
    % ax.XTick = 1:1:length(Fival_compare);
    % set(gca,'XTickLabel',name_all,"LineWidth",1);
    % set(gca,"FontName","Times New Roman","FontSize",12,"LineWidth",1)
    %% 运行时间比较
    % figure(4)
    % color=color_list(randperm(length(color_list)),:);
    % width=0.7; %柱状图宽度
    % for  i=1:length(Fival_compare)-1 
    %    % set(bar(i,Time_compare(i),width),'FaceColor',color(i,:),'EdgeColor',[0,0,0],'LineWidth',1)
    %    set(bar(i,Time_compare(i),width),'EdgeColor',[0,0,0],'LineWidth',1)
    %    hold on
    %    %在柱状图 x,y 基础上 绘制误差 ,low为下误差，high为上误差，LineStyle 误差图样式，'Color' 误差图颜色  
    %    % 'LineWidth', 线宽,'CapSize',误差标注大小
    % %    errorbar(i, y(i), low(i), high(i), 'LineStyle', 'none', 'Color', color(i+3,:), 'LineWidth', 1.5,'CapSize',18);
    % end
    % endofAO=length(Fival_compare)-1;
    % set(bar(i,Time_compare(endofAO),width),'EdgeColor',[0,0,0],'LineWidth',2)
    % hold on
    % % errorbar(i, y(endofAO), low(endofAO), high(endofAO), 'LineStyle', 'none', 'Color', color(endofAO+3,:), 'LineWidth', 2,'CapSize',18);
    % % hold
    % ylabel('Time')
    % ax=gca;
    % ax.XTick = 1:1:length(Fival_compare);
    % set(gca,'XTickLabel',name_all,"LineWidth",1);
    % set(gca,"FontName","Times New Roman","FontSize",12,"LineWidth",1)
    %%


end

