clear;clc;close all
clear gca
addpath('Optimization Algorthims');
addpath('My Optimization Algorthims');
addpath('Apperance');
%% initial
box_pp = 1;  %可选1，或者其他。当等于1，绘制箱型图，否则不绘制
pop_num=30; % Number of search agents 种群数量
Max_iter=500; % Maximum numbef of iterations 最大迭代次数
dim = 30; % 可选 2, 10, 30, 50, 100
iter_count=30;   %遍历次数
func_count=30;  %算法总数
% 1-3 Unimodal Functions
% 4-10 Simple Multimodal Functions
% 11-20 Hybrid Functions
% 21-30 Composition Functions
RESULT=[];   %统计标准差，平均值，最优值等结果 mean std
rank_wilcoxo_RESULT=[];%wilcoxo统计秩和检验结果
rank_friedman_RESULT=[];%friedman排名结果

%% choose the function by its name
for Function_num=1:func_count
    % Function_num=1; % 使用方程的名字，对应 Functions_details 文件
    Function_name=strcat('F',num2str(Function_num));
    [lb,ub,dim,fobj]=Get_Functions_cec2017(Function_num,dim);  %得到具体的方程即目标函数，变量的维度，变量的上下限
    lb=lb*ones(1,dim);
    ub=ub*ones(1,dim);
    %% Initialization of comparison variables for each algorithm
    Time_compare=cat(3);      %算法的运行时间比较
    Fival_compare=cat(3);       %算法的最终目标比较
    curve_compare=cat(3);     %算法的收敛曲线比较
    ranks = [] ;
    rank_wilcoxo_results_compare = [];   %统计秩和检验结果
    results_compare = cat(3);  %统计标准差，平均值，最优值等结果
    box_plot_compare= cat(3);  %统计箱型图结果
    %% names of these algorithms
    name_all=[];     %算法的名称记录
    for indexi=1:iter_count
        %% Calling algorithm
        
        %% the algorithms of myself
        %% 蛇优化算法 SO 
        t1=clock;
        iter=1;
        [fMin_SO,bestX_SO,SO_curve]=SO(pop_num,Max_iter,lb,ub,dim,fobj);     
        t2=clock;
        time_SO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
        Fival_compare(indexi,iter,:)=fMin_SO;
        Time_compare(indexi,iter,:)=time_SO(end);
        curve_compare(indexi,iter,:)=SO_curve;
        name_all{1,iter}='SO';
        ranks(iter,indexi)=fMin_SO;
        all_ranks(indexi, Function_num) = fMin_algorithm;
        iter=iter+1;
        
        %% TCM_SO 算法
        t1=clock;
        [fMin_TCM_SO,bestX_TCM_SO,TCM_SO_curve]=TCM_SO(pop_num,Max_iter,lb,ub,dim,fobj);%TBESO算法
        t2=clock;
        time_TCM_SO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
        Fival_compare(indexi,iter,:)=fMin_TCM_SO;
        Time_compare(indexi,iter,:)=time_TCM_SO(end);
        curve_compare(indexi,iter,:)=TCM_SO_curve;
        name_all{1,iter}='TCM-SO';
        ranks(iter,indexi)=fMin_TCM_SO;
        iter=iter+1;

        %% BDS_SO 算法
        t1=clock;
        [fMin_BDS_SO,bestX_BDS_SO,BDS_SO_curve]=BDS_SO(pop_num,Max_iter,lb,ub,dim,fobj);%TBESO算法
        t2=clock;
        time_BDS_SO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
        Fival_compare(indexi,iter,:)=fMin_BDS_SO;
        Time_compare(indexi,iter,:)=time_BDS_SO(end);
        curve_compare(indexi,iter,:)=BDS_SO_curve;
        name_all{1,iter}='BDS-SO';
        ranks(iter,indexi)=fMin_BDS_SO;
        iter=iter+1;
             
        %% EOBL_SO 算法
        t1=clock;
        [fMin_EOBL_SO,bestX_EOBL_SO,EOBL_SO_curve]=EOBL_SO(pop_num,Max_iter,lb,ub,dim,fobj);%TBESO算法
        t2=clock;
        time_EOBL_SO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
        Fival_compare(indexi,iter,:)=fMin_EOBL_SO;
        Time_compare(indexi,iter,:)=time_EOBL_SO(end);
        curve_compare(indexi,iter,:)=EOBL_SO_curve;
        name_all{1,iter}='EOBL-SO';
        ranks(iter,indexi)=fMin_EOBL_SO;
        iter=iter+1;

        %% TB_SO 算法
        t1=clock;
        [fMin_TB_SO,bestX_TB_SO,TB_SO_curve]=TB_SO(pop_num,Max_iter,lb,ub,dim,fobj);%TBESO算法
        t2=clock;
        time_TB_SO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
        Fival_compare(indexi,iter,:)=fMin_TB_SO;
        Time_compare(indexi,iter,:)=time_TB_SO(end);
        curve_compare(indexi,iter,:)=TB_SO_curve;
        name_all{1,iter}='TB-SO';
        ranks(iter,indexi)=fMin_TB_SO;
        iter=iter+1;

        %% TE_SO 算法
        t1=clock;
        [fMin_TE_SO,bestX_TE_SO,TE_SO_curve]=TE_SO(pop_num,Max_iter,lb,ub,dim,fobj);%TBESO算法
        t2=clock;
        time_TE_SO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
        Fival_compare(indexi,iter,:)=fMin_TE_SO;
        Time_compare(indexi,iter,:)=time_TE_SO(end);
        curve_compare(indexi,iter,:)=TE_SO_curve;
        name_all{1,iter}='TE-SO';
        ranks(iter,indexi)=fMin_TE_SO;
        iter=iter+1;

        %% BE_SO 算法
        t1=clock;
        [fMin_BE_SO,bestX_BE_SO,BE_SO_curve]=BE_SO(pop_num,Max_iter,lb,ub,dim,fobj);%TBESO算法
        t2=clock;
        time_BE_SO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
        Fival_compare(indexi,iter,:)=fMin_BE_SO;
        Time_compare(indexi,iter,:)=time_BE_SO(end);
        curve_compare(indexi,iter,:)=BE_SO_curve;
        name_all{1,iter}='BE-SO';
        ranks(iter,indexi)=fMin_BE_SO;
        iter=iter+1;

        %% TEBSO 算法
        t1=clock;
        [fMin_TBESO,bestX_TBESO,TBESO_curve]=TBESO(pop_num,Max_iter,lb,ub,dim,fobj);%TBESO算法
        t2=clock;
        time_TBESO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
        Fival_compare(indexi,iter,:)=fMin_TBESO;
        Time_compare(indexi,iter,:)=time_TBESO(end);
        curve_compare(indexi,iter,:)=TBESO_curve;
        name_all{1,iter}='TBESO';
        ranks(iter,indexi)=fMin_TBESO;
        iter=iter+1;

    end
    %% wilcoxo统计结果
    for i=2:length(name_all)%统计秩和检验结果
        rs = ranksum(ranks(1,:),ranks(i,:));
        if isnan(rs)  %当z1与z2完全一致时会出现NaN值，这种概率很小，但是要做一个防止报错
            rs=1;
        end
        rank_wilcoxo_results_compare = [rank_wilcoxo_results_compare,rs]; 
    end 
    rank_wilcoxo_RESULT = [rank_wilcoxo_RESULT;rank_wilcoxo_results_compare];  %统计秩和检验结果
    
    %% Friedman 排名


    %%  Storage the datas 数据存储
    disp(Function_name);
    data=[];
    Sheet_name_ms='CEC2017';
    xlRange=strcat('C',num2str(Function_num*3));
    ylRange=strcat('J',num2str(Function_num*3+1));
    range_ms=strcat(xlRange,':',ylRange);
    disp(range_ms);
    for i=1:length(name_all)
        % display([cell2mat(name_all(i)),': ']);
        % box_plot_compare(Function_num,:,i) = Fival_compare(:,i); %统计箱型图结果
        box_plot_compare=[box_plot_compare;Fival_compare(:,i)'];
        result_msmm=[min(Fival_compare(:,i));std(Fival_compare(:,i));mean(Fival_compare(:,i));median(Fival_compare(:,i));max(Fival_compare(:,i))];
        disp([cell2mat(name_all(i)),': 最优值:',num2str(result_msmm(1),4),' 标准差:',num2str(result_msmm(2),4),' 平均值:',num2str(result_msmm(3),4),' 中值:',num2str(result_msmm(4),4),' 最差值:',num2str(result_msmm(5),4)]);
        results_compare = [results_compare,result_msmm];
        data(:,i)=[mean(Fival_compare(:,i));std(Fival_compare(:,i))];
        % mean min max med std
    end
    xlswrite('output_data\mean_std.xlsx',data,Sheet_name_ms,range_ms);%保存平均值和标准差

    %% 绘制箱型图
    if box_pp == 1 
        subplot(5,6,Function_num)  %4行6列
        mycolor = [0.862745098039216,0.827450980392157,0.117647058823529;...
        0.705882352941177,0.266666666666667,0.423529411764706;...
        0.949019607843137,0.650980392156863,0.121568627450980;...
        0.956862745098039,0.572549019607843,0.474509803921569;...
        0.231372549019608,0.490196078431373,0.717647058823529;...
        0.655541222625894,0.122226545135785,0.325468941131211;...
        0.766665984235466,0.955154623456852,0.755161564478523;...
        0.766665984235466,0.955154623456852,0.755161564478523];  %设置一个颜色库
        %% 开始绘图
        %参数依次为数据矩阵、颜色设置、标记符
        box_figure = boxplot(box_plot_compare','color',[0 0 0],'Symbol','o');
        %设置线宽
        set(box_figure,'Linewidth',1.2);
        boxobj = findobj(gca,'Tag','Box');
        for i = 1:length(name_all)   %因为总共有5个算法，这里记者根据自身实际情况更改！
            patch(get(boxobj(i),'XData'),get(boxobj(i),'YData'),mycolor(i,:),'FaceAlpha',0.5,...
                'LineWidth',0.7);
            % patch(get(boxobj(i),'XData'),get(boxobj(i),'YData'),'FaceAlpha',0.5,...
            %     'LineWidth',0.7);
        end
        set(gca,'XTickLabel',name_all);
        title(['F',num2str(Function_num)])
        hold on
    end 
    if box_pp == 1   %保存箱线图
        saveas(gcf,'箱线图.png')
    end

    %% 保存Wilcoxo秩排名
    Sheet_name_wil_rank='CEC2017';
    range_wil='B3:H32';
    xlswrite('output_data\wilcoxo_rank.xlsx',rank_wilcoxo_RESULT,Sheet_name_wil_rank,range_wil);
    RESULT = [RESULT;results_compare];   %统计标准差，平均值，最优值等结果
    
    
end




