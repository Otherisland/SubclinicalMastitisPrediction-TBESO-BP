clear; clc; close all
clear gca
% addpath('Optimization Algorthims');
addpath('My Optimization Algorthims');
addpath('Apperance');

%% Initial settings
box_pp = 0;  % Set to 1 to draw box plots, otherwise set to 0
pop_num = 30; % Number of search agents
Max_iter = 200; % Maximum number of iterations
dim = 30; % Set dimension (2, 10, 30, 50, 100)
iter_count = 30;   % Number of iterations
func_count = 30;  % Total number of functions
algorithms = {'SO', 'TCM_SO', 'BDS_SO', 'EOBL_SO', 'TB_SO', 'TE_SO', 'BE_SO', 'TBESO'};
% Result matrices
RESULT = []; % Matrix to store mean, std, min, median, max
rank_wilcoxo_RESULT = []; % Matrix to store Wilcoxon rank sum test results
rank_friedman_RESULT = []; % Matrix to store Friedman rank results

all_ranks = zeros(length(algorithms), func_count);

%% Loop over functions
for Function_num = 1:func_count
    Function_name = strcat('F', num2str(Function_num));
    [lb, ub, dim, fobj] = Get_Functions_cec2017(Function_num, dim);
    lb = lb * ones(1, dim);
    ub = ub * ones(1, dim);
    
    % Initialize comparison variables for each algorithm
    Time_compare = cat(3);
    Fival_compare = cat(3);
    curve_compare = cat(3);
    ranks = [];
    rank_wilcoxo_results_compare = [];
    results_compare = cat(3);
    box_plot_compare = cat(3);
    
    name_all = [];

    %% Loop over iterations
    for indexi = 1:iter_count
        % Call each algorithm
        
        fMin=Inf*ones(1,length(algorithms));
        for iter = 1:length(algorithms)
            algorithm_name = algorithms{iter};
            t1 = clock;
            % Call the algorithm
            [fMin_algorithm, ~, curve_algorithm] = feval(algorithm_name, pop_num, Max_iter, lb, ub, dim, fobj);
            t2 = clock;
            time_algorithm = t2(end) + t2(end-1) * 60 + t2(end-2) * 3600 - t1(end) - t1(end-1) * 60 - t1(end-2) * 3600;

            % Store results
            Fival_compare(indexi, iter, :) = fMin_algorithm;
            if(fMin_algorithm<fMin(iter))
                fMin(iter)=fMin_algorithm;
            end
            Time_compare(indexi, iter, :) = time_algorithm;
            curve_compare(indexi, iter, :) = curve_algorithm;
            name_all{1, iter} = algorithm_name;
            ranks(iter, indexi) = fMin_algorithm;
            
            % all_ranks(iter, Function_num) = fMin_algorithm;
        end
        all_ranks(:,Function_num)=fMin;
    end

    %% Wilcoxon rank sum test results
    for i = 2:length(name_all)
        rs = ranksum(ranks(1, :), ranks(i, :));
        if isnan(rs)
            rs = 1;
        end
        rank_wilcoxo_results_compare = [rank_wilcoxo_results_compare, rs];
    end
    rank_wilcoxo_RESULT = [rank_wilcoxo_RESULT; rank_wilcoxo_results_compare];

    %% Storage of data
    disp(Function_name);
    data = [];
    Sheet_name_ms = 'CEC2017';
    xlRange = strcat('C', num2str(Function_num * 3));
    ylRange = strcat('J', num2str(Function_num * 3 + 1));
    range_ms = strcat(xlRange, ':', ylRange);
    disp(range_ms);
    for i = 1:length(name_all)
        box_plot_compare = [box_plot_compare; Fival_compare(:, i)'];
        result_msmm = [min(Fival_compare(:, i)); std(Fival_compare(:, i)); mean(Fival_compare(:, i)); median(Fival_compare(:, i)); max(Fival_compare(:, i))];
        disp([cell2mat(name_all(i)), ': 最优值:', num2str(result_msmm(1), 4), ' 标准差:', num2str(result_msmm(2), 4), ' 平均值:', num2str(result_msmm(3), 4), ' 中值:', num2str(result_msmm(4), 4), ' 最差值:', num2str(result_msmm(5), 4)]);
        results_compare = [results_compare, result_msmm];
        data(:, i) = [mean(Fival_compare(:, i)); std(Fival_compare(:, i))];
    end
    xlswrite('output_data\mean_std.xlsx', data, Sheet_name_ms, range_ms);

    RESULT = [RESULT; results_compare];

    %% Draw box plots
    if box_pp == 1
        subplot(5, 6, Function_num)
        mycolor = [0.862745098039216, 0.827450980392157, 0.117647058823529; ...
                   0.705882352941177, 0.266666666666667, 0.423529411764706; ...
                   0.949019607843137, 0.650980392156863, 0.121568627450980; ...
                   0.956862745098039, 0.572549019607843, 0.474509803921569; ...
                   0.231372549019608, 0.490196078431373, 0.717647058823529; ...
                   0.655541222625894, 0.122226545135785, 0.325468941131211; ...
                   0.766665984235466, 0.955154623456852, 0.755161564478523; ...
                   0.766665984235466, 0.955154623456852, 0.755161564478523];
        box_figure = boxplot(box_plot_compare', 'color', [0 0 0], 'Symbol', 'o');
        set(box_figure, 'Linewidth', 1.2);
        boxobj = findobj(gca, 'Tag', 'Box');
        for i = 1:length(name_all)
            patch(get(boxobj(i), 'XData'), get(boxobj(i), 'YData'), mycolor(i, :), 'FaceAlpha', 0.5, 'LineWidth', 0.7);
        end
        set(gca, 'XTickLabel', name_all);
        title(['F', num2str(Function_num)])
        hold on
    end

    %% Save Wilcoxon rank sum results
    Sheet_name_wil_rank = 'CEC2017';
    range_wil = 'B3:H32';
    xlswrite('output_data\wilcoxo_rank.xlsx', rank_wilcoxo_RESULT, Sheet_name_wil_rank, range_wil);

end

    %% Friedman rank results
    [p, tbl, stats] = friedman(all_ranks');
    disp(['Friedman test p-value for ', Function_name, ': ', num2str(p)]);
    rank_friedman_RESULT = stats.meanranks;

    % Save Friedman rank results
    Sheet_name_friedman_rank = 'CEC2017';
    range_friedman = 'C93:J93';
    xlswrite('output_data\mean_std.xlsx', rank_friedman_RESULT, Sheet_name_friedman_rank, range_friedman);
    

