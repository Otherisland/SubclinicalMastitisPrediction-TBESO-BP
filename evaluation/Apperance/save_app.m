% 创建line_styles变量
% 
% currentPath = pwd;
% disp(['Current working directory: ' currentPath]);

%line_styles
line_styles = {'-', '--', ':', '-.'}; %'none'
line_widths=[1,1,1,1,1,1,1,2.5];
% markers={'o','x','+','*','s','d','^','>','v','<','p','h'};%'none','.','_','|'
markers={'*','s','p','^','h','+','x','o'};
color_all= [    
    33/255, 126/255, 195/255;
    138/255, 67/255, 153/255;
    243/255, 204/255, 110/255;
    152/255, 192/255, 100/255;
    84/255, 192/255, 238/255;
    185/255, 77/255, 100/255;
    133/255, 188/255, 255/255;
    217/255, 82/255, 24/255;
    ];  % 明显的红色


% 创建 Apperance 文件夹（如果不存在）
if ~exist('Apperance', 'dir')
    mkdir('Apperance');
end

% 保存 line_styles 和 line_width 到 cec2017/Apperance/style_info.mat 文件
save('Apperance/style_info.mat', 'line_styles', 'line_widths','markers','color_all');
% save('Apperance/style_info.mat', 'line_styles', 'line_widths','markers');

