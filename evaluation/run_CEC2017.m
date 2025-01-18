clc
clear
close all
%%
nPop=50; % ��Ⱥ��

Max_iter=500; % ����������

dim = 100; % ��ѡ 2, 10, 30, 50, 100

%%  ѡ����

Function_name=11; % �������� 1 - 30
[lb,ub,dim,fobj] = Get_Functions_cec2017(Function_name,dim);

%% �����㷨
tic
[Best_score,Best_pos,cg_curve]=WOA(nPop,Max_iter,lb,ub,dim,fobj);

% [Convergence_curve_SO,Best_pos_SO,Best_score_SO]=SO(fobj,POS,SearchAgents_no,Max_iteration,lb,ub,dim,func_num);
% [Convergence_curve_TBESO,Best_pos_TBESO,Best_score_TBESO]=TBESO(fobj,tent_POS,SearchAgents_no,Max_iteration,lb,ub,dim,func_num);
% % [Convergence_curve_,Best_pos,Best_score]=SO(fobj,tent_POS,SearchAgents_no,Max_iteration,lb,ub,dim,func_num);
% 
% [Convergence_curve_Tentchaotic_SO,Best_pos_Tentchaotic_SO,Best_score_Tentchaotic_SO]=T_SO(fobj,tent_POS,SearchAgents_no,Max_iteration,lb,ub,dim,func_num);
% [Convergence_curve_BiSearch_SO,Best_pos_BiSearch_SO,Best_score_BiSearch_SO]=B_SO(fobj,POS,SearchAgents_no,Max_iteration,lb,ub,dim,func_num);
% [Convergence_curve_EOBL_SO,Best_pos_EOBL_SO,Best_score_EOBL_SO]=E_SO(fobj,POS,SearchAgents_no,Max_iteration,lb,ub,dim,func_num);
% 
% [Convergence_curve_TB_SO,Best_pos_TB_SO,Best_score_TB_SO]=TB_SO(fobj,tent_POS,SearchAgents_no,Max_iteration,lb,ub,dim,func_num);
% [Convergence_curve_TE_SO,Best_pos_TE_SO,Best_score_TE_SO]=TE_SO(fobj,tent_POS,SearchAgents_no,Max_iteration,lb,ub,dim,func_num);
% [Convergence_curve_BE_SO,Best_pos_BE_SO,Best_score_BE_SO]=BE_SO(fobj,POS,SearchAgents_no,Max_iteration,lb,ub,dim,func_num);
    
toc

%% plot
figure('Position',[400 200 300 250])
semilogy(cg_curve,'Color','r','Linewidth',1)
%     plot(cg_curve,'Color','r','Linewidth',1)
title(['Convergence curve, Dim=' num2str(dim)])
xlabel('Iteration');
ylabel(['Best score F' num2str(Function_name) ]);
axis tight
grid on
box on
set(gca,'color','none')
legend('WOA')

