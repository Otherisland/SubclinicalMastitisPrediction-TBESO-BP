%%绘制图像
function[] = figure_data(output_train,predictTrainDataset,output_test,predictTestDataset)

%% 分析误差
disp(' ')
disp('训练集误差计算如下: ')
MSE = mean((output_train - predictTrainDataset).^2);
disp(['均方误差MSE = ', num2str(MSE)])
MAE = mean(abs(output_train- predictTrainDataset));
disp(['平均绝对误差MAE = ', num2str(MAE)])
RMSE = sqrt(MSE);
disp(['根均方误差RMSE = ', num2str(RMSE)])
MAPE = mean(abs((output_train - predictTrainDataset)./output_train));
disp(['平均绝对百分比误差MAPE = ', num2str(MAPE*100), '%'])
R = corrcoef(output_train, predictTrainDataset);
R2 = R(1, 2)^2;
disp(['拟合优度R2 = ', num2str(R2)])
disp(' ')
disp('测试集误差计算如下: ')
MSE_test = mean((output_test - predictTestDataset).^2);
disp(['均方误差MSE = ', num2str(MSE_test)])
MAE_test = mean(abs(output_test  - predictTestDataset));
disp(['平均绝对误差MAE = ', num2str(MAE_test)])
RMSE_test = sqrt(MSE_test);
disp(['根均方误差RMSE = ', num2str(RMSE_test)])
MAPE_test = mean(abs((output_test  - predictTestDataset)./output_test ));
disp(['平均绝对百分比误差MAPE = ', num2str(MAPE_test*100), '%'])
R_test = corrcoef(output_test , predictTestDataset);
R2_test = R_test(1, 2)^2;
disp(['拟合优度R2 = ', num2str(R2_test)])

%% 对结果作图
% 训练集
% figure
% plot(output_train, 'b*-', 'LineWidth', 0.8)
% hold on
% plot(output_train, 'ro-', 'LineWidth', 0.8)
% grid on
% xlabel('训练样本序号')
% ylabel('目标')
% legend('实际值', '预测值')
% title({'BP神经网络预测训练集预测值和实际值对比图', ['根均方误差RMSE = ', num2str(RMSE), '拟合优度R2 = ', num2str(R2)]})
% 
% figure
% plot(output_train - predictTrainDataset, 'b*-', 'LineWidth', 0.8)
% grid on
% xlabel('训练样本序号')
% ylabel('预测偏差')
% legend('误差')
% title({'BP神经网络预测训练集预测误差图', ['平均绝对百分比误差MAPE = ', num2str(MAPE*100), '%']})
% 
% % 测试集
% figure
% plot(output_test, 'b*-', 'LineWidth', 0.8)
% hold on
% plot(predictTestDataset, 'ro-', 'LineWidth', 0.8)
% grid on
% xlabel('测试样本序号')
% ylabel('目标')
% legend('实际值', '预测值')
% title({'BP神经网络预测测试集预测值和实际值对比图', ['根均方误差RMSE = ', num2str(RMSE_test), '拟合优度R2 = ', num2str(R2_test)]})
% 
% figure
% plot(output_test - predictTestDataset, 'b*-', 'LineWidth', 0.8)
% grid on
% xlabel('测试样本序号')
% ylabel('预测偏差')
% legend('误差')
% title({'BP神经网络预测测试集预测误差图', ['平均绝对百分比误差MAPE = ', num2str(MAPE_test*100), '%']})
% end