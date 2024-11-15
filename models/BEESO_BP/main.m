%% inilization
warning off
close all
clear
clc

%% load data
data=xlsread('data_1-7_cell_count.xlsx');

%% Split the training and test set
temp = randperm(19591);  
inputTrainDataset = data(temp(1:2000), 1:12)';
outputTrainDataset = data(temp(1:2000), 13)';

inputTestDataset = data(temp(19091:19591), 1:12)';
outputTestDataset = data(temp(19091:19591), 13)';

%% Normalization
[inputn_train, input_ps] = mapminmax(inputTrainDataset, 0, 1);
inputn_test = mapminmax('apply', inputTestDataset, input_ps);
[outputn_train, output_ps] = mapminmax(outputTrainDataset, 0, 1);

%% Determine the neural network architecture
inputnode = length(inputn_train(:, 1));
outputnode = 1;

% Determine the number of neuron nodes in the hidden layer
bound = 3 : 1 : inputnode;
cnt = 1;
for k = bound
    net=newff(inputn_train,outputn_train,k,{'tansig','purelin'},'trainlm');
    net.trainParam.epochs=600;
    net.trainParam.lr=0.01;
    net.trainParam.goal=1e-5;
    net.trainParam.showWindow = 0;
    net=train(net,inputn_train,outputn_train);
    model_out = sim(net, inputn_train);
    e(cnt) = mean((model_out - outputn_train).^2);
    
    cnt = cnt + 1;
end
hiddennode = bound(find(e == min(e), 1));
disp(['The optimal number of hidden layer nodes is: ', num2str(hiddennode), ...
    ', the corresponding normalized RMSE value of the training set is the smallest, as follows. ', ...
    num2str(e((find(e == min(e), 1))))])

%% Calling the optimization algorithm
model_name = sprintf('Currently running, BEESO-BP neural network prediction program.');
disp(model_name);
disp('running... ...')
% tic; 
net=newff(inputn_train,outputn_train,hiddennode,{'tansig','purelin'},'trainlm');
net.trainParam.epochs=5;
net.trainParam.lr=0.01;
net.trainParam.goal=1e-5;
net.trainParam.showWindow=0;
disp(['net.iw{1,1} count is=',num2str(numel(net.iw{1,1}))]);
maxgen=50;
popsize=100;
dim=inputnode*hiddennode+hiddennode+hiddennode*outputnode+outputnode;
lb=-3 * ones(1, dim);
ub=3 * ones(1, dim);
fobj=@(x)func(x,net,inputnode,hiddennode,outputnode,inputn_train,outputn_train);
[curve, optimized_param]=BEESO(fobj,popsize,maxgen,lb,ub,dim, 1);

figure
plot(curve, 'r-', 'LineWidth', 1.0)
grid on
xlabel('Algebra of evolution')
ylabel('Best fitness')
title('Curve of evolution')

%% The model is trained using the optimized parameters

w1=optimized_param(1:inputnode*hiddennode);
B1=optimized_param(inputnode*hiddennode+1:inputnode*hiddennode+hiddennode);
w2=optimized_param(inputnode*hiddennode+hiddennode+1:inputnode*hiddennode+hiddennode+hiddennode*outputnode);
B2=optimized_param(inputnode*hiddennode+hiddennode+hiddennode*outputnode+1:inputnode*hiddennode+hiddennode+hiddennode*outputnode+outputnode);
net.iw{1,1}=reshape(w1,hiddennode,inputnode);
net.lw{2,1}=reshape(w2,outputnode,hiddennode);
net.b{1}=reshape(B1,hiddennode,1);
net.b{2}=reshape(B2,outputnode,1);

%% training
net.trainParam.showWindow=1;
net=train(net,inputn_train,outputn_train);

%% Prediction and de-normalization
model_out1 = sim(net, inputn_train);  % 训练集的归一化预测结果
model_out2 = sim(net, inputn_test);    % 测试集的归一化预测结果
predictTrainDataset = mapminmax('reverse', model_out1, output_ps);  % 反归一化训练集预测结果为原始数量级
predictTestDataset = mapminmax('reverse', model_out2, output_ps);    % 反归一化测试集预测结果为原始数量级
toc;
%% Error of analysis
disp(' ')
disp('The training set error is calculated as follows. ')
MSE = mean((outputTrainDataset - predictTrainDataset).^2);
disp(['MSE = ', num2str(MSE)])
MAE = mean(abs(outputTrainDataset - predictTrainDataset));
disp(['MAE = ', num2str(MAE)])
RMSE = sqrt(MSE);
disp(['RMSE = ', num2str(RMSE)])
MAPE = mean(abs((outputTrainDataset - predictTrainDataset)./outputTrainDataset));
disp(['MAPE = ', num2str(MAPE*100), '%'])
R = corrcoef(outputTrainDataset, predictTrainDataset);
R2 = R(1, 2)^2;
disp(['R2 = ', num2str(R2)])
disp(' ')
disp('The test set error is calculated as follows. ')
MSE_test = mean((outputTestDataset - predictTestDataset).^2);
disp(['MSE = ', num2str(MSE_test)])
MAE_test = mean(abs(outputTestDataset - predictTestDataset));
disp(['MAE = ', num2str(MAE_test)])
RMSE_test = sqrt(MSE_test);
disp(['RMSE = ', num2str(RMSE_test)])
MAPE_test = mean(abs((outputTestDataset - predictTestDataset)./outputTestDataset));
disp(['MAPE = ', num2str(MAPE_test*100), '%'])
R_test = corrcoef(outputTestDataset, predictTestDataset);
R2_test = R_test(1, 2)^2;
disp(['R2 = ', num2str(R2_test)])

