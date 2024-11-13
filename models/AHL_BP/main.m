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

%% Establishing the model
inputnode = length(inputn_train(:, 1));
outputnode = 1;
tic;
net=newff(inputn_train,outputn_train,hiddennode,{'tansig','purelin'},'trainlm');
net.trainParam.epochs=5;
net.trainParam.lr=0.01;
net.trainParam.goal=1e-5;
%%training
net=train(net,inputn_train,outputn_train);

% Prediction and de-normalization
model_out1 = sim(net, inputn_train);
model_out2 = sim(net, inputn_test);
predictTrainDataset = mapminmax('reverse', model_out1, output_ps);
predictTestDataset = mapminmax('reverse', model_out2, output_ps);
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

