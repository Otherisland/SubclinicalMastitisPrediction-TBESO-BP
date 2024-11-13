function[inputn_train,inputn_test,outputn_train,output_ps,outputTrainDataset,outputTestDataset] =load_data(path,n,data_n,mid,row,col)
%% load data
data=xlsread(path);

%% Split the training and test set
temp = randperm(19591);

inputTrainDataset = data(temp(1:2000), 1:12)';
outputTrainDataset = data(temp(1:2000), 13)';

inputTestDataset = data(temp(19091:19591), 1:12)';
outputTestDataset = data(temp(19091 : 19591),13)';
%% Normalization
[inputn_train, input_ps] = mapminmax(inputTrainDataset, 0, 1);
inputn_test = mapminmax('apply', inputTestDataset, input_ps);
[outputn_train, output_ps] = mapminmax(outputTrainDataset, 0, 1);

save("inputn_train.mat","inputn_train");
save("inputn_test.mat","inputn_test");
save("outputn_train.mat","outputn_train");

save("outputTrainDataset.mat","outputTrainDataset");
save("outputTestDataset.mat","outputTestDataset");

end