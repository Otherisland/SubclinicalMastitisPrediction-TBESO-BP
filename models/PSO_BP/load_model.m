%% nilization
warning off
close all
clear
clc;

%% load data
[inputn_train,inputn_test,outputn_train,output_ps,outputTrainDataset,outputTestDataset]=load_data("data_1-7_cell_count.xlsx",19591,2000,1900,13,14);

%traning
[net,predictTrainDataset,predictTestDataset]=PSO_BP(inputn_train,inputn_test,outputn_train,output_ps);
%print the results
figure_data(outputTrainDataset,predictTrainDataset,outputTestDataset,predictTestDataset)