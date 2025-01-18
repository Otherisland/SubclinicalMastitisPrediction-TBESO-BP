%% PSO-BP
function [net,predictTrainDataset,predictTestDataset] =PSO_BP(inputn,inputn_test,outputn,outputps)
%% Determine the neural network architecture
inputnum=length(inputn(:, 1)); 

hiddennum=5;
outputnum=1;

net=newff(inputn,outputn,hiddennum);
 
%% Parameter inilization
tic;
c1 = 1.49445;
c2 = 1.49445;
 
maxgen=50; 
sizepop=100;
 
Vmax=1;Vmin=-1;
popmax=5;
popmin=-5;
 
for i=1:sizepop
    pop(i,:)=inputnum*rands(1,inputnum*hiddennum*2+1);
    V(i,:)=rands(1,inputnum*hiddennum*2+1);
    fitness(i)=fun(pop(i,:),inputnum,hiddennum,outputnum,net,inputn,outputn);
end
 
 
[bestfitness bestindex]=min(fitness);
zbest=pop(bestindex,:);
gbest=pop;
fitnessgbest=fitness;
fitnesszbest=bestfitness;
 
%% Iterative optimization
for i=1:maxgen
    for j=1:sizepop
        V(j,:) = V(j,:) + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:));
        V(j,find(V(j,:)>Vmax))=Vmax;
        V(j,find(V(j,:)<Vmin))=Vmin;
        
        pop(j,:)=pop(j,:)+0.2*V(j,:);
        pop(j,find(pop(j,:)>popmax))=popmax;
        pop(j,find(pop(j,:)<popmin))=popmin;
        
        pos=unidrnd(21);
        if rand>0.95
            pop(j,pos)=5*rands(1,1);
        end
      
        fitness(j)=fun(pop(j,:),inputnum,hiddennum,outputnum,net,inputn,outputn);
    end
    
    for j=1:sizepop
    if fitness(j) < fitnessgbest(j)
        gbest(j,:) = pop(j,:);
        fitnessgbest(j) = fitness(j);
    end
    
    if fitness(j) < fitnesszbest
        zbest = pop(j,:);
        fitnesszbest = fitness(j);
    end
    
    end
    
    yy(i)=fitnesszbest;    
        
end
 
x=zbest;
%%  The optimal initial threshold weight is assigned to the network prediction
w1=x(1:inputnum*hiddennum);
B1=x(inputnum*hiddennum+1:inputnum*hiddennum+hiddennum);
w2=x(inputnum*hiddennum+hiddennum+1:inputnum*hiddennum+hiddennum+hiddennum*outputnum);
B2=x(inputnum*hiddennum+hiddennum+hiddennum*outputnum+1:inputnum*hiddennum+hiddennum+hiddennum*outputnum+outputnum);
 
net.iw{1,1}=reshape(w1,hiddennum,inputnum);
net.lw{2,1}=reshape(w2,outputnum,hiddennum);
net.b{1}=reshape(B1,hiddennum,1);
net.b{2}=B2;
 
%% BP
net.trainParam.epochs=10;
net.trainParam.lr=0.1;
 
[net,per2]=train(net,inputn,outputn);

%% predicting

model_out1 = sim(net, inputn);
model_out2 = sim(net, inputn_test);
predictTrainDataset = mapminmax('reverse', model_out1, outputps);
predictTestDataset = mapminmax('reverse', model_out2, outputps);
toc;

%% restore
save("pso_bp.mat","net");
save("predictTrainDataset.mat","predictTrainDataset");
save("predictTestDataset.mat","predictTestDataset");
save("output_ps.mat","outputps")
end

