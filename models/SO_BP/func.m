function error = func(x,net,inputnode,hiddennode,outputnode,inputn_train,outputn_train)
%子函数用来计算适应度值，前向传播
%提取超参数变量
%x中保存用其他算法确定的权值和阈值，括号里计算的是对应各层权值阈值的索引
%这个括号是数组的意思
%w1=x从1到inputnode*hiddennode，这么多条边
w1=x(1:inputnode*hiddennode);%取到输入层与隐含层连接的权值
%b2=从边到隐藏层那么多个点
B1=x(inputnode*hiddennode+1:inputnode*hiddennode+hiddennode);%隐含层神经元阈值
%w2=从那么多个点到隐藏层到输出层那么多条边
w2=x(inputnode*hiddennode+hiddennode+1:inputnode*hiddennode+hiddennode+hiddennode*outputnode);%取到隐含层与输出层连接的权值
%b2=从上次的边到隐藏层那么多点
B2=x(inputnode*hiddennode+hiddennode+hiddennode*outputnode+1:inputnode*hiddennode+hiddennode+hiddennode*outputnode+outputnode);%输出层神经元阈值
%网络权值赋值
net.iw{1,1}=reshape(w1,hiddennode,inputnode);%隐含层的神经元权值
net.lw{2,1}=reshape(w2,outputnode,hiddennode);%输出层的神经元权值
net.b{1}=reshape(B1,hiddennode,1);%隐含层的神经元阈值
net.b{2}=reshape(B2,outputnode,1); % 输出层的神经元阈值
%网络训练
net=train(net,inputn_train,outputn_train);
% 训练集的归一化仿真值
model_out=sim(net,inputn_train);
error=sqrt(mean((outputn_train - model_out).^2));



