function error = func(x,net,inputnode,hiddennode,outputnode,inputn_train,outputn_train)
%�Ӻ�������������Ӧ��ֵ��ǰ�򴫲�
%��ȡ����������
%x�б����������㷨ȷ����Ȩֵ����ֵ�������������Ƕ�Ӧ����Ȩֵ��ֵ������
%����������������˼
%w1=x��1��inputnode*hiddennode����ô������
w1=x(1:inputnode*hiddennode);%ȡ������������������ӵ�Ȩֵ
%b2=�ӱߵ����ز���ô�����
B1=x(inputnode*hiddennode+1:inputnode*hiddennode+hiddennode);%��������Ԫ��ֵ
%w2=����ô����㵽���ز㵽�������ô������
w2=x(inputnode*hiddennode+hiddennode+1:inputnode*hiddennode+hiddennode+hiddennode*outputnode);%ȡ������������������ӵ�Ȩֵ
%b2=���ϴεıߵ����ز���ô���
B2=x(inputnode*hiddennode+hiddennode+hiddennode*outputnode+1:inputnode*hiddennode+hiddennode+hiddennode*outputnode+outputnode);%�������Ԫ��ֵ
%����Ȩֵ��ֵ
net.iw{1,1}=reshape(w1,hiddennode,inputnode);%���������ԪȨֵ
net.lw{2,1}=reshape(w2,outputnode,hiddennode);%��������ԪȨֵ
net.b{1}=reshape(B1,hiddennode,1);%���������Ԫ��ֵ
net.b{2}=reshape(B2,outputnode,1); % ��������Ԫ��ֵ
%����ѵ��
net=train(net,inputn_train,outputn_train);
% ѵ�����Ĺ�һ������ֵ
model_out=sim(net,inputn_train);
error=sqrt(mean((outputn_train - model_out).^2));



