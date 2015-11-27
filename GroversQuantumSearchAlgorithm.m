nbit = 5; %���ӱ�����Ŀ 
nums = 2^nbit ;%����̬��Ŀ
iter = 140; %������Ŀ
fn = random('norm',1,1,nums,1) ;%cost function 
index = 3; %��Ѱ

P = eye(nums); %phase shift gate
P(1,1) = -1;

H = (1/2)^(1/2)*[1,1;1,-1]; %Hadamard gate
for i = 2:nbit
    H = (1/2)^(1/2)*[H,H;H,-H];
end 

O = eye(nums); %Oracle gate
O(index,index) = -1;

%GroversQuantumSearchAlogrithm
%��ʼ̬
state = (1/nums)^(1/2)*ones(nums,1);
pro = [];
g = H*P*H*O; %grave���� 
for i = 1:iter
    state =   g*state; 
    probability = (state(index))^2;
    pro = [pro,probability];
    fprintf('Search Probability : %f\n',probability);
end 

plot([1:iter],pro);
