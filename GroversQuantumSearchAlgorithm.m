nbit = 5; %量子比特数目 
nums = 2^nbit ;%量子态数目
iter = 140; %迭代数目
fn = random('norm',1,1,nums,1) ;%cost function 
index = 3; %搜寻

P = eye(nums); %phase shift gate
P(1,1) = -1;

H = (1/2)^(1/2)*[1,1;1,-1]; %Hadamard gate
for i = 2:nbit
    H = (1/2)^(1/2)*[H,H;H,-H];
end 

O = eye(nums); %Oracle gate
O(index,index) = -1;

%GroversQuantumSearchAlogrithm
%初始态
state = (1/nums)^(1/2)*ones(nums,1);
pro = [];
g = H*P*H*O; %grave算子 
for i = 1:iter
    state =   g*state; 
    probability = (state(index))^2;
    pro = [pro,probability];
    fprintf('Search Probability : %f\n',probability);
end 

plot([1:iter],pro);
