nbit = 10; %量子比特数目,10bit以上会非常卡
nums = 2^nbit ;%量子态数目
fn = random('norm',1,1,nums,1) ;%cost function 
index = 9; %搜寻的index

P = eye(nums); %phase shift gate
P(1,1) = -1;

H = hadamard(nums);
H = H./norm(H);

O = eye(nums); %Oracle gate
O(index,index) = -1; %实际写的时候要根据CF来书写

%GroversQuantumSearchAlogrithm
%初始态
inistate = (1/nums)^(1/2)*ones(nums,1);
g = H*P*H*O; %grave算子不要写错

m = 1; %initial
sigma = fn(index);
lamada = 6/5; 
Lbbht = 0;
x = randi(nums,1);
ite = 0; %迭代次数
while true;
    ite = ite + 1;
    L = randi(floor(m));
    Lbbht = Lbbht + L;
    state = inistate;
    while L~=0;
        state = g*state;
        L = L-1;
    end
    state = state./norm(state); %强行置1，避免浮点误差
    probability = state.^2'; %概率分布为1
    x = randsrc(1,1,[[1:nums];probability]); %observe Quantum
    fprintf('Lbbht: %d \n',Lbbht);
    if fn(x) == sigma | Lbbht > 4.5*nums^(1/2); %搜寻到值或者迭代次数过多
        fprintf('%d iterations : search index: %d\n',ite,x);
        break;
    else
        m = min(m*lamada,nums^(1/2));
        fprintf('m: %d \n',m);
    end
end

%fprintf('%d iterations : search index: %d',ite,x);