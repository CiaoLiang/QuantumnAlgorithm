nbit = 10; %���ӱ�����Ŀ,10bit���ϻ�ǳ���
nums = 2^nbit ;%����̬��Ŀ
fn = random('norm',1,1,nums,1) ;%cost function 
index = 9; %��Ѱ��index

P = eye(nums); %phase shift gate
P(1,1) = -1;

H = hadamard(nums);
H = H./norm(H);

O = eye(nums); %Oracle gate
O(index,index) = -1; %ʵ��д��ʱ��Ҫ����CF����д

%GroversQuantumSearchAlogrithm
%��ʼ̬
inistate = (1/nums)^(1/2)*ones(nums,1);
g = H*P*H*O; %grave���Ӳ�Ҫд��

m = 1; %initial
sigma = fn(index);
lamada = 6/5; 
Lbbht = 0;
x = randi(nums,1);
ite = 0; %��������
while true;
    ite = ite + 1;
    L = randi(floor(m));
    Lbbht = Lbbht + L;
    state = inistate;
    while L~=0;
        state = g*state;
        L = L-1;
    end
    state = state./norm(state); %ǿ����1�����⸡�����
    probability = state.^2'; %���ʷֲ�Ϊ1
    x = randsrc(1,1,[[1:nums];probability]); %observe Quantum
    fprintf('Lbbht: %d \n',Lbbht);
    if fn(x) == sigma | Lbbht > 4.5*nums^(1/2); %��Ѱ��ֵ���ߵ�����������
        fprintf('%d iterations : search index: %d\n',ite,x);
        break;
    else
        m = min(m*lamada,nums^(1/2));
        fprintf('m: %d \n',m);
    end
end

%fprintf('%d iterations : search index: %d',ite,x);