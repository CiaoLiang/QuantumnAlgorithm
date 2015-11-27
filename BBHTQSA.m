function[x,Lbbht]=BBHTQSA(sigma,nbit,g,fn)
nums = 2^nbit;
m = 1;
ite = 0;
inistate = (1/nums)^(1/2)*ones(nums,1);
lamada = 6/5; 
Lbbht = 0;
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
    %fprintf('Lbbht: %d \n',Lbbht);
    if fn(x) < sigma | Lbbht > 4.5*nums^(1/2); %搜寻到值或者迭代次数过多
        %fprintf('%d iterations : search index: %d\n',ite,x);
        break;
    else
        m = min(m*lamada,nums^(1/2));
        %fprintf('m: %d \n',m);
    end
end