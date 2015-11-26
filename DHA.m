%DHA alogrithm
%Simulation the xmin = arg min{f(x)}
nbit = 12; %the number of Qubit
nums = 2^nbit ;%the number of state
fn = random('norm',1,1,nums,1) ;%cost function 

P = eye(nums); %phase shift gate
P(1,1) = -1;

H = (1/2)^(1/2)*[1,1;1,-1]; %Hadamard gate
for i = 2:nbit
    H = (1/2)^(1/2)*[H,H;H,-H];
end 

inistate = (1/nums)^(1/2)*ones(nums,1);

xi = randi(nums,1);
sigma = fn(xi);
xs = 0;
Ltotal = 0;
Lbbht = 0;
xmin = 0;
while true;
    fprintf('iteration start>>>\n');
    sigma = fn(xi);
    O = eye(nums);
    for i = 1:nums;
        if fn(i)<sigma
            O(i,i) = -1;
        end
    end
    g = H*P*H*O;%grave
    [xs,Lbbht]=BBHTQSA(sigma,nbit,g,fn);%BBHT Quantum Search Algorithm
    Ltotal = Ltotal + Lbbht;
    fprintf('Ltotal:%d\n',Ltotal);
    if fn(xs)>=fn(xi) | Ltotal>22.5*nums^(1/2);
        xmin = xi;
        break;
    else
        xi = xs;
    end
end
fprintf('min f(x):f[%d]:%f\n',xmin,fn(xmin));
fprintf('the min number is :%f\n',min(fn));
