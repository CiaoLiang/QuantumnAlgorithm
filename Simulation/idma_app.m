function [L_e_2,estimate_d_k] = idma_app( L_a_2,c)
% function [L_e_2,estimate_d_k,estimate_c_k] = idma_app( L_a_2,c)

% L_a2 the kth user's a priori info.
% Copyright Apr. 2013, L
% CSE lab, Coumputing Center, Fudan Univ.

c_length = length(c);
LLR_limit = 700;

% Total number of bits
data_length = length ( L_a_2 );
L_total = data_length/c_length;



temp = repmat(c,1,L_total);
L_a_temp = L_a_2 .* temp;
clear temp
for i = 1:L_total
    L_a_k(1,i) = sum(L_a_temp(1,(i-1)*c_length+1:i*c_length),2);
end

L_e_2 = signal_spread( L_a_k,c );

% estimate_c_k = sign(L_e_2);
% estimate_c_k = ((estimate_c_k+1)/2);

L_e_2 = L_e_2 - L_a_2;

for i=1:data_length
    if L_e_2(1,i) > LLR_limit
        L_e_2(1,i) = LLR_limit;
    elseif L_e_2(1,i) < -LLR_limit
        L_e_2(1,i) = -LLR_limit;
    end
end

estimate_d_k = sign(L_a_k);
estimate_d_k = round((estimate_d_k+1)/2);






