function [L_SISO_2,estimate_d_k] = idma_app( L_a2,c)
% L_a2   the kth user's a priori info.  
% Copyright Nov. 2003, babyyong
% CSE lab, Yifu Science Building, Fudan Univ.

c_length = length(c);

% Total number of bits
data_length = length ( L_a2 );
L_total = data_length/c_length;

temp = repmat(c,1,L_total);%temp :1*(clength*ltotal)
L_a_temp = L_a2 .* temp;
clear temp
for i = 1:L_total
    L_a_k(1,i) = sum(L_a_temp(1,(i-1)*c_length+1:i*c_length),2);
end

L_SISO_2 = signal_spread( L_a_k,c ); 
L_SISO_2 = L_SISO_2 - L_a2;

estimate_d_k = sign(L_a_k);
estimate_d_k = (estimate_d_k+1)/2;

[m,n] = size(estimate_d_k);
index = [1,0];
for i = 1:m
    for j = 1:n
        if estimate_d_k(i,j)==0.5
            estimate_d_k(i,j) = index(randi(2));
        end
    end
end