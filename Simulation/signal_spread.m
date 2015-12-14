function output = signal_spread( input,c )
% spread the en_output by spreading sequence "c"
% Copyright Nov. 2003, babyyong
% CSE lab, Yifu Science Building, Fudan Univ.

c_length = length(c);
[K,data_length] = size(input);

for k = 1:K
    temp = input(k,:)'*c;%
    for i = 1:data_length
        output(k,(i-1)*c_length+1:i*c_length) = temp(i,:);
    end
    clear temp
end