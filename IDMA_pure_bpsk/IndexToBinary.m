function [ index ] = IndexToBinary(a,Qubit)
%This Function Change Index to BinarySequence to calculate the cost
%function
%   Input a,nums
%   Output BinarySequence
    index = zeros(1,Qubit);
    for i = 1:Qubit
        index(i) = 2*mod(a,2)-1;
        a = floor(a/2);
    end
end

