%% 
A = [1 3;2 4]*[-1 0; 0 -41]*inv([1 3;2 4])
trueexpA = [1 3;2 4]*[exp(-1) 0; 0 exp(-41)]*inv([1 3;2 4])
cA = eye(2) + A
prevA = zeros(2,2);
i = 2;

while not(prevA == cA)
    prevA = cA;
    cA = cA + (mpower(A,i))./factorial(i);
    if any(cA(abs(cA) >= 1e16))
        k = i
    end
    i = i+1;
end