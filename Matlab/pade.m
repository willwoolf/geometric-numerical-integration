%% pade

A = [
    -121 60
    -160 79
];


function MAT = mypade(n,m,A)
MAT = denominator(n,m,A)\numerator(n,m,A);
end


function MAT = numerator(n,m,A)
MAT = zeros(size(A));
for j = 0:n
    num = factorial(n + m - j)*factorial(n);
    den = factorial(n+m)*factorial(j)*factorial(n-j);
    MAT = MAT + (num/den)*mpower(A,j);
end
end

function MAT = denominator(n,m,A)
MAT = zeros(size(A));
for j = 0:m
    num = factorial(n + m - j)*factorial(m);
    den = factorial(n+m)*factorial(j)*factorial(m-j);
    MAT = MAT + (num/den)*mpower((-A),j);
end
end

