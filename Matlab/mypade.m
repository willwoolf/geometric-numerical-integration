%% pade counterexample

A = [
    -4 1 0;
    2 -1 2;
    2 0 -2
]

firstorderapproximation = inv(eye(3) - A)

eig(firstorderapproximation)

x = [0.75; 0.25; 0.5];

pade_denominator = inv(denominator(1,1,A))
pade_numerator = numerator(1,1,A)
trueexpA = expm(A)

a = -max(max(abs(diag(diag(A)))))
ThrA = A - a*eye(size(A))

pade_threshden = inv(denominator(1,1,ThrA));
pade_threshnum = numerator(1,1,ThrA);

secondpade = denominator(1,1,A)\numerator(1,1,A)

y_true = trueexpA*x
y_brute = secondpade*x
y_optim = padeproduct(1,1,A,x,0)

%% sanity test


N = 7
n = 1:N


samples = 500

Y = zeros(N,1)
for k = 1:N
    E = [];
    for iter = 1:samples
        A = GL(10);
        x = rand(10,1);
        
        y = expm(A)*x;
        ym = positiveSeriesExp(n(k),A,0.4)*x;
        E = [E norm(y-ym,2)/norm(y,2)];
    end
    Y(k) = mean(E);
end

plot(n,Y)


%% order of convergence test

D = 8;
d = zeros(1,D);
d(1) = 1;
for k = 2:D
    d(k) = 2*d(k-1);
end
H = 1./d;

samples = 500;

Y = zeros(1,D)
for r = 1:D
    h = H(r)
    E = []
    for iter = 1:samples
        A = GL(4);
        x = rand(4,1);

%         m_approx = positiveSeriesExp(1,h*A,0.5);
        y = expm(h*A)*x;
%         ym = m_approx*x;
        ym = superPade(h*A,x);

        E = [E norm(y-ym,2)/norm(y,2)];
    end
    Y(r) = mean(E);
end

loglog(H,Y)


%% further pade counterexample

A = [
    -4 1 0;
    2 -1 2;
    2 0 -2
]
x = [3; 1; 2];

S = linspace(0,1,21)
Y = []
for s = linspace(0,1,21)
    y = expm(A)*x;
    ym = positiveSeriesExp(2,A,0.35)*x;
    Y = [Y norm(y-ym,2)];
end
semilogy(S,Y)

% ym = positiveSeriesExp(A,4,s)*x;
% y = expm(A)*x;

%% compare the true exponential to the Pade estimator
N = 10;
A = rand(N);
A = A - diag(diag(A));
for k = 1:N
    A(k,k) = -sum(A(:,k));
end

x = rand(N,1)

h = rand()
eig(A)
eig(inv(eye(N)-h*A))

trueexpA = expm(A);
y = trueexpA*x; % should preserve positivity
y_approx = padeproduct(1,1,A,x,-0.1) % hmm
error = norm(y - y_approx,2)

%% compare superpade

% A = GL(10)
% a = -max(max(abs(diag(diag(A)))))
% thrA = A - a*eye(size(A))
% [d,~] = size(A)
% 
% scale = norm(thrA,1)/((2^1)*(1-a))
% 
% p2 = ceil(log2(sqrt(d)*abs(a)/(1-a)))-1
% inv(denominator(1,1,thrA./p));

A = GL(5);
x = rand(5,1);
y = expm(A)*x
yp = superPade(A,x)


%% Another Monte Carlo

err = []
N = 10
E = zeros(1,N)
for order = 1:N
    M = 10;
    err = zeros(M,1);
    for samples = 1:100
        A = h*rand(M);
        A = A - diag(diag(A));
        for k = 1:M
            A(k,k) = -sum(A(:,k));
        end
        x = rand(M,1);
        
        y = expm(A)*x;

%         [n,m] = nsplit(order);
        
%         ym = padeproduct(n,m,A,x,0);
        ym = positiveSeriesExp(order,A,0.5)*x;
        if any(ym(ym<0))
            sprintf("NOOOOOOOO")
        end
        err(samples) = norm(y-ym,2);
    end
    E(order) = mean(err);
end
X = 1:N
plot(X, E)


%% Monte Carlo Estimator for behaviour of approximations

divs = 501;
startval = 0;
endval = 2;
SCL = linspace(startval,endval,divs);
E = zeros(divs,1);
samples = 500;
for S = 1:divs
    s = SCL(S);
    err = zeros(1,samples);
    for M = 1:samples
        % generate a 10x10 linear system
        N = 10;
        A = GL(N);
        x = rand(N,1);

        y = expm(A)*x;
        % s is used here
%         y_approx = padeproduct(10,10,A,x,s);
%         y_approx = positiveSeriesExp(3,A,s)*x;
        y_approx = superPade(A,x);
        
        % error
        err(M) = norm(y-y_approx, 2);

        % negativity
%         err(M) = size(y_approx(y_approx<0),1)/size(y,1);
    end
    E(S) = mean(err);
end
semilogy(SCL,E)

%% Monte Carlo Estimator for convergence of approximations

divs = 10;
E = zeros(divs,1);
H = E;
samples = 500;
for S = 1:divs
    lambda = 2^(-S);
    err = zeros(1,samples);
    for M = 1:samples
        % generate a 10x10 linear system
        N = 10;
        A = GL(N);
        x = rand(N,1);

        y = expm(lambda*A)*x;
%         y_approx = padeproduct(10,10,A,x,s);
%         y_approx = positiveSeriesExp(2,lambda*A,0)*x;
        y_approx = superPade(lambda*A,x);
%         y_approx = (eye(N) - lambda*A)\x;
        
        % error
        err(M) = norm(y-y_approx, 2)/norm(y, 2);

        % negativity
%         err(M) = size(y_approx(y_approx<0),1)/size(y,1);
    end
    E(S) = mean(err);
    H(S) = lambda;
end
loglog(H,E)

%% Monte Carlo estimator AGAIN: f(order, offset) = error

N = 10;
endval = 5;
divs = 101;
F = zeros(N,divs);
samples = 500;

for n = 1:N
    for S = 1:divs
        s = (S-1)*endval/(divs-1);
        ERR = zeros(samples,1);
        for k = 1:samples
            % make a graph-laplacian random system
            A = GL(10);
            x = rand(10,1);
    
            % evaluate y and the series approximation
            y = expm(A)*x;
            ym = positiveSeriesExp(n,A,s)*x;
            
            % evaluate error
            ERR(k) = norm(y-ym,2);
        end
        F(n,S) = mean(ERR);
    end
end
S = linspace(0,endval,divs)
n = 1:10
mesh(S,n,F)




%% functions

function Y = padeproduct(n,m,A,x,scl)
%% Computes the Pade n,m approximation of [e^A]x
%scl is a scale for the thresholder
a = -scl*max(max(abs(diag(diag(A)))));
ThrA = A - a*eye(size(A));

P = numerator(n,m,ThrA);
Q = denominator(n,m,ThrA);
p = numerator(n,m,a);
q = denominator(n,m,a);
% (q\p)*(Q\P)
Y = (p/q)*(Q\(P*x));
end

function [a,b] = nsplit(k)
%% a+b = k
k_iseven = mod(k+1,2);
a = floor(k/2) + not(k_iseven);
b = floor(k/2);
end

function MAT = numerator(n,m,A)
%% Pade Numerator matrix from formula
MAT = zeros(size(A));
for j = 0:n
    num = factorial(n + m - j)*factorial(n);
    den = factorial(n+m)*factorial(j)*factorial(n-j);
    MAT = MAT + (num/den)*mpower(A,j);
end
end

function MAT = denominator(n,m,A)
%% Pade denominator matrix
MAT = zeros(size(A));
for j = 0:m
    num = factorial(n + m - j)*factorial(m);
    den = factorial(n+m)*factorial(j)*factorial(m-j);
    MAT = MAT + (num/den)*mpower((-A),j);
end
end

function y = superPade(A,x)
%% 1,1 Positivity preserving Pade approximation using magic
a = -max(max(abs(diag(diag(A)))));
[d,~] = size(A);
ThrA = A - a*eye(size(A));

% abs(a) is the 2-norm of thresholded A, which is useful for getting
% stability conditions
% get the power for "squaring"
m = 0;
r = 10;
while r >= 1
    m = m + 1;
    r = sqrt(d)*abs(a)/(2^m);
end

P = numerator(1,1,ThrA./(2^m));
Q = denominator(1,1,ThrA./(2^m));

% also approximate the scalar
p = numerator(1,1,a/(2^m));
q = denominator(1,1,a/(2^m));


z = p/q;
Z = Q\P;
for k = 1:m
    Z = Z*Z;
    z = z*z;
end

y = z*Z*x;
end


function M = T(n,A)
%% truncated exponential
M = eye(size(A));
for k = n:-1:1
    M = M * A./k;
    M = M + eye(size(A));
end
end

function MAT = positiveSeriesExp(n,A,s)
%% Evaluates an approximation of the matrix exponential
% get positive A and diagonal offset
% series approximation
a = -s*max(max(abs(diag(diag(A)))));
ThrA = A - a*eye(size(A));

MAT = T(n,a)*T(n,ThrA);

end

function A = GL(N)
%% generates a graph-laplacian matrix of dimension N
A = rand(N);
A = A - diag(diag(A));
for k = 1:N
    A(k,k) = -sum(A(:,k));
end
end
