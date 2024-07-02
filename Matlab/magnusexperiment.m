%% magnus integrator


init = [0.1, 0.175, 0.15, 1.15, 0.81, 0.5]'
tspan = [0 200];
% opts = odeset(RelTol=1e-12); %12 is as good as it gets
% [T_master,X_master] = ode45(@(t,y) MAPK(t,y)*y, tspan, init, opts);

[T_em2, X_em2] = MagnusEM2(@MAPK, tspan, init, 0.1);
[T_ip2, X_ip2] = MagnusIP2(@MAPK, tspan, init, 0.1);
[T_eb2, X_eb2] = MagnusEB2(@MAPK, tspan, init, 0.0005);

subplot(131)
% plot(T_master, X_master)
plot(T_em2, X_em2)

subplot(132)
plot(T_ip2, X_ip2)

subplot(133)
plot(T_eb2, X_eb2)



%% interpolate data

X_EM2 = interp1(T_em2, X_em2, T_master);
X_IP2 = interp1(T_ip2, X_ip2, T_master);
X_EB2 = interp1(T_eb2, X_eb2, T_master);

% subplot(141)
% plot(T_master, X_master)
% 
% subplot(142)
% plot(T_master, X_EM2)
% 
% subplot(143)
% plot(T_master, X_IP2)
% 
% subplot(144)
% plot(T_master, X_EB2)

%interp1(old time points, data points, new time points) returns resampled
%data

%% error

subplot(131)
plot(T_master, vecnorm((X_EM2-X_master),2,2))

subplot(132)
plot(T_master, vecnorm((X_IP2-X_master),2,2))

subplot(133)
plot(T_master, vecnorm((X_EB2-X_master),2,2))

%% error data

% for each value of the timestep
% % compute the solution of the chosen method at the end time
% % take the norm of the error at the end
% plot global error against timestep
PN = 11;
H = zeros(PN,1);
for i = 1:PN
    H(i) = 1/(2^(i+2));
end

H_special = H./16;

X_true = X_master(end,:)

%begin
%% EM2

E = zeros(size(H)) % vector of global errors
for k = 1:PN
    h=H(k);
    [~, X_emsol] = MagnusEM2(@MAPK, tspan, init, h);
    x_emfin = X_emsol(end,:);
    E(k) = norm(x_emfin' - X_true',2)/norm(X_true',2);
end
loglog(H,E)

%% IP2
% same for IP2
E2 = zeros(size(H))
for k = 1:PN
    h=H(k);
    [~, X_ipsol] = MagnusIP2(@MAPK, tspan, init, h);
    x_ipfin = X_ipsol(end,:);
    E2(k) = norm(x_ipfin' - X_true',2)/norm(X_true',2);
end
loglog(H,E2)

%% EB2
% same for EB2
E3 = zeros(size(H)) % vector of global errors
for k = 1:PN
    h=H(k);
    [~, X_ebsol] = MagnusEB2(@MAPK, tspan, init, h);
    x_ebfin = X_ebsol(end,:);
    E3(k) = norm(x_ebfin' - X_true',2)/norm(X_true',2);
end
loglog(H,E3)

%% IQ2
% same for IQ2
E3 = zeros(size(H)) % vector of global errors
for k = 1:PN
    h=H(k);
    [~, X_iqsol] = MagnusIQ2(@MAPK, tspan, init, h);
    x_iqfin = X_iqsol(end,:);
    E3(k) = norm(x_iqfin' - X_true',2)/norm(X_true',2);
end
loglog(H,E3)

%% functions

function [TOUT, XOUT] = MagnusEM2(ODEMATR, TSPAN, X0, h)
%% Second order Magnus integrator
TOUT = TSPAN(1):h:TSPAN(2);
x = X0;
dim = size(X0,1);
XOUT = zeros(size(TOUT, 2),dim);
XOUT(1,:) = x';

for k = 2:size(TOUT,2)
    t = TOUT(k);
    % two matrix exponentials
    U = (h/2)*ODEMATR(t,x);
    V = expm(U)*x;
    W = h*ODEMATR(t,V);
    x_n = expm(W)*x;
    
    XOUT(k,:) = x_n';

    x = x_n;
end

TOUT = transpose(TOUT);
end

function [TOUT, XOUT] = MagnusIP2(ODEMATR, TSPAN, X0, h)
%% Magnus EM2 using inverse and pade approximations
TOUT = TSPAN(1):h:TSPAN(2);
x = X0;
dim = size(X0,1);
XOUT = zeros(size(TOUT, 2),dim);
XOUT(1,:) = x';

for k = 2:size(TOUT,2)
    t = TOUT(k);
    % inverse (solve) and Pade
    U = (h/2)*ODEMATR(t,x);
    V = (eye(size(U)) - U)\x;
    W = h*ODEMATR(t,V);
    x_n = padeproduct(1,1,W,x,1);
    
    XOUT(k,:) = x_n';

    x = x_n;
end

TOUT = transpose(TOUT);
end


function [TOUT, XOUT] = MagnusEB2(ODEMATR, TSPAN, X0, h)
%% Magnus EM2 using series approximations
TOUT = TSPAN(1):h:TSPAN(2);
x = X0;
dim = size(X0,1);
XOUT = zeros(size(TOUT, 2),dim);
XOUT(1,:) = x';

for k = 2:size(TOUT,2)
    t = TOUT(k);
    % series approximations of exponential, first and second order, with
    % diagonal offset 0
    s = 1;
    U = (h/2)*ODEMATR(t,x);
    V = positiveSeriesExp(1,U,s)*x;
    W = h*ODEMATR(t,V);
    x_n = positiveSeriesExp(2,W,s)*x;
    
    XOUT(k,:) = x_n';

    x = x_n;
end

TOUT = transpose(TOUT);
end

function [TOUT, XOUT] = MagnusIQ2(ODEMATR, TSPAN, X0, h)
%% Magnus EM2 using inverse and positive pade approximation
TOUT = TSPAN(1):h:TSPAN(2);
x = X0;
dim = size(X0,1);
XOUT = zeros(size(TOUT, 2),dim);
XOUT(1,:) = x';

for k = 2:size(TOUT,2)
    t = TOUT(k);
    % inverse (solve) and Pade
    U = (h/2)*ODEMATR(t,x);
    V = (eye(size(U)) - U)\x;
    W = h*ODEMATR(t,V);
    x_n = superPade(W,x);
    
    XOUT(k,:) = x_n';

    x = x_n;
end

TOUT = transpose(TOUT);
end


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


function MAT = positiveSeriesExp(n,A,s)
%% Evaluates an approximation of the matrix exponential
% get positive A and diagonal offset
% series approximation
a = -s*max(max(abs(diag(diag(A)))));
ThrA = A - a*eye(size(A));

MAT = T(n,a)*T(n,ThrA);
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
    M = M * (A./k);
    M = M + eye(size(A));
end
end

function GL = A(y)
GL = [
    -0.04, y(3)*1e4, 0;
    0.04, y(2)*(-3e7), 0;
    0, y(2)*(3e7), 0
];
end

function GL = MAPK(~,y)
a = 0.1;
k1 = 100/3;
k2 = 1/3;
k3 = 50;
k4 = 1/2;
k5 = 10/3;
k6 = 1/10;
k7 = 7/10;
GL = [
    -k7-k1*y(2) 0 0 k2 0 k6;
    0 -k1*y(1) k5 0 0 0;
    0 0 -k3*y(1)-k5 k2 k4 0;
    (1 - a)*k1*y(2) a*k1*y(1) 0 -k2 0 0;
    0 0 k3*y(1) 0 -k4 0;
    k7 0 0 0 0 -k6
];
end

function f = Drv(~,y)
g = A(y);
f = g*y;
end