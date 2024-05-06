%% exponential splitting method on chemical kinetics problem

% [T,Z,X] = threeExponentialMethod(@A, [0 0.3], [1 0 0]');
% [T2, Z2, X2] = twoSubsOneExp(@A, [0 0.3], [1 0 0]');

%gblerr(0.3);
tspan = [12*3600, 72*3600];
init = [9.906e1, 6.624e8, 5.326e11, 1.697e16, 8.725e8, 2.240e8]';
% [TM,XM] = MagnusEM2(@(t,y) B(t,y), tspan, init, 30);
[TS,XS] = threeExponentialMethod(@(t,y) B(t,y), tspan, init, 30);



% [TE, XE] = forwardEulerMethod(@(t,y) B(t,y)*y, tspan, init, 0.01);
% [TC, XC] = clippingMethod(@(t,y) B(t,y)*y, tspan, init, 0.01);
% 
% [TRK, XRK] = ode45(@(t,y) B(t,y)*y, tspan, init);

semilogy(TS,XS)
xlim(tspan)
ylim([10, 1e16])

%% generate the invariants

MC = [1,1,3,2,1,2]*(XC');
ME = [1,1,3,2,1,2]*(XE');
MRK = [1,1,3,2,1,2]*(XRK');

M2C = [0 0 0 0 1 1]*(XC');
M2E = [0 0 0 0 1 1]*(XE');
M2RK = [0 0 0 0 1 1]*(XRK');

%% Linear Test Problem

tspan = [0 4];
xinit = [1];
lambda = -2.5

ltp = @(t,x) lambda*x;

[TE, XE] = forwardEulerMethod(ltp, tspan, xinit, 0.5)
[TI, XI] = backwardEulerMethod(ltp, tspan, xinit, 0.5)
[TC, XC] = clippingMethod(ltp, tspan, xinit, 0.5)

ts = linspace(0,4,50)
xs = exp(lambda.*ts)


%% plotting

semilogy(TE,XE)
title("Euler simulation of a stratospheric reaction model")
xlabel("time t")
ylabel("concentration")
legend("O^{1D}","O","O_3","O_2","NO","NO_2")
ylim([1e-4,1e16])
xlim(tspan)


%% invariants

semilogy(TE,M2E)
title("Euler")
xlabel("time t")
ylabel("invariant mass")
xlim(tspan)



%% ode functions


function [TOUT, YOUT] = forwardEulerMethod(ODEFUNC, TSPAN, Y0, h)
YOUT = [];
TOUT = TSPAN(1):h:TSPAN(2);
y = Y0;

for t = TOUT
    YOUT = [YOUT y]; 
    y = y + h*ODEFUNC(t,y);
end

TOUT = transpose(TOUT);
YOUT = transpose(YOUT);
end


function [TOUT, YOUT] = clippingMethod(ODEFUNC, TSPAN, Y0,h)
YOUT = [];
TOUT = TSPAN(1):h:TSPAN(2);
y = Y0;

for t = TOUT
    YOUT = [YOUT y]; 
    y = y + h*ODEFUNC(t,y);
    y(y<0) = 0;
end

TOUT = transpose(TOUT);
YOUT = transpose(YOUT);
end

function [TOUT, YOUT] = backwardEulerMethod(ODEFUNC, TSPAN, Y0, h)
YOUT = [];
TOUT = TSPAN(1):h:TSPAN(2)
y = Y0;

for t = TOUT
    YOUT = [YOUT y];
    y_init = y + h*ODEFUNC(t,y);
    targetfunc = @(x) x - y - h*ODEFUNC(t+h,x);
    y_iter = fsolve(targetfunc,y_init);
    y = y_iter;
end

TOUT = transpose(TOUT);
YOUT = transpose(YOUT);
end

function [TOUT, YOUT] = exponentialEulerMethod(ODEFUNC, TSPAN, Y0, h)
YOUT = [];
TOUT = TSPAN(1):h:TSPAN(2);
y = Y0;

for t = TOUT
    YOUT = [YOUT y]; 
    y = expm(h*A(t,y))*y
end

TOUT = transpose(TOUT);
YOUT = transpose(YOUT);
end

function [TOUT, ZOUT, XOUT] = threeExponentialMethod(ODEMATR, TSPAN, X0, h)
TOUT = TSPAN(1):h:TSPAN(2);
x = X0;
z = X0;
dim = size(X0,1);
XOUT = zeros(size(TOUT, 2),dim);
ZOUT = XOUT;
XOUT(1,:) = x';
ZOUT(1,:) = z';

for k = 2:size(TOUT,2)
    t = TOUT(k);

    x_hf = expm((h/2)*ODEMATR(t,z))*x;

    z_n = expm(h*ODEMATR(t,x_hf))*z;
    ZOUT(k,:) = z_n'; 

    x_n = expm((h/2)*ODEMATR(t,z_n))*x_hf;
    XOUT(k,:) = x_n';

    x = x_n;
    z = z_n;
end

TOUT = transpose(TOUT);
end

function [TOUT, XOUT] = MagnusEM2(ODEMATR, TSPAN, X0, h)
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


function [TOUT, ZOUT, XOUT] = approximatedExponentials(ODEMATR, TSPAN, X0, h)
TOUT = TSPAN(1):h:TSPAN(2)
x = X0;
z = X0;
dim = size(X0,1);
XOUT = zeros(size(TOUT, 2),dim);
ZOUT = XOUT;
XOUT(1,:) = x';
ZOUT(1,:) = z';

for k = 2:size(TOUT,2)
    x_hf = (eye(3) - (h/2)*ODEMATR(z))\x;

    z_n = (eye(3) - h*ODEMATR(x_hf))\z;
    ZOUT(k,:) = z_n'; 

    x_n = (eye(3) - (h/2)*ODEMATR(z_n))\x_hf;
    XOUT(k,:) = x_n';

    x = x_n;
    z = z_n;
end

TOUT = transpose(TOUT);
end

function [TOUT, ZOUT, XOUT] = twoSubsOneExp(ODEMATR, TSPAN, X0, h)
TOUT = TSPAN(1):h:TSPAN(2)
x = X0;
z = X0;
dim = size(X0,1);
XOUT = zeros(size(TOUT, 2),dim);
ZOUT = XOUT;
XOUT(1,:) = x';
ZOUT(1,:) = z';

for k = 2:size(TOUT,2)
    % approximation
    x_hf = (eye(3) - (h/2)*ODEMATR(z))\x;

    % approximation
    z_n = (eye(3) - h*ODEMATR(x_hf))\z;
    ZOUT(k,:) = z_n'; 
    
    % matrix exponential
    x_n = expm((h/2)*ODEMATR(z_n))*x_hf;
    XOUT(k,:) = x_n';

    x = x_n;
    z = z_n;
end

TOUT = transpose(TOUT);
end

function [TOUT, ZOUT, XOUT] = oneSub(ODEMATR, TSPAN, X0, h)
TOUT = TSPAN(1):h:TSPAN(2)
x = X0;
z = X0;
dim = size(X0,1);
XOUT = zeros(size(TOUT, 2),dim);
ZOUT = XOUT;
XOUT(1,:) = x';
ZOUT(1,:) = z';

for k = 2:size(TOUT,2)
    % approximation
    x_hf = (eye(3) - (h/2)*ODEMATR(z))\x;

    % matrix exponential
    z_n = expm(h*ODEMATR(x_hf))*z;
    ZOUT(k,:) = z_n'; 
    
    % matrix exponential
    x_n = expm((h/2)*ODEMATR(z_n))*x_hf;
    XOUT(k,:) = x_n';

    x = x_n;
    z = z_n;
end

TOUT = transpose(TOUT);
end

function GL = A(y)
GL = [
    -0.04, y(3)*1e4, 0;
    0.04, y(2)*(-3e7), 0;
    0, y(2)*(3e7), 0
];
end

function v = sig(t, TR, TS)
TL = mod((t/3600),24);
if (TR <= TL) & (TL <= TS)
    k = (2*TL - TR - TS)/(TS - TR);
    v = 0.5*(1 + cos(pi*abs(k)*k));
else
    v = 0;
end

end

function STR = B(t,y)
TR = 4.5;
TS = 19.5;
k1 = 2.643*(1e-10)*sig(t,TR,TS)^3;
k2 = 8.018*1e-17;
k3 = 6.12*(1e-4)*sig(t,TR,TS);
k4 = 1.576*1e-15;
k5 = 1.07*(1e-3)*sig(t,TR,TS)^2;
k6 = 7.11*1e-11;
k7 = 1.2*10^-10;
k8 = 6.062*1e-15;
k9 = 1.069*1e-11;
k10 = 1.289*(1e-2)*sig(t,TR,TS);

gamma = k3 + k5 + k4*y(2) + k7*y(1) + k8*y(5);
STR = [
    -(k6 + k7*y(3)), 0,                              k5,                  0,               0,         0;
    k6,              -(k2*y(4) + k4*y(3) + k9*y(6)), k3,                  2*k1,            0,         k10;
    0,               k2*y(4)/3,                      -gamma,              2*k2*y(2)/3,     0,         0;
    k7*y(3)/2,       k4*y(3)+k9*y(6)/2,              (gamma + k7*y(1)/2), -(k1 + k2*y(2)), 0,         k9*y(2)/2;
    0,               0,                              0,                   0,               -k8*y(3),  k10+k9*y(2);        
    0,               0,                              0,                   0,               k8*y(3),   -(k10+k9*y(2))
];
end

function f = Drv(~,y)
g = A(y);
f = g*y;
end

function gblerr(T)
tspan = [0, T];
yinit = [1 0 0]';
DIVS = 18;
D = 1:DIVS;
opts = odeset('RelTol',1e-10);
[TM, XM] = ode45(@(t,y) Drv(t,y), tspan, yinit, opts);
exact = XM(end,:)';

H = 0.3./(2.^(D-1));

ERR = zeros(DIVS,3)

for d = D
    h = 0.3/(2^(d-1));
    [T1, ~, X1] = threeExponentialMethod(@A, tspan, yinit, h);
    [T2, ~, X2] = oneSub(@A, tspan, yinit, h);
    [T3, ~, X3] = twoSubsOneExp(@A, tspan, yinit, h);
    eee = X1(end,:)';
    aee = X2(end,:)';
    aae = X3(end,:)';
    ERR(d, 1) = norm(exact - eee, 2);
    ERR(d, 2) = norm(exact - eae, 2);
    ERR(d, 3) = norm(exact - aae, 2);
end

loglog(H,ERR)
end