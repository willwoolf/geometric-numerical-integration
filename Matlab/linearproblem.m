
TSPAN = [0 100];
INIT = [0 1]';
h = 0.1;

opts = odeset(RelTol=0.5);
[T,Z,X] = threeExponentialMethod(@(y) linear(y), TSPAN, INIT, h);
[T2, X2] = ode45(@(t,y) linear(y)*y, TSPAN, INIT);

subplot(121)
plot(T, X)

subplot(122)
plot(T2, X2)

function [TOUT, ZOUT, XOUT] = threeExponentialMethod(ODEMATR, TSPAN, X0, h)
TOUT = TSPAN(1):h:TSPAN(2);
x = X0;
z = X0;
dim = size(X0,1);
XOUT = zeros(size(TOUT, 2),dim);
ZOUT = XOUT;
XOUT(1,:) = x';
ZOUT(1,:) = z';

ODEMATR(z)

for k = 2:size(TOUT,2)
    x_hf = expm((h/2)*ODEMATR(z))*x;

    z_n = expm(h*ODEMATR(x_hf))*z;
    ZOUT(k,:) = z_n'; 

    x_n = expm((h/2)*ODEMATR(z_n))*x_hf;
    XOUT(k,:) = x_n';

    x = x_n;
    z = z_n;
end
end

function MAT = linear(~)
a = eps^(200);
MAT = [
    -a, 1;
    a, -1
];
end

%