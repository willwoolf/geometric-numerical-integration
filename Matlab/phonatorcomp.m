%% phonation comparison. integrate with ODE45 and with symplectic euler


alpha = 1;
lambda = 0.5;
beta = 3;
omega = 0.3;

TSPAN = [0, 50];
q_init = [1, 3.79]';
p_init = [0, 0]';
X0 = [q_init(1) p_init(1) q_init(2) p_init(2)/alpha];

[TS, QS, PS] = phonatorSymplectic(TSPAN, q_init, p_init, alpha, lambda, beta, omega);

[TR, XR] = ode45(@(t,x) phonator(t,x,alpha, lambda, beta, omega), TSPAN, X0);

F = @(x) x - (1/2)*x^2 + beta*(x + 1/x);
G = @(x) lambda*(x - (1/2)*x^2) + beta*(x + 1/x);

Ham = phonatorHamil(QS,PS,alpha, lambda, beta, omega);

hold on
plot(TS, Ham)


function [TOUT, QOUT, POUT] = phonatorSymplectic(TSPAN, q_init, p_init, alpha, lambda, beta, omega)
%{
    performs the Stormer-Verlet method integration of the phonation model
%}
h = 2^(-5);

f = @(x) 1-x + beta*(1 - 1./x.^2);
g = @(x) lambda*(1-x) + beta*(1 - 1./x.^2);

q_prev = q_init;
p_prev = p_init;

TOUT = TSPAN(1):h:TSPAN(end);
QOUT = zeros(size(TOUT, 2), size(q_init, 1));
QOUT(1,:) = q_init';
POUT = zeros(size(TOUT, 2), size(p_init, 1));
POUT(1,:) = p_init';

for k = 2:size(TOUT, 2)
    q_mid = q_prev + (h/2)*[p_prev(1); p_prev(2)/alpha];
    p_n = p_prev - h*[omega*(q_mid(1) - q_mid(2)) - f(q_mid(1)); omega*(q_mid(2)-q_mid(1)) - g(q_mid(2))];
    q_n = q_mid + (h/2)*[p_n(1); p_n(2)/alpha];
    QOUT(k,:) = q_n';
    POUT(k,:) = p_n';
    q_prev = q_n;
    p_prev = p_n;
end

TOUT = TOUT';

end


function dhdt = phonator(~,x,alpha,lambda,beta,omega)

%take inputs, initialise the h variables
dhdt = zeros(size(x));

h_1 = x(1);
dh_1 = x(2);
h_2 = x(3);
dh_2 = x(4);

%compute the equation

dhdt(1) = dh_1;
dhdt(2) = 1 - h_1 + beta*(1-1/(h_1)^2) + omega*(h_2-h_1);
dhdt(3) = dh_2;
dhdt(4) = (lambda*(1 - h_2) + beta*(1-1/(h_2)^2) + omega*(h_1-h_2));
dhdt(4) = dhdt(4)/alpha;

end

function H = phonatorHamil(QOUT, POUT, alpha, lambda, beta, omega)
F = @(x) x - (1/2)*x.^2 + beta*(x + 1./x);
G = @(x) lambda*(x - (1/2)*x.^2) + beta*(x + 1./x);

sz = size(QOUT,1);
H = zeros(sz, 1);
for k = 1:sz
    H(k) = (1/2)*(POUT(k,1)^2 + (POUT(k,2)^2)/alpha) + (omega/2)*(QOUT(k,1)-QOUT(k,2))^2 - F(QOUT(k,1)) - G(QOUT(k,2));
end
end