%% ODE solvers: forward and backward Euler methods

TSPAN = [0 5];
Y0 = [1];

% [T1,X1] = forwardEulerMethod(@(x)lineartest(x),TSPAN,Y0)
% [T2,X2] = clippingMethod(@(x)lineartest(x),TSPAN,Y0)
% [T2,X2] = implicitMidpointMethod(@(x)pendulum(x),TSPAN,Y0)
% [T3,X3] = backwardEulerMethod(@(x)pendulum(x),TSPAN,Y0)



hold on
% subplot(121)
% plot(X1(:,1),X1(:,2),X2(:,1),X2(:,2));
% subplot(122)
% plot(T1,X1(:,1),T2,X2(:,1));

for N = 2:2:20
    h = 5/N;
    [T,X] = exponentialEulerMethod(@(x) lineartest(x), TSPAN, Y0, h)
    plot(T,X)
end

%plot(T1, X1);
% plot(T2, X2);
% 
% X = linspace(0,5,50);
% Y = exp(-2.5*X);
% 
% plot(X, Y);

% plot(X1(:,1), X1(:,2))
% plot(X2(:,1), X2(:,2))
% plot(X3(:,1), X3(:,2))

function dxdt = oscillator(x)
dxdt = [x(2); -x(1)];
end

function dxdt = lineartest(x)
lambda = -2.5;
dxdt = lambda;
end

%linear test problem is unstable if |xi| < |xi - h l xi|
% |1/(1-hl)| < 1

function dxdt = pendulum(x)
k = 2;
dxdt = [(k^2)*sin(x(1)); x(2)];
end

function [TOUT, YOUT] = forwardEulerMethod(ODEFUNC, TSPAN, Y0)
h = 0.5;
YOUT = [];
TOUT = TSPAN(1):h:TSPAN(2)
y = Y0;

for t = TOUT
    YOUT = [YOUT y]; 
    y = y + h*ODEFUNC(y);
end

TOUT = transpose(TOUT);
YOUT = transpose(YOUT);
end

function [TOUT, YOUT] = clippingMethod(ODEFUNC, TSPAN, Y0)
h = 0.5;
YOUT = [];
TOUT = TSPAN(1):h:TSPAN(2)
y = Y0;

for t = TOUT
    YOUT = [YOUT y]; 
    y = y + h*ODEFUNC(y);
    for i = 1:size(y)
        if y(i) < 0
            y(i) = 0;
        end
    end
end

TOUT = transpose(TOUT);
YOUT = transpose(YOUT);
end

function [TOUT, YOUT] = exponentialEulerMethod(ODEFUNC, TSPAN, Y0, h)
YOUT = [Y0];
TOUT = TSPAN(1):h:TSPAN(2)
y = Y0;

for t = TOUT(2:end)
    y = expm(h*ODEFUNC(y))*y;
    YOUT = [YOUT y]; 
end

TOUT = transpose(TOUT);
YOUT = transpose(YOUT);
end


function [TOUT, YOUT] = backwardEulerMethod(ODEFUNC, TSPAN, Y0)
h = 0.1;
YOUT = [];
TOUT = TSPAN(1):h:TSPAN(2)
y = Y0;

for t = TOUT
    y_init = y + h*ODEFUNC(y);
    targetfunc = @(x) x - y - h*ODEFUNC(x);
    y_iter = fsolve(targetfunc,y_init);
    y = y + h*ODEFUNC(y_iter);
    YOUT = [YOUT y];
end

TOUT = transpose(TOUT);
YOUT = transpose(YOUT);
end


function [TOUT, YOUT] = implicitMidpointMethod(ODEFUNC, TSPAN, Y0)
h = 0.1;
YOUT = [];
TOUT = TSPAN(1):h:TSPAN(2)
y = Y0;

for t = TOUT
    y_init = y + h*ODEFUNC(y);
    targetfunc = @(x) x - y - h*ODEFUNC((x+y)./2);
    y_iter = fsolve(targetfunc,y_init);
    y = y + h*ODEFUNC((y_iter + y)./2);
    YOUT = [YOUT y];
end

TOUT = transpose(TOUT);
YOUT = transpose(YOUT);
end


function [TOUT, YOUT] = trapeziumMethod(ODEFUNC, TSPAN, Y0)
h = 0.01;
YOUT = [];
TOUT = TSPAN(1):h:TSPAN(2)
y = Y0;

for t = TOUT
    y_init = y + h*ODEFUNC(y);
    targetfunc = @(x) x - y - h*(ODEFUNC(x) + ODEFUNC(y))./2;
    y_iter = fsolve(targetfunc,y_init);
    y = y + h*(ODEFUNC(y_iter) + ODEFUNC(y))/2;
    YOUT = [YOUT y];
end

TOUT = transpose(TOUT);
YOUT = transpose(YOUT);
end