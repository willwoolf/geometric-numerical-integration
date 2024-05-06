%% magnus integrator

[T,X] = MagnusEM2(GL,)


function [TOUT, XOUT] = MagnusEM2(ODEMATR, TSPAN, X0, h)
TOUT = TSPAN(1):h:TSPAN(2)
x = X0;
dim = size(X0,1);
XOUT = zeros(size(TOUT, 2),dim);
XOUT(1,:) = x';

for k = 2:size(TOUT,2)
    % two matrix exponentials
    U = (h/2)*A(x);
    V = mypade(2,2,U)*x;
    W = h*A(V);
    x_n = mypade(2,2,W)*x;
    
    XOUT(k,:) = x_n';

    x = x_n;
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

function f = Drv(~,y)
g = A(y);
f = g*y;
end