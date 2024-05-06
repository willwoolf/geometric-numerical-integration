%% Integration of simple pendulum
[T,Y] = eulerSolver(@(q,p)pendulumGradH(q,p), [0 10], [5 0])
%solves the ODE using Euler's method

subplot(121)
plot(Y(:,1),Y(:,2));
subplot(122)
plot(T,Y(:,1));

function gradH = pendulumGradH(q,p)
m = 1;
l = 5;
g = 9.8;
k = sqrt(g/l);

dhdq = (k^2)*q;
dhdp = p;

gradH = [
    dhdq;
    dhdp
];

end

function [tout, yout] = eulerSolver(ODEFUNC,TSPAN,Y0)
h = 0.01;
%note that Y0 should have two columns for a hamiltonian problem in 2*d
y = Y0;
yout = []
tout = (TSPAN(1)):h:TSPAN(2);

J = [0 1; -1 0];

for t = tout
    % euler operation and append to yout
    y(1)
    y(2)
    J*ODEFUNC(y(1),y(2))
    ynew = y + transpose(h.*J*ODEFUNC(y(1), y(2)))
    yout = [yout; ynew];
    y = ynew;
end

tout = transpose(tout);
end

function[tout, yout] = backwardEulerSolver(ODEFUNC,TSPAN,Y0)
h = 0.01;
y = Y0;




end



