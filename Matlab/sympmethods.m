%% Symplectic Euler-VT Solver for arbitrary dimensional problem

% xinit = [0;
%     0;
%     10*cos(2*pi/3 + pi/4);
%     10*sin(2*pi/3 + pi/4);
%     10*cos(4*pi/3 + pi/4);
%     10*sin(4*pi/3 + pi/4);
%     0; 
%     0;
%     0;
%     -100;
%     100;
%     0
% ];


% R1 [ .746156, 0] ; R2 [-0.373078, .238313]; R3 [-0.373078, - .238313] (8)
% V1[ 0, .324677] ; V2 [.764226 , -.162339]; V3 [-[.764226 , -.162339

% xinit = [0.746156;
%     0;
%     -0.373078;
%     0.238313;
%     -0.373078;
%     -0.238313;
%     0;
%     0.324677;
%     0.764226;
%     -0.162339;
%     -0.764226;
%     -0.162339;
% ];

xinit = [1;
    0;
    0;
    0;
    -1;
    0;
    0.3471128135672417;
    0.532726851767674;
    0;
    0;
    -0.3471128135672417;
    -0.532726851767674
];


d = 6;
J = [
    zeros(6), eye(6);
    -eye(6), zeros(6);
]

% xinit = [3.8, 1.8, 0, 0]';
% X0 = [pi/4, 0]';
TSPAN = [0 20];
% opts = odeset('RelTol',1e-12)
% [T,X] = ode45(@(t,x) J*threeBody(t,x), TSPAN, xinit);
[T,X] = Verlet(@(t,x) threeBody(t,x), TSPAN, xinit);

% [T1,X1] = forwardEulerMethod(@(x) pendulum(x), TSPAN, X0);
% [T2,X2] = implicitMidpointMethod(@(x) pendulum(x), TSPAN, X0);
% [T3,X3] = backwardEulerMethod(@(x) pendulum(x), TSPAN, X0);

hold on

% plot(X1(:,1), X1(:,2))
% plot(X2(:,1), X2(:,2))
% plot(X3(:,1), X3(:,2))

plot(X(:,1),X(:,2))
plot(X(:,3),X(:,4),'x')
plot(X(:,5),X(:,6)) 

hold off


%% Functions

function dH = oscillator(x)
dH = [x(1); x(2)];
end

function dH = pendulum(x)
k=-1.2;
dH = [
    (k^2)*sin(x(1));
    x(2)
];
end

function dH = phonator(~, x)
lamda = 0.8;
beta = 3;
omega = 0.3;

q1 = x(1);
q2 = x(2);
p1 = x(3);
p2 = x(4);

f = @(x) 1-x + beta*(1 - 1./x.^2);
g = @(x) lamda*(1-x) + beta*(1 - 1./x.^2);

dH = [
    omega*(q1 - q2) - f(q1);
    omega*(q2 - q1) - g(q2);
    p1;
    p2;
];
end

function Y = phonatorHamilEval(X)
q1 = X(:,1);
q2 = X(:,2);
p1 = X(:,3);
p2 = X(:,4);

lamda = 0.8;
beta = 3;
omega = 0.3;

F = @(x) x - 0.5*x.^2 + beta*(x + 1./x);
G = @(y) lamda*(y - 0.5*y.^2) + beta*(y + 1./y);

Y = 0.5*(p1.^2 + p2.^2) + 0.5*omega*(q1-q2).^2 - F(q1) - G(q2);
end

function dH = threeBody(~,x)
%% Hamiltonian for the three body problem
m1 = 1/3; m2 = m1; m3 = m1;
G = 9.8;

q1 = [x(1); x(2)];
q2 = [x(3); x(4)];
q3 = [x(5); x(6)];
p1 = [x(7); x(8)];
p2 = [x(9); x(10)];
p3 = [x(11); x(12)];

r12 = q1-q2;
d12 = norm(r12);
r23 = q2-q3;
d23 = norm(r23);
r31 = q3-q1;
d31 = norm(r31);

dHdq1 = (G*m2*m1)*(-r12)./(d12^3) + (G*m3*m1)*(r31)./(d31^3);
dHdq2 = (G*m1*m2)*(r12)./(d12^3) + (G*m3*m2)*(-r23)./(d23^3);
dHdq3 = (G*m1*m3)*(-r31)./(d31^3) + (G*m2*m3)*(r23)/(d23^3);
dHdp1 = p1./m1;
dHdp2 = p2./m2;
dHdp3 = p3./m3;

dH = [
    -dHdq1;
    -dHdq2;
    -dHdq3;
    dHdp1;
    dHdp2;
    dHdp3
];

end


%3 body problem
%y3 project

%asymptotic convergence (poincare)

%conjugacy (trapezium and impl midpoint)

% make explicit, then also apply different methods (implicit)
function [TOUT, XOUT] = SymplecticEulerMethod(ODEHAMIL, TSPAN, X0)
%% Symplectic Euler Method (First Order)
%{
SymplecticEulerMethod(ODEHAMIL, TSPAN, X0) integrates the ODE defined by

    X' = J.ODEHAMIL(X)

ODEHAMIL is the Hamiltonian for the ODE. J is the block matrix [0 I; -I 0].
Applies the Symplectic Euler method.
X0 must be a column vector corresponding to initial conditions for the ODE.
TSPAN is the period of integration.
%}
h = 2^(-6);
[m,n] = size(X0);
dim = max(m,n)/2
TOUT = TSPAN(1):h:TSPAN(end);
J = [zeros(dim), eye(dim); -eye(dim), zeros(dim)]
X_init = X0;
XOUT = zeros(size(TOUT, 2), size(X0, 1));
XOUT(1,:) = X0;

for k = 2:size(TOUT, 2)

%perform an estimate using regular Euler.
X_guess = X_init + h*J*ODEHAMIL(TOUT(k), X_init);

%use fsolve to complete the method, applying the regular euler guess as a
%starting point
X = fsolve(@(x) x - X_init - h*J*ODEHAMIL(TOUT(k), [X_init(1:dim,:) ; x(dim+1:2*dim,:)]),X_guess,optimoptions("fsolve","Display","none"));
XOUT(k,:) = X;
X_init = X;
end

TOUT = TOUT';

end

function p_out = Kinetic(ODEHAMIL,t,p)
    H = ODEHAMIL(t, [zeros(6,1); p]);
    p_out = H(7:12);
end

function q_out = Potential(ODEHAMIL,t,q)
    H = ODEHAMIL(t, [q, zeros(6,1)]);
    q_out = H(1:6);
end
function [TOUT, XOUT] = Verlet(ODEHAMIL, TSPAN, X0)
%% Stormer-Verlet Integrator
dim = size(X0)

h = 2^(-10);
[m,n] = size(X0);
dim = max(m,n)/2
TOUT = TSPAN(1):h:TSPAN(end);
J = [zeros(dim), eye(dim); -eye(dim), zeros(dim)]
X_init = X0;
XOUT = zeros(size(TOUT, 2), size(X0, 1));
XOUT(1,:) = X0;

q = X0(1:dim)
p = X0(dim+1:2*dim)

%this method is explicit via. computing the separated hamiltonian
% very poorly written
for k = 2:size(TOUT, 2)
    t = TOUT(k);

    p_half = p - (h/2)*Potential(ODEHAMIL, t, q);
    q_next = q + h*Kinetic(ODEHAMIL, t, p_half);
    p_next = p_half - (h/2)*Potential(ODEHAMIL, t, q_next);
    
    XOUT(k,:) = [q_next' p_next'];

    p = p_next;
    q = q_next;
end

TOUT = TOUT'
end

function [TOUT, YOUT] = forwardEulerMethod(ODEHAMIL, TSPAN, Y0)
%% The classic
h = 0.001;
YOUT = [];
TOUT = TSPAN(1):h:TSPAN(2)
y = Y0;
dim = size(Y0,1)/2;
J = [zeros(dim), eye(dim); -eye(dim), zeros(dim)]

for t = TOUT
    y = y + h*J*ODEHAMIL(y);
    YOUT = [YOUT y]; 
end

TOUT = transpose(TOUT);
YOUT = transpose(YOUT);
end

function [TOUT, YOUT] = implicitMidpointMethod(ODEHAMIL, TSPAN, Y0)
%% Symplectic method!
h = 0.1;
YOUT = [];
TOUT = TSPAN(1):h:TSPAN(2)
y = Y0;
dim = size(Y0,1)/2;
J = [zeros(dim), eye(dim); -eye(dim), zeros(dim)]

for t = TOUT
    y_init = y + h*J*ODEHAMIL(y);
    targetfunc = @(x) x - y - h*J*ODEHAMIL((x+y)./2);
    y_iter = fsolve(targetfunc,y_init);
    y = y + h*J*ODEHAMIL((y_iter + y)./2);
    YOUT = [YOUT y];
end

TOUT = transpose(TOUT);
YOUT = transpose(YOUT);
end


function [TOUT, YOUT] = backwardEulerMethod(ODEHAMIL, TSPAN, Y0)
%% Feels like we only go backwards
h = 0.1;
YOUT = [];
TOUT = TSPAN(1):h:TSPAN(2)
y = Y0;
dim = size(Y0,1)/2;
J = [zeros(dim), eye(dim); -eye(dim), zeros(dim)]

for t = TOUT
    y_init = y + h*J*ODEHAMIL(y);
    targetfunc = @(x) x - y - h*J*ODEHAMIL(x);
    y_iter = fsolve(targetfunc,y_init);
    y = y + h*J*ODEHAMIL(y_iter);
    YOUT = [YOUT y];
end

TOUT = transpose(TOUT);
YOUT = transpose(YOUT);
end