%% long term integration (log timestep) of Robertson

h = 10^-6;
t = 0;
X0 = [1 0 eps]'
XOUT = X0'
ZOUT = X0'
x = X0;
z = X0;
T = [0];
while t < 10^5
    x_hf = expm((h/2)*A(z))*x;

    z_n = expm(h*A(x_hf))*z;
    ZOUT = [ZOUT; z_n']; 

    x_n = expm((h/2)*A(z_n))*x_hf;
    XOUT = [XOUT; x_n'];

    x = x_n;
    z = z_n;
    t = t+h;
    T = [T t];
    h = 1.8*h;
end


function GL = A(y)
GL = [
    -0.04, y(3)*1e4, 0;
    0.04, y(2)*(-3e7), 0;
    0, y(2)*(3e7), 0
];
end