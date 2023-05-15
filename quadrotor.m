function x_next = quadrotor(x,u,dt)

x_next = zeros(12,1);
g = 9.81;
m = 0.650;
l = 0.23;
Ixx = 7.5*10^(-3);
Iyy = Ixx;
Izz = 1.3*10^(-2);

a1 = (Iyy - Izz)/Ixx;
a2 = (Izz - Ixx)/Iyy;
a3 = (Ixx - Iyy)/Izz;

b1 = l/Ixx;
b2 = l/Iyy;
b3 = l/Izz;
ux = cos(x(1))*sin(x(3))*cos(x(5)) + sin(x(1))*sin(x(5));
uy = cos(x(1))*sin(x(3))*sin(x(5)) - sin(x(1))*cos(x(5));

x_next(1) = x(1) + dt*x(2);
x_next(2) = x(2) + dt*(x(4)*x(6)*a1 + b1*u(2));
x_next(3) = x(3) + dt*x(4);
x_next(4) = x(4) + dt*(x(2)*x(6)*a2 + b2*u(3));
x_next(5) = x(5) + dt*x(6);
x_next(6) = x(6) + dt*(x(4)*x(6)*a3 + b3*u(4));
x_next(7) = x(7) + dt*x(8);
x_next(8) = x(8) + dt*((cos(x(1))*cos(x(3))*u(1))/m - g);
x_next(9) = x(9) + dt*x(10);
x_next(10) = x(10) + dt*(ux/m*u(1));
x_next(11) = x(11) + dt*x(12);
x_next(12) = x(12) + dt*(uy/m*u(1));

end

