clc;clear all;close all

%% %% Linear distance from Velocity

load('lineardistancefromvelocity.mat')
load('quatxyzw.mat')
load('angularvelocityfused.mat')
load('controlinputs.mat')

whos 


%% Linear Velocity and Distance
linearvelocity = [slashmatrice210v2slashvelocity(:,1) slashmatrice210v2slashvelocity(:,2) slashmatrice210v2slashvelocity(:,3)];
x_dot = linearvelocity(:,1);    % Linear Velocity in X Direction
y_dot = linearvelocity(:,2);    % Linear Velocity in Y Direction
z_dot = linearvelocity(:,3);    % Linear Velocity in Z Direction

disp('Linear Velocity ')
disp([x_dot, y_dot, z_dot])

n1 = size(linearvelocity,1);         % Number of Measurement
freq_velocity = 50;                 % Frequency of velocity sensor
sample_velocity = 1/freq_velocity;  % Sample time of velocity sensor i.e 0.02sec
t1 = sample_velocity : sample_velocity : n1*sample_velocity;
disp('Duration of t1')
length(t1)

x = cumtrapz(t1,x_dot);  % Position X
y = cumtrapz(t1,y_dot);  % Position Y
z = cumtrapz(t1,z_dot);  % Position Z

disp('Distance')
disp([x, y, z])

%% Convert Quaternions data to Eular Angles
quat = [slashmatrice210v2slashattitude(:,4) slashmatrice210v2slashattitude(:,1) slashmatrice210v2slashattitude(:,2) slashmatrice210v2slashattitude(:,3)];
eulZYX = quat2eul(quat);

%Radians
phi = eulZYX(:,3);      % Roll (rad) (X Axis Rotation)
theta = eulZYX(:,2);    % Pitch (rad) (Y Axis Rotation)
psi = eulZYX(:,1);      % Yaw (rad) (Z Axis Rotation)

disp('Eular Angles (Radians)')
disp([phi, theta, psi])

%Degrees
D = rad2deg([phi theta psi]);
D_phi = D(:,1);
D_theta = D(:,2);
D_psi = D(:,3);

disp('Eular Angles (Degrees)')
disp([D_phi, D_theta, D_psi])

n2 = size(quat,1);           % Number of Measurement
freq_imu = 100;             % Frequency of fused sensor
sample_imu = 1/freq_imu;    % Sample time of fused sensor 
t2 = sample_imu : sample_imu : n2*sample_imu;
disp('Duration of t2')
length(t2)

%% Angular Velocity

AngularVelocity = [slashmatrice210v2slashangularvelocityfused(:,1)   slashmatrice210v2slashangularvelocityfused(:,2)   slashmatrice210v2slashangularvelocityfused(:,3)];
phi_dot = AngularVelocity(:,1);     % Angular Velocity Roll X Direction
theta_dot = AngularVelocity(:,2);   % Angular Velocity Pitch Y Direction 
psi_dot = AngularVelocity(:,3);     % angular Velocity Yaw Z Direction

disp('Angular Velocity')
disp([phi_dot, theta_dot, psi_dot])


n3 = size(AngularVelocity,1);    % Number of Measurement
freq_imu = 100;                 % Frequency of Fused sensor
sample_imu = 1/freq_imu;        % Sample time of Fused sensor 
t3 = sample_imu : sample_imu : n3*sample_imu;
disp('Duration of t3')
length(t3)

n4 = size(controlinputs,1);    % Number of Measurements
freq_rc = 50;                 % Frequency of rc sensor
sample_rc = 1/freq_rc;        % Sample time of rc sensor 
t4 = sample_rc : sample_rc : n4*sample_rc;
disp('Duration of t4')
length(t4)

%%  Plots

%Distance from Linear Velocity
figure()
subplot(3,1,1)
plot(t1,x);xlabel('Time (sec)');ylabel('X (m)');grid on
legend('X distance from velocity')
title('Linear Distance')
subplot(3,1,2)
plot(t1,y);xlabel('Time (sec)');ylabel('Y (m)');grid on
legend('Y distance from velocity')
subplot(3,1,3)
plot(t1,z);xlabel('Time (sec)');ylabel('Z (m)');grid on
legend('Z distance from velocity')

figure()
plot3(x,y,z,'r','linewidth',1.25)
title('Integrated Linear Distance from velocity')
xlabel('X(m)');ylabel('Y(m)');zlabel('Z(m)')
legend('Distance')

%Linear Velocity
figure()
subplot(3,1,1)
plot(t1,x_dot);xlabel('Time (sec)');ylabel('X_dot (m/s)');grid on
legend('X Velocity')
title('Linear Velocity')
subplot(3,1,2)
plot(t1,y_dot);xlabel('Time (sec)');ylabel('Y_dot (m/s)');grid on
legend('Y Velocity')
subplot(3,1,3)
plot(t1,z_dot);xlabel('Time (sec)');ylabel('Z_dot (m/s)');grid on
legend('Z Velocity')

%Eular Angles (Radians)
figure()
subplot(3,1,1)
plot(t2,phi);xlabel('Time (sec)');ylabel('\phi (rad)');grid on
title('Eular Angles (Radians)')
legend('Roll')
subplot(3,1,2)
plot(t2,theta);xlabel('Time (sec)');ylabel('\theta (rad)');grid on
legend('Pitch')
subplot(3,1,3)
plot(t2,psi);xlabel('Time (sec)');ylabel('\psi (rad)');grid on
legend('Yaw')

%Eular Angles (Degrees)
figure()
subplot(3,1,1)
plot(t2,D_phi);xlabel('Time (sec)');ylabel('\phi (degree)');grid on
title('Eular Angles (Degrees)')
legend('Roll')
subplot(3,1,2)
plot(t2,D_theta);xlabel('Time (sec)');ylabel('\theta (degree)');grid on
legend('Pitch')
subplot(3,1,3)
plot(t2,D_psi);xlabel('Time (sec)');ylabel('\psi (degree)');grid on
legend('Yaw')

%Quaternions
figure()
subplot(4,1,1)
plot(t2,quat(:,1));xlabel('Time (sec)');ylabel('magnitude');grid on
title('Quaternions')
legend('w scaler')
subplot(4,1,2)
plot(t2,quat(:,2));xlabel('Time (sec)');ylabel('magnitude');grid on
legend('x scaler')
subplot(4,1,3)
plot(t2,quat(:,3));xlabel('Time (sec)');ylabel('magnitude');grid on
legend('y scaler')
subplot(4,1,4)
plot(t2,quat(:,4));xlabel('Time (sec)');ylabel('magnitude');grid on
legend('z scaler')

%Angular Velocity
figure()
subplot(3,1,1)
plot(t3,phi_dot);xlabel('Time (sec)');ylabel('\phi_dot (m/s^2)');grid on
legend('X Angular Velocity')
title('Angular Velocity')
subplot(3,1,2)
plot(t3,theta_dot);xlabel('Time (sec)');ylabel('\theta_dot (m/s^2)');grid on
legend('Y Angular Velocity')
subplot(3,1,3)
plot(t3,psi_dot);xlabel('Time (sec)');ylabel('\psi_dot (m/s^2)');grid on
legend('Z Angular Velocity')

%Control Inputs
figure()
subplot(4,1,1)
plot(t4,controlinputs(:,1));xlabel('Time (sec)');ylabel('Value');grid on
legend('F Thrust')
title('Control Inputs')
subplot(4,1,2)
plot(t4,controlinputs(:,2));xlabel('Time (sec)');ylabel('Value');grid on
legend('X Moment')
subplot(4,1,3)
plot(t4,controlinputs(:,3));xlabel('Time (sec)');ylabel('Value');grid on
legend('Y Moment')
subplot(4,1,4)
plot(t4,controlinputs(:,3));xlabel('Time (sec)');ylabel('Value');grid on
legend('Z Moment')