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

n1 = size(linearvelocity,1);        % Number of Measurement
freq_tran = 50;                     % Frequency of velocity sensor
sample_tran = 1/freq_tran;          % Sample time of velocity sensor i.e 0.02sec
t1 = sample_tran : sample_tran : n1*sample_tran;
%disp('Duration of t1')
%length(t1)

x = cumtrapz(t1,x_dot);  % Position X
y = cumtrapz(t1,y_dot);  % Position Y
z = cumtrapz(t1,z_dot);  % Position Z

Measurement_tran = [z z_dot x x_dot y y_dot];

disp('Distance')
disp([z, x, y])

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

%n = size(quat,1);           % Number of Measurement
freq_rot = 100;             % Frequency of fused sensor
sample_rot1 = 1/freq_rot;    % Sample time of fused sensor 
%t2 = sample_rot : sample_rot : n*sample_rot;
%disp('Duration of t2')
%length(t2)

Measurement_rot1 = [phi theta psi];

%% Angular Velocity

AngularVelocity = [slashmatrice210v2slashangularvelocityfused(:,1)   slashmatrice210v2slashangularvelocityfused(:,2)   slashmatrice210v2slashangularvelocityfused(:,3)];
phi_dot = AngularVelocity(:,1);     % Angular Velocity Roll X Direction
theta_dot = AngularVelocity(:,2);   % Angular Velocity Pitch Y Direction 
psi_dot = AngularVelocity(:,3);     % angular Velocity Yaw Z Direction

%Measurement_rot = [phi theta psi phi_dot theta_dot psi_dot];

disp('Angular Velocity')
disp([phi_dot, theta_dot, psi_dot])

%n = size(AngularVelocity,1);    % Number of Measurement
freq_rot = 100;                 % Frequency of Fused sensor
sample_rot = 1/freq_rot;        % Sample time of Fused sensor 
%t3 = sample_rot : sample_rot : n*sample_rot;
%disp('Duration of t3')
%length(t3)

Measurement_rot2 = [phi_dot theta_dot psi_dot];

%% Unscented KALMAN Filter (UKF)

% initial values
dt = 0.001;           % sampling period
n = 12;               % Number of states
r = 12;               % Number of measurement
Q = (1e-3)*eye(12);  % Process noice covariance
R = 1*eye(12);        % Measurement noice covariance 
H = eye(12);          % Measurement matrix

x_vector = [pi/8 0.8 pi/8 0.3 pi/8 0.1 0.4 0.2 0.4 0.3 0.2 0.4]'; % initial state
% x_vector = [0.4 0.4 0.2 0.2 0.3 0.4 pi/8 pi/8 pi/8 0.8 0.3 0.1]';

xhat_plus_ukf = x_vector;                                          % initial state estimation (UKF)
p_plus_ukf = diag(x_vector.^2+0.1);                                % initial state uncertainty (UKF)

u = [controlinputs(:,1) controlinputs(:,2) controlinputs(:,3) controlinputs(:,4)];
ux = controlinputs(:,5);
uy = controlinputs(:,6);

%n = size(controlinputs,1);    % Number of Measurements
freq_rc = 50;                 % Frequency of rc sensor
sample_rc = 1/freq_rc;        % Sample time of rc sensor 
%t4 = sample_rc : sample_rc : n*sample_rc;
%disp('Duration of t4')
%length(t4)

% Array for saving data system
yArray = [zeros(r,1)];

% Array for saving data UKF
xhat_minus_ukfArray = [zeros(n,1)];
p_minus_ukfArray = [p_plus_ukf];
xhat_plus_ukfArray = [zeros(n,1)];
K_ukfArray = [zeros(n,r)];
p_plus_ukfArray = [p_plus_ukf];
trace_p_plus_ukfArray = [];

P = p_minus_ukfArray;
x_hat = xhat_minus_ukfArray;

k = 1;
for a = 1 : size(Measurement_tran,1)
    for b = sample_rot1 : sample_rot1 : sample_tran
        y = [Measurement_rot1(k,1); Measurement_rot2(k,1); Measurement_rot1(k,2);...
            Measurement_rot2(k,2); Measurement_rot1(k,3); Measurement_rot2(k,3); Measurement_tran(a,:)'];
        k = k + 1;
    end
    
    [x_hat, P, K_ukf] = UKF_quadrotor(x_hat,u(a,:),dt,y,P,Q,R);
    
    xhat_plus_ukf = x_hat;
    p_plus_ukf = P;
    
    % Saving data system
    yArray = [yArray y];
    
    % Saving data UKF
    xhat_plus_ukfArray = [xhat_plus_ukfArray xhat_plus_ukf];
    K_ukfArray = [K_ukfArray K_ukf];
    p_plus_ukfArray = [p_plus_ukfArray diag(p_plus_ukf)];
    trace_p_plus_ukfArray = [trace_p_plus_ukfArray trace(p_plus_ukf)];

end

%% Plot
clc
t_disp = [0 t];

figure(1)
subplot(3,1,1)
plot(t_disp,yArray(1,:),'b',t_disp,xhat_plus_ukfArray(1,:),'--r')
title('Orientations')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('\phi Roll(rad)');grid on
subplot(3,1,2)
plot(t_disp,yArray(3,:),'b',t_disp,xhat_plus_ukfArray(3,:),'--r')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('\theta Pitch (rad)');grid on
subplot(3,1,3)
plot(t_disp,yArray(5,:),t_disp,xhat_plus_ukfArray(5,:),'--r')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('\psi Yaw (rad)');grid on


figure(2)
subplot(3,1,1)
plot(t_disp,yArray(2,:),t_disp,xhat_plus_ukfArray(2,:),'--r')
title('Angular Velocities')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('\phi dot (rad/s)');grid on
subplot(3,1,2)
plot(t_disp,yArray(4,:),t_disp,xhat_plus_ukfArray(4,:),'--r')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('\theta dot (rad/s)');grid on
subplot(3,1,3)
plot(t_disp,yArray(6,:),t_disp,xhat_plus_ukfArray(6,:),'--r')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('\psi dot (rad/s)');grid on


figure(3)
subplot(3,1,1)
plot(t_disp,yArray(9,:),t_disp,xhat_plus_ukfArray(9,:),'--r')
title('Positions')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('X (m)');grid on
subplot(3,1,2)
plot(t_disp,yArray(11,:),t_disp,xhat_plus_ukfArray(11,:),'--r')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('Y (m)');grid on
subplot(3,1,3)
plot(t_disp,yArray(7,:),t_disp,xhat_plus_ukfArray(7,:),'--r')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('Z (m)');grid on


figure(4)
subplot(3,1,1)
plot(t_disp,yArray(10,:),t_disp,xhat_plus_ukfArray(10,:),'--r')
title('Linear Velocities')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('X dot (m/s)');grid on
subplot(3,1,2)
plot(t_disp,yArray(12,:),t_disp,xhat_plus_ukfArray(12,:),'--r')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('Y dot (m/s)');grid on
subplot(3,1,3)
plot(t_disp,yArray(8,:),t_disp,xhat_plus_ukfArray(8,:),'--r')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('Z dot (m/s)');grid on

figure(5)
subplot(3,1,1)
plot(t_disp,yArray(1,:)-xhat_plus_ukfArray(1,:))
xlabel('Time(sec)','linewidth',2);ylabel('\phi Error(rad)')
title('Estimation error orientation along \phi');grid on
subplot(3,1,2)
plot(t_disp,yArray(3,:)-xhat_plus_ukfArray(3,:))
xlabel('Time(sec)','linewidth',2);ylabel('\theta Error(rad)')
title('Estimation error orientation along \theta');grid on
subplot(3,1,3)
plot(t_disp,yArray(5,:)-xhat_plus_ukfArray(5,:))
xlabel('Time(sec)','linewidth',2);ylabel('\psi (rad) Error(rad)')
title('Estimation error orientation along \psi (rad)');grid on

figure(6)
subplot(3,1,1)
plot(t_disp,yArray(9,:)-xhat_plus_ukfArray(9,:))
xlabel('Time(sec)','linewidth',2);ylabel('Error(m)')
title('Estimation error position X');grid on
subplot(3,1,2)
plot(t_disp,yArray(11,:)-xhat_plus_ukfArray(11,:))
xlabel('Time(sec)','linewidth',2);ylabel('Error(m)')
title('Estimation error position Y');grid on
subplot(3,1,3)
plot(t_disp,yArray(7,:)-xhat_plus_ukfArray(7,:))
xlabel('Time(sec)','linewidth',2);ylabel('Error(m)')
title('Estimation error position Z');grid on

figure(7)
plot3(Measurement_tran(:,3),Measurement_tran(:,5),Measurement_tran(:,1),'r',xhat_plus_ukfArray(9,:),xhat_plus_ukfArray(11,:),xhat_plus_ukfArray(7,:),'b','linewidth',1.25)
legend('Measurement','Estimated');xlabel('X(m)');ylabel('Y(m)');zlabel('Z(m)')
grid on;

figure(8)
plot(Measurement_tran(:,3),Measurement_tran(:,5),'r',xhat_plus_ukfArray(9,:),xhat_plus_ukfArray(11,:),'b','linewidth',1.25)
legend('Measurement','Estimated');xlabel('X(m)');ylabel('Y(m)')
grid on;

