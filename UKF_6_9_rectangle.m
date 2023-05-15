clc;clear all;close all

%% %% Linear distance from Velocity

load('linearvelocity.mat')
load('orientationquaternions.mat')
load('angularvelocity.mat')
load('control_inputs.mat')

whos

%% Linear Velocity and Distance
% linearvelocity = [slashmatrice210v2slashvelocity(:,1) slashmatrice210v2slashvelocity(:,2) slashmatrice210v2slashvelocity(:,3)];
x_dot = linearvelocity(:,1);    % Linear Velocity in X Direction
y_dot = linearvelocity(:,2);    % Linear Velocity in Y Direction
z_dot = linearvelocity(:,3);    % Linear Velocity in Z Direction

disp('Linear Velocity ')
disp([x_dot, y_dot, z_dot])

n1 = size(linearvelocity,1);        % Number of Measurement
freq_tran = 100;                     % Frequency of velocity and position sensor
sample_tran = 1/freq_tran;          % Sample time of velocity and position sensor
t1 = sample_tran : sample_tran : n1*sample_tran;
%disp('Duration of t1')
%length(t1)

%% Artificial rectangle trajectory 6:9 meters, 1 meter altitude

% x = cumtrapz(t1,x_dot);  % Position X
% y = cumtrapz(t1,y_dot);  % Position Y
% z = cumtrapz(t1,z_dot);  % Position Z

x = [linspace(0,0,296)'; linspace(0,9,1266-296)'; 9*ones(1266,1); linspace(9,0,1266)'; zeros(1267,1)];
y = [linspace(0,0,296)'; linspace(0,0,1266-296)'; linspace(0,6,1266)'; linspace(6,6,1266)'; linspace(6,0,1267)'];
z = [linspace(0,1,296)'; linspace(1,1,5065-296)'];

Measurement_tran = [z z_dot x x_dot y y_dot];

disp('Distance')
disp([z, x, y])

%% Convert Quaternions data to Eular Angles

quat = [orientationquaternions(:,4) orientationquaternions(:,1) orientationquaternions(:,2) orientationquaternions(:,3)];
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

disp('Euler Angles (Degrees)')
disp([D_phi, D_theta, D_psi])

%n = size(quat,1);           % Number of Measurement
freq_rot = 400;             % Frequency of quaternions and angular velocities
sample_rot = 1/freq_rot;    % Sample time of fused sensor

Measurement_rot1 = [phi theta psi];

%% Angular Velocity

phi_dot = angularvelocity(:,1);     % Angular Velocity Roll X Direction
theta_dot = angularvelocity(:,2);   % Angular Velocity Pitch Y Direction
psi_dot = angularvelocity(:,3);     % angular Velocity Yaw Z Direction

disp('Angular Velocity')
disp([phi_dot, theta_dot, psi_dot])

Measurement_rot2 = [phi_dot theta_dot psi_dot];

%% Unscented KALMAN Filter (UKF)

% initial values
dt = 0.001;           % sampling period
n = 12;               % Number of states
r = 12;               % Number of measurement
Q = (1e-2)*eye(12);  % Process noice covariance
R = eye(12);        % Measurement noice covariance
H = eye(12);          % Measurement matrix

x_vector = [pi/8 0.8 pi/8 0.3 pi/8 0.1 0.4 0.2 0.4 0.3 0.2 0.4]'; % initial state
% x_vector = [0.4 0.4 0.2 0.2 0.3 0.4 pi/8 pi/8 pi/8 0.8 0.3 0.1]';

xhat_plus_ukf = x_vector;                                          % initial state estimation (UKF)
p_plus_ukf = diag(x_vector.^2+0.1);                                % initial state uncertainty (UKF)

u = [u_1 u_2 u_3 u_4];

freq_ci = 50;                 % Frequency of control inputs
sample_ci = 1/freq_ci;        % Sample time of control inputs

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
for a = 2 : sample_tran/sample_rot*size(Measurement_tran,1)
    
    %     for b = sample_rot1 : sample_rot1 : sample_tran
    %         y = [Measurement_rot1(k,1); Measurement_rot2(k,1); Measurement_rot1(k,2);...
    %             Measurement_rot2(k,2); Measurement_rot1(k,3); Measurement_rot2(k,3); Measurement_tran(a,:)'];
    %         k = k + 1;
    %     end
    
    %% Positions and velocities measures recorded at 100 Hz
    
    if mod(a,sample_tran/sample_rot) == 0
        Measurement_tran_cur = Measurement_tran(floor(a*sample_rot/sample_tran),:); % takes a position and velocity measures each a/4 round number sample as a corresponds to quaternions frequency
    else
        Measurement_tran_cur = Measurement_tran(floor(a*sample_rot/sample_tran+1),:); % keeps the previous measurement value while a/4 is non-integer number
    end
    
    %% Control inputs recorded at 50 Hz
    
    if mod(a,sample_ci/sample_rot) == 0
        u_cur = u(floor(a*sample_rot/sample_ci),:); % same principle as for position and velocity but for a/8
    else
        u_cur = u(floor(a*sample_rot/sample_ci+1),:);
    end
    
    %% Measures used for UKF
    
    y = [Measurement_rot1(a,1); Measurement_rot2(a,1); Measurement_rot1(a,2); Measurement_rot2(a,2);...
        Measurement_rot1(a,3); Measurement_rot2(a,3); Measurement_tran_cur'];
    
    [x_hat, P, K_ukf] = UKF_quadrotor(x_hat,u_cur,dt,y,P,Q,R);
    
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
t_disp = (1:a)*sample_rot;

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

