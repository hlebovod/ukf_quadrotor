%% Quadrotor
clc;clear all;close all
%%

% initial values
dt = 0.01; % sampling period
tf = 10;  % final time
t = 0:dt:tf;
a = length(t);
n = 12;               % Number of states
r = 12;               % Number of measurement
Q = (1e-12)*eye(12);  % Process noice covariance
R = 1*eye(12);        % Measurement noice covariance
% Q = diag([0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001]);  % Process noice covariance
% R = diag([[0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05]]);                         % Measurement noice covariance
H = eye(12);          % Measurement matrix

%x=[roll(1) roll_dot(2) pitch(3) pitch_dot(4) yaw(5) yaw_dot(6) X(7) X_dot(8) Y(9) Y_dot(10) Z(11) Z_dot(12)]
x = [pi/8 0.8 pi/8 0.3 pi/8 0.1 0.4 0.2 0.4 0.3 0.2 0.4]'; % initial state
xhat_plus_ekf = x;                                          % initial state estimation (UKF)
p_plus_ekf = diag(x.^2+0.1);                                % initial state uncertainty (UKF)
W_ukf = ones(2*n,1) / (2*n);                                % UKF weights

% omega1(1)=0.001;
% omega2(1)=0.001;
% omega3(1)=0.001;
% omega4(1)=0.001;

u = [100 1 1 0];

%% Unscented KALMAN Filter (UKF)
clc
% Array for saving data system
xArray = [x];
yArray = [zeros(r,1)];
wArray = [zeros(n,1)];
vArray = [zeros(r,1)];
% Array for saving data UKF
xhat_minus_ekfArray = [zeros(n,1)];
p_minus_ekfArray = [p_plus_ekf];
xhat_plus_ekfArray = [zeros(n,1)];
K_ekfArray = [zeros(n,r)];
p_plus_ekfArray = [p_plus_ekf];
trace_p_plus_ukfArray = [];

P = p_minus_ekfArray;
x_hat = xhat_minus_ekfArray;

for k = 1 : a-1
    
    % =============================
    % Simulate the system equation
    % =============================
    w = sqrt(Q)*randn(n,1);                % w is process noise (white noise with zero mean)
    x = quadrotor(x + w,u,dt);
    
    % =================================
    % Simulate the measurement equation
    % =================================
    v = sqrt(R)*randn(r,1);  % v is measurement noise (white noise with zero mean)
    y = H*x + v;
    
    [x_hat, P, K_ekf] = EKF_quadrotor(x_hat,u,dt,y,P,Q,R);
    
    xhat_plus_ekf = x_hat;
    p_plus_ekf = P;
    
    % Saving data system
    xArray = [xArray x];
    yArray = [yArray y];
    wArray = [wArray w];
    vArray = [vArray v];
    
    % Saving data UKF
    xhat_plus_ekfArray = [xhat_plus_ekfArray xhat_plus_ekf];
    K_ekfArray = [K_ekfArray K_ekf];
    p_plus_ekfArray = [p_plus_ekfArray diag(p_plus_ekf)];
    trace_p_plus_ukfArray = [trace_p_plus_ukfArray trace(p_plus_ekf)];
    
end

%% Plot
clc

figure(1)
subplot(2,1,1)
plot(t,wArray(1,:));xlabel('Time(sec)');ylabel('State Noise');grid on
subplot(2,1,2)
plot(t,vArray(1,:));xlabel('Time(sec)');ylabel('Measurement Noise');grid on


figure(2)
subplot(3,1,1)
plot(t,xArray(1,:),'b',t,xhat_plus_ekfArray(1,:),'--r')
title('Orientations')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('\phi Roll(rad)');grid on
subplot(3,1,2)
plot(t,xArray(3,:),'b',t,xhat_plus_ekfArray(3,:),'--r')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('\theta Pitch (rad)');grid on
subplot(3,1,3)
plot(t,xArray(5,:),t,xhat_plus_ekfArray(5,:),'--r')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('\psi Yaw (rad)');grid on


figure(3)
subplot(3,1,1)
plot(t,xArray(2,:),t,xhat_plus_ekfArray(2,:),'--r')
title('Angular Velocities')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('\phi dot (rad/s)');grid on
subplot(3,1,2)
plot(t,xArray(4,:),t,xhat_plus_ekfArray(4,:),'--r')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('\theta dot (rad/s)');grid on
subplot(3,1,3)
plot(t,xArray(6,:),t,xhat_plus_ekfArray(6,:),'--r')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('\psi dot (rad/s)');grid on


figure(4)
subplot(3,1,1)
plot(t,xArray(9,:),t,xhat_plus_ekfArray(9,:),'--r')
title('Positions')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('X (m)');grid on
subplot(3,1,2)
plot(t,xArray(11,:),t,xhat_plus_ekfArray(11,:),'--r')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('Y (m)');grid on
subplot(3,1,3)
plot(t,xArray(7,:),t,xhat_plus_ekfArray(7,:),'--r')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('Z (m)');grid on


figure(5)
subplot(3,1,1)
plot(t,xArray(10,:),t,xhat_plus_ekfArray(10,:),'--r')
title('Linear Velocities')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('X dot (m/s)');grid on
subplot(3,1,2)
plot(t,xArray(12,:),t,xhat_plus_ekfArray(12,:),'--r')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('Y dot (m/s)');grid on
subplot(3,1,3)
plot(t,xArray(8,:),t,xhat_plus_ekfArray(8,:),'--r')
legend('Real','Estimate')
xlabel('Time (sec)');ylabel('Z dot (m/s)');grid on

figure(6)
subplot(3,1,1)
plot(t,xArray(1,:)-xhat_plus_ekfArray(1,:))
xlabel('Time(sec)','linewidth',2);ylabel('\phi Error(rad)')
title('Estimation error orientation along \phi');grid on
subplot(3,1,2)
plot(t,xArray(3,:)-xhat_plus_ekfArray(3,:))
xlabel('Time(sec)','linewidth',2);ylabel('\theta Error(rad)')
title('Estimation error orientation along \theta');grid on
subplot(3,1,3)
plot(t,xArray(5,:)-xhat_plus_ekfArray(5,:))
xlabel('Time(sec)','linewidth',2);ylabel('\psi (rad) Error(rad)')
title('Estimation error orientation along \psi (rad)');grid on

figure(7)
subplot(3,1,1)
plot(t,xArray(9,:)-xhat_plus_ekfArray(9,:))
xlabel('Time(sec)','linewidth',2);ylabel('Error(m)')
title('Estimation error position X');grid on
subplot(3,1,2)
plot(t,xArray(11,:)-xhat_plus_ekfArray(11,:))
xlabel('Time(sec)','linewidth',2);ylabel('Error(m)')
title('Estimation error position Y');grid on
subplot(3,1,3)
plot(t,xArray(7,:)-xhat_plus_ekfArray(7,:))
xlabel('Time(sec)','linewidth',2);ylabel('Error(m)')
title('Estimation error position Z');grid on

figure(8)
plot3(xArray(9,:),xArray(11,:),xArray(7,:),xhat_plus_ekfArray(9,:),xhat_plus_ekfArray(11,:),xhat_plus_ekfArray(7,:),'r','linewidth',1)
legend('Real','Estimate');xlabel('X(m)');ylabel('Y(m)');zlabel('Z(m)')
grid on;
