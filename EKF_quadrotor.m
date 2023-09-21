function [x,P,K] = EKF_quadrotor(x,u,dt,y,P,Q,R)

fstate = @quadrotor;
% hmeas = @(x,u,dt)[x(:)]; 
% hmeas = @(x,u,dt)[x(1:7); x(9); x(11)]; % with GPS coordinates
hmeas = @(x,u,dt)[x(1:12)]; % without GPS coordinates

% for nonlinear dynamic system:
%           x_k+1 = f(x_k) + w_k
%           y_k   = h(x_k) + v_k
% where w ~ N(0,Q) meaning w is gaussian noise with covariance Q
%       v ~ N(0,R) meaning v is gaussian noise with covariance R
% Inputs:   f: function handle for f(x)
%           x: "a priori" state estimate
%           P: "a priori" estimated state covariance
%           h: fanction handle for h(x)
%           y: current measurement
%           Q: process noise covariance 
%           R: measurement noise covariance
% Output:   x: "a posteriori" state estimate
%           P: "a posteriori" state covariance


[x1,A]=jaccsd(fstate,x,u,dt);    %nonlinear update and linearization at current state
P=A*P*A'+Q;                 %partial update
[y1,H]=jaccsd(hmeas,x1,u,dt);    %nonlinear measurement and linearization
P12=P*H';                   %cross covariance
K=P12*inv(H*P12+R);       %Kalman filter gain
x=x1+K*(y-y1);            %state estimate
P=P-K*P12';               %state covariance matrix

% R=chol(H*P12+R);            %Sk, Cholesky factorization
% U=P12/R;                    %K=U/R'; Faster because of back substitution
% x=x1+U*(R'\(y-y1));         %Back substitution to get state update
% K = U*inv(R');
% P=P-U*U';                   %Covariance update, U*U'=P12/R/R'*P12'=K*P12.

function [z,A]=jaccsd(fun,x,u,dt)
% JACCSD Jacobian through complex step differentiation
% [z J] = jaccsd(f,x)
% z = f(x)
% J = f'(x)

z=fun(x,u,dt);
n=numel(x);
m=numel(z);
A=zeros(m,n);
slack = 1e-4;

for k=1:n
    x1=x;
    x1(k)=x1(k)+slack;
    A(:,k)=(fun(x1,u,dt)-z)/slack;
end