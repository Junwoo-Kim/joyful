% RLS로 등가 외란 파라미터 추정하기 Case 8 - 외란 0, 파라미터 변화
% 2019.08.26 CDSP@mju
% 김준우
% close all
clear all 
clc

global b D1 D2 J w0 Pm L_pos L_neg d FT

%% Definition
t0 = 0;
tf = 25;
b = 950e+3; % critical value constant (W)
D1 = 75e-3*6; % viscous damping constant
D2 = 95; % amortisseur damping constant
J = 550; % combined inertia moment (kg*m^2)
w0 = 2*pi*60; % system angular frequency (rad/s)
Pm = 475e+3; % mechanical power (W)
FT = 5; % Fault Time
d = 0;
%% Gain Design
% sin(del0) = 1
A1 = [0 1 0; 0 -D2/J 1/(J*w0); 0 0 0];
C = [1 0 0];

% sin(del0) = -1
A2 = [0 1 0; 0 -D2/J -1/(J*w0); 0 0 0];

des_poles = [-10 -10 -10]*71;

L_pos = acker(A1', C', des_poles)'
L_neg = acker(A2', C', des_poles)'

%% Simulation

dt = 1e-3; % period

step_size = tf/dt; % sample number

x = zeros(5, step_size); % state
x(:,1) = [0.5107 0 0.5107 0 0]';

substeps = 1;

t = 0:dt:tf; % real time scale

flag = 2; % parameter change flag

omega_dot = zeros(step_size+1, 1);
de = zeros(step_size+1, 1);

% RLS Initializing
A = zeros(4,1);
P = eye(size(A, 1));
P_plot = zeros(length(t), 1);
X = ones(size(A, 1), step_size+1);
lambda = 0.99;

% ANC_RLS Initializing
N = 50;
w = zeros(N, 1);
n1 = zeros(N, 1);
P1 = eye(N);
lambda1 = 0.999;
ANC_e = zeros(step_size+1, 1);
ANC_y = zeros(step_size+1, 1);

for i = 1:step_size
    if flag == 0
        if t(i) >= FT
%             d = 445e+3;
            d = 432e+3;
            flag = 1;
        end
    end
%     if flag == 0
%         if t(i) >= FT
%             D1 = D1/2;
%             D2 = D2/2;
%             J = J/2;
%             flag = 1;
%         end
%     end
    
    for j = 1:substeps
        k1 = SwingModel(x(:, i));
        k2 = SwingModel(x(:, i)+dt/2*k1);
        k3 = SwingModel(x(:, i)+dt/2*k2);
        k4 = SwingModel(x(:, i)+dt*k3);
        x(:, i+1) = x(:, i) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    end
    
    delta = x(1,i);
    omega = x(2,i);
    delta_hat = x(3,i);
    omega_hat = x(4,i);
    de_hat = x(5,i);
    
    omega_dot(i+1) = ((-(b-d)*sin(delta))-(D1*omega*(2*w0+omega))+(Pm-D1*w0*w0))/(J*(w0+omega))-((D2*omega)/J);
    de(i+1) = d - (J*omega*omega_dot(i) + D2*omega^2 + D1*(omega + w0)^2)/sin(delta);
    
    % ANC_RLS
    n1_temp = de(i+1) - de_hat; %delta - delta_hat;
    n1 = [n1_temp; n1(1:N-1)];
    
    K = (P1*n1)/(lambda1+n1'*P1*n1);
    P1 = (P1 - K*n1'*P1)/lambda1;
    ANC_e(i+1) = de_hat - n1'*w;
    w = w + K*(de_hat - n1'*w);
    ANC_y(i+1) = n1'*w;
    

    
    % RLS
    A = [sin(delta) -omega.*omega_dot(i+1) -(w0+omega).^2 -omega.^2]';
    S = de(i+1)*sin(delta);
%     S = de_hat*sin(delta) + ANC_y(i+1)*sin(delta);
%     S = de_hat*sin(delta);
    
    K = (P*A)/(lambda+A'*P*A);
    P = (P - K*A'*P)/lambda;
    X(:,i+1) = X(:,i) + K*(S - A'*X(:,i));
    
    P_plot(i, :) = norm(P, 2);
end
x = x';

delta = x(:,1);
omega = x(:,2);
delta_hat = x(:,3);
omega_hat = x(:,4);
de_hat = x(:,5);
%% ANC Plot
figure(91)
plot(t, ANC_e - de, 'linewidth', 2); grid on;
xlabel('time [sec]');
ylabel('ANC e - de');
title('ANC e - de', 'Fontsize', 13);

figure(92)
plot(t, ANC_y, 'linewidth', 2); grid on;
xlabel('time [sec]');
ylabel('ANC y');
title('ANC y', 'Fontsize', 13);

%% plot
figure(1)
subplot(211)
plot(t, delta*(180/pi), t, delta_hat*(180/pi), '--', 'linewidth', 2); grid on;
% axis([0 tf -2e+3 2.5e+4]);
legend('\delta', '\delta hat', 'location', 'southeast');
xlabel('time [sec]');
ylabel('\delta [deg]');
title('States Estimation', 'Fontsize', 13);
subplot(212)
plot(t, omega*(180/pi), t, omega_hat*(180/pi), '--', 'linewidth', 2); grid on;
% axis([0 tf -50 800]);
legend('\omega', '\omega hat', 'location', 'southeast');
xlabel('time [sec]');
ylabel('\omega [deg/sec]');

figure(2)
plot(t, omega_dot, 'linewidth', 2); grid on;
xlabel('time [sec]');
ylabel('omega dot');
title('omega dot', 'Fontsize', 13);

figure(3)
subplot(211)
plot(t, de, t, de_hat, '--', 'linewidth', 2); grid on;
legend('de', 'dehat', 'location', 'southeast');
xlabel('time [sec]');
ylabel('de & dehat');
title('de & dehat', 'Fontsize', 13);
subplot(212)
plot(t, de - de_hat, 'linewidth', 2); grid on;
xlabel('time [sec]');
ylabel('de - dehat');
title('de, dehat error', 'Fontsize', 13);

% dhat
figure(4)
plot(t, X(1,:)*1e-3, t(length(t)), d*1e-3, '*', 'linewidth', 2); grid on;
axis([0 tf -100 600]);
xlabel('time [sec]');
ylabel(' disturbance ');
title('RLS Estimated d', 'Fontsize', 13);

% Jhat
figure(5)
plot(t, X(2,:), t(length(t)), J, '*', 'linewidth', 2); grid on;
legend('Jhat', 'J', 'location', 'southeast');
axis([0 tf -800 1.6e+3]);
xlabel('time [sec]');
ylabel('J');
title('RLS Estimated J', 'Fontsize', 13);

% D1hat
figure(6)
plot(t, X(3,:), t(length(t)), D1, '*', 'linewidth', 2); grid on;
legend('D1hat', 'D1', 'location', 'southeast');
axis([0 tf 0 0.9]);
xlabel('time [sec]');
ylabel('D1');
title('RLS Estimated D1', 'Fontsize', 13);

% D2hat
figure(7)
plot(t, X(4,:), t(length(t)), D2, '*', 'linewidth', 2); grid on;
legend('D2hat', 'D2', 'location', 'southeast');
axis([0 tf 0 500]);
xlabel('time [sec]');
ylabel('D2');
title('RLS Estimated D2', 'Fontsize', 13);

% P norm
figure(8)
plot(t, P_plot, 'linewidth', 2); grid on;
xlabel('time [sec]');
ylabel('P norm');
title('P norm', 'Fontsize', 13);

LS_d = X(1,length(t))*1e-3
LS_J = X(2,length(t))
LS_D1 = X(3,length(t))
LS_D2 = X(4,length(t))
