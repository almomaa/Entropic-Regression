clearvars; close all; clc;
rng(5)
% Lorenz system ODE function
Lorenz = @(t,x) [10.*(x(2)-x(1));...
                 28.*x(1) - x(1).*x(3) - x(2);...
                 x(1).*x(2) - (8/3).*x(3);];
             
% Now, we solve the forward system:
tao = 0.001;         % time step size
x0  = rand(3,1);    % initial condition
N   = 1500;         % number of observations to obtain

[~,X] = ode45(Lorenz, 0:tao:100, x0); %Transient
x0 = X(end,:)';
[tt,X] = ode45(Lorenz, 0:tao:(N+2)*tao, x0); %Sampled Data


% Now, our data (state variable) stored in the matrix X, and shown in the following figure.
plot3(X(:,1),X(:,2),X(:,3),'.b')
view(20,20)
xlabel('X'); ylabel('Y'); zlabel('Z');
box on; grid on;
axis([-25 25 -25 25 0 50])

% Noise and Outliers
eps1 = 1e-4; % base noise std
eps2 = 0.1;  % outliers std
p    = 0.2;  % corruption probability

X = X + eps1*randn(size(X));

IX = rand(size(X,1),1)<= p;
% IX: Boolian string with the index of corrpted observations
X(IX,:) = X(IX,:) + eps2*randn(size(X(IX,:))); % add outliers noise


% Derivative Estimation (Central difference)
Xdot = (1/(2*tao))*(X(3:end,:)-X(1:end-2,:));
X(1,:) = []; X(end,:) = [];

figure
plot(tt(1:end-2),Xdot,'-')
axis([1 max(tt) -200 200])
xlabel('time ($t$)','interpreter','latex'); 
set(gca,'FontSize',15)
legend('$\dot{x}$','$\dot{y}$','$\dot{z}$',...
       'interpreter','latex','Location','northwest',...
       'FontSize',18)

% Inverse Problem
Phi = polyspace(X,5); % power polynomial expansion of order 5

% Entropic Regression
B = erfit(Phi, Xdot);
disp(B)

% Prediction
[~,Xtrue] = ode45(Lorenz, 0:tao:1.5, x0);
G = @(t,x) (polyspace(x',5)*B)';
[~,Xmodel] = ode45(G, 0:tao:4.5, x0); %Sampled Data

figure
plot3(Xtrue(:,1),Xtrue(:,2),Xtrue(:,3),'ob')
hold on
plot3(Xmodel(:,1),Xmodel(:,2),Xmodel(:,3),'.r')

view(20,20)
xlabel('X'); ylabel('Y'); zlabel('Z');
box on; grid on;
axis([-25 25 -25 25 0 50])
legend('Training Data','Re-produced dynamic','Location','best')
set(gca,'FontSize',15)

figure
[~,Xmodel] = ode45(G, 0:tao:10, x0);
plot3(Xmodel(:,1),Xmodel(:,2),Xmodel(:,3),'-b')
view(20,20)
xlabel('X'); ylabel('Y'); zlabel('Z');
box on; grid on;
axis([-25 25 -25 25 0 50])
set(gca,'FontSize',15)



