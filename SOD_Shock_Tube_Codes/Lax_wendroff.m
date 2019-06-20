%% CFD - Final Project - Lax Wendroff

%% Given Parameters

T_final=0.1644;
L=1;
dx=0.01;
g=1.4;
cfl=0.8;

%% Discretization Parameters

dt=0.001644;
T=0:dt:T_final;
X=0:dx:L;

%% Initialization

Q=zeros(3,size(X,2));
E=zeros(3,size(X,2));
F=zeros(3,size(X,2));
alpha=zeros(1,size(X,2)-1);

Q_h=zeros(3,size(X,2)-1);

%% Initial SOD conditions

Q(1,1:(0.5/dx)+1)=1;
Q(1,(0.5/dx)+2:size(X,2))=0.125;

Q(2,1:(0.5/dx)+1)=0;
Q(2,(0.5/dx)+2:size(X,2))=0;

Q(3,1:(0.5/dx)+1)=2.5;
Q(3,(0.5/dx)+2:size(X,2))=0.25;

E(1,:)=Q(2,:);
E(2,:)=((3-g)/2)*((Q(2,:).^2)./Q(1,:))+(g-1)*Q(3,:);
E(3,:)=(Q(2,:)./Q(1,:)).*(g*Q(3,:)-((g-1)/2)*((Q(2,:).^2)./Q(1,:)));

%% Velocity, Pressure and density from Q and E vectors

u=Q(2,:)./Q(1,:);
rho=Q(1,:);
et=Q(3,:)./Q(1,:);
e=et-u.^2/2;
p=e.*rho*(g-1);

tic;
%% Lax Wendroff scheme

for t=1:size(T,2)
    Q_old=Q;
    
    % Predictor Step
    for l=1:3
        for i=1:size(X,2)-1
            Q_h(l,i)=0.5*(Q_old(l,i)+Q_old(l,i+1))-0.5*(dt/dx)*(E(l,i+1)-E(l,i));
        end
    end
    
    % Recalculating the flux
    E_h(1,:)=Q_h(2,:);
    E_h(2,:)=((3-g)/2)*((Q_h(2,:).^2)./Q_h(1,:))+(g-1)*Q_h(3,:);
    E_h(3,:)=(Q_h(2,:)./Q_h(1,:)).*(g*Q_h(3,:)-((g-1)/2)*((Q_h(2,:).^2)./Q_h(1,:)));
    
    % Corrector Step
    for l=1:3
        for i=2:size(X,2)-1
            Q(l,i)=Q_old(l,i)-(dt/dx)*(E_h(l,i)-E_h(l,i-1));
        end
    end
    
    % Updating the flux
    E(1,:)=Q(2,:);
    E(2,:)=((3-g)/2)*((Q(2,:).^2)./Q(1,:))+(g-1)*Q(3,:);
    E(3,:)=(Q(2,:)./Q(1,:)).*(g*Q(3,:)-((g-1)/2)*((Q(2,:).^2)./Q(1,:)));
    
    % Updating velcity, pressure and density vectors
    u=Q(2,:)./Q(1,:);
    rho=Q(1,:);
    et=Q(3,:)./Q(1,:);
    e=et-u.^2/2;
    p=e.*rho*(g-1);
end

toc;
%% Exact Solution

[p_e,u_e,rho_e]=RiemannExact(1,1,0,0.1,0.125,0,1e-4);

%% Visualization

figure
plot(X,u,'r-','LineWidth',2);
hold on
plot(X,u_e,'k--');
grid on
xlabel('X');
ylabel('Velocity');
title('Velocity - Lax-Wendroff Scheme');
legend('Lax-Wendroff Solution','Exact solution')

figure
plot(X,rho,'g-','LineWidth',2);
hold on
plot(X,rho_e,'k--');
grid on
xlabel('X');
ylabel('Density');
title('Density - Lax-Wendroff Scheme');
legend('Lax-Wendroff Solution','Exact solution')

figure
plot(X,p,'b-','LineWidth',2);
hold on
plot(X,p_e,'k--');
grid on
xlabel('X');
ylabel('Pressure');
title('Pressure - Lax-Wendroff Scheme');
legend('Lax-Wendroff Solution','Exact solution')
