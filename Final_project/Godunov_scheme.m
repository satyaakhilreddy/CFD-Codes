%% CFD - Final Project - Godunov scheme

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

%% Initial SOD conditions

Q(1,1:(0.5/dx))=1;
Q(1,(0.5/dx)+1:size(X,2))=0.125;

Q(2,1:(0.5/dx))=0;
Q(2,(0.5/dx)+1:size(X,2))=0;

Q(3,1:(0.5/dx))=2.5;
Q(3,(0.5/dx)+1:size(X,2))=0.25;

E(1,:)=Q(2,:);
E(2,:)=((3-g)/2)*((Q(2,:).^2)./Q(1,:))+(g-1)*Q(3,:);
E(3,:)=(Q(2,:)./Q(1,:)).*(g*Q(3,:)-((g-1)/2)*((Q(2,:).^2)./Q(1,:)));

% E(1,:)=Q(2,:);
% E(2,:)=((g-3)/2)*Q(1,:).*Q(2,:)+(g-1)*Q(3,:);
% E(3,:)=(Q(2,:)./Q(1,:)).*(g*Q(3,:)-((g-1)/2)*Q(1,:).*Q(2,:));

%% Velocity, Pressure and density from Q and E vectors

u=Q(2,:)./Q(1,:);
rho=Q(1,:);
et=Q(3,:)./Q(1,:);
e=et-u.^2/2;
p=e.*rho*(g-1);

tic
%% Godunov Scheme

for t=1:size(T,2)
%t=0;
%while t<T_final
    Q_old=Q;
    for l=1:3
        % Finding Lax-Friedrich flux
        for i=1:(size(X,2)-1)
            % Velocity of sound in their respective domains
            cl=sqrt((g*p(i))/rho(i));
            cr=sqrt((g*p(i+1))/rho(i+1));
            
            % Finding the maximum eigen value
            A=[abs(u(i)),abs(u(i)+cl),abs(u(i)-cl),abs(u(i+1)),abs(u(i+1)+cr),abs(u(i+1)-cr)];
            alpha(i)=max(A);
            
            % Lax-Friedrich flux
            F(l,i)=0.5*(E(l,i)+E(l,i+1))-0.5*alpha(i)*(Q_old(l,i+1)-Q_old(l,i));
        end
        
        % Time step determination
        %dt=cfl*(dx/max(alpha));
        
        % Time Marching
        for j=2:size(X,2)-1
            Q(l,j)=Q_old(l,j)-(dt/dx)*(F(l,j)-F(l,j-1));
        end
    end
    
    % Updating the flux
    E(1,:)=Q(2,:);
    E(2,:)=((3-g)/2)*((Q(2,:).^2)./Q(1,:))+(g-1)*Q(3,:);
    E(3,:)=(Q(2,:)./Q(1,:)).*(g*Q(3,:)-((g-1)/2)*((Q(2,:).^2)./Q(1,:)));

%     E(1,:)=Q(2,:);
%     E(2,:)=((g-3)/2)*Q(1,:).*Q(2,:)+(g-1)*Q(3,:);
%     E(3,:)=(Q(2,:)./Q(1,:)).*(g*Q(3,:)-((g-1)/2)*Q(1,:).*Q(2,:));
    
    % Updating velcity, pressure and density vectors
    u=Q(2,:)./Q(1,:);
    rho=Q(1,:);
    et=Q(3,:)./Q(1,:);
    e=et-u.^2/2;
    p=e.*rho*(g-1);
    
    %t=t+dt;
end
toc

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
title('Velocity - Gudonov Scheme');
legend('Godunov Solution','Exact solution')

figure
plot(X,rho,'g-','LineWidth',2);
hold on
plot(X,rho_e,'k--');
grid on
xlabel('X');
ylabel('Density');
title('Density - Gudonov Scheme');
legend('Godunov Solution','Exact solution')

figure
plot(X,p,'b-','LineWidth',2);
hold on
plot(X,p_e,'k--');
grid on
xlabel('X');
ylabel('Pressure');
title('Pressure - Gudonov Scheme');
legend('Godunov Solution','Exact solution')