%% CFD - Final Project - Godunov with MUSCL scheme

%% Given Parameters

T_final=0.1644;
L=1;
dx=0.01;
g=1.4;
cfl=0.8;

%% Discretization Parameters

dt=0.00001644;
T=0:dt:T_final;
X=0:dx:L;

%% Initialization

Q=zeros(3,size(X,2));
E=zeros(3,size(X,2));
F=zeros(3,size(X,2));
alpha=zeros(1,size(X,2)-1);

u_l=zeros(1,size(X,2)-1);
u_r=zeros(1,size(X,2)-1);
rho_l=zeros(1,size(X,2)-1);
rho_r=zeros(1,size(X,2)-1);
p_l=zeros(1,size(X,2)-1);
p_r=zeros(1,size(X,2)-1);

Q_l=zeros(1,size(X,2)-1);
Q_r=zeros(1,size(X,2)-1);

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

tic
%% Godunov - MUSCL Scheme

for t=1:size(T,2)
    Q_old=Q;
    
    % Approximations of Q at the i+1/2 boundaries
    for i=1:size(X,2)-1
        if i==1
            Q_l(1,i)=Q(1,i);%+((Q(1,i+1)-Q(1,i))/dx)*(dx/2);
            Q_l(2,i)=Q(2,i);%+((Q(2,i+1)-Q(2,i))/dx)*(dx/2);
            Q_l(3,i)=Q(3,i);%+((Q(3,i+1)-Q(3,i))/dx)*(dx/2);
            
            Q_r(1,i)=Q(1,i+1)-minmod(((Q(1,i+2)-Q(1,i+1))/dx),((Q(1,i+1)-Q(1,i))/dx))*(dx/2);
            Q_r(2,i)=Q(2,i+1)-minmod(((Q(2,i+2)-Q(2,i+1))/dx),((Q(2,i+1)-Q(2,i))/dx))*(dx/2);
            Q_r(3,i)=Q(3,i+1)-minmod(((Q(3,i+2)-Q(3,i+1))/dx),((Q(3,i+1)-Q(3,i))/dx))*(dx/2);
        else
            if i==size(X,2)-1
                Q_l(1,i)=Q(1,i)+minmod(((Q(1,i+1)-Q(1,i))/dx),((Q(1,i)-Q(1,i-1))/dx))*(dx/2);
                Q_l(2,i)=Q(2,i)+minmod(((Q(2,i+1)-Q(2,i))/dx),((Q(2,i)-Q(2,i-1))/dx))*(dx/2);
                Q_l(3,i)=Q(3,i)+minmod(((Q(3,i+1)-Q(3,i))/dx),((Q(3,i)-Q(3,i-1))/dx))*(dx/2);

                Q_r(1,i)=Q(1,i+1);%-((Q(1,i+1)-Q(1,i))/dx)*(dx/2);
                Q_r(2,i)=Q(2,i+1);%-((Q(2,i+1)-Q(2,i))/dx)*(dx/2);
                Q_r(3,i)=Q(3,i+1);%-((Q(3,i+1)-Q(3,i))/dx)*(dx/2);
            else
                Q_l(1,i)=Q(1,i)+minmod(((Q(1,i+1)-Q(1,i))/dx),((Q(1,i)-Q(1,i-1))/dx))*(dx/2);
                Q_l(2,i)=Q(2,i)+minmod(((Q(2,i+1)-Q(2,i))/dx),((Q(2,i)-Q(2,i-1))/dx))*(dx/2);
                Q_l(3,i)=Q(3,i)+minmod(((Q(3,i+1)-Q(3,i))/dx),((Q(3,i)-Q(3,i-1))/dx))*(dx/2);

                Q_r(1,i)=Q(1,i+1)-minmod(((Q(1,i+2)-Q(1,i+1))/dx),((Q(1,i+1)-Q(1,i))/dx))*(dx/2);
                Q_r(2,i)=Q(2,i+1)-minmod(((Q(2,i+2)-Q(2,i+1))/dx),((Q(2,i+1)-Q(2,i))/dx))*(dx/2);
                Q_r(3,i)=Q(3,i+1)-minmod(((Q(3,i+2)-Q(3,i+1))/dx),((Q(3,i+1)-Q(3,i))/dx))*(dx/2);
            end
        end
    end
   
    % Formulation of E
    E_l(1,:)=Q_l(2,:);
    E_l(2,:)=((3-g)/2)*((Q_l(2,:).^2)./Q_l(1,:))+(g-1)*Q_l(3,:);
    E_l(3,:)=(Q_l(2,:)./Q_l(1,:)).*(g*Q_l(3,:)-((g-1)/2)*((Q_l(2,:).^2)./Q_l(1,:)));
     
    E_r(1,:)=Q_r(2,:);
    E_r(2,:)=((3-g)/2)*((Q_r(2,:).^2)./Q_r(1,:))+(g-1)*Q_r(3,:);
    E_r(3,:)=(Q_r(2,:)./Q_r(1,:)).*(g*Q_r(3,:)-((g-1)/2)*((Q_r(2,:).^2)./Q_r(1,:)));
    
    
    for l=1:3
        % Finding Lax-Friedrich flux
        for i=1:(size(X,2)-1)
            % Velocity of sound in their respective domains
            cl=sqrt((g*p_l(i))/rho_l(i));
            cr=sqrt((g*p_r(i))/rho_r(i));
            
            % Finding the maximum eigen value
            A=[abs(u(i)),abs(u(i)+cl),abs(u(i)-cl),abs(u(i+1)),abs(u(i+1)+cr),abs(u(i+1)-cr)];
            alpha(i)=max(A);
            
            % Lax-Friedrich flux
            F(l,i)=0.5*(E_l(l,i)+E_r(l,i))-0.5*alpha(i)*(Q_r(l,i)-Q_l(l,i));
        end
        
        % Time Marching
        for j=2:size(X,2)-1
            Q(l,j)=Q_old(l,j)-(dt/dx)*(F(l,j)-F(l,j-1));
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
title('Velocity - Godunov with MUSCL Scheme');
legend('Godunov with MUSCL Solution','Exact solution')

figure
plot(X,rho,'g-','LineWidth',2);
hold on
plot(X,rho_e,'k--');
grid on
xlabel('X');
ylabel('Density');
title('Density - Godunov with MUSCL Scheme');
legend('Godunov with MUSCL Solution','Exact solution')

figure
plot(X,p,'b-','LineWidth',2);
hold on
plot(X,p_e,'k--');
grid on
xlabel('X');
ylabel('Pressure');
title('Pressure - Godunov with MUSCL Scheme');
legend('Godunov with MUSCL Solution','Exact solution')
