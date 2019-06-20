%% CFD - Final Project - Godunov with MUSCL scheme

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
E(2,:)=((g-3)/2)*Q(1,:).*Q(2,:)+(g-1)*Q(3,:);
E(3,:)=(Q(2,:)./Q(1,:)).*(g*Q(3,:)+((g-1)/2)*Q(1,:).*Q(2,:));

%% Velocity, Pressure and density from Q and E vectors

u=Q(2,:)./Q(1,:);
rho=Q(1,:);
p=(Q(3,:)-0.5*(Q(2,:)./Q(1,:)))*(g-1);

%% Godunov - MUSCL Scheme

for t=1:size(T,2)
    Q_old=Q;
    
    % Approximations of velocity, density and pressure at the i+1/2
    % boundaries
    for i=1:size(X,2)-1
        if i==1
            u_l(i)=u(i)+((u(i+1)-u(i))/dx)*(dx/2);
            rho_l(i)=rho(i)+((rho(i+1)-rho(i))/dx)*(dx/2);
            p_l(i)=p(i)+((p(i+1)-p(i))/dx)*(dx/2);
            
            u_r(i)=u(i+1)-minmod(((u(i+2)-u(i+1))/dx),((u(i+1)-u(i))/dx))*(dx/2);
            rho_r(i)=rho(i+1)-minmod(((rho(i+2)-rho(i+1))/dx),((rho(i+1)-rho(i))/dx))*(dx/2);
            p_r(i)=p(i+1)-minmod(((p(i+2)-p(i+1))/dx),((p(i+1)-p(i))/dx))*(dx/2);
        else
            if i==size(X,2)-1
                u_l(i)=u(i)+minmod(((u(i+1)-u(i))/dx),((u(i)-u(i-1))/dx))*(dx/2);
                rho_l(i)=rho(i)+minmod(((rho(i+1)-rho(i))/dx),((rho(i)-rho(i-1))/dx))*(dx/2);
                p_l(i)=p(i)+minmod(((p(i+1)-p(i))/dx),((p(i)-p(i-1))/dx))*(dx/2);

                u_r(i)=u(i+1)-((u(i+1)-u(i))/dx)*(dx/2);
                rho_r(i)=rho(i+1)-((rho(i+1)-rho(i))/dx)*(dx/2);
                p_r(i)=p(i+1)-((p(i+1)-p(i))/dx)*(dx/2);
            else
                u_l(i)=u(i)+minmod(((u(i+1)-u(i))/dx),((u(i)-u(i-1))/dx))*(dx/2);
                rho_l(i)=rho(i)+minmod(((rho(i+1)-rho(i))/dx),((rho(i)-rho(i-1))/dx))*(dx/2);
                p_l(i)=p(i)+minmod(((p(i+1)-p(i))/dx),((p(i)-p(i-1))/dx))*(dx/2);

                u_r(i)=u(i+1)-minmod(((u(i+2)-u(i+1))/dx),((u(i+1)-u(i))/dx))*(dx/2);
                rho_r(i)=rho(i+1)-minmod(((rho(i+2)-rho(i+1))/dx),((rho(i+1)-rho(i))/dx))*(dx/2);
                p_r(i)=p(i+1)-minmod(((p(i+2)-p(i+1))/dx),((p(i+1)-p(i))/dx))*(dx/2);
            end
        end
    end
    
    % Formulation of Q at the boundaries
    for i=1:size(X,2)-1
        Q_l(l,i)=rho_l(i);
        Q_l(2,i)=rho_l(i)*u_l(i);
        Q_l(3,i)=0.5*rho_l(i)*(u_l(i))^2+p_l(i)/(g-1);

        Q_r(l,i)=rho_r(i);
        Q_r(2,i)=rho_r(i)*u_r(i);
        Q_r(3,i)=0.5*rho_r(i)*(u_r(i))^2+p_r(i)/(g-1);
    end
    
    % Formulation of E
    E_l(1,:)=Q_l(2,:);
    E_l(2,:)=((g-3)/2)*Q_l(1,:).*Q_l(2,:)+(g-1)*Q_l(3,:);
    E_l(3,:)=(Q_l(2,:)./Q_l(1,:)).*(g*Q_l(3,:)+((g-1)/2)*Q_l(1,:).*Q_l(2,:));
     
    E_r(1,:)=Q_r(2,:);
    E_r(2,:)=((g-3)/2)*Q_r(1,:).*Q_r(2,:)+(g-1)*Q_r(3,:);
    E_r(3,:)=(Q_r(2,:)./Q_r(1,:)).*(g*Q_r(3,:)+((g-1)/2)*Q_r(1,:).*Q_r(2,:));
    
    
    for l=1:3
        % Finding Lax-Friedrich flux
        for i=1:(size(X,2)-1)
            % Velocity of sound in their respective domains
            cl=sqrt((g*p_l(i))/rho_l(i));
            cr=sqrt((g*p_r(i))/rho_r(i));
            
            % Finding the maximum eigen value
            A=[abs(u_l(i)),abs(u_l(i)+cl),abs(u_l(i)-cl),abs(u_r(i)),abs(u_r(i)+cr),abs(u_r(i)-cr)];
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
    E(2,:)=((g-3)/2)*Q(1,:).*Q(2,:)+(g-1)*Q(3,:);
    E(3,:)=(Q(2,:)./Q(1,:)).*(g*Q(3,:)+((g-1)/2)*Q(1,:).*Q(2,:));
    
    % Updating velcity, pressure and density vectors
    u=Q(2,:)./Q(1,:);
    rho=Q(1,:);
    p=(Q(3,:)-0.5*(Q(2,:)./Q(1,:)))*(g-1);
end

%% Visualization

figure
plot(X,u,'r-','LineWidth',2);
grid on
xlabel('X');
ylabel('Velocity');
title('Velocity - Gudonov Scheme');

figure
plot(X,rho,'g-','LineWidth',2);
grid on
xlabel('X');
ylabel('Density');
title('Density - Gudonov Scheme');

figure
plot(X,p,'b-','LineWidth',2);
grid on
xlabel('X');
ylabel('Pressure');
title('Pressure - Gudonov Scheme');
