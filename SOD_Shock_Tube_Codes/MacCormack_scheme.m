%% CFD - Final Project - MacCormack scheme

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
Q_str=zeros(3,size(X,2));
alpha=zeros(1,size(X,2)-2);

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
p=(Q(3,:)-0.5*(Q(2,:)./Q(1,:)))*(g-1);

%% MacCormack scheme

for t=1:size(T,2)
%t=0;
%while t<T_final
    Q_old=Q;
    
    % Velocity of sound in their respective domains
%     cl=sqrt((g*p(i))/rho(i));
%     cr=sqrt((g*p(i+1))/rho(i+1));
%     
%     for i=1:size(X,2)-1
%         A=[abs(u(i)),abs(u(i)+cl),abs(u(i)-cl),abs(u(i+1)),abs(u(i+1)+cr),abs(u(i+1)-cr)];
%         alpha(i)=max(A);
%     end
    
%     dt=cfl*(dx/max(alpha));
    
    % Predictor step
    Q_str=Q_old;
    for l=1:3
        for i=2:(size(X,2)-1)
            Q_str(l,i)=Q_old(l,i)-(dt/dx)*(E(l,i+1)-E(l,i));
        end
    end
    
    % Re-calculating the fluxes from Q_str
    E_str(1,:)=Q_str(2,:);
    E_str(2,:)=((3-g)/2)*((Q_str(2,:).^2)./Q_str(1,:))+(g-1)*Q_str(3,:);
    E_str(3,:)=(Q_str(2,:)./Q_str(1,:)).*(g*Q_str(3,:)-((g-1)/2)*((Q_str(2,:).^2)./Q_str(1,:)));
    
    % Corrector Step
    for l=1:3
        for i=2:(size(X,2)-1)
            Q(l,i)=0.5*(Q_old(l,i)+Q_str(l,i)-(dt/dx)*(E_str(l,i)-E_str(l,i-1)));
        end
    end
    
    % Updating the flux
    E(1,:)=Q(2,:);
    E(2,:)=((3-g)/2)*((Q(2,:).^2)./Q(1,:))+(g-1)*Q(3,:);
    E(3,:)=(Q(2,:)./Q(1,:)).*(g*Q(3,:)-((g-1)/2)*((Q(2,:).^2)./Q(1,:)));
    
    % Updating velcity, pressure and density vectors
    u=Q(2,:)./Q(1,:);
    rho=Q(1,:);
    p=(Q(3,:)-0.5*(Q(2,:)./Q(1,:)))*(g-1);
    
%     t=t+dt;
end

%% Visualization

figure
plot(X,u,'r-','LineWidth',2);
grid on
xlabel('X');
ylabel('Velocity');
title('Velocity - MacCormack Scheme');

figure
plot(X,rho,'g-','LineWidth',2);
grid on
xlabel('X');
ylabel('Density');
title('Density - MacCormack Scheme');

figure
plot(X,p,'b-','LineWidth',2);
grid on
xlabel('X');
ylabel('Pressure');
title('Pressure - MacCormack Scheme');
