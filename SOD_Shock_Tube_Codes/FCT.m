%% CFD - Final Project - Flux Corrected Transport

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

E_l=zeros(3,size(X,2)-1);
E_h=zeros(3,size(X,2)-1);
A=zeros(3,size(X,2)-1);

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

%% Flux corrected transport (FCT)

for t=1:size(T,2)
    Q_old=Q;
    
    % Finding the lower and higher order flux values, calculating their
    % difference
    for l=1:3
        E_l(l,1)=E(l,1);
        E_l(l,size(X,2)-1)=E(l,size(X,2)-1);
        E_h(l,1)=E(l,1);
        E_h(l,size(X,2)-1)=E(l,size(X,2)-1);
        
        A(l,1)=E_h(l,1)-E_l(l,1);
        A(l,size(X,2)-1)=E_h(l,size(X,2)-1)-E_l(l,size(X,2)-1);
        for i=2:size(X,2)-2
            % Lower order flux
            E_l(l,i)=E(l,i);
            % Higher order flux
            E_h(l,i)=(7/12)*(E(l,i+1)+E(l,i))-(1/12)*(E(l,i+2)+E(l,i-1));
            
            % Difference of higher and lower order flux
            A(l,i)=E_h(l,i)-E_l(l,i);
        end
    end
    
    % Intermediate solution for Q using lower order scheme
    Q_i=Q_old;
    for l=1:3
        for i=2:size(X,2)-1
            Q_i(l,i)=Q_old(l,i)-(dt/dx)*(E_l(l,i)-E_l(l,i-1));
        end
    end
    
    % Correcting the anti-diffusive flux
    A_c=A;
    for l=1:3
        for i=2:size(X,2)-2
            si=A(l,i)/abs(A(l,i));
            A_c(l,i)=si*max([0,min([abs(A(l,i)),si*((Q_i(l,i+2)-Q_i(l,i+1))/dx),si*((Q_i(l,i)-Q_i(l,i-1))/dx)])]);
        end
    end
    
    % Final solution using corrected flux
    for l=1:3
        for i=2:size(X,2)-1
            Q(l,i)=Q_i(l,i)-(dt/dx)*(A_c(i)-A_c(i-1));
        end
    end
    
    % Updating flux
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
