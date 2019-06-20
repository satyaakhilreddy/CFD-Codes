%% CFD - MEEN 689 -- MID-TERM PROJECT
%  STEADY STATE HEAT CONDUCTION IN 2-D SPACE
%  Method - SOR

clc;
clear all

%% Parameters

%alpha=11.234e-5;
%k=280;
L=0.4;
W=0.3;
T1=40;
T2=0;
T3=10;
T4=0;

Nx=7; % Number of elements in X-direction (Variation of columns in matrix)
Ny=7; % Number of elements in Y-direction (Variation of rows in matrix)
itr_max=500;
w_max=2;
w_d=0.05;
w_min=0.05;

% The PDE is elliptic, therefore the information is transferred
% instantaneously.

%% Non-dimensionalization

% x=L*x_nd
% y=L*x_nd
% T_nd=(T-T3)/(T1-T3)

T1_nd=(T1-T3)/(T1-T3);
T2_nd=(T2-T3)/(T1-T3);
T3_nd=(T3-T3)/(T1-T3);
T4_nd=(T4-T3)/(T1-T3);

del_x=L/(Nx-1);
del_xnd=del_x/L;

del_y=W/(Ny-1);
del_ynd=del_y/L;

%% Analytical Solution

syms m;
x=0:del_x:L;
y=W:-del_y:0;
[X,Y]=meshgrid(x,y);
Ta=((1-cos(m*pi))/(m*pi))*(sinh((m*pi*(W-Y))/L)).*sin((m*pi*X)/L)/(sinh((m*pi*W)/L));
Tb=((1-cos(m*pi))/(m*pi))*(sinh((m*pi*Y)/L)).*sin((m*pi*X)/L)/(sinh((m*pi*W)/L));
T_an=vpa(2*T1*symsum(Ta,m,1,100)+2*T3*symsum(Tb,m,1,100),4);

%% SOR solution method

beta=del_xnd/del_ynd;
% iteration=zeros(1,(w_max-w_min)/w_d+1);
% error=zeros(1,(w_max-w_min)/w_d+1);
% z=0;
% for w=w_min:0.05:w_max
%     z=z+1;
%     
%     % Initial Conditions Setup
%     T_nd=zeros(Ny,Nx);
%     T_nd(1,:)=T3_nd;
%     T_nd(Ny,:)=T1_nd;
%     T_nd(2:Ny-1,1)=T2_nd;
%     T_nd(2:Ny-1,Nx)=T4_nd;
% 
%     T_bar=zeros(Ny,Nx);
%     T_bar(1,:)=T3_nd;
%     T_bar(Ny,:)=T1_nd;
%     T_bar(2:Ny-1,1)=T2_nd;
%     T_bar(2:Ny-1,Nx)=T4_nd;
%     
%     for itr=1:itr_max
%         T_ndold=T_nd;
%         for i=Ny-1:-1:2
%             for j=2:Nx-1
%                 T_bar(i,j)=(1/(2*(1+beta^2)))*(T_bar(i,j-1)+T_ndold(i,j+1)+beta^2*(T_ndold(i-1,j)+T_bar(i+1,j)));
%                 T_nd=T_ndold+w*(T_bar-T_ndold);
%             end
%         end
%         
%         T_nd(1,2:Nx)=T3_nd;
%         T_nd(Ny,2:Nx)=T1_nd;
%         T_nd(:,1)=T2_nd;
%         T_nd(:,Nx)=T4_nd;
%         
%         T=T3+T_nd*(T1-T3);
%         T_old=T3+T_ndold*(T1-T3);
%         
%         if sum(sum((T(2:Ny-1,2:Nx-1)-T_old(2:Ny-1,2:Nx-1)).^2))/sum(sum(T_old(2:Ny-1,2:Nx-1)))<1e-5
%             iteration(z)=itr;
%             error(z)=sum(sum((T(2:Ny-1,2:Nx-1)-T_an(2:Ny-1,2:Nx-1)).^2))/(Nx*Ny);
%             break;
%         end
%         
%         if itr==itr_max
%             iteration(z)=itr;
%             error(z)=sum(sum((T(2:Ny-1,2:Nx-1)-T_an(2:Ny-1,2:Nx-1)).^2))/(Nx*Ny);
%         end
%     end
% end
% 
% figure;
% plot(w_min:w_d:w_max,error,'r-')
% title('Error vs Relaxation Coefficient')
% xlabel('Relaxation Coefficient')
% ylabel('Error')
% set(gca,'XMinorTick','on','YMinorTick','on')
% grid on
% 
% figure
% plot(w_min:w_d:w_max,iteration,'b-')
% title('Iterations vs Relaxation Coefficient')
% xlabel('Relaxation Coefficient')
% ylabel('Iterations')
% set(gca,'XMinorTick','on','YMinorTick','on')
% grid on

% Initial Conditions Setup
T_nd=zeros(Ny,Nx);
T_nd(1,:)=T3_nd;
T_nd(Ny,:)=T1_nd;
T_nd(2:Ny-1,1)=T2_nd;
T_nd(2:Ny-1,Nx)=T4_nd;

T_bar=zeros(Ny,Nx);
T_bar(1,:)=T3_nd;
T_bar(Ny,:)=T1_nd;
T_bar(2:Ny-1,1)=T2_nd;
T_bar(2:Ny-1,Nx)=T4_nd;

w=1.4;
tic
for itr=1:itr_max
    T_ndold=T_nd;
    for i=Ny-1:-1:2
        for j=2:Nx-1
            T_bar(i,j)=(1/(2*(1+beta^2)))*(T_bar(i,j-1)+T_ndold(i,j+1)+beta^2*(T_ndold(i-1,j)+T_bar(i+1,j)));
            T_nd=T_ndold+w*(T_bar-T_ndold);
        end
    end

    T_nd(1,2:Nx)=T3_nd;
    T_nd(Ny,2:Nx)=T1_nd;
    T_nd(:,1)=T2_nd;
    T_nd(:,Nx)=T4_nd;

    T=T3+T_nd*(T1-T3);
    T_old=T3+T_ndold*(T1-T3);
    
    %conv(itr)=sum(sum((T(2:Ny-1,2:Nx-1)-T_old(2:Ny-1,2:Nx-1)).^2))/sum(sum(T_old(2:Ny-1,2:Nx-1)));
    %iter(itr)=itr;

    if sum(sum((T(2:Ny-1,2:Nx-1)-T_old(2:Ny-1,2:Nx-1)).^2))/sum(sum(T_old(2:Ny-1,2:Nx-1)))<1e-5
        err=sum(sum((T(2:Ny-1,2:Nx-1)-T_an(2:Ny-1,2:Nx-1)).^2))/(Nx*Ny);
        break;
    end

    if itr==itr_max
        err=sum(sum((T(2:Ny-1,2:Nx-1)-T_an(2:Ny-1,2:Nx-1)).^2))/(Nx*Ny);
    end
end
toc

% figure;
% plot(iter,conv,'r-')
% title('Convergence Rate')
% xlabel('Iterations')
% ylabel('Convergence')
% grid on

% figure;
% contourf(0:del_x:L,0:del_y:W,flipud(T))
% colorbar
% title('Contour plot for SOR')
% xlabel('X')
% ylabel('Y')
