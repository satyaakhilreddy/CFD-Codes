%% CFD - MEEN 689 -- MID-TERM PROJECT
%  STEADY STATE HEAT CONDUCTION IN 2-D SPACE
%  Method - ADI
%function [T,err]=ADI(Nx,Ny)
clc;

%% Parameters

%alpha=11.234e-5;
%k=280;
L=0.4;
W=0.3;
T1=40;
T2=0;
T3=10;
T4=0;

Nx=10; % Number of elements in X-direction (Variation of columns in matrix)
Ny=10; % Number of elements in Y-direction (Variation of rows in matrix)
itr_max=500;

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

%% Initial Conditions Setup

T_nd=zeros(Ny,Nx);
T_nd(1,2:Nx)=T3_nd;
T_nd(Ny,2:Nx)=T1_nd;
T_nd(:,1)=T2_nd;
T_nd(:,Nx)=T4_nd;

%% ADI method for solution

beta=del_xnd/del_ynd;

ax=zeros(Nx,1);
bx=zeros(Nx,1);
cx=zeros(Nx,1);
dx=zeros(Nx,1);

ay=zeros(Ny,1);
by=zeros(Ny,1);
cy=zeros(Ny,1);
dy=zeros(Ny,1);

T_str=zeros(Ny,Nx);
T_str(1,2:Nx)=T3_nd;
T_str(Ny,2:Nx)=T1_nd;
T_str(:,1)=T2_nd;
T_str(:,Nx)=T4_nd;

tic
for itr=1:itr_max
    T_ndold=T_nd;
    
    % X-sweep 
    for i=Ny-1:-1:2
        for jj=2:Nx-1
            ax(jj)=1;
            bx(jj)=-2*(1+beta^2);
            cx(jj)=1;
            dx(jj)=-beta^2*(T_str(i+1,jj)+T_ndold(i-1,jj));
        end
        
        dx(2)=dx(2)-T_nd(i,1);
        dx(Nx-1)=dx(Nx-1)-T_nd(i,Nx);
        
        dx=TDMA(2,Nx-1,ax,bx,cx,dx);
        
        for p=2:(Nx-1)
            T_str(i,p)=dx(p);
        end
    end
    
    T_str(1,2:Nx)=T3_nd;
    T_str(Ny,2:Nx)=T1_nd;
    T_str(:,1)=T2_nd;
    T_str(:,Nx)=T4_nd;

    
    % Y-sweep
    for j=2:Nx-1
        for ii=2:Ny-1
            ay(ii)=beta^2;
            by(ii)=-2*(1+beta^2);
            cy(ii)=beta^2;
            dy(ii)=-(T_str(ii,j+1)+T_nd(ii,j-1));
        end
        
        dy(2)=dy(2)-beta^2*T_nd(1,j);
        dy(Ny-1)=dy(Ny-1)-beta^2*T_nd(Ny,j);
        
        dy=TDMA(2,Ny-1,ay,by,cy,dy);
        
        for p=2:(Ny-1)
            T_nd(p,j)=dy(p);
        end
    end
    
    T_nd(1,2:Nx)=T3_nd;
    T_nd(Ny,2:Nx)=T1_nd;
    T_nd(:,1)=T2_nd;
    T_nd(:,Nx)=T4_nd;
    
    T=T3+T_nd*(T1-T3);
    T_old=T3+T_ndold*(T1-T3);
    
    if sum(sum((T(2:Ny-1,2:Nx-1)-T_old(2:Ny-1,2:Nx-1)).^2))/sum(sum(T(2:Ny-1,2:Nx-1)))<1e-5
        break;
    end
end
toc
%% Analytical Solution

syms m;
x=0:del_x:L;
y=W:-del_y:0;
[X,Y]=meshgrid(x,y);
Ta=((1-cos(m*pi))/(m*pi))*(sinh((m*pi*(W-Y))/L)).*sin((m*pi*X)/L)/(sinh((m*pi*W)/L));
Tb=((1-cos(m*pi))/(m*pi))*(sinh((m*pi*Y)/L)).*sin((m*pi*X)/L)/(sinh((m*pi*W)/L));
T_an=vpa(2*T1*symsum(Ta,m,1,100)+2*T3*symsum(Tb,m,1,100),4);

%% Error

err=sum(sum((T(2:Ny-1,2:Nx-1)-T_an(2:Ny-1,2:Nx-1)).^2))/(Nx*Ny);

% %% Visualization - T
% 
% figure;
% contourf(0:del_x:L,0:del_y:W,flipud(T))
% colorbar
% title('Contour plot for ADI')
% xlabel('X')
% ylabel('Y')
% 
% %% Visualization - T_an
% 
% figure;
% contourf(0:del_x:L,0:del_y:W,flipud(T_an))
% colorbar
% title('Contour plot for Analytical solution')
% xlabel('X')
% ylabel('Y')

%end
