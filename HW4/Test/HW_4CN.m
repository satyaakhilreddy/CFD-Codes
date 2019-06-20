%% Homework-4 -- Crank Nicholson -- CFD

%% Parameters

L=2;
Ti=100;
Ts=400;
alpha=0.4;
del_t=0.1;
t_start=0;
t_end=1;

%% Modifying the BCs' and other variables according to non-dimensional setting

% x=L*x_nd
% T_nd=(T-Ts)/(Ti-Ts)
% t=(L^2/alpha)*t_nd

tnd_start=(alpha/L^2)*t_start;
tnd_end=(alpha/L^2)*t_end;
del_tnd=(alpha/L^2)*del_t;

% Number of time steps
Nt=round((tnd_end-tnd_start)/del_tnd);

% IC
Ti_nd=(Ti-Ts)/(Ti-Ts);

% BC
Ts_nd=(Ts-Ts)/(Ti-Ts);

%% Setting up of del_xnd for Crank Nicholson

% Since Crank Nicholson scheme is unconditionally stable, we can take any
% value of del_xnd. For uniformity, let us take del_xnd as 0.2
del_xnd=0.2;

% Number of grid points
Nx=round(1/del_xnd)+1;

%% Initial Conditions setup

T=zeros(Nt+1,Nx);
T_nd=zeros(1,Nx);
T_nd(1)=Ts_nd;
T_nd(Nx)=Ts_nd;
for j=2:(Nx-1)
    T_nd(j)=Ti_nd;
end
T(1,:)=T_nd*(Ti-Ts)+Ts;
time=0;

%% Time Marching of the solution

d=del_tnd/(del_xnd)^2;
a=zeros(Nx,1);
b=zeros(Nx,1);
c=zeros(Nx,1);
e=zeros(Nx,1);
for t=1:Nt
    T_nd_old=T_nd;
    
    for j=2:(Nx-1)
        a(j)=-d/2;
        b(j)=(1+d);
        c(j)=-d/2;
        e(j)=(1-d)*T_nd_old(j)+(d/2)*(T_nd_old(j-1)+T_nd_old(j+1));
    end
    
    e(2)=e(2)+(d/2)*T_nd(1);
    e(Nx-1)=e(Nx-1)+(d/2)*T_nd(Nx);
    
    e=TDMA(2,Nx-1,a,b,c,e);
    
    for j=2:(Nx-1)
        T_nd(j)=e(j);
    end
    
    % Updating BCs
    T_nd(1)=Ts_nd;
    T_nd(Nx)=Ts_nd;    
    
    T(t+1,:)=T_nd*(Ti-Ts)+Ts;
end

%% Plot of temperature profile

line={'*-','*:','*-.','*--','o-','o:','o-.','o--','s-','s:','s-.'};
for t=1:Nt+1
    plot(linspace(0,L,Nx),T(t,:),line{t})
    hold on
    legendInfo{t} = ['del_t = ' num2str((t-1)*0.1)];
end
xlabel('Length (ft)');
ylabel('Dimensional T (F)');
title('Crank Nicholson Scheme');
legend(legendInfo)

%% Analytical Solution at t=1hr

syms m;
x=linspace(0,L,Nx);
sum_fun=2*(Ti-Ts)*exp(-1*((m*pi)/L)^2*alpha*1)*(((1-(-1)^m)/(m*pi))*sin((m*pi*x)/L));
T_an=vpa(Ts+symsum(sum_fun,m,1,inf),4);

%% Error 

Err=(1/Nx)*sqrt(sum((T(Nt+1,:)-T_an).^2));