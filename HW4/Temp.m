clc
clear all;
%% Parameters

L=2;
Ti=100;
Ts=400;
alpha=0.4;
del_t=0.0001;
t_start=0;
t_end=1;
%% Analytical Solution at t=1hr

syms m;

%% Error Crank-Nicholson del_x

i=0;
for Nt=3:10
    i=i+1;
    [T,Nx]=Test_code(Nt);
    x=linspace(0,L,Nx);
    sum_fun=2*(Ti-Ts)*exp(-1*((m*pi)/L)^2*alpha*1)*(((1-(-1)^m)/(m*pi))*sin((m*pi*x)/L));
    T_an=vpa(Ts+symsum(sum_fun,m,1,inf),4);
    Er(i)=(1/Nx)*sqrt(sum((T(Nt+1,:)-T_an).^2));
    dt(i)=(t_end-t_start)/Nt;
end