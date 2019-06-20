function [T_an]=ansol(Nx,Ny,L,W,T1,T3)
    syms m;
    del_x=L/(Nx-1);
    del_y=W/(Ny-1);
    
    x=0:del_x:L;
    y=W:-del_y:0;
    [X,Y]=meshgrid(x,y);
    Ta=((1-cos(m*pi))/(m*pi))*(sinh((m*pi*(W-Y))/L)).*sin((m*pi*X)/L)/(sinh((m*pi*W)/L));
    Tb=((1-cos(m*pi))/(m*pi))*(sinh((m*pi*Y)/L)).*sin((m*pi*X)/L)/(sinh((m*pi*W)/L));
    T_an=vpa(2*T1*symsum(Ta,m,1,100)+2*T3*symsum(Tb,m,1,100),4);
end