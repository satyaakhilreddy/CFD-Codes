%% Error Plots

Errx_FTCS=[4.0099 2.286 1.983 1.8422 1.7389];
Errx_CN=[2.535 0.4928 0.1836 0.0908 0.0578];
dx=[1 0.666666667 0.5 0.4 0.333333333];

loglog(dx,Errx_FTCS,'s--')
hold on
loglog(dx,Errx_CN,'*:')
title('Loglog plot for Error vs del_x');
xlabel('del_x');
ylabel('Error');
legend('FTCS','Crank-Nicholson');
mx_FTCS=(log(Errx_FTCS(3))-log(Errx_FTCS(1)))/(log(dx(3))-log(dx(1)));
mx_CN=(log(Errx_CN(2))-log(Errx_CN(1)))/(log(dx(2))-log(dx(1)));

Errt_CN=[0.08337 0.0865 0.0908 0.09677 0.1051 0.1176 0.1372 0.1704 0.23215 0.3666 1.5314];
Errt_FTCS=[1.5348 1.674 1.8422 2.049 2.312 2.654 3.118 3.787 6.1205 38.501 119.1535];
dt=[0.083333333 0.090909091 0.1 0.111111111 0.125 0.142857143 0.166666667 0.2 0.25 0.333333333 0.5];

figure
loglog(dt,Errt_FTCS,'s--')
hold on
loglog(dt,Errt_CN,'*:')
xlim([10^-(1.15) 1]);
title('Loglog plot for Error vs del_t');
xlabel('del_t');
ylabel('Error');
legend('FTCS','Crank-Nicholson');



