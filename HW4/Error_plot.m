%% Error Plots

Errx_FTCS=[2.517 0.4684 0.1634 0.06898 0.03538];
Errx_CN=[2.516 0.4666 0.1618 0.06745 0.0339];
dx=[1 0.666666667 0.5 0.4 0.333333333];

loglog(dx,Errx_FTCS,'s--')
hold on
loglog(dx,Errx_CN,'*:')
title('Loglog plot for Error vs del_x');
xlabel('del_x');
ylabel('Error');
legend('FTCS','Crank-Nicholson');

mx_FTCS=(log(Errx_FTCS(4))-log(Errx_FTCS(3)))/(log(dx(4))-log(dx(3)));
mx_CN=(log(Errx_CN(4))-log(Errx_CN(3)))/(log(dx(4))-log(dx(3)));

Errt_CN=[2.859 0.9761 0.2974 0.09625 0.05315 0.04071 0.03321 0.02789];
Errt_FTCS=[38.5 6.121 3.787 3.119 2.654 2.312 2.05 1.842];
dt=[0.333333333 0.25 0.2 0.166666667 0.142857143 0.125 0.111111111 0.1];

figure
loglog(dt,Errt_FTCS,'s--')
hold on
loglog(dt,Errt_CN,'*:')
xlim([10^-(1.15) 1]);
title('Loglog plot for Error vs del_t');
xlabel('del_t');
ylabel('Error');
legend('FTCS','Crank-Nicholson');

mt_FTCS=(log(Errt_FTCS(4))-log(Errt_FTCS(3)))/(log(dt(4))-log(dt(3)));
mt_CN=(log(Errt_CN(2))-log(Errt_CN(1)))/(log(dt(2))-log(dt(1)));


