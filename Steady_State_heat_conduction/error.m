%% Error 

%% ADI

err_y=[0.08366 0.05479 0.03848 0.02825 0.02165 0.01748];
dely=[0.075 0.06 0.05 0.0429 .0375 0.0333];
m_y=(log(err_y(1))-log(err_y(2)))/(log(dely(1))-log(dely(2)));

err_x=[0.20763 0.12813 0.08761 0.06464 0.04902 0.03999];
delx=[0.1 0.08 0.066666667 0.057142857 0.05 0.044444444];
m_x=(log(err_x(1))-log(err_x(2)))/(log(delx(1))-log(delx(2)));

loglog(dely,err_y)
hold on
loglog(delx,err_x)
title('Loglog plot for Error vs Step size - ADI')
xlabel('Step Size')
ylabel('Error')
legend('del_y','del_x')
grid on

%% SOR

sor_y=[0.060598 0.038573 0.026526 0.0220113 0.018997];
sor_x=[0.283339 0.18636072 0.1322921 0.101574 0.083322 0.0683109];

ms_y=(log(sor_y(1))-log(sor_y(2)))/(log(dely(1))-log(dely(2)));
ms_x=(log(sor_x(1))-log(sor_x(2)))/(log(delx(1))-log(delx(2)));

figure;
loglog(dely(1:5),sor_y)
hold on
loglog(delx,sor_x)
title('Loglog plot for Error vs Step size - SOR')
xlabel('Step Size')
ylabel('Error')
legend('del_y','del_x')
grid on





