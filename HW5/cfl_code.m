%% CFD - HW5 - PLot for different CFL numbers

%% Plot for lower order

del_t=[0.03 0.05 0.06];
x=0:0.05:4;

for i=1:size(del_t,2)
    u1=lower_order(del_t(i));
    plot(x,u1(1+round(2/del_t(i)),:),'-*')
    hold on
end

grid on
title('Upwind scheme for different CFL numbers at t=2s');
xlabel('x');
ylabel('u');
xlim([1 3]);
ylim([-0.1 1.2])
legend('\Deltat=0.03, CFL=0.6','\Deltat=0.05, CFL=1.0','\Deltat=0.06, CFL=1.2')

%% Plot for higher order

figure;
for i=1:size(del_t,2)
    u2=higher_order(del_t(i));
    plot(x,u2(1+round(2/del_t(i)),:),'-*')
    hold on
end

grid on
title('MacCormack scheme for different CFL numbers at t=2s');
xlabel('x');
ylabel('u');
xlim([1 3]);
ylim([-0.1 1.2])
legend('\Deltat=0.03, CFL=0.6','\Deltat=0.05, CFL=1.0','\Deltat=0.06, CFL=1.2')