T_a=[400.0000  316.4250  264.9803  264.9803  316.4250  400.0000];
T_c=[400.0000  320.4278  271.2540  271.2540  320.4278  400.0000];
x=[0    0.4000    0.8000    1.2000    1.6000    2.0000];

T_b=[400.0000  341.3714  291.7214  258.5964  246.9768  258.5964  291.7214  341.3714  400.0000];
T_d=[400.0000  277.6340  418.0988  106.0080  424.5054  106.0080  418.0988  277.6340  400.0000];
y=[0    0.2500    0.5000    0.7500    1.0000    1.2500    1.5000    1.7500    2.0000];

plot(x,T_a,'r*--')
hold  on
plot(y,T_b,'bs-.')
title('CN for different del_x')
xlabel('Length (ft)')
ylabel('Temperature (F)')
legend('del_x=0.4','del_x=0.26')

figure
plot(x,T_c,'r*--')
hold  on
plot(y,T_d,'bs-.')
title('FTCS for different del_x')
xlabel('Length (ft)')
ylabel('Temperature (F)')
legend('del_x=0.4','del_x=0.26')