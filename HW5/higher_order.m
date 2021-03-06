%% CFD - HW5 - Mac-Cormack method

function u=higher_order(del_t)
    %% Number of nodes

    L=4;
    del_x=0.05;
    Nx=L/del_x+1;

    %% Courant Number

    T=6;
    Nt=T/del_t+1;
    cfl=del_t/del_x;

    %% Initialization and IC,BC

    u=zeros(Nt,Nx);

    % Initial Condition

    D=0:0.05:4;
    u(1,1:5)=1;
    u(1,6:26)=1.25-D(6:26);
    u(1,27:81)=0;

    %% MacCormack Method

    u_str=zeros(1,Nx);
    u_str(1,1)=1;
    u_str(1,Nx)=0;

    for i=2:Nt

        % Predictor Step
        for j=2:Nx-1
            u_str(1,j)=u(i-1,j)-(cfl/2)*((u(i-1,j+1))^2-(u(i-1,j))^2);
        end

        % Corrector Step
        for j=2:Nx-1
            u(i,j)=0.5*(u(i-1,j)+u_str(1,j)-(cfl/2)*((u_str(1,j))^2-(u_str(1,j-1))^2));
        end

        % Boundary Conditions
        u(i,1)=1;
        u(i,Nx)=0;
    end
    
    %% Visualization - at different times

    % plot(D,u(1,:),'r-*')
    % hold on
    % plot(D,u(1+round(2/del_t),:),'b-*')
    % hold on
    % plot(D,u(1+round(4/del_t),:),'m-*')
    % hold on
    % plot(D,u(Nt,:),'k-*')
    % grid on
    % ylim([-0.1 1.2])
    % legend('t=0s','t=2s','t=4s','t=6s');
    % title('MacCormack Method for CFL=1')
    % xlabel('x');
    % ylabel('u');

end
