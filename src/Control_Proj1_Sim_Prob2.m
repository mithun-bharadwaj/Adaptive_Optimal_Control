%Control of Robotic Systems course- University of Maryland, College Park
%Code for the paper "Adaptive optimal control for continuous-time linear systems based on policy iteration"_Problem2

global A B K xn un
clc;
x_save=[];
t_save=[];

% System matrices used for simulation purpose
A_nom=[-0.0665 8 0 0; 0 -3.663 3.663 0; -6.86 0 -13.736 -13.736; 0.6 0 0 0];
 A=[-0.0665 11.5 0 0; 0 -2.5 2.5 0; -9.5 0 -13.736 -13.736; 0.6 0 0 0];

B=[0; 0; 13.736; 0];

x0=[0;0.1;0;0]; %initial conditions

[xn,un]=size(B);%size of B. un-column #, xn row #

% Set the weighting matrices for the cost function
Q=eye(4);
R=eye(1);

% Initialize the feedback gain matrix
K=zeros(un,xn);  % Only if A is Hurwitz, K can be set as zero.
P_old=zeros(xn);

T=.05;  %Duration of time for each integration



Dxx=[];XX=[];XU=[];  % Data matrices

X=[x0;kron(x0',x0')';kron(x0,zeros(un,1))]';



P=eye(xn)*10; % Initialize the previous cost matrix
it=0;            % Counter for iterations
p_save=[];       % Track the cost matrices in all the iterations
k_save=[];       % Track the feedback gain matrix in each iterations

[P0,K0]=care(A,B,Q); % Calculate the ideal solution for comparion purpose
k_save=[norm(K-K0)];

while norm(P-P_old)>10^-3   % Stopping criterion for learning
    for i=1:20
    % Simulation the system and at the same time collect online info.
    [t,X]=ode45(@mysys, [it+(i-1)*T,it+i*T],X(end,:));
   
    %Append new data to the data matrices
    Dxx=[Dxx;kron(X(end,1:xn),X(end,1:xn))-kron(X(1,1:xn),X(1,1:xn))];
    XX=[XX;X(end,xn+1:xn+xn^2)-X(1,xn+1:xn+xn^2)];
    XU=[XU;X(end,xn+xn^2+1:end)-X(1,xn+xn^2+1:end)];

    % Keep track of the system trajectories
    x_save=[x_save;X];
    t_save=[t_save;t];
    end
    
   
    P_old=P;                        % Update the previous cost matrix
    
     
    QK=Q+K'*R*K;                    % Update the Qk matrix
    X2=XX*kron(eye(xn),K');         %
    X1=[Dxx,-X2-XU];                % Left-hand side of the key equation
    Y=-XX*QK(:);                    % Right-hand side of the key equation
    
    pp = pinv(X1)*Y;
	
	P = reshape(pp(1:xn*xn), [xn, xn]);
	P = (P +P')/2;
      
     figure(4)
    
    scatter(it,P_old(1,1),'o','filled','b')
    hold on
    scatter(it,P_old(2,3),'+','g','Linewidth',2)
    hold on
    scatter(it,P_old(2,4),'o','r','Linewidth',2)
    hold on
    scatter(it,P_old(3,3),'s','c','Linewidth',2)
    hold on
 
    p_save=[p_save,norm(P-P0)];     % Keep track of the cost matrix
    
    K=inv(R)*transpose(B)*P;          % Get the improved gain matrix
    k_save=[k_save,norm(K-K0)];     % Keep track of the control gains
  
    

    it=it+1     % Update and display the # of iters
end
fprintf('\nNumber of iterations to obtain optimal control solution %d\n',it)

% Plot the trajectories

figure(4)

scatter(it,P0(1,1),'*','b')
hold on
scatter(it,P0(2,3),'*','g')
hold on
scatter(it,P0(2,4),'*','r')
hold on
scatter(it,P0(3,3),'*','c')
hold on
axis([0 it 0 3.5])
xlabel('Time (s)')
title('P Matrix parameters')
legend('P(1,1)','P(2,3)','P(2,4)','P(4,4)','P(1,1)-Optimal','P(2,3)-Optimal','P(2,4)-Optimal','P(4,4)-Optimal')



    
figure(1)
plot([0:length(p_save)-1],p_save,'o',[0:length(p_save)-1],p_save)
legend('||P_i-P^*||')
xlabel('Number of iterations')

figure(2)
plot([0:length(k_save)-1],k_save,'^',[0:length(k_save)-1],k_save)
legend('||K_i-K^*||')
xlabel('Number of iterations')


% Post-learning simulation
[tt,xx]=ode23(@mysys,[t(end) 20],X(end,:)');

% Keep track of the post-learning trajectories
t_final=[t_save;tt];
x_final=[x_save;xx];


figure(3)

plot(t_final,x_final(:,1:4),'Linewidth',2)
axis([0,2,-0.2,0.2])    %axis can be changed to get the complete trejectory
legend('x_1','x_2','x_3','x_4')
xlabel('Time (s)')
title('System States')


% The following nested function gives the dynamics of the sytem.
   function dX=mysys(t,X)
       global A B K xn
        x=X(1:xn);

        u=-K*x;
        dx=A*x+B*u;
        dxx=kron(x',x')';
        dux=kron(x',u')';
        dX=[dx;dxx;dux];
	end

