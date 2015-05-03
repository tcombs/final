%%%%%%% 2 %%%%%%%%%
function f = final_2


function [tg,fg] = getSIR(s_0,i_0,r_0,alpha,beta,lambda,maxTime)

function ff = odeSIR(s,i,r,alpha,beta,lambda,t)
ff = [lambda - lambda*s - beta*s*i ;beta*s*i - alpha*i - lambda*i ; alpha*i - lambda*r];
end

sir_0 = [s_0 i_0 r_0];


[t,sir] = ode45(@(t,sir) odeSIR(sir(1),sir(2),sir(3),alpha,beta,t),[0,maxTime],sir_0);

%return sir
fg = sir;
tg = t;
end

function fh = plotSirVsTime(sir,t,id)
fig=figure;
plot(t,sir(:,1:1),t,sir(:,2:2),t,sir(:,3:3));
title(strcat('SIR vs T ',id));
legend('S','I','R');
xlabel('time');
saveas(fig,strcat('SIRvTIME_',id),'png');
end

function fi = plotSvsI(sir,id)
fig2 = figure;
plot(sir(:,1:1),sir(:,2:2));
title(strcat('S vs I ',id));
xlabel('s(t)');
ylabel('i(t)');
saveas(fig2,strcat('SvI_',id), 'png');
end

%1 
[t1,sir1] = getSIR(0.999,0.001,0.0,1/3,1.05,1/60,200);
sir1(2:2,:)
plotSirVsTime(sir1,t1,'1');

%2 
plotSvsI(sir1,'1');

%Jacobian function
function jac = jacobian(s,i,r,alpha,beta,lambda)
jac = [-lambda-beta*i -beta*s 0;
       beta*i beta*s-alpha-lambda 0;
       0 alpha -lambda;];
end

%reach equalibrium
[t3,sir3] = getSIR(0.999,0.001,0.0,1/3,1.05,1/60,1000);


end