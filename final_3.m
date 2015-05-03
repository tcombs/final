%%%%%%% 2 %%%%%%%%%
function f = final_3


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
plotSirVsTime(sir1,t1,'1');

%2 
plotSvsI(sir1,'1');

%3
%Jacobian function
function jac = jacobian(s,i,r,alpha,beta,lambda)
jac = [-lambda-beta*i -beta*s 0;
       beta*i beta*s-alpha-lambda 0;
       0 alpha -lambda;];
end

%reach equalibrium
[t3,equalibrium] = getSIR(0.999,0.001,0.0,1/3,1.05,1/60,100);
jacob = jacobian(equalibrium(1),equalibrium(2),equalibrium(3),1/3,1.05,1/60)
evs = eig(jacob)
%equalibrium

%4
function init = initialSIR
%pick random s
s =  rand;
range = 1-s;
i = range*rand;
r = 1 - (s + i);

init = [s i r ];
end


for n = 2:6
    %A
    n
    init = initialSIR
    [t,sir] = getSIR(init(1),init(2),init(3),1/3,1.05,1/60,200);
    plotSirVsTime(sir,t,int2str(n));
    plotSvsI(sir,int2str(n));
    

end

%plot all on one graph.

end