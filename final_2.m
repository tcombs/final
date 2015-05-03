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



    init2 = initialSIR
    [t2,sir2] = getSIR(init2(1),init2(2),init2(3),1/3,1.05,1/60,200);
    plotSvsI(sir2,'2');

    init3 = initialSIR
    [t3,sir3] = getSIR(init3(1),init3(2),init3(3),1/3,1.05,1/60,200);
    plotSvsI(sir3,'3');

    init4 = initialSIR
    [t4,sir4] = getSIR(init4(1),init4(2),init4(3),1/3,1.05,1/60,200);
    plotSvsI(sir4,'4');

    init5 = initialSIR
    [t5,sir5] = getSIR(init5(1),init5(2),init5(3),1/3,1.05,1/60,200);
    plotSvsI(sir5,'5');

    init6 = initialSIR
    [t6,sir6] = getSIR(init6(1),init6(2),init6(3),0.1,0.3,3,200);
    plotSvsI(sir6,'6');


%plot all on one graph.
fig2 = figure;
plot(sir1(:,1:1),sir1(:,2:2),sir2(:,1:1),sir2(:,2:2),sir3(:,1:1),sir3(:,2:2),sir4(:,1:1),sir4(:,2:2),sir5(:,1:1),sir5(:,2:2),sir6(:,1:1),sir6(:,2:2));
legend('Start point 1','Start point 2','Start point 3','Start point 4','Start point 5','Start point 6');
title('S vs I ALL 2');
xlabel('s(t)');
ylabel('i(t)');
saveas(fig2,'SvI_All', 'png');
end