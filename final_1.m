%%%%%%% 1 %%%%%%%%%
function f = final_1


function [tg,fg] = getSIR(s_0,i_0,r_0,alpha,beta,R_0,maxTime)

function ff = odeSIR(s,i,r,alpha,beta,t)
ff = [-1*beta*s*i;beta*s*i - alpha * i;alpha * i];
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

%1 a,b
[t1a,sir1a] = getSIR(0.999,0.001,0.0,0.1,0.3,3,200);
[t1b,sir1b] = getSIR(0.999,0.001,0.0,0.1,0.3/4,.75,1000);
plotSirVsTime(sir1a,t1a,'1A');
plotSirVsTime(sir1b,t1b,'1B');

%2 a,b
plotSvsI(sir1a,'1A');
plotSvsI(sir1b,'1B');

%3 a,b
%Jacobian function
function jac = jacobian(s,i,r,alpha,beta)
jac = [-beta*i -beta*s 0;
	   beta*i beta*s-alpha 0;
	   0 alpha 0;];
end

% run sir until we reach equalibrium A
[t3a,sir3a] = getSIR(0.999,0.001,0.0,0.1,0.3,3,1000);
%extract last row of equalibrium
equalibriumA = sir3a(end:end,:)
jacobA = jacobian(equalibriumA(1),equalibriumA(2),equalibriumA(3),.1,.3)
evsA = eig(jacobA)



% run sir until we reach equalibrium B
[t3b,sir3b] = getSIR(0.999,0.001,0.0,0.1,0.3/4,.75,1000);
%extract last row of equalibrium
equalibriumB = sir3a(end:end,:)
jacobB = jacobian(equalibriumB(1),equalibriumB(2),equalibriumB(3),.1,.3/4)
evsB = eig(jacobB)


function init = initialSIR
%pick random s
s =  rand;
range = 1-s;
i = range*rand;
r = 1 - (s + i);

init = [s i r ];
end

%4 Run 5 times with different init conditions
init = initialSIR;



%a, b Figure out the sir models 5 more times

	%A
	init2 = initialSIR
    [ta2,sir2a] = getSIR(init2(1),init2(2),init2(3),0.1,0.3,3,200);
    plotSvsI(sir2a,'2A');
    plotSirVsTime(sir2a,ta2,'2A');

    init3 = initialSIR
    [ta3,sir3a] = getSIR(init3(1),init3(2),init3(3),0.1,0.3,3,200);
    plotSvsI(sir3a,'3A');
    plotSirVsTime(sir3a,ta3,'3A');

    init4 = initialSIR
    [ta4,sir4a] = getSIR(init4(1),init4(2),init4(3),0.1,0.3,3,200);
    plotSvsI(sir4a,'4A');
    plotSirVsTime(sir4a,ta4,'4A');

    init5 = initialSIR
    [ta5,sir5a] = getSIR(init5(1),init5(2),init5(3),0.1,0.3,3,200);
    plotSvsI(sir5a,'5A');
    plotSirVsTime(sir5a,ta5,'5A');

    init6 = initialSIR
    [ta6,sir6a] = getSIR(init6(1),init6(2),init6(3),0.1,0.3,3,200);
    plotSvsI(sir6a,'6A');
    plotSirVsTime(sir6a,ta6,'6A');

    %B
    %[tb,sirb] = getSIR(init(1),init(2),init(3),0.1,0.3/4,.75,1000);
    %plotSirVsTime(sirb,tb,strcat(int2str(n),'B'));
    %plotSvsI(sirb,strcat(int2str(n),'B'));
    %allSIRSB = [allSIRSB;sirb];


    [tb2,sir2b] = getSIR(init2(1),init2(2),init2(3),0.1,0.3/4,.75,1000);
    plotSvsI(sir2b,'2b');
    plotSirVsTime(sir2b,tb2,'2B');

    [tb3,sir3b] = getSIR(init3(1),init3(2),init3(3),0.1,0.3/4,.75,1000);
    plotSvsI(sir3b,'3b');
    plotSirVsTime(sir3b,tb3,'3B');

    [tb4,sir4b] = getSIR(init4(1),init4(2),init4(3),0.1,0.3/4,.75,1000);
    plotSvsI(sir4b,'4b');
    plotSirVsTime(sir4b,tb4,'4B');

    [tb5,sir5b] = getSIR(init5(1),init5(2),init5(3),0.1,0.3/4,.75,1000);
    plotSvsI(sir5b,'5b');
    plotSirVsTime(sir5b,tb5,'5B');

    [tb6,sir6b] = getSIR(init6(1),init6(2),init6(3),0.1,0.3/4,.75,1000);
    plotSvsI(sir6b,'6b');
    plotSirVsTime(sir6b,tb6,'6B');


%plot all in one graph

 

fig2 = figure;
plot(sir1a(:,1:1),sir1a(:,2:2),sir2a(:,1:1),sir2a(:,2:2),sir3a(:,1:1),sir3a(:,2:2),sir4a(:,1:1),sir4a(:,2:2),sir5a(:,1:1),sir5a(:,2:2),sir6a(:,1:1),sir6a(:,2:2));
legend('Start point 1','Start point 2','Start point 3','Start point 4','Start point 5','Start point 6');
title('S vs I ALL A');
xlabel('s(t)');
ylabel('i(t)');
saveas(fig2,'SvI_All_A', 'png');


fig3 = figure;
plot(sir1b(:,1:1),sir1b(:,2:2),sir2b(:,1:1),sir2b(:,2:2),sir3b(:,1:1),sir3b(:,2:2),sir4b(:,1:1),sir4b(:,2:2),sir5b(:,1:1),sir5b(:,2:2),sir6b(:,1:1),sir6b(:,2:2));
legend('Start point 1','Start point 2','Start point 3','Start point 4','Start point 5','Start point 6');
title('S vs I ALL B');
xlabel('s(t)');
ylabel('i(t)');
saveas(fig3,'SvI_All_B', 'png');


end


