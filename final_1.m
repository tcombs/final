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



%a, b Figure out the sir models
for n = 2:6
	%A
	init = initialSIR
    [ta,sira] = getSIR(init(1),init(2),init(3),0.1,0.3,3,200);
    plotSirVsTime(sira,ta,strcat(int2str(n),'A'));
    plotSvsI(sira,strcat(int2str(n),'A'));
    n
    size(sira)
    allSIRSA(:,:,n) = sira

    %B
    [tb,sirb] = getSIR(init(1),init(2),init(3),0.1,0.3/4,.75,1000);
    plotSirVsTime(sirb,tb,strcat(int2str(n),'B'));
    plotSvsI(sirb,strcat(int2str(n),'B'));
    allSIRSB = [allSIRSB;sirb];
end

%plot all on one graph.

 




end


