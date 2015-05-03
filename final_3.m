%%%%%%% 3 %%%%%%%%%
function f = final_3

function rz = rzero(beta,delta,epsilon,gamma)
	rz = (beta*epsilon)/((delta+epsilon)*(delta+gamma))
end


function [tg,fg] = getSEIR(s_0,e_0,i_0,r_0,beta,delta,epsilon,gamma,maxTime)

	function ff = odeSEIR(s,e,i,r,beta,delta,epsilon,gamma,t)
		ff = [-delta*s-beta*s*i+delta ; -(delta+epsilon)*e+beta*s*i ; -(delta+gamma)*i+epsilon*e ; gamma*i - delta * r];
	end
	seir_0 = [s_0 e_0 i_0 r_0];

	[t,seir] = ode45(@(t,seir) odeSEIR(seir(1),seir(2),seir(3),seir(4),beta,delta,epsilon,gamma,t),[0,maxTime],seir_0);

	%return seir
	fg = seir;
	tg = t;
end

betaA = 0.7764
deltaA = 0.9426
epsilonA = 0.9425
gammaA = 0.9455
R0A = rzero(betaA,deltaA,epsilonA,gammaA)



betaB = 0.6049
deltaB = 0.1497
epsilonB = 0.6151
gammaB = 0.1895
R0B = rzero(betaB,deltaB,epsilonB,gammaB)

%get initial conditions
function init = initialSEIR
%pick random s
s =  rand;
range = 1-s;
e = range*rand;
range = range - e;
i = range*rand;
r = 1 - (s + e + i);

init = [s e i r];
end

init = initialSEIR;
[tA,seirA] = getSEIR(init(1),init(2),init(3),init(4),betaA,deltaA,epsilonA,gammaA,100);
[tB,seirB] = getSEIR(init(1),init(2),init(3),init(4),betaB,deltaB,epsilonB,gammaB,100);



%Jacobian
function jac = jacobian(s,e,i,r,beta,delta,epsilon,gamma)
	jac = [-delta-beta*i, 0, -beta*s, 0 ;
		   beta*i,-(delta+epsilon), beta*s 0;
		   0, epsilon, -(delta+gamma), 0;
		   0, 0, gamma, delta];
end

jacobA = jacobian(seirA(1),seirA(2),seirA(3),seirA(4),betaA,deltaA,epsilonA,gammaA)
eigA = eig(jacobA)

jacobB = jacobian(seirB(1),seirB(2),seirB(3),seirB(4),betaB,deltaB,epsilonB,gammaB)
eigB = eig(jacobB)

%plot 3d phase and seir vs time

function fh = plotSeirVsTime(seir,t,id)
	fig=figure;
	plot(t,seir(:,1:1),t,seir(:,2:2),t,seir(:,3:3),t,seir(:,4:4));
	title(strcat('SEIR vs T ',id));
	legend('S','E','I','R');
	xlabel('time');
	saveas(fig,strcat('SEIRvTIME_',id),'png');
end

%3d phase plot
function fi = plot3Dphase(seir,id)
	fig2 = figure;
	plot3(seir(:,1:1),seir(:,2:2),seir(:,3:3));
	title(strcat('S vs E vs I ',id));
	xlabel('s(t)');
	ylabel('e(t)');
	zlabel('i(t)');
	saveas(fig2,strcat('3-Dphase_',id), 'png');
end

plotSeirVsTime(seirA,tA,'1a');
plotSeirVsTime(seirB,tB,'1b');
plot3Dphase(seirA,'1a');
plot3Dphase(seirB,'1b');

%repeat 5 more times


init2 = initialSEIR
[tA2,seirA2] = getSEIR(init2(1),init2(2),init2(3),init2(4),betaA,deltaA,epsilonA,gammaA,100);
[tB2,seirB2] = getSEIR(init2(1),init2(2),init2(3),init2(4),betaB,deltaB,epsilonB,gammaB,100);

init3 = initialSEIR
[tA3,seirA3] = getSEIR(init3(1),init3(2),init3(3),init3(4),betaA,deltaA,epsilonA,gammaA,100);
[tB3,seirB3] = getSEIR(init3(1),init3(2),init3(3),init3(4),betaB,deltaB,epsilonB,gammaB,100);

init4 = initialSEIR
[tA4,seirA4] = getSEIR(init4(1),init4(2),init4(3),init4(4),betaA,deltaA,epsilonA,gammaA,100);
[tB4,seirB4] = getSEIR(init4(1),init4(2),init4(3),init4(4),betaB,deltaB,epsilonB,gammaB,100);

init5 = initialSEIR
[tA5,seirA5] = getSEIR(init5(1),init5(2),init5(3),init5(4),betaA,deltaA,epsilonA,gammaA,100);
[tB5,seirB5] = getSEIR(init5(1),init5(2),init5(3),init5(4),betaB,deltaB,epsilonB,gammaB,100);

init6 = initialSEIR
[tA6,seirA6] = getSEIR(init6(1),init6(2),init6(3),init6(4),betaA,deltaA,epsilonA,gammaA,100);
[tB6,seirB6] = getSEIR(init6(1),init6(2),init6(3),init6(4),betaB,deltaB,epsilonB,gammaB,100);


plotSeirVsTime(seirA2,tA2,'2a');
plotSeirVsTime(seirB2,tB2,'2b');
plot3Dphase(seirA2,'2a');
plot3Dphase(seirB2,'2b');

plotSeirVsTime(seirA3,tA3,'3a');
plotSeirVsTime(seirB3,tB3,'3b');
plot3Dphase(seirA3,'3a');
plot3Dphase(seirB3,'3b');

plotSeirVsTime(seirA4,tA4,'4a');
plotSeirVsTime(seirB4,tB4,'4b');
plot3Dphase(seirA4,'4a');
plot3Dphase(seirB4,'4b');

plotSeirVsTime(seirA5,tA5,'5a');
plotSeirVsTime(seirB5,tB5,'5b');
plot3Dphase(seirA5,'5a');
plot3Dphase(seirB5,'5b');

plotSeirVsTime(seirA6,tA6,'6a');
plotSeirVsTime(seirB6,tB6,'6b');
plot3Dphase(seirA6,'6a');
plot3Dphase(seirB6,'6b');

end