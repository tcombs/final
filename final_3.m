%%%%%%% 3 %%%%%%%%%
function f = final_3

function rz = rzero(beta,delta,epsilon,gamma)
	rz = (beta*epsilon)/((delta+epsilon)*(delta+gamma))
end


function [tg,fg] = getSEIR(s_0,e_0,i_0,r_0,beta,delta,epsilon,gamma,maxTime)

	function ff = odeSEIR(s,e,i,r,beta,delta,epsilon,gamma,t)
		ff = [-delta];
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

end