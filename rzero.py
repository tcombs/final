import random
R0 = .5
beta = 0
delta = 0
epsilon = 0
gamma = 0
while R0 < 1:
	beta = random.random()
	delta = random.random()
	epsilon = random.random()
	gamma = random.random()
	R0 = (beta*epsilon)/((delta+epsilon)*(delta+gamma))

print R0
print beta
print delta
print epsilon
print gamma