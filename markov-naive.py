import sys
import os
import matplotlib.pyplot as plt

class Markov_naive():

	parameters = { 
					"beta_1"  : 0.1,    # infection rate/day of symptomatic
					"beta_2"  : 0.3,   # infection rate/day of asymptomatic
					"beta_4"  : 0.01,  # infection rate/day of isolated
					"gamma_1" : 0.1,   # isolation percent/day of asymptomatic
					"gamma_2" : 0.8,    # isolation percent/day of symptomatic
					"gamma_start" : 10,  # first day to start isolation and medical treatment
					"alpha_1_start" : 14,  # fisrt day to have self-cured
					"alpha_1" : 0.07,       # rate/day of self-cured
					"alpha_2_start" : 12,  # fisrt day to have medical-cured after isolation
					"alpha_2" : 0.1,       # rate/day of medical-cured
					"alpha_3_start" : 10,  # first day to have death
					"alpha_3" : 0.02,      # rate/day of death
					"beta_3_start"  : 3,   # first day to have symptomatic
					"beta_3"  : 0.3        # rate/day of asymptomatic to symptomatic
				}

	initial_population = {
							"N" : 10000, # total
							"A" : 9900,  # common susceptible
							"B" : 100,     # asymptomatic
							"C" : 0,     # symptomatic
							"D" : 0,     # isolated / medical interventional
							"E" : 0,     # recovered / have antibodies
							"F" : 0      # dead
						}


	def __init__(self, parameters=None, initial_population=None):

		if parameters is not None and initial_population is not None:
			self.check_population(initial_population)
			self.check_parameters(parameters)
			if initial_population["B"] == 0 and initial_population["C"] == 0:
				raise ValueError("There should be B or C people in initial population !")
			if initial_population["A"] == 0:
				raise ValueError("Initial population should have susceptible (A) !")
			self.init_population = initial_population.copy()
			self.param = parameters.copy()
		else:
			self.init_population = Markov_naive.initial_population
			self.param = Markov_naive.parameters

		self.initial_population["NewCase"] = 0
		
		self.history = {0:self.init_population}   # {t : population, ...}


	def check_population(self, population):

		t1 = population["A"] + population["B"] + population["C"] + population["D"] + population["E"]
		t2 = population["N"]
		if abs(t2-t1) > 1e-3:
			raise RuntimeError("The population is not self consistent !")


	def check_parameters(self, parameters):

		if parameters["alpha_1"] + parameters["beta_3"] + parameters["gamma_1"] > 1:
			raise ValueError("ERR 1 !")
		if parameters["alpha_2"] + parameters["alpha_3"] > 1:
			raise ValueError("ERR 2 !")
		if parameters["gamma_2"] > 1:
			raise ValueError("ERR 3 !")


	def __gaussian(self, mu, sigma, t):
		return 1/sqrt(2*np.pi*sigma) * np.exp(-(t-mu)**2/sigma**2)

	'''
	def alpha_1(self, t):
		return self.__gaussian(self.param["alpha_1_mu"], self.param["alpha_1_sg"], t)


	def alpha_2(self, t):
		return self.param["alpha_2_hi"] * self.__gaussian(self.param["alpha_2_mu"], self.param["alpha_2_sg"], t)


	def alpha_3(self, t):
		return (1-self.param["alpha_2_hi"]) * self.__gaussian(self.param["alpha_3_mu"], self.param["alpha_3_sg"], t)
	

	def beta_3(self, t):
		return self.__gaussian(self.param["beta_3_mu"], self.param["beta_3_sg"], t)
	'''

	def running(self, days=100):

		for tmp in range(days):

			t = tmp + 1
			pl = self.history[t-1].copy()

			#beta_3 = self.beta_3(t)
			#alpha_1 = self.alpha_1(t)
			#alpha_2 = self.alpha_2(t)
			#alpha_3 = self.alpha_3(t)
			if t >= self.param["gamma_start"]:
				gamma_1 = self.param["gamma_1"]
				gamma_2 = self.param["gamma_2"]
			else:
				gamma_1 = 0
				gamma_2 = 0
			if t >= self.param["alpha_1_start"]:
				alpha_1 = self.param["alpha_1"]
			else:
				alpha_1 = 0
			if t >= self.param["alpha_2_start"] + self.param["gamma_start"]:
				alpha_2 = self.param["alpha_2"]
			else:
				alpha_2 = 0
			if t >= self.param["alpha_3_start"]:
				alpha_3 = self.param["alpha_3"]
			else:
				alpha_3 = 0
			if t >= self.param["beta_3_start"]:
				beta_3 = self.param["beta_3"]
			else:
				beta_3 = 0

			for i in range(100):

				dt = 1/100

				dA = - self.param["beta_1"] * pl["C"] * pl["A"] / pl["N"] \
					- self.param["beta_2"] * pl["B"] * pl["A"] / pl["N"] \
					- self.param["beta_4"] * pl["D"] * pl["A"] / pl["N"]

				dB = - dA - alpha_1 * pl["B"] \
					- gamma_1 * pl["B"] - beta_3 * pl["B"]

				dC = beta_3 * pl["B"] - gamma_2 * pl["C"]

				dD = gamma_1 * pl["B"] + gamma_2 * pl["C"] \
					- alpha_2 * pl["D"] - alpha_3 * pl["D"]

				dE = alpha_1 * pl["B"] + alpha_2 * pl["D"]

				dF = alpha_3 * pl["D"]

				pl["A"] += dA * dt
				pl["B"] += dB * dt
				pl["C"] += dC * dt
				pl["D"] += dD * dt
				pl["E"] += dE * dt
				pl["F"] += dF * dt
				pl["NewCase"] = -dA * dt

			self.history[t] = pl

			if pl["A"] < 0:
				pl["B"] -= (0 - pl["A"])
				pl["A"] = 0


if __name__ == '__main__':
	
	model = Markov_naive()
	model.running(100)

	plh = model.history

	t = list(plh.keys())
	A = []
	BC = []
	newBC = []
	D = []
	E = []
	F = []
	for k,v in plh.items():
		A.append(v["A"])
		BC.append(v["B"]+v["C"]+v["D"])
		newBC.append(v["NewCase"])
		D.append(v["D"])
		E.append(v["E"])
		F.append(v["F"])

	plt.plot(t, A, 'b')
	plt.plot(t, BC, 'r')
	plt.plot(t, D, 'y')
	plt.plot(t, E, 'g')
	plt.plot(t, F, 'k')
	#plt.plot(t, newBC, 'm')
	plt.legend(["susceptible", "infected", "quarantined", "recovered", "dead"])
	#plt.legend(["isolated","new-case"])
	plt.xlabel("Days")
	plt.ylabel("Counts")
	plt.show()




