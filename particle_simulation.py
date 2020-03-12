import sys
import os
import numpy as np
import particle_simulation_config as config
from particle_simulation_config import STAT
from particle_simulation_config import PARA



def A0_to_next(stat_list, stat_time):

	fate = np.random.choice([STAT.A1, STAT.R, STAT.Q], p=[PARA.A0_to_A1, PARA.A0_to_R, PARA.A0_to_Q])

	if fate == STAT.A1:
		time_interval = round(np.random.normal(PARA.mu_A0_to_A1, PARA.sg_A0_to_A1))
	elif fate == STAT.R:
		time_interval = round(np.random.normal(PARA.mu_A0_to_R, PARA.sg_A0_to_R))
	elif fate == STAT.Q:
		time_interval = round(np.random.normal(PARA.mu_A0_to_Q, PARA.sg_A0_to_Q))
	else:
		raise RuntimeError("A0_to_next error")

	time_interval = max(0, time_interval)

	return time_interval, fate


def A1_to_next(stat_list, stat_time):

	assert len(stat_time) > 1, "A1_to_next error !"
	time_A0_to_now = stat_time[-1] - stat_time[0]
	fate = np.random.choice([STAT.S, STAT.R, STAT.Q], p=[PARA.A1_to_S, PARA.A1_to_R, PARA.A1_to_Q])

	if fate == STAT.S:
		time_interval = round(np.random.normal(PARA.mu_A1_to_S, PARA.sg_A1_to_S))
	elif fate == STAT.Q:
		time_interval = round(np.random.normal(PARA.mu_A1_to_Q, PARA.sg_A1_to_Q))
	elif fate == STAT.R:
		time_interval = -1
		while time_interval < 0:
			time_interval = round(np.random.normal(PARA.mu_A0_to_R, PARA.sg_A0_to_R)) - time_A0_to_now
	else:
		raise RuntimeError("A1_to_next error")

	time_interval = max(0, time_interval)

	return time_interval, fate


def S_to_next(stat_list, stat_time):

	assert len(stat_time) > 1, "S_to_next error !"
	time_A0_to_now = stat_time[-1] - stat_time[0]
	fate = np.random.choice([STAT.R, STAT.Q, STAT.D], p=[PARA.S_to_R, PARA.S_to_Q, PARA.S_to_D])

	if fate == STAT.Q:
		time_interval = round(np.random.normal(PARA.mu_S_to_Q, PARA.sg_S_to_Q))
	elif fate == STAT.R:
		time_interval = -1
		while time_interval < 0:
			time_interval = round(np.random.normal(PARA.mu_A0_to_R, PARA.sg_A0_to_R)) - time_A0_to_now
	elif fate == STAT.D:
		time_interval = -1
		while time_interval < 0:
			time_interval = round(np.random.normal(PARA.mu_A0_to_D, PARA.sg_A0_to_D)) - time_A0_to_now
	else:
		raise RuntimeError("S_to_next error")

	time_interval = max(0, time_interval)

	return time_interval, fate


def Q_to_next(stat_list, stat_time):

	assert len(stat_time) > 1, "Q_to_next error !"
	time_A0_to_now = stat_time[-1] - stat_time[0]
	fate = np.random.choice([STAT.R, STAT.D], p=[PARA.Q_to_R, PARA.Q_to_D])

	if fate == STAT.D:
		time_interval = -1
		while time_interval < 0:
			time_interval = round(np.random.normal(PARA.mu_A0_to_D, PARA.sg_A0_to_D)) - time_A0_to_now
		time_interval = round(time_interval * PARA.Q_longer_D)
	elif fate == STAT.R:
		time_interval = -1
		while time_interval < 0:
			time_interval = round(np.random.normal(PARA.mu_A0_to_R, PARA.sg_A0_to_R)) - time_A0_to_now
		time_interval = round(time_interval * PARA.Q_shorter_R)
	else:
		raise RuntimeError("Q_to_next error")

	time_interval = max(0, time_interval)

	return time_interval, fate


def infection_population_density(unhealthy_prt, death_rate, mode="tanh"):

	def sigmoid(x):
		return 2 * (1 / (1 + np.exp(-x)) - 0.5)

	uhr = min(unhealthy_prt, PARA.unhealthy_prt_max) / PARA.unhealthy_prt_max
	dtr = min(death_rate, PARA.death_rate_max) / PARA.death_rate_max

	if mode == "tanh":
		pd = PARA.Pln_den_high - (PARA.Pln_den_high - PARA.Pln_den_low) * \
													np.tanh(2*(uhr + dtr))
	elif mode == "sigmoid":
		pd = PARA.Pln_den_high - (PARA.Pln_den_high - PARA.Pln_den_low) * \
													sigmoid(3*(uhr + dtr))
	else:
		raise ValueError("Unknown mode '%s'" % mode)

	return pd






class particle():

	def __init__(self, pid, status, time_stamp):

		self.id = pid
		self.stat_time = [time_stamp]   # time point of one particle to transfer to different status
		self.stat_list = [status]   	# status list, starts from "A0", the last one is current status
										# There is NO "E" in stat_list !
		#self.time_stamp = time_stamp
		self.future_status = []
		self.set_future_status()

	def one_step(self, time_stamp):
		'''
		if next_stat is not None and next_stat != self.stat_list[-1]:
			# status changed
			self.stat_list.append(next_stat)
			self.stat_time.append(self.time_stamp)
			self.set_future_status(next_stat)
		'''
		# time stamp reach future status
		if time_stamp == self.future_status[1]:
			# change to pre-calculated future
			f_next_stat = self.future_status[0]
			self.stat_list.append(f_next_stat)
			self.stat_time.append(time_stamp)
			self.set_future_status()

	def set_future_status(self):

		current_stat = self.stat_list[-1]
		current_time = self.stat_time[-1]

		if current_stat == STAT.A0:
			# A0 -> A1 or R or Q
			time_interval, next_status = A0_to_next(self.stat_list, self.stat_time)

		elif current_stat == STAT.A1:
			# A1 -> S or R or Q
			time_interval, next_status = A1_to_next(self.stat_list, self.stat_time)

		elif current_stat == STAT.S:
			# S -> R or Q
			time_interval, next_status = S_to_next(self.stat_list, self.stat_time)

		elif current_stat == STAT.Q:
			# Q -> R/D
			time_interval, next_status = Q_to_next(self.stat_list, self.stat_time)

		else:
			# others
			self.future_status = [-1, -1]
			return

		time_point = time_interval + current_time
		self.future_status = [next_status, time_point]





class population():

	pid = 0

	def __init__(self, population_init):

		self.N = 0
		self.evolution_time = 0
		self.population_init = population_init.copy()
		self.particles = []
		self.population_iter = {}
		self.population_density = {}

		self.population_iter[self.evolution_time] = self.population_init.copy()

		for stat, num in self.population_init.items():
			self.N += num
			if stat == STAT.E:
				continue
			for i in range(num):
				p = particle(population.pid, stat, self.evolution_time)
				self.particles.append(p)
				population.pid += 1


	def evolution(self, days):

		daily_S = {}

		for i in range(days):

			prev_plt = self.population_iter[self.evolution_time]
			new_plt = prev_plt.copy()
			if self.evolution_time not in daily_S.keys():
				daily_S[self.evolution_time] = 0
			
			# for groups except for E
			# if future time interval == 0
			for this_particle in self.particles:

				new_plt[this_particle.stat_list[-1]] -= 1

				this_particle.one_step(self.evolution_time)

				new_plt[this_particle.stat_list[-1]] += 1

				if this_particle.stat_list[-1] == STAT.Q \
					and this_particle.stat_time[-1] == self.evolution_time:
					daily_S[self.evolution_time] += 1

			# population density
			unhealthy_prt = (self.N - new_plt[STAT.E] - new_plt[STAT.R]) / self.N
			death_rate = new_plt[STAT.D] / (1e-5 + new_plt[STAT.Q] + new_plt[STAT.R] + new_plt[STAT.D])
			population_density = infection_population_density(unhealthy_prt, death_rate)

			self.evolution_time += 1
			if self.evolution_time not in daily_S.keys():
				daily_S[self.evolution_time] = 0

			# for groups except for E
			# if future time interval > 0
			for this_particle in self.particles:

				new_plt[this_particle.stat_list[-1]] -= 1

				this_particle.one_step(self.evolution_time)

				new_plt[this_particle.stat_list[-1]] += 1

				if this_particle.stat_list[-1] == STAT.Q \
					and this_particle.stat_time[-1] == self.evolution_time:
					daily_S[self.evolution_time] += 1

			# for E group
			E_to_A0 = (prev_plt[STAT.S] * PARA.beta_S \
				+ prev_plt[STAT.A1] * PARA.beta_A1) \
				* population_density
				#* prev_plt[STAT.E] / self.N
			E_to_A0 = int(round(E_to_A0))
			print("Day %d : E_to_A0 = %d" % (self.evolution_time, E_to_A0))
			# add new particles
			for i in range(E_to_A0):
				p = particle(population.pid, STAT.A0, self.evolution_time)
				self.particles.append(p)
				population.pid += 1
			new_plt[STAT.A0] += E_to_A0
			new_plt[STAT.E] -= E_to_A0

			self.population_iter[self.evolution_time] = new_plt.copy()
			self.population_density[self.evolution_time] = population_density

		return daily_S


if __name__ == '__main__':
	
	model = population(config.population_init)
	
	daily_S = model.evolution(150)

	AS = []
	E = []
	dE = []
	A0 = []
	QRD = []
	E_prev = model.population_init[STAT.E]
	for t, pp in model.population_iter.items():
		E.append(pp[STAT.E])
		AS.append(pp[STAT.A1] + pp[STAT.S])
		dE.append(E_prev - pp[STAT.E])
		A0.append(pp[STAT.A0])
		QRD.append(pp[STAT.Q] + pp[STAT.R] + pp[STAT.D])
		E_prev = pp[STAT.E]
		if AS[-1]+E[-1]+A0[-1]+QRD[-1] != model.N:
			raise RuntimeError("Population Error !")

	import matplotlib.pyplot as plt

	plt.subplot(1,3,1)
	plt.plot(list(daily_S.keys()), list(daily_S.values()), 'r.')
	plt.plot(list(daily_S.keys()), dE)
	plt.legend(["Daily_S", "diff(E)"])
	plt.ylabel("Counts")
	
	plt.subplot(1,3,2)
	plt.plot(list(daily_S.keys()), AS)
	plt.plot(list(daily_S.keys()), A0)
	plt.plot(list(daily_S.keys()), QRD)
	plt.legend(["A1+S","A0","Q+R+D"])
	plt.xlabel("Days")

	plt.subplot(1,3,3)
	plt.plot(list(model.population_density.keys()), \
				list(model.population_density.values()))
	plt.legend(["population density"])

	plt.show()








