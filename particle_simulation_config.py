import sys
import os
import numpy as np


class STAT():

	E = 0    # 易感
	A0 = 1   # 不传染潜伏
	A1 = 2   # 传染潜伏
	S = 3    # 传染发病
	Q = 4    # 确诊隔离/医疗介入
	R = 5    # 恢复不易感
	D = 6    # 死亡


population_init = {
	STAT.E : 999900,
	STAT.A0 : 100,
	STAT.A1 : 0,
	STAT.S : 0,
	STAT.Q : 0,
	STAT.R : 0,
	STAT.D : 0
}


class PARA():

	# A1 infection rate
	beta_A1 = 0.3
	# S infection rate
	beta_S = 0.1

	# A0 to A1, R, Q
	A0_to_A1 = 0.9
	A0_to_R = 0.1
	A0_to_Q = 0.0
	# A0 to A1 time distribution
	mu_A0_to_A1 = 2
	sg_A0_to_A1 = 0.5
	# A0 to R time distribution
	mu_A0_to_R = 14
	sg_A0_to_R = 3
	# A0 to Q time distribution
	mu_A0_to_Q = 6
	sg_A0_to_Q = 2
	# A0 to D time distribution
	mu_A0_to_D = 25
	sg_A0_to_D = 10

	# A1 to S, Q, R
	A1_to_S = 0.8
	A1_to_Q = 0.00
	A1_to_R = 0.2
	# A1 to S time distribution
	mu_A1_to_S = 5
	sg_A1_to_S = 3
	# A1 to Q time distribution
	mu_A1_to_Q = 6
	sg_A1_to_Q = 2
	# There is no A1 to R time distribution,
	# because the time are accumulated from A0 status

	# S to Q, R, D
	S_to_Q = 0.8
	S_to_R = (1 - S_to_Q) * 0.96
	S_to_D = (1 - S_to_Q) * 0.04
	# S to Q time distribution
	mu_S_to_Q = 1
	sg_S_to_Q = 0.5
	# There is no S to R time distribution,
	# because the time are accumulated from A0 status.
	
	# Q to R, D
	Q_to_R = 0.99
	Q_to_D = 1 - Q_to_R
	# There is no Q to R/D time distribution,
	# because the time are accumulated from A0 status
	# But, medical treatment will make life longer, and recovery time shorter
	Q_longer_D = 1.5
	Q_shorter_R = 0.9

	# population density (0~1)
	# It has a positive correlation with unhealthy population percent
	# and death rate, in "tanh" formula. tanh(4)=0.999, tanh(2)=0.96
	# Pln_den_high - (Pln_den_high-Pln_den_low) * tanh(2*unhealthy_percent/unhealthy_max + 2*death_rate/death_rate_max)
	Pln_den_high = 1
	Pln_den_low = 0.1
	unhealthy_prt_max = 0.01
	death_rate_max = 0.02











