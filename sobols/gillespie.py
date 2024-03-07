from networkx import connected_watts_strogatz_graph as nxs
from networkx import powerlaw_cluster_graph as nxl
from itertools import starmap, product, combinations
from time import time
from matplotlib import pyplot
from multiprocessing import Pool
import os
from numpy import sqrt, mean, array, matrix, ndarray, argsort, tensordot, sort, ones, zeros
from numpy.linalg import norm
from numpy.random import choice, rand
from sobol import Sobol

states = ['S', 'A', 'I', 'R']
colors = {'S': 'tab:blue', 'A': 'tab:orange', 'I': 'tab:green', 'R': 'tab:red'}
back_reactions = {'S': ['A', 'I'], 'A': [], 'I': [], 'R': []}
forw_reactions = {'S': [],'A': ['S'], 'I': ['S'], 'R': []}
reaction = ['alpha', 'betaA', 'betaI', 'gamma', 'detaA', 'detaI', 'varnu']
interval = {'alpha': (0.001, 0.9), 'betaA': (0.5, 0.9), 'betaI': (0.5, 0.9), 'gamma': (0.001, 0.05), 'detaA': (0.1, 0.51), 'detaI': (0.1, 0.51), 'varnu': (0, 0.02)}
befor = {'alpha': 'A', 'varnu': 'S', 'gamma': 'R', 'detaA': 'A', 'detaI': 'I'}
befor |= {'betaA': ('A', 'S'), 'betaI': ('I', 'S')}
after = {'alpha': 'I', 'varnu': 'R', 'gamma': 'S', 'detaA': 'R', 'detaI': 'R'}
after |= {'betaA': 'A', 'betaI': 'A'}

def update(k, state, elem, x, contact, num):
	for h in contact[k]:
		if (x[h] in back_reactions[x[k]]):
			num[(x[h], x[k])] -= 1
			elem[(x[h], x[k])].remove((h, k))
		if (x[h] in forw_reactions[x[k]]):
			num[(x[k], x[h])] -= 1
			elem[(x[k], x[h])].remove((k, h))
		if (x[h] in back_reactions[state]):
			num[(x[h], state)] += 1
			elem[(x[h], state)] += [(h, k)]
		if (x[h] in forw_reactions[state]):
			num[(state, x[h])] += 1
			elem[(state, x[h])] += [(k, h)]
	num[x[k]] -= 1
	num[state] += 1
	elem[x[k]].remove(k)
	elem[state] += [k]
	x[k] = state

def init(citizens, net):
	dgs = zeros(citizens)
	max = 1
	arg = []
	for (k, h) in net:
		dgs[k] += 1
		dgs[h] += 1
		if (dgs[k] == max):
			arg += [k]
			max = dgs[k]
		if (dgs[k] > max):
			arg = [k]
			max = dgs[k]
		if (dgs[h] == max):
			arg += [h]
			max = dgs[h]
		if (dgs[h] > max):
			arg = [h]
			max = dgs[h]
	temp = array(citizens * ['S'])
	for k in arg: temp[k] = 'A'
	return temp

def stochastic(net, Xparam, days, citizens):
	x = init(citizens, net)
	num = {}
	elem = {}
	for rea in reaction:
		num[befor[rea]] = 0
		elem[befor[rea]] = []
	contact = []
	for k in range(citizens):
		contact += [[]]
		num[x[k]] += 1
		elem[x[k]] += [k]
	for (k, h) in net:
		contact[k] += [h]
		contact[h] += [k]
		if x[h] in back_reactions[x[k]]:
			elem[(x[h], x[k])] += [(h, k)]
			num[(x[h], x[k])] += 1
		if x[k] in back_reactions[x[h]]:
			elem[(x[k], x[h])] += [(k, h)]
			num[(x[k], x[h])] += 1
	y = {'time': [0]}
	for state in states: y[state] = [num[state]]
	t = 0
	while (t < days):
		randomic = rand()
		prob = {}
		rate = 0
		for rea in reaction:
			prob[rea] = num[befor[rea]] * Xparam[rea]
			rate += prob[rea]
		dt = 1 / rate
		t += dt
		for rea in reaction:
			prob[rea] *= dt
		pro = 0
		count = 0
		while (randomic > pro + prob[reaction[count]]):
			pro += prob[reaction[count]]
			count += 1
		rea = reaction[count]
		k = int(float(num[befor[rea]]) * (randomic - pro) / prob[rea])
		if (befor[rea] in states):
			k = elem[befor[rea]][k]
		else:
			h, k = elem[befor[rea]][k]
		update(k, after[rea], elem, x, contact, num)
		for state in states:
			y[state] += [num[state]]
		y['time'] += [t]
		if (num['A'] == num['I'] == 0):
			t = days
			y['time'][- 1] = days
	return y

def gen_net(args):
	assert args == []
	topology = choice(zip([nxs, nxl], ['nxs', 'nxl']))
	citizens = choice([2**k for k in range(8,16)])
	d = 3
	p = rand()
	net = list(topology[0](citizens, d, p).edges())
	return {'topology': topology[1], 'citizens': citizens, 'd': d, 'p': p, 'net': net, 'init': init(citizens, net)}

def save_net(arg, path): pass
	
def load_net(path): pass

def gen_par(args):
	pars = {par: interval[par][0] + rand() * (interval[par][1] - interval[par][0]) for par in reaction}
	pars['betaA'] /= args[0]['d']
	pars['betaI'] /= args[0]['d']
	return pars

def save_par(arg, path): pass

def load_par(path): pass

def gen_sim(args):
	net = args[0]['net']
	citizens = args[0]['citizens']
	days = 16
	param = args[1]
	return stochastic(net, param, days, citizens)

def save_sim(arg, path): pass

def load_sim(path): pass

def QoI(args):
	traj = args[- 1]
	return (array(traj['A']) + array(traj['I'])).max()

if __name__ == '__main__':
	model = [gen_net, gen_par, gen_sim]
	save = [save_net, save_par, save_sim]
	load = [load_net, load_par, load_sim]
	Sobol(QoI, model, save, load, '')([3,3,3])