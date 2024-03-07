from os import mkdir, listdir, rename
from os.path import exists, isdir
from math import log10, floor, ceil
import numpy as np
from itertools import combinations

class Sobol:
    
	def __init__(self, QoI, model, save, load, root):
		self.QoI = QoI
		self.model = list(model)
		self.root = root
		self.save = save
		self.load = load

	def format_path(self, path):
		# aggiungo zeri alle cartelle che sono in path
		# hanno tutte la stessa lunghezza, mia nevrosi
		if not exists(path): mkdir(path)
		format = ceil(log10(len(listdir(path))))
		for file in listdir(path):
			if isdir(path + '/' + file):
				while len(file) < format:
					rename(path + '/' + file, path + '/0' + file)
					file = '0' + file
		return format

	def simulations(self, path, shape, args):
		# intanto formatto path
		step = len(args)
		assert step + len(shape) == len(self.model)
		format = self.format_path(path)
		for k in range(shape[0]):
			# ora devo fare shape[0] simulazioni
			str_fold = str(k)
			if format >= len(str_fold): str_fold = (format - len(str_fold)) * '0' + str_fold
			try: arg = self.load[step](path + '/' + str_fold)
			except: arg = self.moldel[step](args)
			# voglio per ogni simulazione creare
			# un tensore di dimensione shape[1:]
			if len(shape > 0): self.simulations(self, path + '/' + str_fold, shape[1:], args + [arg])
			self.save[step](arg, path + '/' + str_fold)

	def sensitivity(self, path, shape, args):
		step = len(args)
		assert step + len(shape) == len(self.model)
		if len(shape) == 0: x = self.QoI(args)
		else:
			try:
				# se tensore disponibile
				# e di dimensione giusta
				x = np.load(path + '/tensor.npy')
				assert (x.shape >= shape).all()
			except AssertionError:
				# o concatena tensori di
				# dimensione = shape[1:]
				x = []
				format = self.format_path(path)
				for samp in range(shape[0]):
					file = path + '/' + (format - len(str(samp))) * '0' + str(samp)
					x.append(self.sensitivity(file, shape[1:], args + [self.load[step](file)]))
				x = np.array(x)
		np.save(path + '/tensor.npy', x)
		return x

	def __call_(self, shape):
		# prima fai le simulazioni
		self.simulations(self.root, shape, [])
		# poi ti calcoli i tensori
		tensor = self.sensitivity(self.root, shape, [])
		# poi fai order statistics
		steps = len(shape)
		for step in range(steps):
			tensor = tensor.sort(axis=(steps - step))
		# ti fai le medie quadrate
		pick = []
		for r in range(steps + 1): pick += list(combinations(range(steps), r))
		freeze = np.array([(tensor.mean(axis=x)**2).mean() for x in pick])
		freeze.reverse()
		# ora sottrai per varianze
		X = {}
		Y = np.zeros(steps)
		for p, v in zip(pick, freeze):
			temp = v
			for k in range(len(p)):
				for q in combinations(p, k): temp -= X[q]
			X[p] = temp
			for a in p: Y[a] += temp
		# ora se vuoi normalizzare
		Y -= Y[0]
		Y /= np.sqrt(Y[-1])
		return {i: s for i, s in (zip(pick, Y))[1:-1]}