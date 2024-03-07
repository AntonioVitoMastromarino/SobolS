from os import mkdir, listdir, rename
from os.path import exists, isdir
from math import log10, floor, ceil
import numpy as np

class Sobol:
    
	def __init__(self, QoI, model, save, load, root):
		self.QoI = QoI
		self.model = list(model)
		self.root = root
		self.save = save
		self.load = load

	def format_path(self, path):
		if not exists(path): mkdir(path)
		format = ceil(log10(len(listdir(path))))
		for file in listdir(path):
			if isdir(path + '/' + file):
				while len(file) < format:
					rename(path + '/' + file, path + '/0' + file)
					file = '0' + file
		return format

	def simulations(self, path, shape, args):
		step = len(args)
		assert step + len(shape) == len(self.model)
		format = self.format_path(path)
		for k in range(shape[0]):
			str_fold = str(k)
			if format >= len(str_fold): str_fold = (format - len(str_fold)) * '0' + str_fold
			try: arg = self.load[step](path + '/' + str_fold)
			except: arg = self.moldel[step](args)
			if len(shape > 0): self.simulations(self, path + '/' + str_fold, shape[1:], args + [arg])
			self.save[step](arg, path + '/' + str_fold)

	def sensitivity(self, path, shape, args):
		step = len(args)
		assert step + len(shape) == len(self.model)
		if len(shape) == 0: x = self.QoI(args)
		else:
			try:
				x = np.load(path + '/tensor.npy')
				assert (x.shape >= shape).all()
			except AssertionError:
				x = []
				format = self.format_path(path)
				for samp in range(shape[0]):
					file = path + '/' + (format - len(str(samp))) * '0' + str(samp)
					x.append(self.sensitivity(file, shape[1:], args + [self.load[step](file)]))
				x = np.array(x)
		np.save(path + '/tensor.npy', x)
		return x

	def __call_(self, shape):
		self.simulations(self.root, shape, [])
		return self.sensitivity(self.root, shape, [])