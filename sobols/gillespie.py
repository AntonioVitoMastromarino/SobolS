class Gillespie:

	def __init__(self, states, reacts, init_step, next_step, params=None, seed=None):
		#react['befor'] -> react['after']
		self.states = states
		self.reacts = reacts
		self.params = params
		self.acting = {state: {react['befor'] for react in reacts if state in react['befor']} for state in states}
		self.next_step = next_step
		self.init = init_step
		self.elem = {react['befor']: {} for react in self.reacts}
		

	def update(self, k, state, elem, x, contact, num):
		for h in contact[k]:
			if (x[h] in self.back_react[x[k]]):
				num[(x[h], x[k])] -= 1
				elem[(x[h], x[k])].remove((h, k))
			if (x[h] in self.forw_react[x[k]]):
				num[(x[k], x[h])] -= 1
				elem[(x[k], x[h])].remove((k, h))
			if (x[h] in self.back_react[state]):
				num[(x[h], state)] += 1
				elem[(x[h], state)] += [(h, k)]
			if (x[h] in self.forw_react[state]):
				num[(state, x[h])] += 1
				elem[(state, x[h])] += [(k, h)]
		num[x[k]] -= 1
		num[state] += 1
		elem[x[k]].remove(k)
		elem[state] += [k]
		x[k] = state

	def stochastic(self, days, callback=None):
		x = self.init_step
		elem = {elem[react['catal'], react['befor']]: {} for react in self.reacts}
		while (t < days):
			catal, befor, react = self.next_step(x, elem)
			for k, state in zip(befor, react['after']):
				self.update(k, state, elem, x)
			t = callback(catal, befor, react, t)