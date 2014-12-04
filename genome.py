from thinkbayes2 import *
import random
import string

def Weighted_Choice(choices):
	'''
	Obtained from: http://stackoverflow.com/questions/3679694/a-weighted-version-of-random-choice
	Using a list of [choice, weight], pick a random value.
	'''
 	total = sum(w for c, w in choices)
 	r = random.uniform(0, total)
 	upto = 0
 	for c, w in choices:
 		if upto + w > r:
			return c
		upto += w
	assert False, "Shouldn't get here"

def Random_Genome(length_constraints=[100, 1000], weights=[1, 1, 1, 1]):
	'''
	Generate random genomic sequence for test purposes.
	Also allows weights of A, T, C, G in that order for further probability testing
	Based on the random word generator found: http://stackoverflow.com/questions/2030053/random-strings-in-python
	'''
	return ''.join(Weighted_Choice([['A', weights[0]], ['T', weights[1]], ['C', weights[2]], ['G', weights[3]]]) for i in xrange(random.randint(length_constraints[0], length_constraints[1])))

def Read_Genome(genomic_data, length_constraints=[-1, -1], overlap_constraints=[0.2, 0.6]):
	'''
	Returns the forward and backward read of a given section of genomic data
	'''
	# Handle no constraint input
	l = len(genomic_data)
	if(length_constraints[0] < 0):
		length_constraints[0] = l / 5
	if(length_constraints[1] < 0):
		length_constraints[1] = l / 2

	# Deal with poorly chosen constraints
	assert length_constraints[1] >= length_constraints[0], "Minimum length > Maximum length, please fix"
	assert overlap_constraints[1] >= overlap_constraints[0], "Minimum overlap > Maximum overlap, please fix"
	assert length_constraints[1] * (2 - overlap_constraints[0]) < l, "Potential to overread, please decrease length_constraints max or increase overlap_constraints min"

	# Determine, using the constraints, a random forward and backward read
	a = random.randint(length_constraints[0], length_constraints[1])
	b = random.randint(length_constraints[0], length_constraints[1])
	min_length = a if a < b else b
	overlap = int(random.uniform(overlap_constraints[0], overlap_constraints[1]) * min_length)

	# Build indices to use to define string and test later.
	idx = [-1, -1, -1, -1]
	idx[0] = random.randint(0, l - (a + b - overlap))	# Forward read begin
	idx[1] = idx[0] + a 								# Forward read end
	idx[2] = idx[1] - overlap 							# Reverse read begin
	idx[3] = idx[2] + b 								# Reverse read end

	forward_string = genomic_data[idx[0]:idx[1]]
	reverse_string = genomic_data[idx[2]:idx[3]]
	forward = Sequence(forward_string, True)
	reverse = Sequence(reverse_string, False)

	return [forward, reverse, idx]

class Sequence():
	'''
	Wrapper for the sequence list.  Includes original sequence and created sequence for some degree of
	'''
	stringdex = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
	dexstring = 'ATCG'
	def __init__(self, string, forward):
		self.true_sequence = string
		self.read_sequence = self.Read_Sequence(string, forward)

	def Read_Sequence(self, string, forward):
		'''
		string: The original genomic data
		forward: True if forward read, else False.  More accurate data on forward read
		Returns a 
		'''
		l = len(string)
		shift = 0 if forward else l
		read = [0] * l
		for i in range(l):
			read[i] = [0, 0, 0, 0]
			strength = 1
			sig = 0.4

			# Noise due to position in string
			for j in range(4):
				read[i][j] += abs(random.gauss(2 * strength * abs(shift - i) / l, sig))

			# Probability based on the true value
			read[i][self.stringdex[string[i]]] += abs(random.gauss(strength, sig))

			# Affected by nearby units
			if i > 0:
				read[i][self.stringdex[string[i-1]]] += abs(random.gauss(strength * 0.6, sig))
			if i < l - 1:
				read[i][self.stringdex[string[i+1]]] += abs(random.gauss(strength * 0.6, sig))
			
			# Normalize Probabilities
			total = sum(read[i])
			for j in range(4):
				read[i][j] = read[i][j] / total
		return read

	def Read(self):
		'''
		Reads the sequence of "Most Probable" value at each index
		'''
		return ''.join([self.dexstring[x.index(max(x))] for x in self.read_sequence])

	def Read_Seq(self):
		for row in self.read_sequence:
			print(', '.join([str(round(x, 4)) for x in row]))

	def Spoiler(self):
		return self.true_sequence


class Genome(Suite):
	def Likelihood(self, data, hypo):
		'''
		hypo: Starting index of the intended string
		data: [Index, A%, T%, C%, G%] Probabilities of a particular member of an index being A, T, C, or G
		'''
		return 0

if __name__ == '__main__':
	g = Random_Genome([50, 50])
	print(len(g))
	output = Read_Genome(g, [20, 20], [0.5, 0.6])
	print(output[0].Spoiler())
	print(output[0].Read())
	output[0].Read_Seq()


