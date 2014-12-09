from thinkbayes2 import *
import random
import string
import sys

class Settings():
	d = dict()

	def __init__(self):
		self.Set("sigma_correct", 1)
		self.Set("sigma_wrong", 1)
		self.Set("curve_strength", 1)
		self.Set("read_strength", 1)
		self.Set("error_max", 1)
		self.Set("A_weight", 1)
		self.Set("T_weight", 1)
		self.Set("C_weight", 1)
		self.Set("G_weight", 1)

		self.Set("seq_min_len", 100)
		self.Set("seq_max_len", 150)
		self.Set("seq_min_overlap", 0.2)
		self.Set("seq_max_overlap", 0.8)
		self.ParseArgs()

	def Get(self, setting):
		return d[setting]

	def Set(self, setting, value):
		d[setting] = value

	def ParseArgs(self):
		args = sys.argv
		for i in range(len(args)):
			if(args[i][0] == "-" and i + 1 < len(args)):
				self.Set(args[i][1:], float(args[i+1]))

s = Settings()

sigma_correct = 10
sigma_wrong = 10
curve_strength = 1
read_strength = 50
max_val = 50

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
	Wrapper for the sequence list.  Includes original sequence and created sequence with degree of randomness
	'''
	stringdex = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
	dexstring = 'ATCG'
	def __init__(self, string, forward):
		self.true_sequence = string
		self.read_sequence = self.Read_Sequence(string, forward)

	def Read_Sequence(self, string, forward):
		'''
		string: The original genomic data
		forward: True if forward read, else False.  More accurate data on forward read earlier
		Returns a 
		'''
		l = len(string)
		read = ['0'] * l

		shift = 0.0 if forward else 1.0 * l
		for i in range(l):

			# More loss later for forward, less for reversed
			# increasing curve strength decreases the effects earlier on
			# Maxes out at 1 for an "equally likely" scenario of being wrong
			loss = max_val * abs((shift-i)/l)**curve_strength

			# Uses Gaussian Random to check if the value is greater, if fails, then returns a random
			read_truthfully = random.gauss(read_strength, sigma_correct) > random.gauss(loss, sigma_wrong)
			# Might consider making Weighted_Choice always
			read[i] = string[i] if read_truthfully else Weighted_Choice([['A', 1], ['T', 1], ['C', 1], ['G', 1]])
		return ''.join(read)

	def Read(self):
		'''
		Reads the sequence of "Most Probable" value at each index
		'''
		# return ''.join([self.dexstring[x.index(max(x))] for x in self.read_sequence])
		return self.read_sequence

	def Read_Prob(self):
		for row in self.read_sequence:
			print(', '.join([str(round(x, 4)) for x in row]))

	def Read_True(self):
		return self.true_sequence

	def Length(self):
		return len(self.true_sequence)

class Genome_Matcher(Suite):
	def Likelihood(self, data, hypo):
		'''
		hypo: Starting index of the intended string
		data: [Index, reverse_base] The read Base at Index

		Returned probability is the probability that they were matched

		Needs to account for:
		Increasing error in forward string as well as decreasing in reverse
		'''
		index = data[0]
		index_match = hypo + index

		# Deals with the no-match case
		if(index_match >= self.forward_length):
			# Most probably true at later values
			return (1.0 * index / self.min_length)

		forward_base = self.Get_At(index_match)
		reverse_base = data[1]

		if(forward_base == reverse_base):
			# Match Found, determine chance of match accounting for False Positive

			# Accounts for random chance of reading true
			# Forward Probability of reading correctly
			p = self.Prob_True(index_match, True)
			# Reverse Probability of reading correctly
			k = self.Prob_True(index, False)

			positive = p * k + (1-p)*(1-k)/3
			false_positive = (p + k - 2*p*k)/3 + (1-p)*(1-k) * 2/9
			return positive / (positive + false_positive)
		else:
			# Forward Probability of reading correctly
			p = self.Prob_True(index_match, True)
			# Reverse Probability of reading correctly
			k = self.Prob_True(index, False)

			false_negative = 1 - (p * k + (1-p)*(1-k)/3)
			negative = 	1 - ((p + k - 2*p*k)/3 + (1-p)*(1-k) * 2/9)
			return false_negative / (negative + false_negative)

	def Get_At(self, index):
		return self.forward_string[index]

	def Set_Forward(self, forward_string):
		self.forward_string = forward_string

	def Set_Lengths(self, for_l, rev_l):
		self.forward_length = for_l
		self.reverse_length = rev_l
		self.min_length = for_l if for_l < rev_l else rev_l

	def Prob_True(self, index, is_forward):
		length = self.forward_length if is_forward else self.reverse_length
		loss = max_val * (1.0 * index / length)**curve_strength if is_forward else max_val*((length - 1.0 * index)/length)**curve_strength
		pmf_true = MakeNormalPmf(read_strength, sigma_correct, 5, n=101)
		pmf_false = MakeNormalPmf(loss, sigma_wrong, 5, n=101)
		return (1 + 3 * (pmf_true > pmf_false)) * 0.25


'''
Decision Tree:


MATCH:

Assume both are truly (A) [That is, a match exists, happens 0.25 of the time]

	Chances of being each base for Forward
	|    A    |    T    |    C    |    G    |
	|    p    | (1-p)/3 | (1-p)/3 | (1-p)/3 |
	For Reverse
	|    A    |    T    |    C    |    G    |
	|    k    | (1-k)/3 | (1-k)/3 | (1-k)/3 |

Positive Match (P(M^R)):
p * k + (1-p)(1-k)/3

False No Match (P(M^~R)):
1-(p * k + (1-p)(1-k)/3)

NO MATCH:

Assume one is truly (A) and the other is truly (T), no match exists, happens 0.75 of the time

	Chances of being each base for Forward
	|    A    |    T    |    C    |    G    |
	|    p    | (1-p)/3 | (1-p)/3 | (1-p)/3 |
	For Reverse
	|    A    |    T    |    C    |    G    |
	| (1-k)/3 |    k    | (1-k)/3 | (1-k)/3 |

No Match (P(~M^~R)):
1 - (p + k - 2pk)/3 + (1-p)(1-k) * 2/9

False Match (P(~M^R)):
(p + k - 2pk)/3 + (1-p)(1-k) * 2/9

k(1-p)/3 + 2(1-p)(1-k)/9

Probability that Match exists given that match is read:
P(M^R) / [P(M^R) + P(~M^R)]


Probability that a match exists given that data doesn't match:
P(M^~R) / [P(M^~R) + P(~M^~R)]
'''


def Test_Genome():
	g = Random_Genome([50, 50])
	output = Read_Genome(g, [20, 20], [0.5, 0.6])
	print("Overlap: " + str(output[2][1] - output[2][2]))
	print(output[0].Read_True())
	print(output[0].Read())

	print(output[1].Read_True())
	print(output[1].Read())
	#output[0].Read_Seq()

def Reconstruct(forward, reverse, offset):
	'''
	Takes two sequences and attempts to reconstruct them.
	For the overlapping portion: First half is the forward, second is reverse.
	'''
	front = (len(forward.Read()) + offset)/2
	back = front - offset
	first = forward.Read()[0:front]
	last = reverse.Read()[back:]
	return first + last

def main():
	# Defaults
	seq_min_len = 10
	seq_max_len = 15
	seq_min_overlap = 0.3
	seq_max_overlap = 0.7

	args = sys.argv
	if(len(args) > 1):
		seq_min_len = int(args[1])
		seq_max_len = int(args[2])
		seq_min_overlap = float(args[3])
		seq_max_overlap = float(args[4])

	# Define a Random Genomic Sequence owith equally weighted A, T, C, G
	genomic_sequence = Random_Genome([10000, 12000], [1, 1, 1, 1])

	# Get the probabilistic data for each member of the sequence based on reading errors
	forward_data, reverse_data, indices = Read_Genome(genomic_sequence, [seq_min_len, seq_max_len], [seq_min_overlap, seq_max_overlap])

	# Use suite to find the most probable overlap between the two pieces of data
	hypos = range(int(0.5 * seq_min_overlap * forward_data.Length()), forward_data.Length()) # Theoretically, the overlap can occur at any point in the sequence
	suite = Genome_Matcher(hypos)
	suite.Set_Forward(forward_data.Read())
	reverse_read = reverse_data.Read()
	suite.Set_Lengths(len(forward_data.Read()), len(reverse_read))

	# Deconstructs the reverse read string into a dataset and updates
	dataset = [[idx, reverse_read[idx]] for idx in range(len(reverse_read))]
	suite.UpdateSet(dataset)

	print(indices[2] - indices[0])
	suite.Print()
	shift = suite.MostProbable()[0][0]

	print(forward_data.Read())
	print(forward_data.Read_True())
	print(reverse_data.Read())
	print(reverse_data.Read_True())

	print(shift)
	reconstructed = Reconstruct(forward_data, reverse_data, shift)
	print(reconstructed)
	print(genomic_sequence[indices[0]:indices[3]])
	return 0

if __name__ == '__main__':
	main()
