from thinkbayes2 import *
import random
import string
import sys
import os

class Settings():
	d = dict()

	def __init__(self):
		self.Set("sigma_correct", 1.0)
		self.Set("sigma_wrong", 1.0)
		self.Set("curve_strength", 1.0)
		self.Set("read_strength", 4.0)
		self.Set("error_max", 4.0)
		self.Set("A_weight", 1)
		self.Set("T_weight", 1)
		self.Set("C_weight", 1)
		self.Set("G_weight", 1)

		self.Set("min-len", 100)
		self.Set("max-len", 150)
		self.Set("min-overlap", 0.2)
		self.Set("max-overlap", 0.8)
		self.Set("max-runs", -1)
		self.ParseArgs()

	def Get(self, setting):
		if setting not in self.d.keys():
			print("Setting not found")
			return
		return self.d[setting]

	def Set(self, setting, value):
		self.d[setting] = value

	def ParseArgs(self):
		args = sys.argv
		for i in range(len(args)):
			if(args[i][-1] is "=" and i + 1 < len(args)):
				print(args[i] + args[i+1])
				self.Set(args[i][:-1], float(args[i+1]))

print(sys.argv)

s = Settings()

sigma_correct = s.Get('sigma_correct')
sigma_wrong = s.Get('sigma_wrong')
curve_strength = s.Get('curve_strength')
read_strength = s.Get('read_strength')
max_val = s.Get('error_max')

def Prob_True(portion):
	try:
		loss = s.Get('error_max') * portion**s.Get('curve_strength')
	except ValueError:
		loss = 0
	pmf_true = MakeNormalPmf(s.Get('read_strength'), s.Get('sigma_correct'), 5, n=100)
	pmf_false = MakeNormalPmf(loss, s.Get('sigma_wrong'), 5, n=100)
	return pmf_true > pmf_false

def drange(start, stop, steps):
	return [1.0 * x * (stop-start) / steps + start for x in range(steps+1)]

def Prob_True_Graph(err_max, sig_w, read_str, sig_c, cur_str):
	wise = [0] * 101
	count = 0
	for i in drange(0, 1, 100):
		try:
			loss = err_max * i**cur_str
		except ValueError:
			loss = 0
		pmf_true = MakeNormalPmf(read_str,sig_c, 5, n=100)
		pmf_false = MakeNormalPmf(loss, sig_w, 5, n=100)
		wise[count] = [i, pmf_true > pmf_false]
		count = count + 1
	fname = str(err_max) + " " + str(cur_str) + " " + str(read_str) + " " + str(sig_c) +  " " + str(sig_w) + ".txt"
	f = open(fname, 'w')
	for y in wise:
		f.write(str(y[0]) + " " + str(y[1]) + "\n")
	f.close()
	return wise

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
	Takes a string and based on the constraints, returns two reads, a forward and reverse.
	The Sequence object returned has the true data along with the "read" data.
	The last thing returned is a list of the relevant indices [f_beg, f_end, r_beg, r_end]
	'''
	# Handle no constraint input
	l = len(genomic_data)
	if(length_constraints[0] < 0):
		length_constraints[0] = l / 5
	if(length_constraints[1] < 0):
		length_constraints[1] = l / 2

	# Deal with poorly chosen constraints
	assert length_constraints[1] >= length_constraints[0], "Minimum length (" + str(length_constraints[0])+") > Maximum length(" + str(length_constraints[1])+"), please fix"
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

			# Uses Gaussian Random to check if the value is greater, if fails, then returns a random
			# read_truthfully = random.gauss(read_strength, sigma_correct) > random.gauss(loss, sigma_wrong)
			# Might consider making Weighted_Choice always
			portion = 1.0 * i / l if forward else 1 - 1.0 * i/l
			read_true = Prob_True(portion) > random.uniform(0, 1)
			read[i] = string[i] if read_true else Weighted_Choice([['A', 1], ['T', 1], ['C', 1], ['G', 1]])
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
		if(index_match == self.forward_length):
			# Most probably true at later values
			return(1.0 * index / self.min_length)
		if(index_match > self.forward_length):
			return 1

		forward_base = self.Get_At(index_match)
		reverse_base = data[1]

		if(forward_base == reverse_base):
			# Match Found, determine chance of match accounting for False Positive

			# Accounts for random chance of reading true
			# Forward Probability of reading correctly

			f_portion = 1.0 * index_match / self.forward_length
			p = (1 + 3 * Prob_True(f_portion))/4
			# Reverse Probability of reading correctly
			r_portion = 1 - 1.0 * index / self.reverse_length
			k = (1 + 3 * Prob_True(r_portion))/4

			positive = p * k + (1-p)*(1-k)/3
			false_positive = (p + k - 2*p*k)/3 + (1-p)*(1-k) * 2/9
			return positive / (positive + false_positive)
		else:
			# Forward Probability of reading correctly
			f_portion = 1.0 * index_match / self.forward_length
			p = (1 + 3 * Prob_True(f_portion))/4
			# Reverse Probability of reading correctly
			r_portion = 1 - 1.0 * index / self.reverse_length
			k = (1 + 3 * Prob_True(r_portion))/4

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

def ProbGraphs(fname):


def main():
	# Define a Random Genomic Sequence owith equally weighted A, T, C, G
	genomic_sequence = Random_Genome([10000, 12000], [1, 1, 1, 1])

	# Get the probabilistic data for each member of the sequence based on reading errors
	forward_data, reverse_data, indices = Read_Genome(genomic_sequence, 
		[s.Get('min-len'), s.Get('max-len')], 
		[s.Get('min-overlap'), s.Get('max-overlap')])

	val = 0
	while(os.path.isfile("log" + str(val) + ".txt")):
		val = val + 1

	f = open("log" + str(val) + "values.txt", 'w')

	# Use suite to find the most probable overlap between the two pieces of data
	hypos = range(int(0.5 * s.Get('min-overlap') * forward_data.Length()), int((1 - s.Get("min-overlap")) * forward_data.Length())) # Theoretically, the overlap can occur at any point in the sequence
	suite = Genome_Matcher(hypos)
	suite.Set_Forward(forward_data.Read())
	reverse_read = reverse_data.Read()
	suite.Set_Lengths(len(forward_data.Read()), len(reverse_read))

	# Deconstructs the reverse read string into a dataset and updates
	dataset = [[idx, reverse_read[idx]] for idx in range(len(reverse_read))]
	random.shuffle(dataset)

	f.write("Command: " + ' '.join(sys.argv[2:]))
	f.write("History:\n")
	if(s.Get("max-runs") > 1):
		print("Only performing a limited number of tests")
		dataset = dataset[0:int(s.Get("max-runs"))]
	suite.UpdateSetAndLog(dataset, f)
	f.close()

	probs = suite.MostProbable()[0:5]
	for prob in probs:
		print(str(prob[0]) + ": " + str(100 *prob[1])[0:5] + "%")
	shift = probs[0][0]
	print("True Shift: " + str(indices[2] - indices[0]) + " Most Likely: " + str(shift))
	print("Forward Read: " + forward_data.Read())
	print("Forward True: " + forward_data.Read_True())
	print("Reverse Read: " + reverse_data.Read())
	print("Reverse True: " + reverse_data.Read_True())

	reconstructed = Reconstruct(forward_data, reverse_data, shift)

	print("True Overlap")
	print(forward_data.Read())
	print(int(indices[2]-indices[0])* "." + reverse_data.Read())

	print("Probable Overlap")
	print(forward_data.Read())
	print(shift * "." + reverse_data.Read())

	print(reconstructed)

	f = open("log" + str(val) + ".txt", 'w')
	f.write("Command: " + ' '.join(sys.argv))
	f.write("\nForward Read: " + forward_data.Read())
	f.write("\nForward True: " + forward_data.Read_True())
	f.write("\nReverse Read: " + reverse_data.Read())
	f.write("\nReverse True: " + reverse_data.Read_True())
	f.write("\n")

	for i in range(len(probs)):
		f.write("\n" + str(prob[0]) + ": " + str(100 *prob[1])[0:5] + "%")
	f.write("\nTrue Shift: " + str(indices[2] - indices[0]) + " Most Likely: " + str(shift))
	f.write("\nTrue Overlap\n")
	f.write(forward_data.Read() + "\n")
	f.write(int(indices[2]-indices[0]) * " " + reverse_data.Read())
	f.write("\nMost Probable Overlap\n")
	f.write(reverse_data.Read() + "\n")
	f.write(shift * " " + reverse_data.Read())
	f.write("\nReconstructed Sequence: " + reconstructed)
	f.write("\nTrue Sequence:          " + genomic_sequence[indices[0]:indices[3]])
	f.write("\nErrors in the sequence: ")
	f.write(''.join([' ' if reconstructed[x] == genomic_sequence[indices[0] + x] else "x" for x in xrange(min([len(genomic_sequence), len(reconstructed)]))]))
	# print(forward_data.Read_True())
	# print(genomic_sequence[indices[0]:indices[3]])
	f.close()
	return 0

if __name__ == '__main__':
	main()
