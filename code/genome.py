from thinkbayes2 import *
from scipy import stats
import random
import string
import sys
import os
import matplotlib.pyplot as plt

class Argument_Parser():
	'''
	Simple object used to parse and set different command line arguments.
	This is not the primary part of the code, but its stored values are used
	  throughout in order to change the outcome.
	TODO: Make a HELP option that describes what everything does.
	'''
	d = dict()

	def __init__(self):
		'''
		Default Values for the program
		'''

		# Values for the Error function
		self.Set("curve", 1.0)
		self.Set("prob-min", 0.0)
		self.Set("prob-max", 1.0)

		# Unused values for a previous iteration of the code.
		self.Set("sigma-correct", 1.0)
		self.Set("sigma-wrong", 1.0)
		self.Set("read-strength", 4.0)
		self.Set("error_max", 4.0)

		# Weights of each base pair. Currently unused
		self.Set("A-weight", 1)
		self.Set("T-weight", 1)
		self.Set("C-weight", 1)
		self.Set("G-weight", 1)

		# Limits of the matched genome length
		self.Set("length-min", 100)
		self.Set("length-max", 150)
		self.Set("overlap-min", 0.2)
		self.Set("overlap-max", 0.8)
		self.Set("runs-max", -1) # Indicates the max number of trials to run
		self.ParseArgs()

	def Get(self, setting):
		if setting not in self.d.keys():
			raise ValueError('Setting not found: "' + setting + '"')
		return self.d[setting]

	def Set(self, setting, value):
		self.d[setting] = value

	def ParseArgs(self):
		'''
		Iterates through the command line arguments, looks for arguments in the form:
		(setting)= value
		And adds a dictionary value that can be called later.

		At the moment, this only takes in float arguments, so it's just enough to actually run.
		'''
		args = sys.argv
		for i in range(len(args)):
			'''
			Iterates through and finds if the last character is a "="
			Then assigns each of them accordingly.
			'''
			if(args[i][-1] is "=" and i + 1 < len(args)):
				print(args[i] + args[i+1])
				self.Set(args[i][:-1], float(args[i+1]))

print(sys.argv)

s = Argument_Parser()

def Prob_True(portion):
	'''
	Indicates the probability that the returned value correct.
	portion is the fraction through the string ()
	
	Requirements:
	Perfect (1) reading at portion=0
	Lowest reading (minimum 0) at portion = 1

	At all points in between: dProb_True/dPortion < 0 
		(Meaning probability is always decreasing)
	'''

	# 1-(p)^c
	scale = s.Get("prob-max") - s.Get("prob-min")
	y_intercept = s.Get("prob-max")
	if(scale > 1 or scale < 0):
		scale = 1
		y_intercept = 1.0

	return y_intercept - portion ** s.Get('curve') * scale
	# try:
	# 	loss = s.Get('error_max') * portion**s.Get('curve_strength')
	# except ValueError:
	# 	loss = 0
	# pmf_true = MakeNormalPmf(s.Get('read_strength'), s.Get('sigma_correct'), 5, n=100)
	# pmf_false = MakeNormalPmf(loss, s.Get('sigma_wrong'), 5, n=100)
	# return pmf_true > pmf_false

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
	Later, I might change this to extend strings
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
		Returns the read of the sequence without true values
		'''
		# return ''.join([self.dexstring[x.index(max(x))] for x in self.read_sequence])
		return self.read_sequence

	def Read_True(self):
		'''
		Reads
		'''
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
		index        = data[0]
		reverse_base = data[1]
		index_match  = hypo + index

		# Deals with the no-match case because index is out of range
		if(index_match == self.forward_length):
			# Most probably true at later values
			return(1.0 * index / self.min_length)
		if(index_match > self.forward_length):
			return 1

		forward_base = self.forward_string[index_match]

		if(forward_base == reverse_base):
			# Match Found, determine chance of match accounting for False Positive

			# Accounts for random chance of reading correctly because of weighted random
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
			# No match found, what is the possibility that it's a false negative?
			# Forward Probability of reading correctly
			f_portion = 1.0 * index_match / self.forward_length
			p = (1 + 3 * Prob_True(f_portion))/4
			# Reverse Probability of reading correctly
			r_portion = 1 - 1.0 * index / self.reverse_length
			k = (1 + 3 * Prob_True(r_portion))/4

			false_negative = 1 - (p * k + (1-p)*(1-k)/3)
			negative = 	1 - ((p + k - 2*p*k)/3 + (1-p)*(1-k) * 2/9)
			return false_negative / (negative + false_negative)


	def Set_Forward(self, forward_string):
		'''
		Forward Read to compare against
		'''
		self.forward_string = forward_string

	def Set_Lengths(self, for_l, rev_l):
		'''
		Utilized for purposes of updating the probability.
		Lengths are very important to this data because of the portion aspect
		'''
		self.forward_length = for_l
		self.reverse_length = rev_l
		self.min_length = for_l if for_l < rev_l else rev_l

	


'''
Decision Tree:

Index out of Range:
No match will occur guaranteed, but the chance of this happening is index / (length)

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
	# Define a Random Genomic Sequence owith equally weighted A, T, C, G
	genomic_sequence = Random_Genome([10000, 12000], [1, 1, 1, 1])

	# Get the probabilistic data for each member of the sequence based on reading errors
	forward_data, reverse_data, indices = Read_Genome(genomic_sequence, 
		[s.Get('length-min'), s.Get('length-max')], 
		[s.Get('overlap-min'), s.Get('overlap-max')])


	# Use suite to find the most probable overlap between the two pieces of data
	
	# Theoretically, the overlap can occur at any point in the sequence
	hypos = range(int(0.5 * s.Get('overlap-min') * forward_data.Length()), 
		int((1 - s.Get("overlap-min")) * forward_data.Length()))
	suite = Genome_Matcher(hypos)
	suite.Set_Forward(forward_data.Read())
	reverse_read = reverse_data.Read()
	suite.Set_Lengths(len(forward_data.Read()), len(reverse_read))

	# Deconstructs the reverse read string into a dataset and updates
	# Randomizes the order which data is updated to make it interesting
	# Limits the number of runs it does if a setting is made to do so
	dataset = [[idx, reverse_read[idx]] for idx in range(len(reverse_read))]
	random.shuffle(dataset)
	if(s.Get("runs-max") > 1):
		print("Only performing a limited number of tests")
		dataset = dataset[0:int(s.Get("runs-max"))]


	# Finds an unused logfile name by incrementing through the integers.
	val = 0
	while(os.path.isfile("log" + str(val) + ".txt") or 
		  os.path.isfile("log" + str(val) + "summary.txt")):
		val = val + 1

	f = open("log" + str(val) + ".txt", 'w')

	f.write("Arguments: " + ' '.join(sys.argv[1:]) + "\n")
	f.write("History\n")

	# Iterate through the dataset and log the data to be used in graphical analysis by anim.py
	suite.UpdateSetAndLog(dataset, f)

	# Finds a different index based on the most starkingly different value from the line.
	# For some unexplainable reason, this tends to get the actual value much more often
	ex = np.asarray(list(suite.Values()))
	why = np.log([suite[x] for x in suite.Values()])
	slope, intercept, r_value, p_value, std_err = stats.linregress(ex, why)
	adjusted = why - slope * ex
	adjusted_shift = list(adjusted).index(max(adjusted)) + suite.Values()[0]


	probs = suite.MostProbable()[0:5]
	for prob in probs:
		print(str(prob[0]) + ": " + str(100 *prob[1])[0:5] + "%")

	# This is the best shift known by probabilities
	shift = probs[0][0]
	f.write("End History\n")
	f.write(str(indices[2] - indices[0]) + "," + str(shift) + "," + str(adjusted_shift))
	f.close()

	print("True Shift: " + str(indices[2] - indices[0]) + " Most Likely: " + str(shift)
		+ " Adjusted Shift: " + str(adjusted_high))
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

	print("Adjusted Overlap")
	print(forward_data.Read())
	print(adjusted_high * "." + reverse_data.Read())

	print(reconstructed)

	f = open("log" + str(val) + "summary.txt", 'w')
	f.write("Arguments: " + ' '.join(sys.argv[1:]))
	f.write("\nForward Read: " + forward_data.Read())
	f.write("\nForward True: " + forward_data.Read_True())
	f.write("\nReverse Read: " + reverse_data.Read())
	f.write("\nReverse True: " + reverse_data.Read_True())
	f.write("\n")

	for prob in probs:
		f.write(str(prob[0]) + ": " + str(100 *prob[1])[0:5] + "%")
	f.write("\nTrue Shift: " + str(indices[2] - indices[0]) + " Most Likely: " + str(shift)
		+ " Adjusted Shift: " + str(adjusted_high))
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
