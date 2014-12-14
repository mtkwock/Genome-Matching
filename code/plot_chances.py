import numpy as np
import matplotlib.pyplot as plt
import sys

fname = sys.argv[1]

exes = [0] * 101
wise = [0] * 101

f = open(fname, 'r')

for i in range(101):
	split = f.readline().split(" ")
	exes[i] = split[0]
	wise[i] = split[1]

x = np.asarray(exes)
y = np.asarray(wise)

fig, ax = plt.subplots()

ax.set_xlabel("Fraction through data")
ax.set_ylabel("Probability of correctness")
tit = fname[:-4].split(" ") # remove .txt
inte = "S_{err}= |, \sigma_{err}= |, S_{read}= |, \sigma_{err}= |, S_{curve}= ".split("|")

title = "$" + inte[0] + tit[0] + inte[1] + tit[1] + inte[2] + tit[2] + inte[3] + tit[3] + inte[4] + tit[4] + "$"

ax.set_title(title)

line, = ax.plot(x, y)

plt.show()