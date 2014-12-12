"""
A simple example of an animated plot
Obtained from: http://matplotlib.org/examples/animation/simple_anim.html and heavily modified
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import os

def normalize(data):
	total = sum(data)
	return [x / total for x in data]

fig, ax = plt.subplots()
update_text = ax.text(0.02, 0.90, '', transform=ax.transAxes)
ax.set_yscale('log')
ax.set_xlabel('Shift amount')
ax.set_ylabel('Probability')

assert len(sys.argv) > 1, "No filename specified"
assert os.path.isfile(sys.argv[1]), "File does not exist"

f = open(sys.argv[1], 'r')
lines = f.read().split('\n')
history_str = lines[lines.index("History:") + 1:]
xy = [line.split('],') for line in history_str]
for y in range(len(xy)):
	for x in range(len(xy[y])):
		xy[y][x] = ''.join([c for c in xy[y][x] if c not in "[ ]"])

data = []
annotations = []
for time in xy:
	if(len(time) == len(xy[1])):
		data = data + [[[int(v.split(",")[0]), float(v.split(",")[1])] for v in time if len(v.split(",")) > 1]]
		annotations = annotations + [[[v.split(",")[2], v.split(",")[3]] for v in time if len(v.split(",")) > 1]]
	else:
		print(len(time) - len(xy[1]))

# x = np.arange(0, 2*np.pi, 0.01)        # x-array
x = np.asarray([row[0] for row in data[0]])
y = np.asarray(normalize([row[1] for row in data[0] if len(row) > 1]))
line, = ax.plot(x, y)


def animate(i):
	'''
	Iterated function to animate the plot.
	Changes the y-values as well as the text in the top left
	'''
	value = normalize([row[1] for row in data[i] if len(row) > 1])
	line.set_ydata(np.asarray(normalize(value)))  # update the data
	update_text.set_text("Updates: " + str(i+1) + "\nData: " + annotations[i][0][0] + " " + annotations[i][0][1])
	# annotation = plt.annotate('Data Number: ' + str(i), xy=("a","b"))
	return line, update_text

#Init only required for blitting to give a clean slate.
def init():
    line.set_ydata(len(x) * [0])
    return line,

delay = 10000 / len(data)

ani = animation.FuncAnimation(fig, animate, np.arange(0, len(data)-1), init_func=init,
    interval=delay, blit=True)
plt.ylim([10**(-1 * len(data)), 1])
plt.show()