import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 1, 101)
y_c1_s1 = 1 - 1 * x ** 1
y_c2_s1 = 1 - 1 * x ** 2
y_c3_s1 = 1 - 1 * x ** 3
y_c4_s1 = 1 - 1 * x ** 4

plt.plot(x, y_c1_s1, 'r')
plt.plot(x, y_c2_s1, 'b')
plt.plot(x, y_c3_s1, 'g')
plt.plot(x, y_c4_s1, 'k')

y_c1_s2 = 1 - 0.8 * x ** 1
y_c2_s2 = 1 - 0.8 * x ** 2
y_c3_s2 = 1 - 0.8 * x ** 3
y_c4_s2 = 1 - 0.8 * x ** 4

plt.plot(x, y_c1_s2, 'ro')
plt.plot(x, y_c2_s2, 'bo')
plt.plot(x, y_c3_s2, 'go')
plt.plot(x, y_c4_s2, 'ko')

plt.legend(["c = 1 s = 1", "c = 2 s = 1", "c = 3 s = 1", "c = 4 s = 1",
		"c = 1 s = 0.8", "c = 2 s = 0.8", "c = 3 s = 0.8", "c = 4 s = 0.8"], "SouthWest")

plt.xlabel("Fraction through sequence")
plt.ylabel("Probability Correct")
plt.title("Probability that a Base-Pair is Read Correctly")

plt.show()