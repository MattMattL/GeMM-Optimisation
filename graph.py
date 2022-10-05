import matplotlib.pyplot as plt
import re

# Read results from file and decompose
file = open("results.txt", "r")

sizes = []
timeTaken = []
sd = []

for line in file:
	numbers = re.findall("\d+\.?\d+", line)

	sizes.append(float(numbers[0]))
	timeTaken.append(float(numbers[1]))
	sd.append(float(numbers[2]))

# Plot results
figure, plotAxes = plt.subplots()
plotAxes.plot(sizes, timeTaken, linestyle='', marker='o', color='red')

plt.show()

