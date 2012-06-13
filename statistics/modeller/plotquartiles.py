#!/usr/bin/python
#This tool is used to plot the quartile ratio of the input stacks against the
#quartile ratio of the simulated stacks

from optparse import OptionParser

import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
	usage = "plotquartiles.py statistics_file"
	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()
	file = open(args[-1]).readlines()
	x = [np.log10(float(i.split()[2])) for i in file]
	y = [float(i.split()[5])/float(i.split()[4]) for i in file]
	P25err = np.array([float(i.split()[16]) for i in file])
	P25 = np.array([float(i.split()[14]) for i in file])
	P75 = np.array([float(i.split()[15]) for i in file])
	P75err = np.array([float(i.split()[17]) for i in file])
	ysim = P75/P25
	ysimerr = P75err/P25 + P25err*P75/(P25*P25)
	plt.plot(x, y, label="Real Quartile Fraction")
	plt.errorbar(x, ysim, yerr=ysimerr, label="Simulated Quartile Fraction")
	plt.xlabel("log10 Stokes I Intensity")
	plt.legend()
	plt.show()
