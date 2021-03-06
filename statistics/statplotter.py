#!/usr/bin/python
#This is a matplotlib python script to generate a full statisitics plot using
#the statistics table.

import sys
from optparse import OptionParser
import numpy as np
import pyfits
from matplotlib.ticker import *
from pylab import figure,show, setp
from matplotlib import rc

#These lists are used for storing the statistics values read in from the table
bins = []
I = []
P = []
P25 = []
P75 = []
N = []
mean_Poff = []
sig_Poff = []
model_Poff = []
sig_mode = []
P0_ML = []
P0_16 = []
P0_83 = []
P_median = []
Percent = []
stackerr = []
modelerr = []
P2 = []
P3 = []
P4 = []
hline = []


if __name__ == "__main__":
	#rc('font',**{'size':20, 'family':'sans-serif','sans-serif':['Helvetica']})
	rc('text', usetex=True)
	parser = OptionParser("statplotter.py [options] infile")
	(options, args) = parser.parse_args()
	table = open(args[0]).readlines()
	linedict = {}
	for line in table:
		linedict[float(line.split()[2])] = line
	Ilist = linedict.keys()
	Ilist.sort()
	#This filters the upper limit, as we are more interested in the lower fluxes
	Ilist = filter(lambda x: x<1100, Ilist)
	#This reads in the statistics values
	for Ival in Ilist:
		line = linedict[Ival]
		bins.append((float(line.split()[0]), float(line.split()[1])))
		I.append(float(line.split()[2]))
		P.append(float(line.split()[3]))
		P25.append(float(line.split()[4]))
		P75.append(float(line.split()[5]))
		N.append(float(line.split()[6]))
		mean_Poff.append(float(line.split()[7]))
		sig_Poff.append(float(line.split()[8]))
		model_Poff.append(float(line.split()[9]))
		sig_mode.append(float(line.split()[10]))
		P0_ML.append(float(line.split()[11]))
		P0_16.append(float(line.split()[12]))
		P0_83.append(float(line.split()[13]))
		P_median.append(float(line.split()[14]))
		#This percent is based on the most likely P0
		Percent.append(100*float(line.split()[11])/float(line.split()[2]))
		#This percent is based on the observed P
		P2.append(100*float(line.split()[3])/float(line.split()[2]))
		#This percent is based on the 16.5th percentile ML P0
		P3.append(100*float(line.split()[12])/float(line.split()[2]))
		#This percent is based on the 83.5th percentile ML P0
		P4.append(100*float(line.split()[13])/float(line.split()[2]))
		stackerr.append(float(line.split()[8])/(float(line.split()[6])**2))
		modelerr.append(float(line.split()[10])/(float(line.split()[6])**2))
		#This line is a constant percentage polarization
		#hline.append(2.097906)
		hline.append(2.0)
	(xmin, xmax) = (bins[0][0], bins[-1][1])
	fig = figure()
	yprops = dict(rotation=90, horizontalalignment='right', size='small',
	verticalalignment='center', x=-0.01)
	#axprops = dict(yticks=[])
	axprops = dict(yscale='log', xscale='log', autoscalex_on=False)
	#ax1 plots the median, 25th, and 75th percentile PI from the stack vs I
	ax1 = fig.add_axes([0.15, 0.85, 0.8, 0.15], **axprops)
	ax1.plot(I, P75, 'k--', I, P25, 'k--', I, P, 'ks-')
	ax1.set_ylabel(r'$P (mJy/beam)$', **yprops)
	ax1.axis('tight')
	ax1.xaxis.set_visible(False)
	for line in ax1.get_yticklines():
		line.set_markersize(8)
	for line in ax1.yaxis.get_minorticklines():
		line.set_markersize(5)
	axprops['yscale']='linear'
	#ax2 plots the model and stack background median vs I
	ax2 = fig.add_axes([0.15, 0.7, 0.8, 0.15], **axprops)
	#These lines overplot two error-barred lines on the same graph
	ax2.errorbar(I, mean_Poff, yerr=stackerr, fmt='ko-')
	ax2.errorbar(I, model_Poff, yerr=stackerr, fmt='k^--')
	ax2.set_ylabel(r'$offset (mJy/beam)$', **yprops)
	ax2.axis('tight')
	ax2.xaxis.set_visible(False)
	for line in ax2.get_yticklines():
		line.set_markersize(8)
	for line in ax2.yaxis.get_minorticklines():
		line.set_markersize(5)
	axprops['yscale']='log'
	#ax3 plots the model and stack rms vs I
	ax3 = fig.add_axes([0.15, 0.55, 0.8, 0.15], **axprops)
	ax3.plot(I, sig_Poff, 'ko-', I, sig_mode, 'k^--')
	ax3.axis('tight')
	ax3.xaxis.set_visible(False)
	ax3.set_ylabel(r'$\sigma_{PI}(mJy/beam)$', **yprops)
	axprops['yscale']='log'
	for line in ax3.get_yticklines():
		line.set_markersize(8)
	for line in ax3.yaxis.get_minorticklines():
		line.set_markersize(5)
	#ax4 plots the Number of sources vs I
	ax4 = fig.add_axes([0.15, 0.4, 0.8, 0.15], **axprops)
	ax4.plot(I, N, 'ko-')
	ax4.xaxis.set_visible(False)
	ax4.set_ylabel(r'$N_{stack}$', **yprops)
	ax4.axis('tight')
	for line in ax4.get_yticklines():
		line.set_markersize(8)
	for line in ax4.yaxis.get_minorticklines():
		line.set_markersize(5)
	axprops['yscale']='log'
	#ax5 plots stack P, simulation true P0, most likely P0, 16.5th and 83.5th 
	#percentile vs I
	ax5 = fig.add_axes([0.15, 0.25, 0.8, 0.15], **axprops)
	ax5.plot(I, P, color='0.8', marker='o')
	ax5.plot(I, P0_ML, 'ko-', I, P0_16, 'k--', I, P0_83, 'k--')
	#ax5.plot(I, P0_ML, 'ko-', I, P0_16, 'k--', I, P0_83, 'k--', I, P_median, 'r-')
	ax5.xaxis.set_visible(False)
	ax5.set_ylabel(r'$P_0(mJy/beam)$', **yprops)
	ax5.axis('tight')
	for line in ax5.get_yticklines():
		line.set_markersize(8)
	for line in ax5.yaxis.get_minorticklines():
		line.set_markersize(5)
	axprops['yscale']='linear'
	#ax6 plots the Percent polarization based on observed P, P0_ML, P0_16,
	#P0_83, and true simulation P0 vs I
	ax6 = fig.add_axes([0.15, 0.1, 0.8, 0.15], **axprops)
	ax6.plot(I, Percent, 'ko-', I, P2, I, P3, 'k--', I, P4, 'k--', I, hline, 'b--')
	ax6.axis('tight')
	ax6.set_ylabel(r'$\Pi_0', **yprops)
	ax6.set_xlabel(r'Stokes I (mJy/beam)')
	ax6.set_autoscale_on(False)
	ax6.set_ylim([0.0, 6])
	for line in ax6.get_yticklines():
		line.set_markersize(8)
	for line in ax6.yaxis.get_minorticklines():
		line.set_markersize(5)
	show()

