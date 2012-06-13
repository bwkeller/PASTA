#!/usr/bin/python

# Readheader is a python script that uses pyfits to read and print stack FITS 
# header values from a stacked fits image.  Readheader is part of PAST.
#
# Copyright (C) <2009>  <Ben Keller (bwkeller@ucalgary.ca)> 
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
# See doc/LICENSE for full text of the GNU GPL

import sys
from optparse import OptionParser 
from math import sqrt

import pyfits
import numpy as np

def compute_rms(n, bpercent=20):
		"""
		This is an private method used to compute the standard deviation of an 
		outer "frame" 5 pixels thick.
		@type n:			ndarray
		@param n:			The array to compute the frame rms noise for.
		@type bpercent:		int
		@param bpercent:	The percentage of the total image size to make the 
							frame for
		@return:			The value of the rms noise computed over the outer 
							frame.
		"""
		#Default frame thickness
		thickness = 5
		#Enlarge the thickness for bigger images
		if n.shape[0]*bpercent/100 > thickness:
			thickness = n.shape[0]*bpercent/100
		#calculate the size of the box
		bsize = (thickness*2*n.shape[0])+2*((n.shape[1]-2*thickness)*thickness)
		mean = sum(n[:,:thickness].flatten()) + \
		sum(n[:thickness,thickness:(n.shape[1]-thickness)].flatten()) + \
		sum(n[(n.shape[0]-thickness):, 
		thickness:(n.shape[1]-thickness)].flatten()) + \
		sum(n[:,(n.shape[0]-thickness):].flatten())
		mean = mean / (bsize - 1)
		#This computes the square sum of the left edge
		a = sum(np.power((n[:,:thickness]-mean), 2).flatten())
		# The upper edge
		b = sum(np.power((n[:thickness, thickness:(n.shape[1]-thickness)]-mean),
		2).flatten())
		# The right edge
		c = sum(np.power((n[:,:thickness]-mean), 2).flatten())
		# The bottom edge
		d = sum(np.power((n[:,n.shape[0]-thickness:]-mean), 2).flatten())
		return sqrt((a+b+c+d)/(bsize-1))

def compute_bmedian(n, bpercent=20):
		"""
		This is an private method used to compute the median noise of an outer 
		"frame" 5+ pixels thick.
		@type n:			ndarray
		@param n:			The array to compute the frame median noise for.
		@type bpercent:		int
		@param bpercent:	The percentage of the total image size to make the 
							frame for
		@return:			The value of the median noise computed over the									outer frame.
		"""
		#Default frame thickness
		thickness = 5
		#Enlarge the thickness for bigger images
		if n.shape[0]*bpercent/100 > thickness:
			thickness = n.shape[0]*bpercent/100
		#calculate the size of the box
		bsize = (thickness*2*n.shape[0])+2*((n.shape[1]-2*thickness)*thickness)
		#This computes the sum of the left edge
		am = n[:,:thickness].flatten()
		# The upper edge
		bm = n[:thickness,thickness:(n.shape[1]-thickness)].flatten()
		# The right edge
		cm = n[(n.shape[0]-thickness):, 
		thickness:(n.shape[1]-thickness)].flatten()
		# The bottom edge
		dm = n[:,(n.shape[0]-thickness):].flatten()
		return np.median(np.abs(np.concatenate((am, bm, cm, dm))))

if __name__ == "__main__":
	parser = OptionParser("readheader.py [options] infile")
	parser.add_option("-t", action="store_true", dest="tabular", 
	help="Output data in single line format for generating tables")
	(options, args) = parser.parse_args()
	if len(args) > 0:
		infile = pyfits.open(args[0])
	else:
		print "Error Reading input file"
		sys.exit(1)
	inheader = infile[0].header
	units = inheader["BUNIT"]
	min = inheader["STACKMIN"]
	max = inheader["STACKMAX"]
	median = inheader["STACK50"]
	num = inheader["STACKNUM"]
	lowquartile = inheader["STACK25"]
	upquartile = inheader["STACK75"]
	dimensions = infile[0].data.shape
	rms = compute_rms(infile[0].data)
	bmedian = compute_bmedian(infile[0].data)
	maxp = np.maximum.reduce(np.maximum.reduce(infile[0].data))
	minp = np.minimum.reduce(np.minimum.reduce(infile[0].data))
	maxloc = (np.where(infile[0].data == maxp )[0][0], np.where(infile[0].data 
	==	maxp)[1][0])
	minloc = (np.where(infile[0].data == minp )[0][0], np.where(infile[0].data 
	== minp)[1][0])

	if options.tabular:
		print "%(min)12.5E %(max)12.5E %(median)12.5E %(low)12.5E %(up)12.5E %(rms)12.5E %(bmedian)12.5E %(maxp)12.5E %(maxx)3d %(maxy)3d %(minp)12.5E %(minx)3d %(miny)3d %(x)3d %(y)3d %(num)7d %(name)s" % {'name':args[0], 'min':min, 
		'max':max, 'median':median, 'low':lowquartile, 'up':upquartile, 
		'rms':rms, 'bmedian':bmedian, 'maxp':maxp, 'maxx':maxloc[0], 
		'maxy':maxloc[1], 'minp':minp, 'minx':minloc[0], 'miny':minloc[1], 
		'x':dimensions[0], 'y':dimensions[1], 'num':num }
	else:
		print args[0]
		print "\n"
		print "Source Intensity Data"
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
		print "Minimum Stack Intensity:\t"+str(min)+" "+units
		print "Maximum Stack Intensity:\t"+str(max)+" "+units
		print "Median Stack Intensity:\t\t"+str(median)+" "+units
		print "Lowest Quartile Intensity:\t"+str(lowquartile)+" "+units
		print "Highest Quartile Intensity:\t"+str(upquartile)+" "+units
		print "Image Data"
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
		print "Edge RMS Noise:\t\t\t"+str(rms)+" "+units
		print "Edge Median Noise:\t\t"+str(bmedian)+" "+units
		print "Maximum Pixel Intensity:\t"+str(maxp)+" "+units
		print "Minimum Pixel Intensity:\t"+str(minp)+" "+units
		print "Maximum Pixel Location:\t\t"+str(maxloc)
		print "Minimum Pixel Location:\t\t"+str(minloc)
		print "Other Data"
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
		print "Image Dimensions:\t\t"+str(dimensions)
		print "Number of sources:\t\t"+str(num)+"\n"
