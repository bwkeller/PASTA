# Stacky is the class in the PAST suite that is used to generate stack objects.
# Stack statistics can be generated via stacky methods.
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

from math import sqrt

import numpy as np
from scipy import stats

from stampy import Stampy


class Stacky:
	"""
	This class is a datatype for stacks of fits stamps, as well as 
	operations on those stacks.
	"""
	def __init__(self):
		"""
		The constructor of this class merely fills the data with empty values
		"""
		self.data = None
		self.units = None
		self.rejected = 0
		self.is_sorted = False

	def add_stamp(self, stamp, nonoise=True):
		"""
		This method adds another stamp to the stack.
		@type stamp:	Stampy
		@param stamp: The stamp to be added to the stamp stack
		"""
		if stamp.header["flag"]["error_flag"] == 'a':
			#Check to see if the stack is empty
			if self.data == None:
				self.data = stamp.data
				self.units = stamp.header["units"]
				return 'a'
			#Check to see if the stack dimensions match.
			elif self.units != stamp.header["units"]:
				return 'd'
			else:
				try:
					self.data = np.dstack((self.data, stamp.data))
					return 'a'
				except ValueError:
					return 'd'
		else:
			self.rejected += 1
			return stamp.header["flag"]["error_flag"]

	def mean(self):
		"""
		This method returns a stamp of the mean values of the 
		stack.
		@rtype:		numpy array
		@return:	The numpy array containing the mean pixels of the stack
		"""
		#Check to see if there is only on stamp on the stack
		if len(self.data.shape) == 2:
			return self.data
		mean = np.average(self.data, axis=-1)
		return mean

	def median(self):
		"""
		This method returns the a stamp of the median values of the
		stack.
		@rtype:		numpy array
		@return:	The numpy array containing the median pixels of the stack
		"""
		#Check to see if there is only on stamp on the stack
		if len(self.data.shape) == 2:
			return self.data
		avg = np.apply_along_axis(np.median, -1, self.data)
		return avg

	def percentile_centralpix(self, percent):
		"""
		This method is used to find the value of the central pixel in the image
		at a certain percentile along the stack.
		@type percent:		number
		@param percent:		The percentile value to be found
		"""
		#Check to see if there is only on stamp on the stack
		if len(self.data.shape) == 2:
			return self.data[self.data.shape[0]/2.0][self.data.shape[1]/2.0]
		#Pull out the column of pixels at the center
		centralcol = self.data[self.data.shape[0]/2.0][self.data.shape[1]/2.0]
		#Sort the column, and then use linear interpolation to find the value at
		#the given percentile
		centralcol.sort()
		pixval = stats.scoreatpercentile(centralcol, percent)
		return pixval

	def percentile(self, percent):
		"""
		This method returns the a stamp at one percentile value of the
		stack.
		@type percent:	number
		@param percent:	percentile value to be calculated
		@rtype:			numpy array
		@return:		The numpy array containing the median pixels of the 
		stack
		"""
		#Check to see if there is only on stamp on the stack
		if len(self.data.shape) == 2:
			return self.data
		if not self.is_sorted:
			self.sort()
		percentile = np.apply_along_axis(stats.scoreatpercentile, -1, self.data,
		percent)
		return percentile

	def sort(self):
		"""
		This method is used to sort the stack along the z-axis (each pixel 
		column)
		"""
		#Check to see if there is only on stamp on the stack
		if len(self.data.shape) == 2:
			self.is_sorted = True
		else:
			self.data = np.sort(self.data, -1, kind='mergesort')
			self.is_sorted = True
