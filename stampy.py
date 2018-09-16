# Stampy is the class in the PAST suite used to store image data and metadata
# for individual image slices within a stack.  The stack is built from numerous
# stampy objects.
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

import math

import numpy as np
from scipy import ndimage
from math import sqrt

class Heady(dict):
        """
        This class is used a dictionary used for storing the various header flag
        values for a given stamp.
        """
        def __init__(self, name):
                """
                Generate a new header for a stamp.
                @type name:             string
                @param name:    the name for the stamp (generally the name of the 
                                                centered source)
                """
                self["source"] = name #Source Name
                self["stokes"] = '?' #Stokes parameter of the stamp
                self["radec"] = (0,0) #The RA and Declination of the stamp center
                self["noise"] = 0.0 #The rms of the border
                self["galcoords"] = (0,0) #The galactic coordinates of the stamp center
                #The flag dictionary is used to store warnings/errors
                #-has_nans is set if the source has 1 or more NaN pixels
                #-dimension_error is set if the generated stamp dimensions do not match 
                #the dimensions requested (this can occur near the borders of an image)
                #-high_noise is set if the rms in the outer border is larger than the
                #maximum specified noise
                #-near_source is set if there is a likely detection of a nearby source
                #-error_flag is the letter code for the source rejection
                self["flag"] = {"has_nans":False, "dimension_error":False, 
                "high_noise":False,     "near_source":False, "error_flag":'a'}
                self["units"] = "undefined" #The intensity units of the stamp image
                self["bmean"] = 0.0 #The mean intensity of the border
                self["bmedian"] = 0.0 #The median intensity of the border
                self["offset"] = [0,0] #The pixel diff from the center in RA/Dec vs. pix

class Stampy:
        """
        This class is the definition for the image stamp.
        @type dimensions:       number
        @param dimensions:      the length in pixels of the square stamp
        @type fitsfile:         pyfits object
        @param fitsfile:        the FITS file to pull the stamp from
        @type location:         (number, number)
        @param location:        The tuple containing the RA and Declination of the 
        source (in decimal degrees)
        @type name:                     string
        @param name:            The name of the source
        @type real_units:       boolean
        @param real_units:      whether or not to use real-valued units for each pixel
        @type epoch:            string
        @param epoch:           The epoch to use for calculating coordinates 
                                                (B1950 or J2000)
        @type maxnoise:         number
        @param maxnoise:        The maximum noise before the noise flag is flipped.
        @type fitswcs:          AstWCS.WCS
        @param fitswcs:         The WCS for the input fitsfile
        @type savewcs:          boolean
        @param savewcs:         Whether or not to generate a new WCS for the stamp
        """
        def __init__(self, dimensions, fitsfile, location, name="undefined", 
                                real_units=True, epoch="J2000", maxnoise=float('inf'), 
                                fitswcs=None, savewcs=False, galcoords=False):
                self.header = Heady(name)
                #pull out the FITS header from the given file
                self.fitsheader = fitsfile[0].header
                self.header['units'] = self.fitsheader['BUNIT']
                #The reason for the fitswcs argument is the astWCS.WCS method is
                #VERY SLOW.
                #generate a new WCS for the stamp if one isn't given
                if fitswcs is None:
                        fitswcs = wcs.WCS(self.fitsheader, mode="pyfits")
                self.header["radec"] = location
                if galcoords:
                        self.header["galcoords"] = astCoords.convertCoords(epoch, 
                        "GALACTIC", location[0], location[1], fitswcs.getEpoch())
                #This is the ndarray to pull the stamp from
                source = fitsfile[0].data
                #This check is needed for some odd fits files
                if len(fitsfile[0].data.shape) == 4:
                        source = fitsfile[0].data[0][0]
                #pixcentre is the central pixel of the source
                pixcentre = fitswcs.wcs_world2pix(np.atleast_2d(location), 0, ra_dec_order=True)[0]
                print(location, pixcentre)
                if pixcentre[0] < 0 or pixcentre[1] < 0 or np.any(np.isnan(pixcentre)):
                        self.data = np.array([])
                else:
                        self.header["offset"] = [pixcentre[0] -
                            int(round(pixcentre[0])),  pixcentre[1] -
                            int(round(pixcentre[1]))]
                        self.data = source[int(pixcentre[0]-dimensions/2):int(pixcentre[0]+dimensions/2), 
                                int(pixcentre[1]-dimensions/2):int(pixcentre[1]+dimensions/2)]
                #Once again, the astWCS methods are quite slow, disable them by default.
                if savewcs == True:     #generate a new WCS
                        wcsdimensions = dimensions * fitswcs.getPixelSizeDeg()
                        clip = astImages.clipImageSectionWCS(source, fitswcs,
                                location[0], location[1], wcsdimensions)
                        self.fitswcs = clip['wcs']
                else:
                        self.fitswcs = None
                #compute_rms does not work on 0 sized images
                #Check the dimension error flag
                if self.data.shape != (dimensions, dimensions):
                        self.header["flag"]["dimension_error"] = True
                        self.header["flag"]["error_flag"] = 'd'
                else:
                        self.header["bmean"] = self.compute_bmean()
                        self.header["noise"] = self.compute_rms()
                        self.header["bmedian"] = self.compute_bmedian()
                #We currently use bmedian rather than rms in order to 
                #avoid discarding stamps with nearby sources
                #The factor of 0.674443 is used to convert a maximum value for the rms
                #into a maximum value for the median border. (50% of the gaussian falls
                #between -0.674443 and 0.674443 sigma.)
                if self.header["bmedian"] > 0.674443*maxnoise:
                        self.header["flag"]["high_noise"] = True
                        self.header["flag"]["error_flag"] = 'n'
                #This likely means that there is a nearby source, 
                #as the source won't majorly affect the median
                if self.header["noise"] > 0.674443*maxnoise and not self.header["flag"]["high_noise"]:
                        self.header["flag"]["near_source"] = True
                if np.isnan(self.data).any():   #Check for blank pixels
                        self.header["flag"]["has_nans"] = True
                        self.header["flag"]["error_flag"] = 'b'

        def peak(self):
                """
                This returns the peak pixel value of the stamp
                @return:        A tuple containing the maximum pixel 
                                        value and a tuple of its location
                """
                #Check for empty or blanked stamps
                if self.data.size < 2 or self.header["flag"]["has_nans"]: 
                        return (0, (0,0))
                max = self.data.max()
                return (max, (np.where(self.data == max )[0][0], np.where(self.data ==
                                max)[1][0]))
        def min(self):
                """
                This returns the minimum pixel value of the stamp
                @return:        A tuple containing the minimum pixel 
                                        value and a tuple of its location
                """
                if self.data.size < 2 or self.header["flag"]["has_nans"]: 
                        return (0, (0,0))
                min = self.data.min()
                return (min, (np.where(self.data == min)[0][0], np.where(self.data == 
                                min)[1][0]))

        def save(self, filename):
                """
                This saves a fits file based on the internal stamp of the object
                @type filename:         string
                @param filename:        The filename to save the FITS file as
                """
                sliceWCS = self.fitswcs
                sliceWCS.header.update("DATAMIN", np.min(self.data))
                sliceWCS.header.update("DATAMAX", np.max(self.data))
                sliceWCS.updateFromHeader()
                astImages.saveFITS(filename + ".fits", self.data, sliceWCS)

        def compute_bmedian(self, bpercent=20):
                """
                This is an private method used to compute the 
                median noise of an outer "frame" 5+ pixels thick.
                @type bpercent:         int
                @param bpercent:        The percentage of the total image 
                                                        size to make the frame for
                @return:                The value of the median noise 
                                                computed over the outer frame.
                """
                thickness = 5   #Default frame thickness
                #Enlarge the thickness for bigger images
                if self.data.shape[0] * bpercent / 100 > thickness:
                        thickness = int(self.data.shape[0] * bpercent / 100)
                #This computes the sum of the left edge
                am = self.data[:,:thickness].ravel()
                # The upper edge
                bm = self.data[:thickness,thickness:(self.data.shape[1] - thickness)].ravel()
                # The right edge
                cm = self.data[(self.data.shape[0] - thickness):, 
                                thickness:(self.data.shape[1] - thickness)].ravel()
                # The bottom edge
                dm = self.data[:,(self.data.shape[0] - thickness):].ravel()
                return np.median(np.abs(np.concatenate((am, bm, cm, dm))))

        def compute_bmean(self, bpercent=20):
                """
                This is an private method used to compute the 
                mean noise of an outer "frame" 5 pixels thick.
                @type bpercent:         int
                @param bpercent:        The percentage of the total image 
                                                        size to make the frame for
                @return:                The value of the mean noise 
                                                computed over the outer frame.
                """
                thickness = 5   #Default frame thickness
                #Enlarge the thickness for bigger images
                if self.data.shape[0] * bpercent / 100 > thickness:
                        thickness = int(self.data.shape[0] * bpercent / 100)
                #calculate the size of the box
                bsize = (thickness * 2 * self.data.shape[0]) + 2 * \
                                ((self.data.shape[1] - 2 * thickness) * thickness)
                #This computes the sum of the left edge
                am = self.data[:,:thickness].ravel()
                #The upper edge
                bm = self.data[:thickness,thickness:(self.data.shape[1] - thickness)].ravel()
                #The right edge
                cm = self.data[(self.data.shape[0] - thickness):, 
                thickness:(self.data.shape[1] - thickness)].ravel()
                #The bottom edge
                dm = self.data[:,(self.data.shape[0] - thickness):].ravel()
                return np.sum(np.concatenate((am, bm, cm, dm))) / bsize

        def compute_rms(self, bpercent=20, mean=0):
                """
                This is an private method used to compute the 
                rms of an outer "frame" 5 pixels thick.
                @type bpercent:         int
                @param bpercent:        The percentage of the total image 
                                                        size to make the frame for
                @type mean:                     number
                @param mean:            Pass the mean of the box if you want to calculate 
                                                        the standard deviation rather than the rms
                @return:                        The value of the rms noise 
                                                        computed over the outer frame.
                """
                thickness = 5   #Default frame thickness
                #Enlarge the thickness for bigger images
                if self.data.shape[0] * bpercent / 100 > thickness:
                        thickness = int(self.data.shape[0] * bpercent / 100)
                #calculate the size of the box
                bsize = (thickness * 2 * self.data.shape[0]) + 2 * \
                                ((self.data.shape[1] - 2 * thickness) * thickness)
                #This computes the square sum of the left edge
                a = np.power((self.data[:,:thickness] - mean), 2).ravel()
                #The upper edge
                b = np.power((self.data[:thickness, thickness:(self.data.shape[1] - thickness)] - mean), 2).ravel()
                #The right edge
                c = np.power((self.data[:,:thickness] - mean), 2).ravel()
                #The bottom edge
                d = np.power((self.data[:,self.data.shape[0] - thickness:] - mean), 2).ravel()
                try:
                        return np.sqrt(np.sum(np.concatenate((a, b, c, d))) / (bsize - 1))
                except: 
                        print("Error with source @ "+self.header["radec"])

        def regrid(self, factor):
                """
                This method is used to regrid the image to allow >1 pixel accuracy
                @type factor:           number
                @param factor:          the factor to scale the image by
                """
                if not self.header["flag"]["dimension_error"] and not self.header["flag"]["has_nans"]: 
                        self.scale(factor)
                        #calculate the pixel distance to the true center
                        offset = list(map(lambda x: int(round(-x*factor)), self.header['offset']))
                        #shift the image around based on the true center
                        self.data = np.roll(self.data, offset[0], axis=1)
                        self.data = np.roll(self.data, offset[1], axis=0)
                        self.scale(1. / factor)
                        self.data = self.data[1:-1, 1:-1]
                        self.header["bmean"] = self.compute_bmean()
                        self.header["noise"] = self.compute_rms()
                        self.header["bmedian"] = self.compute_bmedian()


        def scale(self, factor):
                """
                This method is used to scale the image stamp using cubic splines in astLib.
                @type factor:           number or (number,number)
                @param factor:          the factor to scale the image by
                """
                if type(factor) == int or type(factor) == float:
                        factor = [factor, factor]       
                self.data = ndimage.zoom(self.data, factor)
                
        def merge(self, stamp):
                """
                This method merges a Q and U stamp into one new stamp, by adding the
                two stamps in quadrature.
                @type stamp:    stampy
                @param stamp:   the stamp to merge into this one
                """
                if stamp.data.shape != self.data.shape:
                        raise ArithmeticException
                self.data = np.add(np.power(self.data, 2), np.power(stamp.data, 2))
                self.data = np.sqrt(self.data)
                self.header["stokes"] = "PI"

        def add(self, value):
                """
                This method adds a constant value offset to each pixel in the stamp.
                @type value:    number
                @param value:   The value to be added to the stamp
                """
                mask = np.ones(self.data.shape)*value
                self.data = self.data + mask

        def multiply(self, value):
                """
                This method multiplies each pixel in the stamp by a constant value.
                @type value:    number
                @param value:   The value to be multiplied by the stamp
                """
                mask = np.ones(self.data.shape)*value
                self.data = self.data * mask

        def add_fakesource(self, peak, FWHM):
                """
                This method adds a seeded 2D gaussian with a specified peak value and 
                FWHM to the position in the source list.
                @type peak:             number
                @param peak:    The peak value for the seeded gaussian
                @type FWHM:             number
                @param FWHM:    The FWHM of the seeded gaussian
                """
                sigma = FWHM/2.35482
                mask = np.zeros(self.data.shape)
                mask[self.data.shape[0]/2, self.data.shape[1]/2] = 1.0
                mask = ndimage.filters.gaussian_filter(mask, sigma, mode='constant')
                mask = ndimage.interpolation.shift(mask, (-1.0*self.header["offset"][0], -1.0*self.header["offset"][1]))
                mask = peak*mask/mask[self.data.shape[0]/2, self.data.shape[1]/2]
                self.data += mask
