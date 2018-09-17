#!/usr/bin/python3

# genstack.py is the primary script of the PAST suite.  genstack takes an input
# sourcelist and one or more FITS files containing the sources, and generates
# a mean and median stacked FITS image set.
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
# See for doc/LICENSE full text of the GNU GPL

import math
from sys import stdout, exit
from os import path
from optparse import OptionParser

import numpy as np

from astropy import wcs
from astropy.io import fits

from stacky import Stacky
from stampy import Stampy

#These lists are used to store the data extracted from the sourcelist
RA = []
Dec = []
Intense = []

#These files are used to store the data files needed by the functions
Stack = Stacky()
noise = None
options = None
args = None
annotation = None
cdelt1 = -20. / 3600
cdelt2 = 20. / 3600
units = None

def extract_data(filename, ra, dec, intense):
    """
    This method pulls data in from a sourcelist and populates the RA, Dec, and
    Intense arrays.
    @type filename:         string
    @param filename:        filename of the source list
    @type ra:                       int
    @param ra:                      column # containing the Right Ascension
    @type dec:                      int
    @param dec:                     column # containing the Declination
    @type intense:          int
    @param intense:         column # containing the Stokes I intensity
    """
    galaxies = open(filename).readlines()
    for i in galaxies:
        # Unrejected is used if the sourcelist is also a noisefile, to filter
        # out rejected sources (which may not be rejected in a different band)
        unrejected = True
        if options.noisefile:
            # The last column is the rejection flag
            if i.split()[12] != 'a':
                unrejected = False
        if unrejected:
            RA.append(float(i.split()[ra]))
            Dec.append(float(i.split()[dec]))
            Intense.append(i.split()[intense])
    print("Done reading infile")

def name_string(ra, dec):
    """
    This method takes a decimal degree value of RA and Dec and generates a
    name string for the source for printing out
    @type ra:               float
    @param ra:              Right Ascension
    @tpye dec:              float
    @param dec:             Declination
    @returns:               string containing the name of the source
    """
    while ra > 360.0:
        ra = ra - 360.0
    while ra < 0.0:
        ra = ra + 360.0
    #Generate the RA section of the source name
    name = "J%(RH)02d%(RM)02d%(RS)02d" % {'RH':math.floor(ra / 15), 'RM':(((ra / 15)
    - math.floor(ra / 15)) * 60),'RS':(((((ra / 15) - math.floor(ra / 15)) * 60)
    - math.floor(((ra / 15) - math.floor(ra / 15)) * 60)) * 60)}
    #Add the Declination sign to the source name
    if dec > 0:
        name = name + "+"
    else:
        name = name + "-"
    dec = abs(dec)
    #Generate the Dec section of the source name
    name = name + "%(DD)02d%(DM)02d%(DS)02d" % {'DD':math.floor(dec),
    'DM':math.floor(60 * (dec - math.floor(dec))),
    'DS':math.floor(60 * ((60 * (dec-math.floor(dec))) -
    math.floor(60 * (dec - math.floor(dec))))) }
    return name

def stack_build(image, size, ra, dec, intense, mywcs, scale, save, maxnoise,
        subtract, fakesource):
    """
    This method loads a stamp from a fits file, and then adds that stamp to
    the stack.  It returns the rejection code character
    @type image:            pyfits.HDUList
    @param image:           pyfits source image to pull stamps from
    @type size:                     int
    @param size:            The size of the stamp to load
    @type ra:                       float
    @param ra:                      Right Ascension of source
    @type dec:                      float
    @param dec:                     Declination of source
    @type intense:          float
    @param intense:         Stokes I intensity of the source
    @type mywcs:                      WCS
    @param mywcs:                     WCS for the source image
    @type scale:            int or None
    @param scale:           The scale factor for regridding
    @type save:                     boolean
    @param save:            Save individual stamps
    @type maxnoise:         float
    @param maxnoise:        The maximum noise for the stamp
    @type subtract:         boolean
    @param subtract:        Subtract border median value from each stamp.
    @type fakesource:       string
    @param fakesource:      The values for the peak, FWHM, and distribution
                                            of a seeded gaussian
    @return:                        The character code of the reject flag
    """
    #This catches stamps already read
    if ra is None:
        return 'x'
    #Apply regridding if requested, and generate a stamp
    if scale is not None:
        #This check allows for noise filtering to be turned off
        if maxnoise is not None:
            stamp = Stampy(size+2, image, (ra, dec), name=name_string(ra, dec),
            fitswcs=mywcs, maxnoise=maxnoise, savewcs=save)
            stamp.regrid(scale)
        else:
            stamp = Stampy(size+2, image, (ra, dec), name=name_string(ra, dec),
            fitswcs=mywcs, savewcs=save)
            stamp.regrid(scale)
    else:
        #This check allows for noise filtering to be turned off
        if maxnoise is not None:
            stamp = Stampy(size, image, (ra, dec), name=name_string(ra, dec),
            fitswcs=mywcs, maxnoise=maxnoise, savewcs=save)
        else:
            stamp = Stampy(size, image, (ra, dec), name=name_string(ra, dec),
            fitswcs=mywcs, savewcs=save)
    if subtract:
        #subtract from to the stamp pixel the border median value
        stamp.add(-1.0*stamp.header["bmedian"])
    #Add a fake 2D Gaussian source to the stamp if selected and stamp accepted
    if fakesource is not None and stamp.header['flag']['error_flag'] == 'a':
        #The peaks will be uniformly distributed between 0 and 2*peak
        if fakesource.split()[-1] == 'u':
            stamp.add_fakesource(np.random.uniform(high =
            2.0*float(fakesource.split()[0])), float(fakesource.split()[1]))
        #The peaks will all be set to the input peak value (delta distribution)
        if fakesource.split()[-1] == 'd':
            stamp.add_fakesource(float(fakesource.split()[0]),
            float(fakesource.split()[1]))
    #Add the stamp to the stack
    reject = Stack.add_stamp(stamp)
    #call the noise file function for the stamp
    save_noise(noise, annotation, stamp, intense, reject)
    if save and reject != 'd' and reject != 'b':
        stamp.save(name_string(ra, dec))
    return reject

def headerize(image):
    """
    This method generates the FITS headers for the final median and mean
    images that are saved.
    @type image:    pyfits image
    @param image:   The FITS image of the stacked results
    @returns:               The pyfits image with the proper header values included
    """
    # These values are used due to a bug in the scale, so my own scale
    # bzero and bscale values must be internally calculated.
    scalingfactor = 0.000001
    maxintvalue = 4294967294
    # This adds the appropriate fits header values
    image[0].header.set('STACKMIN', Stack.percentile_centralpix(0),
                                                    'Minimum value of central pixel')
    image[0].header.set('STACKMAX', Stack.percentile_centralpix(100),
                                                    'Maximum value of central pixel')
    image[0].header.set('STACK50', Stack.percentile_centralpix(50),
                                                    'Median value of central pixel')
    image[0].header.set('STACK25', Stack.percentile_centralpix(25),
                                                    '25th Percentile value of central pixel')
    image[0].header.set('STACK75', Stack.percentile_centralpix(75),
                                                    '75th percentile value of central pixel')
    image[0].header.set('STACKNUM', Stack.data.shape[2],
                                                    'Number of accepted sources used')
    image[0].header.set('BUNIT', units, '')
    image[0].header.set('CRPIX1', 1, '')
    image[0].header.set('CRPIX2', 1, '')
    image[0].header.set('CTYPE1', 'RA-CAR', '')
    image[0].header.set('CTYPE2', 'DEC-CAR', '')
    # CDELT SHOULD be read from input fits file
    image[0].header.set('CDELT1', cdelt1, '')
    image[0].header.set('CDELT2', cdelt2, '')
    image[0].header.set('CRVAL1', -((image[0].data.shape[0]) / 2.0) *
                                                    image[0].header['CDELT1'], '')
    image[0].header.set('CRVAL2', -((image[0].data.shape[1]) / 2.0) *
                                                    image[0].header['CDELT2'], '')
    # These lines are used to scale the image, as some software can only handle
    # int values for the data array
    min = np.minimum.reduce(np.minimum.reduce(image[0].data))
    max = np.maximum.reduce(np.maximum.reduce(image[0].data))
    # These lines are needed to prevent the lowest intensity pixels from being
    # blanks
    min = min - (scalingfactor * abs(min))
    max = max + (scalingfactor * abs(max))
    _zero = (max + min) / 2.
    _scale = (max - min) / (4294967294)
    image[0].scale('int32', '', bzero=_zero, bscale=_scale)
    return image

def save_stacks(name):
    """
    This method saves the mean and median fits files for the stack
    @type name:             string
    @param name:    filename base for FITS images
    """
    meanf = headerize(fits.HDUList(fits.PrimaryHDU(Stack.mean())))
    meanf.writeto(name + "_mean.fits")
    medianf = headerize(fits.HDUList(fits.PrimaryHDU(Stack.median())))
    medianf.writeto(name + "_median.fits")

def save_percentile(name, percent):
    """
    This method saves stamps containing a specific percentile of the stack
    @type name:             string
    @param name:    filename base for FITS images
    @type percent:  number
    @param percent: the percentile to dump out
    """
    percentf = headerize(fits.HDUList(fits.PrimaryHDU( \
    Stack.percentile(percent))))
    percentf.writeto("%s_%8.6f_percent.fits" % (name, percent))


def save_noise(file, ann, stamp, intensity, reject):
    """
    This method generates and saves the noise file and annotation for the stack.
    @type file:                     file
    @param file:            the noise file
    @type ann:                      file
    @param ann:                     the annotation file
    @type stamp:            Stampy
    @param stamp:           the stamp containing the source
    @type intensity:        float
    @param intensity:       the Stokes I intensity for the source
    @type reject:           char
    @param reject:          the rejection flag for the source
    """
    # We don't want any NaNs in the noisefile, so if the central pixel cp is
    # NaN, we set it to 0
    try:
        cp = stamp.data[int(stamp.data.shape[0] / 2)][int(stamp.data.shape[1] / 2)]
        if np.isnan(cp):
            cp = 0
    except IndexError:
        cp = 0
    # Don't write data for stamps that come from outside the image
    if not reject == 'd':
        # Build the output string to write to the noise file
        outstring = "%(ra)3.5f %(dec)3.5f %(I)12.5E %(CP)12.5E %(PI)12.5E %(XPI)4d %(YPI)4d %(Min)12.5E %(XMin)4d %(YMin)4d %(median)12.5E %(noise)12.5E %(rej)s %(name)s\n" % {"ra":float(stamp.header["radec"][0]),
        "dec":float(stamp.header["radec"][1]), "I":float(intensity), "CP":cp,
        "PI":float(stamp.peak()[0]), "XPI":stamp.peak()[1][0],
        "YPI":stamp.peak()[1][1], "Min":float(stamp.min()[0]),
        "XMin":stamp.min()[1][0], "YMin":stamp.min()[1][1],
        "median":float(stamp.header["bmedian"]),
        "noise":float(stamp.header["noise"]), "rej":reject,
        "name":stamp.header['source']}
        file.write(outstring)
        # Build the annotation file
        if not options.annotate:
            # Accepted source
            if reject == 'a':
                colour = "WHITE"
            # Contains blank pixels
            if reject == 'b':
                colour = "YELLOW"
            # Stamp too noisy
            if reject == 'n':
                colour = "BLUE"
            if reject == 'd':
                colour = "RED"
            annstring = ("COLOR " + colour + "\n" + "CIRCLE " +
            str(stamp.header["radec"][0]) + " " +
            str(stamp.header["radec"][1]) + " 0.005\n")
            ann.write(annstring)

if __name__ == "__main__":
    usage = "genstack.py [options] input_list input_image stamp_size"
    # This OptionParser reads CLI arguments
    parser = OptionParser(usage=usage)
    parser.add_option("-r", action="store_true", dest="raw",
    help="Stack a raw list of images")
    parser.add_option("-l", action="store", type="string", dest="locations",
    help="Use a different set of columns in the input file for the RA, Dec, \
    Peak I (Use format 'RA Dec Peak')")
    parser.add_option("-g", action="store", type="int", dest="scale_factor",
    help="Automatically regrid the image for >1 pixel accuracy")
    parser.add_option("-s", action="store_true", dest="save",
    help="Store each stamp as a FITS file WARNING: Regridding may make the \
    astrometry in the stored stamps incorrect")
    parser.add_option("-a", action="store_true", dest="annotate",
    help="Do not generate an annotation file")
    parser.add_option("-n", action="store_true", dest="noisefile",
    help="Read in from a noise file")
    parser.add_option("-m", action="store", type="float", dest="maxnoise",
    help="Maximum noise in mJy")
    parser.add_option("-p", action="store", type="float", dest="percentile",
    help="Write out a specific percentile along with median")
    parser.add_option("-d", action="store_true", dest="dump_percentiles",
    help="Dump out all each slice of the sorted stack (all percentiles)")
    parser.add_option("-b", action="store_true", dest="subtract_background",
    help="Subtracts the median background value from each stamp prior to \
    stacking")
    parser.add_option("-f", action="store", type="string", dest="fakesource",
    help="Seed a Gaussian at the source location (Use format 'Peak FWHM Dist')")
    parser.add_option("-i", action="store", type="string", dest="imglist",
    help="Stack using the list of images stored in this file")
    # Store the CLI arguments for use in other functions
    (options, args) = parser.parse_args()
    if len(args) < 3 and options.imglist is None:
        print("insufficient arguments.  Use --help for options")
        exit(1)
    elif len(args) < 2:
        print("insufficient arguments.  Use --help for options")
        exit(1)
    if path.exists(args[0].split('/')[-1] + args[1].split('/')[-1] +
    "_mean.fits") or path.exists(args[0].split('/')[-1] +
    args[1].split('/')[-1] + "_median.fits"):
        print("outfiles already exists!")
        exit(1)
    #Default columns for input list
    if options.locations is None:
        ra = 8
        dec = 9
        intense = 14
    else:
        ra = int(options.locations.split()[0])
        dec = int(options.locations.split()[1])
        intense = int(options.locations.split()[2])
    if options.noisefile:
        ra = 0
        dec = 1
        intense = 2
    extract_data(args[0], ra, dec, intense)
    noise = open("noise_" + args[0].split('/')[-1] + ".dat", 'w')
    if not options.annotate:
        annotation = open(args[0].split('/')[-1] + ".ann", 'w')
    #These two variables are used to print the % completion
    percent = 0.0
    fraction = 1.0 / len(RA)
    #args[1:-1] are all of the source images
    if options.imglist is None:
        imglist = args[1:-1]
    else:
        imglist = open(options.imglist).readlines()
    for arg in imglist:
        rejects = []
        try:
            image = fits.open(arg)
            # Load and preserve the wcs data once for each source rather than
            # once per stamp, as the WCS constructor is very slow.
            units = image[0].header["BUNIT"]
            mywcs = wcs.WCS(image[0].header)
            try:
                cdelt1 = image[0].header["CDELT1"]
                cdelt2 = image[0].header["CDELT2"]
            except:
                cdelt1 = mywcs.getXPixelSizeDeg()
                cdelt2 = mywcs.getYPixelSizeDeg()
            cen = 0.5*np.array([image[0].data.shape])+0.5
            centre = mywcs.wcs_pix2world(cen, 0, ra_dec_order=True)[0]
            footprint = mywcs.calc_footprint()
            halfdims = np.abs(footprint[0]-footprint[2])
            if centre[0]+halfdims[0] > 360:
                minmax = [centre[0]-halfdims[0]-360.0,
                centre[0]+halfdims[0]-360.0, centre[1]-halfdims[1],
                centre[1]+halfdims[1]]
            else:
                minmax = [centre[0]-halfdims[0], centre[0]+halfdims[0],
                centre[1]-halfdims[1], centre[1]+halfdims[1]]
            # Iterate through the sources
            for i in range(len(RA)):
                if RA[i] is None:
                    continue
                if RA[i] < minmax[0] or RA[i] > minmax[1] or \
                Dec[i] < minmax[2] or Dec[i] > minmax[3]:
                    rejects.append('d')
                else:
                    # This adds a stamp to the stack, and appends the rejection
                    # flag to the rejects list
                    rejects.append(stack_build(image, 2.0 *
                    round(float(args[-1])/ 2.0), RA[i], Dec[i], Intense[i], mywcs,
                    options.scale_factor, options.save, options.maxnoise,
                    options.subtract_background, options.fakesource))
                    # Generate and print a status message showing % completion
                    if rejects[-1] != 'd':
                        percent = percent + fraction
                        stdout.write("Stack generation @ " + "%(P)06f" %
                        {'P':percent * 100} + "%" + "\r")
                        stdout.flush()
        except IOError:
            print(arg + " failed to load")
            #exit(2)
        # This loop is used to throw out stamps already used,
        # to prevent any duplication across images that overlap
        for rej in range(len(rejects)):
            if rejects[rej] != 'd':
                RA[rej] = None
                Dec[rej] = None
                Intense[rej] = None
    try:
        save_stacks(args[0].split('/')[-1] + args[1].split('/')[-1])
    except AttributeError:
        print("No Sources accepted for stack.  No images will be generated.")
    if options.percentile is not None:
        #Write out the percentile stacks
        save_percentile(args[0].split('/')[-1] + args[1].split('/')[-1],
        options.percentile)
    if not options.annotate:
        annotation.close()
    noise.close()
    stdout.write("\n")
    stdout.flush()
    exit(0)
