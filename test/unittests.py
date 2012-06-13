#	PyStack unit tests
#	Ben Keller
#	bwkeller@ucalgary.ca
#	This is a suite of unit testing code for the PyStack suite.  It is mostly
#	low-level coding tests, so most users can ignore this.  It is highly
#	recommended
import unittest
import os
import numpy as np
import pyfits
from astLib import astWCS
from optparse import OptionParser
import genstack
from stampy import Stampy
from stacky import Stacky

class TestStampyStacky(unittest.TestCase):
	def setUp(self):
		self.stack = Stacky()
		self.fits = pyfits.HDUList(pyfits.PrimaryHDU(np.arange(4096).reshape(64,64)))
		self.onesfits = pyfits.HDUList(pyfits.PrimaryHDU(np.ones(4096).reshape(64,64)))
		self.nanfits = pyfits.HDUList(pyfits.PrimaryHDU(np.ones(4096).reshape(64,64)))
		self.nanfits[0].data[0:63][0:63] = np.nan
		for fits in [self.fits, self.onesfits, self.nanfits]:
			fits[0].header.update('BUNIT', 'TestUnits', '')
			fits[0].header.update('CRPIX1', 1, '')
			fits[0].header.update('CRPIX2', 1, '')
			fits[0].header.update('CDELT1', -20.0/3600.0, '')
			fits[0].header.update('CDELT2', 20.0/3600.0, '')
			fits[0].header.update('CRVAL1', -(32*self.fits[0].header['CDELT1']), '')
			fits[0].header.update('CRVAL2', -(32*self.fits[0].header['CDELT2']), '')
			fits[0].header.update('CTYPE1', 'RA-CAR', '')
			fits[0].header.update('CTYPE2', 'DEC-CAR', '')
	def testHeaderValues(self):
		location = (0,0)
		stamp = Stampy(12, self.onesfits, location, name="testname")
		self.assertEqual(stamp.header["units"], "TestUnits")
		self.assertEqual(stamp.header["source"], "testname")
		self.assertEqual(stamp.header["radec"], location)
		self.assertEqual(stamp.header["bmean"], 1)
		self.assertEqual(stamp.header["bmedian"], 1)
	def testFlag(self):
		location = (0, 0)
		badlocation = (64,64)
		stamp = Stampy(12, self.onesfits, location, maxnoise=2)
		self.assert_(not stamp.header["flag"]["highnoise"])
		self.assert_(not stamp.header["flag"]["hasnans"])
		self.assert_(not stamp.header["flag"]["dimensionerror"])
		self.assert_(not stamp.header["flag"]["nearsource"])
		self.assertEqual(stamp.header["flag"]["errorflag"], 'a')
		stamp = Stampy(12, self.fits, location, maxnoise=0)
		self.assert_(stamp.header["flag"]["highnoise"])
		self.assert_(not stamp.header["flag"]["hasnans"])
		self.assert_(not stamp.header["flag"]["dimensionerror"])
		self.assert_(not stamp.header["flag"]["nearsource"])
		self.assertEqual(stamp.header["flag"]["errorflag"], 'n')
		stamp = Stampy(66, self.fits, location)
		self.assert_(not stamp.header["flag"]["highnoise"])
		self.assert_(not stamp.header["flag"]["hasnans"])
		self.assert_(stamp.header["flag"]["dimensionerror"])
		self.assert_(not stamp.header["flag"]["nearsource"])
		self.assertEqual(stamp.header["flag"]["errorflag"], 'd')
		stamp = Stampy(12, self.fits, badlocation)
		self.assert_(not stamp.header["flag"]["highnoise"])
		self.assert_(not stamp.header["flag"]["hasnans"])
		self.assert_(stamp.header["flag"]["dimensionerror"])
		self.assert_(not stamp.header["flag"]["nearsource"])
		self.assertEqual(stamp.header["flag"]["errorflag"], 'd')
		stamp = Stampy(12, self.nanfits, location)
		self.assert_(not stamp.header["flag"]["highnoise"])
		self.assert_(stamp.header["flag"]["hasnans"])
		self.assert_(not stamp.header["flag"]["dimensionerror"])
		self.assert_(not stamp.header["flag"]["nearsource"])
		self.assertEqual(stamp.header["flag"]["errorflag"], 'b')

	def testStackReject(self):
		location = (0,0)
		nanstamp = Stampy(12, self.nanfits, location)
		goodstamp = Stampy(12, self.onesfits, location, maxnoise=2)
		noisestamp = Stampy(12, self.onesfits, location, maxnoise=0)
		dimenstamp = Stampy(66, self.onesfits, location, maxnoise=0)
		self.assertEqual(self.stack.add_stamp(goodstamp), 'a')
		self.assertEqual(self.stack.add_stamp(noisestamp), 'n')
		self.assertEqual(self.stack.add_stamp(dimenstamp), 'd')
		self.assertEqual(self.stack.add_stamp(nanstamp), 'b')

class Testgenstack(unittest.TestCase):
	def setUp(self):
		self.fits = pyfits.HDUList(pyfits.PrimaryHDU(np.ones(4096).reshape(64,64)))
		for fits in [self.fits]:
			fits[0].header.update('BUNIT', 'TestUnits', '')
			fits[0].header.update('CRPIX1', 1, '')
			fits[0].header.update('CRPIX2', 1, '')
			fits[0].header.update('CDELT1', -20.0/3600.0, '')
			fits[0].header.update('CDELT2', 20.0/3600.0, '')
			fits[0].header.update('CRVAL1', -(32*self.fits[0].header['CDELT1']), '')
			fits[0].header.update('CRVAL2', -(32*self.fits[0].header['CDELT2']), '')
			fits[0].header.update('CTYPE1', 'RA-CAR', '')
			fits[0].header.update('CTYPE2', 'DEC-CAR', '')

	def savenoisepos(self, file, ann, stamp, intensity, reject):
		self.assertEqual(stamp.header['source'], "J000012+000300")
		return 0
	def savenoiseneg(self, file, ann, stamp, intensity, reject):
		self.assertEqual(stamp.header['source'], "J000009-005700")
		return 0
	def testNameBuild(self):
		genstack.save_noise = self.savenoisepos
		genstack.stack_build(self.fits, 12, 0.05, 0.05, 1, astWCS.WCS(self.fits[0].header, mode="pyfits"), None, False, 1)
		genstack.save_noise = self.savenoiseneg
		genstack.stack_build(self.fits, 12, 0.04, -0.05, 1, astWCS.WCS(self.fits[0].header, mode="pyfits"), None, False, 1)
	def testSaveStacks(self):
		genstack.save_noise = self.savenoisepos
		genstack.stack_build(self.fits, 12, 0.05, 0.05, 1, astWCS.WCS(self.fits[0].header, mode="pyfits"), None, False, 1)
		genstack.save_noise = self.savenoiseneg
		genstack.stack_build(self.fits, 12, 0.04, -0.05, 1, astWCS.WCS(self.fits[0].header, mode="pyfits"), None, False, 1)
		genstack.units = "TestUnits1"
		genstack.save_stacks("test")
		meanfile = pyfits.open("test_mean.fits")
		medianfile = pyfits.open("test_median.fits")
		self.assertEqual(meanfile[0].data[6,6], 1)
		self.assertEqual(meanfile[0].header["BUNIT"], "TestUnits1")
		self.assertEqual(medianfile[0].header["BUNIT"], "TestUnits1")
		self.assertEqual(meanfile[0].header["CDELT1"], -20.0/3600.0)
		self.assertEqual(medianfile[0].header["CDELT2"], 20.0/3600.0)
		self.assertEqual(meanfile[0].header["CRVAL1"], 120.0/3600.0)
		self.assertEqual(medianfile[0].header["CRVAL2"], -120.0/3600.0) 
		os.remove("test_mean.fits")
		os.remove("test_median.fits")

		
if __name__ == "__main__":
	unittest.main()
