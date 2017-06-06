#MIT License

#Copyright (c) 2017 Kevin M. Lacaille

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.


#Description of the code:

#Cross correlate a set of images to a single image, then save cross correlation product to a file
#This code may be used for unaligned data and for aligned data.

#NOTICE:
#This code has not been optimised nor has it been tested on consumer laptops for large images.


from astropy.io import fits
from scipy.signal import correlate2d
import glob
import time as time  

#make a list of all .fits files
images = glob.glob('/directory/of/your/images/*.fits')

#sort the images by date
#this step is not necessary
#(note: in order for this to work the date must be in format: yearmonthday e.g. 20151202)
images.sort()

#the image to compare all other images to
first_image = fits.open('/directory/to/your/first/image/image.fits')[0].data

#run through all images to cross correlate with first_image
for i in images:
    
	#image to correlate
	f = fits.open(i)
	epoch = f[0].data
	print i
	
	#start timing cross correlation
	start = time.time()

	#compute the cross correlation between epoch and first_image
	cor = correlate2d(epoch, first_image, mode='full', boundary='fill',fillvalue=0)

    #stop timing cross correlation and print duration of computation
	end = time.time()
	print 'cross correlation took: ' + str(int((end - start)/60.0)) + ' minutes and ' + str(int((end-start) - int((end - start)/60.0)*60))+ ' seconds' + '\n'

	#write cross correlation to a file
	data_file = open('/directory/to/place/cross/correlation/catalogues/'+i[:-5]+'_crosscorrelation.txt','w')
	for row in cor:
		for index in row:
			data_file.write(str(index) + '\t')
		data_file.write('\n')	
	data_file.close()

