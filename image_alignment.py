#-
#MIT License
#
#Copyright (c) 2017 Kevin M. Lacaille <kevinlacaille@physics.mcmaster.ca>
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
#

# If this code contributes to a project that leads to a scientific publication, please acknowledge 
# this fact by citing the project. You can acquire the citation from the following journal, and DOI entry: 
# Mairs et al. 2017, ApJ
# DOI: https://www.dropbox.com/s/qaia9w8zxtyuug6/transient_dr_cal.pdf?dl=0

#Description of the code:

#Measure offset between images via a 2D Gaussian fit. Visualize the results.
#This code may be used for unaligned data and for aligned data.

import numpy as np
import pylab as pl
import glob
from scipy.optimize import curve_fit as cf
import time

#define 1D Gaussian
def gauss(x,a,b,c):
	return a*np.exp(-(x-b)**2 / (2*c**2))

#define 2D Gaussian
def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
	xo = float(xo)
	yo = float(yo)    
	a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
	b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
	c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
	g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
	return g.ravel()


#start timing importing data
start1 = time.time()

#import all cross correlated data
fields = glob.glob('directory/to/place/cross/correlation/catalogues/*_crosscorrelation.txt')
fields.sort()
cor_all = []
for i in fields:
	cor_all.append(np.loadtxt(i))
cor_all = np.array(cor_all)

#print time for importing data
end1 = time.time()
print "it took " + str(round(end1-start1)) + " seconds to import all data" + '\n'

#dRA and dDEC from 2D Gaussian fit
dRA = []
dDEC = []
e_dRA = []
e_dDEC = []

#peak cross correlation
peak_cc = []

#how many arcseconds per pixel in your image?
arcsec_pixel = 3.0

#compute the offsets and plot cross correlation
for i in range(len(fields)):
	
	print fields[i]
	print int(np.max(cor_all[i]))

	#append peak cross correlation to array
	peak_cc.append(np.max(cor_all[i]))

	#RA and DEC spaces
	RA_linear = np.arange(-int(len(cor_all[i][0])/2.0), int(len(cor_all[i][0])/2.0)+1)
	DEC_linear = np.arange(-int(len(cor_all[i][0])/2.0), int(len(cor_all[i][0])/2.0)+1)

	#RA-DEC meshgrid
	RA, DEC = np.meshgrid(RA_linear, DEC_linear)


	################
	#2D GAUSSIAN FIT
	################
	
	#limits for fitting (adjust "lower" as needed)
	lower = -5
	upper = -1*lower + 1
	RA_fitting, DEC_fitting = np.meshgrid(RA_linear[np.where(cor_all[i] == np.max(cor_all[i]))[1][0]+lower:np.where(cor_all[i] == np.max(cor_all[i]))[1][0]+upper], DEC_linear[np.where(cor_all[i] == np.max(cor_all[i]))[0][0]+lower:np.where(cor_all[i] == np.max(cor_all[i]))[0][0]+upper])

	#fitting (adjust the initial guess as needed)
	popt, pcov = cf(twoD_Gaussian, (RA_fitting, DEC_fitting), cor_all[i][np.ix_(range(np.where(cor_all[i] == np.max(cor_all[i]))[0][0]+lower,np.where(cor_all[i] == np.max(cor_all[i]))[0][0]+upper),range(np.where(cor_all[i] == np.max(cor_all[i]))[1][0]+lower,np.where(cor_all[i] == np.max(cor_all[i]))[1][0]+upper))].reshape((upper-lower)**2),p0=(np.max(cor_all[i]),0,0,5,5,2,100))

	#fitted data
	data_fitted = twoD_Gaussian((RA, DEC), *popt)

	#append fitted positional offset to arrays
	dRA.append(popt[1])
	dDEC.append(popt[2])
	e_dRA.append(np.sqrt(np.diag(pcov))[1])
	e_dDEC.append(np.sqrt(np.diag(pcov))[2])
	
    #print values of the variables of the 2D Gaussian in terms of arcseconds
    #x offset, x offset uncertainty, y offset, y offset uncertainty
	print round(popt[1]*arcsec_pixel,3), round(np.sqrt(np.diag(pcov))[1],3), round(popt[2]*arcsec_pixel,3), round(np.sqrt(np.diag(pcov))[2],3) + '\n'
	

	##########
	#HEAT MAP
	#########

	#figure
	pl.figure()
	pl.subplots_adjust(hspace=0, wspace=0)
	pl.rc('font', size=13)

	#heat map
	pl.pcolormesh(RA, DEC, cor_all[i])

	#colour bar
	cb = pl.colorbar(cmap='jet')
	cb.set_label('cross correlation')
	pl.clim(vmin=np.min(cor_all[i]),vmax=np.max(cor_all[i]))

	#vertical and horizontal lines going through (0,0)
	pl.axhline(y = 0.5, ls = '--', c = 'k')
	pl.axvline(x = 0.5, ls = '--', c = 'k')

	#contours of 2D Gaussian fit
	pl.contour(RA+0.5, DEC+0.5, data_fitted.reshape(len(RA), len(DEC)), [0.5*data_fitted.max(),0.75*data_fitted.max(),0.85*data_fitted.max(),0.95*data_fitted.max(),0.98*data_fitted.max(),0.99*data_fitted.max()], colors='w')

	#peak of the 2D gaussian fit position
	pl.scatter(popt[1]+0.5,popt[2]+0.5,s=20,facecolor='w',edgecolor='w',alpha=1)

	#limits
	pl.xlim(-10,10)
	pl.ylim(-10,10)

    #labels: units are in terms of pixels
	pl.xlabel(r'RA (pixels)') #this is actually -RA (i.e. the RA is flipped)
	pl.ylabel(r'DEC (pixels)')

    #save figure
    pl.savefig('/directory/to/place/your/figures/'+fields[i][:-4]+'.png',bbox_inches='tight')
	pl.show()
	pl.close()


#convert each positional offset from pixels to arcseconds
#negate the RA to account for backwards map (an artifact from the cross correlation code)
dRA = -arcsec_pixel*np.array(dRA)
dDEC = arcsec_pixel*np.array(dDEC)
e_dRA = arcsec_pixel*np.array(e_dRA)
e_dDEC = arcsec_pixel*np.array(e_dDEC)

#compute the total positional offset and it's uncertainty
offset = np.sqrt(dRA**2 + dDEC**2)
e_offset = 1.0/offset * np.sqrt(dRA**2 * e_dRA**2 + dDEC**2 * e_dDEC**2)

################
# CATALOGUE DATA
################

#write offsets to a file
cat_offset = open('/directory/to/place/offset/catalogue/image_offsets.txt','w')
cat_offset.write('ID\tdRA\te_dRA\tdDEC\te_dDEC\toffset\te_offset\n--\t---\t-----\t----\t------\t------\t--------\n')
for i in range(len(offset)):
    cat_offset.write(fields[i][:-4] + '\t' + str(dRA[i]) + '\t'+ str(e_dRA[i]) + '\t' + str(dDEC[i]) + '\t'+ str(e_dDEC[i]) + '\t' + str(offset[i]) + '\t' + str(e_offset[i]) + '\n')
cat_offset.close()



#########
# FIGURES #change limits and parameters as needed
#########

#bins for histogram (adjust bins as needed)
bins = np.linspace(0,0.16,11)

##########################
# HISTOGRAM OF CORRELATION
##########################
fig = pl.figure()
ax = fig.add_subplot(1,1,1)
pl.subplots_adjust(hspace=0,wspace=0)
pl.rc('font', size=13)
pl.hist(peak_cc,bins=np.logspace(np.log10(min(peak_cc)), np.log10(max(peak_cc)),10))
ax.set_xscale('log')
pl.xlabel('Peak cross correlation')
pl.ylabel('Number')
pl.savefig('/directory/to/place/figures/peak_cc.pdf',bbox_inches='tight')
pl.show()
pl.close()


########################
# CORRELATION VS. OFFSET
########################
pl.figure()
pl.subplots_adjust(hspace=0,wspace=0)
pl.rc('font', size=13)

pl.semilogx(peak_cc, offset, 'bo')

pl.ylabel('Offset (arcsec)')
pl.xlabel('Peak cross correlation')
pl.savefig('/directory/to/place/figures/cc_offset.pdf',bbox_inches='tight')
pl.show()
pl.close()


###########################################
# PERCENT UNCERTAINTY IN OFFSET VS. PEAK CC
###########################################
pl.figure()
pl.subplots_adjust(hspace=0,wspace=0)
pl.rc('font', size=13)
pl.semilogx(peak_cc, e_offset/offset*100, 'bo')
pl.xlabel('Peak cross correlation')
pl.ylabel(r'Offset uncertainty (%)')
pl.savefig('/directory/to/place/figures/offset_uncertainty_percent_cc.pdf',bbox_inches='tight')
pl.show()
pl.close()



##########################################
# PERCENT UNCERTAINTY IN OFFSET VS. OFFSET
##########################################
pl.figure()
pl.subplots_adjust(hspace=0,wspace=0)
pl.rc('font', size=13)
pl.plot(offset, e_offset, 'bo')
pl.xlabel('Offset (arcsec)')
pl.ylabel(r'Offset uncertainty (arcsec)')
pl.savefig('/directory/to/place/figures/offset_uncertainty_offset1.pdf',bbox_inches='tight')
pl.show()
pl.close()

##################################
# UNCERTAINTY IN OFFSET VS. OFFSET
##################################
pl.figure()
pl.subplots_adjust(hspace=0,wspace=0)
pl.rc('font', size=13)
pl.loglog(offset, e_offset/offset*100, 'bo')
pl.plot(np.linspace(1E-3,2E-1,100),3.25*np.linspace(1E-3,2E-1,100)**(-1),'k-',label=r'$\propto$ offset$^{-1}$') #adjust parameters as needed
pl.xlabel('Offset (arcseconds)')
pl.ylabel(r'Offset uncertainty (%)')
pl.xlim(1E-3,2E-1)
pl.legend(loc=1)
pl.savefig('/directory/to/place/figures/offset_uncertainty_offset2.pdf',bbox_inches='tight')
pl.show()
pl.close()

###################
# HISTOGRAM OFFSETS
###################
pl.figure()
pl.subplots_adjust(hspace=0,wspace=0)
pl.rc('font', size=13)
pl.hist(offset,bins=bins)
pl.xlabel('Offset (arcseconds)')
pl.ylabel('Number')
pl.xlim(0,0.16)
pl.savefig('/directory/to/place/figures/offset_distribution.pdf',bbox_inches='tight')
pl.show()
pl.close()

