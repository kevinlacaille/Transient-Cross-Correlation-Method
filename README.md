# Transient-Cross-Correlation-Method
These are the Python scripts I have written to compute the cross correlation between two astronomical images. 
The paper can be found here: https://www.dropbox.com/s/76906qbevx3honj/transient_dr_cal.pdf?dl=0

The procedure goes as follows:

1) Compute the cross correlation between two images (see Figure 13).

2) Measure the offset between the peak of the cross correlation product and (RA,DEC) = (0,0) (equation C2).
	- This “radial offset” is simply a measure of how specially offset one image is from another.
	- Note: The “peak of the cross correlation” is actually an estimated peak of the cross correlation product (see Figure 13). What is done here is the script fits a 2D Gaussian to that peak and then used the position of the peak of the Gaussian as the “peak of the cross correlation”.

3) Correct the image offset by shifting one of the images by this measured radial offset. 
  - This is done in an external data reduction program (e.g. Starlink ORACDR, Starlink makemap)
