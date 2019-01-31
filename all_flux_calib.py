import numpy as np
import star
import os
import pylab as pl

epochs = ['20160115', '20160205', '20160226', '20160318', '20160417']

all_eta = []
all_mean_chi2 = []

for epoch in epochs:

    print epoch

    #create multi-epoch multiplied map, epoch1*epoch2
    star.kappa("mult", epoch+"/"+epoch+"_shifted_mjya2.sdf(-250:-50,100:300)", "../IC348_20151222_850_R1_f200_mjya2.sdf(-250:-50,100:300)", epoch+"/20151222_"+epoch+".sdf")

    #mean and RMS of multi-epoch multiplied map
    star.kappa('stats',epoch+"/20151222_"+epoch+".sdf")
    mean_multi = star.read_starval('stats','mean')[0] #tot: 0.00197281296856839
    sigma_multi = star.read_starval('stats','sigma')[0]

    os.system("rm -f "+epoch+"/20151222_"+epoch+".sdf")


    #create squared map, epoch2*epoch2
    star.kappa("mult", epoch+"/"+epoch+"_shifted_mjya2.sdf(-250:-50,100:300)", epoch+"/"+epoch+"_shifted_mjya2.sdf(-250:-50,100:300)", epoch+"/"+epoch+"_shifted_squared_mjya2.sdf")

    #mean and RMS of squared map
    star.kappa('stats',epoch+"/"+epoch+"_shifted_squared_mjya2.sdf")
    mean_squared = star.read_starval('stats','mean')[0] #tot: 0.165112290522633
    sigma_squared = star.read_starval('stats','sigma')[0]

    os.system("rm -f "+epoch+"/"+epoch+"_shifted_squared_mjya2.sdf")


    #constant to minimize chi^2
    eta = mean_multi / mean_squared
    e_eta = eta * np.sqrt( (sigma_multi/mean_multi)**2 + (sigma_squared/mean_squared)**2 )

    all_eta.append(eta)

    print "eta = " + str(round(eta,5))

    #make chi map
    star.kappa('maths', exp="ia-"+str(eta)+"*ib", ia="../IC348_20151222_850_R1_f200_mjya2.sdf(-250:-50,100:300)", ib=epoch+"/"+epoch+"_shifted_mjya2.sdf(-250:-50,100:300)", out=epoch+"/chi.sdf")

    #make ch^2 map
    star.kappa('mult', epoch+"/chi.sdf", epoch+"/chi.sdf", epoch+"/chi2.sdf")

    os.system("rm -f "+epoch+"/chi.sdf")

    #measure minimized chi^2 map
    star.kappa('stats',epoch+"/chi2.sdf")
    mean_chi2 = star.read_starval('stats','mean')[0]
    sigma_chi2 = star.read_starval('stats','sigma')[0]

    os.system("rm -f "+epoch+"/chi2.sdf")

    all_mean_chi2.append(mean_chi2)

    print "chi^2 = " + str(round(mean_chi2,5))

    #generate flux calibrated map
    star.kappa('cmult', epoch+"/"+epoch+"_shifted_mjya2.sdf", 1-eta, epoch+"/"+epoch+"_shifted_calibrated_mjya2.sdf")

    #generate flux calibrated difference map
    star.kappa('maths', exp="ia-ib", ia="../IC348_20151222_850_R1_f200_mjya2.sdf", ib=epoch+"/"+epoch+"_shifted_calibrated_mjya2.sdf", out=epoch+"/20151222_"+epoch+"_diff.sdf")


    ###RUN TO DOUBLE CHECK ACTUALLY HAVE MINIMIZED CHI^2

    #double check
    all_mean_chi2 = []
    all_sigma_chi2 = []

    #multiplicitive constants to try
    const = eta*np.linspace(-10,10,100)

    #cycle through multiples of eta to see if chi^2 previously got is actually at minimum
    for c in const:

        print c/eta

        #make test chi map
        if c<=0:
            star.kappa('maths',exp="ia"+str(c)+"*ib", ia="../IC348_20151222_850_R1_f200_mjya2.sdf(-250:-50,100:300)", ib=epoch+"/"+epoch+"_shifted_mjya2.sdf(-250:-50,100:300)", out=epoch+"/test_chi.sdf")
        else:
            star.kappa('maths',exp="ia-"+str(c)+"*ib", ia="../IC348_20151222_850_R1_f200_mjya2.sdf(-250:-50,100:300)", ib=epoch+"/"+epoch+"_shifted_mjya2.sdf(-250:-50,100:300)", out=epoch+"/test_chi.sdf")

        #make test chi^2 map
        star.kappa('mult', epoch+"/test_chi.sdf", epoch+"/test_chi.sdf", epoch+"/test_chi2.sdf")

        os.system("rm -f "+epoch+"/test_chi.sdf")

        #measure minimized chi^2 map
        star.kappa('stats',epoch+"/test_chi2.sdf")
        mean_test_chi2 = star.read_starval('stats','mean')[0]
        sigma_test_chi2 = star.read_starval('stats','sigma')[0]

        all_mean_chi2.append(mean_test_chi2)
        all_sigma_chi2.append(sigma_test_chi2)

        os.system("rm -f "+epoch+"/test_chi2.sdf")

        print mean_test_chi2,sigma_test_chi2


    all_mean_chi2 = np.array(all_mean_chi2)
    all_sigma_chi2 = np.array(all_sigma_chi2)


    pl.figure()
    pl.xlabel(r"$\eta$")
    pl.ylabel(r"$\chi^2$")

    pl.plot(const,all_mean_chi2,'k-')
    pl.plot(eta,mean_chi2,'ro')

    pl.axhline(y=mean_chi2,ls='--',c='r')

    pl.show()
    pl.close()


    pl.figure()
    pl.xlabel(r"$\eta$")
    pl.ylabel(r"$\chi^2$")

    pl.plot(const,all_mean_chi2,'k-')
    pl.plot(eta,mean_chi2,'ro')

    pl.axhline(y=mean_chi2,ls='--',c='r')

    pl.ylim(0.0206944,0.0207)
    pl.xlim(-0.03,0.03)

    pl.show()
    pl.close()
    
all_eta = np.array(all_eta)
all_mean_chi2 = np.array(all_mean_chi2)

for i in range(len(epochs)):
    print epochs[i] + '\t' + str(round(all_eta[i],5)) + '\t' + str(round(all_mean_chi2[i],5))
