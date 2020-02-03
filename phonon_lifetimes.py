#!/usr/bin/env python
from __future__ import print_function
import os
from subprocess import check_output,call
from datetime import datetime as datetime   # datetime.datetime.now()

# lammps_pos_to_sum.py bcc 6 3.253 4 4 4 50   # dominiques job; 6sc ; timestep 1fs; jeden 50 schritt ausgegeben; 3.253 alat;
# nun versucht mit ~/Thermodynamics/python_thermodynamics/phonon_lifetimes.py -s bcc -n 6 -dt 50 -a 3.253 -fftpy
#

# /Users/glensk/Dropbox/proj/proj_current/__2017.01_phonon_linewidth_al/__2016.09_phonon_lifetimes_3_nach_elternzeit/tutild_2016_12_02_10sc_300K_4.04_2_hoch_5_steps

# 2*10^5 schritte --> ein qpunkt hat 19MB
# 1*10^6 schritte --> ein qpunkt hat 200MB  --> ein ast geht locker (20 qpunkte sind 4GB)


###########################################
# Timing of this skript a single qpoint for maximal length
###########################################
# started using:
# python -m cProfile -s cumtime ~/Thermodynamics/python_thermodynamics/lammps_pos_to_sum.py -N 10 -a 4.14 -dt 40 -q 8 8 20 -ps -mdstepstocalc_last
###########################################

##########################################################################################
# run argparse first to make it quick
##########################################################################################
import argparse
import textwrap
#import future
#from argparse import ArgumentDefaultsHelpFormatter
class MyDescriptiveError(Exception):
    pass


def create_READMEtxt(directory=False,add=False):
    ''' wiretes a README.txt file '''
    if directory == False:
        directory = os.getcwd()
    # get sha
    pwd = os.getcwd()
    os.chdir(os.environ['scripts'])
    sha = check_output(["git","rev-parse","master"]).decode('utf-8')
    os.chdir(pwd)

    # get time
    time_now = datetime.now()

    # name of RADME
    filepath = directory+'/README_'+time_now.strftime("%Y-%m-%d_%H:%M:%S")+'.txt'

    # write README.txt
    strout=os.path.basename(sys.argv[0])+" "+" ".join(sys.argv[1:])
    with open(filepath, "w") as text_file:
        text_file.write("# using https://github.com/glensk/dotfiles/trunk/scripts\n")
        text_file.write("# to download it: svn checkout https://github.com/glensk/dotfiles/trunk/scripts\n")
        text_file.write("# used sha: "+sha) #+"\n")
        text_file.write("# execution time: "+str(time.time() - start_time)+" seconds.")
        text_file.write("\n")
        if add:
            print('add')
            print(type(add))
            if type(add) == str:
                text_file.write(add+"\n")
            elif type(add) == list:
                for i in add:
                    text_file.write(i+"\n")
        text_file.write("\n")
        text_file.write(strout+"\n")

    print()
    print('written ',filepath)
    print()
    return filepath

def my_exit():
    raise MyDescriptiveError()

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def lorenz(w,w0,a,b):
    ''' w: frequenzabhaengigkeit
        w0: die frequenz
        a hoehe
        b breite'''
    return a/((w-w0)**2.+b**2.)

def lorenzbesser(w,w0,a,b):
    ''' w: frequenzabhaengigkeit
        w0: die frequenz
        a hoehe
        b breite'''
    return a/((w**2-w0**2)**2.+w**2*b**2.)

def lorenzbesserwiki(w,w0,a,G):
    return lorenz_1p_sum_1p(w,w0,a,G)

def lorenz_1p_sum_1p(w,w0,a,G):
    ''' w: frequenzabhaengigkeit
        w0: die frequenz
        a hoehe  (also called the plasma frequency)
        G breite (gamma)
        took this'''
    return (G*a**2)/((w**2-w0**2)**2.+w**2*G**2.) # finite at w=0

def lorenz_1p_sum_2p(w,w0,a,G):
    ''' w: frequenzabhaengigkeit
        w0: die frequenz
        a hoehe  (also called the plasma frequency)
        G breite (gamma)
        took this'''
    lenw = len(w)-1.
    return (G*a**2)/((w**2-w0**2)**2.+w**2*G**2.) + (G*a**2)/(((-w+lenw)**2-w0**2)**2.+(-w+lenw)**2*G**2.)

def get_Gmax(Gmax,x_to_THz):
    if Gmax * x_to_THz > 4.:
        Gmax = 4.0/x_to_THz
    return Gmax


def my_print_params(params,x_to_THz,rt=4):
    for i in params:
        if params[i].vary==False:
            outtf=printred("False")
        else:
            outtf=printgreen("True")
        print(('----',i,str(np.round(params[i].value*x_to_THz,rt)).ljust(6),"(THz)",str(np.round(params[i].min*x_to_THz,rt)).ljust(6) ,str(np.round(params[i].max*x_to_THz,rt)).ljust(6),'vary:',outtf))

##############################################
# lorenz 68
##############################################
def lorenz_68_Hauer2015_1p_0p(w,w1,a1,G1):
    ''' w: frequenzabhaengigkeit
        w1: die frequenz
        G1 breite (gamma)
        a1 hoehe  (also called the plasma frequency); a=(2kT/m) ;
        Smax(70) = a/(G*w1**2)
        ob man hier a oder a**2 schreibt macht keinen unerschied aber
        der fit wird schneller
        whatever is in the numerator (zaehler) does not matter and only scales a
        '''
    #return (a**2)/((w**2-w1**2)**2.+4*w**2*G**2.) # this only halfes the G value (widths)
    #return (G*a)/((w**2-w1**2)**2.+w**2*G**2.)
    #return (G*a**2)/((w**2-w1**2)**2.+w**2*G**2.) # the a**2 instead of a actually makes the fit faster!  # a**2 works without bounds but does not work well with bounds
    #return (G*a**4)/((w**2-w1**2)**2.+w**2*G**2.) # the a**2 instead of a actually makes the fit faster! a**4 works very well with bounds
    return (w1**2.*G1*a1**4)/((w**2-w1**2)**2.+w**2*G1**2.) # this works definitively best with and without constraints

def lorenz_68_Hauer2015_1p_1p(w,w1,a1,G1):
    ''' w: frequenzabhaengigkeit
        w0: die frequenz
        a hoehe  (also called the plasma frequency)
        G breite (gamma)
        took this'''
    lenw = len(w)-1.
    return lorenz_68_Hauer2015_1p_0p(w,w1,a1,G1) + lorenz_68_Hauer2015_1p_0p(-w+lenw,w1,a1,G1)

def lorenz_68_Hauer2015_2p_0p(w,w1,a1,G1,w2,a2,G2):
    return lorenz_68_Hauer2015_1p_0p(w,w1,a1,G1) + lorenz_68_Hauer2015_1p_0p(w,w2,a2,G2)

def lorenz_68_Hauer2015_2p_2p(w,w1,a1,G1,w2,a2,G2):
    lenw = len(w)-1.
    return lorenz_68_Hauer2015_2p_0p(w,w1,a1,G1,w2,a2,G2) + lorenz_68_Hauer2015_2p_0p(-w+lenw,w1,a1,G1,w2,a2,G2)

def lorenz_68_Hauer2015_3p_0p(w,w1,a1,G1,w2,a2,G2,w3,a3,G3):
    return lorenz_68_Hauer2015_1p_0p(w,w1,a1,G1) + lorenz_68_Hauer2015_1p_0p(w,w2,a2,G2) + lorenz_68_Hauer2015_1p_0p(w,w3,a3,G3)


def lorenz_68_Hauer2015_3p_3p(w,w1,a1,G1,w2,a2,G2,w3,a3,G3):
    lenw = len(w)-1.
    return lorenz_68_Hauer2015_2p_0p(w,w1,a1,G1,w2,a2,G2) + lorenz_68_Hauer2015_2p_0p(-w+lenw,w1,a1,G1,w2,a2,G2) + lorenz_68_Hauer2015_1p_1p(w,w3,a3,G3)

def lorenz_68_Hauer2015_4p_0p(w,w1,a1,G1,w2,a2,G2,w3,a3,G3,w4,a4,G4):
    return lorenz_68_Hauer2015_2p_0p(w,w1,a1,G1,w2,a2,G2) +lorenz_68_Hauer2015_2p_0p(w,w3,a3,G3,w4,a4,G4)


def lorenz_68_Hauer2015_4p_4p(w,w1,a1,G1,w2,a2,G2,w3,a3,G3,w4,a4,G4):
    lenw = len(w)-1.
    return lorenz_68_Hauer2015_2p_2p(w,w1,a1,G1,w2,a2,G2) + lorenz_68_Hauer2015_2p_2p(-w+lenw,w1,a1,G1,w2,a2,G2) + lorenz_68_Hauer2015_2p_2p(w,w3,a3,G3,w4,a4,G4) + lorenz_68_Hauer2015_2p_2p(-w+lenw,w3,a3,G3,w4,a4,G4)

##############################################
# lorenz 86
##############################################
def lorenz_86_Hauer2015_1p_0p(w,w0,a,G):
    ''' w: frequenzabhaengigkeit
        w0: die frequenz
        a hoehe  (also called the plasma frequency)
        G breite (gamma)
        hier zu schreiben a**4 mach die groessenordnungen wie beim lorenz_68_Hauer2015
        '''
    return a**4./((w**2.+G**2.)*((w**2.-4.*w0**2)**2.+4.*w**2.*G**2.)) # original
    #return (a)/((w**2+G**2)*((w**2. - 4.*w0**2)**2.+4.*w**2.*G**2.)) # finite at w=0
    #return (a**2)/(((w**2-4.*w0**2)**2.+4*w**2*G**2.)*(w**2+w0**2)) # finite at w=0
    #return (a**2)/(((w**2-w0**2)**2.+w**2*G**2.)*(w**2+w0**2)) # finite at w=0 # works well with bounds but does not work well with bounds
    #return (w0**2.*G*a**4)/(((w**2-4*w0**2)**2.+4*w**2*G**2.)*(w**2+G**2)) # more or less the original 86 Hauer function ... has a too high peak at w=0
    #return ((w1**4.*G1*a1**8.))/(((w**2-w1**2)**2.+w**2*G1**2.)*(w**2+w1**2)) # this works well with and without constraints/bondaries only when we fit to two both lorenzians left and right; when fitting to only one lorenzian this funcion does not work well with bounds (in frequency and lw max)

def lorenz_86_Hauer2015_1p_1p(w,w0,a,G):
    ''' w: frequenzabhaengigkeit
        w0: die frequenz
        a hoehe  (also called the plasma frequency)
        G breite (gamma)
        took this'''
    lenw = len(w)-1.
    return lorenz_86_Hauer2015_1p_0p(w,w0,a,G) + lorenz_86_Hauer2015_1p_0p(-w+lenw,w0,a,G)

def lorenz_86_Hauer2015_2p_0p(w,w1,a1,G1,w2,a2,G2):
    return lorenz_86_Hauer2015_1p_0p(w,w1,a1,G1) + lorenz_86_Hauer2015_1p_0p(w,w2,a2,G2)

def lorenz_86_Hauer2015_3p_0p(w,w1,a1,G1,w2,a2,G2,w3,a3,G3):
    return lorenz_86_Hauer2015_1p_0p(w,w1,a1,G1) + lorenz_86_Hauer2015_1p_0p(w,w2,a2,G2) + lorenz_86_Hauer2015_1p_0p(w,w3,a3,G3)

def lorenz_86_Hauer2015_4p_0p(w,w1,a1,G1,w2,a2,G2,w3,a3,G3,w4,a4,G4):
    return lorenz_86_Hauer2015_1p_0p(w,w1,a1,G1) + lorenz_86_Hauer2015_1p_0p(w,w2,a2,G2) + lorenz_86_Hauer2015_1p_0p(w,w3,a3,G3) + lorenz_86_Hauer2015_1p_0p(w,w4,a4,G4)


def lorenz_86_Hauer2015_2p_2p(w,w1,a1,G1,w2,a2,G2):
    lenw = len(w)-1.
    return lorenz_86_Hauer2015_2p_0p(w,w1,a1,G1,w2,a2,G2) + lorenz_86_Hauer2015_2p_0p(-w+lenw,w1,a1,G1,w2,a2,G2)

def lorenz_86_Hauer2015_3p_3p(w,w1,a1,G1,w2,a2,G2,w3,a3,G3):
    lenw = len(w)-1.
    return lorenz_86_Hauer2015_3p_0p(w,w1,a1,G1,w2,a2,G2,w3,a3,G3) + lorenz_86_Hauer2015_3p_0p(-w+lenw,w1,a1,G1,w2,a2,G2,w3,a3,G3)

def lorenz_86_Hauer2015_4p_4p(w,w1,a1,G1,w2,a2,G2,w3,a3,G3,w4,a4,G4):
    lenw = len(w)-1.
    return lorenz_86_Hauer2015_4p_0p(w,w1,a1,G1,w2,a2,G2,w3,a3,G3,w4,a4,G4) + lorenz_86_Hauer2015_4p_0p(-w+lenw,w1,a1,G1,w2,a2,G2,w3,a3,G3,w4,a4,G4)


##############################################
# michael
##############################################
def michael_mail_1p_0p(w,w0,a,G):
    '''
    whatever is in the numerator (zaehler) does not matter and only scales a
    compare to lorenz_68_Hauer2015_1p_0p
    4.*w**2*G**2. ist fuer HWHM
    1.*w**2*G**2. ist fuer FWHM
    '''
    #return (w0**2.*G*a**4)/((w**2-w0**2)**2.+4.*w**2*G**2.)  # das ist fuer HWHM
    #return (w0**2.*G*a**4)/((w**2-w0**2)**2.+1.*w**2*G**2.)  # das ist fuer FWHM
    #return (w0**2.*G*a**2)/((w**2-w0**2)**2.+w**2*G**2.)  # das ist fuer FWHM (konvergiert nur mit bounds)
    return (w0**2.*G*a**4)/((w**2-w0**2)**2.+w**2*G**2.)     # das ist fuer FWHM (konvergiert ohne bounds)

def michael_mail_1p_1p(w,w0,a,G):
    lenw = len(w)-1.
    return michael_mail_1p_0p(w,w0,a,G) + michael_mail_1p_0p(-w+lenw,w0,a,G)


##############################################
# lorenz 6886
##############################################
def lorenz_6886_Hauer2015_1p_1p(w,w1,a1,G1):
    return lorenz_68_Hauer2015_1p_1p(w,w1,a1,G1)*lorenz_86_Hauer2015_1p_1p(w,w1,a1,G1)

def lorenz_6886_Hauer2015_1p_0p(w,w1,a1,G1):
    return lorenz_68_Hauer2015_1p_0p(w,w1,a1,G1)*lorenz_86_Hauer2015_1p_0p(w,w1,a1,G1)

def lorenz_6886_Hauer2015_2p_0p(w,w1,a1,G1,w2,a2,G2):
    return lorenz_68_Hauer2015_2p_0p(w,w1,a1,G1,w2,a2,G2)*lorenz_86_Hauer2015_2p_0p(w,w1,a1,G1,w2,a2,G2)

def lorenz_6886_Hauer2015_2p_2p(w,w1,a1,G1,w2,a2,G2):
    return lorenz_68_Hauer2015_2p_2p(w,w1,a1,G1,w2,a2,G2)*lorenz_86_Hauer2015_2p_2p(w,w1,a1,G1,w2,a2,G2)

def lorenz_6886_Hauer2015_3p_3p(w,w1,a1,G1,w2,a2,G2,w3,a3,G3):
    return lorenz_68_Hauer2015_3p_3p(w,w1,a1,G1,w2,a2,G2,w3,a3,G3)*lorenz_86_Hauer2015_2p_0p(w,w1,a1,G1,w2,a2,G2,w3,a3,G3)


def lorenz_6886_Hauer2015_4p_4p(w,w1,a1,G1,w2,a2,G2,w3,a3,G3,w4,a4,G4):
    return lorenz_68_Hauer2015_3p_3p(w,w1,a1,G1,w2,a2,G2,w3,a3,G3,w4,a4,G4)*lorenz_86_Hauer2015_2p_0p(w,w1,a1,G1,w2,a2,G2,w3,a3,G3,w4,a4,G4)




def lorenzbesserwiki2(w,w0,a,G,w02,a2,G2):
    ''' w: frequenzabhaengigkeit
        w0: die frequenz
        a hoehe
        b breite'''
    return a**2/((w**2-w0**2)**2.+w**2*G**2.) + a2**2/((w**2-w02**2)**2.+w**2*G2**2.)



##############################################
# parameters for t2
##############################################
def fix_params_long_and_t1_for_t2(params,space_fft_which,path,qpoint,filename,mdsteps,x_to_THz,estimate_w2_a2_G2=False):
    return params


    if space_fft_which == 'space_fft_phase_' and path == "t2_t2_0": # only here crossing
        qpoint_ll0,qpoint_t1_t1_0 = qpoint_map_t2_qpoint_to_L_L_0_and_T1_T1_0_qpoint(qpoint,args)
        #long_estimate = get_from_db_w_a_G(qpoint_ll0,filename,1,mdsteps)
        #t1_estimate   = get_from_db_w_a_G(qpoint_t1_t1_0,filename,1,mdsteps)
        long_estimate = get_from_db_w_a_G(qpoint_ll0,filename,args.peaks,mdsteps)
        t1_estimate   = get_from_db_w_a_G(qpoint_t1_t1_0,filename,args.peaks,mdsteps)
        if type(t1_estimate) != bool:
            #print "**** from db (LL0 T1):",t1_estimate
            params.add('w1', value=t1_estimate[0]/x_to_THz, vary=False)
            params.add('a1', value=t1_estimate[1]/x_to_THz, vary=False)
            params.add('G1', value=t1_estimate[2]/x_to_THz, vary=False)
        if type(long_estimate) != bool:
            #print "**** from db (LL0 Long):",long_estimate
            params.add('w3', value=long_estimate[0]/x_to_THz, vary=False)
            params.add('a3', value=long_estimate[1]/x_to_THz, vary=False)
            params.add('G3', value=long_estimate[2]/x_to_THz, vary=False)

        if estimate_w2_a2_G2:
            f_t1,f_t2,f_l = estimate_freqs_from_at_L_L_0_path(qpoint=qpoint,args=args)
            params.add('w2', value=f_t2/x_to_THz, min=f_t2/x_to_THz*0.5, max = f_t2/x_to_THz*2,vary=True)
            params.add('a2', value=t1_estimate[1]/x_to_THz, min=0.,    max=np.inf)
            params.add('G2', value=t1_estimate[2]/x_to_THz, min=0.,    max=5/x_to_THz)
    return params


def my_little_minimizer(x,y,pin,function,xind,x_to_THz,sigma,verbose):
    values = 3  # 3 means three values on the left and 3 values on the right
    if sigma == 0.:
        parttry = [0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001,0.00005,0.00002,0.00001,0.000005,0.000002,0.000001]
    if sigma > 0.:
        parttry = [0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001,0.00005,0.00002,0.00001,0.000005,0.000002,1e-6,5e-7,2e-7,1e-7,5e-8,2e-8,1e-8,5e-9,2e-9,1e-9,5e-10,2e-10,1e-10,5e-11,2e-11,1e-11,5e-12,2e-12,1e-12,5e-13,2e-13,1e-13,5e-14,2e-14,1e-14]

    #if verbose:
    #    print 'sigma',sigma,'VALUES IN:        ',std,pin*x_to_THz


    for part in parttry:
        for chidx in np.arange(len(pin)):
            #print 'pin',pin
            pin,std = getpin(x,y,pin,function,x_to_THz,sigma=sigma,chidx=chidx,part=part,values=4,verbose=False)
        if verbose:
            print(('sigma',sigma,'part',str(part).ljust(7),std,pin*x_to_THz))
    return pin




def getpin(x,y,pin,function,x_to_THz,sigma,chidx,part=0.01,values=3,verbose=False):
    def getrange(pin,part=0.01,values=values):
        diff = (pin*part)/(values-1)
        min = pin-diff*values
        max = pin+diff*values
        if min < 0.:
            min = 0
            diff = pin/(values+1.)
            #print printred('part'),part,values
            min = pin-diff*values
            max = pin+diff*values
            #diff = (pin*part)/(values-1)
        return np.linspace(min,max,2*values+1)

    def getpuse(pin,chidx,var):
        out = np.array(np.copy(pin))
        out[chidx] = var
        return out
        #if chidx == 0:
        #    return np.array([var,pin[1],pin[2]])
        #if chidx == 1:
        #    return np.array([pin[0],var,pin[2]])
        #if chidx == 2:
        #    return np.array([pin[0],pin[1],var])

    def getout(pin,chidx,var,x,function,sigma):
        #print 'xxx',pin,chidx,var
        puse = np.abs(getpuse(pin,chidx,var))
        #print 'aaa',puse
        if puse[chidx] != var:
            sys.exit("puse[chidx] != var")
        yfit = function(x, *puse)


        diff = y-yfit
        diffabs = np.abs(diff)
        diffsum  = np.sum(diffabs)
        #diffsum  = np.sum(diffabs**2.) # NO, NOT as good for MD convergence
        #diffsum = leastsq_xx = lmfitfit_xx = np.sum(np.abs((y-yfit)**2.))
        #if False:
        #if diffsum*1000. < 324144.304:
        if False:
            np.savetxt('y',y)
            np.savetxt('yfit',yfit)
            np.savetxt('ydiff',diff)
            np.savetxt('ydiffabs',diffabs)
            print('done')
            sys.exit('33ff')
        #print 'diffsum',diffsum*1000.

        ################################################################################
        # using the following definition ends up in (converges to) the same result that lmfit would get
        # BUT : #SEEMS NOT AS STABLE AS WHEN MINIMIZING |y-y'| without squarring
        ################################################################################
        #leastsq = np.abs((y-yfit)**2.)
        #diffsumls  = np.sum(leastsq) #SEEMS NOT AS STABLE AS WHEN MINIMIZING |y-y'| without squaring
        #diffsum = diffsumls

        #diffsum  = np.sum(np.abs((y-yfit)**2.))  converges to the save value eventuall but needs a much longer md run!
        # for Al about converged at 180.000 steps but before slow convergence
        #print diffsum
        if sigma > 0:
            weights = yfit2 = np.copy(yfit)**sigma  # this is ok
            #diffsum = np.sum((diffabs)*yfit2)
            # here we need to use the standard leastsquares procedure (other-
            # wise the minimizer makes a =0)
            #diffsum = np.sum((diffabs**2.)*yfit2) # NOT STABLE
            #diffsum = np.sum(diffabs*yfit2)  # NOT STABLE even quicker so
            diffsum  = np.sum((diffabs**2.)**2.) # is this for sigma == 1?  LMFIT?
            # das diff ist ok so, aber das yfit2 wird einfach zu null und somit "der beste fit"
            #
        #print var,diffsum,abssum2
        #return var,diffsum,abssum2
        return var,diffsum


    def get_var_and_std_best(o):
        aa = np.where(o[:,0]>0.)
        o = o[aa]
        index = np.where(o[:,1] == o[:,1].min())[0][0]
        #print 'index',index
        #print o.shape
        border = False
        if index == 0:
            border = True
        if index == o.shape[0]-1:
            border = True
        return o[index][0],o[index][1],border

    def get_o_out0w_best(pin,chidx,x,function,sigma,part,range_factor,values):
        redo = False
        range_ = getrange(pin[chidx],part=part*range_factor,values=values)
        #print '--> pin',pin
        #print 'chidx',chidx
        #print 'range_',range_
        #print 'r',range_*x_to_THz,pin*x_to_THz,chidx,part,range_factor,values
        #@@print 'r',pin*x_to_THz,chidx,part,range_factor,values
        o = np.zeros((len(range_),2))
        for id,var in enumerate(range_):
            o[id,0],o[id,1] = getout(pin,chidx,var,x,function,sigma)
        varbest,varstd,varusebiggerrange = get_var_and_std_best(o)

        ##print 'out[0]',out[0][0]
        ##print 'out[-1]',out[-1][0]
        ##print 'var0best',var0best,var0std
        ##print 'varwbest',varwbest,varwstd
        #p0 = np.polyfit(out[:,0], out[:,1], 2)
        #pw = np.polyfit(out[:,0], out[:,2], 2)
        #f0 = np.poly1d(p0)
        #fw = np.poly1d(pw)
        #x0 = -f0.c[1]/(2.*f0.c[0])
        #xw = -fw.c[1]/(2.*fw.c[0])
        #if verbose:
        #    print 'x0',x0*x_to_THz,f0(x0)*x_to_THz
        #    print 'xw',xw*x_to_THz,fw(xw)*x_to_THz
        #x0_within = True
        #xw_within = True
        if varusebiggerrange == True:
            import warnings
            warnings.simplefilter('ignore', np.RankWarning)
            p = np.polyfit(o[:,0], o[:,1], 2)
            f = np.poly1d(p)
            varbest_test = -f.c[1]/(2.*f.c[0])
            if varbest_test > varbest:
                check = varbest_test/varbest
            if varbest > varbest_test:
                check = varbest/varbest_test
            #@@print 'x',varbest_test*x_to_THz,f(x)*x_to_THz,check
            if varbest_test > 0. and check < 1.5: # not more than 50 % change
                varbest_test,varstd_test = getout(pin,chidx,varbest_test,x,function,sigma)
                #print '->varbest     ',varbest*x_to_THz,varstd*x_to_THz
                #print '->varbest_test',varbest_test*x_to_THz,varstd_test*x_to_THz
                if varstd_test < varstd:
                    range_ = getrange(pin[chidx],part=part*range_factor,values=values)
                    #print printgreen('changing it!')
                    redo = True
                    varbest,varstd,varusebiggerrange = varbest_test,varstd_test, False

            #x0_within = True


        #if x0 > out[-1][0] or x0 < out[0][0]:
        #    x0_within = False
        #    print 'yes0',x0*x_to_THz
        #    var,a,b = getout(pin,chidx,x0,x,function,sigma)
        #    print 'yes0,v,a,b',np.array([var,a,b])*x_to_THz
        #    if a < var0std:
        #        var0best = var
        #        var0std = a
        #if xw > out[-1][0] or xw < out[0][0]:
        #    xw_within = False
        #    print 'yes1',xw*x_to_THz
        #    var,a,b = getout(pin,chidx,xw,x,function,sigma)
        #    print 'yes1,v,a,b',np.array([var,a,b])*x_to_THz
        #    if b < varwstd:
        #        varwbest = var
        #        varwstd = b
        #if verbose:
        #    print "x0_within",x0_within
        #    print "xw_within",xw_within

        #return o,var0best, varwbest, var0std, varwstd,x0_within,xw_within
        return varbest,varstd,varusebiggerrange,o,redo,range_


    # start with w minimization, get w0
    # get take w0 and get a0 minimization
    # use w0 and a0 and get get G0
    if chidx == 0: add = 'frequency '
    if chidx == 1: add = 'hight     '
    if chidx == 2: add = 'linewidths'
    if verbose > 1:
        print()
        print(( add,'will be adjusted. chidx:',chidx))
        print(( 'pin',pin*x_to_THz))
        print(( 'part',part))
        print(( 'values',values))



    varusebiggerrange = False
    for range_factor in [1.,10.,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11]:
        #print 'range',range_factor
        #print
        #print 'pin',pin*x_to_THz
        #range = getrange(pin[chidx],part=part*range_factor,values=values)
        varbest, varstd, varusebiggerrange,o,redo,range = get_o_out0w_best(pin,chidx,x,function,sigma,part,range_factor,values)
        pin = getpuse(pin,chidx,varbest)
        if redo:
            #print printred('rodo')
            varbest, varstd, varusebiggerrange,o,redo,range = get_o_out0w_best(pin,chidx,x,function,sigma,part,range_factor,values)
            pin = getpuse(pin,chidx,varbest)


        add2 = printgreen(str(varusebiggerrange))
        if varusebiggerrange == True:
            add2 = printred(str(varusebiggerrange))
        if verbose:
            print((add,str(part).ljust(6),'range',str(np.round(range*x_to_THz,4)).ljust(73),'pin',np.round(pin*x_to_THz,4),np.round(varbest*x_to_THz,5),np.round(varstd,5),add2))
        if varusebiggerrange == False:
            break





    if varusebiggerrange == True:
        print(( 'o',range_factor))
        print( o)
        print(( 'len(x)',len(x),'==mdsteps'))
        sys.exit('range seems too small yet')


    pin = getpuse(pin,chidx,varbest)
    return pin,varstd



def fit_curve_fit(function,x,y,p0,bounds,sigma):
    #print "type(bounds in)",type(bounds)
    if type(bounds) == bool:  # bounds == False
        popt0, pcov0 = curve_fit(function, x, y,p0=p0,sigma=sigma) #,maxfev=20000)
    else: # bounds is set to a touple
        popt0, pcov0 = curve_fit(function, x, y,p0=p0,bounds=bounds,sigma=sigma) #,maxfev=20000)
    return popt0, pcov0

def check_if_parameters_converged(pin,pout,d):
    ''' text '''
    check0 = np.abs(1.- pin/pout)

    #for idx,i in enumerate(check0):
    #    print 'check:',i,pin[idx],pout[idx]
    #print 'done '
    check = check0 < d
    verbose = False
    #if verbose:
    #    print "<<",i,popt0*x_to_THz,"   ",check0

    #getout = False
    ##print "lll",len(pin)
    #if len(pin) == 3:

    #    if check[0] == check[1] == check[2] == True:
    #        return True, check0
    #elif len(pin) == 6:
    #    if check[0] == check[1] == check[2] == check[3] == check[4] == check[5] == True:
    #        return True, check0
    #elif len(pin) == 12:
    #    if check[0] == check[1] == check[2] == check[3] == check[4] == check[5] == True:
    #        return True, check0
    #else:
    #    sys.exit('either len 3 or 6 or 12')
    return np.alltrue(check),check0

def fit_several_peaks_by_holding_fix_all_but_one_peak(pin_inner,function,x,y,pin=False,kwfit=False,esto=False,mdsteps=False,verycertainestimate=False,space_fft_which=False,path=False):
    x_to_THz = kwfit[0]
    xind = kwfit[1]
    psmax = kwfit[2]
    verbose = kwfit[3]
    sigma = kwfit[4]
    ############################################################################################
    # be conservative with changes in general, we assume we have a sound GROB fit (good pin)
    # pin = [  883.700137   265.179858    69.21936   1021.962462   343.197157    69.482065] fuer 2 peaks
    ############################################################################################
    if int(len(pin)/3) == 2:  # 2 peaks
        variation = 0.1  # 0.06 was too little
        #variation = 0.33
        if verycertainestimate == True:
            variation = 0.06
    if int(len(pin)/3) == 3:
        variation = 0.05
    if int(len(pin)/3) == 3:
        variation = 0.04
    if int(len(pin)/3) == 4:
        variation = 0.03
    mwl = 1. - variation
    mwr = 1. + variation
    mal = 1. - variation*3.
    mar = 1. + variation*3.
    mGl = 1. - variation*12.
    mGr = 1. + variation*12.

    verbose_inner = False

    modeli = Model(function)
    params = modeli.make_params()
    for peak_idx_vary in np.arange(int(len(pin)/3)):
        if verbose > 2:
            print()
            print((' XX ###############',pin_inner,' peak_idx_vary',peak_idx_vary,"########################## out of",np.arange(int(len(pin)/3))))
            print((' XX verycertainestimate:',verycertainestimate))
        piv = peak_idx_vary
        #print
        #print 'peak_idx_vary',peak_idx_vary
        if verbose_inner:
            print(('pin start:',pin))
        Gval = pin[piv*3+2]
        if verycertainestimate == True:
            Gval = esto[piv*3+2]

        if Gval*x_to_THz > 4.:
            Gval = 2./x_to_THz

        Gmin = pin[piv*3+2]*mGl
        if verycertainestimate == True:
            Gmin = esto[piv*3+2]*mGl
        if Gmin * x_to_THz > 4.:
            Gmin = 0.1/x_to_THz
        if Gmin < 0:
            Gmin = 0

        Gmax = pin[piv*3+2]*mGr
        if verycertainestimate == True:
            Gmax = esto[piv*3+2]*mGr
        if Gmax * x_to_THz > 4.:
            Gmax = 4.0/x_to_THz

        aval = pin[piv*3+1]
        amin = pin[piv*3+1]*mal
        amax = pin[piv*3+1]*mar
        if aval == 0:
            aval = 1
            amin = 0
            amax = np.inf

        if verbose_inner:
            print(( 'a'+str(piv+1),"val:",pin[piv*3+1],mal,mar))
            print(( 'a'+str(piv+1),"min:",pin[piv*3+1]*mal))
            print(( 'a'+str(piv+1),'max:',pin[piv*3+1]*mar))

        wval = pin[piv*3+0]
        wmin = pin[piv*3+0]*mwl
        wmax = pin[piv*3+0]*mwr
        if type(esto) != bool:
            w_unc = 0.15 # 15 %
            if mdsteps >  1000: w_unc = 0.12 # 12 %
            if mdsteps >  3000: w_unc = 0.07 # 12 %
            if mdsteps >  5000: w_unc = 0.04
            if mdsteps >  7000: w_unc = 0.03
            if mdsteps > 10000: w_unc = 0.02
            #print 'esto',esto
            wmin = (esto[piv*3+0]*(1.-w_unc))
            wmax = (esto[piv*3+0]*(1.+w_unc))
            if wval > wmax or wval < wmin:
                wval = esto[piv*3+0]
        if verbose > 2:
            print((' XX wval:'+str(piv+1)+':',wval*x_to_THz,'wmin:',wmin*x_to_THz,'wmax:',wmax*x_to_THz))
            print((' XX a   :',aval*x_to_THz,'amin',amin*x_to_THz,'amax',amax*x_to_THz))
            print((' XX G   :',Gval*x_to_THz,'Gmin',Gmin*x_to_THz,'Gmax',Gmax*x_to_THz))
            print()
            print((' XX wval:'+str(piv+1)+':',wval,'wmin:',wmin,'wmax:',wmax))
            print((' XX a   :',np.round(aval,2),'amin',np.round(amin,2),'amax',np.round(amax,2)))
            print((' XX G   :',np.round(Gval,2),'Gmin',np.round(Gmin,2),'Gmax',np.round(Gmax,2)))
        params.add('w'+str(piv+1), value=wval , min=wmin  , max=wmax)
        params.add('a'+str(piv+1), value=aval , min=amin  , max=amax)
        params.add('G'+str(piv+1), value=Gval , min=Gmin  , max=Gmax)

        psigma = np.array([wval,pin[piv*3+1],Gval])
        if psigma[2]*x_to_THz < 0.3: psigma[2] = 0.5/x_to_THz
        if verbose > 2:
            print(' XX psigma (for weights)',psigma*x_to_THz)
        weights = lorenz_68_Hauer2015_1p_1p(x, *psigma)**sigma  # her an even fctn might be better


        for peak_idx_fix in np.arange(int(len(pin)/3)):
            if peak_idx_vary == peak_idx_fix:
                continue
            pif = peak_idx_fix
            #print 'peak_idx_vary',peak_idx_vary,'peak_idx_fix',peak_idx_fix
            Gval = pin[peak_idx_fix*3+2]
            if Gval*x_to_THz > 4.:
                Gval = 4./x_to_THz

            params.add('w'+str(peak_idx_fix+1), value=pin[peak_idx_fix*3+0], vary=False)
            params.add('a'+str(peak_idx_fix+1), value=pin[peak_idx_fix*3+1], vary=False)
            params.add('G'+str(peak_idx_fix+1), value=Gval                 , vary=False)

            #print 'a'+str(pif+1),"val:",pin[peak_idx_fix*3+1]

        if verbose_inner:
            print('iter',peak_idx_vary)
        if verbose > 2:
            print()
            print(' XX Xwval:'+str(piv+1)+':',wval,'wmin:',wmin,'wmax:',wmax)
            print(' XX Xa   :',np.round(aval,2),'amin',np.round(amin,2),'amax',np.round(amax,2))
            print(' XX XG   :',np.round(Gval,2),'Gmin',np.round(Gmin,2),'Gmax',np.round(Gmax,2))
            if False:
                if peak_idx_vary == 1:
                    #params.add('w'+str(peak_idx_fix+1-1), value=pin[(peak_idx_fix-1)*3+0]*1.02, vary=False)
                    #params.add('w'+str(peak_idx_fix+1), value=pin[peak_idx_fix*3+0]*1.03, vary=False)
                    params.add('w'+str(peak_idx_fix+1), value=1168, vary=False)
                    params.add('G'+str(peak_idx_fix+1), value=553, vary=False)

                    params.add('w'+str(peak_idx_fix+1), value=1230, vary=False)
                    #params.add('a'+str(peak_idx_fix+1), value=3.7629, vary=False)
                    params.add('G'+str(peak_idx_fix+1), value=510, vary=False)

        params = fix_params_long_and_t1_for_t2(params,space_fft_which,path,qpoint,filename,mdsteps,x_to_THz)
        if verbose > 2:
            print('############## zuvor ######################## >>>>>>>>')
            my_print_params(params,x_to_THz,rt=4)
            print('############## zuvor ######################## <<<<<<<<')
        result = modeli.fit(y, params, w=x, weights=weights,method='nelder')
        #if verbose > 2:
        #    my_print_params(params,x_to_THz,rt=4)
        pout = lmfit_result_to_popt(result,xind)
        if verbose > 2:
            #print '############## danach ######################## >>>>>>>>'
            #print(result.fit_report())
            #print '############## danach ######################## <<<<<<<<'
            outorig_y = function(x, *pout)
            xout,yout = ps_to_ps_for_writing(outorig_y,xind,x_to_THz,function,x,pout,args=args)
            parnpsavetxt("out_"+str(peak_idx_vary)+"_"+str(wval),np.transpose(np.array([x*x_to_THz,outorig_y*psmax])),fmt='%.8f %.15f')



        if verbose > 2:
            print(printred(' XX =====>>> pin -0'),np.round(pin*x_to_THz,3))
            print(printred(' XX =====>>> pout-0'),np.round(pout*x_to_THz,3))
            print(printred(' XX =====>>> DIFF  '),np.round((pin-pout)*x_to_THz,3))
        pin = np.copy(pout)
        if verbose_inner:
            print('iter done',peak_idx_vary,"--->",pout)
        #print(result.fit_report())
    if verbose_inner:
        print('whole thing done')

    if verbose > 2:
        print('pout-1->',np.round(pout*x_to_THz,3))
    #! # adapt w's keep everything else fixed THIS DESTROYS 2 2 4 from A
    #! /Users/glensk/Dropbox/proj/proj_current/__2017.03_irina_anharmonic/17_08_03_CrN/DLM_300K_2
    #! psigma = np.copy(pout)
    #! if len(xind) > 0:
    #!     params.add('w1', value=pout[0], min=pout[0]*0.95,max=pout[0]*1.05)
    #!     params.add('a1', value=pout[1], vary=False)
    #!     params.add('G1', value=pout[2], vary=False)
    #!     if psigma[2]*x_to_THz < 0.3: psigma[2] = 0.5/x_to_THz
    #! if len(xind) > 1:
    #!     params.add('w2', value=pout[3], min=pout[3]*0.95,max=pout[3]*1.05)
    #!     params.add('a2', value=pout[4], vary=False)
    #!     params.add('G2', value=pout[5], vary=False)
    #!     if psigma[5]*x_to_THz < 0.3: psigma[2] = 0.5/x_to_THz
    #! if len(xind) > 2:
    #!     params.add('w3', value=pout[6], min=pout[6]*0.95,max=pout[6]*1.05)
    #!     params.add('a3', value=pout[7], vary=False)
    #!     params.add('G3', value=pout[8], vary=False)
    #!     if psigma[8]*x_to_THz < 0.3: psigma[2] = 0.5/x_to_THz
    #! if len(xind) > 3:
    #!     params.add('w4', value=pout[9], min=pout[9]*0.95,max=pout[9]*1.05)
    #!     params.add('a4', value=pout[10],vary=False)
    #!     params.add('G4', value=pout[11],vary=False)
    #!     if psigma[11]*x_to_THz < 0.3: psigma[2] = 0.5/x_to_THz
    #! weights = function(x, *psigma)**sigma
    #! result = modeli.fit(y, params, w=x, weights=weights)
    #! pout = lmfit_result_to_popt(result,xind)

    #! if verbose > 2:
    #!     print 'pout-2->',np.round(pout*x_to_THz,3)
    #! # adapt w's keep everything else fixed
    return pout

def fit_2_peaks_fix_one_peak(model,x,y,pin,fixleft=False,fixright=False,plotit=False,kwfit=False,fixleft_values=False):
    if type(fixleft) != bool:
        sys.exit('fixleft needs to be True or False')
    if type(fixright) != bool:
        sys.exit('fixright needs to be True or False')
    if fixleft == fixright:
        sys.exit('fixleft == fixright')
    x_to_THz = kwfit[0]
    xind = kwfit[1]
    psmax = kwfit[2]
    verbose = kwfit[3]

    ###########################################################
    # be conservative with changes in general
    ###########################################################
    mwl = 0.9
    mwr = 1.1
    mal = 0.9
    mar = 1.1
    mGl = 0.9
    mGr = 1.1
    ratio_act = False

    ###########################################################
    # check initial fit for ratio
    ###########################################################
    if False:
        outorig_y = model(x, *pin)
        xout,yout = ps_to_ps_for_writing(outorig_y,xind,x_to_THz,model,x,pin,args=args)
        ratio = y.max()/yout.max()
        if 0.2 < ratio < 3.0:
            mal = 0.666
            mar = 1.5
        if 0.6666 < ratio < 1.5:  # be conservative if close to the truth
            mal = 0.9
            mar = 1.1
        if ratio <= 0.2 or ratio >= 3.:
            ratio_act = True
            #print '222: (y*psmax).max()',(y*psmax).max()
            #print '222: outorig_y.max()',(yout*psmax).max()
            #parnpsavetxt("INIT_B_ACT",np.transpose(np.array([xout,yout*psmax])),fmt='%.8f %.5f')
            #sys.exit()

        #print printblue('222: ratio: '+str(ratio))
        # welcher peak muss angepasst werden?
        #print '222 xind',xind*x_to_THz
        xind_peak_THz_curr = xout[np.where(yout==yout.max())[0][0]]
        #print '222 xind_peak_THz_curr',xind_peak_THz_curr
        xind_peak_THz = (find_nearest(xind*x_to_THz, xind_peak_THz_curr))
        #print '222 xind_peak_THz',xind_peak_THz
        xind_peak_idx = np.where(xind*x_to_THz==xind_peak_THz)[0][0]
        #print 'xind_peak_idx',xind_peak_idx





    #print 'popt in    :',np.around(pin*x_to_THz,decimals=2)
    #print 'fixleft values:',fixleft_values
    #print 'popt in THz:',np.around(pin*x_to_THz,decimals=2)
    #gmodel = Model(lorenz_68_Hauer2015_2p_0p)
    #gmodel = Model(lorenz_68_Hauer2015_2p_2p)
    #print 'modelxx',model
    #np.savetxt("yyy",y)  # ist auf 1 skaliert
    #sys.exit()
    modeli = Model(model)
    if fixleft:
        wval = pin[0]
        Gval = pin[2]
        #print 'wval1;',wval,wval*x_to_THz
        #print 'Gval1;',Gval,Gval*x_to_THz
        if type(fixleft_values) != bool:
            wval = fixleft_values[0][0]/x_to_THz
            Gval = fixleft_values[1][0]/x_to_THz
            #print 'wval1;',wval,wval*x_to_THz
            #print 'Gval1;',Gval,Gval*x_to_THz
            pin[0] = wval
            pin[2] = Gval
        #print 'wval2;',wval,wval*x_to_THz
        #print 'Gval2;',Gval,Gval*x_to_THz
        #sys.exit()
        modeli.set_param_hint('w1', value=pin[0], vary=False)
        modeli.set_param_hint('a1', value=pin[1], vary=False)
        modeli.set_param_hint('G1', value=pin[2], vary=False)

        Gmax = pin[5]*5.0
        if Gmax * x_to_THz > 4.:
            Gmax = 4./x_to_THz
        modeli.set_param_hint('w2',  value=pin[3], min=pin[3]*mwl,max=pin[3]*mwr)
        modeli.set_param_hint('a2',  value=pin[4], min=pin[4]*mal,max=pin[4]*mar)
        modeli.set_param_hint('G2',  value=pin[5], min=pin[5]*mGl,max=pin[4]*mGr) # 0.8 is too much a constraint when few mdsteps
        psigma = np.array([pin[3],pin[4],pin[5]])
        if psigma[2]*x_to_THz < 0.3: psigma[2] = 0.5/x_to_THz
        #print 'psigma',psigma,psigma*x_to_THz

        if ratio_act and xind_peak_idx == 1:
            if verbose > 1:
                print(printblue('222: ratio (change right peak): '+str(ratio)))
            modeli.set_param_hint('w2',  value=pin[3], vary=False)
            modeli.set_param_hint('G2',  value=pin[5], vary=False) #
            if ratio < 0.2:
                modeli.set_param_hint('a2',  value=pin[4]*ratio, min=0,max=pin[4]/10.)
            if ratio > 8.:
                modeli.set_param_hint('a2',  value=pin[4]*ratio, min=0,max=pin[4]*10.)



    if fixright:
        Gmax = pin[2]*5.0
        Gmin = pin[2]*0.3
        if Gmax * x_to_THz > 4.:
            Gmax = 4./x_to_THz
        modeli.set_param_hint('w1',  value=pin[0], min=pin[0]*mwl,max=pin[0]*mwr)
        modeli.set_param_hint('a1',  value=pin[1], min=pin[1]*mal,max=pin[1]*mar)
        modeli.set_param_hint('G1',  value=pin[2], min=pin[2]*mGl,max=pin[2]*mGr) #

        psigma = np.array([pin[0],pin[1],pin[2]])
        if psigma[2]*x_to_THz < 0.3: psigma[2] = 0.5/x_to_THz

        modeli.set_param_hint('w2', value=pin[3], vary=False)
        modeli.set_param_hint('a2', value=pin[4], vary=False)
        modeli.set_param_hint('G2', value=pin[5], vary=False)

        if ratio_act and xind_peak_idx == 0:
            if verbose > 1:
                print(printblue('222: ratio (change left peak): '+str(ratio)))
            modeli.set_param_hint('w1',  value=pin[0], vary=False)
            modeli.set_param_hint('G1',  value=pin[2], vary=False) #
            if ratio < 0.2:
                modeli.set_param_hint('a1',  value=pin[1]*ratio, min=0,max=pin[1]/10.)
            if ratio > 8.:
                modeli.set_param_hint('a1',  value=pin[1]*ratio, min=0,max=pin[1]*10.)


    #print 'modeli.param_hints :',modeli.param_hints
    #print "param names        :",modeli.param_names
    #print 'independent_vars   :',modeli.independent_vars
    #print "modeli.func        :",modeli.func
    weights = lorenz_68_Hauer2015_1p_1p(x, *psigma)
    p=np.copy(pin)
    #print 'p',p
    result_in = modeli.fit(y, w=x, weights=weights, w1=p[0], a1=p[1], G1=p[2], w2=p[3], a2=p[4], G2=p[5])
    #sys.exit()
    popt = lmfit_result_to_popt(result_in,np.zeros(2))
    if verbose > 1:
        print('fix_l:',fixleft,'fix_r:',fixright,np.round(popt*x_to_THz,4))
    #sys.exit()
    #result = Model(gaussian).fit(y, x=x, amp=5, cen=5, wid=1)
    #print(result.fit_report())
    #print "----------",result.params['w01']
    #print 'aaaaaaaaaaa',result.params.pretty_print()
    #print 'bbbbbbbbbbb',result.best_values
    #d = result.best_values
    #for i in d:
    #    print i, d[i]
    #popt=np.abs(np.array([d.get('w01'),d.get('a1'),d.get('G1'),d.get('w02'),d.get('a2'),d.get('G2')]))
    #popt = np.array([result.best_values[k] for k in result.best_values])
    #if verbose:
    #    #if True:
    #    print 'popt >>>>>>:',np.around(popt*x_to_THz,decimals=2)

    #if plotit:
    #    plt.plot(x, y, 'bo')
    #    plt.plot(x, result.init_fit, 'k--')
    #    plt.plot(x, result.best_fit, 'r-')
    #    plt.xlim([0,2000])
    #    plt.ion()
    #    plt.show()
    ##return np.around(popt*x_to_THz,decimals=5)  # in THz
    return popt

def lmfit_result_to_popt(result,xind,indices=['w','a','G'],addindex=True):
    popt0=[]
    for i in np.arange(len(xind))+1:
        for j in indices:
            var=j
            if addindex:
                var = j+str(i)
            out=result.best_values.get(var)
            #print 'i:',i,var,out #np.round(out,2)#,'\t',np.round(out*x_to_THz,2)
            popt0.append(out)
    popt0 = np.array(popt0)
    return popt0

def fit_iterative(function,x,y,sigma=False,maxvew=10,strr=False,kwargs=False,args=False):
    ''' function is the defined fuction e.g. lorenz where lorenz == lorenz(w,w0,G,a)
    x and y are the xvalues and yvalues
    d is an array or a list containing the relative differences for the variable to converge to
    xind is an array with [[xpos1, ppos2, ...], x_to_THz]
    maxvew is the maximal number of iterations
    sigma can be either False or a int
            strr muss xa (whole powerspektrum left/right)
            oder xl (only left part of ps) sein!
    '''
    if args.verbose:
        print()
        print("************************* fit_iterative ************************")
    acc=14
    if np.round(y[0],acc) != np.round(y[-1],acc):
        print('y[0] ',y[0])
        print('y[-1]',y[-1])
        print('diff',y[0]-y[-1])
        print('diff',np.round(y[0],acc)-np.round(y[-1],acc))
        print(y[0],y[1],y[2])
        print(y[-1],y[-2],y[-3])
        print()
        print(y[1])
        print(y[-1])
        print()
        print(y[:4])
        print(y[-4:])
        sys.exit("PROBLEM 1, FUNCTION IS NOT SYMMETRIC")
    if np.round(y[1],acc) != np.round(y[-2],acc):
        print('y[1] ',y[1])
        print('y[-2]',y[-2])
        print('diff',y[1]-y[-2])
        print('diff',np.round(y[1],acc)-np.round(y[-2],acc))
        print()
        print(y[:4])
        print(y[-4:])
        print()
        print('y[0:3]',y[0],y[1],y[2])
        print('y[neg]',y[-1],y[-2],y[-3])
        print()
        print('y[0:3]',y[0:3])
        print('y[neg]',y[-3:][::-1])
        print()
        print(y[1])
        print(y[-1])
        print()
        print(y[:4])
        print(y[-4:])
        print('---')
        sys.exit("PROBLEM 2, FUNCTION IS NOT SYMMETRIC")
    xind = kwargs[0]
    qpointstr = kwargs[1]
    xl = kwargs[2]
    yl = kwargs[3]
    xa = kwargs[4]
    ya = kwargs[5]
    verbose = kwargs[6]
    mdsteps = kwargs[7]
    check_qpoint_for_crossing = kwargs[8]
    qpoint  = kwargs[9]
    psmax = kwargs[10]
    x_to_THz = kwargs[11]
    space_fft_which = kwargs[12]

    path = qpoint_get_l_0_0_or_l_l_0_or_l_l_l_or_t_0_0_or_t_t_t_or_t1_t1_t1_or_t2_t2_t2_from_qpoint(qpoint,args)
    crossing_fq_lt = kwargs[-3]
    filename = kwargs[-2]
    take = kwargs[-1]
    if len(xind) == 1:
        if xind[0] == 0:
            take2 = str(len(xind))+'p_'+str(len(xind))+'p'
            function = eval(take+take2)
            return False,False,function

    #print "check_qpoint_for_crossing",check_qpoint_for_crossing,qpoint
    if strr != "xa" and strr != "xl":
        sys.eixt("strr has to be xl or xa")
    if type(sigma) == bool and sigma == True:
        sigma=1.

    #print "check_qpoint_for_crossing:",check_qpoint_for_crossing
    #print "crossing_fq_lt:",crossing_fq_lt
    #sys.exit()
    if type(sigma) != bool and type(sigma) != int:
        print("type(sigma):",type(sigma))
        sys.exit('sigma has to be either False(==default) or an int')


    ################################################################
    # when several peaks: add bounds in general
    # also when one peak! otherwise it might not converge
    ################################################################
    if len(xind) < 1:
        sys.exit("len(xind) needs to be at least one")


    ######################################################################
    # INIT GROB FREQUENCY: more or less fix the frequency
    # INIT GROB LIFETIME : from 0 to 3
    # INIT GROB a parameter: to get this is the aim of the game
    ######################################################################
    do_estimate_from_db = True # do this only for G (not for a)  # THIS IS GENERALLY TRUE
    write_for_debugging = args.write_for_debugging
    if do_estimate_from_db == True:
        xind = np.array([6.43,15.17,16.72])/x_to_THz           # 4 2 0  # CrN ASD 300K 222
        xind = np.array([9.6,15.2,16.7])/x_to_THz           # 2 2 4  # CrN ASD 300K 222
        xind = np.array([6.5,9.5,15.45,17.0])/x_to_THz      # 2 2 4  # CrN DLM 1000K 222
        xind = np.array([6.2,9.4,15.2,16.8])/x_to_THz      # 2 2 4  # CrN
        xind = np.array([3.8,14.9])/x_to_THz               # 6 1 0   333 TiN 1000K
        xind = np.array([7.,15.16,16.12])/x_to_THz         # 6 3 0   333 TiN 1000K
        xind = np.array([5.38,15.5])/x_to_THz              # 7 5 0       TiN
        xind = np.array([7.85,15.16,16.16])/x_to_THz       # 3 3 6  # 1000K
        xind = np.array([7.4,7.6,15.2,15.7])/x_to_THz      # 3 3 6  # 300K
        xind = np.array([7.3,8.8,15.8,16.9])/x_to_THz      # 2 2 4  # 1000K 222 TiN
        xind = np.array([7.4,8.2,8.55])/x_to_THz      # 3 3 8  # 300K Al
        xind = np.array([5.1,8.8])/x_to_THz      # 4 4 8  # 900K Al
        xind = np.array([2.83,0.036,1.2])/x_to_THz           # qp=8 8 8 for dominiques Ti bcc
        xind = np.array([5.3])/x_to_THz           # qp=8 8 8 for dominiques Ti bcc
    #if do_estimate_from_db == False:
    #    esto = xind
    #    esto_THz = xind*x_to_THz

    verycertainestimate = False
    qlen = qpoint_to_length_qpoint(qpoint[0],qpoint[1],qpoint[2],args,get_lll_double_length=True)
    xind_init = np.copy(xind)

    ############################################################
    # load sigma0 from db if available
    ############################################################
    filenamesigma0 = filename.split("sigma_")[0]+"sigma_0"
    from_db = get_from_db_w_a_G(qpoint,filenamesigma0,args.peaks,mdsteps)
    ## change from db manually
    #from_db = np.array([2.83,0.036,1.2])
    #from_db = np.array([2.65,0.034,1.45])
    dp = np.repeat(np.array([[0.01,0.01,0.01]]), len(xind),axis=0).flatten()

    if strr == 'xa':
        take2 = str(len(xind))+'p_'+str(len(xind))+'p'
        function = eval(take+take2)

    if strr == 'xl':
        take2 = str(len(xind))+'p_0p'
        function = eval(take+take2)
    gmodel = Model(function)

    if args.verbose:
        print("**** xind (1)        :",xind,"-> (USED FOR GROB):",np.round(xind*x_to_THz,2),"(THz)")
        print('**** len(xind)       :',len(xind))
        print('**** x_to_THz        :',x_to_THz)
        print('**** qlen            :',qlen)
        print('**** qpoint          :',qpoint)
        print('**** path            :',path)
        print('**** psmax           :',psmax)
        print('**** take            :',take)
        print('**** function        :',function)
        print('**** mdsteps         :',mdsteps)
        print('**** sigma           :',sigma)
        print('**** strr            :',strr)
        print('**** space_fft_which :',space_fft_which)
        print('**** crossing        :',check_qpoint_for_crossing)
        print('**** write_for_debugg:',args.write_for_debugging)
        print('**** LOADED DB       :',from_db)
        print('**** ANMERKUNG       : if LOADED from DB freqs significantly different than the xind && md dependence has significant wiggles, try better the one from xind!')






    haveestimatesigma0 = False
    if type(from_db) != bool:
        esto_THz = np.copy(from_db)

        def check_if_xind_and_from_db_are_same(esto_THz,xind,x_to_THz):
            #print 'ee',esto_THz
            #print 'ee',esto_THz[0::3]
            #print 'ex',xind*x_to_THz
            freqs_db = esto_THz[0::3]
            freqs_xin = xind*x_to_THz
            if len(freqs_db == freqs_xin):
                aa = np.isclose(freqs_db,freqs_xin,rtol=0.05)
                #print aa,type(aa),aa.all()
                return aa.all()
            else:
                return False


    # a) check if freqs dfifferent from xind and from_db (this can be done all the time)
        same = check_if_xind_and_from_db_are_same(esto_THz,xind,x_to_THz)
        #same = True
    # b) if a) shows differences, check if convergence as function of mdsteps exist.
        if same == False:
            esto_THz = False
            do_estimate_from_db = False
            if mdsteps > 2000:
                print(printred('**** LOADED DB       : seems not OK!, mdsteps '+str(mdsteps)+" sigma "+str(sigma)),"xind",xind*x_to_THz,"from DB:",esto_THz)
            #print_and_save_lifetimesgood(args = args,printtoscreen=False,return_mdconvergence=[qpoint,sigma])
    # c) if b) is yes -> take xind only and start from there!
        else:
            esto = esto_THz/x_to_THz
            popt0 = np.copy(esto)
            verycertainestimate = True
            haveestimatesigma0 = True



    #########################################################################
    # have not the exact point in the DB for sigma0 ... DO ESTIMATE FROM DB
    #########################################################################
    if haveestimatesigma0 == False:
        esto_THz = False
        esto     = False
        if args.verbose:
            print('**** do_est_from_db  :',do_estimate_from_db)
            if do_estimate_from_db:
                print('**** filename        :',filename)
                # filename e.g.:    lorenz_68_Hauer2015_xa_sigma_1
                #                   lorenz_6886_Hauer2015_xa_sigma_1
            if write_for_debugging:
                print(printred('**** write_for_debugging ****'))


        if do_estimate_from_db:
            esto_THz = get_estimate_w_a_G(qpoint,filename,mdsteps,args.verbose,args=args,add_filename=filename,x_to_THz=x_to_THz)
            #print 'es01',esto_THz
            if type(esto_THz) != bool:
                #esto_THz[0] = 5.03
                #if take == 'lorenz_6886_Hauer2015_':
                #    esto_THz[2]=0.11
                if np.count_nonzero(esto_THz) == 0 or np.count_nonzero(esto_THz) != len(esto_THz):
                    esto_THz = False
                if len(xind)*3 != len(esto_THz):
                    esto_THz = False
                if type(esto_THz) != bool:
                    esto = esto_THz/x_to_THz
                    xind = esto[::3]
        if args.verbose: # no matter if it is type bool or not
            print('**** GOT FROM DB     :',esto_THz) # is already in THz
            #if type(esto_THz) == bool:
            #    sys.exit('222222333')

        if type(esto) != bool:
            popt0 = np.copy(esto)
            verycertainestimate = True
            haveestimatesigma0 = True
        else:
            haveestimatesigma0 = False


    #########################################################################
    # still no estimate from DB. DO ESTIMATE FROM DB with lorenz_68_Hauer2015_xa_sigma_1
    #########################################################################
    if haveestimatesigma0 == False:
        esto_2     = False
        esto_THz_2 = False
        if do_estimate_from_db and type(esto_THz) == bool and filename == "lorenz_6886_Hauer2015_xa_sigma_1":
            esto_THz_2 = get_estimate_w_a_G(qpoint,"lorenz_68_Hauer2015_xa_sigma_1",mdsteps,args.verbose,args=args,add_filename=filename,x_to_THz=x_to_THz)
            #print 'es02',esto_THz
            if np.count_nonzero(esto_THz_2) == 0 or np.count_nonzero(esto_THz_2) != len(esto_THz_2):
                esto_THz_2 = False
            if len(xind)*3 != len(esto_THz):
                esto_THz = False
            if type(esto_THz_2) != bool:
                esto_2 = esto_THz_2/x_to_THz
                xind = esto_2[::3]

        if type(esto_2) != bool:
            popt0 = np.copy(esto_2)
            verycertainestimate = True
            haveestimatesigma0 = True
        else:
            haveestimatesigma0 = False



    #print("+++++++++++5 write_for_debugging:",write_for_debugging)
    #print("+++++++++++5 haveestimatesigma0 :",haveestimatesigma0)

    #########################################################################
    # if still no estimate from DB, try to do the fit
    #########################################################################
    if haveestimatesigma0 == False:
        ######################################
        # for t2_t2_0 get estimated freqs
        ######################################

        if args.verbose > 1:
            print('**** KKK space_fft_which',space_fft_which)
        show = False
        if mdsteps == args.mdstepstocalc_all.max(): show = True  # only when last mdstep
        if type(args.mdsteps) != bool: show = True
        #if show:
        if False:
            if type(esto_THz) != bool:
                print('**** esto_THz (db)   :',np.round(esto_THz,4))
            else:
                print('**** esto_THz (db)   :',esto_THz)
        #sys.exit()
        if args.verbose:
            if type(esto_THz) != bool:
                if args.verbose: # > 1:
                    print('**** take estimate db:',esto_THz)
                    print("**** xind (2)        :",xind,"-> (USED FOR GROB):",np.round(xind*x_to_THz,2),"(THz)")

        if args.verbose:
            if type(esto_THz_2) != bool:
                print('**** esto_THz_2 (db) :',np.round(esto_THz_2,4))
            else:
                print('**** esto_THz_2 (db) :',esto_THz_2)
        #sys.exit('345677')

        #print 'tya',type(esto)
        if mdsteps > 20000 and type(esto) == np.ndarray:
            verycertainestimate = True
        if args.verbose:
            print('**** verycertainestimate:',verycertainestimate)
        if args.verbose > 1:
            print("**** xind             ---> :",np.round(xind,3),"(USED FOR GROB):",np.round(xind*x_to_THz,3),"(THz)")
        if args.verbose > 1:
            print('**** len(xind)        ---> :',len(xind))
        #print("+++++++++++5 write_for_debugging:",write_for_debugging)
        if write_for_debugging:
            #print 'x:',x
            #print 'y:',y
            print('**** writing INIT_A_ORIG_POWERSPECTRUM',(y*psmax).max())
            parnpsavetxt("INIT_A_ORIG_POWERSPECTRUM_y_noscale_is_psmax__we_fit_to_y",np.transpose(np.array([x*x_to_THz,y*psmax])),fmt='%.8f %.5f')
            #parnpsavetxt("INIT_A_ORIG_POWERSPECTRUM_y_scaled_to_1",np.transpose(np.array([x*x_to_THz,y])),fmt='%.8f %.5f')

        #########################################################################
        # hier haben wir ein belastbares xind welches die richtige laenge hat
        #########################################################################

        if write_for_debugging and type(esto) != bool:
            #print 'x:',x
            #print 'y:',y
            print('**** writing INIT_B_DB:',do_estimate_from_db)
            out = function(x, *esto)
            #parnpsavetxt("INIT_DB_EST",np.transpose(np.array([x*x_to_THz,y*psmax])),fmt='%.8f %.5f')
            parnpsavetxt("INIT_B_DB",np.transpose(np.array([x*x_to_THz,out*psmax])),fmt='%.8f %.5f')

        if len(xind) == 3:
            gmodel = Model(lorenz_68_Hauer2015_3p_3p)
            function = eval("lorenz_68_Hauer2015_3p_3p")
            if args.verbose > 1:
                print('**** function        :',function)
        #sys.exit('have a good estimate?')



        #######################################################################
        # try to get a sensible estimate if nothing is in the DB
        #######################################################################
        params = gmodel.make_params()
        if args.verbose > 3:
            print('params0-----------------------')
            my_print_params(params,x_to_THz,rt=4)

        G_est = (qlen*2.)/x_to_THz     # das ist nur eine Hausnummer!  (stimmt nur wenn ein peak)
        if G_est == 0:
            G_est = 0.5/x_to_THz
        amin=0.03
        if len(xind) ==1:
            amin=0.003
        if len(xind) > 0:
            #print 'it seems len(xind) > 0 is ok'
            #params.add('w1', value=xind[0],         min=xind[0]*0.999,  max=xind[0]*1.001)
            params.add('w1', value=xind[0],         min=0,vary=False)
            params.add('a1', value=0.02/x_to_THz,   min=amin,           max=np.inf)
            params.add('G1', value=G_est,           min=0,              max=3.5/x_to_THz)
        if len(xind) > 1:
            #params.add('w2', value=xind[1],         min=xind[1]*0.999,  max=xind[1]*1.001)
            params.add('w2', value=xind[1],         min=0,vary=False)
            params.add('a2', value=0.02/x_to_THz,   min=amin,              max=np.inf)
            params.add('G2', value=G_est,           min=0,              max=3.5/x_to_THz)
        if len(xind) > 2:
            #params.add('w3', value=xind[2],         min=xind[2]*0.999,  max=xind[2]*1.001)
            params.add('w3', value=xind[2],         min=0,vary=False)
            params.add('a3', value=0.02/x_to_THz,   min=amin,              max=np.inf)
            params.add('G3', value=G_est,           min=0,              max=3.5/x_to_THz)
        if len(xind) > 3:
            #params.add('w4', value=xind[3],         min=xind[3]*0.999,  max=xind[3]*1.001)
            params.add('w4', value=xind[3],         min=0,vary=False)
            params.add('a4', value=0.02/x_to_THz,   min=amin,              max=np.inf)
            params.add('G4', value=G_est,           min=0,              max=3.5/x_to_THz)



        #result = gmodel.fit(y, params, w=x)
        ##print(result.fit_report())
        #popt0 = lmfit_result_to_popt(result,xind)
        #outorig_y = function(x, *popt0)
        #xout,yout = ps_to_ps_for_writing(outorig_y,xind,x_to_THz,function,x,popt0)
        #if True:
        #    print 'writing INIT_B_GROB',np.round(popt0*x_to_THz,3)
        #    parnpsavetxt("INIT_B_GROB",np.transpose(np.array([xout,yout*psmax])),fmt='%.8f %.15f')

        if args.verbose > 3:
            print('params1-----------------------')
            my_print_params(params,x_to_THz,rt=4)
            #print 'xind',xind
            #print 'esto',esto
        if type(esto_THz) != bool:
            if len(xind) > 0:
                G_max = get_Gmax(esto[2]*1.2,x_to_THz)
                params.add('w1', value=esto[0], vary=False)
                params.add('a1', value=esto[1], min=esto[1]/10.,    max=esto[1]*10.1) # 10.1 so taht min != max
                params.add('G1', value=esto[2], min=esto[2]*0.8,    max=G_max)
                if check_qpoint_for_crossing:
                    params.add('G1', value=esto[2],min=esto[2]*0.1,max=G_max)
                if verycertainestimate:
                    params.add('G1', value=esto[2],min=esto[2]*0.8,max=esto[2]*1.2)

            if len(xind) > 1:
                G_max = get_Gmax(esto[5]*1.2,x_to_THz)
                params.add('w2', value=esto[3], vary=False)
                params.add('a2', value=esto[4], min=esto[4]/10.,    max=esto[4]*10.1)
                params.add('G2', value=esto[5], min=esto[5]*0.8,    max=G_max)
            if len(xind) > 2:
                G_max = get_Gmax(esto[8]*1.2,x_to_THz)
                params.add('w3', value=esto[6], vary=False)
                params.add('a3', value=esto[7], min=esto[7]/10.,    max=esto[7]*10.1)
                params.add('G3', value=esto[8], min=esto[8]*0.8,    max=G_max)
            if len(xind) > 3:
                G_max = get_Gmax(esto[11]*1.2,x_to_THz)
                params.add('w4', value=esto[9], vary=False)
                params.add('a4', value=esto[10], min=esto[10]/10.,    max=esto[10]*10.1)
                params.add('G4', value=esto[11], min=esto[11]*0.8,    max=G_max)

        if args.verbose > 3:
            print('params2-----------------------')
            my_print_params(params,x_to_THz,rt=4)

        if type(esto_THz) == bool and type(esto_THz_2) != bool:
            if len(xind) > 0:
                G_max = get_Gmax(esto_2[2]*1.2,x_to_THz)
                params.add('w1', value=esto_2[0], vary=False)
                print('now a1:',1./x_to_THz)
                params.add('a1', value=1./x_to_THz, min=0.,  max=np.inf) # 10.1 so taht min != max
                params.add('G1', value=esto_2[2], min=esto_2[2]*0.8,    max=G_max)
                if check_qpoint_for_crossing:
                    #print "CROSSING"
                    params.add('G1', value=esto_2[2], min=esto_2[2]*0.3,    max=G_max)

            if len(xind) > 1:
                G_max = get_Gmax(esto[5]*1.2,x_to_THz)
                params.add('w2', value=esto_2[3], vary=False)
                params.add('a2', value=1./x_to_THz, min=0.,    max=np.inf)
                params.add('G2', value=esto_2[5], min=esto_2[5]*0.8,    max=G_max)
            if len(xind) > 2:
                G_max = get_Gmax(esto[8]*1.2,x_to_THz)
                params.add('w3', value=esto_2[6], vary=False)
                params.add('a3', value=1./x_to_THz, min=0.,    max=np.inf)
                params.add('G3', value=esto_2[8], min=esto_2[8]*0.8,    max=G_max)
            if len(xind) > 3:
                G_max = get_Gmax(esto[11]*1.2,x_to_THz)
                params.add('w4', value=esto_2[9], vary=False)
                params.add('a4', value=1./x_to_THz, min=0.,    max=np.inf)
                params.add('G4', value=esto_2[11], min=esto_2[11]*0.8,    max=G_max)






        params = fix_params_long_and_t1_for_t2(params,space_fft_which,path,qpoint,filename,mdsteps,x_to_THz,estimate_w2_a2_G2=True)

        if args.verbose > 2:
            print('params esto2-----------------------')
            my_print_params(params,x_to_THz,rt=4)
        #sys.exit()
        #print('y')
        #print(y)
        #print('x')
        #print(x)
        #np.savetxt('y.dat',y)
        #np.savetxt('x.dat',x)
        #print('ymin',y.min())
        result = gmodel.fit(y, params, w=x)
        if args.verbose > 2:
            print('params esto2 fit_report()-----------------------')
            print((result.fit_report()))
        popt0 = lmfit_result_to_popt(result,xind)
        outorig_y = function(x, *popt0)
        xout,yout = ps_to_ps_for_writing(outorig_y,xind,x_to_THz,function,x,popt0,args=args)
        if args.verbose > 1:
            print("**** ----> INIT_C   ",np.round(popt0*x_to_THz,4))
        if write_for_debugging:
            print('writing INIT_C_GROB ')#,np.round(popt0*x_to_THz,3)
            #parnpsavetxt("INIT_C_GROB",np.transpose(np.array([xout,yout*psmax])),fmt='%.8f %.15f')
            parnpsavetxt("INIT_C_GROB",np.transpose(np.array([x*x_to_THz,outorig_y*psmax])),fmt='%.8f %.15f')

        ratio = y.max()/yout.max()  # get also y.mean
        if ratio <= 0.12 or ratio >= 3.5:  # ratio of 3 is ok when there is a big scatter! also check if the mean is scattered as this!
            if args.verbose > 1:
                print(printred('111: ratio: '+str(ratio)+" mdsteps: "+str(mdsteps)+" for > 1 origfunction is larger;"))
                print("y           max",y.max())
                print("yout (fit)  max",yout.max())
                print("popt0, hmmmm. for 6886 it is necessary to use weighting around maximum for all fits")

            if ratio >= 3.5:
                params['a1'].max = np.inf
                params['a1'].value = 1./x_to_THz
                if len(xind) > 1:
                    params['a2'].max = np.inf
                    params['a2'].value = 1./x_to_THz
                if len(xind) > 2:
                    params['a3'].max = np.inf
                    params['a3'].value = 1./x_to_THz
                if len(xind) > 3:
                    params['a4'].max = np.inf
                    params['a3'].value = 1./x_to_THz
                #jif len(xind) == 1:
                #G_max = get_Gmax(esto[2]*1.2,x_to_THz)
                #params.add('w1', value=esto[0], vary=False)
                #params.add('a1', value=1.     , min=0,    max=np.inf) # 10.1 so taht min != max
                #params.add('G1', value=esto[2], min=esto[2]*0.8,    max=G_max)

                if args.verbose > 3:
                    print('params33-----------------------')
                    my_print_params(params,x_to_THz,rt=4)
                result = gmodel.fit(y, params, w=x)
                #print(result.fit_report())
                popt0 = lmfit_result_to_popt(result,xind)
                outorig_y = function(x, *popt0)
                xout,yout = ps_to_ps_for_writing(outorig_y,xind,x_to_THz,function,x,popt0,args=args)
                if args.verbose > 1:
                    print("**** ----> INIT_CC  ",np.round(popt0*x_to_THz,4))
                ratio = y.max()/yout.max()  # get also y.mean
                if ratio <= 0.2 or ratio >= 3.5:  # ratio of 3 is ok when there is a big scatter! also check if the mean is scattered as this!
                    if args.verbose > 1:
                        print(printred('111: ratio: '+str(ratio)+" mdsteps: "+str(mdsteps)))
                        print("y           max",y.max())
                        print("yout (fit)  max",yout.max())

        #xind_peak_THz_curr = xout[np.where(yout==yout.max())[0][0]]
        #if args.verbose > 1:
        #    print 'xind_peak_THz_curr',xind_peak_THz_curr
        #xind_peak_THz = (find_nearest(xind*x_to_THz, xind_peak_THz_curr))
        #if args.verbose > 1:
        #    print 'xind_peak_THz',xind_peak_THz
        #xind_peak_idx = np.where(xind*x_to_THz==xind_peak_THz)[0][0]
        #if args.verbose > 1:
        #    print 'xind_peak_idx',xind_peak_idx


        ###########################################################################
        # if a is out of bounds, MAKE THIS ALSO IF a IS NOT OUT OF BOUNDS
        # keep everything but a
        # often it only makes sense to vary a and G
        ###########################################################################
        #if ratio <= 0.2 or ratio >= 3. and type(esto_THz) != bool:
        if True:
            # adapt a
            if len(xind) > 0:
                params.add('w1', value=popt0[0], vary=False)
                params.add('a1', value=popt0[1], min=0,              max=np.inf)
                params.add('G1', value=popt0[2], vary=False)
            if len(xind) > 1:
                params.add('w2', value=popt0[3], vary=False)
                params.add('a2', value=popt0[4], min=0,              max=np.inf)
                params.add('G2', value=popt0[5], vary=False)
            if len(xind) > 2:
                params.add('w3', value=popt0[6], vary=False)
                params.add('a3', value=popt0[7], min=0,              max=np.inf)
                params.add('G3', value=popt0[8], vary=False)
            if len(xind) > 3:
                params.add('w4', value=popt0[9], vary=False)
                params.add('a4', value=popt0[10],min=0,              max=np.inf)
                params.add('G4', value=popt0[11],vary=False)
            #params = fix_params_long_and_t1_for_t2(params,space_fft_which,path,qpoint,filename,mdsteps,x_to_THz)
            result = gmodel.fit(y, params, w=x)
            popt0 = lmfit_result_to_popt(result,xind)
            if args.verbose > 1:
                print("**** ----> INIT_D1   ",np.round(popt0*x_to_THz,4))

            # adapt G
            if len(xind) > 0:
                Gmin = popt0[2]*0.099
                Gmax = get_Gmax(popt0[2]*3.,x_to_THz)
                params.add('w1', value=popt0[0], vary=False)
                params.add('a1', value=popt0[1], vary=False)
                params.add('G1', value=popt0[2], min=Gmin , max=Gmax)
                if verycertainestimate:
                    Gmin = popt0[2]*0.8
                    Gmax = get_Gmax(popt0[2]*1.2,x_to_THz)
                    params.add('G1', value=popt0[2], min=Gmin , max=Gmax)
            if len(xind) > 1:
                Gmin = popt0[5]*0.333
                Gmax = get_Gmax(popt0[5]*3.,x_to_THz)
                params.add('w2', value=popt0[3], vary=False)
                params.add('a2', value=popt0[4], vary=False)
                params.add('G2', value=popt0[5], min=Gmin , max=Gmax)
                if verycertainestimate:
                    Gmincertain = popt0[5]*0.8
                    Gmaxcertain = get_Gmax(popt0[5]*1.2,x_to_THz)
                    params.add('G2', value=popt0[5], min=Gmincertain, max=Gmaxcertain)
            if len(xind) > 2:
                Gmin = popt0[8]*0.333
                Gmax = get_Gmax(popt0[8]*3.,x_to_THz)
                params.add('w3', value=popt0[6], vary=False)
                params.add('a3', value=popt0[7], vary=False)
                params.add('G3', value=popt0[8], min=Gmin , max=Gmax)
                if verycertainestimate:
                    Gmincertain = popt0[8]*0.8
                    Gmaxcertain = get_Gmax(popt0[8]*1.2,x_to_THz)
                    params.add('G3', value=popt0[8], min=Gmincertain, max=Gmaxcertain)
            if len(xind) > 3:
                Gmin = popt0[11]*0.333
                Gmax = get_Gmax(popt0[11]*3.,x_to_THz)
                params.add('w4', value=popt0[9], vary=False)
                params.add('a4', value=popt0[10],vary=False)
                params.add('G4', value=popt0[11], min=Gmin, max=Gmax)
                if verycertainestimate:
                    Gmincertain = popt0[11]*0.8
                    Gmaxcertain = get_Gmax(popt0[11]*1.2,x_to_THz)
                    params.add('G4', value=popt0[11], min=Gmincertain, max=Gmaxcertain)
            #params = fix_params_long_and_t1_for_t2(params,space_fft_which,path,qpoint,filename,mdsteps,x_to_THz)
            result = gmodel.fit(y, params, w=x)
            popt0 = lmfit_result_to_popt(result,xind)
            if args.verbose > 1:
                print("**** ----> INIT_D2   ",np.round(popt0*x_to_THz,4))

            # adapt a
            if len(xind) > 0:
                params.add('w1', value=popt0[0], vary=False)
                params.add('a1', value=popt0[1], min=0,              max=np.inf)
                params.add('G1', value=popt0[2], vary=False)
            if len(xind) > 1:
                params.add('w2', value=popt0[3], vary=False)
                params.add('a2', value=popt0[4], min=0,              max=np.inf)
                params.add('G2', value=popt0[5], vary=False)
            if len(xind) > 2:
                params.add('w3', value=popt0[6], vary=False)
                params.add('a3', value=popt0[7], min=0,              max=np.inf)
                params.add('G3', value=popt0[8], vary=False)
            if len(xind) > 3:
                params.add('w4', value=popt0[9], vary=False)
                params.add('a4', value=popt0[10],min=0,              max=np.inf)
                params.add('G4', value=popt0[11],vary=False)
            params = fix_params_long_and_t1_for_t2(params,space_fft_which,path,qpoint,filename,mdsteps,x_to_THz)
            result = gmodel.fit(y, params, w=x)
            popt0 = lmfit_result_to_popt(result,xind)
            if args.verbose > 1:
                print("**** ----> INIT_D3   ",np.round(popt0*x_to_THz,4),"len(xind)",len(xind))



            #sys.exit('yola789-')
            # adapt G
            if len(xind) > 0:
                Gmin = popt0[2]*0.1
                Gmax = get_Gmax(popt0[2]*3.,x_to_THz)
                params.add('w1', value=popt0[0], vary=False)
                params.add('a1', value=popt0[1], vary=False)
                params.add('G1', value=popt0[2], min=Gmin          , max=Gmax)
                if verycertainestimate:
                    Gmin = popt0[2]*0.8
                    Gmax = get_Gmax(popt0[2]*1.2,x_to_THz)
                    params.add('G1', value=popt0[2], min=Gmin , max=Gmax)
            if len(xind) > 1:
                Gmin = popt0[5]*0.333
                Gmax = get_Gmax(popt0[5]*3.,x_to_THz)
                params.add('w2', value=popt0[3], vary=False)
                params.add('a2', value=popt0[4], vary=False)
                params.add('G2', value=popt0[5], min=Gmin          , max=Gmax       )
                if verycertainestimate:
                    Gmin = popt0[5]*0.8
                    Gmax = get_Gmax(popt0[5]*1.2,x_to_THz)
                    params.add('G2', value=popt0[5], min=Gmin , max=Gmax)
            if len(xind) > 2:
                Gmin = popt0[8]*0.333
                Gmax = get_Gmax(popt0[8]*3.,x_to_THz)
                params.add('w3', value=popt0[6], vary=False)
                params.add('a3', value=popt0[7], vary=False)
                params.add('G3', value=popt0[8], min=Gmin          , max=Gmax       )
            if len(xind) > 3:
                Gmin = popt0[11]*0.333
                Gmax = get_Gmax(popt0[11]*3.,x_to_THz)
                params.add('w4', value=popt0[9], vary=False)
                params.add('a4', value=popt0[10],vary=False)
                params.add('G4', value=popt0[11], min=Gmin           , max=Gmax      )
            params = fix_params_long_and_t1_for_t2(params,space_fft_which,path,qpoint,filename,mdsteps,x_to_THz)
            result = gmodel.fit(y, params, w=x)
            popt0 = lmfit_result_to_popt(result,xind)

            # adapt a and G
            if len(xind) > 0:
                params.add('w1', value=popt0[0], vary=False)
                params.add('a1', value=popt0[1], min=popt0[1]/10.,max=popt0[1]*10)
                params.add('G1', value=popt0[2], min=popt0[2]/10.,max=popt0[2]*10)
                if verycertainestimate:
                    Gmin = popt0[2]*0.8
                    Gmax = get_Gmax(popt0[2]*1.2,x_to_THz)
                    params.add('G1', value=popt0[2], min=Gmin , max=Gmax)
            if len(xind) > 1:
                params.add('w2', value=popt0[3], vary=False)
                params.add('a2', value=popt0[4], vary=False)
                params.add('G2', value=popt0[5], vary=False)
                if verycertainestimate:
                    Gmin = popt0[5]*0.8
                    Gmax = get_Gmax(popt0[5]*1.2,x_to_THz)
                    params.add('G2', value=popt0[5], min=Gmin , max=Gmax)
            if len(xind) > 2:
                params.add('w3', value=popt0[6], vary=False)
                params.add('a3', value=popt0[7], vary=False)
                params.add('G3', value=popt0[8], vary=False)
            if len(xind) > 3:
                params.add('w4', value=popt0[9], vary=False)
                params.add('a4', value=popt0[10],vary=False)
                params.add('G4', value=popt0[11],vary=False)
            params = fix_params_long_and_t1_for_t2(params,space_fft_which,path,qpoint,filename,mdsteps,x_to_THz)
            result = gmodel.fit(y, params, w=x)
            popt0 = lmfit_result_to_popt(result,xind)
            if args.verbose > 1:
                print("**** ----> INIT_D4   ",np.round(popt0*x_to_THz,4),"len(xind)",len(xind))



            # adapt a
            if len(xind) > 0:
                params.add('w1', value=popt0[0], vary=False)
                params.add('a1', value=popt0[1], min=0,              max=np.inf)
                params.add('G1', value=popt0[2], vary=False)
            if len(xind) > 1:
                params.add('w2', value=popt0[3], vary=False)
                params.add('a2', value=popt0[4], min=0,              max=np.inf)
                params.add('G2', value=popt0[5], vary=False)
            if len(xind) > 2:
                params.add('w3', value=popt0[6], vary=False)
                params.add('a3', value=popt0[7], min=0,              max=np.inf)
                params.add('G3', value=popt0[8], vary=False)
            if len(xind) > 3:
                params.add('w4', value=popt0[9], vary=False)
                params.add('a4', value=popt0[10],min=0,              max=np.inf)
                params.add('G4', value=popt0[11],vary=False)
            params = fix_params_long_and_t1_for_t2(params,space_fft_which,path,qpoint,filename,mdsteps,x_to_THz)
            result = gmodel.fit(y, params, w=x)
            popt0 = lmfit_result_to_popt(result,xind)

            # adapt G
            if len(xind) > 0:
                Gmin = popt0[2]*0.1
                Gmax = get_Gmax(popt0[2]*3.,x_to_THz)
                params.add('w1', value=popt0[0], vary=False)
                params.add('a1', value=popt0[1], vary=False)
                params.add('G1', value=popt0[2], min=Gmin          , max=Gmax       )
            if len(xind) > 1:
                Gmin = popt0[5]*0.333
                Gmax = get_Gmax(popt0[5]*3.,x_to_THz)
                params.add('w2', value=popt0[3], vary=False)
                params.add('a2', value=popt0[4], vary=False)
                params.add('G2', value=popt0[5], min=Gmin          , max=Gmax       )
            if len(xind) > 2:
                Gmin = popt0[8]*0.333
                Gmax = get_Gmax(popt0[8]*3.,x_to_THz)
                params.add('w3', value=popt0[6], vary=False)
                params.add('a3', value=popt0[7], vary=False)
                params.add('G3', value=popt0[8], min=Gmin          , max=Gmax       )
            if len(xind) > 3:
                Gmin = popt0[11]*0.333
                Gmax = get_Gmax(popt0[11]*3.,x_to_THz)
                params.add('w4', value=popt0[9], vary=False)
                params.add('a4', value=popt0[10],vary=False)
                params.add('G4', value=popt0[11], min=Gmin           , max=Gmax      )
            params = fix_params_long_and_t1_for_t2(params,space_fft_which,path,qpoint,filename,mdsteps,x_to_THz)
            result = gmodel.fit(y, params, w=x)
            popt0 = lmfit_result_to_popt(result,xind)

            # adapt a
            if len(xind) > 0:
                params.add('w1', value=popt0[0], vary=False)
                params.add('a1', value=popt0[1], min=0,              max=np.inf)
                params.add('G1', value=popt0[2], vary=False)
            if len(xind) > 1:
                params.add('w2', value=popt0[3], vary=False)
                params.add('a2', value=popt0[4], min=0,              max=np.inf)
                params.add('G2', value=popt0[5], vary=False)
            if len(xind) > 2:
                params.add('w3', value=popt0[6], vary=False)
                params.add('a3', value=popt0[7], min=0,              max=np.inf)
                params.add('G3', value=popt0[8], vary=False)
            if len(xind) > 3:
                params.add('w4', value=popt0[9], vary=False)
                params.add('a4', value=popt0[10],min=0,              max=np.inf)
                params.add('G4', value=popt0[11],vary=False)
            params = fix_params_long_and_t1_for_t2(params,space_fft_which,path,qpoint,filename,mdsteps,x_to_THz)
            result = gmodel.fit(y, params, w=x)
            popt0 = lmfit_result_to_popt(result,xind)

            # adapt G
            if len(xind) > 0:
                Gmin = popt0[2]*0.1
                Gmax = get_Gmax(popt0[2]*3.,x_to_THz)
                params.add('w1', value=popt0[0], vary=False)
                params.add('a1', value=popt0[1], vary=False)
                params.add('G1', value=popt0[2], min=Gmin          , max=Gmax       )
            if len(xind) > 1:
                Gmin = popt0[5]*0.333
                Gmax = get_Gmax(popt0[5]*3.,x_to_THz)
                params.add('w2', value=popt0[3], vary=False)
                params.add('a2', value=popt0[4], vary=False)
                params.add('G2', value=popt0[5], min=Gmin          , max=Gmax       )
            if len(xind) > 2:
                Gmin = popt0[8]*0.333
                Gmax = get_Gmax(popt0[8]*3.,x_to_THz)
                params.add('w3', value=popt0[6], vary=False)
                params.add('a3', value=popt0[7], vary=False)
                params.add('G3', value=popt0[8], min=Gmin          , max=Gmax       )
            if len(xind) > 3:
                Gmin = popt0[11]*0.333
                Gmax = get_Gmax(popt0[11]*3.,x_to_THz)
                params.add('w4', value=popt0[9], vary=False)
                params.add('a4', value=popt0[10],vary=False)
                params.add('G4', value=popt0[11], min=Gmin           , max=Gmax      )
            params = fix_params_long_and_t1_for_t2(params,space_fft_which,path,qpoint,filename,mdsteps,x_to_THz)
            if args.verbose > 2:
                print('############## vor INIT_D fit ######################## >>>>>>>>')
                my_print_params(params,x_to_THz,rt=4)
                print('############## vor INIT_D fit ######################## <<<<<<<<')
            result = gmodel.fit(y, params, w=x)
            popt0 = lmfit_result_to_popt(result,xind)


            outorig_y = function(x, *popt0)
            xout,yout = ps_to_ps_for_writing(outorig_y,xind,x_to_THz,function,x,popt0,args=args)
            if args.verbose > 1:
                print("**** ----> INIT_D   ",np.round(popt0*x_to_THz,4))

            if write_for_debugging:
                print('writing INIT_D_GROB ')#,np.round(popt0*x_to_THz,3)
                #parnpsavetxt("INIT_D_GROB",np.transpose(np.array([xout,yout*psmax])),fmt='%.8f %.15f')
                parnpsavetxt("INIT_D_GROB",np.transpose(np.array([x*x_to_THz,outorig_y*psmax])),fmt='%.8f %.15f')
            #sys.exit('hebe hier auf')

            ratio = y.max()/yout.max()
            if ratio <= 0.12 or ratio >= 4.5:

                popt0 = my_little_minimizer(x,y,popt0,function,xind,x_to_THz,sigma=0.,verbose=args.verbose)
                outorig_y = function(x, *popt0)
                xout,yout = ps_to_ps_for_writing(outorig_y,xind,x_to_THz,function,x,popt0,args=args)

            ratio = y.max()/yout.max()
            if ratio <= 0.12 or ratio >= 4.5:
                if mdsteps > 3000:
                    print(printred('222: ratio: '+str(ratio)+' mdsteps: '+str(mdsteps)+"  len(xind)"+str(len(xind))))
                #sys.exit('bad fit!')
                return np.zeros(int(len(xind)*3)),0,function #,pcov0


        if False: # only for debugging
            print('xmax',x.max())
            print('xout.max',xout.max())
            import matplotlib.pyplot as plt
            plt.plot(x*x_to_THz, y, 'b.-')
            plt.plot(x*x_to_THz, result_GROB.best_fit, 'r-')
            plt.plot(xout, yout, 'k--')
            plt.xlim([7.,10.])
            #plt.ion()
            plt.show()
            sys.exit()


    if args.verbose > 1:
        print("**** ----> INIT_EXT ",np.round(popt0*x_to_THz,4))



    ###########################################################################
    # INIT_B_MIDDLE  (minimize one peak after the other iterative)
    # too many peaks at once the minimizer can not handle
    # therefore:
    #   - if one peak   -> probably nothing to do but to use now the sigma
    #   - if two peaks  -> fit iteratively one peak and then the other
    #   - if four peaks -> fit iteratively one peak and then the other
    #                   -> currently only for one case, we should not bother
    ###########################################################################

    if True:
        # IS ALREADY VERY CLOSE TO DB
        #print 'popt0',popt0
        #print 'xind',xind
        #print 'sigma',sigma
        #print("x:",x,len(y),type(x))
        #print("y:",y,len(y),type(y))
        #np.savetxt("xout.dat",x)
        #np.savetxt("yout.dat",y)
        # curr curr curr
        exp_fitber = False
        if exp_fitber:
            cut_x_min = 3.5 # THz
            cut_x_max = 7   # THz

            cxmin = int(cut_x_min/x_to_THz)
            cxmax = int(cut_x_max/x_to_THz)
            x = x[cxmin:cxmax]
            y = y[cxmin:cxmax]
        #print("x:",x,len(y),type(x))
        #print("y:",y,len(y),type(y))
        #print("cxmin",cxmin)
        #print("cxmax",cxmax)
        #sys.exit()
        pin_stable = my_little_minimizer(x,y,popt0,function,xind,x_to_THz,sigma=0.,verbose=args.verbose)
        #sys.exit()
        if sigma == 0:
            ######################################################################
            # for sigma == 0 it turnes out htat my_little_minimizer converges to
            # the same result for large number of mdsteps
            # it converges however much quicker with the number of mdsteps
            ######################################################################
            if args.verbose > 1:
                print("**** ----> INIT_OUT ",np.round(pin_stable*x_to_THz,4))
            if True:
                print(printred('#################### very stable sigma = 0. ###############'))
                print(mdsteps,'pin_stable for weight',pin_stable*x_to_THz)
                print(printred('#################### very stable sigma = 0. ###############'))
                #sys.exit()
            return np.abs(pin_stable),0,function

        # wenn mann weightsfix benutzt um fuer sigma = 1.0 zu fitten geht dies schneller als mit
        # mit (function(x, *pin))**sigma ... konvergiert aber zum gleichen
        # dies ist auch ein gutes verfahren ...ist auch ok fuer 1BZ (space fft phase)
        # schlaegt aber fehl fuer z.b. 12 12 0
        weightsfix = (function(x, *pin_stable))**sigma

        #if args.verbose:
        if True:
            print(printred('#################### very stable sigma = 0. ###############'))
            print(mdsteps,'pin_stable for weight',pin_stable*x_to_THz)
            print(printred('#################### very stable sigma = 0. ###############'))
        #sys.exit('ka1234334')

    if False:
        pin_quest = my_little_minimizer(x,y,pin_stable,function,xind,x_to_THz,sigma=1.,verbose=args.verbose)
        if True:
            print(printred('####################             sigma = 1. ###############'))
            print(mdsteps,'pin_quest            ',pin_quest*x_to_THz)
            print(printred('####################             sigma = 1. ###############'))





    # in case no convergence is reached, the grob fit is taken!
    pin_start = np.copy(popt0)  # this is done for every iteration
    imax = 0
    for i in np.arange(100):  # sometimmes 30 are necessary, for small number of steps 70 might be necessary
        kwfit=[x_to_THz,xind,psmax,args.verbose,sigma]
        pinkeep = np.copy(popt0)
        if i == 0:
            pin = np.copy(popt0)
        else:
            pin = np.copy(pin_inner)

        if args.verbose > 1:
            print("**** ----> INIT_FS",i,np.round(pin*x_to_THz,4))
            if write_for_debugging:
                print('**** writing INIT_X_WD:')
                out = function(x, *pin)
                #parnpsavetxt("INIT_DB_EST",np.transpose(np.array([x*x_to_THz,y*psmax])),fmt='%.8f %.5f')
                parnpsavetxt("INIT_X_WD"+str(i),np.transpose(np.array([x*x_to_THz,out*psmax])),fmt='%.8f %.5f')


        pin_start_i = np.copy(pin)  # this is done for every iteration
        popt0sigmaa = np.copy(popt0)
        len_xind_2_normal = False
        if len(xind) == 2: len_xind_2_normal = True
        if len(xind) == 2 and check_qpoint_for_crossing and qlen >= 0.75 and qlen <= 0.85:
            len_xind_2_normal = False

        if len(xind) == 1:         # by definition there is no crossing
            for varyone in np.arange(9):
                w_unc = 0.15 # 15 %
                ####################################################################################
                # Here it turn out not to be good to vary 3 variables at once, vary one at a time
                ####################################################################################
                try:
                      params
                except NameError:
                    params = gmodel.make_params()

                params.add('w1', value=pin[0], vary=False)
                params.add('a1', value=pin[1], vary=False)
                params.add('G1', value=pin[2], vary=False)
                if varyone == 0:
                    params.add('w1', value=pin[0], min=pin[0]*w_unc, max=pin[0]*(1.+w_unc))
                if varyone == 1:
                    params.add('a1', value=pin[1], min=pin[1]*0.666, max=pin[1]*1.5)
                if varyone == 2:
                    params.add('G1', value=pin[2], min=pin[2]*0.666, max=pin[2]*1.5)
                # in very few cases sigma can get too small for low q vecs
                if popt0sigmaa[2]*x_to_THz < 0.4: popt0sigmaa[2] = 1.0/x_to_THz
                imax += 1
                if args.verbose > 2:
                    print('############## zuvor ######################## >>>>>>>>')
                    my_print_params(params,x_to_THz,rt=4)
                    print('############## zuvor ######################## <<<<<<<<')

                if True:
                    #print 'sigma',sigma
                    weights = (function(x, *pin))**sigma
                    #print 'pin weight',pin*x_to_THz
                    result = gmodel.fit(y, params, w=x, weights=weights,method='leastsq') #'nelder')
                    pin = lmfit_result_to_popt(result,xind)
                    #print 'pin out',pin*x_to_THz


                #@ lmfit funktioniert einfach nitcht gut!! def residual(params, x, y):
                #@ lmfit funktioniert einfach nitcht gut!!     #print 'in--------'
                #@ lmfit funktioniert einfach nitcht gut!!     #my_print_params(params,1,rt=4)
                #@ lmfit funktioniert einfach nitcht gut!!     #print 'in--------'
                #@ lmfit funktioniert einfach nitcht gut!!     w1 = params['w1']
                #@ lmfit funktioniert einfach nitcht gut!!     a1 = params['a1']
                #@ lmfit funktioniert einfach nitcht gut!!     G1 = params['G1']
                #@ lmfit funktioniert einfach nitcht gut!!     model = lorenz_68_Hauer2015_1p_1p(x,w1,a1,G1)
                #@ lmfit funktioniert einfach nitcht gut!!     return np.abs((y-model))
                #@ lmfit funktioniert einfach nitcht gut!! from lmfit import minimize,report_fit
                #@ lmfit funktioniert einfach nitcht gut!! out = minimize(residual, params, args=(x, y))
                #@ lmfit funktioniert einfach nitcht gut!! #print 'out',out
                #@ lmfit funktioniert einfach nitcht gut!! print report_fit(out)



                #yf = function(x, *pin)
                ##np.savetxt('xy',np.transpose((x,y)))
                ##np.savetxt('xyf',np.transpose((x,yf)))
                #abssum  = np.sum(np.abs((y-yf)))    # THIS SHOULD BE THE ACTUAL WEIGHT
                #abssum2 = np.sum(np.abs((y-yf))*weights)
                #leastsq = np.sum((y-yf)*(y-yf))
                #sys.exit('1231233')
                getout,check0 = check_if_parameters_converged(pin_start_i,pin,dp)
                #print i,varyone,'pin',pin*x_to_THz,getout,check0 #,'weights',weights*x_to_THz
                pin_inner = np.copy(pin)
                if args.verbose > 2:
                    #print 'args.verbose',args.verbose
                    print('############## danach ######################## >>>>>>>>')
                    #print(result.fit_report())
                    print('pin inner',pin*x_to_THz)
                    #print "--->>>",printred(str(abssum)),printred(str(abssum2))
                    #print "--->>>",printblue(str(leastsq))
                    print('############## danach ######################## <<<<<<<<')

            #sys.exit('12345ee')
        if len(xind) > 1:
            pin_inner = np.copy(pin)
            for i_inner in np.arange(10):
                pin_inner_start_i = np.copy(pin_inner)
                #print
                #print "**** ----> INIT_ES",i,i_inner,np.round(pin_inner*x_to_THz,3)
                imax += 1
                pin_inner = fit_several_peaks_by_holding_fix_all_but_one_peak(pin_inner,function,x,y,pin=pin_inner,kwfit=kwfit,esto=esto,mdsteps=mdsteps,verycertainestimate=verycertainestimate,space_fft_which=space_fft_which,path=path)
                #print "**** ----> INIT_B ",i,i_inner,pin_inner
                getout_inner,check0 = check_if_parameters_converged(pin_inner_start_i,pin_inner,dp)
                if args.verbose > 1:
                    print("**** ----> INIT_E ",i,i_inner,np.round(pin_inner*x_to_THz,3),getout_inner)
                if getout_inner == True:
                    pin = pin_inner
                    if args.verbose > 1:
                        print("**** ----> INIT_EF",i,i_inner,np.round(pin*x_to_THz,3),getout_inner)
                    break
                #sys.exit('6yt5')
        #print ' **pin_start_i',pin_start_i
        #print ' ** ppp pin',pin
        getout,check0 = check_if_parameters_converged(pin_start_i,pin,dp)

        if getout== True:
            break
    if args.verbose > 1: # in case no convergence is reached, the grob fit is taken!
        print("**** ----> INIT_F ",i,np.round(pin*x_to_THz,3),getout)

    if write_for_debugging:
        outorig_y = function(x, *pin)
        xout,yout = ps_to_ps_for_writing(outorig_y,xind,x_to_THz,function,x,pin,args=args)
        print('writing INIT_F      ')#,np.round(popt0*x_to_THz,3)
        parnpsavetxt("INIT_F",np.transpose(np.array([x*x_to_THz,outorig_y*psmax])),fmt='%.8f %.15f')
        parnpsavetxt("INIT_FF",np.transpose(np.array([xout,yout*psmax])),fmt='%.8f %.15f')

    #sys.exit('6yt5oof')
    ############################################################################
    # and a final triel (with weighting according to sigma)
    # with all variables let to minimize then compare std/or chisquare
    ############################################################################
    #print 'pin final',pin*x_to_THz
    #pin_qqq = my_little_minimizer(x,y,pin,function,xind,x_to_THz,sigma,args.verbose)
    if False:
        obrained_pin           = np.copy(pin)
        obrained_result        = np.copy(result)
        obrained_result_chisqr = np.copy(result.chisqr)
        print('result 1.',obrained_result_chisqr)
        if len(xind) > 0:
            params.add('w1', value=pin[0], min=pin[0]*0.95, max=pin[0]*1.05)
            params.add('a1', value=pin[1], min=pin[1]*0.95, max=pin[1]*1.05)
            params.add('G1', value=pin[2], min=pin[2]*0.95, max=pin[2]*1.05)
        if len(xind) > 1:
            params.add('w2', value=pin[3], min=pin[3]*0.95, max=pin[3]*1.05)
            params.add('a2', value=pin[4], min=pin[4]*0.95, max=pin[4]*1.05)
            params.add('G2', value=pin[5], min=pin[5]*0.95, max=pin[5]*1.05)
        if len(xind) > 2:
            params.add('w3', value=pin[6], min=pin[6]*0.95, max=pin[6]*1.05)
            params.add('a3', value=pin[7], min=pin[7]*0.95, max=pin[7]*1.05)
            params.add('G3', value=pin[8], min=pin[8]*0.95, max=pin[8]*1.05)
        if len(xind) > 3:
            params.add('w4', value=pin[9], min=pin[9]*0.95, max=pin[9]*1.05)
            params.add('a4', value=pin[10],min=pin[10]*0.95, max=pin[10]*1.05)
            params.add('G4', value=pin[11],min=pin[11]*0.95, max=pin[11]*1.05)

        weights = (function(x, *pin))**sigma
        result = gmodel.fit(y, params, w=x, weights=weights)
        pin = lmfit_result_to_popt(result,xind)
        check_result_chisqr = np.copy(result.chisqr)
        print('result 2.',check_result_chisqr)
        if args.verbose > 1:
            print("**** ----> INIT_G ",i,np.round(pin*x_to_THz,3),getout)
        if write_for_debugging:
            outorig_y = function(x, *pin)
            xout,yout = ps_to_ps_for_writing(outorig_y,xind,x_to_THz,function,x,pin,args=args)
            print('writing INIT_G      ')#,np.round(popt0*x_to_THz,3)
            parnpsavetxt("INIT_G",np.transpose(np.array([xout,yout*psmax])),fmt='%.8f %.15f')

    if getout == True: # INIT_END
        return np.abs(pin),imax,function #,np.abs(pcov0)
    else:
        #return np.zeros(len(dp)),imax #,pcov0 # global name 'dp' is not defined
        return np.zeros(int(len(xind)*3)),imax,function #,pcov0

def number_of_atoms_to_struct(atoms,args):
    fcc = []
    bcc = []
    for i in np.arange(100)+1:
        #print i,"fcc:",4*(i**3),"\tbcc:",2*(i**3)
        fcc.append(4*(i**3))
        bcc.append(2*(i**3))
    #print fcc
    #if args.verbose > 1:
    #    print("atoms",atoms)
    for i,k in enumerate(np.arange(100)+1):
        #print "fcc["+str(i)+"] in bcc:",fcc[i] in bcc
        pass
    if atoms in fcc:
        sc = N = int(round(float((atoms/float(4))**(1/3.)),0))
        return "fcc",sc
    if atoms in bcc:
        sc = N = int(round(float((atoms/float(2))**(1/3.)),0))
        return "bcc",sc

def get_inputvariables_from_calculation(positionsfilename = "positions.*",infilefilename = "in_file_dynamics.in", inN=False,inalat=False,indt=False,inqvec=False,verbose=False,args=False):
    ''' get the lammps input variables
        gets:
            - atoms (number of atoms e.g. 32)
            - sc = N (supercellsize e.g. 2)
            - alat (in Angstrom e.g. 4.14)
            - supercelllength (in Angstrom e.g. 8.28)
            - dt (timestep of dumping of positions)

    '''
    positionsfile = glob.glob(positionsfilename)
    if len(positionsfile) != 1:
        #print "ERROR: positionsfile:",positionsfile
        return False

    infile = glob.glob(infilefilename)  # in_file_dynamics.in
    #print "infile:",infile
    if len(infile) != 1:
        print(("ERROR: infile:",infile))
        return False

    #if not args.notverbose:
    print("JOB             : LAMMPS")
    atoms1 = os.popen('wc -l '+positionsfile[0] + "| awk '{print $1}'").read()
    atoms = int(atoms1)-8
    atoms2 = int(os.popen('head -n 2 '+positionsfile[0] + "| tail -1 | awk '{print $1}'").read().rstrip())
    #print("atoms1:",atoms1)
    #print("atoms2:",atoms2)
    if atoms != atoms2:
        print("positionsfile:",positionsfile[0])
        print("ERROR: atoms: (wc -l positionsfile -8)",atoms,"atoms2 (head -n 2 positionsfile):",atoms2,"from",positionsfile[0])
        return False


    supercelllength = float(os.popen('head -n 4 '+positionsfile[0] + "| tail -1 | awk '{print $2}'").read().rstrip())
    if verbose:
        print("supercelllength :",supercelllength)
    #sc = N = int(round(float((atoms/float(args.usestruct))**(1/3.)),0))
    args.structure, args.supercell = number_of_atoms_to_struct(atoms,args)
    sc = N = args.supercell
    alat = a = supercelllength/N
    if verbose:
        print("alat            :",alat)

    # check if structure and atoms make sense
    if type(args.usestruct) != bool:
        atoms3 = args.usestruct * sc**3.
        if atoms != atoms3:
            print("args:",args)
            print("usestruct:",args.usestruct, "4=fcc, 2=bcc use option: -s {bcc,fcc}")
            print("ERROR: atoms:",atoms,"atoms3:",atoms3, "maybe you have bcc instad of fcc or viec versa? define this!")
            return False

    stepslammps = int(os.popen('grep "^run " '+infile[0] + "| awk '{print $2}'").read().rstrip())
    timestepfs = int(os.popen('grep "^timestep " '+infile[0] + "| awk '{print $2*1000}'").read().rstrip())
    dt = dump = os.popen('grep "^dump dump1" '+infile[0] + "| awk '{print $5}'").read().rstrip()
    try:
        dt = dump = int(os.popen('grep "^dump dump1" '+infile[0] + "| awk '{print $5}'").read().rstrip())*timestepfs
    except ValueError:
        dt = False

    if verbose:
        print("timestepfs      :",timestepfs)
        print("atoms           :",atoms,type(atoms))
        print("sc              :",sc)
        print("stepslammps     :",stepslammps)
        print("dump            :",dump)
    stepswritten = stepslammps/dt+1
    if verbose:
        print("stepswritten    :",stepswritten)
    lineswritten = stepswritten*atoms
    if verbose:
        print("lineswritten    :",lineswritten)

    #########################################
    # check against invariables
    #########################################
    if type(args.folder) == bool:  # if False --> check
        if type(inN) != bool:
            if inN != N:
                print("args.folder:",args.folder)
                print("inN:",inN)
                print("N  :",N)
                sys.exit("ERROR: inN not equal N!")
    if type(inalat) != bool:
        if inalat != alat:
            print("(inalat      ):",inalat)
            print("(alat        ):",alat)
            print("(os.getcwd() ):",os.getcwd())
            if args.folder == False:
                sys.exit("ERROR: inalat not equal alat!")
            else:
                print("now changing from inalat:",inalat,"to alat  :",alat)
                args.alat = alat
                inalat = alat
    if type(indt) != bool:
        if indt != dt:
            print("indt:",indt)
            print("dt  :",dt)
            sys.exit("ERROR: indt not equal dt!")

    #########################################
    # check against invariables
    #########################################
    qvec = inqvec
    if type(inqvec) == bool:
        qvec = [ 'all', 'all', 'all' ]


    if args.supercell == False: args.supercell = N
    elif args.supercell != N: sys.exit("args.supercell, N: "+str(args.supercell)+" "+str(N))


    if args.dt == False: args.dt = dt
    elif args.dt != dt: sys.exit("args.dt, dt: "+str(args.dt)+" "+str(dt))

    if args.alat == False: args.alat = alat
    elif args.alat != alat: sys.exit("args.alat, alat: "+str(args.alat)+" "+str(alat))

    if args.qvec == False: args.qvec = qvec
    return

def get_inputvariables_from_c_skript_md_long_tox_job():
    if os.path.isfile("out_info.txt"):
        print("JOB             : C SKRIPT MD_LONG_TOX.c")
    else:
        return


    sc = int(float(os.popen("grep \"^N              :\" out_info.txt| awk '{print $3}'").read().rstrip()))
    alat = float(os.popen("grep \"^alat_lattice   :\" out_info.txt| awk '{print $3}'").read().rstrip())
    dt = int(float(os.popen("grep \"^dt             :\" out_info.txt| awk '{print $3}'").read().rstrip())*10**15)
    l = int(float(os.popen("grep \"^l              :\" out_info.txt| awk '{print $3}'").read().rstrip()))
    #print('sc:',sc)
    #print('alat:',alat)
    if l != 1:
        print('dt (grep fromskript):',dt)
        print('l  (grep fromskript):',l)
    args.supercell = sc
    args.dt = (dt*l)/10**12
    args.alat = alat
    return

def get_inputvariables_from_vasp_job():
    vaspinput = False
    if os.path.isfile("INCAR") and os.path.isfile("POSCAR"): # and os.path.isfile("POSITIONs"):
        vaspinput = True
    if vaspinput:
        print("JOB             : VASP")
    else:
        print("JOB             : not VASP")
        return
    try:
        atoms = int(float(os.popen("head -n 7 POSCAR | tail -1 | xargs -n 1 | awk '{sum+=$1} END {print sum}'").read().rstrip()))
        #print('atoms111',atoms)
        if atoms == 0:
            atoms = int(float(os.popen("head -n 6 POSCAR | tail -1 | xargs -n 1 | awk '{sum+=$1} END {print sum}'").read().rstrip()))
    except ValueError or atoms == 0:
        atoms = int(float(os.popen("head -n 6 POSCAR | tail -1 | xargs -n 1 | awk '{sum+=$1} END {print sum}'").read().rstrip()))
        #print('atoms222',atoms)

    dt = int(float(os.popen("grep POTIM INCAR | awk '{print $3}'").read().rstrip()))
    scale = float(os.popen("head -n 2 POSCAR | tail -1 | awk '{print $1}'").read().rstrip())
    sc11 = float(os.popen("head -n 3 POSCAR | tail -1 | awk '{print $1}'").read().rstrip())
    #print "atoms:",atoms
    #print "scale:",scale
    #print "sc11:",sc11
    args.dt = dt
    args.alat = scale*sc11

    #print('atoms',atoms)
    args.structure, args.supercell = number_of_atoms_to_struct(atoms,args)
    if args.structure == "fcc": args.usestruct = 4;
    if args.structure == "bcc": args.usestruct = 2;


    args.alat = scale*sc11/args.supercell
    return

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def help(p = None):
    string = '''
    - use this skript without options if in a lammps job (in_file_dynamics.in and positions*.dat available!)
      otherwise define required options manually!

    - required input is: (if no means of parsing from in_file_dynamics.in and positions*)
                    - N = sc = n
                    - dt
                    - a = alat
                    - q = qvec

    - examples:
        phonon_lifetimes.py (use without options if in a lammps job)
        phonon_lifetimes.py -N 10 -a 4.13 -dt 40 -q l 0 0
        phonon_lifetimes.py -N 3 -a 4.14 -dt 1 -q lt 0 0 -lt
        phonon_lifetimes.py -N 3 -a 4.14 -dt 1 -q l l l -ps -fftpy
        phonon_lifetimes.py -sc 10 -dt 40 -a 4.14 -q all 0 0
        phonon_lifetimes.py -s bcc -n 6 -dt 50 -a 3.253 -fftpy -pm -ps --space_fft_phase False (for dominiques bcc and omega job's)
        phonon_lifetimes.py -sc 3 -a 4.237 -atoms 216 -p 2 -h 1 -fftpy (TiN 3x3x3sc)

    - how to convert trj_lammps.out -> trj_lammpsnew.out:
             grep "^1 " trj_lammps.out | awk \'{print $2,$3,$4}\' > trj_lammpsnew.out

    '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter) #ArgumentDefaultsHelpFormatter)

    p.add_argument('-s',   '--structure',choices=[ 'fcc', 'bcc' ],
       help='which structure was calculated? Currently only fcc and bcc are supported. Simple cubic should work too.', type=str, default='fcc')
    p.add_argument("--space_fft_exp", type=str2bool, nargs='?',const=True, default=True,help="evaluate phonon frequencies/lifetimes by standard equation? default = True")
    p.add_argument("--space_fft_phase", type=str2bool, nargs='?',const=True, default=True,help="evaluate phonon frequencies/lifetimes in first BZ (space_fft_phase)? default = True")
    p.add_argument("--sigma0", type=str2bool, nargs='?',const=True, default=True,help="evaluate phonon frequencies/lifetimes by fitting to unweighted lorenzian (sigma = 0)? default = True")
    p.add_argument("--sigma1", type=str2bool, nargs='?',const=True, default=True,help="evaluate phonon frequencies/lifetimes by fitting to unweighted lorenzian (sigma = 0)? default = True")
    p.add_argument('-f',   '--folder',
       help='which folder to evaluate; expansion symbols are allowed but enclose in quotes e.g. "run_*" or "run_[1-5]"', type=str, default=False)
    p.add_argument('-u',   '--usestruct',help=argparse.SUPPRESS, type=int, default=4) # usestruct is 4 for fcc (4 atoms in cubic cell) and 2 for bcc (2 atoms in cubic cell)
    p.add_argument('-lalodiff',   '--separation_acoustic_optical',help='minimal distance in Thz of optical and acoustic branches (this helps to find the peaks in the power spectrum)' , type=float, default=3) # usestruct is 4 for fcc (4 atoms in cubic cell) and 2 for bcc (2 atoms in cubic cell)
    #p.add_argument('-at','--atoms', help="number of atoms",type=int, default=99) # usestruct is 4 for fcc (4 atoms in cubic cell) and 2 for bcc (2 atoms in cubic cell)
    p.add_argument('-n', '-sc','--supercell' , required=False,
       help='supercell: 2x2x2 -> 2', type=int, default=False)
    p.add_argument('-nm', '-scm','--supercellmultiply' , required=False,
       help='multiply atomic positions of the supercell: make 4x4x4 from 2x2x2 supercell', type=int, default=False)
    p.add_argument('-p', '--peaks' , required=False,
       help='amount of peaks to search in the powerspektrum', type=int, default=1)
    #p.add_argument('-at', '--atomms' , required=False,
    #   help='number of atoms', type=int, default=False)
    p.add_argument('--N', required=False,
       help=argparse.SUPPRESS, type=int, default=False)
    p.add_argument('--db', required=False,
            help=argparse.SUPPRESS, type=str, default='sql.db')
    #p.add_argument('--db', required=False,
    #        help=argparse.SUPPRESS, type=str, default=":memory:")
    p.add_argument('-dt',   '--dt' , required=False,
       help='timestep of your calculation in femto seconds.', type=float, default=False)
    p.add_argument('-a',   '--alat' , required=False,
       help='alat of your qubic cell in angstromr; e.g. 4.13', type=float, default=False)
    #p.add_argument('-ka', '--ka',action='append',nargs=3,metavar=('q1','q2','q3'),required=False,type=str,help='help:')


    p.add_argument('-ll',   '--ll' , required=False, action='append', type=list,
            help=argparse.SUPPRESS, default=False)
    p.add_argument('-mdstepstocalc_all',   '--mdstepstocalc_all' , required=False, action='append', type=int,
            help=argparse.SUPPRESS, default=False)
    p.add_argument('-mdstepstocalc_all_ps',   '--mdstepstocalc_all_ps' , required=False, action='append', type=int,
            help=argparse.SUPPRESS, default=False)
    p.add_argument('-mdstepstocalc_last','--mdstepstocalc_last', action='store_true', default=False,
        help='calculate powerspectrum only for all (maximum number of) steps and not intermediate step values to estimate the MD error;')
    p.add_argument('--mdsteps','-mdsteps' , required=False, type=int, nargs='+',default=False,
       help='calculate powerspectrum for particular mdsteplength; e.g. -mdsteps 30000 60000')
    p.add_argument('-mdstepstocalc_all_show',   '--mdstepstocalc_all_show' , required=False, action='store_true', default=True, help=argparse.SUPPRESS)


    p.add_argument('-q',   '--qvec' , required=False, action='append',nargs=3, type=str,
            help=textwrap.dedent('''\
       q1 q2 q3;
       Use 3 0 0      (three integers) to specivy a particular q-vector;
       use l 0 0      speify all in Longitudinal [ i 0 0 ] direction;
       use t 0 0      speify all in Transversal  [ 2N i 0 ] direction;
       use lt 0 0     speify all in Long and Trans   [ 2N i 0 ] direction;

       use l l l      speify all in Longitudinal [ i i i ] direction;
       use t t t      speify all in Transversal  [ ? ? ? ] direction;

       use l l 0      speify all in Longitudinal [ i i 0 ] direction;
       use t1 t1 0    speify all in Transversal1  [ i i 2N ] direction;
       use t2 t2 0    speify all in Transversal1  [ 2N+i 2N-i 0 ] direction;

       use all 0 0    if one of qvec is "all" then all branches are evaluated in Longitudinal and Transversal direction;
       use all l l    if one of qvec is "all" then all branches are evaluated in Longitudinal and Transversal direction;
       '''))
    p.add_argument('-qgar', '-garqp','--add_garching_qp','-add_gar_ap','-qg', action='store_true', default=False,
        help='add qpoints for higher brillouin zones; default = False')


    p.add_argument('-sql',   '--sql' , required=False, action='append',nargs='*', type=str,
            help='show the entries of the sql.tmp database and exit')

    p.add_argument('-pd','--pos_direct', action='store_true', default=False,
        help='set this flag if your input POSITIONs are in direct coords; this is only used in -fftpy')
    p.add_argument('-pm','--pos_modulo', action='store_true', default=False,
        help='set this flag if your input POSITIONs are are not mapped to original cell (go outside of supercell);this is only used if -fftpy')
    p.add_argument('-cl','--create_lammps_inputfile', action='store_true', default=False,
        help='create in_file_dynamics.in including defined qpoints')
    p.add_argument('-space_fft','--space_fft', action='store_true', default=False,
        help='get space_fft_x_x_x.npy from xaa_, xab_, ... folder')

    p.add_argument('-lt','--lifetimes', action='store_true', default=False,
            help='print lifetimes in the specified q-vecotor direction (and write to file); (units: THz );')
    p.add_argument('-wfd','--write_for_debugging', action='store_true', default=False,
            help='write INIT_{A,B,C,D,E,F} files')
    p.add_argument('-wv','--write_var_files', action='store_true', default=False,
            help='write var_{w,a,G} which show fits to the lifetimes as fct. of mdsteps')


    p.add_argument('-freqs','--freqs', action='store_true', default=False,
            help='print frequencies in the specified q-vecotor direction (and write to file); (units: THz);')

    #        help='average several lifetimes/frequencies files); (units: THz );')
    p.add_argument('--lifetimesaverage','-la','-ltsum', '-psav' , required=False, type=str, nargs='*',default=False,
       help=textwrap.dedent('''\
               average several lifetimes/frequencies files use
               --lifetimesaverage "subfolder*/ps_summary/" _lorenz_68_Hauer2015_xa_sigma_1__1__
               or
               --lifetimesaverage "subfolder*/ps_summary/" _lorenz_68_Hauer2015_xa_sigma_1__2__
               or
               --lifetimesaverage "subfolder*/ps_summary/" _good__1__
               or
               --lifetimesaverage "subfolder*/ps_summary/" _good__2__
               BUT remember to put the path in "" !
               '''))

    p.add_argument('-showeq','--showeq', action='store_true', default=False,
        help='print equivalent qpoint to screen;')
    p.add_argument('-showlammps','--showlammpsinput', action='store_true', default=False,
        help='print input for lammps inputfile for this paritcular qpoint(s) to screen;')

    p.add_argument('-fftpy','--make_space_fft_from_py', action='store_true', default=False,
        help='create space_fft_x_x_x.npy file using pandas read csv;')
    p.add_argument('-fftlammpslog','--make_space_fft_from_lammpslog', action='store_true', default=False,
        help='create space_fft_x_x_x.npy lammps.log file using numpy genfromtxt;')
    p.add_argument('-fftc','--make_space_fft_from_c', action='store_true', default=False,
        help='create space_fft_x_x_x.npy file using c++ skript;')
    p.add_argument('-ps','--make_power_spectrum', action='store_true', default=False,
        help='create ps_x_x_x_MDSTEPS.dat file and evaluate error due to MD/smoothing/background;')
    p.add_argument('-write_full_ps','-write_ps_full','-write_ps','--write_full_ps','--write_ps_full', action='store_true', default=False,
        help='Write the full powerspectrum to hd. Files can be relatively larg, 15-150 MB (default = False)')
    p.add_argument('-write_fitted_ps','-write_ps_fitted','--write_fitted_ps','--write_ps_fitted', action='store_true', default=False,
        help='Write the fitted powerspectrum to hd. Files can be relatively larg (default = False)')
    p.add_argument('-ser', '-seriell','--seriell','-serial', action='store_true', default=False,
        help='calculate everything seriell instead of parallell(==delfault); default = False')
    p.add_argument('-write_smooth_in', '-write_si','--write_smooth_in','-si' , action='store_true', default=False,
        help='Write the smoothed powerspectrum which is used to find the linewidths.')
    p.add_argument('-i',   '--lammps_infile',
       help='name of the lammps inputfile. can also be read from file infile.infilefilename', type=str, default='in_file_dynamcs.in')


    p.add_argument('-smoothing',   '--smoothing' , required=False,
       help='wite to hd a particular smoothing ; e.g. -mdsteps 30000 -smoothing 0.001', type=float, default=False)
    p.add_argument('-e',   '--exit' ,  action='store_true', default=False,
       help='exit after collecting job information')
    p.add_argument('-ee',   '--exclude_existing' ,  action='store_true', default=False,
       help='exclude existing space_fft_xxx files from qpoints to be used')
    p.add_argument('-v','--verbose',
            help='verbose', action='count', default=False)
    p.add_argument('-nv','-quiet','--notverbose',
            help='dont write anything but errors/warnings', action='count', default=False)
    p.add_argument('-atoms','--atoms', required=False,
       help="number of atoms in the system", type=int, default=False)

    return p
p = help()  # this gives the possibility to change some __init__ settings
args = p.parse_args()

##########################################################################################
# argparse first done
##########################################################################################
from scipy.optimize import curve_fit
from scipy import optimize
from scipy.integrate import simps

def importt(args):
    ''' import statements
    # timing: 3.73059201241 for pandas
    # timing: 1.33269906044 for scipy.fftpack
    '''
    global interp1d
    from scipy.interpolate import interp1d
    global Model
    global minimize
    from lmfit import Model,minimize
    global warnings
    global time
    import time
    anfang= time.time()
    start = time.time()
    global fft
    global ifft
    from scipy.fftpack import fft,ifft
    end = time.time()
    if args.verbose > 1:
        print("scipy.fftpack     imported in",str((end-start)),"sec.")
    start = time.time()
    global sys
    import sys
    global sqlite3
    import sqlite3
    end = time.time()
    if args.verbose > 1:
        print("sys               imported in",str((end-start)),"sec.")
    start = time.time()
    global os
    import os
    end = time.time()
    if args.verbose > 1:
        print("os                imported in",str((end-start)),"sec.")
    start = time.time()
    global math
    import math
    end = time.time()
    if args.verbose > 1:
        print("math              imported in",str((end-start)),"sec.")
    start = time.time()
    global glob
    import glob
    end = time.time()
    if args.verbose > 1:
        print("glob              imported in",str((end-start)),"sec.")
    start = time.time()
    global np
    import numpy as np
    end = time.time()
    if args.verbose > 1:
        print("numpy             imported in",str((end-start)),"sec.")
    start = time.time()
    global permutations
    from itertools import permutations
    global multiprocessing
    import multiprocessing
    end = time.time()
    if args.verbose > 1:
        print("itertools         imported in",str((end-start)),"sec.")
    start = time.time()
    #from pandas import set_option
    end = time.time()
    if args.verbose > 1:
        print("pandas set_option imported in",str((end-start)),"sec.")
    if args.make_space_fft_from_lammpslog or args.make_space_fft_from_py:
        start = time.time()
        global read_csv
        global set_option
        from pandas import read_csv,set_option
        set_option('max_columns', 0)
        end = time.time()
        if args.verbose > 1:
            print("pandas read csv   imported in",str((end-start)),"sec.")
    if args.verbose > 1:
        print("-----------------------------------------------------")
    #import pyfftw
    #import h5py
    #import timeit

    end = time.time()
    # OPTIONS
    if args.verbose > 1:
        print("all python module imported in",str((end-anfang)),"sec.")
        print("-----------------------------------------------------")
    return

importt(args)
np.set_printoptions(suppress=True)   # display arrays withou 000000
np.set_printoptions(precision=6,linewidth=150)

def lifetimes_freqs_summary_read_save(save=False,lt_lorz=False,fq_lorz=False,lt_good=False,fq_good=False):
    file_fq_lorz = "freqs_lorz.dat"
    file_lt_lorz = "lifet_lorz.dat"
    file_fq_good = "freqs_good.dat"
    file_lt_good = "lifet_good.dat"

    if save == True:
        fmt='%.0f %.0f %.0f    %.3f %.3f %.3f %.3f %.0f'
        if fq_lorz.shape[0] > 0:
            #print "fq_lorz-->",fq_lorz.shape,fq_lorz.shape[0]
            #print "fq_lorz.shape",fq_lorz.shape
            #print "fq_lorz",fq_lorz
            #print
            #print "fq... -->"
            #print fq_lorz[:,:8]
            np.savetxt(file_fq_lorz,fq_lorz,fmt=fmt)
        if lt_lorz.shape[0] > 0:
            np.savetxt(file_lt_lorz,lt_lorz,fmt=fmt)
        if fq_good.shape[0] > 0:
            print("fq_good",fq_good)
            np.savetxt(file_fq_good,fq_good,fmt=fmt)
        if lt_good.shape[0] > 0:
            print("lt_good",lt_good)
            np.savetxt(file_lt_good,lt_good,fmt=fmt)
        return

    if os.path.exists(file_fq_lorz):
        fq_lorz = np.loadtxt(file_fq_lorz)
    else:
        #fq_lorz = np.zeros((1,6))
        #fq_lorz = []
        #fq_lorz = np.array([])
        fq_lorz = np.empty((0,8), float)

    if os.path.exists(file_lt_lorz):
        lt_lorz = np.loadtxt(file_lt_lorz)
    else:
        #lt_lorz = np.zeros((1,6))
        lt_lorz = np.empty((0,8), float)

    if os.path.exists(file_fq_good):
        fq_good = np.loadtxt(file_fq_good)
    else:
        #fq_good = np.zeros((1,6))
        fq_good = np.empty((0,8), float)

    if os.path.exists(file_lt_good):
        lt_good = np.loadtxt(file_lt_good)
    else:
        #lt_good = np.zeros((1,6))
        lt_good = np.empty((0,8), float)

    return fq_lorz,lt_lorz,fq_good,lt_good

def lifetimes_freqs_summary_get_q(qp1,qp2,qp3,fq_lorz,lt_lorz,fq_good,lt_good):
    qp1=int(qp1)
    qp2=int(qp2)
    qp3=int(qp3)
    # print "--> q1,q2,q3:",qp1,qp2,qp3
    #print "--> q1,q2,q3:",qp1,qp2,qp3,"||||",fq_lorz.shape
    #print "qp1:",qp1,type(qp1)
    #print "qp2:",qp1,type(qp2)
    #print "qp3:",qp1,type(qp3)
    #print "fq_lorz:",fq_lorz,"||",fq_lorz[0]
    out1, out2, out3, out4 = False, False, False, False
    #for zeile in fq_lorz:
    #    print "-->> zeile:",zeile
    #    for i in zeile:
    #        print i
    #print "x:",fq_lorz,len(fq_lorz)
    #print "y:",len(fq_lorz.shape)
    for i in fq_lorz:
        #print "i1;:",i
        if i[0] == qp1 and i[1] == qp2 and i[2] == qp3:
            out1 = i
            #if len(fq_lorz.shape) == 1:
            #    out1 = [i]
    for i in lt_lorz:
        if i[0] == qp1 and i[1] == qp2 and i[2] == qp3:
            out2 = i
            #if len(lt_lorz.shape) == 1:
            #    out2 = [i]
    for i in fq_good:
        if i[0] == qp1 and i[1] == qp2 and i[2] == qp3:
            out3 = i
            #if len(fq_good.shape) == 1:
            #    out3 = [i]
    for i in lt_good:
        if i[0] == qp1 and i[1] == qp2 and i[2] == qp3:
            out4 = i
            #if len(lt_good.shape) == 1:
            #    out4 = [i]

    #print "--------9999999:",out2 # out2 kann auch bool sein (False)
    #print "--------88888888",out2.shape
    #print out1,out2,out3,out4
    return out1,out2,out3,out4

def lifetimes_freqs_summary_change(data_saved,line_new):
    qp1=int(line_new[0])
    qp2=int(line_new[1])
    qp3=int(line_new[2])
    #print "qp--->:",qp1,qp2,qp3,type(qp1)
    line_found = False
    line_idx = False
    #print "data-saved:",data_saved
    for idx,line in enumerate(data_saved):
        #print "line:",line
        #print int(line[0]),qp1
        #print int(line[1]),qp2
        #print int(line[2]),qp3
        if int(line[0]) == qp1 and int(line[1]) == qp2 and int(line[2]) == qp3:
            line_found = True
            line_idx = idx
            #print "??",line_found,idx,line_idx


    #print "line_found?",line_found,line_idx
    if line_found == True:
        #print "CHANGE 1"
        data_saved[line_idx] = line_new

    if line_found == False:
        #print "append line (0)",data_saved,"||",line_new,"==",data_saved.shape[0]
        if data_saved.shape[0] == 0:  # just add the line
            #data_saved = np.hstack((data_saved, np.array([line_new])))
            #print "APPEND first time"
            data_saved = np.append(data_saved, np.array([line_new]), axis=0)

        #data_saved.append(line_new)  # works for lists: []
        #data_saved = np.vstack((data_saved, np.array([4,5,6])))
        #data_saved = np.vstack((data_saved, line_new))
        else:  # change particular line
            #print "APPEND"
            data_saved = np.vstack((data_saved, line_new))
        #print "append line (1)",data_saved,"||",line_new

    return data_saved

##########################################################################################
# related to fitting peaks
##########################################################################################
def indexes(y, thres=0.3, min_dist=1):
    """Peak detection routine.

    Finds the numeric index of the peaks in *y* by taking its first order difference. By using
    *thres* and *min_dist* parameters, it is possible to reduce the number of
    detected peaks. *y* must be signed.

    Parameters
    ----------
    y : ndarray (signed)
        1D amplitude data to search for peaks.
    thres : float between [0., 1.]
        Normalized threshold. Only the peaks with amplitude higher than the
        threshold will be detected.
    min_dist : int
        Minimum distance between each detected peak. The peak with the highest
        amplitude is preferred to satisfy this constraint.

    Returns
    -------
    ndarray
        Array containing the numeric indexes of the peaks that were detected
    """
    if isinstance(y, np.ndarray) and np.issubdtype(y.dtype, np.unsignedinteger):
        raise ValueError("y must be signed")

    thres = thres * (np.max(y) - np.min(y)) + np.min(y)
    min_dist = int(min_dist)

    # compute first order difference
    dy = np.diff(y)

    # propagate left and right values successively to fill all plateau pixels (0-value)
    zeros,=np.where(dy == 0)

    # check if the singal is totally flat
    if len(zeros) == len(y) - 1:
        return np.array([])

    while len(zeros):
        # add pixels 2 by 2 to propagate left and right value onto the zero-value pixel
        zerosr = np.hstack([dy[1:], 0.])
        zerosl = np.hstack([0., dy[:-1]])

        # replace 0 with right value if non zero
        dy[zeros]=zerosr[zeros]
        zeros,=np.where(dy == 0)

        # replace 0 with left value if non zero
        dy[zeros]=zerosl[zeros]
        zeros,=np.where(dy == 0)

    # find the peaks by using the first order difference
    peaks = np.where((np.hstack([dy, 0.]) < 0.)
                     & (np.hstack([0., dy]) > 0.)
                     & (y > thres))[0]

    # handle multiple peaks, respecting the minimum distance
    if peaks.size > 1 and min_dist > 1:
        highest = peaks[np.argsort(y[peaks])][::-1]
        rem = np.ones(y.size, dtype=bool)
        rem[peaks] = False

        for peak in highest:
            if not rem[peak]:
                sl = slice(max(0, peak - min_dist), peak + min_dist + 1)
                rem[sl] = True
                rem[peak] = False

        peaks = np.arange(y.size)[~rem]

    return peaks

def get_idx_of_peaks(cb,peaks_amount):
    ''' this still can be made a lot faster by first checking the maximal and minimal thres
        before checking all 10 '''
    orders = 8
    orders = 6
    orders = 4
    verbose = False
    nn=30

    def get_thres_min_thres_max_or_result(thres_in,peaks_amount,loop=False,verbose=False,nn=10):
        if verbose:
            print()
            print()
            print("thres_in:",thres_in)
        out = np.zeros((2,len(thres_in)))
        loopbreak = False
        for idx,thres in enumerate(thres_in):
            peak_idx = indexes(cb, thres=thres, min_dist=nn) # 0.001 finds more peaks than 0.1
            if len(peak_idx) < peaks_amount and loop == False:
                if verbose:
                    print("###################!!!!!!!!!!!!!!!!!! thres:",idx,thres,len(peak_idx),peak_idx)
                out[0,idx] = thres
                out[1,idx] = round(len(peak_idx),3)
                loopbreak = idx
                break
            else:
                if verbose:
                    print(idx,"thres:",thres,len(peak_idx)) # ,peak_idx
            out[0,idx] = thres
            out[1,idx] = round(len(peak_idx),3)
        if verbose:
            print("out[0]:",out[0])
            print("out[1]:",out[1])
            print("loopbreak:",loopbreak,type(loopbreak))
        if type(loopbreak) is not bool and loop == False:
            minidx = (np.where(out[0] == 0.)[0]).min()
            #print minidx
            out = out[:,:minidx]
            #print "out[0]:",out[0]
            #print "out[1]:",out[1]
        out_idx       = np.where(out[1] >= peaks_amount)[0]
        out_thres     = out[0,out_idx]
        out_num_peaks = out[1,out_idx]
        idx_min = out_idx[-1]
        #print "out_idx  :",out_idx[-1]
        #print "out_thres:",out_thres,out_thres[-1]
        #print "out_num_p:",out_num_peaks,out_num_peaks[-1]
        if verbose:
            print("idx_min:",idx_min)
        thres_min = out[0,idx_min]
        thres_max = out[0,idx_min+1]
        #if loop == False:
        #    thres_max = out[0,-1]+out[0,0]
        #    p_2 = 0
        #else:
        #    thres_max = out[0,idx_min+1]
        #    p_2 = out[1,idx_min+1]
        if loop:
            thres_step_next = thres_min    # in loop which goes through orders
        else:
            thres_step_next = out[0,0]/10. # in loop which diveds range by 10
        p_1, p_2 = out[1,idx_min], out[1,idx_min+1]
        if verbose:
            print("---> p_1, p_2:",p_1, p_2,"---> thres_min, thres_max:",thres_min, thres_max,"step_next:",thres_step_next)
        return p_1,p_2,thres_min,thres_max,thres_step_next


    thres_in = [10./10**i for i in np.arange(1,orders)[::-1]]
    p_1,p_2,thres_min,thres_max,step = get_thres_min_thres_max_or_result(thres_in,peaks_amount,loop=True,verbose=verbose,nn=nn)
    #print "k1  p_1, p_2:",p_1, p_2, "---> thres_min, thres_max:",thres_min, thres_max
    if p_1 == peaks_amount: return indexes(cb, thres=thres_min, min_dist=3)
    if p_2 == peaks_amount: return indexes(cb, thres=thres_max, min_dist=3)


    thres_in = np.arange(thres_min,thres_max+step,step)
    p_1,p_2,thres_min,thres_max,step = get_thres_min_thres_max_or_result(thres_in,peaks_amount,verbose=verbose,nn=nn)
    #print "k2  p_1, p_2:",p_1, p_2, "---> thres_min, thres_max:",thres_min, thres_max
    if p_1 == peaks_amount: return indexes(cb, thres=thres_min, min_dist=3)
    if p_2 == peaks_amount: return indexes(cb, thres=thres_max, min_dist=3)

    thres_in = np.arange(thres_min,thres_max+step,step)
    p_1,p_2,thres_min,thres_max,step = get_thres_min_thres_max_or_result(thres_in,peaks_amount,verbose=verbose,nn=nn)
    #print "k3  p_1, p_2:",p_1, p_2, "---> thres_min, thres_max:",thres_min, thres_max
    if p_1 == peaks_amount: return indexes(cb, thres=thres_min, min_dist=3)
    if p_2 == peaks_amount: return indexes(cb, thres=thres_max, min_dist=3)

    thres_in = np.arange(thres_min,thres_max+step,step)
    p_1,p_2,thres_min,thres_max,step = get_thres_min_thres_max_or_result(thres_in,peaks_amount,verbose=verbose,nn=nn)
    #print "k3  p_1, p_2:",p_1, p_2, "---> thres_min, thres_max:",thres_min, thres_max
    if p_1 == peaks_amount: return indexes(cb, thres=thres_min, min_dist=3)
    if p_2 == peaks_amount: return indexes(cb, thres=thres_max, min_dist=3)

    return False

##########################################################################################
# sqlite3 stuff to store lifetimes/frequencies
##########################################################################################
def create_db(db):
    if db != ":memory:":
        if os.path.exists(db):
            #print "datbase exists already!"
            return

    #if db != ":memory:":
    #if not os.path.exists(db):

    print("creating db!!")
    ##############################
    # create db
    ##############################

    conn = sqlite3.connect(db)

    c = conn.cursor()

    #c.execute('''CREATE TABLE lifetimes
    #             (id INTEGER PRIMARY KEY, ltfq text, q real, value real, a real, b real)''')
    print("making db lifetimes!")
    c.execute('''CREATE TABLE lifetimes
                 (id INTEGER PRIMARY KEY, ltfq text, q1 int,q2 int,q3 int, mdsteps int,peakidx int, valuefq real, valuelt real)''')
    print("making db mdstepsall!")
    c.execute('''CREATE TABLE mdstepsall
                 (id INTEGER PRIMARY KEY, lll text, mdsteps int)''')

    # Save (commit) the changes
    conn.commit()

    # We can also close the connection if we are done with it.
    # Just be sure any changes have been committed or they will be lost.
    conn.close()
    #print "created databse",db
    return

def insert_to_db_mdsteps_all(c,conn,mdsteps_all):
    for mdsteps in mdsteps_all:
        c.execute('''INSERT INTO mdstepsall(id, lll, mdsteps) VALUES(NULL,?,?)''', ("kkk",int(mdsteps)))
    conn.commit()
    return

def insert_to_db(c,conn,ltfq,q1,q2,q3,mdsteps,peakidx,valuefq,valuelt): #,errup,errdn,errmd):
    verbose = False
    # if not os.path.exists(db):
    #    create_db(db)
    # if verbose:
    #    print "db:",db
    #    print "ltfq:",ltfq
    #    print "q1:",q1
    #    print "value:",value
    #    print "errup:",errup
    #    print "errdn:",errdn

    #if os.path.exists(db):
    #print 'adding data'
    #conn = sqlite3.connect(db)
    #c = conn.cursor()

    #c.execute('''INSERT INTO lifetimes(ltfq, a, b) VALUES(?,?,?)''', ('fq',1.0,9.81))
    #c.execute("INSERT INTO lifetimes VALUES (NULL,'lt',0.1,0.5,0.6)")

    for i in np.arange(10):
        goout=True
        #print 'iiiiiiiii',i

        try:
            c.execute('''INSERT INTO lifetimes(id, ltfq, q1,q2,q3,mdsteps,peakidx, valuefq,valuelt) VALUES(NULL,?,?,?,?,?,?,?,?)''', (ltfq,q1,q2,q3,mdsteps,peakidx,valuefq,valuelt)) #,errup,errdn,errmd))
        except (sqlite3.OperationalError,sqlite3.DatabaseError) as error:
            time.sleep(1)
            if i > 7:
                print(('WARNING: sleep INSERT i:'+str(i)+" "+str(q1)+"_"+str(q2)+"_"+str(q3)+" mdsteps: "+str(mdsteps)))
            goout=False

        if goout == True:
            break




    #try:
    #    c.execute('''INSERT INTO lifetimes(id, ltfq, q1,q2,q3,mdsteps, value,a, b,c) VALUES(NULL,?,?,?,?,?,?,?,?,?)''', (ltfq,q1,q2,q3,mdsteps,value,errup,errdn,errmd))
    #except sqlite3.OperationalError:
    #    time.sleep(1)
    #    print_warning('sleep insert 1 '+str(q1)+"_"+str(q2)+"_"+str(q3)+"_"+str(mdsteps))
    #    c.execute('''INSERT INTO lifetimes(id, ltfq, q1,q2,q3,mdsteps, value,a, b,c) VALUES(NULL,?,?,?,?,?,?,?,?,?)''', (ltfq,q1,q2,q3,mdsteps,value,errup,errdn,errmd))

    # Save (commit) the changes
    conn.commit()
    #conn.close()
    return

def get_db_mdstepsall(db):
    out = []
    conn = sqlite3.connect(db)
    c = conn.cursor()
    # if nothing is written yet
    try:
        c.execute('SELECT * FROM mdstepsall ORDER BY id')
    except sqlite3.OperationalError:
        return False
    #print "ddd:",print_db_c_to_screen(c)
    #print "eee:",print_db_c_to_mdstepsall(c)

    for row in c.execute('SELECT * FROM mdstepsall ORDER BY id'):
        #print row,"\t||",row[0],row[1],"||"#,row.rowid
        out.append(row[2])
    conn.close()
    #print "out:",out
    #print "lne:",len(out)
    if len(out) < 5:
        return False
    else:
        return np.array(out)


def print_db_to_screen(db):
    conn = sqlite3.connect(db)
    c = conn.cursor()
    for row in c.execute('SELECT * FROM lifetimes ORDER BY id'):
        print(row) #,"\t||",row[0],row[1],"||"#,row.rowid
    conn.close()
    return

def print_db_c_to_screen(c):
    for row in c.execute('SELECT * FROM lifetimes ORDER BY id'):
        print(row,"\t||",row[0],row[1],"||")#,row.rowid
    return

def print_db_c_to_screen_mdstepsall(c):
    for row in c.execute('SELECT * FROM mdstepsall RDER BY id'):
        print(row,"\t||",row[0],row[1],"||")#,row.rowid
    return

def get_row_id(c,conn,ltfq='lt',q1=False,q2=False,q3=False,mdsteps=False,peakidx=False):
    try:
        myiterator = c.execute('SELECT * FROM lifetimes ORDER BY id')
    except sqlite3.DatabaseError:
        time.sleep(1)
        print_warning('sleep get_row_id 1 '+str(q1)+"_"+str(q2)+"_"+str(q3)+"_"+str(mdsteps))
        myiterator = c.execute('SELECT * FROM lifetimes ORDER BY id')

    #for row in c.execute('SELECT * FROM lifetimes ORDER BY id'):
    #print "yooo"
    for row in myiterator:
        #print 'rooowww',row,"\t||",row[0],row[1]#,"||"#,row.rowid
        if row[1] == ltfq:
            if row[2] == q1 and row[3]==q2 and row[4]==q3:
                if row[5]==mdsteps and row[6]==peakidx:
                    #conn.close()
                    conn.commit()
                    return row

    #conn.close()
    return False

def sql_db_get_data(c,conn,kw,q1,q2,q3,mdsteps):
    ''' (6, u'lt_good', 16, 0, 0, 10000, 1.593, 0.0, 0.0, 0.0) '''
    for row in c.execute('SELECT * FROM lifetimes ORDER BY id'):
        #print "rrr:",row[1],kw
        if row[1] == kw and row[2]==q1 and row[3]==q2 and row[4]==q3 and row[5]==mdsteps:
            return row
    return False

def sql_db_get_mdsteps_max(c,conn):
    ''' (6, u'lt_good', 16, 0, 0, 10000, 1.593, 0.0, 0.0, 0.0) '''
    mdstepsmax = 0
    for row in c.execute('SELECT * FROM lifetimes ORDER BY id'):
        #print "row[5]:",row[5]
        if row[5] > mdstepsmax:
            mdstepsmax = row[5]
    return mdstepsmax


def old_unused():
    print("ka")
    #def get_row_id_unopened_database(db,ltfq='lt',q1=False,q2=False,q3=False):
    #    misses mdsteps
    #    conn = sqlite3.connect(db)
    #    c = conn.cursor()
    #    for row in c.execute('SELECT * FROM lifetimes ORDER BY id'):
    #        #print row,"\t||",row[0],row[1],"||"#,row.rowid
    #        if row[1] == ltfq:
    #            if row[2] == q1 and row[3]==q2 and row[4]==q3:
    #                conn.close()
    #                return row
    #
    #    conn.close()
    #    return False

    #def change_or_add_value_to_db_which_is_unopened(db,ltfq,q1,q2,q3,mdsteps,value,errup,errdn,errmd):
    #   misses mdsteps
    #    if db != ":memory:":
    #        if not os.path.exists(db):
    #            create_db(db)
    #
    #    verbose = False
    #    if verbose:
    #        print "db:",db
    #        print "ltfq:",ltfq
    #        print "value:",value
    #        print "errup:",errup
    #        print "errdn:",errdn
    #        # make sure that your q has three digits
    #        print "q1:",q1
    #
    #    digits = 3
    #    #q = round(q,digits)
    #    q1 = int(q1)
    #    q2 = int(q2)
    #    q3 = int(q3)
    #    #print "v1:",value
    #    value = round(value,digits)
    #    #print "v2:",value
    #    #value = format(value, '.3f')
    #    #print "v3:",value
    #    errup = round(errup,digits)
    #    errdn = round(errdn,digits)
    #    errmd = round(errmd,digits)
    #
    #    if verbose:
    #        print "q1:",q1
    #        print 'update data, looking for ltfq=',ltfq,"and q1=",q1
    #    sys.exit('the next line wont work')
    #    row = get_row_id(db,ltfq=ltfq,q1=q1,q2=q2,q3=q3,mdsteps=mdsteps)
    #    #print "type:",type(row)
    #    if type(row) != bool:  # renew the line
    #        #print "found row"
    #        id=row[0]
    #        #print "row:",row,"id:",id
    #        conn = sqlite3.connect(db)
    #        c = conn.cursor()
    #        #print "UPDATE!! id:",id
    #        c.execute('''UPDATE lifetimes SET value=?,a = ?,b=?,c=? WHERE id = ? ''',
    #            (value,errup,errdn,errmd,id))
    #        ##c.execute('''UPDATE lifetimes SET a=99.8 WHERE ltfq="fq"''')
    #        #c.execute('''UPDATE lifetimes SET b = ? WHERE ltfq = ? ''',
    #        #         (88.4, 'fq'))
    #        conn.commit()
    #        conn.close()
    #    else:
    #        #print "did not find row!, adding it!"
    #        insert_to_db(db,ltfq,q1,q2,q3,mdsteps,value,errup,errdn,errmd)
    #    return
    return

def add_value_to_db_tmp(qp,mdsteps,peakidx,ltfq,fq,lt,aa):
    if ltfq == 0 and fq == 0 and aa == 0:
        pass
    else:
        global pppx_add_db
        pppx_add_db.append([qp[0],qp[1],qp[2],mdsteps,peakidx,ltfq,fq,lt,aa])
    return


def change_or_add_value_to_db(c,conn,qp,mdsteps,peakidx,ltfq,  valuefq,valuelt): #, errup,errdn,errmd):
    ''' this assumes that the db exists and is open and will be closed externally
        called by:
        change_or_add_value_to_db(c,conn,qpoint,'fq_good',*fq_good_1)
    '''
    mdsteps = int(mdsteps)
    verbose = False
    if verbose:
        print("ltfq:",ltfq)
        #print "value:",value
        #print "errup:",errup
        #print "errdn:",errdn
        print("mdsteps:",mdsteps)
        print("peakidx:",peakidx)
        # make sure that your q has three digits

    digits = 3
    #q = round(q,digits)
    q1 = int(qp[0])
    q2 = int(qp[1])
    q3 = int(qp[2])
    #print "v1:",value
    valuefq = round(valuefq,digits)   # lifetime or freq
    valuelt = round(valuelt,digits)   # lifetime or freq
    #print "v2:",value
    #value = format(value, '.3f')
    #print "v3:",value
    #errup = round(errup,digits)
    #errdn = round(errdn,digits)
    #errmd = round(errmd,digits)

    if verbose:
        print("q1:",q1)
        print("q2:",q2)
        print("q3:",q3)
        print('update data, looking for ltfq=',ltfq,"and q1=",q1)
    row = get_row_id(c,conn,ltfq=ltfq,q1=q1,q2=q2,q3=q3,mdsteps=mdsteps,peakidx=peakidx)  # return either False or rowname
    #print "type:",type(row),row
    #sys.exit()
    if type(row) != bool:  # renew the line
        #print "found row"
        id=row[0]
        #print "row:",row,"id:",id
        #conn = sqlite3.connect(db)
        #c = conn.cursor()
        #print "UPDATE!! id:",id
        for i in np.arange(10):
            goout=True
            #print 'iiiiiiiii',i

            try:
                c.execute('''UPDATE lifetimes SET valuefq=?, valuelt=?,a=?,b=?,c=? WHERE id = ? ''',(valuefq,valuelt,id)) #,errup,errdn,errmd,id))
            except (sqlite3.OperationalError,sqlite3.DatabaseError) as error:
                time.sleep(1)
                if i > 7:
                    print(('WARNING: sleep change_or_add i:'+str(i)+" "+str(q1)+"_"+str(q2)+"_"+str(q3)+" mdsteps: "+str(mdsteps)))
                goout=False

            if goout == True:
                break



        ##c.execute('''UPDATE lifetimes SET a=99.8 WHERE ltfq="fq"''')
        #c.execute('''UPDATE lifetimes SET b = ? WHERE ltfq = ? ''',
        #         (88.4, 'fq'))
        conn.commit()
        #conn.close()
    else:
        #print "did not find row!, adding it!"
        insert_to_db(c,conn,ltfq,q1,q2,q3,mdsteps,peakidx,valuefq,valuelt) #,errup,errdn,errmd)
    return


##########################################################################################
# general functions and helperfunctions ##################################################
##########################################################################################
def printblue(text):
    return ('\x1b[4;34;40m' + text + '\x1b[0m')
def printorange(text):
    return ('\x1b[4;33;40m' + text + '\x1b[0m')
def printgreen(text):
    return ('\x1b[4;32;40m' + text + '\x1b[0m')
def printred(text):
    return ('\x1b[4;31;40m' + text + '\x1b[0m')

def parnpsavetxt(filename,data,fmt=False):
    #np.savetxt(filenameout,ps[:len(ps)/factor_psout])  # saves only half the powerspektrum (left peak)
    if type(fmt) == bool:
        p = multiprocessing.Process(target=np.savetxt, args=[filename,data])
    else:
        p = multiprocessing.Process(target=np.savetxt, args=[filename,data,fmt])
    p.start()
    return

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def get_lifetimessummary_from_several_folders(base="run_*",lifetimestring="2.0"):
    if args.lifetimesaverage == []:
        pass
        #print "[]",args.lifetimesaverage
    else:
        base = args.lifetimesaverage[0]
        #print "[[[",args.lifetimesaverage
    alatstring = str(args.alat)
    lifetimestring = str(args.dt)
    np.set_printoptions(precision=4)
    np.set_printoptions(suppress=True)
    searchforall_dispdir = ["l_0_0", "t_0_0", "l_l_0", "t1_t1_0", "t2_t2_0" , "l_l_l", "t_t_t", 'lg_0_0','lg_lg_0','lg_lg_lg','tg_0_0','tg2_tg2_0' ]
    forfiles= [ "lifetimesgood_", "freqsgood_" ]
    forfiles=[args.lifetimesaverage[1]]
    for idxdispdir,dispdir in enumerate(searchforall_dispdir):
        for forefile in forfiles:
            for fqlt in [ 'fq','lt']:
                #print
                #print
                #print
                #filename = fqlt+'*'+forefile+"*"+dispdir+"_alat"+alatstring+"_dt"+lifetimestring+"_THz.dat"
                filename = fqlt+forefile+dispdir+"_alat"+alatstring+"_dt"+lifetimestring+"_THz.dat"
                searchstring = base+filename
                found=glob.glob(searchstring)
                #if len(found) == 0:
                if len(found) <1:
                    print("resulting base:",base)
                    print("     forfiles:",forfiles)
                    print("     dispdir:",dispdir)
                    print("resulting filename:",filename)
                    print("-->resulting searchstring:",searchstring)
                    print("")
                    print("--- this are the details ---")
                    print("forefile:",forefile)
                    print("dispdir:",dispdir)
                    print("alatstring:",alatstring)
                    print("--> lifetimestring:",lifetimestring)
                    for i in np.arange(len(searchstring)-1)+1:
                        searchstring= searchstring[:-i]+'*'
                        found=glob.glob(searchstring)
                        #print 'i',i,searchstring,"fournd?",len(found)
                    continue
                    #sys.exit("____________ no files found ___________________")
                if False:
                    if len(found) >=1:
                        for i in found:
                            print('found',i)
                        #print "found!!!!!!!!!!!!!",len(found),found
                if len(found) >=1:
                    #print "found!!!!!!!!!!!!!",len(found),found
                    if idxdispdir == 0:
                        for idx,i in enumerate(found):
                            print("$$",idx,i)
                    #print 'filename',filename
                    #sys.exit()
                    for idx,i in enumerate(found):
                        #data=np.loadtxt(i)
                        data=np.genfromtxt(i)
                        #print "idx:",idx,i,data
                        if idx == 0:
                            #print "idx:",idx,"data.shape:",data.shape,len(found),len(data.shape)
                            if len(data.shape) == 2:
                                datacumulated = np.zeros((len(found),data.shape[0],data.shape[1]))
                            if len(data.shape) == 1:
                                datacumulated = np.zeros((len(found),data.shape[0]))
                            #print "tt;",datacumulated.shape
                        print('i',i)
                        print('data',data.shape)
                        print(data)
                        datacumulated[idx]=data
                    #print "mean ++++++++++++++++++++++++++++",dispdir
                    #print "shape:",datacumulated.shape,len(datacumulated.shape)
                    mean = np.mean(datacumulated,axis=0)
                    max = np.max(datacumulated,axis=0)
                    min = np.min(datacumulated,axis=0)
                    #print "kk:",len(datacumulated.shape)
                    if len(datacumulated.shape) > 2: # more than one line
                      mean[:,5] = 0
                      max[:,5] = 0
                      min[:,5] = 0
                    else:  # only one line
                      mean[5] = 0
                      max[5] = 0
                      min[5] = 0
                    #print mean
                    #print "std  ++++++++++++++++++++++++++++",dispdir
                    std = np.std(datacumulated,axis=0)
                    #print std
                    #print "ste  ++++++++++++++++++++++++++++",dispdir,"################",len(found)
                    ste = np.std(datacumulated,axis=0)/np.sqrt(len(found))
                    #print ste


                    #print "ste  part +++++++++++++++++++++++",dispdir
                    #toget = np.std(datacumulated,axis=0)/np.sqrt(len(found))
                    #print "len(ste.shape):",len(ste.shape)
                    if len(ste.shape) == 2:
                        #print "------> mean[:,2],ste[:,1]:",mean[:,2],ste[:,1]
                        #mean[:,2] = mean[:,6]+ (ste[:,1]+std[:,1])/2.  # error up
                        #mean[:,3] = mean[:,7]+ (ste[:,1]+std[:,1])/2.  # error down
                        mean[:,2] = max[:,1]-mean[:,1] # error up
                        mean[:,3] = -(min[:,1]-mean[:,1])# error down
                    if len(ste.shape) == 1:
                        #print "------>   mean[2],ste[1]  :",mean[2],ste[1]
                        #mean[2] = mean[2]+ (ste[1]+std[1])/2.       # error up
                        #mean[3] = mean[3]+ (ste[1]+std[1])/2.       # error down
                        mean[2] = max[1]-mean[1] # error up
                        mean[3] = -(min[1]-mean[1]) # error down
                    #print "mean  again  +++++++++++++++++++++",searchfor
                    #print mean
                    if os.path.isdir('ps_summary_SUM/') != True:
                        os.makedirs('ps_summary_SUM')
                    np.savetxt('ps_summary_SUM/'+filename,mean,fmt="%.5f")
                    if len(ste.shape) == 1: # more than one line
                        #np.savetxt(filename,np.transpose(mean),fmt="%.5f")
                        #print "mean:",mean
                        np.savetxt('ps_summary_SUM/'+filename,np.array([mean]),fmt="%.5f")
                    print()
        pass
    print('written to ps_summary_SUM/')
    return

def qpointstring(q1,q2=False,q3=False):
    ''' creates string for the corresponding qpoint for files to write out'''
    #print "q1:",q1,type(q1)
    #print "q2:",q2,type(q2)
    #print "q3:",q3,type(q3)
    qpointstr=False
    if type(q2)==bool and type(q3)==bool:
        #print "1Xx"
        if type(q1)==np.ndarray or type(q1) == list:
            #print "2Xx"
            if len(q1) == 3:
                qpointstr = str(q1[0])+"_"+str(q1[1])+"_"+str(q1[2])
    else:
        qpointstr = str(q1)+"_"+str(q2)+"_"+str(q3)
    #print qpointstr
    #sys.exit()
    return qpointstr

def qpointstring_to_qpoint(qpointstring):
    results = qpointstring.split("_")
    #print 'res',results
    return [int(i) for i in results]


def get_all_qpoints(listeqpoints,args=False,verbose=False):
    '''
    [ 1 0 0 ] # particular qpoint
    [ 3 3 3 ] # particular qpoint

    [ l 0 0 ] # 100-L, get all long in this direction
    [ l l 0 ] # 110-L, get all long in this direction
    [ l l l ] # 111-L, get all long in this direction

    [ t 0 0 ] # 100-T, get all trans in this direction
    [ t t t ] # 111-T, get all trans in this direction
    [ t1 t1 0 ] # 110-T1, get all trans1 in this direction
    [ t2 t2 0 ] # 110-T2, get all trans2 in this direction
    '''
    N = args.supercell
    qpoints_all =[]
    appropratestring = [ 'l', 'l2', 'lgar', 't', 't1', 't2', 'all', 'lt', 'tnew', 'lg', 'tg', 'tg1','tg2']
    appropratestring = [ 'l', 'l2', 'lgar', 't', 't1', 't2', 'all', 'lt', 'tnew', 'lg', 'tg', 'tg1','tg2','gar']
    #print listeqpoints,type(listeqpoints)
    #verbose = False
    if verbose:
        print("listeqpoints",listeqpoints)
        #sys.exit()
    # DEFAULT QPOINTS
    if listeqpoints is None:
        listeqpoints = [['all','all','all']]
    if args.exit == True and listeqpoints is None:
        listeqpoints = [['all','all','all']]
    #for N1N2N3 in listeqpoints:
    #    print N1N2N3
    #sys.exit()
    if verbose:
        if not args.notverbose:
            print("qpoints_all  :     [ qpoint      ] qpstr     symeq     Nr. symeqtot     Reduced wave vector")
            print("-------------------------------------------------------------------------------------------")
    for N1N2N3 in listeqpoints:
        #print "iiiiiiiii:",N1N2N3
        N1 = N1N2N3[0]
        N2 = N1N2N3[1]
        N3 = N1N2N3[2]

        #print "type(N1):",N1,type(N1)
        #print "str(N1) :",N1,int(N1)
        #print "type(N2):",N2,type(N2)
        #print "str(N2) :",N2,int(N2)
        #print "type(N3):",N3,type(N3)
        liste = [N1,N2,N3]
        for idx,i in enumerate(liste):
            #print "i:",i,type(i)
            if type(i) != str and type(i) != int and type(i) != float:
                sys.exit("ERROR: qpoint needs to be a string, integer or float!")


            try:
                liste[idx]=int(i)
            except ValueError:
                if i not in appropratestring:
                    print("listeqpoints:",listeqpoints)
                    print("liste:",liste)
                    print("appropratestring:",appropratestring)
                    print("i    :",i)
                    sys.exit("ERROR: qpoint as string should one of \"l,t,t1,t2,l3\" but is \""+i+"\"")

        N1 = liste[0]
        N2 = liste[1]
        N3 = liste[2]

        if type(N1) != str and type(N2) != str and type(N3) != str:
            #qpoints_all = [[N1,N2,N3]] # we have an particular qpoint
            qpoints_all.append([N1,N2,N3])

        # in case however we had a string in there
        N = int(N)  # supercellsize
        all = False
        l100 = False
        t100 = False

        t100new = False
        l110 = False
        t1_110 = False
        t2_110 = False
        l111 = False
        t111 = False
        #gar = False
        #gar = True
        if N1 == 'all' or N2 == 'all' or N3 == 'all': all = True            # [ All ]
        #if N1 == 'gar' or N2 == 'gar' or N3 == 'gar': gar = True            # [ gar ]


        # [ 1 0 0 ]
        if (N1 == 'l' and N2 == 0. and N3 == 0.): l100 = True               # [ N 0 0 ] (L)
        if (N1 == 't' and N2 == 0. and N3 == 0.): t100 = True               # [ N 0 0 ] (T)
        if (N1 == 'tnew' and N2 == 0. and N3 == 0.): t100new = True               # [ N 0 0 ] (T)
        if (N1 == 'lt' and N2 == 0. and N3 == 0.): l100 = True;t100 = True  # [ N 0 0 ] (L+T)

        lg100 = False  # == 3BZ
        tg100 = False  # == 3BZ
        if (N1 == 'lg' and N2 == 0. and N3 == 0.): lg100 = True               # [ N 0 0 ] (L)
        if (N1 == 'tg' and N2 == 0. and N3 == 0.): tg100 = True               # [ N 0 0 ] (L)

        l2100 = False # == 2BZ
        t2100 = False # == 2BZ
        if (N1 == 'l2' and N2 == 0. and N3 == 0.): l2100 = True               # [ N 0 0 ] (L)
        if (N1 == 't2' and N2 == 0. and N3 == 0.): t2100 = True               # [ N 0 0 ] (L)

        # [ 1 1 0 ]
        if (N1 == 'l' and N2 == 'l' and N3 == 0.): l110 = True               # [ N N 0 ] (L)
        if (N1 == 't1' and N2 == 't1' and N3 == 0.): t1_110 = True               # [ N N 0 ] (L)
        if (N1 == 't2' and N2 == 't2' and N3 == 0.): t2_110 = True               # [ N N 0 ] (L)

        lg110 = False
        #tg2_110 = False  # same as the one taken t2_t2_0
        if (N1 == 'lg' and N2 == 'lg' and N3 == 0.): lg110 = True               # [ N N 0 ] (L)
        #if (N1 == 'tg2' and N2 == 'tg2' and N3 == 0.): tg2_110 = True               # [ N N 0 ] (L)

        # [ 1 1 1 ]
        if (N1 == 'l' and N2 == 'l' and N3 == 'l'): l111 = True               # [ N N N ] (L)
        if (N1 == 't' and N2 == 't' and N3 == 't'): t111 = True               # [ N N N ] (L)
        if (N1 == 'lt' and N2 == 'lt' and N3 == 'lt'): l111 = True;t111 = True # [ N N N ] (L)

        lg111 = False
        tg111 = False
        if (N1 == 'lg' and N2 == 'lg' and N3 == 'lg'): lg111 = True               # [ N N N ] (L)
        if (N1 == 'tg' and N2 == 'tg' and N3 == 'tg'): tg111 = True               # [ N N N ] (L)

        def check_if_to_add(args,qpoint,qpointsall,verbose=False):
            ''' add or not '''
            #print "args.exclude_existing:",args.exclude_existing
            if args.exclude_existing == False:
                #print "in False"
                if qpoint not in qpoints_all:
                    if verbose: print(qpoint,args.exclude_existing,"->",True,len(qpointsall))
                    return True
                else:
                    if verbose: print(qpoint,args.exclude_existing,"->",False,len(qpointsall))
                    return False
            elif args.exclude_existing == True:
                check1 = "space_fft_"+qpointstring(qpoint)+".npy"
                check = os.path.isfile(check1)
                # if file exists, dont add
                if check:
                    if verbose: print(qpoint,args.exclude_existing,"->",check,"->",False,len(qpointsall))
                    return False
                else:
                    if verbose: print(qpoint,args.exclude_existing,"->",check,"->",True,len(qpointsall))
                    return True
            else:
                sys.exit("Error 789")
            return


        symeq= []




        def qpoints_all_append(qp,add="",verbose=verbose):
            cr = check_qpoints_for_crossing(qp, args.supercell, args.structure)
            co=""
            if cr:
                co='\u2713'
                co="CROSSING"
            #print "cr:",cr
            qpoints_all.append(qp)
            i = qp  # i[0], i[1], i[2]
            ############ get length of qpoint
            #print "qp:",qp
            ja = 3
            symeq.append(equivalent_qs(i).shape[0])
            addstr = qpoint_get_l_0_0_or_l_l_0_or_l_l_l_or_t_0_0_or_t_t_t_or_t1_t1_t1_or_t2_t2_t2_from_qpoint(i,args)
            if verbose:
                if not args.notverbose:
                    #print "    "+add.ljust(13),"[",str(i[0]).ljust(ja),str(i[1]).ljust(ja),str(i[2]).ljust(ja),"]",qpointstring(i).ljust(10),"("+str(equivalent_qs(i).shape[0])+")".ljust(5),'\t',str(len(qpoints_all)).ljust(4),"("+str(np.sum(symeq))+")".ljust(3),'\t',"",str(qpoint_to_length_qpoint(q1=i[0],q2=i[1],q3=i[2],N=N,struct=args.structure)),'\t',co
                    print("    "+add.ljust(13),"[",str(i[0]).ljust(ja),str(i[1]).ljust(ja),str(i[2]).ljust(ja),"]",qpointstring(i).ljust(10),"("+str(equivalent_qs(i).shape[0])+")".ljust(4),'\t',str(len(qpoints_all)-1).ljust(4),"("+str(np.sum(symeq))+")".ljust(3),'\t',"",str(qpoint_to_length_qpoint(q1=i[0],q2=i[1],q3=i[2],args=args)).ljust(6),co,'\t',addstr)
            return
        #print "agg:",args.showeq
            #if args.showlammpsinput or args.create_lammps_inputfile: print equivalent_qs(i,show_lammps_inputfile=True,write_lammps_inputfile=args.create_lammps_inputfile,N=args.supercell,alat=args.alat,dt=args.dt)

        #print "all:",all
        #standard = True
        #if gar == True:
        #    standard = False
        #if args.add_garching_qp:
        #    gar = True
        #    standard = True
        if True:
            # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 0 0 ] >>>>>>>>>>>>>>>>>>>>>>
            if l100 or all == True:                                         # [ N 0 0 ] (L)
                for i in np.arange(N)+1:
                    qp = [i,0,0]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N 0 0 (L)")
                #if verbose: print
                if verbose and not args.notverbose: print()
            if t100 or all == True:                                         # [ N 0 0 ] (T)
                for i in np.arange(N)+1:
                    qp = [2*N,i,0]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N 0 0 (T)")
                #if verbose: print
                if verbose and not args.notverbose: print()
            #if t100new or all == True:                                         # [ N 0 0 ] (T)
            #    for i in np.arange(N)+1:
            #        qpoints_all.append([2*N,i,2*N])
            # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 0 0 ] >>>>>>>>>>>>>>>>>>>>>>

            # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 1 0 ] >>>>>>>>>>>>>>>>>>>>>>
            if args.structure == 'fcc': NN = N
            if args.structure == 'bcc': NN = N/2
            if l110 == True or all == True:      # [ N N 0 ] (L)
                for i in np.arange(NN)+1:
                    qp = [i,i,0]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N 0 (L)")
                if verbose and not args.notverbose: print()
            if t1_110 == True or all == True:    # [ N N 0 ] (T1)
                for i in np.arange(NN)+1:
                    qp = [2*N+i,2*N-i,0]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N 0 (T1)")
                if verbose and not args.notverbose: print()
            if t2_110 == True or all == True:    # [ N N 0 ] (T2)
                for i in np.arange(NN)+1:
                    qp = [i,i,2*N]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N 0 (T2)")
                if verbose and not args.notverbose: print()
            # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 1 0 ] >>>>>>>>>>>>>>>>>>>>>>


            # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 1 1 ] >>>>>>>>>>>>>>>>>>>>>>
            if l111 == True or all == True:     # [ N N N ] (L)
                if args.structure == 'fcc':
                    for i in np.arange(N/2)+1:
                        i = int(i)
                        qp = [i,i,i]
                        if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N N (L)")
                    qpadd_LO_TO = [N,N,N]
                    if check_if_to_add(args,qpadd_LO_TO,qpoints_all): qpoints_all_append(qpadd_LO_TO,"N N N (LO_TO)")
                    if verbose and not args.notverbose: print()

                if args.structure == 'bcc':
                    for i in np.arange(N)+1:
                        qp = [i,i,i]
                        if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N N (L)")
                    if verbose and not args.notverbose: print()
            if t111 == True or all == True:     # [ N N N ] (T)
                if args.structure == 'fcc':
                    #for i in (np.arange(N/2)+N/2)[::-1]:
                    #print 'aaa',(np.arange(N/2)+N/2),'NNN',N,(np.arange(N/2)+N/2)[::-1]
                    #print 'ka1',np.arange(N)
                    #print 'ka2',np.arange(N/2)
                    #print 'ka3',np.arange(N/2)+N/2
                    #print 'N',N,'ka4',np.arange(int(N/2.))+1
                    for i in np.arange(int(N/2.))+1:
                        #print "iii:",i,"NNN",N
                        qp = [i,i,2*N-i]
                        if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N N (T)")
                    #if verbose: print
                    if verbose and not args.notverbose: print()
                if args.structure == 'bcc':
                    for i in (np.arange(N)+1): #[::-1]:
                        qp = [i,i,2*N-i]
                        if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N N (T)")
                    #if verbose: print
                    if verbose and not args.notverbose: print()
            # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 1 1 ] >>>>>>>>>>>>>>>>>>>>>>

        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ gar ] >>>>>>>>>>>>>>>>>>>>>>


        #print "gar1:",gar,len(qpoints_all)
        #if gar == True:
        if True:
            if lg100 or all == True:                                         # [ N 0 0 ] (L)
                for i in np.arange(N)+1:                    # [ N 0 0 ] (L)
                    qp = [2*N+i,0,0]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N 0 0 (Lgar)")
                if verbose and not args.notverbose: print()
            if tg100 or all == True:                                         # [ N 0 0 ] (L)
                for i in np.arange(N)+1:                    # [ N 0 0 ] (T)
                    qp = [2*N,2*N,i]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N 0 0 (Tgar)")
                #if verbose: print
                if verbose and not args.notverbose: print()
            if l2100 or all == True:                                         # [ N 0 0 ] (L)
                for i in (np.arange(N+1)+N)[::-1]:                    # [ N 0 0 ] (L)
                    qp = [i,0,0]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N 0 0 (2 BZ)")
                if verbose and not args.notverbose: print()


            if lg110 == True or all == True:      # [ N N 0 ] (L)
                for i in np.arange(N)+1:                    # [ N N 0 ] (L)
                    qp = [2*N+i,2*N+i,0]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N 0 (Lgar)")
                #if verbose: print
                if verbose and not args.notverbose: print()
            #if tg2_110 == True or all == True:    # [ N N 0 ] (T2)
            #    for i in np.arange(N)+1:                    # [ N N 0 ] (T2)
            #        qp = [i,i,-2*N]
            #        if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N 0 (T2gar)")
            #    #if verbose: print
            #    if verbose and not args.notverbose: print
            #if tg1_110 == True or all == True:    # [ N N 0 ] (T2)
            #    for i in np.arange(N)+1:                    # [ N N 0 ] (T1)
            #        qp = [2*N+i,2*N-i,0]
            #        if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N 0 (T1gar)")
            #    #if verbose: print
            #    if verbose and not args.notverbose: print
            if lg111 == True or all == True:     # [ N N N ] (L)
                if args.structure == 'fcc':
                    for i in np.arange(N/2)+1:                  # [ N N N ] (L)
                        i = int(i)
                        qp = [N+i,N+i,N+i]
                        if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N N (Lgar)")
                    if verbose and not args.notverbose: print()
                if args.structure == 'bcc':
                    for i in np.arange(N)+1:                  # [ N N N ] (L)
                        qp = [N+i,N+i,N+i]
                        if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N N (Lgar)")

                    if verbose and not args.notverbose: print()

            if tg111 == True or all == True:     # [ N N N ] (L)
                if args.structure == 'fcc':
                    for i in np.arange(N/2)+1:                  # [ N N N ] (T)
                        i = int(i)
                        qp = [N-i,N-i,N+i]
                        if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N N (Tgar)")
                    if verbose and not args.notverbose: print()
                    #if verbose: print
                #if args.structure == 'bcc':
                #    for i in np.arange(N)+1:                  # [ N N N ] (T)
                #        qp = [N-i,N-i,N+i]
                #        if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N N (Tgar)")
        #print "gar2:",gar,len(qpoints_all)

        #if len(qpoints_all) == 0:
        #    sys.exit("ERROR: this path is not known!")
    if args.showeq:
        if verbose:
            print("-------------------------------------------------------------------------------------------")
        for idx,i in enumerate(qpoints_all):
            if args.showeq: print(equivalent_qs(i,show_without_ipi=True))
        if verbose:
            print("-------------------------------------------------------------------------------------------")
    return qpoints_all

def alltrue(items):
    print('items',items) #, end=' ')
    if all(item == True for item in items):
        return True
    else:
        return False

def qpoint_map_t2_qpoint_to_L_L_0_and_T1_T1_0_qpoint(qpoint,args):
    N = args.supercell
    if qpoint[0] == qpoint[1] and qpoint[0]<=N and qpoint[2] == 2*N:   # [ T2 T2 0 ]
        pass
        #return "t2_t2_0"
    else:
        sys.exit("THIS is not t2_t2_0")
    #qp = [2*N+i,2*N-i,0]
    i = qpoint[0]
    t10 = N*2 + qpoint[0]
    t20 = N*2 - qpoint[0]
    t30 = 0
    return np.array([i,i,0]),np.array([t10,t20,t30])



def qpoint_get_l_0_0_or_l_l_0_or_l_l_l_or_t_0_0_or_t_t_t_or_t1_t1_t1_or_t2_t2_t2_from_qpoint(qpoint,args):
    qpointlong = np.sort(qpoint)[::-1]
    if qpoint[1] == qpoint[2] == 0:                 # [ L 0 0 ]  1BZ und hoehere BZ
        return "l_0_0"
    if qpoint[0] == qpoint[1] and qpoint[2] == 0:   # [ L L 0 ]  1BZ und hoehere BZ
        return "l_l_0"
    if qpoint[0] == qpoint[1] == qpoint[2]:     # [ L L L ]  1BZ und hoehere BZ
        return "l_l_l"

    # it is apparently a transversal branch
    N = args.supercell
    if qpoint[0] == 2*N and qpoint[1]<=N and qpoint[2] == 0:   # [ T 0 0 ]
        return "t_0_0"
    if qpoint[0] == 2*N and qpoint[1]==2*N and qpoint[2] <= N:   # [ Tgar 0 0 ]
        return "t_0_0"
    if qpoint[0] > 2*N and qpoint[0] <=3*N and qpoint[1]<2*N and qpoint[0]+qpoint[1]==4*N and qpoint[2] == 0:   # [ T1 T1 0 ]
        return "t1_t1_0"
    if qpoint[1] == qpoint[2] and qpoint[1]<=N and qpoint[0] == 2*N:   # [ T2 T2 0 ]
        return "t2_t2_0"
    if qpoint[0] == qpoint[1] and qpoint[0]<=N and qpoint[2] == 2*N:   # [ T2 T2 0 ]
        return "t2_t2_0"
    if qpoint[0] == qpoint[1] and qpoint[0]<=N and qpoint[2] == -2*N:   # [ T2gar T2gar 0 ]
        return "t2_t2_0"
    if qpoint[0] == qpoint[1] and qpoint[0]<N and qpoint[2] > N:   # [ T T T ]
        return "t_t_t"
    return False

def qpoint_map_longitudinal_to_transversal(qpoint,args=False):
    if qpoint_is_longitudinal(qpoint,passhigherBZ=True,args=args):
        #print 'only map longitudinal qpoints'
        N = args.supercell
        qr = qpoint_get_l_0_0_or_l_l_0_or_l_l_l_or_t_0_0_or_t_t_t_or_t1_t1_t1_or_t2_t2_t2_from_qpoint(qpoint,args)
        if qr == 'l_0_0':
        #if qpoint[1] == qpoint[2] == 0:  # [ L 0 0 ]
            i = qpoint[0]
            qp = [2*N,i,0]
            return [qp]
        #if qpoint[0] == qpoint[1] and qpoint[2] == 0:  # [ L L 0 ]
        if qr == 'l_l_0':
            i = qpoint[0]
            qp1 = [2*N+i,2*N-i,0] # (T1)
            qp2 = [i,i,2*N]  # (T2)
            return [qp1,qp2]
        #if qpoint[0] == qpoint[1] == qpoint[2]:  # [ L L L ]
        if qr == 'l_l_l':
            qlen = qpoint_to_length_qpoint(qpoint[0],qpoint[1],qpoint[2],args,get_lll_double_length=False,verbose=False)
            #print 'qlen',qlen
            gothrough =  get_all_qpoints([['t','t','t']],args=args,verbose=False)
            #print 'gothrough',gothrough
            for i in gothrough:
                qlen2 = qpoint_to_length_qpoint(i[0],i[1],i[2],args,get_lll_double_length=False,verbose=False)
                #print 'qlen2',qlen2
                if qlen == qlen2:
                    #print('jjjjiii',[i])
                    return [i]
            if False: # old
                if args.structure == 'fcc':
                    gothrough = (np.arange(N/2)+N/2)[::-1]
                if args.structure == 'bcc':
                    gothrough = (np.arange(N)+1)
                print('gothrough',gothrough)
                for idx,i in enumerate(gothrough):
                    qp = [i,i,2*N-i]  # here we get all qpoint of this cell
                    print('qpoint-->',qpoint, 'idx',idx,'i',i,'all',(np.arange(N/2)+N/2)[::-1],'qp-->',qp,idx+1==qpoint[0])
                    if idx+1 == qpoint[0]:
                        print(('iiiiiii',[qp]))
                        return [qp]

            #i = qpoint[0]
            #qp = [i,i,2*N-i]
            #return [qp]
            #return np.array([0,0,0]) # in case of no clue

    else:
        #return 'qpoint is transversal xxyy'
        return False

def qpoint_is_longitudinal(qpoint,passhigherBZ=False,args=False):
    ishigherBZ = False
    if qpoint[0] > args.supercell:
        ishigherBZ = True
    if qpoint[1] > args.supercell:
        ishigherBZ = True
    if qpoint[2] > args.supercell:
        ishigherBZ = True
    #print 'ishigherBZ',ishigherBZ
    if passhigherBZ == False and ishigherBZ == True:
        return False
    if qpoint[1] == qpoint[2] == 0:  # [ L 0 0 ]
        return True
    if qpoint[0] == qpoint[1] and qpoint[2] == 0:  # [ L L 0 ]
        return True
    if qpoint[0] == qpoint[1] == qpoint[2]:  # [ L L L ]
        return True
    return False

def qpoint_to_length_qpoint(q1,q2,q3,args,get_lll_double_length=False,verbose=False):
    verbose = False
    N = args.supercell
    roundto = 4
    qo = qpoint_get_l_0_0_or_l_l_0_or_l_l_l_or_t_0_0_or_t_t_t_or_t1_t1_t1_or_t2_t2_t2_from_qpoint(np.array([q1,q2,q3]),args)
    # in a sc of 5x5x5 the qpoints at [400] and [600] are similar; also similar are [100] and [900]
    # so everything repeats ater sc*2 == 2*N
    # the lenght of qpoint [400] has the same length as [600]
    repeats_after = 2.

    if args.structure == 'fcc' and qo == 't_t_t':
        repeats_after = 1.
    if args.structure == 'fcc' and qo == 'l_l_l':
        repeats_after = 1.

    if args.structure == 'bcc' and qo == 'l_l_0':
        repeats_after = 1.
    if args.structure == 'bcc' and qo == 't1_t1_0':
        repeats_after = 1.
    if args.structure == 'bcc' and qo == 't2_t2_0':
        repeats_after = 1.

    ka0 = np.array([int(q1),int(q2),int(q3)])%(repeats_after*N)
    l = length = ka0.max()/float(N)



    if verbose:
        print('args.structure',args.structure)
        print('qo',qo)
        print("repeats_after",repeats_after)
        print("N",N)
        print("repeats_after*N",repeats_after*N)
        print("ka0",ka0)
        print('l',l)
        print('supercell:',N)
        print('--------ka> ',q1,q2,q3)
        print('ka  ',ka0)
        print('l1   ',l)

    if l > repeats_after/2.:
        l = length = repeats_after - l
    if verbose:
        print('l2   ',l)
    #ka1 = ka0[np.nonzero(ka0)[0]]/float(N)
    #if verbose:
    #    print 'ka1',ka1
    if verbose:
        if False:
            sys.exit('outt334456')
    #return round(ka1[0],roundto)
    #sys.exit('33')
    return round(l,roundto)
    ##@@ #ka2=ka%2
    ##@@ # ka%N makes higher BZ to 1st BZ this is necessary to get the right length in the 1BZ
    ##@@ ka2=ka%N
    ##@@ if verbose:
    ##@@     print 'ka2 ',ka2
    ##@@ # hmmm.. is ka22 really makes sense ..... for high symmetry directions this should be ok
    ##@@ # since any nonzero number would be the same
    ##@@ ka22 = ka2[np.nonzero(ka2)[0]]
    ##@@ if verbose:
    ##@@     print 'ka22',ka22

    ##@@ if verbose:
    ##@@     sys.exit('outt33445')
    ##@@ return round(ka22[0],roundto)

    ##@@ # ok, ka3 = ka%1 does not make sens to me
    ##@@ ka3=ka%1
    ##@@ if verbose:
    ##@@     print 'ka3',ka3
    ##@@ ka4 = ka3[np.nonzero(ka3)[0]]
    ##@@ if verbose:
    ##@@     print 'ka4',ka4
    ##@@     print 'len(ka4)',len(ka4)
    ##@@ if len(ka4) == 0:
    ##@@     qlen = 0
    ##@@ else:
    ##@@     qlen = ka22.min()  # this is however not working for ttt path
    ##@@ if verbose:
    ##@@     print 'qlen1',qlen
    ##@@ ############################################
    ##@@ # qlen correction for fcc [t t t]
    ##@@ ############################################
    ##@@ if q1==q2 and q3 != 0 and q3 != q1 and q3 != 2*N and args.structure=='fcc':
    ##@@     if len(ka4) > 0:
    ##@@         qlen = ka4.min()
    ##@@ if verbose:
    ##@@     print 'qlen2',qlen

    ##@@ ############################################
    ##@@ # for [l l l] and [t t t]
    ##@@ ############################################
    ##@@ if get_lll_double_length:
    ##@@     if q1==q2 and q3 != 0 and q3 != 2*N:
    ##@@         if round(qlen,roundto) > 1.:
    ##@@             qlen = qlen%1
    ##@@         return round(qlen*2,roundto)
    ##@@ #if in1==in2==in3=='t':
    ##@@ #    qlen = ka4.min()
    ##@@ #@    #elif struct=='bcc':
    ##@@ #@    #    print "kkk:",ka,ka%2
    ##@@ #print 'qlen3',qlen
    ##@@ if round(qlen,roundto) > 1.:
    ##@@     qlen = qlen%1
    ##@@ if verbose:
    ##@@     print 'out',round(qlen,roundto)
    ##@@ if verbose:
    ##@@     sys.exit('outt3344')
    ##@@ return round(qlen,roundto)

def get_fq_lt_md_error(lt_or_fq_series):
    idxmax = np.where(lt_or_fq_series == lt_or_fq_series.max())[0].max()

    lifetimestocheckerrormd = np.zeros(len(lt_or_fq_series))
    for ind,yy in enumerate(lt_or_fq_series):
        #print ind,y,"   ",lt_or_fq_series[:ind+1][-6:]
        considerlt = lt_or_fq_series[:ind+1][-5:]
        lifetimestocheckerrormd[ind] = considerlt.max() - considerlt.min()

    idxmax = np.where(lifetimestocheckerrormd == lifetimestocheckerrormd.max())[0].max()
    for ind,yy in enumerate(lifetimestocheckerrormd):
        if ind <= idxmax:
            lifetimestocheckerrormd[ind] = lifetimestocheckerrormd[idxmax]
    errorup   = lifetimestocheckerrormd

    def make_error_down_never_larger_as_value(series,error_md):
        out = np.zeros(len(series))
        if len(series) != len(error_md):
            return np.ones(len(series))*999.

        for idx,idy in enumerate(series):
            if series[idx] > error_md[idx]:
                out[idx] = error_md[idx]
            else:
                out[idx] = series[idx]
        return out

    errordown = make_error_down_never_larger_as_value(lt_or_fq_series,errorup)
    return errorup, errordown

def get_qpoint_mapping_index_dispdir(args,dispdir_all):
    qpoints_all_forplotting = np.array(get_all_qpoints([['all', 'all', 'gar']],args))  # dont change this, it is correct
    #for i in qpoints_all_forplotting:
    #    print 'iopuy',i
    #sys.exit()
    list_map_qpoint_to_id_dispdir = np.zeros((len(qpoints_all_forplotting),6))
    id=0
    for id_dispdir,dispdir in enumerate(dispdir_all):               # schleife ueber [ l 0 0 ], [ t 0 0], ...
        if args.verbose > 2:
            print('id_dispdir',id_dispdir,'out of len(dispdir_all)',len(dispdir_all)-1)
        in1, in2, in3 = dispdir[0],dispdir[1],dispdir[2]            # in1,in2,in3 = t1, t1, 0 oder l, 0, 0
        qpoints_dispdir = get_all_qpoints([[in1, in2, in3]],args)   # dont change this, it is correct
        for qp in qpoints_dispdir:                                  # qp = [1,0,0] oder [10,10,32] oder ...
            #print 'shapexx1',list_map_qpoint_to_id_dispdir.shape,'id',id,' --> qp:',qp[0],qp[1],qp[2],'<-- dispdir',dispdir,'id_dispdir:',id_dispdir, 'id',id #,'len(dispdir_all)',dispdir_all  #'len(qpoints_dispdir)',len(qpoints_dispdir)
            if id >= list_map_qpoint_to_id_dispdir.shape[0]:
                z = np.zeros((1, 6), dtype=list_map_qpoint_to_id_dispdir.dtype)
                list_map_qpoint_to_id_dispdir = np.concatenate((list_map_qpoint_to_id_dispdir,z), axis=0)
            list_map_qpoint_to_id_dispdir[id,0]=qp[0]
            list_map_qpoint_to_id_dispdir[id,1]=qp[1]
            list_map_qpoint_to_id_dispdir[id,2]=qp[2]
            list_map_qpoint_to_id_dispdir[id,3]=id_dispdir
            list_map_qpoint_to_id_dispdir[id,4]=id
            list_map_qpoint_to_id_dispdir[id,5]=0       # is the number of peaks, will later needs to be changed, THIS CAN JUST SAVE THE MAXIMUM NUMBER OF PEAKS WHICH SCHOULD BE ENOUGH
            id+=1
    return list_map_qpoint_to_id_dispdir

def get_qpoint_to_dispdir_idx(qpoint,list_map_qpoint_to_id_dispdir):
    verbose = False
    if verbose:
        print('qpoint:',qpoint)
    for idx,i in enumerate(list_map_qpoint_to_id_dispdir):
        if verbose:
            print('idx',idx,'/',len(list_map_qpoint_to_id_dispdir)-1,'i',i,"::",i[0],qpoint[0])
        if i[0]==qpoint[0] and i[1]==qpoint[1] and i[2]==qpoint[2]:
            if verbose:
                print('yo',i,"-->",int(i[3]))
            return int(i[3]),idx,int(i[5])
    return False

def get_dispdir_to_peaks_max():
    pass

def list_flatten(l, ltypes=(list, tuple)):
    ltype = type(l)
    l = list(l)
    #print "ltype:",ltype
    #print "l    :",l
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)


def get_fq_lt_summary_line(f,lastline,qp,printtoscreen=True):
    roundto=3
    ka = str(round(lastline[0],4)                        ).ljust(8)+"\t"+\
         str(round(lastline[1],roundto))+"\t"+\
         str(round(lastline[2],roundto))+"\t"+\
         str(round(lastline[3],roundto))+"\t"+\
         str(round(lastline[4],7))+"\t"+\
         str(999999999)+"___"+str(qp[0])+"_"+str(qp[1])+"_"+str(qp[2])+"\t"+\
         str(round(lastline[6],roundto))+"\t"+\
         str(round(lastline[7],roundto))+"\t"+\
         str(round(lastline[8],roundto))+"\t"+\
         "\n"
    ktoscreen = "\t["+str(qp[0])+" "+str(qp[1])+" "+str(qp[2])+"]\n"
    if printtoscreen:
        print(ka.rstrip()+ktoscreen.rstrip())
    f.write(ka)
    return

def model_lin(x,a,m):
    return m*x + a

def square_model(x,b):
    return x**b

def model_square_shifted(x,a,shift,max):
    #return a+ b*x +c*x**2
    return a*(x-shift)**2. + max

def model_constant(x,out):
    return out

def get_from_db_w_a_G(qpoint,func,peaks,mdsteps):
    '''
    ## sigma is saved in the func
    print n gives:
                                                                  freq   lt     aa
    [3, 0, 0, 6800, 1, 'lorenz_68_Hauer2015_xa_sigma_0_bounds_T', 7.802, 0.546, 0.1591]
    [3, 0, 0, 6800, 2, 'lorenz_68_Hauer2015_xa_sigma_0_bounds_T', 16.176, 0.426, 0.1928]
    [3, 0, 0, 6100, 1, 'lorenz_68_Hauer2015_xa_sigma_0_bounds_T', 7.713, 0.49, 0.1679]
    [3, 0, 0, 6100, 2, 'lorenz_68_Hauer2015_xa_sigma_0_bounds_T', 16.23, 0.427, 0.2082]
    '''
    if False:
        print('qpoint',qpoint)
        print('func',func)
        print('peaks',peaks,'---',list(range(peaks)))
        print('mdsteps',mdsteps)

    the_filename="sql.tmp"
    if not os.path.isfile(the_filename):
        return False
    out = np.zeros(peaks*3)

    with open(the_filename, 'rb') as f:
        try:
            all = pickle.load(f, encoding='latin1')            # geht ueber alle eintraege
        except TypeError:
            all = pickle.load(f)            # geht ueber alle eintraege

        for peakidx,peak in enumerate(np.arange(peaks)+1):
            #print 'peakidx,peak',peakidx,peak
            for idn,n in enumerate(all):
                #print 'xyp1',n
                if n[0] == qpoint[0] and n[1] == qpoint[1] and n[2] == qpoint[2] and n[3] == mdsteps and n[4] == peak and n[5] == func:
                    out[peakidx*3+0] = n[6]
                    out[peakidx*3+1] = n[8]
                    out[peakidx*3+2] = n[7]
                    break
                    #return np.array([n[6],n[8],n[7]])
                    #print 'xyp',n
                    #listout.append([n[3],n[4],n[5],n[6],n[7],n[8]])
                    #listout.append([n[3],n[4],n[6],n[7],n[8]]) # hier ist sigma 1 und 0 drin
    #print 'out1',out,"nonzero",np.count_nonzero(out),"len(out)",len(out)
    if np.count_nonzero(out) == 0 or np.count_nonzero(out) != len(out):
        out = False
    #print 'out2',out
    return out

def get_estimate_w_a_G(qpoint,func,for_mdsteps,verbose=False,args=False,add_filename='',x_to_THz=1.):
    '''
    print n gives:
                                                                  freq   lt     aa
    [3, 0, 0, 6800, 1, 'lorenz_68_Hauer2015_xa_sigma_0_bounds_T', 7.802, 0.546, 0.1591]
    [3, 0, 0, 6800, 2, 'lorenz_68_Hauer2015_xa_sigma_0_bounds_T', 16.176, 0.426, 0.1928]
    [3, 0, 0, 6100, 1, 'lorenz_68_Hauer2015_xa_sigma_0_bounds_T', 7.713, 0.49, 0.1679]
    [3, 0, 0, 6100, 2, 'lorenz_68_Hauer2015_xa_sigma_0_bounds_T', 16.23, 0.427, 0.2082]
    [3, 0, 0, 6000, 1, 'lorenz_68_Hauer2015_xa_sigma_0_bounds_T', 7.73, 0.616, 0.1677]
    [3, 0, 0, 6000, 2, 'lorenz_68_Hauer2015_xa_sigma_0_bounds_T', 16.127, 0.39, 0.2088]
    [3, 0, 0, 5400, 1, 'lorenz_68_Hauer2015_xa_sigma_0_bounds_T', 7.778, 0.565, 0.1985]
    [3, 0, 0, 5400, 2, 'lorenz_68_Hauer2015_xa_sigma_0_bounds_T', 15.926, 0.629, 0.2511]
    [3, 0, 0, 5000, 1, 'lorenz_68_Hauer2015_xa_sigma_0_bounds_T', 7.677, 0.656, 0.1974]
    [3, 0, 0, 5000, 2, 'lorenz_68_Hauer2015_xa_sigma_0_bounds_T', 16.2, 0.467, 0.2472]
    [3, 0, 0, 4700, 1, 'lorenz_68_Hauer2015_xa_sigma_0_bounds_T', 7.645, 0.457, 0.2123]
    [3, 0, 0, 4700, 2, 'lorenz_68_Hauer2015_xa_sigma_0_bounds_T', 16.186, 0.533, 0.2677]
    [3, 0, 0, 4000, 1, 'lorenz_68_Hauer2015_xa_sigma_0_bounds_T', 7.75, 0.726, 0.2562]
    [3, 0, 0, 4000, 2, 'lorenz_68_Hauer2015_xa_sigma_0_bounds_T', 16.266, 0.588, 0.309]
    for lorenz_68_Hauer2015 G follows the square model however
    for lorenz_86_Hauer2015 G follows rather the linear model
    func ist z.b. im sql.tmp: lorenz_68_Hauer2015_xa_sigma_1_bounds_T
    func ist z.b. im sql.tmp: lorenz_68_Hauer2015_xa_sigma_1
    '''
    #verbose = 5
    #print 'kk---------------------------------------------------------------------------------------',verbose

    #verbose = False
    the_filename="sql.tmp"
    if not os.path.isfile(the_filename):
        return False

    listout = []
    if verbose > 3:
        print('#db1# qpoint  :',qpoint)
        print('#db2# funct   :',func)

    with open(the_filename, 'rb') as f:
        try:
            all = pickle.load(f, encoding='latin1')            # geht ueber alle eintraege
        except TypeError:
            all = pickle.load(f)            # geht ueber alle eintraege
        for idn,n in enumerate(all):
            if verbose > 4:
                print('init!!!',n)
            #print 'xyp1',n
            if n[0] == qpoint[0] and n[1] == qpoint[1] and n[2] == qpoint[2] and n[5] == func:
                #print 'xyp',n
                #listout.append([n[3],n[4],n[5],n[6],n[7],n[8]])
                listout.append([n[3],n[4],n[6],n[7],n[8]]) # hier ist sigma 1 und 0 drin

    if verbose > 3:
        print('len(listout) 00 :',len(listout))
        for i in listout:
            print("00",i)

    if len(listout) <= 1:  # nothing found
        listout = []
        with open(the_filename, 'rb') as f:
            try:
                all = pickle.load(f, encoding='latin1')            # geht ueber alle eintraege
            except TypeError:
                all = pickle.load(f)            # geht ueber alle eintraege
            for idn,n in enumerate(all):
                if n[0] == qpoint[0] and n[1] == qpoint[1] and n[2] == qpoint[2] and n[5] == func:
                    #print 'xyp2',n
                    #listout.append([n[3],n[4],n[5],n[6],n[7],n[8]])
                    listout.append([n[3],n[4],n[6],n[7],n[8]]) # hier ist sigma 1 und 0 drin
        if verbose > 3:
            print('len(listout) 01 :',len(listout))
            for i in listout:
                print("01",i)

    if len(listout) <= 1 and args.peaks == 1:  # nothing found , take from good
        listout = []
        with open(the_filename, 'rb') as f:
            try:
                all = pickle.load(f, encoding='latin1')            # geht ueber alle eintraege
            except TypeError:
                all = pickle.load(f)            # geht ueber alle eintraege
            for idn,n in enumerate(all):
                if n[0] == qpoint[0] and n[1] == qpoint[1] and n[2] == qpoint[2] and n[5] != 'good':
                    #print 'xyp3',n,n[7]
                    #print 'xyp4',n[8]
                    #print 'nnn good',n
                    #listout.append([n[3],n[4],n[5],n[6],n[7],n[8]])
                    listout.append([n[3],n[4],n[6],n[7],n[8]]) # hier ist sigma 1 und 0 drin
        if verbose > 3:
            print('len(listout) 02 :',len(listout))
            for i in listout:
                print("02",i)

    if len(listout) == 0:
        return False

    if verbose > 3:
        print('#db3# len(listout):',len(listout))

    allpeaks = np.sort(np.unique(np.array(listout)[:,1]))
    if verbose > 3:
        print('#db4# allpeaks:',allpeaks)
        print('#db5# len(listout):',len(listout))

    if len(listout) > 0:
        if verbose > 3:
            print('#db6# allpeaks:',allpeaks,'len(allpeaks)',len(allpeaks))
            for i in allpeaks.astype(int):
                out = np.array(listout)[np.where(np.array(listout)[:,1]==i)]
                if verbose > 3:
                    print('#db7#','i',i)
                    print('#db8#')
                    print(out)
                    #for outt in out:
                    #    print '#db#',outt
                    print('#db9#','--')

    if len(listout) == len(allpeaks):
        outout = []
        for i in allpeaks:
            out = np.array(listout)[np.where(np.array(listout)[:,1]==i)]
            if verbose > 3:
                print('out 88',out)
            outout.append(out[0,2:])
            if verbose > 3:
                print('#db10# -> out',i,outout)
        outout = np.array(outout)
        outout[:,[1, 2]] = outout[:,[2, 1]]
        if verbose > 3:
            print('#db11# -> out',i,outout)
        return np.array(outout).flatten()


    results_out = []
    for i in allpeaks.astype(int):
        out = np.array(listout)[np.where(np.array(listout)[:,1]==i)]


        ######################################################################
        # get rid of all datapoints which have a frequency of 0 or above 30
        ######################################################################
        out1 = out[out[:,2]<30.]
        out2 = out1[out1[:,2]>0.]
        if len(out) != len(out2):
            #print 'out 0:',out
            #print 'out 1:',out1
            #print 'out 2:',out2
            out=out2

        #####################################
        # get out (db)
        #####################################
        out = out[out[:,0].argsort()]
        if len(out) == 0:
            return False

        #####################################
        # get freq
        #####################################
        x = out[:,0]
        y = out[:,2]
        if verbose > 3:
            print("#db13#")
            print("#db14#",i,'DATA:')
            print("#db15#","out:")
            print(out)
            print("#db16#",'i:',i,'x:',x,len(x))
            print("#db16#",'i:',i,'y:',y,len(y))
            print("#db16#",'i: len(out)',len(out))


        const = Model(model_constant)
        weights = square_model(x,1.5)
        weights = (weights/weights.max())*10.
        #print 'i:',i,'x:',x,len(x)
        #print 'i:',i,'y:',y,len(y)
        #print 'i:',i,'weights:',weights,len(weights)
        #np.savetxt("xy",np.transpose((x,y)))
        #np.savetxt("xweights",np.transpose((x,weights)))
        #print 'i:',i,'y.mean():',y.mean()
        result = const.fit(y, x=x, weights=weights,out=y.mean())
        popt0 = lmfit_result_to_popt(result,np.zeros(1),indices=['out'],addindex=False)
        if verbose > 3:
            print('#db16x1 --> old:',result,"new:",popt0[-1])
        result_w = model_constant(for_mdsteps,*popt0)
        if len(out) > 10:
            if True:
                #print 'in'
                square_shifted = Model(model_square_shifted)
                square_shifted.set_param_hint('shift' ,value=x.max(), vary=False)
                weights = square_model(x,1.5)
                weights = (weights/weights.max())*10.
                #np.savetxt("xweightsnew",np.transpose((x,weights)))
                #print 'weights',weights
                result = square_shifted.fit(y, x=x, weights=weights, a=1., shift=x.max(), max=1)
                popt0 = lmfit_result_to_popt(result,np.zeros(1),indices=['a','shift','max'],addindex=False)
                if verbose > 3:
                    print('#db16x2 --> old:',result_w,"new:",popt0[-1])
                #result_w = model_constant(for_mdsteps,*popt0)


        #sys.exit()

        if verbose > 3:
        #if True:
            #print 'FF: popt0 lmfit peak('+str(int(i))+')',popt0,'y.mean():',y.mean()
            print("#db17 fit freq --------------> :",i,'FF: popt0 lmfit peak('+str(int(i))+') FF (THz) ----> ',result_w)
            print("")


        results_out.append(result_w)
        if verbose > 3:
            print('#db# XX 1 ---> results_out:',results_out)

        if args.write_var_files: # for debugging
            if True:
                a = np.transpose(np.array([x,y]))
                b = a[a[:,0].argsort()]
                np.savetxt("var_w_"+str(int(i))+"_"+add_filename+'_aa.dat',b)
                a = np.transpose(np.array([x,np.ones(len(x))*result_w]) )
                b = a[a[:,0].argsort()]
                np.savetxt("var_w_"+str(int(i))+"_"+add_filename+'_lm.dat',b)

        #####################################
        # get aa  !! aa is index 4, GG is index 3 !!!
        #####################################
        x = out[:,0]
        y = out[:,4]

        square_shifted = Model(model_square_shifted)
        square_shifted.set_param_hint('shift' ,value=x.max(), vary=False)
        weights = square_model(x,1.5)

        if len(x) > 1:
            result = square_shifted.fit(y, x=x, weights=weights, a=1., shift=x.max(), max=1)
            popt0 = lmfit_result_to_popt(result,np.zeros(1),indices=['a','shift','max'],addindex=False)
            result_a = model_square_shifted(for_mdsteps,*popt0)
        else:
            result_a = y[0]
        #print 'GG: popt0 lmfit',popt0
        results_out.append(result_a)
        if verbose > 3:
            print('#db18 fit aa: popt0 lmfit peak('+str(int(i))+')',popt0)
            print('#db# XX 2 ---> results_out:',results_out)

        if args.write_var_files: # for debugging
            if True:
                a = np.transpose(np.array([x,y]))
                b = a[a[:,0].argsort()]
                np.savetxt("var_a_"+str(int(i))+"_"+add_filename+'_aa.dat',b)
                a = np.transpose(np.array([x,np.ones(len(x))*result_w]) )
                b = a[a[:,0].argsort()]
                np.savetxt("var_a_"+str(int(i))+"_"+add_filename+'_lm.dat',b)

        #####################################
        # get aa  !! aa is index 4, GG is index 3 !!!
            #@ import matplotlib.pyplot as plt
            #@ plt.plot(x, y*x_to_THz, 'o', label='a Original data', markersize=3)
            #@ plt.plot(x, result.best_fit*x_to_THz, 'g--',label='a lmfit')
            #@ plt.legend()
            #@ #plt.show()
            #@ plt.savefig("var_a_"+str(int(i))+"_"+str(for_mdsteps)+"_"+add_filename+'.pdf', transparent=True,dpi=600,bbox_inches='tight', pad_inches=0)
            #@ #plt.clf()


        #####################################
        # get G
        #####################################
        x = out[:,0]
        y = out[:,3]
        #print 'x1',x
        #print 'y1',y

        ind_corr = y<4.
        y = y[ind_corr]
        x = x[ind_corr]

        #print 'weights',weights
        #print 'x',x
        #print 'y',y
        #xg20000i = x>20000
        #xg20000 = x[xg20000i]
        #yg20000 = y[xg20000i]
        #print 'xg20000',xg20000

        #print 'y',y
        #print 'y<4',y<4.
        #print 'x2',x
        #print 'y2',y
        if len(x) == 0:
            return False
        elif len(x) == 1:
            result_G = y[0]
        else:
            # return a*(x-shift)**2. + max
            #square_shifted = Model(model_square_shifted)
            #square_shifted.set_param_hint('shift' ,value=x.max(), vary=False)
            const = Model(model_constant)
            weights = square_model(x,1.8)
            #np.savetxt("xy",np.transpose((x,y)))
            weights = (weights/weights.max())
            #np.savetxt("xweightsnew",np.transpose((x,weights)))
            #print 'weights',weights
            #result = square_shifted.fit(y, x=x, weights=weights, a=1., shift=x.max(), max=1)
            #print "--------"
            #print(result.fit_report())
            #res_quad = result.chisqr
            #popt0 = lmfit_result_to_popt(result,np.zeros(1),indices=['a','shift','max'],addindex=False)
            result = const.fit(y, x=x, weights=weights,out=y.mean())
            popt0 = lmfit_result_to_popt(result,np.zeros(1),indices=['out'],addindex=False)
            if verbose > 3:
                print(i,'#db20 GG(1): popt0 lmfit',popt0)
            #result_G = model_square_shifted(for_mdsteps,*popt0)
            result_G = model_constant(for_mdsteps,*popt0)
        if verbose > 3:
            print(i,"#db20 GG(2) result_G:",result_G)
        ###############################
        #### append G #################
        results_out.append(result_G)
        if verbose > 3:
            print(i,"#db20 GG(3) from lmfit:",result_G)
        #sys.exit('123er4')
        #result_G_all = model_square_shifted(x,*popt0)
        #if verbose > 3:
        #    print i,"#db20 GG(4)",i,np.transpose(np.array([x,result_G_all]))
        #    print '#db# XX 3 ---> results_out:',results_out

        #print 'result.residual G quad',res_quad
        #@ if True:
        #@     lin_model = Model(model_lin)
        #@     weights = square_model(x,1.5)
        #@     result_lin = lin_model.fit(y, x=x, weights=weights, a=1., m=1)
        #@     res_lin = result_lin.chisqr
        #@     #print printblue('peak '+str(i)+' result.residual G quad/lin '+str(np.round(res_quad/res_lin,2)))
        #@     if res_quad/res_lin >= 1.:
        #@         popt0_lin = lmfit_result_to_popt(result_lin,np.zeros(1),indices=['a','m'],addindex=False)
        #@         result_G_lin = model_lin(for_mdsteps,*popt0_lin)
        #@         #print 'result_G',result_G
        #@         #print 'result_G_lin',result_G_lin
        #@         result_G = np.copy(result_G_lin)


        if args.write_var_files: # for debugging
            if True:
                a = np.transpose(np.array([x/1000.,y]))
                b = a[a[:,0].argsort()]
                if verbose > 3:
                    print('saving',"var_G_"+str(int(i))+"_"+add_filename+'_aa.dat')
                np.savetxt("var_G_"+str(int(i))+"_"+add_filename+'_aa.dat',b)
                #a = np.transpose(np.array([x,np.ones(len(x))*result_w]) )
                #a = np.transpose(np.array([x,result.best_fit*x_to_THz]))
                a = np.transpose(np.array([x/1000.,result.best_fit]))
                b = a[a[:,0].argsort()]
                if verbose > 3:
                    print('saving',"var_G_"+str(int(i))+"_"+add_filename+'_lm.dat')
                np.savetxt("var_G_"+str(int(i))+"_"+add_filename+'_lm.dat',b)

                #plt.plot(x, y*x_to_THz, 'o', label='G Original data', markersize=3)
                #plt.plot(x, result.best_fit*x_to_THz, 'g--',label='G lmfit')
                #plt.legend()
                ##plt.show()
                #plt.savefig("var_G_"+str(int(i))+"_"+str(for_mdsteps)+"_"+add_filename+'.pdf', transparent=True,dpi=600,bbox_inches='tight', pad_inches=0)
                ##plt.clf()





        #plt.plot(x, y, 'bo')
        #plt.plot(x, result.init_fit, 'k--')
        #plt.plot(x, result.best_fit, 'r-')
        #plt.xlim([0,2000])
        #plt.ion()
        if verbose > 3:
            print(i,"#db21")
    if verbose > 3:
        print(i,"#db22")

    return np.abs(np.array(results_out))



def estimate_freqs_from_at_L_L_0_path(qpoint,args):
    ''' returns frequency and estimated separation to other freqs
    --> in the end only l l 0 branch (crossign) is important the ohters are not
    1. get freq from L 0 0 point
    2. get freq from T 0 0 point
    '''
    verbose = False
    #print
    i = qpoint
    path = qpoint_get_l_0_0_or_l_l_0_or_l_l_l_or_t_0_0_or_t_t_t_or_t1_t1_t1_or_t2_t2_t2_from_qpoint(qpoint,args)

    #print 'qpoint::',qpoint,'path',path
    if type(path) == bool:
        return False
    stepsmaxL = 0
    stepsmaxT = 0
    freqoutL = 0
    freqoutT = 0
    qlen = qpoint_to_length_qpoint(q1=i[0],q2=i[1],q3=i[2],args=args)


    the_filename="sql.tmp"
    if not os.path.isfile(the_filename):
        freqoutL = 9.0 # for Al
        freqoutT = 5.8 # for Al
    else:
        with open(the_filename, 'rb') as f:  # laed die sql database
            try:
                all = pickle.load(f, encoding='latin1')            # geht ueber alle eintraege
            except TypeError:
                all = pickle.load(f)            # geht ueber alle eintraege

        for idn,n in enumerate(all):   # geht bei redo shorter nur bis [14 14 32] == T2 (LLL is missing completely)
            #print n[0],n[1],n[2],n[0]==2*args.supercell,n[1]==args.supercell,n[2]== 0
            if n[0]==args.supercell and n[1]==n[2]== 0:
                if n[3] >= stepsmaxL:
                    stepsmaxL = n[3]
                    freqoutL = n[6]
            if n[0]==2*args.supercell and n[1]==args.supercell and n[2]== 0:
                if n[3] >= stepsmaxT:
                    stepsmaxT = n[3]
                    freqoutT = n[6]

    if freqoutL == 0:
        freqoutL = 9.0 # for Al
    if freqoutT == 0:
        freqoutT = 5.8 # for Al

    if path=='l_l_l' or path=='t_t_t':
        qlen = qlen*2

    outfreqL = freqoutL*np.sin(qlen*np.pi/2)
    outfreqT = freqoutT*np.sin(qlen*np.pi/2)
    outfreqLcorr = 0
    outfreqTcorr = 0
    # THe correction seems to make everything worse its really good for al as it is
    corr_T2_T1 = freqoutL/10 # THZ in the middle
    if path == 'l_l_0' or path == 't1_t1_0' or path == 't2_t2_0':
        outfreqLcorr = -corr_T2_T1*np.sin(qlen*np.pi)  # T2
        #outfreqTcorr = corr_T2_T1*np.sin(qlen*np.pi)  # T1
    outfreqLONG = outfreqT + 3.8*np.sin(qlen*np.pi)
    if verbose:
        print("freq L/T along [L 0 0]",freqoutL,freqoutT,'qlen',qlen)
        print('outfreq HIGH (or T2 when along [L L 0]):',outfreqL,'corr',outfreqLcorr)
        print('outfreq LOW  (or T1 when along [L L 0]):',outfreqT,'corr',outfreqTcorr)
        print('outfreq LONG (is HIGH WHEN NOT [L L 0]):',outfreqLONG) #,'corr',outfreqLONG
    return outfreqT+outfreqTcorr, outfreqL+outfreqLcorr,outfreqLONG


def irinasumlifetimeweighted(qp,lt):
    weight=0
    if qp==[1,0,0]:
        weight=6
    if qp==[2,0,0]:
        weight=3
    if qp==[1,1,0]:
        weight=12
    if qp==[1,1,1]:
        weight=4
    if qp==[4,1,0]:
        weight=6
    if qp==[4,2,0]:
        weight=3
    if qp==[5,3,0]:
        weight=12
    if qp==[1,1,4]:
        weight=12
    if qp==[1,1,3]:
        weight=4
    #print 'qp',qp,'weight',weight,'lt',lt,(1./lt)*weight
    if weight > 0:
        #print 'lt',lt
        #print 'tt',1./lt
        import warnings
        warnings.simplefilter('ignore', RuntimeWarning)
        x = (1./lt)
        x[x == np.inf] = 0
        x[x == -np.inf] = 0
        #print 'xx',x
        return x,x*weight
    else:
        return lt*0,lt*0



def print_and_save_lifetimesgood(args = False,verbose=True,printtoscreen=True,the_filename="sql.tmp",return_mdconvergence=False):
    ''' keyword is
            - lifetimesgood l 0 0
            - lifetimesgood t2 t2 t2 mev
            - lifetimesgood l l l meV
            - freqsgood t1 t1 0 meV
    '''
    import shutil
    if printtoscreen:
        print("-------------------------------------------------------")
        print("---------------- def print_and_save_lifetimesgood() ---")
        print("-------------------------------------------------------")
    #print 'return_mdconvergence::==',return_mdconvergence
    dispdir_all = [['l','0','0'],['t','0','0'],['l','l','0'],['t1','t1','0'],['t2','t2','0'],['l','l','l'],['t','t','t'], ['lg','0','0'],['tg','0','0'],['lg','lg','0'],['lg','lg','lg'],['tg','tg','tg'], ['l2','0','0'],['l2','l2','l2']]
    dispdir_all_peaks = [[]]*len(dispdir_all)  # [[1], [1], [1], [1], [1, 2], [1], [1]]
                                                #necessary to save weather the dispdir has 1 peak everywhere; or wheather it has 2 peaks everywhere or only 2 peaks for certain qpoints and 1 peak otherwise
    list_map_qpoint_to_id_dispdir = get_qpoint_mapping_index_dispdir(args,dispdir_all)  # 0 [ 1.  0.  0.  0.  0.]          (== for qpoint [1 0 0] == [l 0 0])
                                                                                        # 1 [ 2.  0.  0.  0.  1.]          (== for qpoint [1 0 0] == [l 0 0])
                                                                                        # .....
                                                                                        # 80 [  1.   1.   1.   5.  80.]    (== for qpoint [1 1 1] == [l l l])
                                                                                        # .....
                                                                                        # 94 [  9.   9.  23.   6.  94.]    (== for qpoint [9 9 23] == [t t t])
                                                                                        # 95 [  8.   8.  24.   6.  95.]    (== for qpoint [8 8 24] == [t t t])
    #for io in list_map_qpoint_to_id_dispdir:
    #    print io

    #print
    #sys.exit('kaiop')
    get_mdstepstocalc(space_fft=False) #,args=args)
    mdstepstocalc_in = args.mdstepstocalc_all
    if printtoscreen:
        print("mdstepstocalc_in:",mdstepstocalc_in,"elements:",len(mdstepstocalc_in),type(mdstepstocalc_in),mdstepstocalc_in.shape)
    mdstepstocalc_all = []                              # list: [ 1000, 2000, 5000, ..., 150000 ]
    fittypes_all = []                                   # list: [ 'good', 'lorenz_68_Hauer2015_xa_sigma_1_bounds_T', ...]
    #for idx,i in enumerate(list_map_qpoint_to_id_dispdir):
    #    print "in:",idx,i
    #########################################################
    # get mdstepstocalc_all
    # get fittypes_all
    # get dispdir_all_peaks
    #########################################################
    newlist = []
    dispdir_all_peaks_max = 0
    lista=[]
    listb=[]

    with open(the_filename, 'rb') as f:  # laed die sql database
        try:
            all = pickle.load(f, encoding='latin1')            # geht ueber alle eintraege
        except TypeError:
            all = pickle.load(f)            # geht ueber alle eintraege
        #print 'ya',type(all),len(all)
        #sys.exit()
        #for idn,n in enumerate(all):   # geht bei redo shorter nur bis [14 14 32] == T2 (LLL is missing completely)
        #    print "n:",n  # n: [1, 0, 0, 1000, 1, 'good', 0.825, 0.05]
        #sys.exit('iopyui')

        for idn,n in enumerate(all):   # geht bei redo shorter nur bis [14 14 32] == T2 (LLL is missing completely)
            if args.verbose > 2:
                print()
                print('all idx:',idn,n)
                print("all idx:",idn,"||","qp:",n[0],n[1],n[2],"mdstep:",n[3],"peaks:",n[4],"ltfq:",n[5],"fq:",n[6],"lt:",n[7])
            #print 'idn:',idn,'/',len(all),"n:",n  # n: [1, 0, 0, 1000, 1, 'good', 0.825, 0.05]
            #sys.exit()
            mdstepstocalc_all.append(n[3])
            fittypes_all.append(n[5])
            if False:
                print("n[0]",n[0])
                print("n[1]",n[1])
                print("n[2]",n[2])
                print("list_map_qpoint_to_id_dispdir",list_map_qpoint_to_id_dispdir)
            id_dispdir,id_qp,peaks_amount = get_qpoint_to_dispdir_idx([n[0],n[1],n[2]],list_map_qpoint_to_id_dispdir)
            if args.verbose > 2:
                print('all idx:',idn,id_dispdir,id_qp,peaks_amount)
            #print '@@@----> n[0],n[1],n[2]:',n[0],n[1],n[2],'-->',id_dispdir,id_qp,peaks_amount
            if list_map_qpoint_to_id_dispdir[id_qp,5] < n[4]:
                list_map_qpoint_to_id_dispdir[id_qp,5]=n[4]
            #if n[4] != 1:
            #    print 'id_qp:',id_qp,'n:',n
            #print "dispdir_all_peaks:",dispdir_all_peaks[id_dispdir]==[]
            if dispdir_all_peaks[id_dispdir] == []:
                dispdir_all_peaks[id_dispdir] = [n[4]]
            else:
                dispdir_all_peaks[id_dispdir] = [dispdir_all_peaks[id_dispdir],n[4]]
                mylist = list_flatten(dispdir_all_peaks[id_dispdir])
                dispdir_all_peaks[id_dispdir]= list(sorted(set(mylist)))
            #print '-->',dispdir_all_peaks
            if n[4] > dispdir_all_peaks_max:
                dispdir_all_peaks_max = n[4]

            removegood1216 = False
            if removegood1216:   # remove good for [t2,t2,0]
                qii=np.array([12,13,14,15,16])

                #if n[1]==n[2] and n[2] == 32 and n[0] in qii and n[5]=='good':
                if n[0]==n[1] and n[2] == 32 and n[5] == 'good' and n[0] in qii:
                    print('idn',idn,'n:',n,'id_qp:',id_qp)
                else:
                    newlist.append(n)

            removemdstep40000 = False
            if removemdstep40000:
                if n[3] == 1000000:
                    print('idn',idn,'n:',n,'id_qp:',id_qp)
                else:
                    newlist.append(n)

            remove4qp = False
            if remove4qp:
                if n[0] == 3 and n[1] == 3 and n[2] == 6 and n[4]==4:
                    print(n)
                    continue
                else:
                    newlist.append(n)

            remove4 = False
            if remove4:
                if n[4]==4:
                    print(n)
                    continue
                else:
                    newlist.append(n)

            remove224 = False  # remove from qp 2 2 4 the 4th peak  ( peak numbering starts with 1) (1,2,3,4 are 4 peaks)
            if remove224:
                if n[0]==n[1]==2:
                    print('kk',n)
                if n[0]==n[1]==2 and n[4]==4:
                    print(n)
                    continue
                else:
                    newlist.append(n)


            #if n[0]==n[1] and n[2] == 32 and n[3] == 3750000 and n[0] in np.array([12]):

            #if n[0]==n[1] and n[2] == 32 and n[3] == 3750000 and n[0] in np.array([12]):

            #if n[0]==n[1] and n[2] == 32 and n[3] == 3750000 and n[0] in np.array([12]):
            #    print 'idn',idn,'n:',n,'id_qp:',id_qp
            #if n[0]==n[1] and n[2] == 32 and n[3] == 15000000 and n[0] in np.array([14]):
            #    print 'idn',idn,'n:',n,'id_qp:',id_qp
            #if n[0]==n[1] and n[2] == 32 and n[0] in np.array([12,13,14,15,16]) and n[5]=='lorenz_68_Hauer2015_xa_sigma_1_bounds_T' and n[4]== 1 and 3750000==n[3]:
            #   lista.append(n)
            #if n[0]==n[1] and n[2] == 32 and n[0] in np.array([12,13,14,15,16]) and n[5]=='lorenz_68_Hauer2015_xa_sigma_1_bounds_T' and n[4]== 2 and 3750000==n[3]:
            #    listb.append(n)
            #    print 'idn',idn,'n:',n,'id_qp:',id_qp

            #if n[0] == 1 and n[1] == 0 and n[2]== 0 and n[5] == 'lorenz_68_Hauer2015_xa_sigma_1_bounds_T' and n[3]==44990: #np.hhhmdstepstocalc_all.max():
            #    print 'nnnnn',n

    #print
    #print "lista:"
    #for i in lista:
    #    print i
    #print
    #print "listb:"
    #for i in listb:
    #    print i
    #sys.exit()

    #sys.exit('ka22344')
    if printtoscreen:
        print('len(newlist)',len(newlist))
    if removegood1216 == True and removemdstep40000 == True: sys.exit('never make both to true')
    if removegood1216 or removemdstep40000 or remove4qp or remove4 or remove224:
        if len(all) != len(newlist):   # remove good for [t2,t2,0]
            with open(the_filename, 'wb') as f:
                pickle.dump(newlist, f)
        sys.exit()

    #for idx,i in enumerate(list_map_qpoint_to_id_dispdir):
    #    print "out:",idx,i
    fittypes_all = np.unique(fittypes_all)
    #dispdir_all_peaks_max = np.array(list_map_qpoint_to_id_dispdir[:,5]).max().astype(int)  WRONG
    mdstepstocalc_all = np.sort(np.unique(mdstepstocalc_all))
    #print "mdstepstocalc_all      (1):",mdstepstocalc_all,"elements:",len(mdstepstocalc_all),type(mdstepstocalc_all),mdstepstocalc_all.shape
    if printtoscreen:
        print("mdstepstocalc_all (ps)    :",mdstepstocalc_all*0.001*args.dt)
        print("dispdir_all_peaks_max     :",dispdir_all_peaks_max)

    if type(mdstepstocalc_in) == np.ndarray:
        if len(mdstepstocalc_all) > len(mdstepstocalc_in):
            mdstepstocalc_all = mdstepstocalc_in
            #pass
    if printtoscreen:
        print("mdstepstocalc_all      (2):",mdstepstocalc_all,"elements:",len(mdstepstocalc_all),type(mdstepstocalc_all),mdstepstocalc_all.shape)
        print("dispdir_all_peaks_max:",dispdir_all_peaks_max)
    #sys.exit()

    #############################################################################
    # make sorted list (sorted by fittypes_all, dispdir_all, mdstepstocalc_all)
    #############################################################################
    fq_fittype_qpoint_mdstep    = np.zeros((len(fittypes_all),len(list_map_qpoint_to_id_dispdir),len(mdstepstocalc_all),dispdir_all_peaks_max))
    aa_fittype_qpoint_mdstep    = np.zeros((len(fittypes_all),len(list_map_qpoint_to_id_dispdir),len(mdstepstocalc_all),dispdir_all_peaks_max))
    lt_fittype_qpoint_mdstep    = np.zeros((len(fittypes_all),len(list_map_qpoint_to_id_dispdir),len(mdstepstocalc_all),dispdir_all_peaks_max))
    peaks_fittype_qpoint_mdstep = np.zeros((len(fittypes_all),len(list_map_qpoint_to_id_dispdir),len(mdstepstocalc_all),dispdir_all_peaks_max))
    print()

    #sys.exit('ka223')
    with open(the_filename, 'rb') as f:
        for idn,n in enumerate(all):  # hier haben wir schon zeilen wie: [32, 10, 0, 9000000, 1, 'lorenz_86_Hauer2015_xl_sigma_1_bounds_T', 4.109, 1.09]
            ##################################################################
            # do this only for mdstes which are in the main list, not others
            ##################################################################

            # remove lorenz_86_Hauer
            if 'lorenz_86_Hauer2015' in n[5]:
                continue

            #remove sigma 0 runs
            if 'sigma_0_' in n[5]:
                continue

            #remove sigma 2 runs
            if 'sigma_2_' in n[5]:
                continue

            # remove unnecessary mdsteps
            if n[3] not in mdstepstocalc_all:
                continue

            # get rid of xl, only take xa
            if len(n[5].split('_xl_sigma')) == 2:
                continue

            #print 'n[5]',n[5],n[5].split('_xl_sigma'),len(n[5].split('_xl_sigma'))
            #print 'nnn',n[5],n[5]

            #print
            #print n[3],"n:",n,n[3] in mdstepstocalc_all
            #id_md = np.where(mdstepstocalc_all==n[3])
            #print 'id_md',id_md,type(id_md)
            #id_md = np.where(mdstepstocalc_all==n[3])[0]
            #print 'id_md',id_md,type(id_md)
            #id_md = np.where(mdstepstocalc_all==n[3])[0][0]
            #print 'id_md',id_md,type(id_md)

            # id_md (array([1]),) <type 'tuple'>
            # id_md [1] <type 'numpy.ndarray'>
            # id_md 1 <type 'numpy.int64'>
            #sys.exit()
            id_md = np.where(mdstepstocalc_all==n[3])[0][0]
            id_ft = np.where(fittypes_all==n[5])[0][0]
            id_dispdir,id_qp,peaks_amount = get_qpoint_to_dispdir_idx([n[0],n[1],n[2]],list_map_qpoint_to_id_dispdir)
            id_peak = int(n[4]-1)   # one peak equals 1 but the index has to be zero

            #peaks_amount = dispdir_all_peaks[id_dispdir]
            #if id_qp == 79 and id_md == 12:
            #    print 'idn',idn,'n:',n,'id_qp:',id_qp,"--------->>>>>>","id_ft,id_qp,id_md,id_peak:",id_ft,id_qp,id_md,id_peak

            #@@% if n[0]==n[1]==12 and n[2]==32 and n[5]=='good':
            #@@%     print 'idn',idn,'n:',n,'id_qp:',id_qp,"--------->>>>>>","id_ft,id_qp,id_md,id_peak:",id_ft,id_qp,id_md,id_peak
            #@@%     print '88888888888'
            #@@%     sys.exit()
            peaks_fittype_qpoint_mdstep[id_ft,id_qp,id_md] = n[4]
            #try:
            #    fq_fittype_qpoint_mdstep[id_ft,id_qp,id_md,id_peak] = n[6]
            #except IndexError:
            #    print "+++",id_ft,id_qp,id_md,id_peak,'n:',n
            #    sys.exit()

            fq_fittype_qpoint_mdstep[id_ft,id_qp,id_md,id_peak] = n[6]
            lt_fittype_qpoint_mdstep[id_ft,id_qp,id_md,id_peak] = n[7]
            aa_fittype_qpoint_mdstep[id_ft,id_qp,id_md,id_peak] = n[8]

    #############################################################################
    # check which mdsteps are not there for every qpoint (do this by checking where zero
    #############################################################################
    if printtoscreen:
        print("fq_fittype_qpoint_mdstep.shape:",fq_fittype_qpoint_mdstep.shape)
        print("fittypes_all.shape            :",fittypes_all.shape,fittypes_all)
        print("len(fittypes_all)             :",len(fittypes_all))
        print("len(dispdir_all)              :",len(dispdir_all))
        print("mdstepstocalc_all.shape       :",mdstepstocalc_all.shape)
    #print
    #print ",fq_fittype_qpoint_mdstep.shape:",fq_fittype_qpoint_mdstep[:,:,15,0]
    #print
    #print ",fq_fittype_qpoint_mdstep.shape:",fq_fittype_qpoint_mdstep[:,:,16,0]
    #print
    #print np.where(fq_fittype_qpoint_mdstep[:,:,:,0]==0)
    #print
    #print np.where(fq_fittype_qpoint_mdstep[:,:,:,0]==0)[2]
    #sys.exit()
    #for id_ft,ft in enumerate(fittypes_all):
    #    for id_dispdir,dispdir in enumerate(dispdir_all):



    dir = "ps_md_convergence"
    if not os.path.exists(dir):
        os.makedirs(dir)
    dir = "ps_summary"
    if not os.path.exists(dir):
        os.makedirs(dir)
    faktor=args.dt*0.001

    if printtoscreen:
        print("dispdir_all:",dispdir_all)
        print("dispdir_all_peaks:",dispdir_all_peaks)
        print("dispdir_all:",dispdir_all)
        print()
    #sys.exit()
    laufidx =1
    irinanoweight = np.zeros(len(mdstepstocalc_all))
    irinaweight = np.zeros(len(mdstepstocalc_all))
    for id_ft,ft in enumerate(fittypes_all):
        print('id_ft,ft',id_ft,ft)
        #if id_ft != 1:
        #    continue
        for id_dispdir,dispdir in enumerate(dispdir_all):                                     # schleife ueber [ l 0 0 ], [ t 0 0], ...
            print('id_dispdir,dispdir',id_dispdir,dispdir)
            in1, in2, in3 = dispdir[0],dispdir[1],dispdir[2]            # is 'l','0','0' or 't2','t2','0'
            qpoints_dispdir = get_all_qpoints([[in1, in2, in3]],args)   # dont change this, it is correct
                                                                        # hier ist qpoints_all_forplotting nur die qunkte eines zweiges!!


            #for id_peak_first in dispdir_all_peaks[id_dispdir]:         # THIS IS THE FIRST  NUMBER OF THE PEAK (either 1 or 2 for al) (TOTAL NUMBER OF PEAKS in THIS DIRECTION [t2 t2 0] has 2 [l l 0] has 1
            #    for id_peak_second in np.arange(id_peak_first)+1:       # THIS IS THE SECOND NUMBER OF THE PEAK ([1 1 32] has [1] and [16 16 32] has [1, 2])
            #        #print 'dispdir:',dispdir,'id_dispdir:',id_dispdir,"dispdir_all_peaks[id_dispdir]:",dispdir_all_peaks[id_dispdir],"id_peak:",id_peak_first,'id_peak_second:',id_peak_second
            #                                                                        # if whole direction has only one peak          -> peak_1_0
            #                                                                        # if several peaks:
            #                                                                        # for all qpoints which have only one peak      -> peak_2_0
            #                                                                        # for all qpoints which have tow peaks(loop)    -> peak_2_1 und peak_2_2
            if id_ft == 2 and id_dispdir == 5:
                print("qpoints_dispdir",qpoints_dispdir)
            if id_ft == 2 and id_dispdir == 6:
                print("qpoints_dispdir",qpoints_dispdir)

            #for id_peak in dispdir_all_peaks[id_dispdir]:
            for idx_peak,id_peak in enumerate(dispdir_all_peaks[id_dispdir]):   # for al: id_peak is mainly 1 but 2 for [t2,t2,0]
                #print '===='
                if len(dispdir_all_peaks[id_dispdir]) > 2:
                    print ('this has to be thought about')
                #if len(dispdir_all_peaks[id_dispdir]) == 2 and id_peak == 1:
                #    idx_peak = 0

                ############################################################
                # for [t2, t2, 0] get the correct peak
                ############################################################
                id_peak_name = id_peak
                changed=False

                ## change the peaks for aluminium, dont do ths for TiN
                #@ AL if dispdir[0] == 't2' and changed==False and id_peak == 1 and ft != 'good':
                #@ AL     id_peak_name = 2
                #@ AL     changed=True
                #@ AL if dispdir[0] == 't2' and changed==False and id_peak == 2 and ft != 'good':
                #@ AL     id_peak_name = 1
                #@ AL     changed=True

                #if dispdir[0] == 't2':
                #    print
                #print 'ft',ft,ft[-1]
                fname = str(ft)+"__"+str(id_peak_name)+"__"+str(in1)+"_"+str(in2)+"_"+str(in3)+"_alat"+str(args.alat)+"_dt"+str(args.dt)+"_"+str("THz")+'.dat'
                in1_exp = np.copy(in1)
                in2_exp = np.copy(in2)
                in3_exp = np.copy(in3)
                if in1=='lg': in1_exp = 'l'
                if in2=='lg': in2_exp = 'l'
                if in3=='lg': in3_exp = 'l'
                if in1=='tg': in1_exp = 't'
                if in2=='tg': in2_exp = 't'
                if in3=='tg': in3_exp = 't'
                fname_exp = str(ft)+"__"+str(id_peak_name)+"__"+str(in1_exp)+"_"+str(in2_exp)+"_"+str(in3_exp)+"_alat"+str(args.alat)+"_dt"+str(args.dt)+"_"+str("THz")+'.dat'
                #print 'fname',fname,'id_peak:',id_peak,dispdir
                file_fq = open("ps_summary/fq_"+fname, 'w')
                file_lt = open("ps_summary/lt_"+fname, 'w')
                file_fq.write("#q\twidth\terr+\terr-\tsigma\t999999999___qpoint\terr_unt\terr_ob\terr_md\n")
                file_lt.write("#q\twidth\terr+\terr-\tsigma\t999999999___qpoint\terr_unt\terr_ob\terr_md\n")


                for iabc,qp in enumerate(qpoints_dispdir):
                    if id_ft in [2,5] and id_dispdir in [5,6]:  # id_dispdir 5 & 6 sind lll und ttt
                        print("qpoints_dispdirA",id_dispdir,qpoints_dispdir,"iabc,qp",iabc,qp,"type(return_mdconvergence)",type(return_mdconvergence),in1,in2,in3)
                    #if type(return_mdconvergence) != bool and return_mdconvergence[0] != qp and str(return_mdconvergence[1]) != ft[-1]:
                    #    continue
                    if type(return_mdconvergence) != bool:
                        if return_mdconvergence[0] != qp or str(return_mdconvergence[1]) != ft[-1]:
                            continue
                    #print 'iabcde,qp',iabc,qp
                    idx_peak_take = idx_peak
                    id_peak_tak = id_peak
                    # hier gehe ich zur zeit noch ueber alle qpunkte des zweiges
                    # hier sollte ich eigentlich nur die qpunkte nehmen welche
                    #id_qp = qpoints_all_list.index(str(qp[0])+"_"+str(qp[1])+"_"+str(qp[2]))
                    id_dispdir,id_qp,peaks_amount = get_qpoint_to_dispdir_idx([qp[0],qp[1],qp[2]],list_map_qpoint_to_id_dispdir)   # peaks amount can be 1 or 2 for al
                    qlen = qpoint_to_length_qpoint(qp[0],qp[1],qp[2],args)
                    #print 'iabcde,qp',iabc,qp
                    #print 'dispdir:',dispdir,'id_dispdir:',id_dispdir,"dispdir_all_peaks[id_dispdir]:",dispdir_all_peaks[id_dispdir]
                    fq = np.zeros((4,len(mdstepstocalc_all)))
                    aa = np.zeros((4,len(mdstepstocalc_all)))
                    lt = np.zeros((4,len(mdstepstocalc_all)))
                    fq[0] = lt[0] = aa[0] = mdstepstocalc_all*faktor  # time in ps
                    #print "mdstepstocalc_all       :",mdstepstocalc_all
                    #print "mdstepstocalc_all*faktor:",mdstepstocalc_all*faktor

                    #if ft == 'lorenz_68_Hauer2015_xa_sigma_1_bounds_T' and in1=='t2' and id_peak > peaks_amount: # and qp[0] in np.array([12,13,14,15,16]):
                    #    print '   bbb id_peak',id_peak,'pa:',peaks_amount,'idx_peak',idx_peak,qp,fq[1][-1]
                    #    idx_peak_take = 0

                    if id_peak > peaks_amount:
                        idx_peak_take = 0

                    #print 'kkk id_peak',id_peak,'pa:',peaks_amount,'idx_peak',idx_peak,"->0",qp
                    #if peaks_amount == 1:
                    #    idx_peak = 0
                    ##@@ if peaks_amount == 2 and id_peak == 2:
                    ##@@     idx_peak = 1
                    #if qp[0]==qp[1]==14 and qp[2] == 32: # and n[3] == 15000000 and n[0] in np.array([14]):
                    #    print 'dispdir:',dispdir,'id_dispdir:',id_dispdir,"dispdir_all_peaks[id_dispdir]:",dispdir_all_peaks[id_dispdir],"len(dispdir_all_peaks[id_dispdir]):",len(dispdir_all_peaks[id_dispdir]),"id_peak:",id_peak,'idx_peak',idx_peak,"qp:",qp,'peaks_amount',peaks_amount #_first,'id_peak_second:',id_peak_second
                    fq[1] = fq_fittype_qpoint_mdstep[id_ft,id_qp,:,idx_peak_take]    # fq values
                    fq[2],fq[3] = get_fq_lt_md_error(fq[1])                     # fq error

                    aa[1] = aa_fittype_qpoint_mdstep[id_ft,id_qp,:,idx_peak_take]    # fq values
                    aa[2],aa[3] = get_fq_lt_md_error(aa[1])                     # fq error

                    #if ft == 'lorenz_68_Hauer2015_xa_sigma_1_bounds_T' and in1=='t2': # and qp[0] in np.array([12,13,14,15,16]):
                    #    print 'kkk id_peak',id_peak,'pa:',peaks_amount,'idx_peak',idx_peak,'idx_peak_take:',idx_peak_take,qp,'fq[1][-1]:',fq[1][-1]

                    if id_ft in [2,5] and id_dispdir in [5,6]:
                        print("qpoints_dispdirB",id_dispdir,qpoints_dispdir,"iabc,qp",iabc,qp,"type(return_mdconvergence)",type(return_mdconvergence),in1,in2,in3,'fq',fq[1][-1]) #,'lt',lt[-1])

                    if fq[1][-1] == 0:
                        #print 'fq0:dispdir:',dispdir,'id_dispdir:',id_dispdir,"dispdir_all_peaks[id_dispdir]:",dispdir_all_peaks[id_dispdir],"len(dispdir_all_peaks[id_dispdir]):",len(dispdir_all_peaks[id_dispdir]),"id_peak:",id_peak,'idx_peak',idx_peak,"qp:",qp,'peaks_amount',peaks_amount #_first,'id_peak_second:',id_peak_second
                        #print 'fq0: dispdir:',dispdir,"qp:",qp,'peaks_amount',peaks_amount #_first,'id_peak_second:',id_peak_second
                        continue
                        #sys.exit()
                    lt[1] = lt_fittype_qpoint_mdstep[id_ft,id_qp,:,idx_peak_take]    # lt values
                    lt[2],lt[3] = get_fq_lt_md_error(lt[1])                     # lt error

                    if id_ft in [2,5] and id_dispdir in [5,6]:
                        print("qpoints_dispdirC",id_dispdir,qpoints_dispdir,"iabc,qp",iabc,qp,"type(return_mdconvergence)",type(return_mdconvergence),in1,in2,in3,'fq',fq[1][-1],'lt',lt[1][-1])

                    #@ if qp[0]==qp[1]==14 and qp[2] == 32: # and n[3] == 15000000 and n[0] in np.array([14]):
                    #@     print
                    #@     print 'id_ft:',id_ft,"lt[2],lt[3]",lt[2][-1],lt[3][-1]
                    #@     if lt[2][-1] > 0.1 or lt[3][-1] > 0.1:
                    #@         print 'fff',fq[1]
                    #@         print 'lll',lt[1]
                    #@         #np.savetxt('aaa'
                    #@         for oo in np.arange(len(lt[1])):
                    #@             print oo,int(fq[0][oo]/faktor),int(fq[0][oo]),fq[1][oo],lt[1][oo],lt[2][oo],lt[3][oo]
                    #@         #print 'lt[1]:',lt[1]
                    #@         #print 'lt[2]:',lt[2]
                    #@ #if idx_peak != 0:
                    #@ #    print 'qp',qp,fq_fittype_qpoint_mdstep[id_ft,id_qp,:]
                    #@ #    sys.exit()

                    ###############################################
                    # save ps_md_convergence
                    # here we would need to remove the old peaks ones which are not there anymore
                    ###############################################
                    if True:
                        #ooo=np.transpose(fq)
                        #print 'ka',np.where(ooo[:,1]==0)
                        #for ippp,ppp in enumerate(np.transpose(fq)):
                        #    print ippp,ppp
                        #sys.exit()
                        namee=str(qp[0])+"_"+str(qp[1])+"_"+str(qp[2])+"_"+ft+"_"+str(in1)+"_"+str(in2)+"_"+str(in3)
                        namee=fname
                        namee=fname+"_"+str(qp[0])+"_"+str(qp[1])+"_"+str(qp[2])+".dat"
                        where = np.where(fq[1]==0)[0]
                        ###############################################
                        # where do we have frequencies which are 0?
                        ###############################################
                        #if len(where) > 0:
                            #print "namee:",namee
                            #print fq[0]
                            #print fq[1]
                            #print where,len(where)
                        #print namee
                        #print np.transpose(fq)
                        np.savetxt("ps_md_convergence/fq_"+namee,np.transpose(fq),fmt='%8.3f %9.3f %9.3f %9.3f')
                        np.savetxt("ps_md_convergence/lt_"+namee,np.transpose(lt),fmt='%8.3f %9.3f %9.3f %9.3f')
                        np.savetxt("ps_md_convergence/aa_"+namee,np.transpose(aa),fmt='%8.3f %9.3f %9.3f %9.3f')

                    #####################################
                    # some analysis
                    #####################################
                    #if qp[0]== 1 and qp[1]==0 and qp[2] == 0: # and n[3] == 15000000 and n[0] in np.array([14]):
                    #    print iabc,'dispdir:',dispdir,'id_dispdir:',id_dispdir,"dispdir_all_peaks[id_dispdir]:",dispdir_all_peaks[id_dispdir],"len(dispdir_all_peaks[id_dispdir]):",len(dispdir_all_peaks[id_dispdir]),"id_peak:",id_peak,'idx_peak',idx_peak,"qp:",qp,'peaks_amount',peaks_amount,'fq[1]',fq[1][-1],'lt[1]',lt[1][-1],namee
                    #if ft == 'lorenz_68_Hauer2015_xa_sigma_1_bounds_T' and in1=='t2': # and qp[0] in np.array([12,13,14,15,16]):
                    #    print 'bbb id_peak',id_peak,'pa:',peaks_amount,'idx_peak',idx_peak,'idx_peak_take:',idx_peak_take,qp,'fq[1][-1]:',fq[1][-1]
                    #####################################
                    # for irina
                    #####################################
                    if ft == 'space_fft_phase_lorenz_68_Hauer2015_xa_sigma_1':
                        #print laufidx,str(dispdir).ljust(17),'id_d:',id_dispdir,"ddap:",dispdir_all_peaks[id_dispdir],"len:",len(dispdir_all_peaks[id_dispdir]),"id_peak:",id_peak,'idx_peak',idx_peak,"qp:",qp,'pks',peaks_amount,'fq[1]',fq[1][-1],'lt[1]',lt[1][-1],lt[1]
                        laufidx+=1
                        tmp1,tmp2= irinasumlifetimeweighted(qp,lt[1])
                        if printtoscreen:
                            pass
                            #print 'tmp1',tmp1
                            #print 'tmp2',tmp2
                        irinanoweight += tmp1
                        irinaweight += tmp2

                    #################################################################################################
                    # write lt_space_fft_phase_lorenz_68_Hauer2015_xa_sigma_1__1__t_t_t_alat4.14_dt1_THz.dat or so
                    #################################################################################################
                    lastlinefq = [qlen,fq[1][-1],fq[2][-1],fq[3][-1],0,0,0,0,0]
                    if qlen == 0.0 and fq[1][-1] < 0.02 and lt[1][-1] < 0.02:
                        pass
                    else:
                        get_fq_lt_summary_line(file_fq,lastlinefq,qp,printtoscreen=False)

                    lastlinelt = [qlen,lt[1][-1],lt[2][-1],lt[3][-1],0,0,0,0,0]
                    if qlen == 0.0 and fq[1][-1] < 0.02 and lt[1][-1] < 0.02:
                        pass
                    else:
                        get_fq_lt_summary_line(file_lt,lastlinelt,qp,printtoscreen=False)
                    if id_ft in [2,5] and id_dispdir in [5,6]:
                        print("qpoints_dispdir",qpoints_dispdir,"iabc,qp",iabc,qp,"type(return_mdconvergence)",type(return_mdconvergence),in1,in2,in3,'lastlinelt',lastlinelt)
                        print("ab",file_lt)
                file_fq.close()
                file_lt.close()
                print("file_fq done",file_fq)


                #print('in1',in1,in2,in3)
                if in1=='lg' or in1=='tg': # and in2=='0'==in3:
                    if not os.path.exists("ps_summary_exp"):
                        os.makedirs("ps_summary_exp")
                    shutil.copyfile("ps_summary/lt_"+fname,"ps_summary_exp/lt_"+fname_exp)
                    shutil.copyfile("ps_summary/fq_"+fname,"ps_summary_exp/fq_"+fname_exp)
                if in1=='t1' or in1=='t2': # and in2=='0'==in3:
                    if not os.path.exists("ps_summary_exp"):
                        os.makedirs("ps_summary_exp")
                    shutil.copyfile("ps_summary/lt_"+fname,"ps_summary_exp/lt_"+fname_exp)
                    shutil.copyfile("ps_summary/fq_"+fname,"ps_summary_exp/fq_"+fname_exp)


                #sys.exit()
    sys.exit()
    if printtoscreen:
        print()
    #print 'irinanoweight',irinanoweight
    #print 'irinaweight',irinaweight/150.
    #print "kk",np.transpose([mdstepstocalc_all/1000.,irinanoweight])
    sys.exit()
    np.savetxt("ps_md_convergence/irina_average_phononlifetime_noweight.dat",np.transpose([mdstepstocalc_all/1000.,irinanoweight]),fmt='%8.3f %9.3f')
    np.savetxt("ps_md_convergence/irina_average_phononlifetime_weight.dat",np.transpose([mdstepstocalc_all/1000.,irinaweight/150.]),fmt='%8.3f %9.3f')

    sys.exit()

    if printtoscreen:
        print('done')
    return

def equivalent_qs(q,show_without_ipi=False,show_lammps_inputfile=False,write_lammps_inputfile=False,N=False,alat=False,dt=False,show_ipi_and_without_ipi=False):
    ''' - finds all the equivalent q vectors for the given q vector
        - q can be a list or a numpy array e.b. [2,0,0] or np.array([2,0,0])
    '''
    #print "type(q):",type(q)
    if type(q) == list:
        q = np.array(q)
    if type(q) != np.ndarray:
        sys.exit("q needs to be a list or a numpy array!")
    if len(q) != 3:
        sys.exit("q needs to be of length 3")

    ergebnis=[];
    dum1=np.array(list(permutations(q)))
    for i in range(0,dum1.shape[0]):
        for m in range(8):
            dum4=m;
            dum31=(-1)**(dum4%2);
            dum4=(dum4-(dum4%2))/2;
            dum32=(-1)**(dum4%2);
            dum4=(dum4-(dum4%2))/2;
            dum33=(-1)**(dum4%2);
            dum3=(dum31, dum32, dum33)
            ergebnis.append((dum3*dum1[i]))
    a=np.array(ergebnis)
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, idx = np.unique(b, return_index=True)

    equivalents = a[idx]
    #print equivalents
    #print equivalents*1j*2*np.pi

    # fixedq: [ 1 0 0 ] and not general  (== old)
    #if q[1] != 0:
    #    sys.exit("currently for [N,0,0] implemented")
    #if q[2] != 0:
    #    sys.exit("currently for [N,0,0] implemented")
    #N = q[0]
    ## aus qvector [N,0,0] wird:
    #equivalents_fixq = np.array([[N,0,0],[0,N,0],[0,0,N],[0,0,-N],[0,-N,0],[-N,0,0]])
    #return equivalents_fixq*1j*2*np.pi
    #return equivalents*1j*2*np.pi   # has + - entries

    ########################################################
    # remove all i = -i
    ########################################################
    out = []
    for i in list(equivalents):
        i=list(i)
        j=[i[0]*-1,i[1]*-1,i[2]*-1]

        #print "OUT:",out
        #print "i:",i,type(i),"j:",j,"(i in out):",i in out,"(j in out):",j in out,"OUT:",out,type(i)
        if (i in out) or (j in out):
            #print "exists"
            continue
        else:
            if len(out) == 0:
                out = [i]
            else:
                out.append(i)
            #print i,"-------------",out #i in out,"_______________",out
        #print
    out= np.array(out)
    if show_without_ipi == True:
        return out


    if show_lammps_inputfile == True or write_lammps_inputfile:
        if type(N) == bool:
            sys.exit('ERROR: you need to specify N if show_lammps_inputfile')
        print(out.shape)  # z.b. 12,3
        outstr = []
        maxvalue = 0
        qs_computes = []
        if write_lammps_inputfile:
            runline = os.popen('grep "^run" in_file_dynamics.in').read().rstrip()
            os.popen("sed -i 's|^run.*||' in_file_dynamics.in").read().rstrip()
            print("runline:",runline)
            target = open("in_file_dynamics.in", 'a')
        for idx,sym_qp in enumerate(out):
            # x
            if sym_qp[0] == 0: o1 = ""
            if sym_qp[0] != 0:
                o1 = str(sym_qp[0])+"*x+";
                if maxvalue < sym_qp[0]: maxvalue = sym_qp[0]
            # y
            if sym_qp[1] == 0: o2 = ""
            if sym_qp[1] != 0:
                o2 = str(sym_qp[1])+"*y+"
                if maxvalue < sym_qp[1]: maxvalue = sym_qp[1]

            # z
            if sym_qp[2] == 0: o3 = ""
            if sym_qp[2] != 0:
                o3 = str(sym_qp[2])+"*z+"
                if maxvalue < sym_qp[2]: maxvalue = sym_qp[2]

            #print idx,sym_qp
            #for xyzidx,xyz in enumerate(sym_qp):

            #if out[idx][idx] == 0:
            #    outstr.append(" ")
            #else:
            #    outstr.append(str(out[idx][idx]))
            out=(o1+o2+o3).replace("+-", "-")[:-1]
            #out2='(('+out+')/${alat}*2.0*PI/'+str(maxvalue)+')'
            #out2='(('+out+')/${alat}*2.0*PI/'+str(N)+')'
            out2='(('+out+')/'+str(alat)+'*2.0*PI/'+str(N)+')'
            #print idx,sym_qp

            ########################################################################
            # schreibe variablen
            ########################################################################
            realcom = [ 'r','i' ]
            for idx2,i in enumerate([ 'cos', 'sin' ]):
                variablename = qpointstring(q)+'_'+str(idx+1)+str(realcom[idx2])
                write1 = "variable "+variablename+" atom "+i+out2
                print(write1)
                if write_lammps_inputfile:
                    target.write(write1+"\n")
                qs_computes.append(variablename)
        ########################################################################
        # schreibe fuer jede variable ein compute
        ########################################################################
        for variablename in qs_computes:
            write2 = "compute c"+variablename+" all reduce sum v_"+variablename
            if write_lammps_inputfile:
                target.write(write2+"\n")
        for variablename in qs_computes:
            write3 = "variable "+variablename+"_ equal c_c"+variablename
            if write_lammps_inputfile:
                target.write(write3+"\n")
        print()
        write4 = 'fix '+qpointstring(q)+' all print '+str(dt)+' "'+ ' '.join(['${'+s+'_}' for s in qs_computes]) + '" append space_fft_'+qpointstring(q)+".dat screen no"
        if write_lammps_inputfile:
            target.write(write4+"\n")
            target.write(runline+"\n")
            target.close()
        print()
        print('thermo_style custom '+ ' '.join(['c_c'+s for s in qs_computes]))
            #print idx,sym_qp,"outstr:","cos"+out2
            #print "           outstr:","sin"+out2
        #print "cos("+outstr[0]+"*x+",outstr[1],"*y+",outstr[2],"*z)"
        return
    if show_ipi_and_without_ipi:
        return out*1j*2*np.pi,out
    else:
        return out*1j*2*np.pi

def print_warning(text):
    print(bcolors.OKGREEN + "##########################################################################################################################################################################################" +bcolors.ENDC)
    print(bcolors.OKGREEN + "##########################################################################################################################################################################################" +bcolors.ENDC)
    print(bcolors.WARNING + "##########################################################################################################################################################################################" +bcolors.ENDC)
    print(bcolors.FAIL + "##########################################################################################################################################################################################" +bcolors.ENDC)
    print(bcolors.FAIL + "##########################################################################################################################################################################################" +bcolors.ENDC)
    print(bcolors.FAIL + text +bcolors.ENDC)
    print(bcolors.FAIL + "##########################################################################################################################################################################################" +bcolors.ENDC)
    print(bcolors.FAIL + "##########################################################################################################################################################################################" +bcolors.ENDC)
    print(bcolors.WARNING + "##########################################################################################################################################################################################" +bcolors.ENDC)
    print(bcolors.OKGREEN + "##########################################################################################################################################################################################" +bcolors.ENDC)
    print(bcolors.OKGREEN + "##########################################################################################################################################################################################" +bcolors.ENDC)

def mysmooth(l,N):
    if (N % 2 == 0): #even
        sys.exit("N has to be odd")
    if N <=2:
        sys.exit("make N > 2")

    d = (N-1)/2
    #print "d:",d
    #print "l:",l
    #print
    out = np.copy(l)
    for idx,i in enumerate(l):
        neg=idx-d
        if neg<0:
            neg=0
        #print idx,i,"window:",l[neg:idx+d+1]
        out[idx] = np.mean(l[neg:idx+d+1])
        #if idx >= d and idx <= len(l)-d:
        #    print idx,i,"window:",l[idx-d:idx+d+1]
        #else:
        #    print idx,i,"window:"
    return out

##########################################################################################
# functions to create space_fft ##########################################################
##########################################################################################
# needs read csv
def get_space_fft_from_positions(filename, qpoints_all = False, mdsteps_per_chunk = 1, scale_data = 1,chunkitornot = False,args=args,filenames=False):
    ''' reads a file in chunks to save memory and creates the sum (for a certain list of
        qpoints) which can be directly fourier transformed.

    - per chunk, amount of lines given by 'args.atoms' (=atoms_in_supercell) times
      mdsteps_per_chunk

    - chunkitornot = {False,'chunk','full'}  where full reads full file and chunk reads in in chunks.

    - example:
        filename = pos_lammps
        qpoint = np.array([1,0,0])
        args.atoms        = 32      # == number of atoms; for a 2x2x2 fcc supercell
        mdsteps_per_chunk = 100     # this would read 32*100 lines per chunk
        scale_data = 4.13*2         # 4.13 Angstrom * 2; this is necessary to get DIRECT coords
    '''
    print("######################################################")
    print("# starting get_space_fft_from_positions()")
    print("######################################################")
    #qpoints_all=np.array([[1,0,0],[2,0,0],[3,0,0],[4,0,0]])
    #qpoints_all=np.array([[2,0,0]])
    #qpoints_all=np.array([[2,0,0],[3,0,0]])
    print('filename (1):',filename)
    if os.path.isfile(filename) != True:
        sys.exit("ERROR:"+str(filenmae)+" does not exist")
    #print "filename:",filename
    #print "scale_data:",scale_data,type(scale_data)
    if type(args.atoms) == bool: sys.exit("ERROR: args.atoms needs to be a number!")


    ############################################################
    # start doint this for al list of qpoints
    # this would only work only when 2*10**8/amountqpoint (say 20) = 1*10**7 is the max number of timesteps
    ############################################################
    qs_all = dict()
    sum_all = dict()
    sum_all_NEW = dict()

    ############################################################
    # start reading in the positions file
    ############################################################
    import time
    start1 = time.time()
    start = time.time()
    #print "fn:",filenames
    if filename not in filenames:
        sys.exit('ERROR: dont know this filename:'+str(filename))



    print("start reader using c engine!")
    if filename == 'dum': scale_data = args.supercell
    if filename == 'posmichael': scale_data = args.supercell

    filesize = os.path.getsize(filename)
    print("filesize             :",filesize,'Byte')
    print("filesize             :",filesize/(1024**2),"MB (if < 10 GB currently no splitting necessary)")

    print('....                 : determining length of file '+printred("(may take a while)"))
    def file_len(fname):
        with open(fname) as f:
            for i, l in enumerate(f):
                if i == 100000:
                    print("counting lines ... 100000       out of ?")
                if i == 1000000:
                    print("counting lines ... 1000000      out of ? ")
                if i == 10000000:
                    print("counting lines ... 10000000     out of ?")
                if i == 100000000:
                    print("counting lines ... 100000000    out of ?")
                pass
        return i + 1
    len_file = file_len(filename)
    print("len_file (lines)     :",len_file)
    len_steps = len_file/args.atoms
    print("len_steps            :",len_steps,'!!!!!!!!!!!')



    if len_file%args.atoms !=0 and args.supercellmultiply == False:
        sys.exit('lines of file is not divisible by number of atoms without remainder!! (1)')
    elif len_file%args.atoms !=0 and args.supercellmultiply != False:
        if (len_file*args.supercellmultiply**3)%args.atoms != 0:
            print("len_file will be multiplied by",args.supercellmultiply,"^3 ADD THIS IN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    chunksize=args.atoms*mdsteps_per_chunk
    if chunksize < 1: sys.exit("ERROR: chunksize needs to be greater 0!")
    print("args.atoms           :",args.atoms)
    print("mdsteps_per_chunk    :",mdsteps_per_chunk)
    print("resuling chunksize(1):",chunksize,'(==lines for',mdsteps_per_chunk,'steps)')

            #sys.exit('lines of file is not divisible by number of atoms without remainder!! (2)')
    # cat POSITIONs| awk '{print $1,$2,$3}' > tmp_xyz
    # mv tmp_xyz > POSITIONs
    # THIS IS STILL NECESSARY!!!
    if os.path.isfile("tmp_xyzz"):
        os.remove("tmp_xyzz")
    if os.path.exists('POSITIONs'):
        print('POSITIONs            : get columns')
        import subprocess
        columns = int(subprocess.check_output('head -1 '+'POSITIONs'+' | wc -w', shell=True,\
                stderr=subprocess.STDOUT))
        print('POSITIONs columns    :',columns,type(columns))

        ####################################################################################
        # this is not necessary with delim_whitespace=True in pandas.read_csv
        ####################################################################################
        #if columns == 3:
        #    print('POSITIONs columns    : has already 3 columns -> awk is still necessary to be able to read_csv... may take a while')
        #    os.popen("cat POSITIONs| awk '{print $1,$2,$3}' > tmp_xyzz; mv tmp_xyzz POSITIONs").read()
        #    print("os.popen POSITIONs done")
        #elif columns == 6:
        #    print('POSITIONs columns    : has 6 columns -> awk is still necessary to be able to read_csv... may take a while')
        #    os.popen("cat POSITIONs| awk '{print $1,$2,$3,$4,$5,$6}' > tmp_xyzz; mv tmp_xyzz POSITIONs").read()
        #    print("os.popen POSITIONs done")
        #else:
        #    sys.exit('can not read POSITONs')
    else:
        sys.exit('no POSITIONs file found')

    if filesize < 1*10**10: # 10 GB
        print("filesize             : smaller 1*10**10 (10 GB) -> no splitting necessary -> setting chunksize to None")
        chunksize=None
    else:
        print("WARNING: this POSITIONS/pos/trj_lammpsnew.out/trj_lammpsnew.out_noJumps_small file ist greater than 10GB, consider the c skript using -fftc since this might take a while.")
    print("---> read_csv        : from filename",filename,printred("(may take a while)"))
    #reader = read_csv(filename, sep=' ', header=None,engine='c',chunksize=chunksize) #,usecols=[1,2,3]) # 3.7 sec
    reader = read_csv(filename,header=None,delim_whitespace=True)  # works for /Users/glensk/Dropbox/Albert/proj/proj_current/__2017.01_phonon_linewidth_al/check_62_LA_POSITIONs/generate_positios/4x4x4sc
    print("---> read_csv        : from filename",filename,"finished!")
    print("resuling chunksize(2):",chunksize)
    #print(reader)
    #sys.exit('kba')



    ############################################################
    # look for existing qpoint files
    ############################################################
    maxavail = np.zeros(len(qpoints_all))
    print("qpoints_all          :",qpoints_all)
    for idx,qpoint in enumerate(qpoints_all):
        qpointstr = qpointstring(qpoint)
        #print "qpoint   :",qpoint,"qpointstr:",qpointstr,"mdsteps_per_chunk:",str(mdsteps_per_chunk)
        #qpointfiles = glob.glob("sum_all_new_"+qpointstr+"__*__"+str(mdsteps_per_chunk)+".npy")
        qpointfiles = glob.glob("space_fft_"+qpointstr+"__*__"+str(mdsteps_per_chunk)+".npy")

        maxavailqpoint = 0
        for j in qpointfiles: # for particular qpoint
            #check = int(j.split("sum_all_new_"+qpointstr+"__")[1].split("__"+str(mdsteps_per_chunk)+".npy")[0])
            check = int(j.split("space_fft_"+qpointstr+"__")[1].split("__"+str(mdsteps_per_chunk)+".npy")[0])
            if check > maxavailqpoint:
                maxavailqpoint=check
        maxavail[idx] = maxavailqpoint
        #print qpointstr, maxavailqpoint,maxavail,int(maxavail.min())
    maxavail = int(maxavail.min())
    print("maxavail             :",maxavail,"maxavail checks if some chunks for particular nuber of mdsteps were performed already")
    if maxavail == 0:
        maxavail = -1;  # just so that nothing goes wrong on step 0
    print("maxavail             :",maxavail,"(is never changed after here; check if inbetweeen files exist)")

    ############################################################
    ############################################################
    # start iterating over infile
    ############################################################
    ############################################################
    print('args.pos_direct      :',args.pos_direct)
    if args.pos_direct == True:
        scale_data = 1.
    print("starting iterating over infile (to get to relative coordinates)... deviding by",scale_data,"to get relative coords; so your largest number in e.g. POSITIONS should be around",scale_data)
    print("for dum: scale data = N;")
    print("POSITIONs files have coords in Angstrom --> for POSITIONs: scale_data = N*alat")
    print('scale_data           :',scale_data)
    print('chunksize            :',chunksize)

    nextcheck = mdsteps_per_chunk
    xxx_ideal_equilibrium = np.zeros((args.atoms,3))
    xxx_ideal_equilibrium = False

    ####################################################################################
    # This loop could be easily parellized to make it quicker
    ####################################################################################
    for chunk,pos in enumerate(reader):
        print('pos')
        ########################################
        # case were we have chunks
        ########################################
        if chunksize != None:   # just to check if there were file made before (intermediate steps)
            if chunk < maxavail:
                if chunk == mdsteps_per_chunk or chunk == nextcheck:
                    nextcheck = chunk + mdsteps_per_chunk
                continue
            # remove nan colums when there was whitespace in the file
            xxx = pos.as_matrix()
            print('xxx (11)')
            print(xxx[:3])





            xxx = xxx[:,~np.all(np.isnan(xxx), axis=0)]
            xxx = xxx[:,:3] # (8000000, 3) is correct; (8000000, 6) includes forces
            #print 'xxx2'
            #print xxx


        ########################################
        # case were we dont have chunks
        ########################################
        if chunksize == None and chunk > 0:
            continue
        if chunksize == None and chunk == 0:
            xxx = reader.as_matrix()

        print('xxx (12)')
        print(xxx[:3])
        #print 'xxx3',xxx.shape
        #print xxx
        xxx = xxx[:,:3] # (8000000, 3) is correct; (8000000, 6) includes forces
        print('chunk',chunk,"len(xxx)",len(xxx))
        #print 'xxx4',xxx.shape
        #print xxx
        #sys.exit('2345443')


        ############################################################################
        # multiply atoms if necessary
        ###########################################################################
        # awk '{a=16.28;print $1,$2,$3,$1+a,$2,$3,$1,$2+a,$3,$1+a,$2+a,$3,$1,$2,$3+a,$1+a,$2,$3+a,$1,$2+a,$3+a,$1+a,$2+a,$3+a}' POSITIONs > tmp
        #print 'xxx1'
        print(xxx)
        print('shape',xxx.shape)
        if args.supercellmultiply:
            if xxx.shape[0] == len_file:
                print("xxx.shape[0] == len_file, this is ok when args.supercellmultiply!")
            else:
                sys.exit('xxx.shape[0] 1= len_file, this is a problem')

            ##################### this here is all made for args.supercellmultiply == 2
            if args.supercellmultiply != 2:
                sys.exit('currently not')
            xxx_tmp = np.zeros((len_file*8,3))
            xxx_tmp_tmp = np.zeros((args.supercellmultiply**3,3))
            ai = args.alat*4. #(args.supercell #*args.supercellmultiply
            #print('ai',ai)
            print('############################# supercell multiply ################')
            print('############################# supercell multiply ################')
            print('############################# supercell multiply ################')
            import datetime
            for idxab in np.arange(len(xxx)):
                if idxab%100000 == 0:
                    print('idxab',idxab,'outof',len(xxx),datetime.datetime.now())
                #print('aaa',xxx[idxab])
                #print('aab',np.array([xxx[idxab][0]          ,xxx[idxab][1]          ,xxx[idxab][2]]))

                xxx_tmp_tmp[0] = np.array([xxx[idxab][0]   ,xxx[idxab][1]   ,xxx[idxab][2]])
                xxx_tmp_tmp[1] = np.array([xxx[idxab][0]+ai,xxx[idxab][1]   ,xxx[idxab][2]])
                xxx_tmp_tmp[2] = np.array([xxx[idxab][0]   ,xxx[idxab][1]+ai,xxx[idxab][2]])
                xxx_tmp_tmp[3] = np.array([xxx[idxab][0]+ai,xxx[idxab][1]+ai,xxx[idxab][2]])

                xxx_tmp_tmp[4] = np.array([xxx[idxab][0]   ,xxx[idxab][1]   ,xxx[idxab][2]+ai])
                xxx_tmp_tmp[5] = np.array([xxx[idxab][0]+ai,xxx[idxab][1]   ,xxx[idxab][2]+ai])
                xxx_tmp_tmp[6] = np.array([xxx[idxab][0]   ,xxx[idxab][1]+ai,xxx[idxab][2]+ai])
                xxx_tmp_tmp[7] = np.array([xxx[idxab][0]+ai,xxx[idxab][1]+ai,xxx[idxab][2]+ai])
                #print('bbb',xxx_tmp_tmp)
                # aus zeile 1 werden zeilen 1-8
                # aus zeile 2 werden zeilen 9-17
                xxx_tmp[idxab*8:idxab*8+8] = xxx_tmp_tmp
                #print()
                #print(xxx_tmp[:30])
                #if idxab == 3:
                #    sys.exit('kkk')
            xxx = xxx_tmp
        #sys.exit('hhh33')
        print('done111')
        xxx = xxx/float(scale_data)
        xxx = np.reshape(xxx,((-1,args.atoms,3)))   # shape: (4000, 256, 3) for last step ...
        if args.pos_modulo == True:
            xxx = np.remainder(xxx,1.)
        #print 'xxx2',xxx
                                                    # ... for last step it could be not 3 but 2 oder 1 (32,3)
        ##############################################
        # check what numbers/positions we have
        ##############################################
        if chunk == 0:
            print("xxx.max()            :",xxx.max())
            if xxx.max() > 1.01:
                #print xxx
                #print 'xxx.shape',xxx.shape
                #print 'now xxy'
                #print xxy
                #print 'xxy.shape',xxy.shape
                #print 'save 1'
                #np.savetxt('kk1.dat',xxx)
                #print 'save 2'
                #np.savetxt('kk2.dat',xxy)
                print_warning('xxx.max() should not be above 1.0')
                sys.exit('xxx.max() should not be above 1.0')
            if xxx.max() < 0.7:
                print_warning('xxx.max() is below 0.7 ---> check scale_data')
                sys.exit('xxx.max() is below 0.7 ---> check scale_data')
            print("xxx.max()            :",xxx.max(),'(is not > 1.0 and < 0.7 and therefore ok!)')
            print('xxx.shape            :',xxx.shape)


        ############################### get transversal to first BZ
        ############################### get transversal to first BZ CURRENT
        print('xxx args.space_fft_phase:',args.space_fft_phase)
        if args.space_fft_phase == True:
            print('############################# make space_fft_phase ################')
            print('############################# make space_fft_phase ################')
            print('############################# make space_fft_phase ################')
            print('make xxx_ideal ...')
            N=args.supercell
            if type(xxx_ideal_equilibrium) == bool:
                xxx_ideal_equilibrium = np.round([(xxx[0]*2*N+0.5)%(2*N)-0.5])/2.
                print('have now xxx_ideal_equilibrium positions ...')
            print('make xxx_dr ...')
            xxx_dr = ((xxx*2*N+0.5)%(2*N)-0.5)/2. - xxx_ideal_equilibrium

            def get_sum_from_sigma_qs(xxx,xxx_dr,qpoint,qpoint_is_longitudinal_):
                '''
                    the magnitude of the sigma does not matter at all
                    the direction of the sigma [ 1 0 0 ] vs [ -1 0 0] ??
                '''
                print('########################################################')
                print('# get_sum_from_sigma_qs(xxx,xxx_dr,qpoint)')
                print('########################################################')

                #qs_all = equivalent_qs(qpoint)
                #qs_all_n = equivalent_qs(qpoint,show_without_ipi=True)
                #qs_all_a,qs_all_b = equivalent_qs(qpoint,show_ipi_and_without_ipi=True)
                qs_all,sigma = equivalent_qs(qpoint,show_ipi_and_without_ipi=True)

                print('qpoint   -->',qpoint)
                print('sigma    -->')
                print(sigma)
                print('qs_all   -->')
                print(qs_all)
                print('info     -->',len(qs_all),xxx.shape,xxx.shape[0],xxx_dr.shape)

                # for longitudinal use qs_all[idx] == sigma[idx]
                #        qpint   :   qpoint == sigma
                # z.B. [ 1 0 0 ] : [ 1 0 0 ]
                # z.B. [ 2 0 0 ] : [ 2 0 0 ]
                # z.B. [ 3 0 0 ] : [ 3 0 0 ]
                # z.B. [ 3 3 0 ] : [ 3 3 0 ]
                # z.B. [ 4 4 4 ] : [ 4 4 4 ]
                # for transversal use longitudinal qpoints but go through other sigmas!
                #        qpint   :   all qpoints
                # z.B. [ 1 0 0 ] : [1 0 0] [ 0 1 0 ], [ 0 0 1]
                # z.B. [ 0 1 0 ] : [ 1 0 0 ], [ 0 0 1]
                # z.B. [ 1 1 1 ] : [ 1 1 -1 ], [ 1 -1 1 ],[ 1 -1 -1 ]
                # if True: # MAKE Longitudinal (all longitudinal can be saved as it is)
                # if True: # MAKE Longitudinal (all longitudinal can be saved as it is)
                # if True: # MAKE Longitudinal (all longitudinal can be saved as it is)
                ps = np.zeros((len(qs_all),xxx.shape[0]),dtype=np.complex64)
                writedown = False
                #qpointstr_ = qpointstring(qpoint)
                #if os.path.isfile('space_fft_phase_'+qpointstr_+'.npy'):
                #    print('space_fft_phase_'+qpointstr_+'.npy does already exist')
                #else:
                #    print('space_fft_phase_'+qpointstr_+'.npy does NOT YET exist')
                if True:
                    for idx,qs in enumerate(qs_all):
                        qpointstr = qpointstring(sigma[idx])
                        print('(LONGITUDINAL) idx',idx+1,'/',len(qs_all),'(BEDINGUNG for LONGITUDINAL: qs_all[idx] == sigma[idx] == ',sigma[idx],')')
                        writedown = True
                        if idx == 0:
                            print('XXXXXXXXX LONGITUDINAL XXXXXXXXXXXXXXXX')
                            print('XXXXXXXXX LONGITUDINAL XXXXXXXXXXXXXXXX')
                        xxx_phase = np.sum(sigma[idx]*xxx_dr,axis=2)
                        #print 'idx,qs',idx,qs,qpoint,qs_all_n[idx]
                        #print
                        #xxx_amplitude=xxx_phase*np.exp(np.sum(xxx*N*qs,axis=2));
                        xxx_amplitude=xxx_phase*np.exp(np.sum(xxx*qs_all[idx],axis=2)); # version lange gueltig until 2018_04_03
                        #xxx_amplitude=xxx_phase*np.exp(np.sum(xxx_ideal_equilibrium*qs_all[idx],axis=2));
                        ps[idx] = np.sum(xxx_amplitude,axis=1)
                        if False: # only for debugging
                            bla=ps[idx]
                            out=np.abs(fft(bla))**2.;
                            print('qpointstr',qpointstr)
                            print(printred("space_fft_phase_L_"+qpointstr))
                            np.savetxt("space_fft_phase_L_"+qpointstr,out)

                if writedown:
                    qpointstr = qpointstring(qpoint)
                    np.save("space_fft_phase_"+qpointstr,ps)
                    #np.save("space_fft_phase_"+qpointstr+"_L",ps)
                    print("(LONGITUDINAL) space_fft_phase_"+qpointstr,'written!')
                    #print "(LONGITUDINAL) space_fft_phase_"+qpointstr+"_L",'written!'
                print('XXXXXXXXX LONGITUDINAL DONE XXXXXXXXXXXXXXXX')
                print('XXXXXXXXX LONGITUDINAL DONE XXXXXXXXXXXXXXXX')

                # if True: # MAKE Transversal (here it would be best to give it the name from higher BZ)
                # if True: # MAKE Transversal (here it would be best to give it the name from higher BZ)
                if True: # MAKE Transversal (here it would be best to give it the name from higher BZ)
                    print('XXXXXXXXX TRANSVERSAL XXXXXXXXXXXXXXXX')
                    print('XXXXXXXXX TRANSVERSAL XXXXXXXXXXXXXXXX')
                    qpt = qpoint_map_longitudinal_to_transversal(qpoint,args=args)
                    if qpt is None:
                        qpt = np.array([[0,0,0]]) # in case of no clue
                    print('qpt',qpt)
                    #sys.exit('out123')
                    for tqidx,tq in enumerate(qpt): # no effect when just one transversal branch
                        if os.path.isfile('space_fft_phase_'+qpointstring(tq)+'.npy'):
                            print("The qpointfile",qpointstring(tq)+".npy"," exists")
                            continue
                        T = False
                        T1 = False
                        T2 = False
                        if len(qpt) == 1:  # normal transversal mode (L00 and LLL)
                            ps = np.zeros((len(qs_all)*(len(qs_all)-1),xxx.shape[0]),dtype=np.complex64)
                            T = True
                            textout = 'T'
                        if len(qpt) == 2: # tow transversal modes
                            textout = ''
                            if (np.array(tq) == np.array(qpt[0])).all():
                                T1 = True
                                ps = np.zeros((6,xxx.shape[0]),dtype=np.complex64)
                                textout = 'T1'
                            else:                       # T2
                                T2 = True
                                ps = np.zeros((len(qs_all)*(len(qs_all)-1)-6,xxx.shape[0]),dtype=np.complex64)
                                textout = 'T2'
                        #sys.exit('kaa')
                        print('textout',textout)
                        idxnext=0
                        for idx,qs in enumerate(qs_all):
                            sigma_sub = np.delete(sigma, idx,0)
                            qs_all_sub = np.delete(qs_all,idx,0)
                            #print 'sigma_sub -->'
                            #print sigma_sub
                            #print 'qs_all_sub --->'
                            #print qs_all_sub
                            for idx_sub,qs_sub in enumerate(qs_all_sub):
                                qs_sub_showarray = ((qs_sub/(2*np.pi)).imag).astype(int)
                                #print 'idx,qs',idx,qs,'idx_sub,qs_sub',idx_sub,qs_sub_showarray
                                goon = False
                                #print '(TRANSVERSAL) '+textout+' idx',idxnext+1,'/',len(ps),' ---> qpoint',qpoint,'qpt',qpt,'tq',tq,textout,np.abs(sigma[idx]),np.abs(sigma_sub[idx_sub])
                                #print '...',(np.abs(sigma[idx]) == np.abs(sigma_sub[idx_sub])).all()
                                if T:
                                    goon = True
                                if T1 and (np.abs(sigma[idx]) == np.abs(sigma_sub[idx_sub])).all():
                                    #print 'check t1',(np.abs(sigma[idx]) == np.abs(sigma_sub[idx_sub]))
                                    #print '--...',(np.abs(sigma[idx]) == np.abs(sigma_sub[idx_sub])).all
                                    goon = True
                                if T2 and not (np.abs(sigma[idx]) == np.abs(sigma_sub[idx_sub])).all():
                                    goon = True
                                    #print 'check t2'

                                if goon:
                                    #print '(TRANSVERSAL) '+textout+' idx',idxnext+1,'/',len(ps),' ---> qpoint',qpoint,'qpt',qpt,'tq',tq,textout,np.abs(sigma[idx]),np.abs(sigma_sub[idx_sub]),'goon',goon
                                    print('(TRANSVERSAL) '+textout+' idx',idxnext+1,'/',len(ps),' ---> qpoint',qpoint,'tq',tq,textout,np.abs(sigma[idx]),sigma_sub[idx_sub],'goon',goon)
                                    #print '(TRANSVERSAL) '+textout+' ---> idx,qs,qs_all[idx]',idx,qs,qs_all[idx]
                                    #print '(TRANSVERSAL) '+textout+' ---> tq',tq,'len(qpt)',len(qpt),'qpt',qpt
                                    #print '(TRANSVERSAL) '+textout+' ---> idxnext',idxnext,"out of",ps.shape
                                    #print '(TRANSVERSAL) '+textout+' ---> sigma[idx]',sigma[idx]
                                    #print '(TRANSVERSAL) '+textout+' ---> idx_sub',idx_sub,'sigma_sub[idx_sub] ',sigma_sub[idx_sub]
                                    #print 'qs_all_sub[idx_sub]',qs_all_sub[idx_sub]
                                    xxx_phase = np.sum(sigma_sub[idx_sub]*xxx_dr,axis=2)
                                    xxx_amplitude=xxx_phase*np.exp(np.sum(xxx*qs_all[idx],axis=2)); # same as long ... lange gueltig until 2018_04_03
                                    #xxx_amplitude=xxx_phase*np.exp(np.sum(xxx_ideal_equilibrium*qs_all[idx],axis=2)); # same as long
                                    ps[idxnext] = np.sum(xxx_amplitude,axis=1)
                                    idxnext+=1

                                    if False: # only for debugging
                                        bla=ps[idxnext-1]
                                        out=np.abs(fft(bla))**2.;
                                        print('sigma[idx]',sigma[idx])
                                        print('sigma_sub[idx_sub]',sigma_sub[idx_sub])
                                        qpointstr = qpointstring(sigma[idx])
                                        print('qpointstr',qpointstr)
                                        sigmastr = qpointstring(sigma_sub[idx_sub])
                                        print('sigmastr',sigmastr)
                                        print(printred("space_fft_phase_"+textout+"_"+qpointstr+"__"+sigmastr))
                                        np.savetxt("space_fft_phase_"+textout+"_"+qpointstr+"__"+sigmastr,out)
                                        print('no duplicates in T2 found and no duplicates in T1 found')
                        qpointstr = qpointstring(tq)
                        np.save("space_fft_phase_"+qpointstr,ps)
                        print("(TRANSVERSAL) "+textout+" space_fft_phase_"+qpointstr,'written!')
                    print('XXXXXXXXX TRANSVERSAL DONE XXXXXXXXXXXXXXXX')
                    print('XXXXXXXXX TRANSVERSAL DONE XXXXXXXXXXXXXXXX')
                    print()

                return ps


            ##############################################################
            # HERE WE CAN EASILY GO THROUGH ALL LONGITUDINAL WAVE VECTORS
            ##############################################################
            #for qpoint in np.array([[1,0,0],[2,0,0],[3,0,0],[4,0,0]]):
            #for qpoint in np.array([[1,1,1],[2,2,2]]):
            #for qpoint in np.array([[1,1,0],[2,2,0],[3,3,0],[4,4,0]]):
            #for qpoint in np.array([[1,0,0]]):

            ####################################################################
            # qpoint_is_longitudinal (2nd line) filters all longitudinal qpoints
            ####################################################################
            #for qpoint in np.array([[4,4,0]]):
            #for qpoint in np.array([[1,0,0]]):
            if True:
                #for qpoint in np.array([[1,1,0]]): # for 1 1 8
                for qpoint in qpoints_all:
                    qpointstr_ = qpointstring(qpoint)
                    if os.path.isfile('space_fft_phase_'+qpointstr_+'.npy'):
                        print('space_fft_phase_'+qpointstr_+'.npy does already exist')
                    else:
                        print('space_fft_phase_'+qpointstr_+'.npy does NOT YET exist')
                        qpoint_is_longitudinal_ = qpoint_is_longitudinal(qpoint,passhigherBZ=False,args=args)
                        print("qpoint_is_longitudinal? (or maybe out of first BZ? or first BZ alredy calculated?)",qpoint_is_longitudinal_)
                        print("in case you want to calculate a transversal qpoint only, you need to delete firest the corresponding longitudinal space_fft_phase qpoint, the one selected transveral qpoint will be calculated and additionally the deleted longitudinal")
                        if qpoint_is_longitudinal_:
                            #print 'qpoint --> (make in first BZ)',qpoint
                            #####################################################################################################
                            # --> THIS CREATES THE space_fft_phase_xxx.npy FILE
                            # --> THIS CREATES THE space_fft_phase_xxx.npy FILE
                            ps = get_sum_from_sigma_qs(xxx,xxx_dr,qpoint,qpoint_is_longitudinal_)   # --> THIS CREATES THE space_fft_phase_xxx.npy FILE
                            #####################################################################################################
            #sys.exit('doneeee')

        print('############################# make space_fft (high bz trans) ################')
        print('############################# make space_fft (high bz trans) ################')
        print('############################# make space_fft (high bz trans) ################')

        mdstepbegin = chunk*mdsteps_per_chunk
        mdstepend = chunk*mdsteps_per_chunk+xxx.shape[0]
        print("chunk (=index):"+str(chunk)+"  mdstepsbegin:"+str(mdstepbegin)+"  mdstepsend:"+str(mdstepend),"nextcheck:",nextcheck,"len_steps:",len_steps)

        ##################################################################################
        # calculate for intermediate qpoints the space_fft
        ##################################################################################
        print("this qpoints will be calculated")
        #for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste [ 1,0,0] [2,0,0] ] [3,0,0]
        #    print(idx+1,"/",len(qpoints_all),qpoint) #"qpoint:",qpoint,len(qpoints_all),

        print()
        print("now lets calculate them")
        for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste [ 1,0,0] [2,0,0] ] [3,0,0]
            print(idx+1,"/",len(qpoints_all),qpoint) #"qpoint:",qpoint,len(qpoints_all),
            qpointstr = qpointstring(qpoint)
            if os.path.isfile('space_fft_'+qpointstr+'.npy'):
                print('space_fft_'+qpointstr+'.npy does already exist')
            else:
                print('space_fft_'+qpointstr+'.npy does NOT YET exist',idx+1,"/",len(qpoints_all),qpoint) #"qpoint:",qpoint,len(qpoints_all),
                #print idx+1,"/",len(qpoints_all),"qpoint:",qpoint,len(qpoints_all),
                qs_all[idx] = equivalent_qs(qpoint) # [1 0 0] [0 1 0] [0 0 1]
                sum_all[idx] = np.zeros((len(qs_all[idx]),mdstepend-mdstepbegin),dtype=np.complex64)  # np.complex_32 would also be an option!
                if chunksize == None:
                    print('will be written since no chunks')

                for ind_qs,qs in enumerate(qs_all[idx]):  # fuer alle symmetrieequivalenten qpoints also fuer [1 0 0] haben wir [1 0 0] [0 1 0] [0 0 1]
                    #print "qs:",qs
                    # calculate via exp when dumping every chunk (amount number of steps)
                    ### MAIN !!!!!!!!!!!!!!!!!!!  (inner sum to make one column per md step, outer sum to make one complex number per mdstep)

                    # check q und -q: also in /Users/glensk/Dropbox/proj/proj_current/__2017.01_phonon_linewidth_al/__2017.01_phonon_lifetimes_4_ab2017/check_32_q_minus_q_auswirkung_timeinversion
                    #np.save("xxx",xxx)
                    #sys.exit()
                    # z.B. qs = np.array([0,0,1])*2*np.pi*1j
                    # array([ 0.+0.j        ,  0.+0.j        ,  0.+6.28318531j])
                    # z.B. qs = np.array([0,0,-1])*2*np.pi*1j
                    # array([ 0.+0.j        ,  0.+0.j        , -0.-6.28318531j])
                    #
                    #
                    # In [19]: qs = np.array([0,0,1])*2*np.pi*1j
                    #
                    # In [20]: ampl = np.sum(np.exp(np.sum((xxx*qs),axis=2)),axis=1);print ampl;np.save("100ampl",ampl)
                    # [ 0.89848075-0.16619753j  0.90832258-0.17126346j  0.91808158-0.176336j   ...,
                    #  -0.15524963-0.10834658j -0.15839359-0.1109407j  -0.16157893-0.11330545j]
                    #
                    # In [21]: qs = np.array([0,0,-1])*2*np.pi*1j
                    #
                    # In [22]: ampl = np.sum(np.exp(np.sum((xxx*qs),axis=2)),axis=1);print ampl;np.save("100amplmin",ampl)
                    # [ 0.89848075+0.16619753j  0.90832258+0.17126346j  0.91808158+0.176336j   ...,
                    #  -0.15524963+0.10834658j -0.15839359+0.1109407j  -0.16157893+0.11330545j]
                    #
                    #print "qpoint:",qpoint
                    #print "qs_all[idx][ind_qs]:",qs_all[idx][ind_qs]
                    #sys.exit()

                    #sum_all_NEW[idx][ind_qs] = np.sum(xxx_ph*np.exp(np.sum((xxx*float(N)*qs_all[idx][ind_qs]),axis=2)),axis=1)
                    #sum_all_NEW[idx][ind_qs] = xxx_ph*np.exp(np.sum((xxx*float(N)*qs_all[idx][ind_qs]),axis=2))
                    #sum_all_NEW[idx][ind_qs] = np.exp(np.sum((xxx*float(N)*qs_all[idx][ind_qs]),axis=2))

                    #kkk = np.exp(np.sum((xxx*float(N)*qs_all[idx][ind_qs]),axis=2))
                    #print "kkk.shap:",kkk.shape,kkk
                    #print
                    #print "xxx.phase.shape:",xxx_phase.shape,xxx_phase
                    #print
                    #ab = xxx_phase*kkk
                    #print "---->",ab.shape,ab
                    #print
                    #ac = np.sum(ab,axis=1)
                    #print ".....> 2:",ac.shape,ac
                    #sys.exit()
                    #sum_all_NEW[idx] = np.sum(xxx_phase*np.exp(np.sum((xxx*float(N)*qs_all[idx][ind_qs]),axis=2)),axis=1)

                    eulerforexp = False  # --> gibt beides exakt die gleiche streuamplitude
                    if eulerforexp == False:
                        sum_all[idx][ind_qs] = np.sum(np.exp(np.sum((xxx*qs_all[idx][ind_qs]),axis=2)),axis=1)

                    # calculate via exp when dumping only at the very end
                    # sum_all[idx][ind_qs,mdstepbegin : mdstepend] = np.sum(np.exp(np.sum((xxx*qs_all[idx][ind_qs]),axis=2)),axis=1
                    # OUTPUT TO SCREEN WHEN NECESSARY
                    #print "sum_all:",chunk,idx,ind_qs,sum_all[idx][ind_qs].shape,"<-- anzahl mdsteps"
                    #pass

                    # calculate via eulers formula (sin, cos) when dumping every chunk (amount number of steps)
                    ### MAIN !!!!!!!!!!!!!!!!!!!  (first sum to make one column per md step, second sum to make one complex number per mdstep)
                    if eulerforexp == True:
                        a = np.sum((xxx*qs_all[idx][ind_qs]),axis=2)  # np.sum is over all atoms of one particular mdstep to give one column (per atom one number (complex))
                        sum_all[idx][ind_qs].real = np.sum(np.cos(a.imag),axis=1)  # correct ralpart due to euler
                        sum_all[idx][ind_qs].imag = np.sum(np.sin(a.imag),axis=1)  # correct imagpart to euler
                    #print "realp. :",chunk,idx,ind_qs,sum_all[idx][ind_qs].shape,"<-- anzahl mdsteps"
                    #print "imagp. :",chunk,idx,ind_qs,sum_all[idx][ind_qs].shape,"<-- anzahl mdsteps"
                    #print
                if chunksize == None:
                    print('writeit')
                    qs_all[idx] = equivalent_qs(qpoint)
                    qpointstr = qpointstring(qpoint)
                    np.save("space_fft_"+qpointstr,sum_all[idx])
        #print "----------------",chunk,"start",sum_all[idx].shape
        #print sum_all
        #print "----------------",chunk,"fin",sum_all[idx].shape
        #if chunk == 2:
        #    sys.exit()

        ##################################################################################
        # save intermediate qpoins
        # this is also done when we have no chunks
        ##################################################################################
        if chunksize != None:   # just to check if there were file made before (intermediate steps)
            for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
                qs_all[idx] = equivalent_qs(qpoint)
                qpointstr = qpointstring(qpoint)
                #np.save("sum_all_new_"+qpointstr+"__"+str(chunk)+"__"+str(mdsteps_per_chunk),sum_all[idx])
                np.save("space_fft_"+qpointstr+"__"+str(chunk)+"__"+str(mdsteps_per_chunk),sum_all[idx])

        ##################################################################################
        # for every 10-100 chunks, put all qpoints together and remove the last ones.
        # e.g. chunk  0-9  put all together to 9
        # e.g. chunk 10-19 put all together to 19
        # ONLY WHEN CHUNKS
        ##################################################################################
        if chunksize != None:   # just to check if there were file made before (intermediate steps)
            if chunk == mdsteps_per_chunk or chunk == nextcheck:

                loadall = np.arange(nextcheck-mdsteps_per_chunk,nextcheck+1)
                #print loadall
                for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
                    qpointstr = qpointstring(qpoint)
                    #a = np.load("sum_all_new_"+qpointstr+"__"+str(loadall[0])+"__"+str(mdsteps_per_chunk)+".npy")
                    a = np.load("space_fft_"+qpointstr+"__"+str(loadall[0])+"__"+str(mdsteps_per_chunk)+".npy")
                    for b in loadall[1:]:
                        #b = np.load("sum_all_new_"+qpointstr+"__"+str(b)+"__"+str(mdsteps_per_chunk)+".npy")
                        b = np.load("space_fft_"+qpointstr+"__"+str(b)+"__"+str(mdsteps_per_chunk)+".npy")
                        a = np.concatenate((a, b), axis=1)
                    # this will overwrite the last chunkfile
                    #np.save("sum_all_new_"+qpointstr+"__"+str(chunk)+"__"+str(mdsteps_per_chunk),a)
                    np.save("space_fft_"+qpointstr+"__"+str(chunk)+"__"+str(mdsteps_per_chunk),a)
                    for b in loadall[:-1]:  # delete all but last cunkfile
                        os.remove("space_fft_"+qpointstr+"__"+str(b)+"__"+str(mdsteps_per_chunk)+".npy")

                nextcheck = chunk + mdsteps_per_chunk
        ##################################################################################
        # loop over intermediate qpoins done
        ##################################################################################

    ##################################################################################
    # get last steps
    # ONLY WHEN CHUNKS
    ##################################################################################
    if chunksize != None:   # just to check if there were file made before (intermediate steps)
        if chunk < nextcheck:
            loadall = np.arange(nextcheck-mdsteps_per_chunk,chunk+1)
            print(loadall)
            for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
                qpointstr = qpointstring(qpoint)
                #a = np.load("sum_all_new_"+qpointstr+"__"+str(loadall[0])+"__"+str(mdsteps_per_chunk)+".npy")
                a = np.load("space_fft_"+qpointstr+"__"+str(loadall[0])+"__"+str(mdsteps_per_chunk)+".npy")
                for b in loadall[1:]:
                    #b = np.load("sum_all_new_"+qpointstr+"__"+str(b)+"__"+str(mdsteps_per_chunk)+".npy")
                    b = np.load("space_fft_"+qpointstr+"__"+str(b)+"__"+str(mdsteps_per_chunk)+".npy")
                    a = np.concatenate((a, b), axis=1)
                # this will overwrite the last chunkfile
                #np.save("sum_all_new_"+qpointstr+"__"+str(chunk)+"__"+str(mdsteps_per_chunk),a)
                np.save("space_fft_"+qpointstr+"__"+str(chunk)+"__"+str(mdsteps_per_chunk),a)
                for b in loadall[:-1]:  # delete all but last cunkfile
                    os.remove("space_fft_"+qpointstr+"__"+str(b)+"__"+str(mdsteps_per_chunk)+".npy")


    print("DONE! mdstepend:",mdstepend,"------------------------------------------------")
    end = time.time()
    print("#1) read in data :",bcolors.FAIL + str((end - start)) +" sec"+ bcolors.ENDC)
    print("ranaming files:")
    ##################################################################################
    # rename space_fft_6_2_0__x__1000.npy to space_fft_6_2_0.npy
    ##################################################################################
    #sys.exit()
    if chunksize != None:   # just to check if there were file made before (intermediate steps)
        for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
            qpointstr = qpointstring(qpoint)
            #file = glob.glob("sum_all_new_"+qpointstr+"__*")
            file = glob.glob("space_fft_"+qpointstr+"__*")
            print(len(file),file)
            if len(file) == 1:
                #os.rename(file[0], 'sum_all_new_'+qpointstr+'.npy')
                os.rename(file[0], 'space_fft_'+qpointstr+'.npy')

    #for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
    #    qpointstr = qpointstring(qpoint)
    #    sum_all[idx] = sum_all[idx][:,:mdstepend]
    #    #np.save("sum_all_new_"+str(qpoint[0])+str(qpoint[1])+str(qpoint[2]),sum_all)
    #    np.save("sum_all_new_"+qpointstr,sum_all[idx][:,:mdstepend])
    #    #np.savetxt("sum_all.txt",sum_all)
    #    print "sum_all[idx].shape:",sum_all[idx].shape
    #    #for ind_qs,qs in enumerate(qs_all):
    #    #    np.savetxt('sum'+str(qpoint[0])+str(qpoint[1])+str(qpoint[2])+str(ind_qs)+'_'+str(mdstepend)+'.out',sum)
    return

def get_space_fft_prepare_lammps_log(filename):
    ''' this reads directly from lammps.log; comments above and below the output should be no problem '''
    print("starting get_space_fft_from_lammps_log ...")
    #if type(qpoint) == bool:
    #    sys.exit("ERROR qpoint is bool")
    print("filename     :",filename)
    skipfirst=False
    from itertools import islice
    with open(filename) as myfile:
        head = list(islice(myfile, 2000))
        for idx,i in enumerate(head):
            if "Memory usage per processor =" in i:
                skipfirst = idx+2
    skipfirst = int(skipfirst)

    print("skipfirst    :",skipfirst)
    if type(skipfirst)== bool:
        sys.exit("ERROR: did not find skipfirst")


    skiplast = os.popen('tail -200 '+filename+' | grep -n "Loop time of" | sed \'s|:.*||\'').read()
    skiplast = skiplast.rstrip()
    print("skiplast1    :",skiplast) #,":",type(skiplast)
    dothisall = True

    if skipfirst == 0 and skiplast == "" and os.path.isfile(filename+'.head') and os.path.isfile(filename+'.tail'):
        dothisall = False


    if dothisall:
        if skiplast == "":
            sys.exit("ERROR: skiplast not found")
        skiplast = 200 - int(skiplast) + 1
        print("skiplast2    :"+str(skiplast)+":",type(skiplast))
        if type(skiplast) == bool:
            sys.exit("ERROR: did not find skiplast")

        oben_chars = os.popen('head -'+str(skipfirst)+' '+filename+' | wc -c').read()  # writes log.lammps.head
        unten_chars = os.popen('tail -'+str(skiplast)+' '+filename+' | wc -c').read()  # writes log.lammps.head
        #print "oben_chars1  :",oben_chars
        #print "unten_chars1 :",unten_chars
        print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        count_oben=oben_chars.rstrip('\n')
        count_unten=unten_chars.rstrip('\n')
        print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        print("count_oben   :",count_oben)
        print("count_unten  :",count_unten)


        # write head (wird ueberschrieben)
        if skipfirst > 0:
            print("write ",filename,".head")
            os.system('head -'+str(skipfirst)+' '+filename+' > '+filename+'.head')  # writes log.lammps.head
        # write tail (wird ueberschrieben)
            print("write ",filename,".tail")
        os.system('tail -'+str(skiplast)+' '+filename+' > '+filename+'.tail')  # writes log.lammps.head



        # remove from the top (works great)
        # this was taken from: http://stackoverflow.com/questions/17330188/remove-first-n-lines-of-a-file-in-place-in-unix-command-line
        print("remove from the top ...")
        os.system('dd if="'+filename+'" bs="'+str(count_oben)+'" skip=1 of="'+filename+'" conv=notrunc')
        print("remove from the top truncate ...")
        os.system('truncate -s "-'+str(count_oben)+'" "'+filename+'"')

        # remove from the bottom
        # dd if=/dev/null of=log.lammps bs=1 seek=$(echo $(stat --format=%s log.lammps ) - $( tail -n1 log.lammps | wc -c) | bc ) #which removes from the bottom
        #os.system('dd if=/dev/null of='+filename+' bs=1 seek=$(echo $(stat --format=%s log.lammps ) - $( '+str(count_unten)') | bc )') #which removes from the bottom
        #for i in np.arange(skiplast):
        #    os.system('dd if=/dev/null of=log.lammps bs=1 seek=$(echo $(stat --format=%s log.lammps ) - $( tail -n1 log.lammps | wc -c) | bc )') #which removes from the bottom
        #
        # faster without the loop
        print("remove from the bottom ...")
        os.system('dd if=/dev/null of=log.lammps bs=1 seek=$(echo $(stat --format=%s log.lammps ) - $( tail -n'+str(skiplast)+' log.lammps | wc -c) | bc )') #which removes from the bottom

    print("get qpoints_read_in ...")
    qpoints_read_in = os.popen('tail -1 '+filename+'.head').read()
    qpoints_read_in = qpoints_read_in.rstrip()
    #print qpoints_read_in
    header = qpoints_read_in
    header = header.split(" ")

    print("nun muss gruperit werden um die space_fft_x_x_x zu speichern.",len(header))
    print(header)

    list_of_qpoints = []
    for i in header:
        if len(i.split("c_c")) == 2: # nur wenn es das gibt wird es interesant
            qpoint_string_in="_".join((i.split("c_c")[1]).split("_")[0:3])  # gives the ['10', '10', '20'] for every qpoint
            list_of_qpoints.append(qpoint_string_in)

    #print "list_of_qpoints:",list_of_qpoints
    fromm_too = []
    fromm_too_name = []
    for idx,i in enumerate(list_of_qpoints):
        #print idx,i
        if idx==0:
            current = i
            fromm = idx
            fromm_too_name.append(i)
            #print "from",idx
        if idx > 0:
            if i == current:
                continue
            else:
                current = i
                too = idx #-1
                #print "to",idx #-1
                fromm_too.append([fromm,too])
                fromm = idx
                fromm_too_name.append(i)
                #print "from",idx

    #print "to",len(list_of_qpoints) #-1
    fromm_too.append([fromm,len(list_of_qpoints)])

    #print
    #for idx,i in enumerate(list_of_qpoints):
    #    print idx,i
    #for idx,i in enumerate(fromm_too):
    #    print idx,i,fromm_too_name[idx]

    #for i in fromm_too:
    #    print "-->",i[0],i[1],"qp:",header[i[0]],"lll:",list_of_qpoints[i[0]]
    #    filename="space_fft_"+list_of_qpoints[i[0]]+".npy"
    #    filecontent=sum_all_new[:,i[0]:i[1]]
    #    np.savetxt(filename,filecontent)
    #print "lll:|",len(list_of_qpoints)
    #print "kkk:|",len(header)
    #print "???:|",header[98]
    for i in np.arange(len(list_of_qpoints),len(header)):
        #print i,header[i]
        fromm_too.append([i,i+1])
        fromm_too_name.append(header[i])
    #print
    #for idx,i in enumerate(fromm_too):
    #    print idx,i,fromm_too_name[idx]

    return fromm_too, fromm_too_name

# needs read csv
def get_space_fft_from_lammps_log(filename):
    ''' this reads directly from lammps.log; comments above and below the output should be no problem '''

    from_to, from_to_name = get_space_fft_prepare_lammps_log(filename)
    print()
    for idx,i in enumerate(from_to):
        print(idx,i,from_to_name[idx])
    print()
    print("pandas read csv ... c engine (15 GB log.lammps takes 29.15% of the cmmc's memory)")
    reader = read_csv(filename, sep='\s+',header=None,engine='c')
    print("save as matrix")
    reader = reader.as_matrix()
    print("save space_fft ... lines",reader.shape)
    lines = reader.shape[0]
    print("ll:",lines)
    print()
    for idx,i in enumerate(from_to):
        print(idx,i,from_to_name[idx])
        filename="space_fft_"+from_to_name[idx]+".npy"
        filecontent=np.transpose(reader[:,i[0]:i[1]])
        columns = filecontent.shape[0]
        print("cc:",columns,filecontent.shape)
        if columns > 1:
            print("columns > 1")
            sum_all_new = np.zeros((columns/2,lines),dtype=np.complex64)  # 4*10^8 seems to be max for cluster
            for j in np.arange(columns/2):
                sum_all_new[j] = filecontent[2*j] +1j * filecontent[2*j+1]  #out[:,0] - 1j * out[:,1]
                np.save(filename,sum_all_new)
        else:
            print("columns = 1")
            np.save(filename,filecontent)
            #sum_all_new[5,chunk*chunksize:chunk*chunksize+out.shape[0]] = out[:,0] - 1j * out[:,1]
        #np.save(filename,filecontent.astype(np.complex64))


    if False:
        sys.exit()
        #print sum_all_new,sum_all_new.shape
        #print
        #tempin = sum_all_new[:,-1]

        ### analysis of hov many datapots should be excluded
        ### -> REUSLT: NOT excluding the first 4 points results in an error of 0.4Kelvin
        ###
        ##@cutoffmax = 40
        ##@out = np.zeros((cutoffmax,3))

        ##@for cutoff in np.arange(cutoffmax):
        ##@    temp = np.copy(tempin[cutoff:])


        ##@    #np.savetxt("temp.dat",temp)
        ##@    temp_average = np.copy(temp)
        ##@    temp_std = np.copy(temp)
        ##@    temp_ste = np.copy(temp)
        ##@    print "Temperature average is",temp.mean(),"Kelvin"
        ##@    for idx,i in enumerate(temp):
        ##@        #print idx,i,temp[:(idx+1)].mean()
        ##@        temp_average[idx]=temp[:(idx+1)].mean()
        ##@        temp_std[idx]=temp[:(idx+1)].std()
        ##@        temp_ste[idx]=temp[:(idx+1)].std()/np.sqrt(idx+1)
        ##@        #if idx == 3:
        ##@        #    sys.exit()
        ##@    #np.savetxt("temp_average.dat",temp_average)
        ##@    #np.savetxt("temp_std.dat",temp_std)
        ##@    #np.savetxt("temp_ste.dat",temp_ste)
        ##@    print cutoff,temp_ste[-1]
        ##@    out[cutoff][0] = cutoff
        ##@    out[cutoff][1] = temp.mean()
        ##@    out[cutoff][2] = temp_ste[-1]
        ##@np.savetxt("out.dat",out)
        ##@print
        ##@print

        #cutoff=4
        #temp = np.copy(tempin[cutoff:])


        ##np.savetxt("temp.dat",temp)
        #temp_average = np.copy(temp)
        #temp_std = np.copy(temp)
        #temp_ste = np.copy(temp)
        #for idx,i in enumerate(temp):
        #    #print idx,i,temp[:(idx+1)].mean()
        #    temp_average[idx]=temp[:(idx+1)].mean()
        #    temp_std[idx]=temp[:(idx+1)].std()
        #    temp_ste[idx]=temp[:(idx+1)].std()/np.sqrt(idx+1)
        #    #if idx == 3:
        #    #    sys.exit()
        #np.savetxt("temp_average.dat",temp_average)
        #np.savetxt("temp_std.dat",temp_std)
        #np.savetxt("temp_ste.dat",temp_ste)
        #print "Temperature average is",temp.mean(),"Kelvin +-",temp_ste[-1]
        #print



        #sys.exit()
        ##sum_all_new = np.zeros((columns,2*10**8),dtype=np.complex64)  # 4*10^8 seems to be max for cluster
    return

def get_space_fft_from_lammps_log_old_but_in_chunks(filename, columns = 8, qpoint = False, has_temperature_column = True):
    ''' - this is for the [N,0,0] qpoint
        - the result (sum_all_new) can however be evaluated without knowledge of N
          if lammps could finc c_1 c_2 this could be automated (at least the skipfirst)
          can python also find the end (of a possibly huge file?)  make it fixed first
    '''
    print("starting lammppace_fft_from_lammps_log_old_but_in_chunks ...")
    if type(qpoint) == bool:
        sys.exit("ERROR qpoint is bool")
    skipfirst = None
    skiplast = None

    skipfirst = None
    skiplast = 27
    print("fi:",filename)
    from itertools import islice
    with open(filename) as myfile:
        head = list(islice(myfile, 2000))
        for idx,i in enumerate(head):
            #if "c_1 c_2" in i:
            if "Memory usage per processor =" in i:
                #print idx,i
                skipfirst = idx+1

    print("skipfirst",skipfirst)

    checkline = False
    # chck line i+1 to know how many elements
    if skipfirst != None:
        with open(filename) as myfile:
            head = list(islice(myfile, 200))
            for idx,i in enumerate(head):
                if idx == skipfirst+1:
                    checkline = i

    if type(checkline) != bool:
        #print "checkline:",checkline
        checkline = checkline.split(" ") #.remove("' '")
        checkline[-1] = checkline[-1].strip()  # remove '\n'
        checkline = [_f for _f in checkline if _f]
        #print "checkline:",checkline,
        elements = len(checkline)
        if has_temperature_column == True:
            elements = elements - 1
        # chck line i+1 to know how many elements



    print("read data from ",filename," at some point, if this file gets too large, a reading in chunks may become necessary")
    #print "with 108082514 lines this gives a memory Erroro --> chunk!"

    ###################################
    # new approach read in chunks
    ###################################
    print("read data from lammpsfiel with frequencies done!")
    #sum_all_new = np.zeros((columns/2,2*10**8),dtype=np.complex64)  # 4*10^8 seems to be max for cluster
    sum_all_new = np.zeros((columns,2*10**8),dtype=np.complex64)  # 4*10^8 seems to be max for cluster
    chunksize=100000
    chunksize=3;


    #reader = read csv(filename, sep=r"\s*",header=None,skiprows=skipfirst,chunksize=chunksize,error_bad_lines=False)
    #reader = read csv(filename, sep=r"\s*",header=None,skiprows=skipfirst,chunksize=chunksize)
    #reader = read csv(filename, sep=r"\s*",header=None,skiprows=skipfirst,skipfooter=30)

    # reader for log_for_sum_all
    print("columns:",columns)
    print("columns/2:",columns/2,np.arange(columns/2))
    print("sum_all_new.shape",sum_all_new.shape)
    #sys.exit()
    for chunk,pos in enumerate(reader):
        #print "chunk    :",chunk
        #print "pos      :",pos
        #print "pos.shape:",pos.shape
        #print "pos.as_matrix():"#,pos.as_matrix()
        #print "pos"


        #!#sys.exit()
        out = pos.as_matrix()
        #print "read line nr:",chunk*chunksize

        #print "kkk:",out
        #print "out.shape:",out.shape,"out[:,0]:",out[:,0]
        #print "chunk:",chunk,"chunksize:",chunksize,"out.shape:",out.shape,"out.shape[0]:",out.shape[0],"|||",chunk*chunksize,chunk*chunksize+out.shape[0]
        #
        #
        # sum_all[ind_qs,mdstepbegin : mdstepend] = np.sum(np.exp(np.sum((xxx*qs_all[ind_qs]),axis=2)),axis=1)
        # np.exp(x+iy) wird in lammps eingesetzt wobei wir aber ueber alle positions gehen.
        # cos(x/${alat}/${N}*2.0*${N}*PI)        == realteil
        # sin(x/${alat}/${N}*2.0*${N}*PI)        == imaginaerteil
        #       wobei /${alat}/${N} das scaling (in direkte coordinaten) ist
        #       *2.0*${N}*PI ist das equivalent_qs (bekommt nur N und macht daraus die qpunkte)
        #

        # ths should be the general case probably
        for idx in np.arange(columns/2):  #  array([0, 1, 2, 3]) for dom and  array([0, 1, 2]) for fcc [N 0 0]
            #print "idx:",idx
            sum_all_new[idx,chunk*chunksize:chunk*chunksize+out.shape[0]]         = out[:,idx*2] + 1j * out[:,idx*2+1]
            #print "idx+col:",idx+columns/2
            sum_all_new[idx+columns/2,chunk*chunksize:chunk*chunksize+out.shape[0]] = out[:,idx*2] - 1j * out[:,idx*2+1]

        # the next 3 lines word for fcc
        #sum_all_new[0,chunk*chunksize:chunk*chunksize+out.shape[0]] = out[:,0] + 1j * out[:,1]
        #sum_all_new[1,chunk*chunksize:chunk*chunksize+out.shape[0]] = out[:,2] + 1j * out[:,3]
        #sum_all_new[2,chunk*chunksize:chunk*chunksize+out.shape[0]] = out[:,4] + 1j * out[:,5]
        ## (for dominique) sum_all_new[3,chunk*chunksize:chunk*chunksize+out.shape[0]] = out[:,6] + 1j * out[:,7]

        #####################################################################################################
        ### die naechsten 3 sind eigentlich redundant und sollten eigentlich
        ### erst spaeter bei der berechnung des powerspektrums
        ### beruecksichtigt werden.
        ### DIS GILT FUER DIE LONGITUDINALEN, JEDOCH AUCH FUER DIE TRANSVERSALEN???
        #####################################################################################################
        #sum_all_new[5,chunk*chunksize:chunk*chunksize+out.shape[0]] = out[:,0] - 1j * out[:,1]
        #sum_all_new[4,chunk*chunksize:chunk*chunksize+out.shape[0]] = out[:,2] - 1j * out[:,3]
        #sum_all_new[3,chunk*chunksize:chunk*chunksize+out.shape[0]] = out[:,4] - 1j * out[:,5]

    #print "chunklast!:",chunk,"out.shape[0]",out.shape[0],"--------->",chunk*chunksize+out.shape[0]
    sum_all_new = sum_all_new[:,:(chunk*chunksize+out.shape[0])]
    #print "sum_all_new:"
    #print sum_all_new
    #sys.exit()

    if abs(sum_all_new[:,0].max()) > 10**-13:
        print("first line:",sum_all_new[:,0])
        #sys.exit("first line seems to be too high")
    #qs1 = data[...,0] + 1j * data[...,1]
    #qs2 = data[...,2] + 1j * data[...,3]
    #qs3 = data[...,4] + 1j * data[...,5]
    #qs4 = data[...,4] + 1j * -data[...,5]
    #qs5 = data[...,2] + 1j * -data[...,3]
    #qs6 = data[...,0] + 1j * -data[...,1]

    #print "data.shape:",data.shape
    #print data
    #print "-----------"
    ##sum_all_new = np.zeros((6,data.shape[0]),dtype=np.complex64)
    #i=0
    #for vorzeichen in np.array([+1.0,-1.0]):
    #    for qspair in np.array([[0,1],[2,3],[4,5]]):
    #        #print "vz:",vorzeichen,"qspair:",qspair
    #        #qs = data[...,qspair[0]] + 1j * vorzeichen * data[...,qspair[1]]
    #        #qs = 1j*Data[...,1]; result += Data[...,0]
    #        qs = 1j * vorzeichen * data[...,qspair[1]];qs += data[...,qspair[0]]  # saves memory, no intermediate result
    #        print "qs:",qs,qs.shape,i
    #        sum_all_new[i] = qs
    #        i=i+1
    #print "sum_all_new[0]:",sum_all_new[0]
    #np.save("sum_all_new_zuvor",sum_all_new)
    #sum_all_new = sum_all_new[~np.isnan(sum_all_new).any(axis=1)]
    #print "sum_all_new[0]:",sum_all_new[0]

    print("np.save(sum_all_new)")
    qpointstr = qpointstring(qpoint)
    #np.save("sum_all_new_"+str(qpoint[0])+str(qpoint[1])+str(qpoint[2]),sum_all)
    #np.save("sum_all_new_"+qpointstr,sum_all_new)
    np.save("space_fft_"+qpointstr,sum_all_new)

    #np.save("sum_all_new",sum_all_new)
    return sum_all_new

def get_space_fft_from_xaa_xab_folder(qpstring):
    ''' qpstring is a list of all qpoints which are to be evaluated
    e.b. [ '1_0_0', '2_0_0', ...] '''
    #print "qpstring:",qpstring
    #print "len(qpstring):",len(qpstring),type(qpstring)
    if type(qpstring) == str:
        qpstring = [ qpstring ]
    #def sum_them_all(j):  # where j is one the qpstirng
    foldername="x[a-z][a-z]_"
    filenames_vor_qp="sum_all_new_"
    filenames_nach_qp=".dat"

    folder=sorted(glob.glob(foldername))

    files=sorted(glob.glob(folder[0]+"/"+filenames_vor_qp+"*"+filenames_nach_qp))
    #qpstring= []
    ##print "xx:",qpstring
    #print "----------------folder-------------------"
    #print "len(folder):",len(folder),qpstring

    #for num,i in enumerate(folder):
    #    print num,i
    #print
    #print "----------------files for folder[0]------"
    #for i in files:
    #    print i
    #for i in files:
    #    qpstring.append(i.split(filenames_vor_qp)[1].split(filenames_nach_qp)[0])
    #print qpstring

    #######################################################
    ## start filtering for specific branch
    #######################################################
    ## currently only evaluate [1 1 1] L branch
    #qpstringsubset = []
    #for i in qpstring:
    #    #print "possible:",i
    #    alle = i.split("_")
    #    #print alle,alle[0]==alle[1]==alle[2]

    #    # 1 1 1 (L)
    #    #if alle[0]==alle[1]==alle[2]:
    #    #    qpstringsubset.append(i)

    #    # 1 1 0 (L)
    #    #if alle[0]==alle[1] and alle[2] == '0':
    #    #    qpstringsubset.append(i)

    #    ## 1 1 0 (T1)
    #    #if alle[0]==alle[1] and alle[2] == '20':
    #    #    qpstringsubset.append(i)


    #    ## 1 0 0 (T)
    #    #if alle[0]=='20' and (int(alle[1]) <= 10 and int(alle[1]) >=1) and alle[2]=='0':
    #    #    qpstringsubset.append(i)

    #    #if alle[0]==alle[2]=='20':  # 1 0 0 Tnew
    #    #    qpstringsubset.append(i)
    #    print "alle:",alle[0],alle[1],alle[2],i,"==============",int(alle[0])==20+int(i)
    #    if int(alle[0])==20+int(i) and int(alle[1])==20 -int(i) and alle[2]=='0':
    #        qpstringsubset.append(i)


    #print
    #print "subset:",qpstringsubset
    #qpstring=qpstringsubset
    #print
    #print "out:",qpstring
    #print
    #for idx,i in enumerate(qpstring):
    #    print "in--->>>:",idx,i

    ######################################################
    # start filtering for specific branch ... DONE
    ######################################################

    #theList= range(6)
    #N=2
    #print theList
    #subList = [theList[n:n+N] for n in range(0, len(theList), N)]
    #print subList
    #
    ##print grouper(2,np.arange(6))


    #from multiprocessing import Process
    #def doit(j):
    #    print "all stuff to do with j",j
    #return
    for idxj,j in enumerate(qpstring):  # this needs to be parllelized
        #print "j:",j
        #sys.exit()
    #def sum_them_all(j):
        if os.path.isfile('space_fft_'+j+'.npy'):
            print('space_fft_'+j,"does already exist...continue")
            continue
        else:
            print('space_fft_'+j+'.npy',"does not exist --> creating it!")


        #print "idxj:",idxj,"j:",j
        #if j != "1_0_0":
        #    continue
        for idxi,i in enumerate(folder):
            #strout="idxj:",idxj,"j:",j,"         idxi:",idxi,"i:",i
            file=i+"/"+filenames_vor_qp+j+filenames_nach_qp
            #print "importing "+file
            data=np.loadtxt(file)
            #print file,"data.shape:",data.shape,len(data.shape)
            if len(data.shape) == 1:
                print(data,type(data))
                data = np.array([data])
                #print "--.>",data,type(data)
            #print file,"data.shape:",data.shape,data.shape[0],data.shape[1]
            #print "importing "+file+" done!"

            shape = data.shape
            lines=shape[0]
            columns=shape[1]
            theList=list(range(columns))
            N=2
            groupby = [theList[n:n+N] for n in range(0, len(theList), N)]
            sum_all_curr = np.zeros((columns/2,lines),dtype=np.complex64)
            if idxi==0: #  and idxj==0:
                print("IDXI +++ 0000000000000000000000")
                sum_all_new = np.zeros((columns/2,lines),dtype=np.complex64)
                #print "shape:",shape
                #print "groupby:",groupby

            for idk,k in enumerate(groupby):
                #print idk,k
                #print data[:,k[0]],data[:,k[0]].shape
                #print data[:,k[1]],data[:,k[1]].shape
                #print sum_all_curr[:,idk].shape
                #print
                #print
                sum_all_curr[idk].real = data[:,k[0]]
                sum_all_curr[idk].imag = data[:,k[1]]
            if idxi==0: # and idxj==0:
                print("idxj:",idxj,"j:",j,"idxi:",idxi,"i:",i,file+" b aaa fft.shape:",sum_all_new.shape,"curr.shape:",sum_all_curr.shape)
                sum_all_new = sum_all_curr
                print("idxj:",idxj,"j:",j,"idxi:",idxi,"i:",i,file+" a aaa fft.shape:",sum_all_new.shape,"curr.shape:",sum_all_curr.shape)
            else:
                print("idxj:",idxj,"j:",j,"idxi:",idxi,"i:",i,file+" b bbb fft.shape:",sum_all_new.shape,"curr.shape:",sum_all_curr.shape)
                sum_all_new = np.hstack((sum_all_new,sum_all_curr))
                print("idxj:",idxj,"j:",j,"idxi:",idxi,"i:",i,file+" a bbb fft.shape:",sum_all_new.shape,"curr.shape:",sum_all_curr.shape)
                #print "saa:",sum_all_new
            print()
        np.save("space_fft_"+j,sum_all_new)
    return


##########################################################################################
# functions to create space_fft from lammps outputfiles ##################################
##########################################################################################
def lammps_check_if_trj_lammps_and_trj_lammpsnew_equal(f1,f2,rm=False):
    if f1 == f2:
        print("both files equal therefore return")
        return True
    c1 = os.popen('tail -1 '+f1+" | awk '{print $NF}'").read().rstrip().rstrip()
    c2 = os.popen('tail -1 '+f2+" | awk '{print $NF}'").read().rstrip().rstrip()
    c3 = os.popen('tail -1 '+f1+" | awk '{print $(NF-1)}'").read().rstrip().rstrip()
    c4 = os.popen('tail -1 '+f2+" | awk '{print $(NF-1)}'").read().rstrip().rstrip()
    print("c1,c2:",c1,c2)
    print("c3,c4:",c3,c4)
    if c1 == c2 and c3==c4:
        print("c1 == c2 and c3 == c4")
        if type(rm) != bool:
            if rm == f1 or rm == f2:
                if rm == f1:
                    print("removing:",rm,"since same content in",f2)
                if rm == f2:
                    print("removing:",rm,"since same content in",f1)
                os.remove(rm)
        return True

    else:
        print(f1,"and",f2,"are not equal!")
        return False

def lammps_check_if_trj_lammps_and_trj_lammpsnew_equal_also_check_head(f1,f2,rm=False):
    if f1 == f2:
        print("both files equal therefore return")
        return True
    c1 = os.popen('tail -1 '+f1+" | awk '{print $NF}'").read().rstrip().rstrip()
    c2 = os.popen('tail -1 '+f2+" | awk '{print $NF}'").read().rstrip().rstrip()
    c3 = os.popen('tail -1 '+f1+" | awk '{print $(NF-1)}'").read().rstrip().rstrip()
    c4 = os.popen('tail -1 '+f2+" | awk '{print $(NF-1)}'").read().rstrip().rstrip()

    c5 = os.popen('head -1 '+f1+" | awk '{print $(NF-1)}'").read().rstrip().rstrip()
    c6 = os.popen('head -1 '+f2+" | awk '{print $(NF-1)}'").read().rstrip().rstrip()
    c7 = os.popen('head -1 '+f1+" | awk '{print $(NF)}'").read().rstrip().rstrip()
    c8 = os.popen('head -1 '+f2+" | awk '{print $(NF)}'").read().rstrip().rstrip()
    print("c1,c2:",c1,c2)
    print("c3,c4:",c3,c4)
    if c1 == c2 and c3==c4 and c5 == c6 and c7 == c8:
        print("c1 == c2 and c3 == c4 and c5 == c6 and c7 == c8")
        if type(rm) != bool:
            if rm == f1 or rm == f2:
                if rm == f1:
                    print("removing:",rm,"since same content in",f2)
                if rm == f2:
                    print("removing:",rm,"since same content in",f1)
                os.remove(rm)
        return True

    else:
        print(f1,"and",f2,"are not equal!")
        return False

def lammps_split_lammps_trajectory_in_xaa_files(filenamelammps="trj_lammps.out",filenamepos="trj_lammpsnew.out",positionsfilename = "positions.*",infilefilename = "in_file_dynamics.in",linesokformemory = 800000000,qpoints_all=False,args=False,verbose=False):
    ''' splits lammps trajectory file in xaa, xab, xac ...'''
    print(bcolors.OKGREEN + "lammps_split_lammps_trajectory_in_xaa_files" +bcolors.ENDC)
    print()
    print()
    print("das splitten und so weiter sollte alles auf einen knoten submittet werden!")
    print("daher: exit after grep und split!")
    print()
    print()
    positionsfile = glob.glob(positionsfilename)
    if len(positionsfile) != 1:
        print("positionsfile:",positionsfile)
        sys.exit("Not one positionsfile but "+str(len(positionsfile)))

    print("positionsfile:",positionsfile)

    infile = glob.glob(infilefilename)
    print("infile:",infile)
    if len(infile) != 1:
        print("infile:",infile)
        sys.exit("Not one infile but "+str(len(infile)))

    atoms1 = os.popen('wc -l '+positionsfile[0] + "| awk '{print $1}'").read()
    atoms = int(atoms1)-8
    atoms2 = int(os.popen('head -n 2 '+positionsfile[0] + "| tail -1 | awk '{print $1}'").read().rstrip())
    print("atoms2:",atoms2)
    if atoms != atoms2:
        print("atoms:",atoms,"atoms2:",atoms2)
        sys.exit("atoms not atoms2")

    supercelllength = float(os.popen('head -n 4 '+positionsfile[0] + "| tail -1 | awk '{print $2}'").read().rstrip())
    print("supercelllength:",supercelllength)
    print("structure:",args.structure,"usestruct:",args.usestruct)
    sc = N = round(float((atoms/args.usestruct)**(1/3.)),0)
    alat = supercelllength/N
    print("alat:",alat)
    stepslammps = int(os.popen('grep "^run " '+infile[0] + "| awk '{print $2}'").read().rstrip())
    dump = int(os.popen('grep "^dump dump1" '+infile[0] + "| awk '{print $5}'").read().rstrip())

    print("atoms:",atoms,type(atoms))
    print("sc:",sc)
    print("stepslammps:",stepslammps)
    print("dump:",dump)
    stepswritten = stepslammps/dump+1
    print("stepswritten    :",stepswritten)
    lineswritten = stepswritten*atoms
    print("lineswritten    :",lineswritten,"in trj_lammpsnew.out")
    print("linesokformemory:",linesokformemory)
    nsteps = linesokformemory / atoms
    nsteps = int(round(linesokformemory / atoms/10,0)*10)
    print("nsteps:",nsteps,"per xaa file")
    from string import ascii_lowercase
    stepsremain = stepswritten

    print("######################################################################################")
    print("# grep (now sed)")
    print("# it would be best to submit this to the que on one core (grep and sed)")
    print("# it is assumed thtat this was already done by the que, if not, comment in")
    print("######################################################################################")
    if os.path.isfile(filenamelammps):
        print(filenamelammps,"exists!")
        if not os.path.isfile(filenamepos):
            print(filenamepos,"does not exist! --> comment in the next lines!!")
            #print "here we have to ensure that we are in a screen session"
            #import socket
            #hostname=(socket.gethostname())
            #if hostname == 'cmmc001' or hostname == 'cmmc002':
            #    checkscreen = os.popen('echo $STY').read().rstrip()
            #    if checkscreen == "":
            #        onscreen = False
            #        print "onscreen:",onscreen
            #        sys.exit("ERROR: YOU WANT TO GREP A POSSIBLY LARGE FILE and should do this in a screen session!")
            #    else:
            #        onscreen = True
            #        print "onscreen:",onscreen,os.getcwd()
            #        print "10 GB take approx. 2min 20 sec ..."
            #        print "100 GB take approx. 1400 sec = 23 min ..."
            #        #os.popen('grep "^1 " '+filenamelammps + "| awk '{print $2,$3,$4}' > trj_lammpsnew.out").read()
            #        print "hier nun der submit fuer sed und split welcher dann wieder lammps_pos_to_sum startet und am besten gleich noch die c++ jobs startet"
            #        print "sed und split brauchen nur die nsteps, das ist einfach"
            #        print "hier sollte nur ein ping gegeben werden wenn sed/split/ gebraucht wird da sonst alles nicht da"
            #        os.popen("sed -n 's/^1 //p' "+filenamelammps + " > trj_lammpsnew.out").read()
            #        if os.path.isfile(filenamelammps) and os.path.isfile(filenamepos):
            #            lammps_check_if_trj_lammps_and_trj_lammpsnew_equal(filenamelammps,filenamepos,rm=filenamelammps)
    else:
        print(filenamelammps,"does not exist anymore, presumably grep already done!")

    if os.path.isfile(filenamelammps) and os.path.isfile(filenamepos):
        lammps_check_if_trj_lammps_and_trj_lammpsnew_equal(filenamelammps,filenamepos,rm=filenamelammps)

    print()
    print()
    splitanzahl=-1
    print("      stepstot","\t","steps  ","\t","stepsremain")
    print("---------------------------------------------------")
    lastfile = False
    lastamoutoflines = False
    xaa_filenames = []
    xaa_steps = []
    xaa_lines = []

    # make here an array
    idlast=0
    for idx,c in enumerate(ascii_lowercase):
        for idy,d in enumerate(ascii_lowercase):
            nstepscurr = nsteps
            #print "(stepsremain - nsteps):",(stepsremain - nsteps)
            if (stepsremain - nsteps) > 0:
                stepsremain = stepsremain - nsteps
                exittrue = False
            else:
                nstepscurr = stepsremain
                stepsremain = 0
                exittrue = True
            splitanzahl = splitanzahl+1
            lines = nstepscurr*atoms
            print(idlast,"x"+c+d,stepswritten,"\t",nstepscurr,"\t",stepsremain,splitanzahl,lines,exittrue)
            xaa_filenames.append("x"+c+d)
            xaa_steps.append(nstepscurr)
            xaa_lines.append(lines)
            lastfile = "x"+c+d
            lastamoutoflines = lines
            idlast = idlast + 1
            if exittrue: break
        if exittrue: break

    #for i in np.arange(20)+1:
    #    atoms = 4*i**3
    #    print  i,atoms,linesokformemory / atoms, int(round(linesokformemory / atoms/10,0)*10)

    if splitanzahl == 0:
        if os.path.isfile(filenamepos) and not os.path.isfile('xaa'):
            os.rename(filenamepos, "xaa")
        if os.path.isfile(filenamepos) and os.path.isfile('xaa'):
            print("move one to other if those are equal")
            lammps_check_if_trj_lammps_and_trj_lammpsnew_equal_also_check_head(filenamepos,'xaa',rm=filenamepos)
        if not os.path.isfile(filenamepos) and os.path.isfile('xaa'):
            pass # this is also ok
        if not os.path.isfile(filenamepos) and not os.path.isfile('xaa'):
            sys.exit("where are the positions?")
        xaa_filenames = [ filenamepos ]
        xaa_filenames = [ 'xaa' ]
        xaa_steps = [ stepswritten ]
        xaa_lines = [ lineswritten ]

    print("splitanzahl:",splitanzahl)
    print("xaa_filenames:",xaa_filenames)
    print("xaa_steps:",xaa_steps)
    print("xaa_lines:",xaa_lines)
    print()
    print()
    print()
    print("######################################################################################")
    print("# split")
    print("######################################################################################")
    if splitanzahl > 0:
        print("splitanzhal > 0, in fact it is:",splitanzahl)
        if os.path.isfile("xaa"):
            print("xaa file already exists, splitting presumable already started/done!")
        else:
            print("splitting by:",str(int(nsteps*atoms)),"since xaa does not exist! --> splitting")
            print("hmmm ...unexpected... this should have been done by the cluster")

        #    checkscreen = os.popen('echo $STY').read().rstrip()
        #    if checkscreen == "":
        #        onscreen = False
        #        print "onscreen:",onscreen
        #        #sys.exit("ERROR: YOU WANT TO SPLIT A POSSIBLY LARGE FILE and should do this in a screen session!")
        #    else:
        #        onscreen = True
        #        print "onscreen:",onscreen
        #    os.popen('split -l '+str(int(nsteps*atoms))+" "+filenamepos).read()

    print("lastfile:",lastfile)
    print("lastamoutoflines:",lastamoutoflines)
    if splitanzahl > 0 and os.path.isfile("xaa") and os.path.isfile(lastfile) and os.path.isfile(filenamepos):
        lammps_check_if_trj_lammps_and_trj_lammpsnew_equal(filenamepos,lastfile,rm=filenamepos)

    print("######################################################################################")
    print("now check if last file .."+xaa_filenames[-1]+".. has the correct number of lines:")
    print("######################################################################################")

    if os.path.isfile(xaa_filenames[-1]+"_wcl"):
        print(xaa_filenames[-1]+"_wcl","exists, therefore correct number of lines!")

    if not os.path.isfile(xaa_filenames[-1]+"_wcl") and os.path.isfile(xaa_filenames[-1]):
        print("it seems that,",xaa_filenames[-1]+"_wcl","does not exist and ",xaa_filenames[-1],"exists!")
        wclast = int(os.popen('wc -l '+xaa_filenames[-1] + "| awk '{print $1}'").read().rstrip())
        print("wclast:",wclast,type(wclast))
        print("xaa_lines[-1]:",xaa_lines[-1],type(xaa_lines[-1]))
        if wclast != xaa_lines[-1]:
            target = open(xaa_filenames[-1]+"_wclerror", 'w')
            target.write("wc -l gives "+str(wclast)+" lines\n")
            target.close()
            sys.exit("wclast not last xaa_lines")
        else:
            print("writing that number of lines is ok",xaa_filenames[-1]+"_wcl")
            target = open(xaa_filenames[-1]+"_wcl", 'w')
            target.write("wc -l gives "+str(wclast)+" lines\n")
            target.close()
            print("######################################################################################")
            print("seems to have correct number of lines!")
            print("######################################################################################")

    if not os.path.isfile(xaa_filenames[-1]+"_wcl") and not os.path.isfile(xaa_filenames[-1]):
        print("######################################################################################")
        print("the last file, ",lastfile," does not exist! I will therefore exit here")
        print("######################################################################################")
        sys.exit()

    print("######################################################################################")
    print("now check if xaa_ folder exist, if not, create them and send to cluster or not")
    print("######################################################################################")
    xaa_folder= glob.glob("xa?_")
    print("len(xaa_folder):",len(xaa_folder))
    if len(xaa_folder) == 0 and check_for_all_space_fft_and_rm_xaa_files_folders(sc,args) == False:
        print("##################################################################################")
        print("create xaa_ folders for c++ skript (even one folder should be submitted to que)")
        print("##################################################################################")
        print("all xaa_filenames:",xaa_filenames)
        avail_xaa_filenames = []
        avail_xaa_steps = []
        for idxx,op in enumerate(xaa_filenames):
            #print "xaa_file:",op
            if os.path.isfile(op) == True:
                avail_xaa_filenames.append(xaa_filenames[idxx])
                avail_xaa_steps.append(xaa_steps[idxx])
        print("avail_xaa_filenames:",avail_xaa_filenames,len(avail_xaa_filenames))
        if len(avail_xaa_filenames) > 0:
            #print "I am here!!!"
            #sys.exit()
            submit_c_skript_parallel(avail_xaa_filenames,avail_xaa_steps,sc,alat,atoms,executable="/cmmc/u/aglen/sascha/source/test",args=args)
        else:
            print("##################################################################################")
            print("no xa files available!")
            print("##################################################################################")
    else:
        print("xaa folder exist already.")
    return xaa_filenames,xaa_steps

def submit_c_skript_parallel(xaa_filenames,xaa_steps,sc,alat,atoms,executable="/cmmc/u/aglen/sascha/source/test",qsub="/u/aglen/submit.lammps.lifetimejob.c_sascha.cmfe40.sh",args=False):
    print()
    print("##############################################################################")
    print("now applying C++ skript (will do this even for one xaa folder)")
    print("this should also create a joblist to submit s++ skript to cluster or directly")
    print("do so if necessary")
    print("this should be started in parallel in any case! -> submit to cluster")
    print("##############################################################################")
    qpoints_all = get_all_qpoints(args.qvec,args)

    hier = os.getcwd()
    for idx,filename in enumerate(xaa_filenames):
        print(idx,filename)
        os.chdir(hier)
        lammps_create_c_SETTINGS_file(infile=filename,atoms=atoms,steps=xaa_steps[idx],scale=alat*sc,qpoints_all=qpoints_all,executable=executable)
        os.chdir(filename+"_")
        print("cwd:",os.getcwd())
        print("here execture qsub")
        print(os.popen('qsub '+qsub).read().rstrip())   # this print statement submits the job
    return

def lammps_create_c_SETTINGS_file(filename="SETTINGS",infile=False,atoms=False,steps=False,scale=False,qpoints_all=False,executable=False):
    ''' infile: xaa'''
    if type(infile) == bool:
        sys.exit("specify infile")
    if type(atoms) == bool:
        sys.exit("specify atoms")
    if type(scale) == bool:
        sys.exit("specify scale")
    if type(infile) == bool:
        sys.exit("specify infile")
    if type(qpoints_all) == bool:
        sys.exit("specify qpoints_all")
    if type(executable) == bool:
        sys.exit("specify executable")

    if not os.path.isfile(infile): sys.exit("file "+infile+" does not exist!")
    if not os.path.isfile(executable): sys.exit("file "+executable+" does not exist!")

    if not os.path.exists(infile+"_"): os.makedirs(infile+"_")

    settingsfile = infile+"_"+"/SETTINGS"
    if os.path.isfile(settingsfile):
        os.remove(settingsfile)

    target = open(settingsfile, 'w')
    target.write("NATOMS = "+str(atoms)+"\n")
    target.write("NSTEPS = "+str(steps)+"\n")
    for idx,i in enumerate(qpoints_all):
        target.write("Q"+str(idx+1)+" = "+str(i[0])+" "+str(i[1])+" "+str(i[2])+"\n")
    target.write("SCALE = "+str(scale)+"\n")
    target.write("INFILE = ../"+infile+"\n")
    target.close()
    import shutil
    shutil.copy2(executable, infile+"_/test")

#############
# UNUSED
#############
def get_space_fft_from_positions_fast_read_at_once_UNUSED(filename, qpoints_all = False, scale_data = 1., N=False, usestruct=False,args=args):
    ''' dont use engine='python' in read csv which is 6 times slower,
        - this module is aimed for speed, therefore loads the positions as whole into memory
        - chunking however should also work with c engine
        - when checked, cmmc001 had 56GB of free ram
        - 10x10x10sc (4000atoms) with 10^6steps (every 40 wirtten) did not fit into memory
          also it should, I guess it was close ;)

        then, go qpoint by qpoint and do the calculation(s)
        when
    '''
    print("#################################### try_lammpsnew.out #########################")
    #filename = 'trj_lammpsnew.outfirst3'
    #filename = 'trj_lammpsnew.outfirst300'     # 0.86 vs 2.2 sec
    #filename = 'trj_lammpsnew.outfirst3000'     # 3.7 vs 22.0 sec (head -12000000) = 12*10^6 lines of 10*10^8
                                                # size in memory: 0.288 GB
    #filename = 'trj_lammpsnew.outfirst30000'   # 36.7sec; size in memory: 2.88GB (1/10 of whole file)
    #filename = 'trj_lammpsnew.out'             # lines: 1000004000; 23GB; so should be 12GB in memory;
    filename = 'trj_lammpsnew.out'

    import time
    start1 = time.time()
    start = time.time()
    #df = read csv(filename, sep=' ', header=None,engine='c',) # 3.7 sec
    df = read_csv(filename, sep=' ', header=None,engine='c',chunksize=4000*2) # 3.7 sec
    for chunk,pos in enumerate(df):
        print("chunk:",chunk)
        xx = pos/float(scale_data)
        print("pos  :",xx)
        print(xx.as_matrix())
        sys.exit()
    #df = read csv(filename, sep=' ', header=None,engine='c', dtype={'0': np.float32, '1': np.float32, '2': np.float32}) # 3.7 sec
    print(sys.getsizeof(df)*1e-09,"GB", "df == pandas read in")#,df.nbytes
    print(df[0].dtype) # float64
    end = time.time()
    print("#1) read in data :",bcolors.FAIL + str((end - start)) +" sec"+ bcolors.ENDC)

    start = time.time()
    df = df/float(scale_data)
    print(sys.getsizeof(df)*1e-09,"GB", "df == pandas")
    df = df.as_matrix()
    print(sys.getsizeof(df)*1e-09,"GB", "df == nupy",df.nbytes*1e-09,"in GB")
    z = np.copy(df)
    print(sys.getsizeof(z)*1e-09,"GB", "df == numpy copy",z.nbytes*1e-09,"in GB")
    #z = z.astype(complex64)
    z = np.complex64(z)
    print(sys.getsizeof(z)*1e-09,"GB", "df == numpy copy",z.nbytes*1e-09,"in GB")



    qs_all = equivalent_qs([1,1,1])
    qs_all = equivalent_qs([1,0,0])
    print("00000000000:",df.shape)
    xxx = df*qs_all[0]
    print(sys.getsizeof(xxx)*1e-09,"GB")
    dfqs = np.exp(np.sum(df*qs_all[0],axis=1)).reshape(-1, N**3*args.usestruct).sum(axis=1) # the reshape is quick, it is probably the np.exp which takes time
    print(sys.getsizeof(dfqs)*1e-09,"GB")
    print("------>",dfqs.dtype)
    #aaa[4000:8000].sum()
    #aaa[8000:12000].sum()

    end = time.time()
    print("#2) np.exp(np.sum()):",bcolors.FAIL + str((end - start)) +" sec"+ bcolors.ENDC)
    start = time.time()
    #dfqs = dfqs.reshape(-1, N**3*args.usestruct).sum(axis=1)
    end = time.time()
    print("#3) reshape.sum:",bcolors.FAIL + str((end - start)) +" sec"+ bcolors.ENDC)
    #aaa = np.sum(np.exp(np.sum(xxx*qs_all[0],axis=1)))
    print()
    print(dfqs[:10])
    print(dfqs.shape)
    end1 = time.time()
    print("#5) total:",bcolors.HEADER+  str((end1 - start1)) +" sec"+ bcolors.ENDC)
    return dfqs

# needs read csv
def get_space_fft_from_hd5f_UNUSED(filename='dump_h5md.h5',qpoint=[1,0,0],scale_data = False):
    print("#################################### hd5f #########################")
    import time
    start1 = time.time()

    filename = 'dump_h5md.h5'
    start = time.time()
    f = h5py.File(filename, 'r')
    positions = f['particles/all/position']
    end = time.time()
    print("#1) hdf5 read:",bcolors.FAIL + str((end - start)) +" sec"+ bcolors.ENDC)
    start = time.time()
    print(positions['value'].shape)
    positions = positions['value'][:]/float(scale_data)
    print("----->",positions.shape)
    end = time.time()
    print("#2) tonumpy:",bcolors.FAIL + str((end - start)) +" sec"+ bcolors.ENDC)
    print("positions.shape:",positions.shape)
    print(sys.getsizeof(positions)*1e-09,"GB", "df == numpy copy",positions.nbytes*1e-09,"in GB")


    start = time.time()
    out = np.zeros(positions.shape[0])
    end = time.time()
    print("#3) init empty out:",bcolors.FAIL + str((end - start)) +" sec"+ bcolors.ENDC)
    print("sh:",out.shape)
    qs_all = equivalent_qs(qpoint)
    print("qs_all[]:",qs_all[0])
    start = time.time()
    for step,pos in enumerate(positions):
        #xxx = positions['value'][step]/float(scale_data)
        #out[step] = np.sum(np.exp(np.sum(xxx*qs_all[0],axis=1)))
        out[step] = np.sum(np.exp(np.sum(pos*qs_all[0],axis=1)))
    end = time.time()
    print("#4) exp:",bcolors.FAIL + str((end - start)) +" sec"+ bcolors.ENDC)
    print(out.shape)
    print(out[:10])
    end1 = time.time()
    print("#5) total:",bcolors.HEADER+  str((end1 - start1)) +" sec"+ bcolors.ENDC)



    return


##########################################################################################
# functions related to creation of the powerspectrum #####################################
##########################################################################################
def get_mdstepstocalc(space_fft=False): #, args = False):
    ''' just determines which mdsteps you want to run over
    0) if args.mdstepstocalc_all is defined return that NO!!! wait until 3) if args.mdsteps
    0) if args.mdstepstocalc_all is in sql db get it
    1) if args.mdstepstocalc_all not defined in db sql to that
    2) if args.mdstepstocalc_all_ps not defined to that
    3) return args.mdsteps if necessary
    4) return args.mdstepstocalc_last if necessary
    '''
    #print "aaaa",args.mdstepstocalc_all,type(args.mdstepstocalc_all)
    #if type(args.mdstepstocalc_all) == np.ndarray:  DONT DO THAT if args.mdsteps are definec
    #    return args.mdstepstocalc_all

    #if type(args.mdstepstocalc_all) != np.ndarray:
    # I DONE USE NOW SQL
    #    args.mdstepstocalc_all = get_db_mdstepsall(args.db)
    #if type(args.mdsteps) != bool:
    #    return args.mdsteps
    #print "bbba",args.mdstepstocalc_all,type(args.mdstepstocalc_all)
    #######################################################################
    # make loading before calculating since it does not need the space_fft
    #######################################################################
    global args
    #print 'args ? dt',args.dt   # hier kommt immer das richtige dt raus
    if type(args.mdstepstocalc_all) != np.ndarray and os.path.isfile("sql.mdsteps") == True:
        if args.verbose > 1:
            print('loading from file .... since',type(args.mdstepstocalc_all))
        args.mdstepstocalc_all = np.loadtxt("sql.mdsteps")
        #print 'after loading .... type:',type(args.mdstepstocalc_all)


    faktor=args.dt*0.001
    if type(args.mdstepstocalc_all) != np.ndarray and type(space_fft) != bool:
        print("calculating ...")
        #if args.verbose > 2:
            #print "making mdstepstocalc ..."
        #####################################################
        # MDSTEPSTOCALC Manual
        #####################################################
        mdstepstocalc_all = np.array([1000,5000,9000,10000,20000,50000,90000,100000,120000,140000,160000,180000,200000,500000,\
                900000,1000000,120000,124000,126000,128000,2000000,5000000,9000000,10000000,20000000,50000000,90000000,100000000])
        #####################################################
        # MDSTEPSTOCALC last order of magnitude to check for error
        #####################################################
        mdstepstocalc_all = np.arange(11)*0
        verb = True
        verb = False
        if verb:
            print(mdstepstocalc_all,type(mdstepstocalc_all))
        steps_tocheck_md_error = int(space_fft.shape[1]/10)
        if verb:
            print(steps_tocheck_md_error,type(steps_tocheck_md_error))
        for idx,i in enumerate((np.arange(10)+1)[::-1]):
            add = int(steps_tocheck_md_error*10 - i*steps_tocheck_md_error)
            #print "add:",add,type(add)
            mdstepstocalc_all[i]=int(add)
            #np.append(mdstepstocalc_all,add)
        mdstepstocalc_all = np.trim_zeros(mdstepstocalc_all)

        # appendlast:
        # das /10*10 macht aus 20000{1,9} 200000  aber aus 200010 macht es 200010
        mdstepstocalc_all = np.sort(np.append(mdstepstocalc_all,[space_fft.shape[1]/10*10]))

        #####################################################
        # MDSTEPSTOCALC schoenheitskorrekturen der zahlen
        #####################################################
        if verb:
            print("mdstepstocalc_all:",mdstepstocalc_all)
        mdstepstocalc_all = np.sort(np.unique(mdstepstocalc_all[np.where(mdstepstocalc_all <= space_fft.shape[1])[0]]))

        if verb:
            print("mdstepstocalc_all:",mdstepstocalc_all)
        #if type(mdstepstocalc) == bool:
        mdstepstocalc = mdstepstocalc_all
        if verb:
            print("mdstepstocalc (in the function 1):",mdstepstocalc,type(mdstepstocalc))

        #print "2 mdstepstocalc_all:",mdstepstocalc_all
        mdstepsmax=mdstepstocalc_all[-1]
        #print "2 mdstepstocalc:",mdstepstocalc
        #print mdstepstocalc.min(),-1.*(int(math.log10(mdstepstocalc.min()))-1)
        ordmag = int(math.log10(mdstepstocalc.min()))
        #print mdstepstocalc.min(),10**ordmag
        factchange=10**ordmag
        factchange=10**ordmag

        if verb:
            print("mdstepstocalc (in the function 2.0):",mdstepstocalc,type(mdstepstocalc))

        mdstepstocalc = np.unique(mdstepstocalc/factchange*factchange) # remove elements wich are too close
        if verb:
            print("mdstepstocalc (in the function 2.1):",mdstepstocalc,type(mdstepstocalc))
        mdstepstocalc[-1]=mdstepsmax
        mdstepstocalc = mdstepstocalc[np.where(mdstepstocalc > 999)[0]]
        mdstepstocalc_all = mdstepstocalc_all[np.where(mdstepstocalc_all > 999)[0]]
        if mdstepstocalc_all[0] >= 20000:
            mdstepstocalc_all = np.insert(mdstepstocalc_all,0,10000)
        if mdstepstocalc_all[0] >= 10000:
            mdstepstocalc_all = np.insert(mdstepstocalc_all,0,5000)
        if mdstepstocalc_all[0] >= 5000:
            mdstepstocalc_all = np.insert(mdstepstocalc_all,0,1000)
        #print "--> mdstepstocalc:", mdstepstocalc# ,mdstepstocalc[np.where(mdstepstocalc > 999)[0]]



        #print 'yy',mdstepstocalc
        #print "mdstepstocalc_all xx:",mdstepstocalc_all
        #sys.exit()
        time = mdstepstocalc*faktor
        #print "timeaa:",mdstepstocalc*faktor

        if verb:
            print("mdstepstocalc (in the function 3):",mdstepstocalc,type(mdstepstocalc))

        def add_times(time,max,step,stepver=[8,6,4,2,1]):
            verbose=False
            if verbose:
                print()
                print("time in ps:",time)
                print("time in steps:",time/faktor)
                print("time[0]:",time[0],"max:",max,"step:",step)
            # time[0] = 1000, step = 10000
            # time[0] = 1000, max = 100, step = 10
            if time[0] <= step:
                return time
            elif time[0] > step:
                maxall = np.array([time[0],max]).min()
                if verbose:
                    print("maxall:",maxall)
                    print("mdstepstocalc in:",time)
                for i in stepver:
                    if i*step < max:
                        if verbose:
                            print("i*step:",i*step)
                        if i*step < time.max():
                            time = np.concatenate((np.array([i*step]),time))
                #print "mdstepstocalc out:",time
                return np.unique(time)


        if verb:
            print("timeaa:",mdstepstocalc*faktor,mdstepstocalc)  # this is correct (*faktor is in ps)
        mdstepstocalc = add_times(mdstepstocalc*faktor,10000,2000)/faktor # add points 10-100 ps

        def add_to_mdsteps(mdstepstocalc,add):
            #print 'kk',mdstepstocalc,add
            if mdstepstocalc.max() >= add:
                if add not in mdstepstocalc:
                    mdstepstocalc = np.sort(np.append(mdstepstocalc,add))
            #print 'oo',mdstepstocalc,add
            return mdstepstocalc
        #verb = True
        mdstepstocalc = add_to_mdsteps(mdstepstocalc,15000)
        mdstepstocalc = add_to_mdsteps(mdstepstocalc,20000)

        if verb:
            print('mm1',mdstepstocalc)
            #print "timebb:",mdstepstocalc*faktor,mdstepstocalc
        mdstepstocalc = add_times(mdstepstocalc*faktor,1000,100)/faktor # add 100-1000 ps
        if verb:
            print('mm2',mdstepstocalc)
            #print "timecc:",mdstepstocalc*faktor,mdstepstocalc
        mdstepstocalc = add_times(mdstepstocalc*faktor,100,10,stepver=[1,2,3,4,5,6,7,8,9][::-1])/faktor # add points 10-100 ps
        if verb:
            print('mm3',mdstepstocalc)
            #print "timecc:",mdstepstocalc*faktor,mdstepstocalc
            #print "timedd:",mdstepstocalc*faktor,mdstepstocalc
        mdstepstocalc = add_times(mdstepstocalc*faktor,10,1)/faktor # add points 10-100 ps
        if verb:
            print('mm4',mdstepstocalc)
            #print "timecc:",mdstepstocalc*faktor,mdstepstocalc
            #print "timeee1:",mdstepstocalc*faktor
            #print "timeee2:",mdstepstocalc
        #sys.exit()
        # take only those where we have more than 400 structures
        mdstepstocalc = mdstepstocalc[np.where(mdstepstocalc >= 400)[0]]
        if verb:
            #print "timeee1--->:",mdstepstocalc*faktor
            print("timeee2--->:",mdstepstocalc)

        #sys.exit()
        args.mdstepstocalc_all = mdstepstocalc.astype(int)
        #print "#####",args.mdstepstocalc_all
        print("inserting into db ...",args.mdstepstocalc_all)
        #insert_to_db_mdsteps_all(c,conn,args.mdstepstocalc_all)

        #sys.exit('77')
    if type(args.mdstepstocalc_all) == np.ndarray and os.path.isfile("sql.mdsteps") != True:
        print("saving....")
        np.savetxt("sql.mdsteps",args.mdstepstocalc_all.astype(int))


    if type(args.mdstepstocalc_all_ps) != np.ndarray and type(args.mdstepstocalc_all) == np.ndarray:
        args.mdstepstocalc_all_ps = (args.mdstepstocalc_all*faktor).astype(int)

    if type(args.mdsteps) != bool:
        return np.array(args.mdsteps)

    if type(args.mdstepstocalc_last) != bool:
        return np.array([args.mdstepstocalc_all[-1].astype(int)])

    if type(args.mdstepstocalc_all) == np.ndarray:
        #print "it was already known ..."
        return args.mdstepstocalc_all.astype(int)
    else:
        return False


def space_fft_to_powerspectrum(        qpoint, cute_timeinversion = True, idxprint =0, idxprintmax = 0,idxx=0,idxx_sum=0,space_fft_which=False):
    ''' calculates the power_spectrum for a certain qpoint (from space_fft)
    - c is the sqlite3 databse which has to be open and is closed afterwards externally
    - space_fft contins for a single q point all (symmetry) equivalent q points
    - output of this skript is the power_spectrum for the corresponding q-point
    - qpoint is just for saving the power_spectrum
    '''
    global args
    qpointstr = qpointstring(qpoint)
    if args.verbose:
        print("qpoint:",qpoint)
        print("execute_timeinversion:",execute_timeinversion)
        print("idxprint:",idxprint)
        print("idxprintmax:",idxprintmax)
    dt = args.dt
    check_qpoint_for_crossing = check_qpoints_for_crossing(qpoint, args.supercell, args.structure,space_fft_which=space_fft_which)
    if args.verbose:
        print("check_qpoint_for_crossing:",check_qpoint_for_crossing)
        print("space_fft_which 778899   :",space_fft_which)
        #sys.exit('778899')
        #print "dt:",dt
        #print "----------------------------------- all args ----------------------------"
        #print "args 99:",args
        #print "----------------------------------- all args ----------------------------"
    loadfile = space_fft_which+qpointstring(qpoint)+".npy"
    if os.path.isfile(loadfile) != True:
        print(printred("file "+loadfile+" does not exist!"))
        return

    # this does not ned to be loaded when all power_spectrum_files exist!
    #print printblue('loading '+loadfile)
    path = qpoint_get_l_0_0_or_l_l_0_or_l_l_l_or_t_0_0_or_t_t_t_or_t1_t1_t1_or_t2_t2_t2_from_qpoint(qpoint,args)
    #print printblue(str(idxx+1)+"/"+str(idxx_sum)+" np.load("+loadfile+"); Timeinv. "+str(execute_timeinversion)+"; args.seriell "+str(args.seriell)+"; crossing "+str(check_qpoint_for_crossing)+"; space_fft_which: "+space_fft_which+" "+path)
    space_fft = np.load(loadfile)
    #print "np.load(space_fft) done!",type(space_fft),space_fft.shape
    #if args.verbose: print "execute_timeinversion:",execute_timeinversion
    if execute_timeinversion != True:
        print_warning("execute_timeinversion is False!!!!")

    if dt == False:
        sys.exit("ERROR: need dt!")

    #####################################################################################
    # get mdstepstocalc (das koennte eigentlich schon ganz am anfang gemacht werden, bis auf space_fft)
    #####################################################################################
    #print '??? 0',args.mdstepstocalc_all,"args dt?",args.dt
    mdstepstocalc = get_mdstepstocalc(space_fft=space_fft) #,args=args)
    #print '??? 1',mdstepstocalc
    if args.mdstepstocalc_all_show == True:
        args.mdstepstocalc_all_show = False  # just show this information once
        if args.verbose:
            print("-----> args.mdstepstocalc_all  (before get_mdstepstocalc)  :")

        amd = args.mdstepstocalc_all.astype(int)
        ams = args.mdstepstocalc_all_ps.astype(int)
        if idxx < 1:
            print("args.mdstepstocalc_all   : [",amd[0],amd[1],"...",amd[-2],amd[-1],"]",len(args.mdstepstocalc_all),'idxprint',idxprint)
            print("args.mdstepstocalc_all_ps: [",ams[0],ams[1],"...",ams[-2],ams[-1],"]")
            if len(mdstepstocalc) <= 5:
                print("mdstepstocalc                    :",mdstepstocalc.astype(int))   # DONT REMOVE THIS< THIS IS REALLY CALCULATED
                print("mdstepstocalc (ps)               :",mdstepstocalc*0.001*args.dt) # DONT REMOVE THIS< THIS IS REALLY CALC
                #sys.exit()
            if len(mdstepstocalc) > 5:
                amd=mdstepstocalc.astype(int)
                ams=(mdstepstocalc*0.001*args.dt).astype(int)
                print("mdstepstocalc     : [",amd[0],amd[1],"...",amd[-2],amd[-1],"]",len(mdstepstocalc))
                print("mdstepstocalc (ps): [",ams[0],ams[1],"...",ams[-2],ams[-1],"]")
                #sys.exit()

    if not args.notverbose:
        #print printblue(str(idxx+1)+"/"+str(idxx_sum)+" np.load("+loadfile+"); executre Timeinversion is "+str(execute_timeinversion)+"; args.seriell is "+str(args.seriell)+"; crossing "+str(check_qpoint_for_crossing)+";")
        path = qpoint_get_l_0_0_or_l_l_0_or_l_l_l_or_t_0_0_or_t_t_t_or_t1_t1_t1_or_t2_t2_t2_from_qpoint(qpoint,args)
        q1 = qpoint[0]
        q2 = qpoint[1]
        q3 = qpoint[2]
        length = qpoint_to_length_qpoint(q1,q2,q3,args)
        add = ""
        if space_fft_which == "space_fft_":
            add = "      "
        #print printblue(str(idxx+1)+"/"+str(idxx_sum)+" np.load("+loadfile+"); Timeinv. "+str(execute_timeinversion)+"; args.seriell "+str(args.seriell)+"; crossing "+str(check_qpoint_for_crossing)+"; "),
        print(printgreen(add+str(path)+" --> "+str(length)))
    if os.path.isfile(loadfile) != True:
        sys.exit("ERROR 44: :"+loadfile+" does not exist! use either -fftpy option (slow but implemented) or use the -fftc optein (fast but currently only working with lammps jobs) to get it!")
    #print '??? 2',args.mdstepstocalc_all
    ################################################################################
    # calculate the lifeimes
    ################################################################################
    get_lifetimes_for_different_timesteps_ser_par(mdstepstocalc,qpoint,space_fft,idxprint=idxx,idxprintmax=idxx_sum,args=args,space_fft_which=space_fft_which)
    #print '??? 3',args.mdstepstocalc_all
    return

def get_lifetimes_for_different_timesteps_ser_par(mdstepstocalc,qpoint,space_fft,idxprint=998,idxprintmax=998,args=False,space_fft_which=False):
    ''' If I am not mistaken, this is done for every q-point'''
    qpointstr = qpointstring(qpoint)
    jobsa = []
    serial = False

    #print "done7xx"
    if args.seriell: # seriell !!!!!!!!!!!!!!!!!!!
        start = time.time()
        for idx,mdsteps in enumerate(mdstepstocalc[::-1]):
            #print "filenameout:",filenameout
            get_lifetime_for_particular_mdsteps(mdsteps,mdstepstocalc,space_fft,qpoint,args=args,space_fft_which=space_fft_which)

    if not args.seriell: # parallel  !!!!!!!!!!!!!!!!!
        #print "parallel"
        def currentfunc_helper(args):
            return get_lifetime_for_particular_mdsteps(*args)

        # print printblue('not serial!')
        start = time.time()
        manager = multiprocessing.Manager()
        return_dict = manager.dict()
        #print "rd in:",return_dict
        jobs = []
        for idx,i in enumerate(mdstepstocalc):
            arguments = (i,mdstepstocalc[::-1],space_fft,qpoint,args,space_fft_which)
            p = multiprocessing.Process(target=currentfunc_helper, args=[arguments])
            jobs.append(p)
            p.start()
        for job in jobs:
            job.join()
    return


def space_fft_to_powerspectrum_ser_par(qpoints_all, execute_timeinversion = True, args = False,idxx=0,idxx_sum=0,space_fft_which='space_fft_'):
    ######################################################################
    # if seriell
    ######################################################################
    sys.exit('ho889')
    if args.seriell:  # !! seriell
        for idx,i in enumerate(qpoints_all):
            space_fft_to_powerspectrum(i,execute_timeinversion,idx,len(qpoints_all),idxx,idxx_sum,space_fft_which=space_fft_which)
    ######################################################################
    # if parallel, split to chunks of 40, otherwise problems using osx
    ######################################################################
    if not args.seriell:  # !!! parallel
        def chunks(l, n):
            """Yield successive n-sized chunks from l."""
            for i in range(0, len(l), n):
                yield l[i:i + n]
        for qpoints_all_chunk in chunks(qpoints_all,40):
            jobs = []
            for idx,i in enumerate(qpoints_all_chunk):
                def space_fft_to_powerspectrum_helper(args):
                    space_fft_to_powerspectrum(*args)
                arguments = (i,execute_timeinversion,idx,len(qpoints_all_chunk),idxx,idxx_sum,space_fft_which)
                p = multiprocessing.Process(target=space_fft_to_powerspectrum_helper, args=[arguments])
                jobs.append(p)
                p.start()

            for job in jobs: # just to wait until job is done
                job.join()
    return
# DELETE SINCE DEFINED TWICE!
#def ps_to_ps_for_writing(ps,xind,x_to_THz,function=False,x=False,popt=False,psmax=1):
#    ''' xind are the maxima of the function and is an array of one or several numbers '''
#    #if type(x) != bool:
#    #    print 'xx',x,len(x)
#    #    #ps_alt = eval(take+take2)(x, *popt)
#    #    ps_alt = function(x, *popt)
#    #    print 'ps_alt',ps_alt,len(ps_alt)
#    #    print 'lenps',len(ps)
#
#    THz_max = 20
#    THz_max = 10  # for bcc Ti
#    x_to_THz_full = x_to_THz
#    THz_max_test = np.array(xind).max()*x_to_THz_full
#    #print "xind:",xind,"THz_max_test:",THz_max_test
#    ################################################################
#    # Here we define the max frequency to be twice the peak frequency
#    ################################################################
#    if THz_max_test*2 < THz_max:
#        THz_max =  THz_max_test*2
#    if type(xind) != bool:
#        if len(xind) > 1:
#            THz_max = 23
#    #print "len(xind):",len(xind)
#    #print "THz_max:",THz_max
#
#
#    factor_psout = 2  # 2 is a save bet since power spectrum is doubled
#    #np.savetxt('kka',ps)
#    #yout = ps[:len(ps)/factor_psout]  # das ist um den zweiten teil (zweite seite) des powerspektrums
#                                      # wegzubekommen; bei [t2,t2,0] kommt aber schon nur die haelfte rein
#                                      # hmmm schein generell nicht mehr notwendig zu zein.
#    yout = ps
#    #print '   1 yout.max()',yout.max()
#    #np.savetxt('kkb',yout)
#    xout = np.arange(len(yout))*x_to_THz_full
#    #print xout,len(xout)
#    xout1 = xout[np.where(xout < THz_max)[0]]
#    yout1 = yout[:len(xout1)]
#    #print '   2 yout1.max()',yout1.max()
#    #print 'yout.max():',yout1.max()
#    #print 'yout[-1]  :',yout1[-1]
#    #print 'ratio     :',yout1.max()/yout1[-1]
#    #print 'xout1.max()',xout1.max()
#    #print "ps1:",ps.max()
#    minsteps = 500
#    #print 'type(x)',type(x)
#    #print 'type(function)',type(function)
#    #print 'type(popt)',type(popt)
#    #print 'len(xout1)',len(xout1)
#    #print 'minsteps:',minsteps
#    #print 'popt',popt
#    if type(x) != bool and type(function) != bool and type(popt) != bool and len(xout1) < minsteps:
#        lt_min = np.array(popt[2:][::3]*x_to_THz).min()
#        #print 'lt_min:',lt_min
#        if lt_min < 0.2:
#            minsteps = 4000
#            minsteps = 40000  # war mal fuer Irinas CrN notwendig
#        #print 'lenx',len(xout1),len(yout1),THz_max,'ll',len(ps),len(xout1)
#        for i in [0.1,0.01, 0.001]:
#            xout1 = np.arange(0,THz_max,i)
#            if len(xout1) > minsteps:
#                break
#        yout1 = function(xout1/x_to_THz, *popt)*psmax
#        #print 'xout2.max()',xout1.max()
#        #print "ps2:",yout1.max()
#        #np.savetxt('kkk',yout1)
#        #sys.exit()
#    return xout1,yout1
#

def make_powerspectrum_for_particular_mdstep_length(space_fft,mdsteps,args):
    power_spectrum=np.zeros(int(mdsteps))
    serial = True
    if serial:  # serial
        #print printblue('fftserial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        for i in np.arange(len(space_fft)): # laeuft ueber alle symmetrieequivalten qs (6 stueck oder 3 stucek i == [ 0, 1, 2])
            #startt = time.time()
            ### MAIN
            #@ amplitude = (space_fft[i][:mdsteps])[:mdsteps]
            #@ #print 'in ampl.shape:',ampl.shape,ampl  # real and imag
            #@ ##############################################################################
            #@ ######### old definition MAIN (WORKING)
            #@ ##############################################################################
            #@ fourier_transformierte = fft(amplitude)
            #@ #print "old ft.shape:",ft.shape,ft       # real and imag
            #@ a =np.abs(fourier_transformierte)**2.                       # only real

            if args.verbose > 4:
                print('yostage 333',i)
            a =np.abs(fft((space_fft[i][:int(mdsteps)])[:int(mdsteps)]))**2.  # put into one command to save (hopefully) memory
            if args.verbose > 4:
                print('yostage 3334',i)
            # timeinversion
            #aa = np.roll(a[::-1],1)
            #a += aa
            ##############################################################################
            ######### new definition MAIN
            ##############################################################################
            #@ampl = (space_fft[i][:mdsteps])[:mdsteps]
            #@amplka=np.empty(len(ampl),dtype=np.complex64)
            #@amplka.real=ampl.real
            #@amplka.imag=-ampl.imag

            #@intensity = np.abs(ampl)**2.            # only real
            #@np.savetxt("i"+str(i),intensity)
            #@intensity = ampl*amplka
            #@#print "new intensity.shape:",intensity.shape,intensity
            #@b = fft(intensity)                      # real and imag
            #@bb = fft(b)                      # real and imag
            #@np.savetxt("bb.real"+str(i),bb.real)
            #@power_spectrumb+=b

            #@#print "a.shape:",a.shape,a
            #@#print "b.shape:",b.shape,b
            #@#print

            # WRONG a =np.abs(fft(((space_fft[i][:mdsteps])[:mdsteps])**2))  # put into one command to save (hopefully) memory
            # DOES NOT WORK a =fft(np.abs((space_fft[i][:mdsteps])[:mdsteps])**2)  # put into one command to save (hopefully) memory
            #b =fft(np.abs(((space_fft[i][:mdsteps])[:mdsteps]))**2.)
            #print "n1a:",a.shape,type(a),a[:2]
            #print "n1b:",b.shape,type(b),b[:2]

            print_single_qpoints = False  # done for dominique for bcc
            if print_single_qpoints == True and mdsteps == mdstepstocalc.max():
                tmpfilename = "ps"+qpointstring(qpoint)+"__"+str(i)+"__"+str(mdsteps)+".dat"
                parnpsavetxt(tmpfilename,a)

            power_spectrum+=a
        #endp = time.time()
        #print "fft1 ",str((endp-startp)),"sec."

    if args.verbose > 4:
        print('yostage 234')

    if not serial: # parallel
        ## DONT! FOR SMALL JOBS THIS IS A BIT FASTER approx 15% but for larger jobs this makes
        ## everything slower by approx 17.5%,  maybe due to using more memory?
        sys.exit()
        ## DONT! FOR SMALL JOBS THIS IS A BIT FASTER approx 15% but for larger jobs this makes
        ## everything slower,  maybe due to using more memory?
        #print printblue('fftparallel !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        startp = time.time()
        def currentfuncc_helper(args):
            return thefunc(*args)
        def thefunc(i,space_fft,mdsteps,return_dictt):
            calc = np.abs(fft((space_fft[i][:mdsteps])[:mdsteps]))**2.
            return_dictt[i] = calc
            return calc

        manager = multiprocessing.Manager()
        return_dictt = manager.dict()
        jobs = []
        for i in np.arange(len(space_fft)): # laeuft ueber alle symmetrieequivalten qs (6 stueck)
            arguments = (i,space_fft,mdsteps,return_dictt)
            p = multiprocessing.Process(target=currentfuncc_helper, args=[arguments])
            jobs.append(p)
            p.start()
        for job in jobs:
            job.join()
        #print "#################33"
        X=np.zeros(mdsteps)
        d = return_dictt
        for k, v in list(d.items()):
            #print k, v
            power_spectrum+=v
        #print "don#################33"
        #print d
        #for k,v in d: X[k] = v
        #print type(d.values())
        #power_spectrum=X
        return_dictt = None
        #endp = time.time()
        #print "fft1 ",str((endp-startp)),"sec."
        #print



    #
    #np.savetxt("psb_pos.real",power_spectrumb.real)

    #power_spectrumb += np.roll(power_spectrumb[::-1],1)
    #print power_spectrumb
    #np.savetxt("ps_pos",power_spectrum)
    if args.verbose > 4:
        print('yostage 344')

    execute_timeinversion = True
    if execute_timeinversion:
        #a = np.roll(a[::-1],1)
        # this can only be done as done in the first definition
        #np.savetxt("ps_neg",np.roll(power_spectrum[::-1],1))
        power_spectrum += np.roll(power_spectrum[::-1],1)
        print_single_qpoints = False # for dominique
        if print_single_qpoints == True and mdsteps == mdstepstocalc.max():
            tmpfilename = "ps"+qpointstring(qpoint)+"__"+str(i)+"-__"+str(mdsteps)+".dat"
            parnpsavetxt(tmpfilename,a)
        #power_spectrum+=a
    #end = time.time()
    #print "fft ",mdsteps,i,str((end-startt)),"sec."

    if args.verbose > 4:
        print('yostage 345')
    ##########################################################################
    # write out power_spectrum (simple list, only y without x values;
    # length of mdsteps/2 which changes in this loop)
    ##########################################################################

    #power_spectrum[:2] = 0;  # power_spectrum is just a list of values (only y without x)
    #power_spectrum[-2:] = 0;  # power_spectrum is just a list of values (only y without x)

    #np.savetxt('ps_orig',power_spectrum)
    #np.savetxt('ps_orig_xy',np.transpose(np.array([np.arange(len(power_spectrum)),power_spectrum])))
    #sys.exit()

    ##########################################################################
    # make powerspectrum symmetric (necessery although execute_timeinversion is used!)
    ##########################################################################
    power_spectrum = power_spectrum+power_spectrum[::-1]
    if power_spectrum[2] < 0.:
        sys.exit('ps is < 0!!')
    return power_spectrum


##def load_or_make_powerspectrum_for_particular_mdstep_length(space_fft,space_fft_which,qpointstr,mdsteps,args,power_spectrum_filename=False):
##    print('power_spectrum_filename',power_spectrum_filename)
##    ############################# a) load if possible ##########################
##    if power_spectrum_filename != False:
##        if os.path.isfile(power_spectrum_filename):
##            print(printblue('loading from (1)'+power_spectrum_filename))
##            power_spectrum = np.load(power_spectrum_filename)
##        else:
##            sys.exit('this does not exist 77')
##            #return False
##
##
##    ############################# a) create if necessary ##########################
##    if power_spectrum_filename == False:
##        power_spectrum_filename = "ps_save/ps_"+space_fft_which+qpointstr+"_"+str(mdsteps)+".npy"
##        if os.path.isfile(power_spectrum_filename):
##            print(printblue('loading from (2)'+power_spectrum_filename))
##            power_spectrum = np.load(power_spectrum_filename)
##        else:
##            print('does not exist',power_spectrum_filename)
##            if type(space_fft) == bool:
##                print('creating powerspectrum! (1)',power_spectrum_filename)
##
##                loadfile = space_fft_which+qpointstring(qpoint)+".npy"
##                if os.path.isfile(loadfile) != True:
##                    print(printred("file "+loadfile+" does not exist!"))
##                    return False
##
##                # this does not ned to be loaded when all power_spectrum_files exist!
##                space_fft = np.load(loadfile)
##                #print "np.load(space_fft) done!",type(space_fft),space_fft.shape
##                power_spectrum = make_powerspectrum_for_particular_mdstep_length(space_fft,mdsteps,args)
##            else:
##                print('ss',space_fft.shape)
##                power_spectrum = make_powerspectrum_for_particular_mdstep_length(space_fft,mdsteps,args)
##    if not os.path.isfile(power_spectrum_filename):
##        np.save(power_spectrum_filename,power_spectrum)
##    return power_spectrum


def get_lifetime_for_particular_mdsteps(mdsteps,mdstepstocalc,space_fft,qpoint,args=False,space_fft_which=False):
    ''' at this stage the space_fft is already loaded '''
    qpointstr = qpointstring(qpoint)
    power_spectrumb=np.empty(int(mdsteps),dtype=np.complex64)
    if args.verbose > 4:
        print('yostage 123')
    ##########################################################################
    # get power spectrum (needs to stay all the time with both peaks until end for fitting)
    ##########################################################################

    path = qpoint_get_l_0_0_or_l_l_0_or_l_l_l_or_t_0_0_or_t_t_t_or_t1_t1_t1_or_t2_t2_t2_from_qpoint(qpoint,args)
    print(printblue(str(idxx+1)+"/"+str(idxx_sum)+" np.load("+loadfile+"); Timeinv. "+str(execute_timeinversion)+"; args.seriell "+str(args.seriell)+"; crossing "+str(check_qpoint_for_crossing)+"; space_fft_which: "+space_fft_which))
    space_fft = np.load(loadfile)
    #print "np.load(space_fft) done!",type(space_fft),space_fft.shape
    #if args.verbose: print "execute_timeinversion:",execute_timeinversion
    if execute_timeinversion != True:
        print_warning("execute_timeinversion is False!!!!")

    if dt == False:
        sys.exit("ERROR: need dt!")

    #####################################################################################
    # get mdstepstocalc (das koennte eigentlich schon ganz am anfang gemacht werden, bis auf space_fft)
    #####################################################################################
    #print '??? 0',args.mdstepstocalc_all,"args dt?",args.dt
    mdstepstocalc = get_mdstepstocalc(space_fft=space_fft) #,args=args)
    #print '??? 1',args.mdstepstocalc_all
    #sys.exit()
    if args.mdstepstocalc_all_show == True:
        args.mdstepstocalc_all_show = False  # just show this information once
        if args.verbose:
            print("-----> args.mdstepstocalc_all  (before get_mdstepstocalc)  :")

        #print "space_fft.shape",space_fft.shape
        #print "space_fft.shape",space_fft.shape[1]
        #if args.verbose:
        amd = args.mdstepstocalc_all.astype(int)
        ams = args.mdstepstocalc_all_ps.astype(int)
        if idxx < 1:
            print("args.mdstepstocalc_all   : [",amd[0],amd[1],"...",amd[-2],amd[-1],"]",len(args.mdstepstocalc_all),'idxprint',idxprint)
            print("args.mdstepstocalc_all_ps: [",ams[0],ams[1],"...",ams[-2],ams[-1],"]")
            if len(mdstepstocalc) <= 5:
                print("mdstepstocalc                    :",mdstepstocalc.astype(int))   # DONT REMOVE THIS< THIS IS REALLY CALCULATED
                print("mdstepstocalc (ps)               :",mdstepstocalc*0.001*args.dt) # DONT REMOVE THIS< THIS IS REALLY CALC
                #sys.exit()
            if len(mdstepstocalc) > 5:
                amd=mdstepstocalc.astype(int)
                ams=(mdstepstocalc*0.001*args.dt).astype(int)
                print("mdstepstocalc     : [",amd[0],amd[1],"...",amd[-2],amd[-1],"]",len(mdstepstocalc))
                print("mdstepstocalc (ps): [",ams[0],ams[1],"...",ams[-2],ams[-1],"]")
                #sys.exit()

    if not args.notverbose:
        #print printblue(str(idxx+1)+"/"+str(idxx_sum)+" np.load("+loadfile+"); executre Timeinversion is "+str(execute_timeinversion)+"; args.seriell is "+str(args.seriell)+"; crossing "+str(check_qpoint_for_crossing)+";")
        path = qpoint_get_l_0_0_or_l_l_0_or_l_l_l_or_t_0_0_or_t_t_t_or_t1_t1_t1_or_t2_t2_t2_from_qpoint(qpoint,args)
        q1 = qpoint[0]
        q2 = qpoint[1]
        q3 = qpoint[2]
        length = qpoint_to_length_qpoint(q1,q2,q3,args)
        add = ""
        if space_fft_which == "space_fft_":
            add = "      "
        #print printblue(str(idxx+1)+"/"+str(idxx_sum)+" np.load("+loadfile+"); Timeinv. "+str(execute_timeinversion)+"; args.seriell "+str(args.seriell)+"; crossing "+str(check_qpoint_for_crossing)+"; "),
        print(printgreen(add+str(path)+" --> "+str(length)))
    if os.path.isfile(loadfile) != True:
        sys.exit("ERROR 44: :"+loadfile+" does not exist! use either -fftpy option (slow but implemented) or use the -fftc optein (fast but currently only working with lammps jobs) to get it!")
    #print '??? 2',args.mdstepstocalc_all
    ################################################################################
    # calculate the lifeimes
    ################################################################################
    get_lifetimes_for_different_timesteps_ser_par(mdstepstocalc,qpoint,space_fft,idxprint=idxx,idxprintmax=idxx_sum,args=args,space_fft_which=space_fft_which)
    #print '??? 3',args.mdstepstocalc_all
    return

def get_lifetimes_for_different_timesteps_ser_par(mdstepstocalc,qpoint,space_fft,idxprint=998,idxprintmax=998,args=False,space_fft_which=False):
    ''' If I am not mistaken, this is done for every q-point'''
    qpointstr = qpointstring(qpoint)
    jobsa = []
    serial = False

    #print "done7xx"
    if args.seriell: # seriell !!!!!!!!!!!!!!!!!!!
        start = time.time()
        for idx,mdsteps in enumerate(mdstepstocalc[::-1]):
            #print "filenameout:",filenameout
            get_lifetime_for_particular_mdsteps(mdsteps,mdstepstocalc,space_fft,qpoint,args=args,space_fft_which=space_fft_which)

    if not args.seriell: # parallel  !!!!!!!!!!!!!!!!!
        #print "parallel"
        def currentfunc_helper(args):
            return get_lifetime_for_particular_mdsteps(*args)

        # print printblue('not serial!')
        start = time.time()
        manager = multiprocessing.Manager()
        return_dict = manager.dict()
        #print "rd in:",return_dict
        jobs = []
        for idx,i in enumerate(mdstepstocalc):
            arguments = (i,mdstepstocalc[::-1],space_fft,qpoint,args,space_fft_which)
            p = multiprocessing.Process(target=currentfunc_helper, args=[arguments])
            jobs.append(p)
            p.start()
        for job in jobs:
            job.join()
    return


def space_fft_to_powerspectrum_ser_par(qpoints_all, execute_timeinversion = True, args = False,idxx=0,idxx_sum=0,space_fft_which='space_fft_'):
    ######################################################################
    # if seriell
    ######################################################################
    if args.seriell:  # !! seriell
        for idx,i in enumerate(qpoints_all):
            space_fft_to_powerspectrum(i,execute_timeinversion,idx,len(qpoints_all),idxx,idxx_sum,space_fft_which=space_fft_which)
    ######################################################################
    # if parallel, split to chunks of 40, otherwise problems using osx
    ######################################################################
    if not args.seriell:  # !!! parallel
        def chunks(l, n):
            """Yield successive n-sized chunks from l."""
            for i in range(0, len(l), n):
                yield l[i:i + n]
        for qpoints_all_chunk in chunks(qpoints_all,40):
            jobs = []
            for idx,i in enumerate(qpoints_all_chunk):
                def space_fft_to_powerspectrum_helper(args):
                    space_fft_to_powerspectrum(*args)
                arguments = (i,execute_timeinversion,idx,len(qpoints_all_chunk),idxx,idxx_sum,space_fft_which)
                p = multiprocessing.Process(target=space_fft_to_powerspectrum_helper, args=[arguments])
                jobs.append(p)
                p.start()

            for job in jobs: # just to wait until job is done
                job.join()
    return

def ps_to_ps_for_writing(ps,xind,x_to_THz,function=False,x=False,popt=False,psmax=1,args=False):
    ''' xind are the maxima of the function and is an array of one or several numbers '''
    #if type(x) != bool:
    #    print 'xx',x,len(x)
    #    #ps_alt = eval(take+take2)(x, *popt)
    #    ps_alt = function(x, *popt)
    #    print 'ps_alt',ps_alt,len(ps_alt)
    #    print 'lenps',len(ps)

    THz_max = 20
    if args.structure ==  'bcc':
        THz_max = 10
    x_to_THz_full = x_to_THz
    THz_max_test = np.array(xind).max()*x_to_THz_full
    #print "xind:",xind,"THz_max_test:",THz_max_test
    ################################################################
    # Here we define the max frequency to be twice the peak frequency
    ################################################################
    if THz_max_test*2 < THz_max:
        THz_max =  THz_max_test*2
    if type(xind) != bool:
        if len(xind) > 1:
            THz_max = 23
    #print "len(xind):",len(xind)
    #print "THz_max:",THz_max


    factor_psout = 2  # 2 is a save bet since power spectrum is doubled
    #np.savetxt('kka',ps)
    #yout = ps[:len(ps)/factor_psout]  # das ist um den zweiten teil (zweite seite) des powerspektrums
                                      # wegzubekommen; bei [t2,t2,0] kommt aber schon nur die haelfte rein
                                      # hmmm schein generell nicht mehr notwendig zu zein.
    yout = ps
    #print '   1 yout.max()',yout.max()
    #np.savetxt('kkb',yout)
    xout = np.arange(len(yout))*x_to_THz_full
    #print xout,len(xout)
    xout1 = xout[np.where(xout < THz_max)[0]]
    yout1 = yout[:len(xout1)]
    #print '   2 yout1.max()',yout1.max()
    #print 'yout.max():',yout1.max()
    #print 'yout[-1]  :',yout1[-1]
    #print 'ratio     :',yout1.max()/yout1[-1]
    #print 'xout1.max()',xout1.max()
    #print "ps1:",ps.max()
    minsteps = 500
    #print 'type(x)',type(x)
    #print 'type(function)',type(function)
    #print 'type(popt)',type(popt)
    #print 'len(xout1)',len(xout1)
    #print 'minsteps:',minsteps
    #print 'popt',popt
    if type(x) != bool and type(function) != bool and type(popt) != bool and len(xout1) < minsteps:
        lt_min = np.array(popt[2:][::3]*x_to_THz).min()
        #print 'lt_min:',lt_min
        if lt_min < 0.2:
            minsteps = 4000
            minsteps = 40000  # war mal fuer Irinas CrN notwendig
        #print 'lenx',len(xout1),len(yout1),THz_max,'ll',len(ps),len(xout1)
        for i in [0.1,0.01, 0.001]:
            xout1 = np.arange(0,THz_max,i)
            if len(xout1) > minsteps:
                break
        yout1 = function(xout1/x_to_THz, *popt)*psmax
        #print 'xout2.max()',xout1.max()
        #print "ps2:",yout1.max()
        #np.savetxt('kkk',yout1)
        #sys.exit()
    return xout1,yout1


def make_powerspectrum_for_particular_mdstep_length(space_fft,mdsteps,args):
    power_spectrum=np.zeros(int(mdsteps))
    serial = True
    if serial:  # serial
        #print printblue('fftserial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        for i in np.arange(len(space_fft)): # laeuft ueber alle symmetrieequivalten qs (6 stueck oder 3 stucek i == [ 0, 1, 2])
            #startt = time.time()
            ### MAIN
            #@ amplitude = (space_fft[i][:mdsteps])[:mdsteps]
            #@ #print 'in ampl.shape:',ampl.shape,ampl  # real and imag
            #@ ##############################################################################
            #@ ######### old definition MAIN (WORKING)
            #@ ##############################################################################
            #@ fourier_transformierte = fft(amplitude)
            #@ #print "old ft.shape:",ft.shape,ft       # real and imag
            #@ a =np.abs(fourier_transformierte)**2.                       # only real

            if args.verbose > 4:
                print('yostage 333',i)
            a =np.abs(fft((space_fft[i][:int(mdsteps)])[:int(mdsteps)]))**2.  # put into one command to save (hopefully) memory
            if args.verbose > 4:
                print('yostage 3334',i)
            # timeinversion
            #aa = np.roll(a[::-1],1)
            #a += aa
            ##############################################################################
            ######### new definition MAIN
            ##############################################################################
            #@ampl = (space_fft[i][:mdsteps])[:mdsteps]
            #@amplka=np.empty(len(ampl),dtype=np.complex64)
            #@amplka.real=ampl.real
            #@amplka.imag=-ampl.imag

            #@intensity = np.abs(ampl)**2.            # only real
            #@np.savetxt("i"+str(i),intensity)
            #@intensity = ampl*amplka
            #@#print "new intensity.shape:",intensity.shape,intensity
            #@b = fft(intensity)                      # real and imag
            #@bb = fft(b)                      # real and imag
            #@np.savetxt("bb.real"+str(i),bb.real)
            #@power_spectrumb+=b

            #@#print "a.shape:",a.shape,a
            #@#print "b.shape:",b.shape,b
            #@#print

            # WRONG a =np.abs(fft(((space_fft[i][:mdsteps])[:mdsteps])**2))  # put into one command to save (hopefully) memory
            # DOES NOT WORK a =fft(np.abs((space_fft[i][:mdsteps])[:mdsteps])**2)  # put into one command to save (hopefully) memory
            #b =fft(np.abs(((space_fft[i][:mdsteps])[:mdsteps]))**2.)
            #print "n1a:",a.shape,type(a),a[:2]
            #print "n1b:",b.shape,type(b),b[:2]

            print_single_qpoints = False  # done for dominique for bcc
            if print_single_qpoints == True and mdsteps == mdstepstocalc.max():
                tmpfilename = "ps"+qpointstring(qpoint)+"__"+str(i)+"__"+str(mdsteps)+".dat"
                parnpsavetxt(tmpfilename,a)

            power_spectrum+=a
        #endp = time.time()
        #print "fft1 ",str((endp-startp)),"sec."

    if args.verbose > 4:
        print('yostage 234')

    if not serial: # parallel
        ## DONT! FOR SMALL JOBS THIS IS A BIT FASTER approx 15% but for larger jobs this makes
        ## everything slower by approx 17.5%,  maybe due to using more memory?
        sys.exit()
        ## DONT! FOR SMALL JOBS THIS IS A BIT FASTER approx 15% but for larger jobs this makes
        ## everything slower,  maybe due to using more memory?
        #print printblue('fftparallel !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        startp = time.time()
        def currentfuncc_helper(args):
            return thefunc(*args)
        def thefunc(i,space_fft,mdsteps,return_dictt):
            calc = np.abs(fft((space_fft[i][:mdsteps])[:mdsteps]))**2.
            return_dictt[i] = calc
            return calc

        manager = multiprocessing.Manager()
        return_dictt = manager.dict()
        jobs = []
        for i in np.arange(len(space_fft)): # laeuft ueber alle symmetrieequivalten qs (6 stueck)
            arguments = (i,space_fft,mdsteps,return_dictt)
            p = multiprocessing.Process(target=currentfuncc_helper, args=[arguments])
            jobs.append(p)
            p.start()
        for job in jobs:
            job.join()
        #print "#################33"
        X=np.zeros(mdsteps)
        d = return_dictt
        for k, v in list(d.items()):
            #print k, v
            power_spectrum+=v
        #print "don#################33"
        #print d
        #for k,v in d: X[k] = v
        #print type(d.values())
        #power_spectrum=X
        return_dictt = None
        #endp = time.time()
        #print "fft1 ",str((endp-startp)),"sec."
        #print



    #
    #np.savetxt("psb_pos.real",power_spectrumb.real)

    #power_spectrumb += np.roll(power_spectrumb[::-1],1)
    #print power_spectrumb
    #np.savetxt("ps_pos",power_spectrum)
    if args.verbose > 4:
        print('yostage 344')

    execute_timeinversion = True
    if execute_timeinversion:
        #a = np.roll(a[::-1],1)
        # this can only be done as done in the first definition
        #np.savetxt("ps_neg",np.roll(power_spectrum[::-1],1))
        power_spectrum += np.roll(power_spectrum[::-1],1)
        print_single_qpoints = False # for dominique
        if print_single_qpoints == True and mdsteps == mdstepstocalc.max():
            tmpfilename = "ps"+qpointstring(qpoint)+"__"+str(i)+"-__"+str(mdsteps)+".dat"
            parnpsavetxt(tmpfilename,a)
        #power_spectrum+=a
    #end = time.time()
    #print "fft ",mdsteps,i,str((end-startt)),"sec."

    if args.verbose > 4:
        print('yostage 345')
    ##########################################################################
    # write out power_spectrum (simple list, only y without x values;
    # length of mdsteps/2 which changes in this loop)
    ##########################################################################

    #power_spectrum[:2] = 0;  # power_spectrum is just a list of values (only y without x)
    #power_spectrum[-2:] = 0;  # power_spectrum is just a list of values (only y without x)

    #np.savetxt('ps_orig',power_spectrum)
    #np.savetxt('ps_orig_xy',np.transpose(np.array([np.arange(len(power_spectrum)),power_spectrum])))
    #sys.exit()

    ##########################################################################
    # make powerspectrum symmetric (necessery although execute_timeinversion is used!)
    ##########################################################################
    power_spectrum=power_spectrum+power_spectrum[::-1]
    return power_spectrum


def load_or_make_powerspectrum_for_particular_mdstep_length(space_fft,space_fft_which,qpointstr,mdsteps,args,power_spectrum_filename=False):
    #print('power_spectrum_filename',power_spectrum_filename)
    ############################# a) load if possible ##########################
    if type(power_spectrum_filename) != bool:
        if os.path.isfile(power_spectrum_filename):
            power_spectrum = np.load(power_spectrum_filename)
            if args.verbose > 1 or mdsteps == args.mdstepstocalc_all.max():
                print(printblue('loading 1 from '+power_spectrum_filename),power_spectrum.shape)
            return power_spectrum
        else:
            return False


    ############################# a) create if necessary ##########################
    if power_spectrum_filename == False:
        power_spectrum_filename = "ps_save/ps_"+space_fft_which+qpointstr+"_"+str(mdsteps)+".npy"
        #print 'power_spectrum_filename 2',power_spectrum_filename
        #sys.exit()
        if os.path.isfile(power_spectrum_filename):
            power_spectrum = np.load(power_spectrum_filename)
            if args.verbose > 1 or mdsteps == args.mdstepstocalc_all.max():
                print(printblue('loading 2 from '+power_spectrum_filename),power_spectrum.shape)
            # dont know what this is for #np.savetxt(power_spectrum_filename+".dat",power_spectrum)
        else:
            #print
            #print 'power_spectrum_filename (dne) will be created...',power_spectrum_filename,qpointstr,type(space_fft)
            if type(space_fft) == bool:
                #print 'qpointstring',qpointstr
                #print qpointstring_to_qpoint(qpointstr)
                loadfile = space_fft_which+qpointstr+".npy"
                #print 'loadfile:',loadfile
                #print 'creating powerspectrum! (2)',loadfile,printred("-->"),power_spectrum_filename
                if os.path.isfile(loadfile) != True:
                    print(printred("file "+loadfile+" does not exist!"))
                    return False

                # this does not ned to be loaded when all power_spectrum_files exist!
                space_fft = np.load(loadfile)
                #print "np.load(space_fft) done!",type(space_fft),space_fft.shape
                power_spectrum = make_powerspectrum_for_particular_mdstep_length(space_fft,mdsteps,args)
                #np.savetxt(loadfile+".dat",power_spectrum)
            else:
                #print 'ss1',space_fft.shape,"all"
                power_spectrum = make_powerspectrum_for_particular_mdstep_length(space_fft,mdsteps,args)
                #print 'ss2',len(power_spectrum),mdsteps
        if not os.path.isfile(power_spectrum_filename):
            np.save(power_spectrum_filename,power_spectrum)
    return power_spectrum


def get_lifetime_for_particular_mdsteps(mdsteps,mdstepstocalc,space_fft,qpoint,args=False,space_fft_which=False):
    ''' at this stage the space_fft is already loaded '''
    qpointstr = qpointstring(qpoint)
    power_spectrumb=np.empty(int(mdsteps),dtype=np.complex64)
    if args.verbose > 4:
        print('yostage 123')
    ##########################################################################
    # get power spectrum (needs to stay all the time with both peaks until end for fitting)
    ##########################################################################

    path = qpoint_get_l_0_0_or_l_l_0_or_l_l_l_or_t_0_0_or_t_t_t_or_t1_t1_t1_or_t2_t2_t2_from_qpoint(qpoint,args)
    #print 'path',path
    #print 'space_fft_which',space_fft_which
    saveorig=False
    ##################################################################################################################
    # *** if this is True, this aparently loads _onepeak for t2_t2_t2
    # *** e.g. try: rm */*4_4_8*;
    ##################################################################################################################
    loadonepeak = False
    if space_fft_which == 'space_fft_phase_' and path == "t2_t2_0" and loadonepeak == True:
        #if path == "t2_t2_0":
        #print 'if1'
        power_spectrum_filename = "ps_save/ps_"+space_fft_which+qpointstr+"_"+str(mdsteps)+"_onepeak.npy"
        power_spectrum = load_or_make_powerspectrum_for_particular_mdstep_length(space_fft,space_fft_which,qpointstr,mdsteps,args,power_spectrum_filename = power_spectrum_filename)
        #print("type(power_spectrum)",type(power_spectrum))
        if type(power_spectrum) == bool:
            power_spectrum_orig = load_or_make_powerspectrum_for_particular_mdstep_length(False,space_fft_which,qpointstr,mdsteps,args)
            #print 'get long'
            #np.savetxt('psorig',power_spectrum_orig)
            qpoint_ll0, qpoint_t1_t1_0 = qpoint_map_t2_qpoint_to_L_L_0_and_T1_T1_0_qpoint(qpoint,args)
            power_spectrum_long = load_or_make_powerspectrum_for_particular_mdstep_length(False,space_fft_which,qpointstring(qpoint_ll0),mdsteps,args)
            #print 'get trans'
            #np.savetxt('psl',power_spectrum_long)
            power_spectrum_t1 = load_or_make_powerspectrum_for_particular_mdstep_length(False,space_fft_which,qpointstring(qpoint_t1_t1_0),mdsteps,args)
            #np.savetxt('pst',power_spectrum_t1)
            power_spectrum = np.copy(power_spectrum_orig - power_spectrum_long-power_spectrum_t1)
            #np.savetxt('psout',power_spectrum)
        if not os.path.isfile(power_spectrum_filename):
            #print('saveing x')
            np.save(power_spectrum_filename,power_spectrum)
        #print 'done'
        #print power_spectrum.shape
        #print power_spectrum_long.shape
        #print power_spectrum_t1.shape
        #np.savetxt('psall',power_spectrum)
        #np.savetxt('psl',power_spectrum_long)
        #np.savetxt('pst',power_spectrum_t1)
        #np.savetxt('psout',power_spectrum-power_spectrum_long-power_spectrum_t1)
        #sys.exit()



    ##################################################################################################################
    # *** if this is set to true, this aparently loads for t2_t2_t2 the whole spectrum (more close to experiment)
    ##################################################################################################################
    elif space_fft_which == 'space_fft_' and path == "t2_t2_0" and loadonepeak == False:
        power_spectrum_onepeak_filename = "ps_save/ps_"+space_fft_which+qpointstr+"_"+str(mdsteps)+"_onepeak.npy"
        # lade das volle "onepeak" powerspectrum
        power_spectrum = load_or_make_powerspectrum_for_particular_mdstep_length(space_fft,space_fft_which,qpointstr,mdsteps,args,power_spectrum_filename = power_spectrum_onepeak_filename)

        vvv = True
        if vvv:
            print('type(pow)',type(power_spectrum))
            print('ps ..',power_spectrum)
            print('saving it')
            #np.savetxt('ps',power_spectrum)
            #print 'ps1',power_spectrum
            #sys.exit('123')
            #print('ps oben',power_spectrum)
            #sys.exit()
        #if type(power_spectrum) != bool:
        #    if mdsteps == args.mdstepstocalc_all.max():
        #        print(printblue("WRITING0 : "+"ps_full/"+filenameout_full).rstrip('\n'))
        #    x_to_THz_full = 1000./(len(power_spectrum)*args.dt)
        #    xout,yout = ps_to_ps_for_writing(ps,xind,x_to_THz_full,args=args)
        #    parnpsavetxt("ps_full/"+filenameout_full,np.transpose(np.array([xout,yout])),fmt='%.8f %.5f')
        #if type(power_spectrum) == bool:
        #print('XXX ps:',power_spectrum)
        if type(power_spectrum) == bool:
            saveorig=True
            filenameout_full = "ps_"+space_fft_which+qpointstring(qpoint)+"_"+str(mdsteps)+".dat.full"
            power_spectrum_orig = load_or_make_powerspectrum_for_particular_mdstep_length(False,space_fft_which,qpointstr,mdsteps,args)
            #np.savetxt('ps_orig',power_spectrum_orig)
            #sys.exit('123')
            qpoint_ll0, qpoint_t1_t1_0 = qpoint_map_t2_qpoint_to_L_L_0_and_T1_T1_0_qpoint(qpoint,args)
            #print 'qpoint_ll0',qpoint_ll0,'qpoint_t1_t1_0',qpoint_t1_t1_0
            #np.savetxt('psorig',power_spectrum_orig)
            power_spectrum_long = load_or_make_powerspectrum_for_particular_mdstep_length(False,space_fft_which,qpointstring(qpoint_ll0),mdsteps,args)
            #np.savetxt('ps_orig_long',power_spectrum_long)
            #sys.exit('123')
            #np.savetxt('psl',power_spectrum_long)
            #print 'done2'


            #power_spectrum_t1 = load_or_make_powerspectrum_for_particular_mdstep_length(False,space_fft_which,qpointstring(qpoint_t1_t1_0),mdsteps,args)
            #np.savetxt('pst',power_spectrum_t1)
            #np.savetxt('pst1BZ',power_spectrum_t1_1BZ)
            longmax_x = np.where(power_spectrum_long == power_spectrum_long.max())[0][0]
            longmax_y = power_spectrum_long[longmax_x]
            origmax_y = power_spectrum_orig[longmax_x]
            #print('longmax_x',longmax_x)
            #print('longmax_y',longmax_y)
            #print('origmax_y',origmax_y)
            power_spectrum_long = (power_spectrum_long/longmax_y)*origmax_y
            #np.savetxt('xx_orig_long',power_spectrum_long)
            #np.savetxt('xx_orig',power_spectrum_long)
            #sys.exit('123')

            power_spectrum = np.copy(power_spectrum_orig - power_spectrum_long) #-power_spectrum_t1)
            filenameout_full_min_long = "ps_"+space_fft_which+qpointstring(qpoint)+"_"+str(mdsteps)+".dat.full_min_long"
            power_spectrum_full_min_long = np.copy(power_spectrum)
            #np.savetxt('xx_orig_minus_long',power_spectrum)
            #sys.exit('123')
            #np.savetxt('psnew',power_spectrum)
            #sys.exit('123')


            ## This is a good move for 4 4 8  -> reduces 2 peaks to one peak!
            ## This is a wrong move for 1 1 8, 2 2 8, 3 3 8  --> have only one peak
            ## Therefore make this only for qlen > 0.75
            if False:
                qlen = qpoint_to_length_qpoint(qpoint[0],qpoint[1],qpoint[2],args)
                if qlen > 0.75:
                    #print 'qlen',qlen,'now get t1 from first BZ == space_fft_phase_'
                    # this is not perfect but ok
                    print('now loading space_fft_phase_'+qpointstring(qpoint_t1_t1_0)+"_"+str(mdsteps))
                    power_spectrum_t1_1BZ = load_or_make_powerspectrum_for_particular_mdstep_length(False,'space_fft_phase_',qpointstring(qpoint_t1_t1_0),mdsteps,args)
                    #np.savetxt("power_spectrum_t1_1BZ",power_spectrum_t1_1BZ)
                    if type(power_spectrum_t1_1BZ) == bool:
                        return False
                    power_spectrum = np.copy(power_spectrum - power_spectrum_t1_1BZ) #-power_spectrum_t1)
                    #if power_spectrum[2] < 0.:
                    #    sys.exit('smaller 0 is ps!!')
            if not os.path.isfile(power_spectrum_onepeak_filename):
                np.save(power_spectrum_onepeak_filename,power_spectrum)
            #np.savetxt('psnew',power_spectrum)
            #sys.exit('34rttu')

    else:
        #print(printblue('else no going mc'))
        power_spectrum = load_or_make_powerspectrum_for_particular_mdstep_length(space_fft,space_fft_which,qpointstr,mdsteps,args)





    #np.savetxt("ps",power_spectrum)
    #sys.exit()
    #np.savetxt("psb_pos_neg.real",power_spectrumb.real)
    #np.savetxt("psb_pos_neg.imag",power_spectrumb.imag)
    #sys.exit('done23e435;w')

    ps=power_spectrum
    filenameout = "ps_"+space_fft_which+qpointstring(qpoint)+"_"+str(mdsteps)+".dat"
    if ps[2] < 0.:
        print('see psxxx')
        np.savetxt('psxxx',ps)
        sys.exit('ps < 0')

    ##########################################################################
    # WRITE out power_spectrum full
    # x_ind_from full
    ##########################################################################
    check_qpoint_for_crossing = check_qpoints_for_crossing(qpoint, args.supercell, args.structure,space_fft_which=space_fft_which)
    #if args.verbose:
    #    print "check_qpoint_for_crossing:",check_qpoint_for_crossing
    #    print "space_fft_which 556677   :",space_fft_which
    #    #sys.exit('556677')
    get_xind_from_full = False
    xind        = False     # depends on mdsteps
    xind_in_THz = False     # does not depend on mdsteps
    #if mdsteps == args.mdstepstocalc_all.max() and
    writefull = False
    if args.mdstepstocalc_all.max() <= 50000 and mdsteps <= 50000:
        writefull = True
    if mdsteps == args.mdstepstocalc_all.max():
        writefull = True

    if args.write_full_ps:
        writefull = True
    #print 'writefull',writefull
    #sys.exit()
    #if writefull:
    needit=False
    if not os.path.isfile("ps_full/"+filenameout):
        needit = True
    if saveorig and not os.path.isfile("ps_full/"+filenameout_full):
        needit = True
    if saveorig and not os.path.isfile("ps_full/"+filenameout_full_min_long):
        needit = True

    #print('XXX needit:',needit)
    if needit:
        #print('needit')
        if not os.path.isfile("ps_full/"+filenameout):
            #    print printblue("EXISTING0: "+"ps_full/"+filenameout).rstrip('\n')
            x_to_THz_full = 1000./(args.dt*float(mdsteps))
            x_to_THz = np.copy(x_to_THz_full)
            #print "x_to_THz_full:",x_to_THz_full
            if check_qpoint_for_crossing == True:
                xind = get_ind_of_x_peaks(ps,peaks=args.peaks*2,x_to_THz=x_to_THz_full,sep_in_THz_min=1.0,verbose=False,qpoint=qpoint,space_fft_which=space_fft_which)
                print('xind @ full ps (unsmoothed + CROSSING):',np.round(xind*x_to_THz_full,2),'(mdsteps=',mdsteps,')')
            else:
                xind = get_ind_of_x_peaks(ps,peaks=args.peaks  ,x_to_THz=x_to_THz_full,sep_in_THz_min=args.separation_acoustic_optical,qpoint=qpoint,space_fft_which=space_fft_which)


            xind_in_THz = np.copy(xind*x_to_THz_full)
            get_xind_from_full = True
            #print 'xind @ full ps (unsmoothed):',np.round(xind*x_to_THz_full,2),'(mdsteps=',mdsteps,')'
            #print 'xind @ full 99 (unsmoothed):',xind,'x_to_THz',x_to_THz

            #if os.path.isfile("ps_full/"+filenameout):
            #    print printblue("EXISTING0: "+"ps_full/"+filenameout).rstrip('\n')
            #else:
            if mdsteps == args.mdstepstocalc_all.max():
                print(printblue("WRITING0 : "+"ps_full/"+filenameout).rstrip('\n'))
            x_to_THz_full = 1000./(len(ps)*args.dt)
            xout,yout = ps_to_ps_for_writing(ps,xind,x_to_THz_full,args=args)
            parnpsavetxt("ps_full/"+filenameout,np.transpose(np.array([xout,yout])),fmt='%.8f %.5f')

            if saveorig and not os.path.isfile("ps_full/"+filenameout_full):
                print(printblue("WRITING0 : "+"ps_full/"+filenameout_full).rstrip('\n'))
                xout,yout = ps_to_ps_for_writing(power_spectrum_orig,xind,x_to_THz_full,args=args)
                parnpsavetxt("ps_full/"+filenameout_full,np.transpose(np.array([xout,yout])),fmt='%.8f %.5f')

            if saveorig and not os.path.isfile("ps_full/"+filenameout_full_min_long):
                print(printblue("WRITING0 : "+"ps_full/"+filenameout_full_min_long).rstrip('\n'))
                xout,yout = ps_to_ps_for_writing(power_spectrum_full_min_long,xind,x_to_THz_full,args=args)
                parnpsavetxt("ps_full/"+filenameout_full_min_long,np.transpose(np.array([xout,yout])),fmt='%.8f %.5f')


    #sys.exit('writefulloutxxxx')
    #print 'ka',get_xind_from_full,args.verbose
    if get_xind_from_full == True and args.verbose:
        print('xind from full:',np.round(xind_in_THz,2),'(mdsteps=',mdsteps,')')
    #sys.exit()
    smooth=True   # this seems closer to Experiment
    smooth=False  # this seems more accurate from Theoretical standpoint
    ps = power_spectrum
    for i in np.arange(100):
        psmax = power_spectrum.max()
        counnnt = np.where(ps==psmax)[0]
        #print(i,"counnnt",counnnt)
        if counnnt[0] == 0:
            ps[0] = ps[i+1]
            ps[-1] = ps[i+1]
            #print(i,"xxxcounnnt",counnnt)
        else:
            break
        #if type(counnnt) == bool:
        #    break
        #if counnnt == 0:
        #    break
    #sys.exit()
    #############################################################################
    # HIER KANN ps zwar was verschoben sein, aber xy ist genau richtig!
    # beim powerspectrum_to_powerspectrum_2d_sparse muss man sicherstellen das man
    # mindestens 3 punkte (oder besser 5) hat. Wenn schon das psin (das rohe powerspektrum)
    # nur 6 punkte hat, dann sollte man da am bestn nichts tun (was gleichzusetzen ist mit
    # einer erhoehung von mindestens=5000 zu hoeher).
    #############################################################################
    pslenfull=len(ps)
    #print 'pslenfull',pslenfull
    #xy,x,ps,x_to_THz = powerspectrum_to_powerspectrum_2d_sparse(ps,smooth=smooth,mindestens=50000)
    xy,x,ps,x_to_THz = powerspectrum_to_powerspectrum_2d_sparse(ps,smooth=smooth,mindestens=1000)
    pslencurr=len(ps)
    #print 'pslencurr',pslencurr
    x_to_THz = 1000./(args.dt*float(mdsteps))*pslenfull/pslencurr
    if type(xind_in_THz) != bool:
        xind = xind_in_THz/x_to_THz
    #print 'ppp full',pslenfull
    #print 'ppp curr',pslencurr
    #print "x_to_THz (oben):",x_to_THz


    writesmooth = False
    if mdsteps == args.mdstepstocalc_all.max() or type(args.mdsteps) != bool or args.write_smooth_in:
        writesmooth = True
    if writesmooth == True:
        if os.path.isfile("ps_smooth/"+filenameout): # and False:
            print(printblue("EXISTING2: ps_smooth/"+filenameout),'(mdsteps=',mdsteps,')','compacted by',pslencurr/pslenfull)
        else:
            print(printblue("WRITING2 : ps_smooth/"+filenameout),'(mdsteps=',mdsteps,')')
            #parnpsavetxt("ps_smooth/"+filenameout,np.transpose(np.array([x*x_to_THz,ps])),fmt='%.8f %.5f')
            xout,yout = ps_to_ps_for_writing(ps,x,x_to_THz,args=args)
            parnpsavetxt("ps_smooth/"+filenameout,np.transpose(np.array([xout,yout])),fmt='%.8f %.5f')

    #sys.exit('645rt45')
    #sys.exit()
    ##########################################################################
    # WRITE pssparse if necessary
    ##########################################################################
    writeps_smooth_in = False
    #if args.write_smooth_in: writeps_smooth_in = True
    #if args.write_full_ps and mdsteps == mdstepstocalc.max(): writeps_smooth_in = True
    if check_qpoint_for_crossing == True and mdsteps == mdstepstocalc.max(): writeps_smooth_in = True
    if (args.write_smooth_in or args.write_full_ps) == True and mdsteps == mdstepstocalc.max():
        if pslenfull != pslencurr:
            if not args.notverbose:
                print(printblue("WRITING1: "+filenameout+".smooth_in ... (this is written at the last step when mdsteps == mdstepstocalc.max()) "+qpointstring(qpoint)+"  --->"+str(mdsteps)))
            #print "saveunten"
            #np.savetxt("psoutxyunten",xy)  # --> xy scheint ok zu sein
            parnpsavetxt("ps_full/"+filenameout+".smooth_in",xy[:xy.shape[0]/2])

    #ps=ps/ps.max()*.9; # full with peaks left and right
    ##########################################################################
    # CROSSING
    # WRITE out power_spectrum_cut
    # (needs to stay all the time with both peaks until end for fitting)
    ##########################################################################
    #if False:
    q1 = float(qpointstr.split("_")[0])
    q3 = float(qpointstr.split("_")[2])
    if check_qpoint_for_crossing == True:
        q1_q3=q1/(q3/2.)
        if q1_q3 > 0.77:
            if mdsteps == mdstepstocalc.max():
                checkwriteps = mdsteps
            else:
                checkwriteps = False
            cut_the_ps = True
            if cut_the_ps:
                pscut = get_max_min_max_coords_idx_of_full_ps_when_bandcrossing(xy,x,ps,qpointstr,checkwriteps,verbose=False)  # ps kommt hier skaliert raus auf 0.9 (wobei nur der andere peak auch hoecher sein kann)
                #np.savetxt("ps1",ps)
                pscut2 = pscut/pscut.max()*psmax # danach geht ps hoch auch ca. 10**12
                if mdsteps == mdstepstocalc.max():
                    fn = "ps_smooth/"+filenameout+".cut"
                    if os.path.isfile(fn) and False:
                        print(printblue("EXISTING3: "+fn))
                    else:
                        print(printblue("WRITING3 : "+fn))
                        #parnpsavetxt(fn,np.transpose(np.array([x*x_to_THz,pscut2])),fmt='%.8f %.5f')
                        xout,yout = ps_to_ps_for_writing(pscut2,x,x_to_THz,args=args)
                        parnpsavetxt(fn,np.transpose(np.array([xout,yout])),fmt='%.8f %.5f')

    ######################################################
    # good: auswertung !!!!!!!!!!!!!!!
    ######################################################
    if check_qpoint_for_crossing == True:
        if q1_q3 > 0.77: goodps = pscut2
        else: goodps = ps
    else: goodps = ps

    if True:
        sigmagood, fqgood,ltgood,ltgoodmin,ltgoodmax,xmaxwritegood,goodsmoothfunc = get_goodsmoothing_for_ps(goodps,args.dt,mdsteps,allowed_func = 0, allowed_der = 0,args=args,stringadd=['good',str(mdsteps)])
        outstr = '     fq: '+str(fqgood)+"   lt:"+str(ltgood)

        ########################################################
        # good: save to db
        ########################################################
        if True:
            global pppx_add_db
            peakidx=1
            add_value_to_db_tmp(qpoint,mdsteps,peakidx,space_fft_which+'good',fqgood,ltgood,0.0)
            #print qpoint,mdsteps,peakidx,space_fft_which+'good',fqgood,ltgood,0.0

        ########################################################
        # good: write if necessary
        ########################################################
        if args.write_fitted_ps == True and mdsteps == args.mdstepstocalc_all.max():
            filenameoutgood = filenameout+"_good"+str(sigmagood)+".dat"
            if os.path.isfile("ps_fitted/"+filenameoutgood):
                print(printblue("EXISTING4: ps_fitted/"+filenameoutgood+outstr))
            else:
                print(printblue("WRITING4 : ps_fitted/"+filenameoutgood+outstr))
                xout,yout = ps_to_ps_for_writing(goodsmoothfunc,x,x_to_THz,args=args)
                parnpsavetxt("ps_fitted/"+filenameoutgood,np.transpose(np.array([xout,yout])),fmt='%.8f %.5f')

    #sys.exit()
    #######################################################
    # get estimate for the frequency
    #######################################################
    #x_to_THz = 1000./(args.dt*float(mdsteps))*pslenfull/pslencurr
    peaks_lookfor = args.peaks

    if check_qpoint_for_crossing == True: peaks_lookfor = args.peaks*2
    qlen = qpoint_to_length_qpoint(qpoint[0],qpoint[1],qpoint[2],args,get_lll_double_length=True)
    if qlen == 0.0:
        peaks_lookfor = 1
    #if space_fft_which == 'space_fft_phase_':
    #    path = qpoint_get_l_0_0_or_l_l_0_or_l_l_l_or_t_0_0_or_t_t_t_or_t1_t1_t1_or_t2_t2_t2_from_qpoint(qpoint,args)
    #    if path == "t2_t2_0":
    #        peaks_lookfor = args.peaks*3

    if args.verbose > 1:
        print()
        print("&&&&&&&&&&&&&&&&&&&&&&&&& xind get      &&&&&&&&&&&&&&&&&&&&&&&&")
        print("&&&& mdsteps                             :",mdsteps)
        print("&&&& pslenfull                           :",pslenfull)
        print("&&&& pslencurr                           :",pslencurr)
        print("&&&& pslenfull/pslencurr                 :",pslenfull/pslencurr)
        print("&&&& xind 000                            :",xind)
        print("&&&& x_to_THz (unten)                    :",x_to_THz)
        print("&&&& check_qpoint_for_crossing           :",check_qpoint_for_crossing)
        print("&&&& args.separation_acoustic_optical    :",args.separation_acoustic_optical)
        print('&&&& peaks look for                      :',peaks_lookfor)


    #sys.exit()
    #######################################################
    # get estimate for the frequency from db
    #######################################################
    estimate_from_db = False
    if mdsteps != args.mdstepstocalc_all.max() and os.path.isfile('sql.tmp'):
        filename_for_db = False
        estimate_from_db = get_estimate_w_a_G(qpoint,filename_for_db,mdsteps,verbose=False,args=args,add_filename=filename,x_to_THz=x_to_THz)
        if type(estimate_from_db) != bool:
            xind = estimate_from_db[0::3]/x_to_THz


            if len(xind) != peaks_lookfor:
                xind = False

        if args.verbose:
            if type(xind) != bool:
                print('xind $ db xind:',np.round(xind*x_to_THz,2))
            else:
                print('xind $ db xind:',xind)


    if type(xind) == bool:
        if check_qpoint_for_crossing != True:
            if args.verbose:
                print('get_ind_of_x_peaks ....')
            xind = get_ind_of_x_peaks(ps,peaks=args.peaks  ,x_to_THz=x_to_THz,sep_in_THz_min=args.separation_acoustic_optical,qpoint=qpoint,space_fft_which=space_fft_which)

        if check_qpoint_for_crossing == True:
            # sep_in_THz_min = 0.1 is not enout for /cmmc/u/aglen/TiN/222_300K_SUM 2 2 4
            xind = get_ind_of_x_peaks(ps,peaks=peaks_lookfor,x_to_THz=x_to_THz,sep_in_THz_min=0.3,qpoint=qpoint,space_fft_which=space_fft_which)
            print('&&&& xind @ ccc                          :',xind*x_to_THz)
            #if type(xind) == bool and q1_q3 < 0.77: # second peak was not found since too close to the first one
            if type(xind) == bool:
                print('xind was bool',xind)
                xind = get_ind_of_x_peaks(ps,peaks=args.peaks,x_to_THz=x_to_THz,sep_in_THz_min=1.0,qpoint=qpoint,space_fft_which=space_fft_which)
                xind = np.array([0.98*xind[0],1.005*xind[0]])
                print('CRUDE ESTIMATE OF FREQ since found only one peak qpoint:',qpoint,"mdsteps:",mdsteps,"xind:",xind*x_to_THz)


    #qlen = qpoint_to_length_qpoint(qpoint[0],qpoint[1],qpoint[2],args,get_lll_double_length=True)
    if qlen == 0.0: # and len(xind) == 2:
        #print('qlen',qlen)
        #print("lne(xind):",len(xind))
        #print("xind:",xind)
        xind = get_ind_of_x_peaks(ps,peaks=1,x_to_THz=x_to_THz,sep_in_THz_min=args.separation_acoustic_optical,qpoint=qpoint,space_fft_which=space_fft_which)
        #print("xind:",xind)
        #sys.exit()

    if args.verbose:
        print('xind @ xind   :',np.round(xind*x_to_THz,2),'sep_in_THz_min:',args.separation_acoustic_optical,'args.peaks',args.peaks,'xind:',xind,'type(xind):',type(xind))
        #sys.exit('111111223')

    #sys.exit('nach xind')

    def fit_new_lorenz(ps_in,x_to_THz,xind,verbose,psmax,qpoint,mdsteps,args,space_fft_which):
        qpointstr = qpointstring(qpoint)
        check_qpoint_for_crossing = check_qpoints_for_crossing(qpoint, args.supercell, args.structure,space_fft_which)
        if args.verbose:
            print("check_qpoint_for_crossing:",check_qpoint_for_crossing)
            print("space_fft_which 445566   :",space_fft_which)
            #sys.exit('445566')
        #peaks_lookfor = args.peaks
        #if check_qpoint_for_crossing == True: peaks_lookfor = args.peaks*2
        #print "xind:",xind
        #np.savetxt("psin",ps_in)
        xind = np.array(xind)
        THz_to_x = 1./x_to_THz
        psmax=ps_in.max()
        #print "psmax (1):",psmax
        ps = ps_in/psmax   # ohne diesen schritt haut das fitting nicht hin
        #print "psmax (2):",ps.max(),ps.max()/psmax
        xvalues=len(ps)
        yall = ps
        xall = np.arange(xvalues)

        yl = ps[0:int(xvalues/2)]   # only take one peak
        xl = np.arange(len(yl))
        xr = xl+xl[-1]+1
        yr = ps[int(xvalues/2):]


        x = xa = xall
        y = ya = yall
        bounds=False
        max_fq=0
        min_fq=0
        max_lt=0
        pl=False
        pr=False
        if args.verbose > 1:
            print('xind in fit_new_lorenz',np.round(xind*x_to_THz,2),'mdsteps',mdsteps)
            print('len(xind) in fit_new_lorenz',len(xind))

        if len(xind) == 1:
            pl="1p_0p"
            pa="1p_1p"
            dp = [0.001,0.1,0.005]
            min_fq = xind[0]*0.85
            max_fq = xind[0]*1.5   # especially for the short md runs of 2,3,4 ps and so on it is necessary to be more flexible
            max_lt = 4. * THz_to_x #  4 THz
            bounds=((min_fq,0,0),(max_fq,np.inf,max_lt)) # geht nach redifinition von a in lorenzhauer
            #bounds=((min_fq,0,0),(max_fq,np.inf,3000)) # geht nach redifinition von a in lorenzhauer

        if len(xind) == 2:
            pl="2p_0p"
            pa="2p_2p"
            dp = [0.001,0.1,0.005,0.001,0.1,0.005]
            min_fq_l = xind[0]*0.85
            max_fq_l = xind[0]*1.3
            min_fq_r = xind[1]*0.85
            max_fq_r = xind[1]*1.3
            max_lt = 5. * THz_to_x #  5 THz
            bounds=((min_fq_l,0,0,min_fq_r,0,0),(max_fq_l,np.inf,max_lt,max_fq_r,np.inf,max_lt))

        if len(xind) == 4:
            pl="4p_0p"
            pa="4p_4p"
            dp = [0.001,0.1,0.005,0.001,0.1,0.005,0.001,0.1,0.005,0.001,0.1,0.005]
            min_fq_ll = xind[0]*0.85
            max_fq_ll = xind[0]*1.3
            min_fq_l  = xind[1]*0.85
            max_fq_l  = xind[1]*1.3
            min_fq_r  = xind[2]*0.85
            max_fq_r  = xind[2]*1.3
            min_fq_rr = xind[3]*0.85
            max_fq_rr = xind[3]*1.3
            max_lt = 5. * THz_to_x #  5 THz
            bounds=((min_fq_ll,0,0, min_fq_l,0,0, min_fq_r,0,0, min_fq_rr,0,0),(max_fq_ll,np.inf,5,max_fq_l,np.inf,5, max_fq_r,np.inf,5,max_fq_rr,np.inf,5))

        if args.verbose > 3:
            #np.savetxt("ps_in",np.transpose(np.array([xall*100,ps_in])))
            print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
            print("xall       :",len(xall),"[",xall[:3],"...",xall[-3:])
            print("xl         :",len(xl),"[",xl[:3],"...",xl[-3:])
            print("xr         :",len(xr),"[",xr[:3],"...",xr[-3:])
            print()
            print("yall       :",len(yall),"[",yall[:3],"...",yall[-3:])
            print("yl         :",len(yl),"[",yl[:3],"...",yl[-3:])
            print("yr         :",len(yr),"[",yr[:3],"...",yr[-3:])
            print("xind       :",xind)
            print("xind [THz] :",xind*x_to_THz)
            print("len(xind)  :",len(xind))
            print("ps[xind]   :",ps[xind.astype(int)])
            print("x_to_THz   :",x_to_THz)
            print("psmax      :",psmax)
            print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")

        def ppplot(popt,x_to_THz,im,bounds=False):
            rt=2 # number of digits
            bs="F"
            if type(bounds) != bool:
                bs="T"
            #return str(np.round(popt*x_to_THz,rt))+"  "+str(np.round(popt,rt))+" im:"+str(im)
            if type(popt) == bool:
                popt = np.array([0,0,0])

            out=""
            for peakidx in np.arange(len(popt)/3):
                peakidx = int(peakidx)
                string0  = str(np.round(popt[0+peakidx*3]*x_to_THz,rt)).ljust(6)
                string1  = str(np.round(popt[1+peakidx*3]*x_to_THz,rt)).ljust(8)
                string2  = str(np.round(popt[2+peakidx*3]*x_to_THz,rt)).ljust(6)
                string0o = str(np.round(popt[0+peakidx*3],rt)).ljust(6)
                string1o = str(np.round(popt[1+peakidx*3],rt)).ljust(8)
                string2o = str(np.round(popt[2+peakidx*3],rt)).ljust(6)
                string3=string4=string5=""
                string3o=string4o=string5o=""

                out1=printred("["+string0+"  "+string1+"  "+string2+"]") #+" im:"+str(im)+"   B:"+bs #+" ["+string0o+"  "+string1o+"  "+string2o+"]"

                if peakidx==0:
                    out +=out1
                else:
                    out+="\n"+out1
            if len(popt) == 3:
                return out
            else:
                return "\n"+out

        def makefit(take=False,strr="",x=x,y=y,sigma=False, kwargs=False,b=False,args=False):
            '''
            strr muss xa (whole powerspektrum left/right)
            oder xl (only left part of ps) sein!
            xind is an array with [[xpos1, ppos2, ...], x_to_THz]
            kwargs = [peaks, qpointstr, xl,yl,xa,ya,verbose,mdsteps,check_qpoint_for_crossing,qpoint,psmax]
            kwargs = [xind, qpointstr, xl,yl,xa,ya,verbose,mdsteps,check_qpoint_for_crossing,qpoint,psmax,x_to_THz]
            '''
            if strr != "xa" and strr != "xl":
                sys.exit("strr has to be xl or xa")


            xind = kwargs[0]
            qpointstr = kwargs[1]
            xl = kwargs[2]
            yl = kwargs[3]
            xa = kwargs[4]
            ya = kwargs[5]
            verbose = kwargs[6]
            mdsteps = kwargs[7]
            check_qpoint_for_crossing = kwargs[8]
            qpoint  = kwargs[9]
            psmax = kwargs[10]
            x_to_THz = kwargs[11]
            space_fft_which = kwargs[12]


            peaks = len(xind)
            if peaks < 1:
                sys.exit('no peak found here')

            if strr == 'xa':
                x = xa
                y = ya
                take2 = str(len(xind))+'p_'+str(len(xind))+'p'

            if strr == 'xl':
                x = xl
                y = yl
                take2 = str(len(xind))+'p_0p'

            if type(sigma)==bool:
                if sigma==False:
                    strsigma="F"
                if sigma==True:
                    strsigma="T"
            elif type(sigma) == int or type(sigma)==float:
                strsigma=str(sigma)
            else:
                sys.exit("sigma format wrong")

            bs="F"
            if type(b) != bool:
                bs="T"
            filename=space_fft_which+take+strr+"_sigma_"+strsigma #+"_bounds_"+bs
            #print "xind:------------>",xind,xind.min()
            #print "dp:",dp

            ##################################################
            # if crossing get left peak from [l l 0]
            ##################################################
            check_qpoint_for_crossing = kwargs[8]
            lower_peaks_fq_lt=False

            if args.verbose > 1:
                print("lower_peaks_fq_lt:",lower_peaks_fq_lt)
            #print 'kwargs1:',kwargs
            kwargs.append(lower_peaks_fq_lt)
            kwargs.append(filename)
            kwargs.append(take)
            #print 'kwargs1:',kwargs



            #print 'mdsteps:',mdsteps
            #print 'xind:',xind
            #print '===',np.ones(len(xind)*3)*0.000001
            #print 'take+take2',take+take2
            try:
                #print 'x:',x
                #print 'y:',y
                #print 'sigma',sigma
                popt,im,function = fit_iterative(eval(take+take2),x,y,sigma=sigma,strr=strr,kwargs=kwargs,args=args)
            except RuntimeError:
                popt = np.ones(len(xind)*3)*0.000001
                im=0
                #function = False

            if type(popt) == bool or type(im) == bool or type(function) == bool:
                popt = np.ones(len(xind)*3)*0.000001
                im=0
                #function = False

            #####################################################
            # check if boundary conditions were hit
            #####################################################
            #print('pppp6',popt*x_to_THz)
            setpopt0=False
            for peakidx in np.arange(len(popt)/3):
                #print "peakidx:",peakidx
                addtoindex=peakidx*3
                #print('peakidx',peakidx)
                peakidx = int(peakidx)
                popt_fq  = round(float(popt[0+peakidx*3]),3)
                popt_aa  = round(float(popt[1+peakidx*3]),3)
                popt_lt  = round(float(popt[2+peakidx*3]),3)
                #if type(b) != bool:
                #    b_fq_min = round(float(b[0][0+peakidx*3]),3)
                #    b_fq_max = round(float(b[1][0+peakidx*3]),3)
                #    b_lt_min = round(float(b[0][2+peakidx*3]),3)
                #    b_lt_max = round(float(b[1][2+peakidx*3]),3)
                #else:
                b_fq_min = 0.
                b_fq_max = 24./x_to_THz
                b_lt_min = 0.
                b_lt_max = 5./x_to_THz
                if mdsteps < 12000:
                    b_lt_max = 7./x_to_THz

                b_lt_max = 8./x_to_THz   # for 12 12 0 in 4sc al necessary
                #if check_qpoint_for_crossing:
                #    b_lt_max = 4.6/x_to_THz

                popt_fq_THz  = round(popt_fq*x_to_THz,3)
                popt_aa_THz  = round(popt_aa*x_to_THz,4)
                popt_lt_THz  = round(popt_lt*x_to_THz,3)
                b_fq_min_THz = round(b_fq_min*x_to_THz,3)
                b_fq_max_THz = round(b_fq_max*x_to_THz,3)
                b_lt_min_THz = round(b_lt_min*x_to_THz,3)
                b_lt_max_THz = round(b_lt_max*x_to_THz,3)
                #print "b:",b,len(b[0])
                #print "p:",popt,len(popt)
                #print
                #print "b_fq_min",b_fq_min
                #print "popt_fq:",popt_fq
                #print "b_fq_max",b_fq_max
                #print
                #print "b_lt_min",b_lt_min
                #print "popt_lt:",popt_lt
                #print "b_lt_max",b_lt_max
                if popt_fq >= b_fq_max or popt_fq <= b_fq_min:
                    setpopt0=True
                    if mdsteps*0.001*args.dt >= 10 or args.verbose > 2: # if greater then 10 ps
                        print("ERROR (fq): qp:",qpointstr,'mdsteps:',mdsteps,filename," popt_fq reached boundary;",b_fq_min_THz,"<?",popt_fq_THz,"<?",b_fq_max_THz)
                        print('type(b):',type(b),np.round(b,decimals=2))
                        print('estimate: xind:',xind,'(ps)',xind*x_to_THz)
                        print('result: popt (ps):',popt*x_to_THz)

                if popt_lt >= b_lt_max or popt_lt <= b_lt_min:
                    setpopt0=True
                    #if args.verbose > 1:
                    if mdsteps*0.001*args.dt >= 10 or args.verbose > 2: # if greater then 10 ps
                        print("ERROR (lt): qp:",qpointstr,'mdsteps:',mdsteps," popt_lt reached boundary;",b_lt_min_THz,"<?",popt_lt_THz,"<?",b_lt_max_THz)
                        print('type(b):',type(b),b)
                        print('estimate: xind:',xind,'(ps)',xind*x_to_THz)
                        print('result: popt (ps):',popt*x_to_THz)
                #####################################################
                # and write to db
                # setpopt0 wird zur zeit nicht angewendet sonder alles einfach geschrieben, das ist gar nicht so schlecht
                #####################################################
                global pppx_add_db
                add_value_to_db_tmp(qpoint,mdsteps,peakidx+1,filename,popt_fq_THz,popt_lt_THz,popt_aa_THz)
                #sys.exit('hooh7')


            if setpopt0==True:
                #print filename.ljust(30),p(popt,x_to_THz,im,bounds=b)
                popt=np.ones(popt.shape)*0.000001

            save = False
            if mdsteps == args.mdstepstocalc_all.max(): save = True  # only when last mdstep
            if type(args.mdsteps) != bool: save = True

            #if args.write_fitted_ps == True:
            #if True:
            #print('type(args.mdsteps)',type(args.mdsteps))
            #print('save',save)
            #print('pppp7',popt*x_to_THz)
            if save: #args.verbose: # or save:
                add = ""
                #if space_fft_which == "space_fft_":
                if True:
                    add = "      "

                    print(printred("**** ==>"),"ps_fitted/ps"+qpointstr+"_"+str(mdsteps)+".dat_"+filename+add,ppplot(popt,x_to_THz,im,bounds=b))
            #if save and type(popt) != bool: # and setpopt0 != True:
            #if True:
            #print("mdsteps-->",mdsteps)
            #print("mdstepsmax-->",args.mdstepstocalc_all.max())
            if args.write_fitted_ps == True or type(args.mdsteps) != bool or mdsteps == args.mdstepstocalc_all.max():
                #yout = eval(take+take2)(x, *popt)
                yout = function(x, *popt)
                #print 'xind:',xind
                #print "MAX:",yout.max()
                #xoutt,youtt = ps_to_ps_for_writing(yout*psmax,xind,x_to_THz,eval(take+take2),x,popt,psmax)
                xoutt,youtt = ps_to_ps_for_writing(yout*psmax,xind,x_to_THz,function,x,popt,psmax,args=args)
                np.savetxt("ps_fitted/ps"+qpointstr+"_"+str(mdsteps)+".dat_"+filename,np.transpose(np.array([xoutt,youtt])))
                #print 'xoutt[:10]',xoutt[:10],xoutt[-10:]
                #print "^^^^^^^^^^^^^^^^^^ saved ^^^^^^^^^^^^^^^^^^^^^^^^^^^"
                #print
                #sys.exit()
            return



        verbose = True
        #peaks=[xind,x_to_THz,psmax]
        #print "x_to_THz:",x_to_THz,xind
        #print 'xindooo',xind

        kwargs = [xind, qpointstr, xl,yl,xa,ya,verbose,mdsteps,check_qpoint_for_crossing,qpoint,psmax,x_to_THz,space_fft_which]
        dofit = True
        if dofit:
            #sys.exit('exitsys123')
            take1="lorenz_68_Hauer2015_"
            if args.sigma0 == True:
                makefit(take = take1,strr="xa",kwargs=kwargs,sigma=0,b=bounds,args=args) # sigma 1 macht besseren fit bei hohen q's wo hintergrund hoch wird
            if args.sigma1 == True:
                makefit(take = take1,strr="xa",kwargs=kwargs,sigma=1,b=bounds,args=args) # sigma 1 macht besseren fit bei hohen q's wo hintergrund hoch wird
            #makefit(take = take1,strr="xa",kwargs=kwargs,sigma=2,b=bounds,args=args)
            #take2="lorenz_6886_Hauer2015_"   # whith bound it is not working for one peak
            #makefit(take = take2,strr="xa",kwargs=kwargs,sigma=0,b=bounds,args=args)
            #makefit(take = take2,strr="xa",kwargs=kwargs,sigma=1,b=bounds,args=args)   # der 6886 ist fuer kleine q zu hoch bei w=0THz
            #makefit(take = take2,strr="xa",kwargs=kwargs,sigma=2,b=bounds,args=args)

        # all with sigma=False just dont geht the maximum right
        #makefit(take = take1,strr="xl",kwargs=kwargs,sigma=False)
        #makefit(take = take1,strr="xl",kwargs=kwargs,sigma=False,b=bounds)
        #makefit(take = take1,strr="xa",kwargs=kwargs,sigma=False)
        #makefit(take = take1,strr="xa",kwargs=kwargs,sigma=False,b=bounds)
        # ohne gewichtutake1+ng (sigma): ist der fit besser der xindr einen peak fittet (fuer 16 0 0)

        #makefit(take = take1,strr="xl",kwargs=kwargs,sigma=1)
        #makefit(take = take1,strr="xl",kwargs=kwargs,sigma=1,b=bounds)
        #makefit(take = take1,strr="xa",kwargs=kwargs,sigma=1)
        # mit gewichtuntake1+g (sigma=1) ist der fit besser der xindide peaks beruecksichtigt (16 0 0)
        # fuer (16 0 0)take1+ bester fit: gewichtung von sigma=1xindnd 1p_1p
        #makefit(take = take1,strr="xl",kwargs=kwargs,sigma=2)
            #makefit(take = take1,strr="xl",kwargs=kwargs,sigma=2,b=bounds)
        #makefit(take = take1,strr="xa",kwargs=kwargs,sigma=2)
        #makefit(take = take1,strr="xa",kwargs=kwargs,sigma=2,b=bounds)
        # die gewichtung macht grundsaetzlich die linienbreite ixinder kleiner
        # speziell fuer (16 0 0) ist durch die starke assymmetrixinddie gewichtung von sigma=2 besser als sigma=1
        take2="lorenz_86_Hauer2015_"   # whith bound it is not working for one peak
        # DONT USE SOMEHOW CAN HIT HAFT THE LIFETIMEBUNDDARY AND NOT GO ON
        #makefit(take = take2,strr="xl",kwargs=kwargs,sigma=False,b=bounds)
        #makefit(take = take2,strr="xa",kwargs=kwargs,sigma=False,b=bounds)
            #makefit(take = take2,strr="xl",kwargs=kwargs,sigma=1    ,b=bounds)
            #makefit(take = take2,strr="xa",kwargs=kwargs,sigma=1    ,b=bounds)
        #makefit(take = take2,strr="xl",kwargs=kwargs,sigma=2    ,b=bounds)
        #makefit(take = take2,strr="xa",kwargs=kwargs,sigma=2    ,b=bounds)

        #@ # include later
        #@ makefit(take = take2,strr="xl",kwargs=kwargs,sigma=False)
        #@ makefit(take = take2,strr="xa",kwargs=kwargs,sigma=False)
        #@ makefit(take = take2,strr="xl",kwargs=kwargs,sigma=1    )
        #@ makefit(take = take2,strr="xa",kwargs=kwargs,sigma=1    )
        #@ makefit(take = take2,strr="xl",kwargs=kwargs,sigma=2    )
        #@ makefit(take = take2,strr="xa",kwargs=kwargs,sigma=2    )
        #@ # (16  0 0) best: lorenz_86_Hauer2015_1p_1p_xa_1        xind
        #@ # (14  0 0) best: lorenz_86_Hauer2015_1p_1p_xa_1        xind
        #@ # (48 16 0) best: lorenz_86_Hauer2015_1p_1p_xa_1        xind
        #@ # (14  0 0)  muss: sigma=1  (sigma=F ist schlecht)      xind
        #@ # (48 16 0)  muss: sigma=1  (sigma=F ist sehr!! schlechtxind
        #@ # ob bei (14 0 0) sigma=1 oder sigma=2 besser ist kann mxind nicht unterscheiden aber sigma=F ist schlecht


        return
    #print 'type(xind)',type(xind)
    if type(xind) == list or type(xind) == np.array or type(xind) == np.ndarray:
        fit_new_lorenz(ps,x_to_THz,xind,verbose=False,psmax=psmax,qpoint=qpoint,mdsteps=mdsteps,args=args,space_fft_which=space_fft_which)


    #@ if check_qpoint_for_crossing == True:
    #@     sigmagood = sigmamin
    #@     fqgood = fqmin
    #@     ltgood = ltmin
    #@     ltgoodmin,ltgoodmax = ltminmin,ltminmax
    #@ #############################################################
    #@ # write here the smoothigfunc, further down the lifetimes
    #@ #############################################################
    #@ finaldata_min  = [sigmamin, fqmin,ltmin,ltminmin,ltminmax,int(mdsteps)]
    #@ finaldata_good = [sigmagood, fqgood,ltgood,ltgoodmin,ltgoodmax,int(mdsteps)]
    #@ finaldata_lorenz = [sigmagood, freq_lorenz,lifetime_lorenz,ltgoodmin,ltgoodmax,int(mdsteps)]
    #@ #print "finaldata_min :",finaldata_min
    #@ #print "finaldata_good:",finaldata_good
    #@ if mdsteps == mdstepstocalc.max():
    #@     filenameoutlorenz = filenameout+".smooth_out_lorenz.dat"

    #@     # adapt to 2d which is smaller
    #@     if not args.notverbose:
    #@         print printblue("WRITING3: ps_smooth/"+filenameout+" and ps_smooth/"+filenameoutgood+" and ps_smooth/"+filenameoutlorenz)
    #@     o = x[::2]/2
    #@     on = x[:len(fit_lorenz)]
    #@     parnpsavetxt("ps_smooth/"+filenameoutlorenz,   np.transpose(np.array([on,fit_lorenz]))[:xmaxwritemin],fmt='%.1f %.5f')

    #@ if type(return_dict) != bool:
    #@     return_dict[mdsteps] = [finaldata_good,finaldata_min,filenameout,finaldata_lorenz]
    #@ return finaldata_good,finaldata_min,filenameout,finaldata_lorenz
    return


def smoothing_power_spectrum(power_spectrum,sigma,scale=True):
    ''' smoothes a gaussian by sigma '''
    steps = power_spectrum.shape[0]                 # 1000
    nu = np.arange(0,steps)/float(steps)            # [ 0., 0.001, ..., 0.999 ]
    kern = np.exp(-(np.mod(nu+0.5,1)-0.5)**2/(2*sigma**2))  # Gaussian distribution in 1-D
    kern = kern/np.sum(kern)
    ps = np.real(ifft(fft(power_spectrum)*fft(kern)))
    if scale:
        ps=ps/np.max(ps)*.9;
    return ps

def get_uncertainty_smoothed_powerspectrum(power_spectrum_smooth,dt,sigma = "",args=False,stringadd='',verbose=False,color='blue'):
    ''' - dt is the sampled timestep in femtoseconds
        - power_spectrum_smooth is the full power spectrum having both peaks (left and right)
        - sigma is only for plotting to screen
    '''
    verbose = args.verbose
    data = power_spectrum_smooth
    #np.savetxt("ka"+str(sigm)+".dat",data)  #dies schein schon oversmoothed zu sein
    #np.savetxt("ka.dat",data)  #dies schein schon oversmoothed zu sein
    #sys.exit()

    xvalues=data.shape[0]
    #print "xvalues all:",xvalues
    datahalf = data[0:int(xvalues/2)]   # only take one peak
    #np.savetxt("kb.dat",datahalf)
    #print "datahalf:",datahalf

    def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return array[idx]


    # x index of absolute maximum
    #np.savetxt("data",data)
    #sys.exit()
    x_ind_max = np.where(datahalf == datahalf.max())[0][0]  # 963  (9.63 THz)
    y_val_max = datahalf[x_ind_max]
    #if args.verbose:
    #    print "x_ind_max 77:",x_ind_max

    if x_ind_max == 0:
        return 0 ,0, 0 ,xvalues, 0,0,0,99999999999999999999999999

    def printit(x_ind_max,y_val_max,x_ind_leftdata_min,y_val_leftdata_min,y_val_leftdata_max,y_lifetime_max,y_lifetime_min):
        print("############################################################")
        print("# analyse the hight of the background")
        print("#                          .---------------------- y_val_max:",y_val_max)
        print("#                         .|.")
        print("#                        . | .")
        print("#                       .  |  .")
        print("#                     .    |   .")
        print("#                   .      |    .")
        print("#                 .        |     .")
        print("#               .          |      .")
        print("#.---------- . ------------|------ .------------- y_val_leftdata_max:",y_val_leftdata_max)
        print("#  .      .                |        .")
        print("#     .------------------------------- . -------- y_val_leftdata_min:",y_val_leftdata_min,"(therefore linewidth measured at y = y_lifetime_max:",y_lifetime_max," y_lifetime_min:",y_lifetime_min,")")
        print("#     |                    |               .  .   .")
        print("############################################################")
        print("#     x--------------------|----------------------- x_ind_leftdata_min:",x_ind_leftdata_min)
        print("#")
        print("#                          x                        x_ind_max:",x_ind_max)
        print("#x                                                  x_val_leftdata_max:")
        print("############################################################")
        return

    y_val_leftdata_min = datahalf[:x_ind_max].min()
    x_ind_leftdata_min = np.where(datahalf[:x_ind_max] == y_val_leftdata_min)[0][0]

    if x_ind_leftdata_min == 0:
        y_val_leftdata_max = datahalf[0]
    else:
        y_val_leftdata_max = datahalf[:x_ind_leftdata_min].max()

    #############################################################
    ## maximal linewidths if background not substracted
    ## NOT substracting the background is not something the experimentalists can do!!!
    ## THEREFORE: get another max .... TODO
    #############################################################
    y_lifetime_max              = datahalf.max()/2.
    y_lifetime_max              = datahalf.max()/2.+y_val_leftdata_min/2  # Therefore changed definition
    y_ind_max_over2_left_max    = find_nearest(datahalf[0:x_ind_max], y_lifetime_max)
    y_ind_max_over2_right_max   = find_nearest(datahalf[x_ind_max:-1], y_lifetime_max)
    x_ind_left_over2_max        = np.where(datahalf[0:x_ind_max] == y_ind_max_over2_left_max)[0][0]
    x_ind_right_over2_max       = np.where(datahalf == y_ind_max_over2_right_max)[0][0]

    #############################################################
    ## minimal linewidths if maximal background is substracted
    #############################################################
    y_lifetime_min = datahalf.max()/2.+y_val_leftdata_max/2.
    y_ind_max_over2_left_min  = find_nearest(datahalf[0:x_ind_max], y_lifetime_min)
    y_ind_max_over2_right_min = find_nearest(datahalf[x_ind_max:-1], y_lifetime_min)
    x_ind_left_over2_min = np.where(datahalf[0:x_ind_max] == y_ind_max_over2_left_min)[0][0]
    x_ind_right_over2_min = np.where(datahalf == y_ind_max_over2_right_min)[0][0]

    if args.verbose > 3:
        printit(x_ind_max,y_val_max,x_ind_leftdata_min,y_val_leftdata_min,y_val_leftdata_max, y_lifetime_max,y_lifetime_min)


    #############################################################
    ## "normal/intermediate"  linewidths if background is substracted
    #############################################################
    y_lifetime = datahalf.max()/2.+y_val_leftdata_min/2.
    #print "datahalf.max()/2     :",datahalf.max()/2.
    #print "y_val_leftdata_min   :",y_val_leftdata_min
    #print "y_lifetime           :",y_lifetime
    #sys.exit()

    y_lifetime_third = datahalf.max()/3.
        #return freq,lifetime, lifetimemin,xvalues, r+l,ll+lr+rl+rr, x_ind_right_over1000

    #############################################################
    # dis ist nur schoen (und korrekt) bei einem relativ glatten gaussian
    #############################################################
    # x index of left
    #print "x_ind_max:",x_ind_max
    #print "datahalf[0:x_ind_max]",datahalf[0:x_ind_max]
    #print "datahalf.max()/2.:",datahalf.max()/2.
    y_ind_max_over2_left  = find_nearest(datahalf[0:x_ind_max], y_lifetime)
    y_ind_max_over3_left  = find_nearest(datahalf[0:x_ind_max], y_lifetime_third)
    #print "y_ind_max_over2_left",y_ind_max_over2_left
    #print "y_ind_max_over3_left",y_ind_max_over3_left
    x_ind_left_over2 = np.where(datahalf[0:x_ind_max] == y_ind_max_over2_left)[0][0]
    x_ind_left_over3 = np.where(datahalf[0:x_ind_max] == y_ind_max_over3_left)[0][0]

    y_ind_max_over2_right = find_nearest(datahalf[x_ind_max:-1], y_lifetime)
    y_ind_max_over3_right = find_nearest(datahalf[x_ind_max:-1], y_lifetime_third)

    y_ind_max_over2_right = find_nearest(datahalf[x_ind_max:-1], y_lifetime)
    y_ind_max_over3_right = find_nearest(datahalf[x_ind_max:-1], y_lifetime_third)


    #print "y_ind_max_over2_right",y_ind_max_over2_right
    x_ind_right_over2 = np.where(datahalf == y_ind_max_over2_right)[0][0]
    x_ind_right_over3 = np.where(datahalf == y_ind_max_over3_right)[0][0]
    #sys.exit()
    #print "x_ind_left_over2",x_ind_left_over2,"y:",y_ind_max_over2_left
    #print "x_ind_right_over2",x_ind_right_over2, "y:",y_ind_max_over2_right


    l = len(np.where(np.diff(datahalf[x_ind_left_over2:x_ind_max])<0)[0])    # using x_ind_left_over2 insted of x_ind_left_over3 makes a TA2 realistic errorbar
    #print l
    #print "bbb"
    r = len(np.where(np.diff(datahalf[x_ind_max:x_ind_right_over2])>0)[0])   # using x_ind_left_over2 insted of x_ind_left_over3 makes a TA2 realistic errorbar
    #print "r+l:",r+l
    #sys.exit()

    #############################################################
    # checke die ableitung des peaks
    #############################################################
    ###############################
    # erstmal die linke haelfte
    ###############################
    #print x_ind_left_over3,x_ind_max,x_ind_right_over3
    ll = lr = 0.
    #print "x_ind_max:",x_ind_max
    if args.verbose > 3:
        print("x_ind_left_over2,y_ind_max_over2_left:",x_ind_left_over2,y_ind_max_over2_left)
        print("x_ind_left_over3,y_ind_max_over3_left:",x_ind_left_over3,y_ind_max_over3_left)

    #if (x_ind_max - x_ind_left_over3) > 1:
    #    datatake = datahalf[x_ind_left_over3:x_ind_max]
    if (x_ind_max - x_ind_left_over2) > 1:
        datatake = datahalf[x_ind_left_over2:x_ind_max]
        ableitung1 = np.diff(datatake)  # 95:[97] funktioniert k --> ableitung ist nur ein wert
        #print ableitung1
        ableitung1_ind_max = np.where(ableitung1 == ableitung1.max())[0][0]
        ableitung1l = ableitung1[:ableitung1_ind_max]
        ableitung1r = ableitung1[ableitung1_ind_max:]
        ll = len(np.where(np.diff(ableitung1l)<0)[0])
        lr = len(np.where(np.diff(ableitung1r)>0)[0])

        #np.savetxt("ka.dat",datahalf)
        #np.savetxt("kb.dat",datatake)
        #np.savetxt("kc.dat",ableitung1)
        #np.savetxt("kd.dat",ableitung1l)
        #np.savetxt("ke.dat",ableitung1r)
        #sys.exit()
    rl = rr = 0.
    #if (x_ind_right_over3-x_ind_max) > 1:
    #    datatake = datahalf[x_ind_max:x_ind_right_over3]
    if (x_ind_right_over2-x_ind_max) > 1:
        datatake = datahalf[x_ind_max:x_ind_right_over2]
        ableitung1 = np.diff(datatake)  # 95:[97] funktioniert k --> ableitung ist nur ein wert
        #np.savetxt("ka.dat",datatake)
        #np.savetxt("kb.dat",ableitung1)
        #print ableitung1
        ableitung1_ind_max = np.where(ableitung1 == ableitung1.min())[0][0]
        ableitung1l = ableitung1[:ableitung1_ind_max]
        ableitung1r = ableitung1[ableitung1_ind_max:]

        #np.savetxt("kc.dat",ableitung1l)
        #np.savetxt("kd.dat",ableitung1r)


        rl = len(np.where(np.diff(ableitung1l)>0)[0])
        rr = len(np.where(np.diff(ableitung1r)<0)[0])
        #print rl,rr
    #rl = rr = 0.
    #print ll,rr

    #############################################################
    # dis ist zuverlaessiger bei viel zu dicht gesampleten
    #############################################################
    #print
    #print "bbb:",y_lifetime
    #if len(np.where(datahalf>y_lifetime)[0]) == 0:
    #    np.savetxt("kkk.dat",datahalf)
    #print "ccc:",len(datahalf)
    #print "kkk:",len(np.where(datahalf>y_lifetime)[0])
    #if len(np.where(datahalf>y_lifetime)[0]) == 0:
    #    np.savetxt("kkk.dat",datahalf)
    #    print "now 0"

    #print "aaa:",np.where(datahalf>y_lifetime)[0][0]
    #np.savetxt("datahalf",datahalf)
    #print "y_lifetime:",y_lifetime
    absmin = np.where(datahalf>y_lifetime)[0][0]     # nur die spitze des gaussians (linke seite)
    absmax = np.where(datahalf>y_lifetime)[0][-1]    # nur die spitze des gaussians (rechte seite)

    if len(np.where(datahalf>y_lifetime)[0]) == 0:
        np.savetxt("kkk.dat",datahalf)
        os._exit(1)
        sys.exit("now 0")

    #mult = dt*10 / float(xvalues)  ## Wrong!  see example 6x6x6sc 6_0_0 qpoint between dt={10,20}
    mult = 1000./(args.dt*float(xvalues))
    #print "fx:",float(xvalues)
    #print "dh:",len(datahalf)
    #print "mult:",mult
    #print "x_ind_max:",x_ind_max
    freq = x_ind_max * mult
    lifetime     = (x_ind_right_over2 - x_ind_left_over2) * mult
    lifetimemax  = (x_ind_right_over2_max - x_ind_left_over2_max) * mult
    diff_lifetimemax = np.abs(lifetimemax - lifetime)
    lifetimemin0 = (x_ind_right_over2_min - x_ind_left_over2_min) * mult
    lifetimemin1 = (absmax - absmin) * mult          # top of gaussian
    lifetimemin = np.array([lifetimemin0, lifetimemin1]).min()
    diff_lifetimemin = np.abs(lifetimemin - lifetime)





    ######################################################################################
    # now new fitting approach to lorenzian distribution
    ######################################################################################
    lifetimemax = lifetime + diff_lifetimemax
    lifetimemin = lifetime - diff_lifetimemin

    #print lifetimemin,lifetime, lifetimemax
    #vor_freq = x_ind_max/float(xvalues)
    #vor_lifetime = (x_ind_right_over2 - x_ind_left_over2)/float(xvalues)
    #vor_lifetimemin = (absmax - absmin)/float(xvalues)
    ##print "vor_freq:",vor_freq
    ###print "vor_lifetime:",vor_lifetime
    ###print "vor_lifetimemin:",vor_lifetimemin
    #print "Smoothed function:"
    #if args.verbose and verbose:
    if args.verbose > 2:
        if stringadd[0] == 'min':
            stringadd[0] = 'min '
        #ljust adds zeros (0) to the end
        #ka = printblue(stringadd)
        if color == 'blue':
            print("Freq "+printblue(stringadd[0])+" "+stringadd[1].ljust(10)+":",str(round(freq,2)).ljust(4,'0'),"THz; width: ",str(round(lifetime,2)).ljust(4,'0'),"THz", "||", x_ind_left_over2,x_ind_max,x_ind_right_over2,"||", "sigma:",str(sigma),"r+l",r+l,"rr,ll",ll+lr+rl+rr,"ll:",ll,"rr:",rr,"lr:",lr,"rl:",rl)
        if color == 'red':
            print(printred("Freq ")+printred(stringadd[0])+" "+stringadd[1].ljust(10)+":",str(round(freq,2)).ljust(4,'0'),printred("THz; width: "),printred(str(round(lifetime,2)).ljust(4,'0')),"THz", "||", x_ind_left_over2,x_ind_max,x_ind_right_over2,"||", "sigma:",str(sigma),"r+l",r+l,"rr,ll",ll+lr+rl+rr,"ll:",ll,"rr:",rr,"lr:",lr,"rl:",rl)
            #"lifetime min:",lifetimemin,\
    #print "indizes:",x_ind_left_over2,"<-",x_ind_max,"->",x_ind_right_over2 #,"from" ,xvalues,"xvalues"

    # x cutoff when writing files to harddisk
    y_ind_max_over1000_right = find_nearest(datahalf[x_ind_max:-1], datahalf.max()/1000.)
    x_ind_right_over1000 = np.where(datahalf == y_ind_max_over1000_right)[0][0]
    #print "xxx:",x_ind_right_over1000
    #print "--->>",freq,lifetime, lifetimemin, lifetimemax ,xvalues, r+l,ll+lr+rl+rr, x_ind_right_over1000
    return freq,lifetime, lifetimemin, lifetimemax ,xvalues, r+l,ll+lr+rl+rr, x_ind_right_over1000


def fit_lorenzian_to_powerspektrum_old_fit_one_peak(power_spectrum_full_left_right_peak,dt,freqestimate=False,lifetimeestimate=False,stringadd="",check_qpoint_for_crossing = False,verbose=False,printwarning=False):
    vb=verbose
    vbresult=False
    #vb = True
    #np.savetxt("tmp.dat",power_spectrum_full_left_right_peak)

    if type(freqestimate) == bool:
        freq = 4.  # THz
        sys.exit("this function currently needs a relatively good estimate!")
    else:
        freq=freqestimate
    if type(lifetimeestimate) == bool:
        lifetime = 1.
        sys.exit("this function currently needs a relatively good estimate!")
    else:
        lifetime=lifetimeestimate

    #np.savetxt("ka"+str(sigm)+".dat",data)  #dies schein schon oversmoothed zu sein
    #np.savetxt("ka.dat",data)  #dies schein schon oversmoothed zu sein

    xvalues=power_spectrum_full_left_right_peak.shape[0]
    mult = 1000./(dt*float(xvalues))
    #print "xvalues all:",xvalues
    y = power_spectrum_full_left_right_peak[0:xvalues/2]   # only take one peak
    x=np.arange(len(y));
    #print "y:",y
    #print "x:",x


    if vb:
        print()
        print("################### crossing? ",check_qpoint_for_crossing)
    #if vb:
    #    print "init freq    :",freq,freq/mult
    #    print "init lifetime:",lifetime,lifetime/mult
    # wenn wir mal doch den berech zum fitten einschraenken wollen
    #bereich_links  = int(p0_freq) - 10*int(p0_life)
    #bereich_rechts = int(p0_freq) + 10*int(p0_life)
    #if bereich_links < 0:
    #    bereich_links = 0
    #if bereich_rechts > len(y):
    #    bereich_rechts = len(y)
    #np.savetxt("y.dat",y)

    maxfev=200000
    conv_crit = 10**-4
    p0_freq = np.abs(freq/mult)
    p0_hoeh = 100000000000
    p0_hoeh = 1
    if lifetime >= 3.:
        lifetime = 3.
    p0_life = np.abs(lifetime/mult)
    p0_life_max_crossing = np.abs(3.05/mult)

    p0in=np.array([p0_freq,p0_hoeh,p0_life])

    #############################################################
    # erster grober fit um erstmal popt0[1] zu bekommen
    #############################################################
    fd1 = 0.05  # 5 %
    ld1 = 0.05  # 5 %
    fd1 = 0.01  # 5 %
    ld1 = 0.01  # 5 %
    bounds1=((p0_freq - fd1*p0_freq, 0,p0_life-ld1*p0_life),(p0_freq + fd1*p0_freq,np.inf,p0_life+ld1*p0_life))
    fd2 = 0.3
    ld2 = 0.9
    bounds=((p0_freq - fd2*p0_freq, 0,p0_life-ld2*p0_life),(p0_freq + fd2*p0_freq,np.inf,p0_life+ld2*p0_life))
    if check_qpoint_for_crossing:
        bounds=((p0_freq - fd2*p0_freq, 0,p0_life-ld2*p0_life),(p0_freq + fd2*p0_freq,np.inf,p0_life_max_crossing))

    if type(freqestimate) == bool or type(lifetimeestimate) == bool:
        bounds=((0, 0,0),(20/mult,np.inf,4))

    if vb:
        #print "-------------- start -----------"
        print("---> poptin  ->",p0in*mult, len(y),lifetime,"-->",stringadd) #,"sigma:",sigma,"str:",stringadd #, "mult:",mult
    # THIS IS ONLY TO GET A FIRST ESTIMATE
    popt0, pcov0 = curve_fit(lorenzbesserwiki, x, y,p0=p0in,bounds=bounds1)
    yfit=lorenzbesserwiki(x, popt0[0], popt0[1], popt0[2])
    if y.max()/yfit.max() > 10. or y.max()/yfit.max() < 0.1:
        for i in np.arange(20):
            p0_hoeh = 10**i
            #print "i,p0;",i,p0_hoeh
            popt0, pcov0 = curve_fit(lorenzbesserwiki, x, y,p0=np.array([p0in[0],p0_hoeh,p0in[2]]),bounds=bounds1)
            yfit=lorenzbesserwiki(x, popt0[0], popt0[1], popt0[2])
            if vb:
                print("i,p0,hoeh:",i,p0_hoeh,y.max()/yfit.max())
            if y.max()/yfit.max() <= 10. and y.max()/yfit.max() > 0.1:
                break


    popt = popt0
    popt[2]=abs(popt[2])
    sigmaa=lorenzbesserwiki(x, popt[0], popt[1], popt[2])
    if vb:
        print("---> poptout ->",popt*mult,y.max()/yfit.max())

    #############################################################
    # genauer fit
    #############################################################
    for i in np.arange(10):
        if i == 0:
            poptold = popt
        if vb:
            print("---> poptin :",i,popt*mult,popt)
        sigmaa=lorenzbesserwiki(x, popt[0], popt[1], popt[2])
        #sigmaa=lorenzbesserwiki(x, popt[0], popt[1], popt[2])**2
        #sigmaa=lorenzbesserwiki(x, popt[0], popt[1], popt[2])**4
        #sigmaa=lorenzbesserwiki(x, popt[0], popt[1], popt[2])**10
        #sigmaa=lorenzbesserwiki(x, popt[0], popt[1], popt[2])**100
        #sigmaa=lorenzbesserwiki(x, popt[0], popt[1], popt[2])**1000
        #sigmaa=lorenzbesserwiki(x, popt[0], popt[1], popt[2])**100
        #print "write !!!!!!!! sigma"
        #np.savetxt("sig.dat",1./sigmaa)
        popt, pcov = curve_fit(lorenzbesserwiki, x, y,sigma=1./sigmaa,p0=popt,bounds=bounds) #,maxfev=maxfev) #,p0=[813,10**9,1])
        popt[2]=abs(popt[2])
        #perr = np.sqrt(np.diag(pcov))
        #print "perr:",perr
        #if vb:
        #    print "---> poptout:",i,popt*mult,abs(1. - popt[0]/poptold[0]),abs(1. - popt[2]/poptold[2])
        #if abs(1. - popt[0]/poptold[0]) <= 10**-4:
        #    print "1<2?",1,np.abs(1. - popt[0]/poptold[0])
        #    print "1<2?",2,10**-4
        #    sys.exit()
        if abs(1. - popt[0]/poptold[0]) <= conv_crit and abs(1. - popt[2]/poptold[2]) <= conv_crit:
            if vb:
                print("done popt:",popt*mult,popt)
            break
        poptold = popt
        #if i == 4:
        #    sys.exit()
    if vb:
        print("---> poptout->:",popt*mult,abs(1. - popt[0]/poptold[0]),abs(1. - popt[2]/poptold[2]))
        #print "freq:",freq,"freq_new:",popt[0]*mult
        #print "lifetime:",lifetime,"lifetimenew:",popt[2]*mult
    #if vb:
    #    print "write kx.dat using popt:",popt*mult
    #    np.savetxt
    #    np.savetxt("ki"+str(len(power_spectrum_full_left_right_peak))+".dat",np.transpose((x, y)))
    #    np.savetxt("kx"+str(len(power_spectrum_full_left_right_peak))+".dat",np.transpose((x, lorenzbesserwiki(x, popt[0], popt[1], popt[2]))))

    freq = np.abs(popt[0])*mult
    lifetime = popt[2]*mult
    ################################################
    # checks
    ################################################
    if np.abs(1. - freqestimate/freq) > 0.1:
        if printwarning:
            print("freqestimate:",freqestimate)
            print("freq:",freq)
            print("res->",freqestimate/freq,np.abs(1. - freqestimate/freq))
            print_warning("freq/freqestimate seems wrong")
    if np.abs(1. - lifetimeestimate/lifetime) > 2:
        if printwarning:
            print("lifetimeestimate:",lifetimeestimate)
            print("lifetime:",lifetime)
            print("res->",freqestimate/freq,np.abs(1. - freqestimate/freq))
            print_warning("lifetime/lifetimeestimate seems wrong")
    yfit=lorenzbesserwiki(x, popt[0], popt[1], popt[2])
    if y.max()/yfit.max() > 10. or y.max()/yfit.max() < 0.1:
        if printwarning:
            print_warning("y seems wrong")
    #os._exit(1)

    #lifetimemax = lifetime + diff_lifetimemax
    #lifetimemin = lifetime - diff_lifetimemin
    if vb:
        print("################### crossing? ",check_qpoint_for_crossing)
        #print "max rauschen",np.where(y == y.max())[0][0]*mult,np.where(y == y.max())[0][0]
    #np.savetxt("ppp.dat",power_spectrum_full_left_right_peak)
    od = 4
    if vb:
        print("++++++++++++++++",round(freq,od),round(lifetime,od))
    return round(freq,od),round(lifetime,od),yfit


def fit_lorenzian_to_powerspektrum(power_spectrum_full_left_right_peak,dt,freqestimate=False,lifetimeestimate=False,stringadd="",check_qpoint_for_crossing = False,verbose=False,printwarning=False,take='lorenzbesserwiki'):
    #power_spectrum_full_left_right_peak = power_spectrum_full_left_right_peak/power_spectrum_full_left_right_peak.max()
    vb=verbose
    vbresult=False
    vb = False
    vbb = False
    psmax = power_spectrum_full_left_right_peak.max()

    if type(freqestimate) == bool:
        freq = 4.  # THz
        sys.exit("this function currently needs a relatively good estimate!")
    else:
        freq=freqestimate

    if type(lifetimeestimate) == bool:
        lifetime = 1.
        sys.exit("this function currently needs a relatively good estimate!")
    else:
        lifetime=lifetimeestimate

    y_all = power_spectrum_full_left_right_peak
    xvalues=len(y_all)
    x_all = np.arange(xvalues)

    mult = 1000./(dt*float(xvalues))
    y = power_spectrum_full_left_right_peak[0:xvalues/2]   # only take one peak
    x=np.arange(len(y));

    xr = x+x[-1]+1
    yr = power_spectrum_full_left_right_peak[xvalues/2:]

    if vb:
        print("x_all:",len(x_all),x_all)
        print("x    :",len(x),x)
        print("xr   :",len(xr),xr)

        print()
        print("y_all:",len(y_all),y_all)
        print("y    :",len(y),y)
        print("yr   :",len(yr),yr)

        print()
        print("################### crossing? ",check_qpoint_for_crossing)

    maxfev=200000
    conv_crit = 10**-4
    p0_freq = np.abs(freq/mult)
    p0_hoeh = 100000000000
    p0_hoeh = 1
    if lifetime >= 3.:
        lifetime = 3.
    p0_life = np.abs(lifetime/mult)
    p0_life_max_crossing = np.abs(3.05/mult)

    p0in=np.array([p0_freq,p0_hoeh,p0_life])

    #############################################################################
    # 1. erster grober (nur linker peak) fit um erstmal popt0[1] zu bekommen
    #############################################################################
    fd1 = 0.05  # 5 %
    ld1 = 0.05  # 5 %
    fd1 = 0.01  # 5 %
    ld1 = 0.01  # 5 %
    bounds1=((p0_freq - fd1*p0_freq, 0,p0_life-ld1*p0_life),(p0_freq + fd1*p0_freq,np.inf,p0_life+ld1*p0_life))
    fd2 = 0.3
    ld2 = 0.9
    bounds=((p0_freq - fd2*p0_freq, 0,p0_life-ld2*p0_life),(p0_freq + fd2*p0_freq,np.inf,p0_life+ld2*p0_life))
    if check_qpoint_for_crossing:
        bounds=((p0_freq - fd2*p0_freq, 0,p0_life-ld2*p0_life),(p0_freq + fd2*p0_freq,np.inf,p0_life_max_crossing))

    if type(freqestimate) == bool or type(lifetimeestimate) == bool:
        print("freqsestimate == bool or lifetimeestimate == bool --> no estimates!")
        bounds=((0, 0,0),(20/mult,np.inf,4))

    if vb:
        #print "-------------- start -----------"
        print("---> poptin_aa  ->",p0in*mult, len(y),lifetime,"-->",stringadd) #,"sigma:",sigma,"str:",stringadd #, "mult:",mult
    # THIS IS ONLY TO GET A FIRST ESTIMATE
    popt0, pcov0 = curve_fit(lorenzbesserwiki, x, y,p0=p0in,bounds=bounds1)
    yfit=lorenzbesserwiki(x, popt0[0], popt0[1], popt0[2])
    if y.max()/yfit.max() > 10. or y.max()/yfit.max() < 0.1:
        for i in np.arange(18): # dont go above 18, otherwise overflow at mac
            p0_hoeh = 10**i
            #print "i,p0;",i,p0_hoeh
            popt0, pcov0 = curve_fit(lorenzbesserwiki, x, y,p0=np.array([p0in[0],p0_hoeh,p0in[2]]),bounds=bounds1)
            yfit=lorenzbesserwiki(x, popt0[0], popt0[1], popt0[2])
            if vb:
                print("i,p0,hoeh:",i,p0_hoeh,y.max()/yfit.max())
            if y.max()/yfit.max() <= 10. and y.max()/yfit.max() > 0.1:
                break


    popt = popt0
    popt[2]=abs(popt[2])
    sigmaa=lorenzbesserwiki(x, popt[0], popt[1], popt[2])
    if vb:
        print("---> (0) poptout ->",popt*mult,y.max()/yfit.max())
    #############################################################
    # 2. genauer fit (nur linker peak)
    #############################################################
    for i in np.arange(10):
        if i == 0:
            poptold = popt
        if vb:
            print("---> poptin_bb i:",i,"popt*mult:",popt*mult,"popt:",popt)
            print("bounds:",bounds)
        sigmaa=lorenzbesserwiki(x, popt[0], popt[1], popt[2])
        #sigmaa=lorenzbesserwiki(x, popt[0], popt[1], popt[2])**2
        #sigmaa=lorenzbesserwiki(x, popt[0], popt[1], popt[2])**4
        #sigmaa=lorenzbesserwiki(x, popt[0], popt[1], popt[2])**10
        #sigmaa=lorenzbesserwiki(x, popt[0], popt[1], popt[2])**100
        #sigmaa=lorenzbesserwiki(x, popt[0], popt[1], popt[2])**1000
        #sigmaa=lorenzbesserwiki(x, popt[0], popt[1], popt[2])**100
        #print "write !!!!!!!! sigma"
        #np.savetxt("sig.dat",1./sigmaa)
        popt, pcov = curve_fit(lorenzbesserwiki, x, y,sigma=1./sigmaa,p0=popt,bounds=bounds) #,maxfev=maxfev) #,p0=[813,10**9,1])
        popt[2]=abs(popt[2])
        #perr = np.sqrt(np.diag(pcov))
        #print "perr:",perr
        #if vb:
        #    print "---> poptout:",i,popt*mult,abs(1. - popt[0]/poptold[0]),abs(1. - popt[2]/poptold[2])
        #if abs(1. - popt[0]/poptold[0]) <= 10**-4:
        #    print "1<2?",1,np.abs(1. - popt[0]/poptold[0])
        #    print "1<2?",2,10**-4
        #    sys.exit()
        if abs(1. - popt[0]/poptold[0]) <= conv_crit and abs(1. - popt[2]/poptold[2]) <= conv_crit:
            if vb:
                print("done popt (genauer fit):",popt*mult,popt)
            break
        poptold = popt
        #if i == 4:
        #    sys.exit()
    if vb:
        print("---> (1) poptout->:",popt*mult,abs(1. - popt[0]/poptold[0]),abs(1. - popt[2]/poptold[2]))

    freq = np.abs(popt[0])*mult
    lifetime = popt[2]*mult
    popt_one_peak = popt


    ################################################
    # 3. check fit (nur linker peak)
    ################################################
    if np.abs(1. - freqestimate/freq) > 0.1:
        if printwarning:
            print("freqestimate:",freqestimate)
            print("freq:",freq)
            print("res->",freqestimate/freq,np.abs(1. - freqestimate/freq))
            print_warning("freq/freqestimate seems wrong")
    if np.abs(1. - lifetimeestimate/lifetime) > 2:
        if printwarning:
            print("lifetimeestimate:",lifetimeestimate)
            print("lifetime:",lifetime)
            print("res->",freqestimate/freq,np.abs(1. - freqestimate/freq))
            print_warning("lifetime/lifetimeestimate seems wrong")
    yfit=lorenzbesserwiki(x, popt[0], popt[1], popt[2])
    if y.max()/yfit.max() > 10. or y.max()/yfit.max() < 0.1:
        if printwarning:
            print_warning("y seems wrong")

    od = 4
    if vbb:
        print("+++++++++one_peak+++++++ freq / lifetime",round(freq,od),round(lifetime,od))


    ################################################
    # 4. erste tests 2 fitten 2 peaks
    ################################################
    if False: # only for debugging / testing
        #yfit_all=lorenzbesserwiki(np.arange(xvalues), popt[0], popt[1], popt[2])
        x_all_lr =np.arange(2*xvalues)-xvalues
        x_all_lr2=np.arange(4*xvalues)-2*xvalues
        yfit_all=lorenzbesserwiki(x_all_lr, popt[0], popt[1], popt[2])
        print("####################### checks #############################")
        print("x:",x)
        print("x_all_lr:",x_all_lr)
        delta = xvalues - popt[0]
        #yfit_rechts=lorenzbesserwiki(2*popt[0]-x, popt[0], popt[1], popt[2])
        def lorenzbesserwiki_ne(w,w0,a,b):
            ''' w: frequenzabhaengigkeit
                w0: die frequenz
                a hoehe  (also called the plasma frequency)
                b breite (gamma)
                took this'''
                #
                # S=4*w0^2*b/pi./((w.^2-w0.^2).^2+4*b.^2*w.^2);
            #w = -w + 18505+18500
            w = -w + 37500 # (len(xvalues))
            return (a**2)/((w**2-w0**2)**2.+w**2*b**2.) # finite at w=0


        sigmaa=lorenzbesserwiki(x, popt[0], popt[1], popt[2])
        popt, pcov = curve_fit(lorenzbesserwiki, x, y,sigma=1./sigmaa) #,p0=popt,bounds=bounds) #,maxfev=maxfev) #,p0=[813,10**9,1])
        print("final no bound (0):",popt)
        print("final no bound (0):",popt*mult)



        print("xx:",2*popt[0]-x_all_lr)
        #np.savetxt("tmp_in.dat",power_spectrum_full_left_right_peak)
        #np.savetxt("tmp_l.dat",np.transpose((x,y)))
        #np.savetxt("tmp_r.dat",np.transpose((xr,yr)))
        #yfit_l=lorenzbesserwiki(x_all_lr, popt[0], popt[1], popt[2])
        #np.savetxt("tmp_l_fit.dat",np.transpose((x_all_lr2,yfit_l)))
        yfit_l=lorenzbesserwiki(x_all_lr2, popt[0], popt[1], popt[2])
        np.savetxt("tmp_l_fit.dat",np.transpose((x_all_lr2,yfit_l)))


        #yfit_t=lorenzbesserwiki(x_all_lr2, popt[0]+10000, popt[1], popt[2])  # wrong symmetric
        #np.savetxt("tmp_t_fit.dat",np.transpose((x_all_lr2,yfit_t)))

        #np.savetxt("tmp_out_all.dat",np.transpose((x_all_lr,yfit_all)))

        ######## working
        #yfit_rechts=lorenzbesserwiki(2*popt[0]-x_all_lr2-3500, popt[0], popt[1], popt[2])
        #np.savetxt("tmp_t_fit.dat",np.transpose((x_all_lr2+len(x),yfit_rechts)))
        print("lennnnn:",len(x_all))
        yfit_r=lorenzbesserwiki_ne(x_all_lr2, popt[0], popt[1], popt[2])
        np.savetxt("tmp_r_fit.dat",np.transpose((x_all_lr2,yfit_r)))

    ################################################
    # 5. fitten 2 peaks
    ################################################
    def lorenzbesserwiki_lr(w,w0,a,b):
        ''' w: frequenzabhaengigkeit
            w0: die frequenz
            a hoehe  (also called the plasma frequency)
            b breite (gamma)
            took this'''
        lenw = len(w)-1.
        return (a**2)/((w**2-w0**2)**2.+w**2*b**2.) + (a**2)/(((-w+lenw)**2-w0**2)**2.+(-w+lenw)**2*b**2.)

    #np.savetxt("tmp_in.dat",power_spectrum_full_left_right_peak)
    #yfit_l= lorenzbesserwiki(x_all, popt[0], popt[1], popt[2])
    #np.savetxt("tmp_l.dat",np.transpose((x_all,yfit_l)))
    for i in np.arange(10):
        if i == 0:
            poptold = popt
        if vb:
            print("---> poptin_cc i:",i,"popt*mult:",popt*mult,"popt:",popt)
        sigmaa=lorenzbesserwiki_lr(x_all, popt[0], popt[1], popt[2])
        #sigmaa=lorenzbesserwiki(x, popt[0], popt[1], popt[2])**2
        #sigmaa=lorenzbesserwiki(x, popt[0], popt[1], popt[2])**4
        #print "bounds:",bounds
        popt, pcov = curve_fit(lorenzbesserwiki_lr, x_all, y_all,sigma=1./sigmaa,p0=popt,bounds=bounds) #,maxfev=maxfev) #,p0=[813,10**9,1])
        popt[2]=abs(popt[2])
        if abs(1. - popt[0]/poptold[0]) <= conv_crit and abs(1. - popt[2]/poptold[2]) <= conv_crit:
            if vb:
                print("done popt 33:",popt*mult,popt)
            break
        poptold = popt
    if vb:
        print("---> (2) poptout->:",popt,popt*mult,abs(1. - popt[0]/poptold[0]),abs(1. - popt[2]/poptold[2]))
        #print "freq:",freq,"freq_new:",popt[0]*mult
        #print "lifetime:",lifetime,"lifetimenew:",popt[2]*mult

    #yfit_lr= lorenzbesserwiki_lr(x_all, popt[0], popt[1], popt[2])
    #np.savetxt("tmp_lr_fit.dat",np.transpose((x_all,yfit_lr)))

    #print "popt;",popt
    #np.savetxt("tmp_in.dat",power_spectrum_full_left_right_peak)
    #yfit_lrn= lorenzbesserwiki_lr(x_all, popt[0], popt[1], popt[2])
    #np.savetxt("tmp_lrn_fit.dat",np.transpose((x_all,yfit_lrn)))

    #print "final no bound (1):",popt
    #print "final no bound (1):",popt*mult

    freq = np.abs(popt[0])*mult
    lifetime = popt[2]*mult
    popt_two_peak = popt
    if vbb:
        print("+++++++++two_peak+++++++ freq / lifetime",round(freq,od),round(lifetime,od))

    #######################################################################################################
    # 6. if band crossing (currently not, does not look good in camparision to one_peak fit which is good)
    #######################################################################################################
    #if check_qpoint_for_crossing:
    check_qpoint_for_crossing_out = False
    if check_qpoint_for_crossing_out:
        #print "lenx ",len(x_all)
        #print "lenyfit",len(yfit_lrn)
        #print "lenorig",len(power_spectrum_full_left_right_peak)
        #rest = power_spectrum_full_left_right_peak-yfit_lrn
        #np.savetxt("tmp_rest.dat",np.transpose((x_all,rest)))

        def lorenzbesserwiki_lr_4_peaks(w,w0,a,b,w02,a2,b2):
            ''' w: frequenzabhaengigkeit
                w0: die frequenz
                a hoehe  (also called the plasma frequency)
                b breite (gamma)
                took this'''
            lenw = len(w)-1.
            return (a**2)/((w**2-w0**2)**2.+w**2*b**2.) + (a**2)/(((-w+lenw)**2-w0**2)**2.+(-w+lenw)**2*b**2.)+(a2**2)/((w**2-w02**2)**2.+w**2*b2**2.) + (a2**2)/(((-w+lenw)**2-w02**2)**2.+(-w+lenw)**2*b2**2.)

        #bounds=((p0_freq - fd2*p0_freq, 0,p0_life-ld2*p0_life),(p0_freq + fd2*p0_freq,np.inf,p0_life_max_crossing))
        print("popt;",popt)
        bounds =(( 0.99*popt[0],0.99*popt[1],0.99*popt[2],0.1*popt[0],0.0   ,0.1*popt[2]),\
                (  1.01*popt[0],1.01*popt[1],1.01*popt[2],10.*popt[0],np.inf,10.*popt[2]))
        print("bounds:",bounds)
        #popt, pcov = curve_fit(lorenzbesserwiki_lr_4_peaks, x_all, y_all,p0=[popt[0],popt[1],popt[2],13212,44344156,4575]) #,bounds=bounds) #,maxfev=maxfev) #,p0=[813,10**9,1])
        ####################
        # first shot
        ####################
        popt, pcov = curve_fit(lorenzbesserwiki_lr_4_peaks, x_all, y_all, p0=[popt[0],popt[1],popt[2],popt[0],popt[1],popt[2]],bounds=bounds)


        #yfit_4= lorenzbesserwiki_lr_4_peaks(x_all, popt[0], popt[1], popt[2],popt[3],popt[4],popt[5])
        #print "y_fit_4:",yfit_4

        ####################
        # now accurate
        ####################
        bounds =(( 0.95*popt[0],0.90*popt[1],0.90*popt[2],0.5*popt[3],0.1*popt[4],0.5*popt[5]),\
                (  1.10*popt[0],1.10*popt[1],1.10*popt[2],2.0*popt[3],10.*popt[4],2.0*popt[5]))
        for i in np.arange(10):
            if i == 0:
                poptold = popt
            if vb:
                print("---> poptin_cc i:",i,"popt*mult:",popt*mult,"popt:",popt)
            sigmaa=lorenzbesserwiki_lr_4_peaks(x_all, *popt)
            #sigmaa=lorenzbesserwiki_lr_4_peaks(x_all, *popt)**2
            #print "bounds:",bounds
            popt, pcov = curve_fit(lorenzbesserwiki_lr_4_peaks, x_all, y_all,sigma=1./sigmaa,p0=popt,bounds=bounds) #,maxfev=maxfev) #,p0=[813,10**9,1])
            popt[2]=abs(popt[2])
            #perr = np.sqrt(np.diag(pcov))
            #print "perr:",perr
            #if vb:
            #    print "---> poptout:",i,popt*mult,abs(1. - popt[0]/poptold[0]),abs(1. - popt[2]/poptold[2])
            #if abs(1. - popt[0]/poptold[0]) <= 10**-4:
            #    print "1<2?",1,np.abs(1. - popt[0]/poptold[0])
            #    print "1<2?",2,10**-4
            #    sys.exit()
            if abs(1. - popt[0]/poptold[0]) <= conv_crit and abs(1. - popt[2]/poptold[2]) <= conv_crit:
                if vb:
                    print("done popt 33:",popt*mult,popt)
                break
            poptold = popt
        if vb:
            print("---> (3) poptout->:",popt,popt*mult,abs(1. - popt[0]/poptold[0]),abs(1. - popt[2]/poptold[2]))
            #print "freq:",freq,"freq_new:",popt[0]*mult
        freq = np.abs(popt[0])*mult
        lifetime = popt[2]*mult
        freq2 = np.abs(popt[3])*mult
        lifetime2 = popt[5]*mult
        if vbb:
            print("++++++++four_peak++++++++ freq / lifetime",round(freq,od),round(lifetime,od),round(freq2,od),round(lifetime2,od))



    ####################
    # write stuff
    ####################
    if True:
        #np.savetxt("tmp_in.dat",power_spectrum_full_left_right_peak)
        yfit_all=lorenzbesserwiki(x_all, *popt_one_peak)
        #np.savetxt("tmp_one.dat",np.transpose((x_all,yfit_all)))
        yfit_all=lorenzbesserwiki_lr(x_all, *popt_two_peak)
        #np.savetxt("tmp_two.dat",np.transpose((x_all,yfit_all)))
        if check_qpoint_for_crossing_out:
            yfit_4= lorenzbesserwiki_lr_4_peaks(x_all, *popt)
            #np.savetxt("tmp_4.dat",np.transpose((x_all,yfit_4)))

    #sys.exit()
    return round(freq,od),round(lifetime,od),yfit



def find_nearest(array,value,returnindex=False):
    ''' find value in array
    '''
    idx = (np.abs(array-value)).argmin()
    if returnindex:
        return idx
    return array[idx]

def get_xind_min_between_two_xind(ps,xind1,xind2):
    min = ps[xind1:xind2].min()
    #print "min;",min
    out = np.where(ps[xind1:xind2] == min)[0][0]
    #print "out:",out,out+xind1
    return out+xind1

def get_ind_of_x_peaks(ps,peaks,x_to_THz,sep_in_THz_min=False,verbose=False,qpoint=False,space_fft_which=False):
    ''' finds the x indizes of the correspondingg peaks
    ps: the 1d-array powerspectrum
    for al a sep_in_THz_min of 1 works for all data
    a sep_in_THz_min of 2 or 3 would increase the speed for some qpoints but not work for others
    '''
    verbose = False
    if verbose:
        print('--> qp',qpoint,type(qpoint))
        print('--> peaks',peaks)
        print('--> sep_in_THz_min',sep_in_THz_min)
    if type(qpoint) != bool and space_fft_which == 'space_fft_phase_':
        path = qpoint_get_l_0_0_or_l_l_0_or_l_l_l_or_t_0_0_or_t_t_t_or_t1_t1_t1_or_t2_t2_t2_from_qpoint(qpoint,args)
        if verbose:
            print('--> path',path,space_fft_which)
        if path == 't2_t2_0':
            f1,f2,flong = estimate_freqs_from_at_L_L_0_path(qpoint=qpoint,args=args)
            if verbose:
                print('--> f1,f2',f1,f2,abs(f2-f1),abs(f2-f1)*0.7)
            if peaks == 2:
                #print 'f2',f2,f1,abs(f2-f1)
                if sep_in_THz_min < abs(f2-f1)*0.7:
                    sep_in_THz_min = abs(f2-f1)*0.7
                #print '-->',np.array([int(f2/x_to_THz),int(f1/x_to_THz)])
                #sys.exit()
                #return np.array([int(f2/x_to_THz),int(f1/x_to_THz)])
    if verbose:
        print('--> sep_in_THz_min',sep_in_THz_min)
    verbose = False
    peaks= int(peaks)
    #print('len(ps)',len(ps))
    #print('len(ps)/2',len(ps)/2)
    #print('int(len(ps)/2)',int(len(ps)/2))
    #np.savetxt("kkk.dat",ps)
    y = ps[:int(len(ps)/2)]
    if verbose:
        print('get_ind_of_x_peaks:',"peaks:",peaks,"len(ps):",len(ps))
        print('get_ind_of_x_peaks:','check symmetry')
        print('get_ind_of_x_peaks:','ps[0],ps[-1]',ps[0],ps[-1])
        print('get_ind_of_x_peaks:','ps[10],ps[-1-10]',ps[10],ps[-1-10])
        if len(ps) > 1000:
            print('get_ind_of_x_peaks:',ps[1000],ps[-1-1000])
        if len(ps) > 10000:
            print('get_ind_of_x_peaks:',ps[10000],ps[-1-10000])
    #print "x_to_THz:",x_to_THz
    #print "sep_in_THz_min:",sep_in_THz_min
    if peaks== 1:
        #return np.where(ps == ps.max())[0]  # this has tow maxima, seft peak and right peak
        out = np.where(y == y.max())[0]
        if verbose:
            print(out,out.shape,type(out))
        if len(out) == 1:
            return np.array([out[0]])  # return an array
        else:
            print("problem 33, out has not only one maximum",out)
            return False
    if peaks>= 1:
        #print "sep_in_THz_min;",sep_in_THz_min
        #np.savetxt('yy',y)
        threslist=np.array([0.0001,0.001,0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9,1.0,2.0,10.0,100.0])
        thresmin=threslist.min()
        thresmax=threslist.max()
        for thres in threslist:
            out = indexes(y, thres=thres, min_dist=sep_in_THz_min/x_to_THz) # 0.001 finds more peaks than 0.1
            #print '----> thres a',thres,'out',out, 'len(out)',len(out),'peaks',peaks
            #print '----> thres a',thres,'out',np.round(out*x_to_THz,2), 'len(out)',len(out),'peaks',peaks
            if len(out) > peaks:
                thresmin = thres
            if len(out) < peaks:
                thresmax = thres

            if verbose:
                print("thres b:",thres,out*x_to_THz,out.shape,type(out))
            if len(out) == peaks:
                return np.array(out)
            if peaks == 2 and len(out[out*x_to_THz>10]) == 1: # this is the optical peak
                opticalpeak = out[out*x_to_THz>10][0]
                maxpeak = np.where(ps==ps.max())[0][0]
                #print 'out',out,"-->",maxpeak,'opt',opticalpeak
                if opticalpeak != maxpeak:
                    return np.array([maxpeak,opticalpeak])

            if len(out) < peaks:
                if verbose:
                    print("thresmin: c ",thresmin)
                    print("thresmax: c",thresmax)
                break
        #sys.exit()
        ###############################
        # here it did not find the minimum
        ###############################
        #print 'thresmin d',thresmin,type(thresmin)
        #print 'thresmax d',thresmax,type(thresmax)
        if thresmin == thresmax:
            return False

        threslist = np.arange(thresmin,thresmax,thresmin/10.)
        #print 'threslist',threslist
        thresmin=threslist.min()
        thresmax=threslist.max()
        for thres in threslist:
            out = indexes(y, thres=thres, min_dist=sep_in_THz_min/x_to_THz) # 0.001 finds more peaks than 0.1
            if len(out) > peaks:
                thresmin = thres
            if len(out) < peaks:
                thresmax = thres

            if verbose:
                print("thres:",thres,out*x_to_THz,out.shape,type(out)) #,len(out[out*x_to_THz>10])
            if len(out) == peaks:
                return np.array(out)
            if len(out) < peaks:
                #print "thresmin:",thresmin
                #print "thresmax:",thresmax
                break
    return False

def get_max_min_max_coords_idx_of_full_ps_when_bandcrossing(xy,x,ps,qpointstr,checkwriteps,verbose=False):
    ''' expecially for the 900K runs we also want to find the minimum of the smoothed curve
        a) find min of smoothed curve
        b) make averages of smoothed curve and find the true minimum around the minimum (target would be to have say 100 values in the range of +-1% around minimum
    '''
    ps=ps/ps.max()*.9; # full with peaks left and right
    #verbose=True
    if verbose:
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    psout=ps
    #print "qpointsss:",qpointstr.split("_")
    #print "qpointsss:",qpointstr.split("_")[0]
    #print "qpointsss:",qpointstr.split("_")[2]
    q1 = float(qpointstr.split("_")[0])
    q3 = float(qpointstr.split("_")[2])
    q1_q3=q1/(q3/2.)
    #if q1/(q3/2.) >
    #np.savetxt("ps",ps)
    #print "33:",ps.shape,x.shape,xy.shape
    ps=ps[:len(ps)/2]
    x=x[:len(x)/2]
    #print "x:",x.min(),x.max()
    xy=xy[:len(xy)/2]
    #print "44:",ps.shape,x.shape,xy.shape
    #np.savetxt("aa2.dat",ps)
    #np.savetxt("aa3.dat",pss)
    #np.savetxt("aa4.dat",xy)
    def get_max_min_max_of_function(ps,verbose=False,above_y_set=1,q1_q3=False):
        #np.savetxt("ps",ps)
        #verbose = True
        xmax = np.where(ps==ps.max())[0][0]
        if verbose:
            print()
            print("GET_MAX_MIN_MAX_OF_FUNCTION:")
            print("INNNN xmax:",xmax)
        ymax = ps[xmax]
        if verbose:
            print("xmax,ymax:",xmax,ymax)
        # now go to left and right of smoothfunc
        yminleft  = ps[:xmax].min()
        yminright = ps[xmax:].min()
        yminright = ps[-1]
        if verbose:
            print("yminleft :",yminleft)
            print("yminright:",yminright)
        xminleft  = np.where(ps[:xmax]==yminleft )[0][0]
        #xminrigth = np.where(ps[xmax:]==yminright)[0][0]+xmax
        xminrigth = len(ps)

        if verbose:
            print("xminleft :",xminleft)
            print("xminrigth:",xminrigth)
            print("now make function which takes all values above")
        above = np.array([0.8,.7,.6,.5,.4,.3,0.2,0.1])

        if verbose:
            print("above_y_set:",above_y_set)
        if above_y_set == 1:
            above = (np.arange(9)/10.)[np.arange(9)/10.>np.array([yminleft,yminright]).min()][::-1]
        else:
            above = (np.arange(18)/20.)[np.arange(18)/20.>np.array([yminleft,yminright]).min()][::-1]

        #above = above[np.where(above >= 0.15)] we need  teiler in about 0.1 for e.g. /cmmc/ptmp/aglen/Understand_phonon_lifetimes/check_3_sc_convergence_300K/nvt_tutild_pot_cutoff_6/sc_7/ps5_5_14
        def get_teiler_daten(ps,mindestens):
            lenps = len(ps)
            teiler = 1
            for i in np.arange(16)+1:
                tryteiler = 1*10**i
                #print i,teiler,lenps/teiler,(lenps/teiler)>=5000
                if ((lenps/tryteiler)>=mindestens) == True:
                    #print i,tryteiler,"Ja"
                    teiler = tryteiler
                else:
                    #print i,tryteiler,"NEIN"
                    break
            #print "teiler:",teiler
            N = teiler
            return N

        N = get_teiler_daten(ps,mindestens=500)
        if verbose:
            print("N==teiler:",N,"len(ps) aendert sich nicht durch das smoothing:",len(ps))

        #print "NNNin:",N,len(ps)
        if N > 1:
            ps = mysmooth(ps,N+1)  # mysmooth(y,N*3+1)  ---> dieser ist also immer ungerade!
        #print "now smoothed len(ps):",len(ps)

        if verbose:
            print("N:",N,len(ps))
        #np.savetxt("ps2",ps2)
        if verbose:
            print("above:",above)
            print(len(above)/2)
            #la = len(above)/2
            #von=int(la*1./3.)
            #bis=int(la+la*1./3.)
            #print above[von:bis]
            #print above[:von][::-1]
            #print "get middle 30 percent:",above[2la:2./3.,la+la

        for above_y in above:
            if verbose:
                print()
                print("------------------------")
                print("above_y:",above_y,"<<<<<<<<<<<<<<--------")
                print("------------------------")
            above_y_out = above_y
            allv = np.where(ps[xminleft:xminrigth]>above_y)[0]
            if len(allv) <=1:
                continue
            diff = np.diff(allv)
            diffmax = diff.max()
            if verbose:
                print("len allv:",len(allv), "diffmax>1?:",diffmax) #,"allv:",allv,diff,"len diff:",len(diff),"--------> diffmax:",diffmax
            if diffmax > 1:
                #print above_y,"allv:",allv,len(allv)
                firsthoeckerleft  = xminleft + allv[0]
                secondhoeckerright= xminleft + allv[-1]
                #print "firsthoeckerleft  :",firsthoeckerleft
                #print "secondhoeckerright:",secondhoeckerright
                xmin = np.where(ps[firsthoeckerleft:secondhoeckerright]==ps[firsthoeckerleft:secondhoeckerright].min())[0][0]
                xmin = xmin+firsthoeckerleft
                #if verbose:
                #    print above_y,"----------------------------> xmin:",xmin
                #break
                #print "xmin:",xmin
                #print "diffmax:",diffmax
                #print "allv[0]:",allv[0]
                #print "diffmax idx:",np.where(diff==diffmax)[0][0]
                #print "firsthoeckerleft  :",firsthoeckerleft
                #print "secondhoeckerright:",secondhoeckerright

                xmaxleft = np.where(ps[:xmin]==ps[:xmin].max())[0][0]
                y_xmaxleft = ps[xmaxleft]
                #print "CROSSING --> xmaxleft:",xmaxleft,y_xmaxleft
                y_xmin = ps[xmin]
                #print "CROSSING --> xmin:",xmin,y_xmin
                xmaxright = np.where(ps[xmin:]==ps[xmin:].max())[0][0]+xmin
                y_xmaxright = ps[xmaxright]
                #print "CROSSING --> xmaxright:",xmaxright,y_xmaxright
                if verbose:
                    print("ka1:",(xmin-xmaxleft))
                    print("ka2:",xminrigth)
                    print("xmin:",xmin)
                    print("xmaxleft:",xmaxleft)
                    print("xminrigth:",xminrigth)
                    print("(xmin-xmaxleft)/float(xminrigth)")
                aaa = (xmin-xmaxleft)/float(xminrigth)
                bbb = (xmaxright-xmaxleft)/float(xminrigth)
                ccc = (xmin-xmaxleft)/float(len(ps))
                ddd = (xmaxright-xmin)/float(len(ps))
                ccc = np.array([ccc,ddd]).min()
                if verbose:
                    print(above_y,"----------------------------> leftpeak,min,rightpeak (X):",xmaxleft ,xmin,xmaxright,"|a>0.002?:",aaa,"|b",bbb,"c>0.002?:",ccc)
                    print(above_y,"----------------------------> leftpeak,min,rightpeak (Y):",y_xmaxleft,y_xmin,y_xmaxright,"|")
                    print()
                #if aaa > 0.002 and ccc > 0.002:  # a sollte gleich c sein
                if q1_q3 > 0.8:
                    eee = 1. - float(xmin)/float(xmaxright)
                else:
                    eee = 1.
                #print "q1_q3:",q1_q3
                #print "aaa:",aaa,q1_q3
                #print "eee:",eee


                if aaa > 0.002 and ccc > 0.001 and eee > 0.1:  # a sollte gleich c sein
                    # a gut ########################################################
                    # a ist 0.0299818009697 bei 8_8_20 (t2_t2_0)    # 300K
                    # a ist 0.0876429458276 bei 9_9_20              # 300K
                    # a ist 0.0044179368235 bei 4_4_8 (t2)          # 900K (DFT GGA)
                    # a ist 0.00444444444444 bei 2_2_4 (t2)         # 300K (DFT CrN)
                    #
                    # a ungut #######################################################
                    # a is 0.0068 (== gut)
                    # ein zu kleines a ist 0.000721783607492
                    # ein zu kleines a ist 0.000792169219702
                    # ein zu kleines a ist 0.011 bei 9_9_20 900K /check_4_temperatuer_effects/900K_LDA_4.08

                    # c gut ########################################################
                    # c gut @ 900K 4_4_8    0.0031
                    # c gut @ 300K 10_10_20 0.114 (smothed 0.105)
                    # c gut @ 300K  8_8_20  0.049 (smoothed 0.038)
                    # c gut @ 300K 2_2_4    0.0014 (300K CrN)
                    #                                                   0.00066 kann aber auch gut sein....(bei wenigen daten)
                    # c ungut ########################################################
                    # c ungut @ 900K  10_10_20  0.055(=gut)  (smoothied 0.00069=ungut)
                    #break
                    return xmaxleft ,xmin,xmaxright,xminleft,xminrigth,yminleft,len(ps),y_xmaxleft,y_xmin,y_xmaxright,above_y

                #if aaa < 0.02 or ccc < 0.01:
                #    # checke erst ob zwischen max1 und max2 evtl. ein anderer wert hoecher ist als
                #HIER HIER HIER
                #@@@if above_y == 0.1:
                #@@@    # checke erst ob zwischen max1 und max2 evtl. ein anderer wert hoecher ist als
                #@@@    ps[:xmaxleft

        #print undefinedvariable_raises_unboundlocalerror
        return
        #return xmaxleft ,xmin,xmaxright,xminleft,xminrigth,yminleft,len(ps),y_xmaxleft,y_xmin,y_xmaxright,above_y

    #########################################
    # eigentliches powerspektrum
    #########################################
    #xmaxleft1,xmin1,xmaxright1,xminleft1,xminrigth1,yminleft1,lenps1 = get_max_min_max_of_function(ps,verbose=False)
    #print "&&&&&&&&&&&&&&&&&&&&&&&&&&&&> leftpeak,min,rightpeak (X):",xmaxleft1,xmin1,xmaxright1,"|a",(xmin1-xmaxleft1)/float(xminrigth1),"|b",(xmaxright1-xmaxleft1)/float(xminrigth1)

    #########################################
    # gesmoothtes powerspektrum
    # lets use this minimum in general since it is better
    #########################################

    #print "x:",x.min(),x.max()
    #print "------------------"
    #np.savetxt("pss",ps)
    #print "------------------"
    try:
        verboseshow=False
        xmaxleft,xmin,xmaxright,xminleft,xminrigth,yminleft,lenps,y_xmaxleft,y_xmin,y_xmaxright,above_y = get_max_min_max_of_function(ps,verbose=verboseshow,q1_q3=q1_q3)
        if verboseshow:
            f=1000.
            print("len(ps):",len(ps)*2,'xxx1:',xmaxleft*f/len(ps),xmin*f/len(ps),xmaxright*f/len(ps))
    except (UnboundLocalError, TypeError):
        #print "HHHHHHHHHHHHHHHHHHHHHHHHH"
        try:
            xmaxleft,xmin,xmaxright,xminleft,xminrigth,yminleft,lenps,y_xmaxleft,y_xmin,y_xmaxright,above_y = get_max_min_max_of_function(ps,verbose=verboseshow,above_y_set=2,q1_q3=q1_q3)
            if verboseshow:
                f=1000.
                print("--> len(ps):",len(ps)*2,'xxx1:',xmaxleft*f/len(ps),xmin*f/len(ps),xmaxright*f/len(ps))
        except (UnboundLocalError, TypeError):
            #print "YYYYYYYYYYYYYYYYYYYYY"
            #print_warning("IT SEMMS THERE IS ONLY ONE MAXIMUM for ps lenghts of "+str(len(psout)))
            print("IT SEMMS THERE IS ONLY ONE MAXIMUM for ps lenghts of "+str(len(psout)))
            #np.savetxt("psout",psout)
            return psout/psout.max()*.9

    #if verbose:
    xmaxleft2 = x[xmaxleft]
    xmin2 = x[xmin]
    xmaxright2 = x[xmaxright]
    #if verbose:
    if False:
        #np.savetxt("ps448_"+str(len(ps)),ps)
        print("left, min, right:",xmaxleft2,xmin2,xmaxright2,"--> above_y:",above_y,"--> len",len(ps),"--> totallen:",len(ps)*2.)
        print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&> leftpeak,min,rightpeak (X):",xmaxleft2,xmin2,xmaxright2,"|a",(xmin-xmaxleft)/float(xminrigth),"|b",(xmaxright-xmaxleft)/float(xminrigth))
        print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&> leftpeak,min,rightpeak (Y):",y_xmaxleft,y_xmin,y_xmaxright,"<<<- diese Y sind egal da nur die X'es gebraucht werden und das ps zurechtzuschneiden")
        np.savetxt("pstest",ps)

    #print "x:",x.min(),x.max()
    #np.savetxt("ps2",ps2)
    #np.savetxt("ps3",mysmooth(ps2,101))
    #np.savetxt("ps",ps)
    if xmaxleft > 2:  # > 1 macht fehler weil dann kein a.min
        #print "xml:",xmaxleft,xmaxleft
        a = ps[:xmaxleft][2:]  # first 2 values are 0!!
        #np.savetxt("a",a)
        #print "a:",a
        #print "a.min():",a.min()

        xminfueruntergrund = np.where(a==a.min())[0][0]
        y = a[xminfueruntergrund]
        if False:
            #print "a:",a
            print("xminfueruntergrund:",xminfueruntergrund,"y:",y)


        # change ps to psout!!! to keep the length!
        psout[0:xmin] = y  # linke seite des powerspektrums setze auf y
        psout[-xmin:] = y  # rechte seite des powerspektrums setze auf y
        psout[0:xmin] = 0  # linke seite des powerspektrums setze auf 0 besser fuer lorenzfit
        psout[-xmin:] = 0  # rechte seite des powerspektrums setze auf 0 besser fuer lorenzfit
    #print "22;",x.shape,ps.shape
    #def xyplot(x,y):
    #    return np.transpose(np.array([x,y]))
    if verbose:
        print("psmax :",psout.max())
    #print "psmax2:",ps[:len(ps)/2.max()
    #ps=ps/ps[:len(ps)/2].max()*.9; # full with peaks left and right
    psout=psout/psout.max()*.9; # full with peaks left and right

    if checkwriteps:
        if verbose:
            print("che",checkwriteps)
        #np.savetxt("ps"+qpointstr+"_"+str(checkwriteps)+".dat.cut",xyplot(x,ps))
    #pssave=pssave/pssave[:len(pssave)/2].max()*.9; # full with peaks left and right
    #ps=ps/ps[:len(ps)/2].max()*.9; # full with peaks left and right
    # funktion spiegeln
    #out = np.concatenate([psout,ps[::-1]])
    return psout

def get_goodsmoothing_for_ps_with_crossing(ps):
    min_smoothing = 1./len(ps)
    order_magnitude = -1.*(int(math.log10(min_smoothing))-1)
    #print order_magnitude
    order_magnitude = np.arange(1,order_magnitude+1)[::-1]
    #print "ord:",order_magnitude
    for ido,i in enumerate(order_magnitude):
        #print "------>",ido,i,10**(-i)
        order_magnitude[ido] = 10**(-i)
    loop_smoothing = order_magnitude
    #print "2:",loop_smoothing

    #print "STARTING loop 0 -> (logarithmic loop) 1e-07 .. 1e-01" #,loop_smoothing
    working = []
    minimum = []
    x_max1list = []
    x_max2list = []
    minimumalt = []
    for smoothingfact in loop_smoothing[::-1]:   # [::-1] heisst einfach nur die loopsmoothing andersrum lesen, alle werte bleiben erhaleten
    #for smoothingfact in loop_smoothing:
        print("smoothingfact:",smoothingfact)
        if len(working) >= 2: continue   # to orders of magnitude are enough, the third one often fails
        out = smoothing_power_spectrum(ps,smoothingfact)  # ps and out1 are the full power spectrum (left and right peak)
        #np.savetxt("smooth"+str(smoothingfact)+".dat",out)
        if out[0] > 0.5 : continue
        working.append(True)

        # debug line crossing
        #np.diff(out)
        #outsave=np.diff(out)/np.max(np.diff(out))*.9;
        #np.savetxt("pssmoothed"+str(smoothingfact)+"diff.dat",out)

        outhalfright = out[:len(out)/2]
        outhalfleft = np.copy(outhalfright)
        #np.savetxt("newpssmoothed"+str(smoothingfact)+".dat",out1)
        x_max1 = np.where(outhalfright == outhalfright.max())[0][0]    # das ist aber fuer den [ 8 8 20 ] schon der zweite maximalwert
        #print "X_MAXXXXX1:",x_max1


        outhalfright[0:x_max1] = 0
        outhalfleft[x_max1:] = 0

        #np.savetxt("newpssmoothed"+str(smoothingfact)+".dat",outhalfright)
        out1diffright = np.diff(outhalfright)
        out1diffright[x_max1-1] = 0
        #np.savetxt("diff"+str(smoothingfact)+".dat",out1diffright)
        out1diffmaxright = np.where(out1diffright == out1diffright.max())[0][0]

        #np.savetxt("pssmoothed"+str(smoothingfact)+".dat",outhalfleft)
        out1diffleft = np.diff(outhalfleft)
        out1diffleft[x_max1-1] = 0
        #np.savetxt("pssmootheddiff"+str(smoothingfact)+".dat",out1diffleft)
        out1diffmaxleft = np.where(out1diffleft== out1diffleft.max())[0][0]

        out1diffleft[0:out1diffmaxleft] = 0
        out1diffmaxleftright = np.where(out1diffleft == out1diffleft.min())[0][0]
        out1diffleft[0:out1diffmaxleftright] = 0
        out1diffleft = out1diffleft*10000
        lmmm = np.where(out1diffleft == out1diffleft.min())[0][0]
        out1diffleft[0:lmmm] = -1
        lrrr = np.where(out1diffleft == out1diffleft.max())[0][0]
        out1diffleft[lrrr:] = 1

        lminimum = find_nearest(out1diffleft, 0.,returnindex=True)
        minimumalt.append(lminimum)
        #print "lll;",lminimum," this is important for [ 8 8 20 ]

        #out1diffmaxleftrightright = np.where(out1diffleft== out1diffleft.max())[0][0]
        #out1diffleft[0:out1diffmaxleftrightright] = 1

        #np.savetxt("pssmootheddiff"+str(smoothingfact)+".dat",out1diffleft)


        #print "outdiffmaxright:",out1diffmaxright,"out1diffmaxleft:",out1diffmaxleft

        tmp = np.copy(outhalfright)
        tmp[0:out1diffmaxright] = 0
        x_max2right = np.where(tmp == tmp.max())[0][0]

        tmp = np.copy(outhalfright)
        tmp[0:x_max1] = 1
        #np.savetxt("outoutka"+str(smoothingfact)+".dat",tmp)
        tmp[x_max2right:] = 1
        #np.savetxt("outoutkb"+str(smoothingfact)+".dat",tmp)
        x_minimumright = np.where(tmp == tmp.min())[0][0]   # [ 10 10 20 ]: 339748, 345544, nehme dass kleinere

        print(smoothingfact,"x_min:", x_minimumright,"l_min:", lminimum)
        if x_minimumright < x_max1:
            x_minimumright = len(ps)
        print(smoothingfact,"x_max1:",x_max1,out1diffmaxright,"x_max2right:",x_max2right,round(float(x_max2right)/float(x_max1),2),"x_minimumright:",x_minimumright)

        minimum.append(x_minimumright)
        #print "fuer den transversalen brauchen wir den hoeheren (rechten) gaussian."
        #print "ps muss entweder auf 0 gesetzt werden beim minimum, oder besser, man kopiert den hintergrund bis zum minimum."
        #print " bis zu 1/3 von x_max1 ist alles nur hintergrund, dass kann kopiert werden"
        #template_len = len(template)
        x_max1list.append(x_max1)
        x_max2list.append(x_max2right)

    print("minialt:",minimumalt)
    minimum = np.array(minimum).min()
    minimumalt = np.array(minimumalt).min()
    x_max1 = np.array(x_max1list).min()
    x_max2 = np.array(x_max2list).min()
    if x_max1 == x_max2: minimum = minimumalt
    print()
    print("MINIMUM:",minimum)
    print("X_MAX1 :",x_max1)
    print("X_MAX2 :",x_max2)
    print()

    #minimum = 339748
    #x_max1 = 258439
    template = np.copy(ps[:x_max1/3+2])
    psout = np.copy(ps)
    #np.savetxt("psin",ps)
    print("x_max1:",x_max1)
    print("minimum:",minimum)
    for x in np.arange(100):
        ab  = x * x_max1/3
        bis = x * x_max1/3 + x_max1/3
        ab_  = len(psout)-ab    # spiegelverkehrt
        bis_ = len(psout)-bis   # spiegelverkehrt
        if bis > minimum:
            bis = minimum
            bis_ = len(psout)-minimum
        print("ia:",x,"ab:",ab,"bis:",bis,"checks:",psout[ab:bis].shape,template.shape,"lenK:|",len(psout[ab:bis]))
        psout[ab:bis+1] = template[:len(psout[ab:bis+1])]   # das minus x ist wichtig weil: np.arange(10)[0:4] == array([0, 1, 2, 3]) und np.arange(10)[5:10] ==  array([5, 6, 7, 8, 9])
        print("ib:",x,"ab:",bis_,"bis:",ab_,"lenK:",len(psout[bis_-1:ab_-1]),len(template[:len(psout[ab:bis])][::-1]),len(template[:len(psout[bis_:ab_])][::-1]))
        psout[bis_-1:ab_] = template[:len(psout[ab:bis+1])][::-1]
        print("done")
        if bis >= minimum:
            break
    #np.savetxt("psout",psout)
    return psout

def get_goodsmoothing_for_ps(ps,dt,mdsteps,allowed_func = 0, allowed_der = False,args=False,stringadd=''):
    ''' if allowed_der = False: this will not be checked
        if allowed_der = 0 or another integer it is made sure that this is fulfilled
        an dieser stelle ist ps schon gesmoothed um fakrot 10 oder 100 oder mehr
    '''
    #faktor = 1./ps.max()*0.9;
    #ps=ps*faktor; # full with peaks left and right
    #print "args.sm:",args.smoothing
    #print "11111111:",ps.max()
    if type(args.smoothing) != bool:  # this is only run when explicitly specified in args. inputparameters, default is FALSE
        out1 = smoothing_power_spectrum(ps,args.smoothing)  # ps is the full power spectrum (left and right peak)
        #np.savetxt("out1",out1)
        #sys.exit()
        f,lt,ltmin,ltmax,xv,rl,rrll, xcutoff= get_uncertainty_smoothed_powerspectrum(out1,dt,args.smoothing,args=args,stringadd=stringadd)
        od = 4
        return args.smoothing, round(f,od),round(lt,od),round(ltmin,od),round(ltmax,od),xcutoff,out1   # out1 is the smoothed function

    #if args.verbose:
    #    print printblue("###############################  "+stringadd+"   ######### mdsteps: "+str(mdsteps))
    if type(allowed_func) != int:
        sys.exit("ERROR:allowed_func has to be an integer")
    if type(allowed_der) != int:
        allowed_der = 9**40
        #print allowed_der

    #######################################################
    # determine smoothingfactors
    #######################################################
    # determine smoothingfactors
    min_smoothing = 1./len(ps)
    #print "lenps:",len(ps)
    order_magnitude = -1.*(int(math.log10(min_smoothing))-1)
    #print order_magnitude
    order_magnitude = np.arange(1,order_magnitude+1)[::-1]
    #print "ord:",order_magnitude
    for ido,i in enumerate(order_magnitude):
        #print "------>",ido,i,10**(-i)
        order_magnitude[ido] = 10**(-i)
    loop_smoothing = order_magnitude
    #print "loop_smoothing:",loop_smoothing

    #if args.verbose: print "STARTING loop 0 -> (logarithmic loop) 1e-07 .. 1e-01" #,loop_smoothing
    bestworking=False
    #np.savetxt("ps.dat",ps)
    new = True
    if new:
        tmparray = np.zeros((len(loop_smoothing),8))
    #print tmparray
    #print
    #np.savetxt("ps",ps)
    for idxtmp,smoothingfact in enumerate(loop_smoothing): # do them all and only afterwards take the broadest which is working, otherwise a local minimum might screw up the result
        #print idxtmp,"smooth:",smoothingfact,"from:",loop_smoothing
        out1 = smoothing_power_spectrum(ps,smoothingfact)  # ps and out1 are the full power spectrum (left and right peak)
        #np.savetxt("out1_"+str(smoothingfact),out1)
        #sys.exit()
        #print "lenout1:",len(out1)
        #if args.write_full_ps:
        #    np.savetxt("out_loop1_"+str(smoothingfact),out1)


        # debug line crossing
        #np.savetxt("pssmoothed"+str(smoothingfact)+".dat",out1)
        #np.diff(out1)
        #outsave=np.diff(out1)/np.max(np.diff(out1))*.9;
        #np.savetxt("pssmoothed"+str(smoothingfact)+"diff.dat",outsave)
        #np.savetxt("kkb.dat",out1)
        try:
            f,lt,ltmin,ltmax,xv,rl,rrll,xcutoff = get_uncertainty_smoothed_powerspectrum(out1,dt,smoothingfact,args=args,stringadd=stringadd)
        except ValueError:
            pass
        if not new:
            if rl <= allowed_func and rrll <= allowed_der:
                bestworking = smoothingfact
                break
        if new:
            tmparray[idxtmp] = f ,lt,ltmin,ltmax,xv,rl,rrll,xcutoff
            #print tmparray
            #print
    if new:
        printtmparray = False
        if printtmparray:
            print(tmparray)
        rlall = tmparray[:,5]
        rrllall = tmparray[:,6]
        bestworking = False
        for idx in np.arange(len(loop_smoothing))[::-1]:
            if printtmparray:
                print(idx,loop_smoothing[idx],rlall[idx],rrllall[idx],bestworking)
            if rlall[idx] <= allowed_func and rrllall[idx] <= allowed_der:
                bestworking = loop_smoothing[idx]
                smoothingfact = bestworking
            else:
                break

    #print "(mdsteps)",mdsteps,"bestworking:",bestworking
    ##sys.exit()
    #for idxtmp,smoothingfact in enumerate(loop_smoothing[::-1]):
    #    print "idxtmp:",idxtmp,smoothingfact

    #sys.exit()
    #print "rlall",rlall
    #print "rrllall",rrllall
    #sys.exit()
    loop_smoothing = smoothingfact/10.*np.array([1,2,3,4,5,6,7,8,9,10])  # redo 10 for easier coding, 10 will definitively work.
    loop_smoothing_next = smoothingfact/100.*np.array([1,2,3,4,5,6,7,8,9,10])  # redo 10 for easier coding, 10 will definitively work.
    if type(bestworking) == bool:
        print("no bestworking found! maybe band crossing??? check this!!!!!! in this case you would have 2 gaussians and you'd need to get the linewidths of both")
        bestworking = 1.0
    #if args.verbose: print "DONE loop 1 ->","best:",bestworking,"try:",loop_smoothing

    for smoothingfact in loop_smoothing:
        out1 = smoothing_power_spectrum(ps,smoothingfact)  # ps is the full power spectrum (left and right peak)
        f,lt,ltmin,ltmax,xv,rl,rrll, xcutoff= get_uncertainty_smoothed_powerspectrum(out1,dt,smoothingfact,args=args,stringadd=stringadd)
        if rl <= allowed_func and rrll <= allowed_der:
            bestworking = smoothingfact
            break
    if args.verbose > 2:
        print("-----------------------------------------------------------------------------------")
    #print "loop_smoothing:",loop_smoothing
    delta = loop_smoothing[0]/10.
    #print "delta:",delta
    #if args.verbose: print "best:",bestworking
    #print np.arange(1,11)[::-1]
    loop_smoothing = bestworking - delta * np.arange(0,11)[::-1]  # include 0 so bestworking will definitively work (easier code)
    loop_smoothing = np.trim_zeros(loop_smoothing)
    #order_magnitude = -1.*(int(math.log10(min_smoothing))-1)
    loop_smoothing = loop_smoothing[loop_smoothing > 0]
    #if args.verbose: print "DONE loop 2 ->","best:",bestworking,"try:",loop_smoothing
    for smoothingfact in loop_smoothing:
        out1 = smoothing_power_spectrum(ps,smoothingfact)  # ps is the full power spectrum (left and right peak)
        f,lt,ltmin,ltmax,xv,rl, rrll, xcutoff= get_uncertainty_smoothed_powerspectrum(out1,dt,smoothingfact,args=args,stringadd=stringadd)
        if rl <= allowed_func and rrll <= allowed_der:
            bestworking = smoothingfact
            break
    if args.verbose > 2:
        print("-----------------------------------------------------------------------------------")
    #if args.verbose:
    #    print "DONE loop 3 ->",smoothingfact,"best:",bestworking

    #    print "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"
    #    print "bestworking:",bestworking
    sigma = bestworking
    #print "22222222:",ps.max()
    out1 = smoothing_power_spectrum(ps,sigma,False)  # ps is the full power spectrum (left and right peak)
    fqgood, ltgood, ltmin, ltmax, xv, rl, rrll, xcutoff= get_uncertainty_smoothed_powerspectrum(out1,dt,sigma,args=args,stringadd=stringadd,verbose=True,color='red')
    #print "wriet out1",mdsteps
    od = 4
    #print "33333333:",ps.max(),out1.max()
    #print "len(ps):",len(ps)
    #print "len(out1):",len(out1)
    return sigma, round(fqgood,od),round(ltgood,od),round(ltmin,od),round(ltmax,od),xcutoff,out1   # out1 is the smoothed function

def powerspectrum_to_powerspectrum_2d_sparse(ps,smooth=True,mindestens=10000,verbose=False):
    ''' print "len of ps:",len(ps)
        print " immer so fitten dass mindestens 5000 werte drin bleiben "
        print " immer so fitten dass das ps (out) eine gerade anzahl an eintraegen hat -> das macht aber mysmooth nicht!"
        x_to_THz = 1000./(args.dt*float(mdsteps))*pslenfull/pslencurr  # sollte 0.01 nicht ueberschreiten
        x_to_THz = 1000./(args.dt*float(mdsteps))*len(ps)/pslencurr  # sollte 0.01 nicht ueberschreiten
    '''
    #verbose = True
    if len(ps) <= mindestens:
        x = np.arange(len(ps))
        y = ps
        xy=np.transpose(np.array([x,y]))
        x_to_THz = 1000./(len(x)*args.dt)
        return xy,x,ps,x_to_THz

    psmax=ps.max()  # notewendig da ansonsten der l 0 0 zweig reskaliert wird
    tmp = ps[:int(len(ps)/2)]
    anzahl_punkte_ueber_FWHM = len(np.where(tmp > tmp.max()/2.)[0])
    if verbose:
        print(" xx len(ps):",len(ps))
        print(" xx ps.max:",tmp.max())
        print(" xx anz punkte > halbps.max():",np.where(tmp > tmp.max()/2.)[0])
        print(" xx anz punkte               :",len(np.where(tmp > tmp.max()/2.)[0]))
        print(' xx anzahl_punkte_ueber_FWHM :',anzahl_punkte_ueber_FWHM)
    if anzahl_punkte_ueber_FWHM <= 10:
        x = np.arange(len(ps))
        y = ps
        xy=np.transpose(np.array([x,y]))
        x_to_THz = 1000./(len(x)*args.dt)
        return xy,x,ps,x_to_THz

    def get_teiler_daten(ps,mindestens):
        lenps = len(ps)
        teiler = 1
        for i in np.arange(16)+1:
            tryteiler = 1*10**i
            #print i,teiler,lenps/teiler,(lenps/teiler)>=5000
            if ((lenps/tryteiler)>=mindestens) == True:
                #print i,tryteiler,"Ja"
                teiler = tryteiler
            else:
                #print i,tryteiler,"NEIN"
                break
        #print "teiler:",teiler
        N = teiler
        return N

    N = get_teiler_daten(ps,mindestens)
    if verbose:
        print(" xx N = get_teiler_daten:",N)
    if N < 2:
        x = np.arange(len(ps))
        y = ps
        xy=np.transpose(np.array([x,y]))
        x_to_THz = 1000./(len(x)*args.dt)
        return xy,x,ps,x_to_THz



    #print "N:",N
    if verbose:
        print(' xx len(ps) OLD a:',len(ps))
    x = np.arange(len(ps))[::N]+N/2
    if verbose:
        print(' xx len(x)  NEW a:',len(x))
        print(' xx x.min()',x.min(),'x.max()',x.max())
        print()

    for i in np.arange(20):
        x_to_THz = 1000./(len(x)*args.dt)
        if verbose:
            print(' xx i',i,'x_to_THz:',str(x_to_THz).ljust(10),"x_to_THz > 0.02?",x_to_THz>0.02,'N',N,"N>2?",N>2)
        if x_to_THz > 0.02 and N > 2:  # we want to resolve a lifetime change of 0.02 or better (0.06 ist zu grob)  # if N == 1 dont go on
            #print "N:",N,mindestens,"len(x):",len(x),"x_to_THz:",1000./(len(x)*args.dt)
            #print "i:",i,"Nold:",N,
            N = N / 10
            #print "Nnew:",N
            x = np.arange(len(ps))[::N]+N/2
            #try:
            #    x = np.arange(len(ps))[::N]+N/2
            #except ValueError:
            #    print "N:",N
            #    printwarning("N = "+str(N))
            #    sys.exit()
        else:
            break
    x = x/N
    if verbose:
        print(' xx len(x)  NEW b:',len(x))
    #print "N fin:",N,mindestens,"len(x):",len(x),"x_to_THz:",x_to_THz
    #print "x:",x[:3],x[-3:],len(x)
    #print "ps 11 (ps.max):",ps.max(),ps.shape

    # das mean nimmt hier immer den mean aus einem bestimmten bereich;
    # dies skaliert im endeffekt auch die daten
    y = np.mean(ps.reshape(-1, N), 1)    # hier wird irgendwie skaliet!! das will ich aber nicht
    #print "ps 22 ( y.max):",y.max(),y.shape
    y = y/y.max() * psmax
    if verbose:
        print(' xx len(y)  NEW b:',len(y))
    #@ print "ps 33 ( y.max):",y.max(),y.shape
    #@ def func(y, a):
    #@         return y*a
    #@ #import scipy.optimize as optimization
    #@ #print optimization.curve_fit(func, y, ps, x0, sigma)

    #@ def lin(y,a):
    #@     return y*a
    #@ # ps has a length of 3750000
    #@ #  y has a length of 37500
    #@ print "len(ps):",len(ps)
    #@ print "len(x ):",len(x)
    #@ print "len(y ):",len(y)
    #@ print "len(ps[x]):",len(ps[x]),ps[x].max()
    #@ np.savetxt("psx",ps[x])
    #@ popt0, pcov0 = curve_fit(lin, y, ps[x])
    #@ print "popt0:",popt0
    #@ #sys.exit()
    if smooth:
        print("smooth is TRUE:",N) # N = z.b. 100
        y = mysmooth(y,N+1)  # mysmooth(y,N*3+1)  ---> dieser ist also immer ungerade!
    #y = y/y.max()*0.9
    #y = y/y.max()*psmax  # bring it to old hight
    #print "x.sh:",x.shape
    #print "y.sh:",y.shape
    xy=np.transpose(np.array([x,y]))

    #print "xy:",xy[-3:]

    #ps=ps/ps[:len(ps)/2].max()*.9; # full with peaks left and right
    #y=xy[:,1];y=y/y[:len(y)/2].max()*.9;xy[:,1]=y;
    #return np.transpose(np.array([x,y])),x,y
    #np.savetxt("ps_ka",ps)
    #sys.exit()
    #return xy,x,ps

    #x_to_THz = 1000./(args.dt*float(len(ps)))*len(ps)/len(y) # sollte 0.01 nicht ueberschreiten
    #xx_to_THz = 1000./(args.dt)/len(y) # sollte 0.01 nicht ueberschreiten
    #xxx_to_THz = 1000./(len(y)*args.dt) # sollte 0.01 nicht ueberschreiten
    #print "x_to_THZ(infunc):",x_to_THz,xx_to_THz,xxx_to_THz
    #print "len(ps):",len(ps),len(y),len(ps)/len(y)
    return xy,x,y,x_to_THz

def check_qpoints_for_crossing(qpoint, N, structure,space_fft_which="space_fft_"):
    '''
    just gives true for fcc t1 for qpoint [[ 10 10 20 ], [ 9 9 20 ], [ 8 8 20 ]]
    NOW: add True when t2 for space_fft_phase
    it seems we got rid of crossing
    '''
    return False
    #print 'space_fft_which:',space_fft_which
    if space_fft_which == 'space_fft_phase_':
        return False
        #print 'in2'
        #if qpoint[0] == qpoint[1] and qpoint[2] == 2*N:
        #    #print 'in3'
        #    return True
        #else:
        #    return False
    #print "iiiiiiiiiii",qpoint,N
    returnvalue = False
    #crossingfrom=float(8)/float(10) # there seems to be crossing in DFT GGA 300K [3 3 8]
    #crossingfrom=float(7.5)/float(10)
    crossingfrom=1./np.sqrt(2)  # 0.707 seems wrong, only one peak in /cmmc/ptmp/aglen/Understand_phonon_lifetimes/check_3_sc_convergence/nvt_tutild_pot_cutoff_6/sc_7 for qpoint [5 5 14]
    crossingfrom=0.79     # at 900K we cant distinguish [7 7 18] 0.777 in /cmmc/ptmp/aglen/Understand_phonon_lifetimes/check_3_sc_convergence/nvt_tutild_pot_cutoff_6/sc_9
    crossingfrom=0.75     # at 300K we can distinguish [ 3 3 8] ==0.75 in /cmmc/ptmp/aglen/Understand_phonon_lifetimes/check_3_sc_convergence_300K/nvt_tutild_pot_cutoff_6/sc_4
    crossingfrom=0.72     # at 300K we can distinguis [ 5 5 14] ==0.714 in /cmmc/ptmp/aglen/Understand_phonon_lifetimes/check_3_sc_convergence_300K/nvt_tutild_pot_cutoff_6/sc_7

    if qpoint[0]==qpoint[1] and qpoint[2]==2*N:  # bedingutn t1
        if structure == 'fcc': # bedingung fcc
            if float(qpoint[0])/float(N)>=crossingfrom:   # nur [ {8,9,10} {8,9,10} 20]
                returnvalue = True
    elif qpoint[0]==qpoint[1] and qpoint[2]==-2*N:  # bedingutn t1
        if structure == 'fcc': # bedingung fcc
            if float(qpoint[0])/float(N)>=crossingfrom:   # nur [ {8,9,10} {8,9,10} 20]
                returnvalue = True
    else:
        returnvalue = False
    return returnvalue

def get_correct_input_variables(args):
    #print "in structure     :",args.structure
    #print "in N (supercell) :",args.supercell
    #print "in alat          :",args.alat
    #print "in dt [fs]       :",args.dt
    #print "in usestruct     :",args.usestruct
    #print "in qvec          :",args.qvec
    #print "---------------------------------"
    if args.structure == "fcc": args.usestruct = 4;
    if args.structure == "bcc": args.usestruct = 2;
    if os.path.isfile("struct.bcc"):
        args.structure = "bcc"
        args.usestruct = 2
    #print "ar:",args.lammps_infile
    #print glob.glob("infile.*")
    ##if os.path.isfile(glob.glob("infile.*")):
    ##    args.infile = glob.glob("infile.*")[0]
    ##    print 'ka',args.infile

    #N = args.supercell
    #alat = args.alat
    #dt = args.dt

    ########################################
    # in case of LAMMPS job
    ########################################
    get_inputvariables_from_calculation(positionsfilename = "positions.*",infilefilename = "in_file_dynamics.in", inN=args.supercell, inalat=args.alat, indt=args.dt, inqvec=args.qvec, verbose=args.verbose,args=args)
    ########################################
    # in case of VASP job
    ########################################
    if args.supercell == False or args.alat == False or args.dt == False:
        get_inputvariables_from_vasp_job()
    ########################################
    # in case of md_long_tox (LA) c skript
    ########################################
    if args.supercell == args.alat == args.dt == False:
        get_inputvariables_from_c_skript_md_long_tox_job()

    #if args.structure == "fcc": args.usestruct = 4;
    #if args.structure == "bcc": args.usestruct = 2;
    if not args.notverbose:
        print("structure       :",args.structure)
        print("N (supercell)   :",args.supercell)

    if args.supercellmultiply:
        args.supercell = args.supercell*args.supercellmultiply
        if not args.notverbose:
            print("N (multiply )   :",args.supercellmultiply)
            print("N (supercell final):",args.supercell)

    if type(args.atoms) == bool:
        args.atoms = args.usestruct * (args.supercell**3)


    if not args.notverbose:
        print("alat            :",args.alat)
        print("dt [fs]         :",args.dt)
        print("usestruct       :",args.usestruct)
        print("number of atoms :",args.atoms)
        print("qvec            :",args.qvec)  # qvec            : [['l', '0', '0'], ['t', '0', '0']]

    if args.dt > 50:
        sys.exit('Error: your timestep dt seems too large...')
    if args.dt < 0.1:
        sys.exit('Error: your timestep dt seems quite small, is it corect?')
    #sys.exit('ka')
    if type(args.qvec) == bool:
        sys.exit("ERROR: specify qvec!")
    if type(args.dt) == bool:
        sys.exit("ERROR: dt! This is also necessary for writing the output correctly.")
    # now print all qvec
    #qpoints_all  :     [ qpoint      ] qpstr     symeq     Nr. symeqtot     Reduced wave vector
    # -------------------------------------------------------------------------------------------
    #     N 0 0 (L)     [ 1   0   0   ] 1_0_0      (3)     	1    (3)   	 0.0625
    #     N 0 0 (L)     [ 2   0   0   ] 2_0_0      (3)     	2    (6)   	 0.125
    #     N 0 0 (L)     [ 3   0   0   ] 3_0_0      (3)     	3    (9)   	 0.1875
    #     N 0 0 (L)     [ 4   0   0   ] 4_0_0      (3)     	4    (12)   	 0.25
    #     N 0 0 (L)     [ 5   0   0   ] 5_0_0      (3)     	5    (15)   	 0.3125
    #     N 0 0 (L)     [ 6   0   0   ] 6_0_0      (3)     	6    (18)   	 0.375
    #     N 0 0 (L)     [ 7   0   0   ] 7_0_0      (3)     	7    (21)   	 0.4375
    #     N 0 0 (L)     [ 8   0   0   ] 8_0_0      (3)     	8    (24)   	 0.5
    #     N 0 0 (L)     [ 9   0   0   ] 9_0_0      (3)     	9    (27)   	 0.5625
    #     N 0 0 (L)     [ 10  0   0   ] 10_0_0     (3)     	10   (30)   	 0.625
    #     N 0 0 (L)     [ 11  0   0   ] 11_0_0     (3)     	11   (33)   	 0.6875
    #     N 0 0 (L)     [ 12  0   0   ] 12_0_0     (3)     	12   (36)   	 0.75
    #     N 0 0 (L)     [ 13  0   0   ] 13_0_0     (3)     	13   (39)   	 0.8125
    #     N 0 0 (L)     [ 14  0   0   ] 14_0_0     (3)     	14   (42)   	 0.875
    #     N 0 0 (L)     [ 15  0   0   ] 15_0_0     (3)     	15   (45)   	 0.9375
    print('aa',args.qvec)
    qpoints_all = get_all_qpoints(args.qvec,args,verbose=True)


    if args.create_lammps_inputfile:
        print("commenot out writing of positions")
        os.popen("sed -i 's|^dump dump1 all xyz|#dump dump1 all xyz|' in_file_dynamics.in").read().rstrip()
    return qpoints_all

######################################################################
# in case of lammps job
######################################################################
def check_for_lammps_job():
    ''' if lammps inputfile is found, lammps job is evaluated '''
    print(bcolors.FAIL + "check_from_space_fft" +bcolors.ENDC)
    check_for_lammps_infile = glob.glob("in_file_dynamics.in")
    #print "check_for_lammps_infile:",check_for_lammps_infile
    if not len(check_for_lammps_infile) == 1:
        print("This does not look like a lammsp run!")
        return
    else:
        print(" -> 1 I assume this is a lammps job since I found the in_file_dynamics.in file")
        if not check_from_space_fft(args):
            check_from_xaa_folder(args)
    print(bcolors.FAIL + "check_from_space_fft DONE" +bcolors.ENDC)
    return

def check_for_all_ps_in_sql_db_and_create_if_necessary(N,args):
    ''' this currently creates only the ps files and and now also checks if those are
    necessary to create '''
    #print bcolors.OKGREEN + "check_for_all_ps_and_create_if_necessary" +bcolors.ENDC
    print(bcolors.OKGREEN + "check_for_all_ps_in_sql_db_and_create_if_necessary" +bcolors.ENDC)
    ####################################################################
    # check for ps xxx files
    ####################################################################
    qpoints_all = get_all_qpoints(args.qvec,args)
    #for i in qpoints_all:
    #    print "hier:",i
    all_ps_files_there = True
    create_db(args.db)
    conn = sqlite3.connect(args.db)
    c = conn.cursor()
    for i in qpoints_all:
        #ps_file = glob.glob("ps"+qpointstring(i)+"_*.dat.lifetimesgood.dat")
        ##print ps_file,len(ps_file)
        #if len(ps_file) < 1:
        #    all_ps_files_there = False
        out = get_row_id(c,conn,ltfq='fq_good',q1=i[0],q2=i[1],q3=i[2],mdsteps=mdsteps)
        #print out,type(out)
        if type(out) != tuple:
            print("missing qpoint:",i)
            all_ps_files_there = False
    conn.close()
    print("all_ps_files_there:",all_ps_files_there)

    if not all_ps_files_there:
        print("need to create ps_files! (or at least one!)")
        sys.exit()
    return

def check_for_all_space_fft_and_rm_xaa_files_folders(N,args):
    print(bcolors.OKGREEN + "check_for_all_space_fft_and_rm_xaa_files_folders" +bcolors.ENDC)
    qpoints_all = get_all_qpoints(args.qvec,args)
    remove_all = True
    #print "qpall:",qpoints_all,type(qpoints_all)
    for i in qpoints_all:
        #print "iii;",i,type(i)
        #print qpointstring(i)
        if not os.path.isfile("space_fft_"+qpointstring(i)+".npy"):
            remove_all  = False
            print(" -> some (or all) space_fft_files are missing ...")
            return False

    #if remove_all == True:
    #    print " -> seems all space_fft_xxx folder are availalbe! -> after removing xaa.. -> get ps"
    #    # get all xaa files:
    #    import shutil
    #    # DONT REMOVE THOSE USER HAS TO DECIDE IF HE WANTS THIS
    #    #xaa_files = glob.glob("x[a-z]?")
    #    #for i in xaa_files:
    #    #    print "removing",i
    #    #    os.remove(i)
    #    #xaa_files = glob.glob("x[a-z]?_wcl")
    #    #for i in xaa_files:
    #    #    print "removing",i
    #    #    os.remove(i)
    #    #xaa_folder= glob.glob("x[a-z]?_")
    #    #for i in xaa_folder:
    #    #    print "removing",i
    #    #    shutil.rmtree(i)
    return remove_all

def get_space_fft_from_xaa_xab_folder_parallel(N,args):
    '''
    ########################################
    # make space_fft
    ########################################'''
    print(bcolors.OKGREEN + "get_space_fft_from_xaa_xab_folder_parallel" +bcolors.ENDC)
    qpoints_all = get_all_qpoints(args.qvec,args)
    num_cores = multiprocessing.cpu_count()
    print("num_cores++:",num_cores)
    #print "qpstring:",qpstring
    #print "qpstring:",qpoints_all
    #sys.exit() # unten
    jobs = []
    for i in qpoints_all:
        print("get_space_fft_from_xaa_xab_folder for",i)
        p = multiprocessing.Process(target=get_space_fft_from_xaa_xab_folder, args=[qpointstring(i)])
        jobs.append(p)
        p.start()
    for job in jobs: # this is just to wait until job is finished
        job.join()
    return True

def check_for_xaa_folder_and_make_everything_else(N,args):
    ''' at this stage we assume we do this since no space_fft files exist. '''
    print(bcolors.OKGREEN + "check_for_xaa_folder_and_make_everything_else" +bcolors.ENDC)
    xaa_file = glob.glob("x??_wcl")
    #print "len(xaa_file):",len(xaa_file)
    if len(xaa_file) == 1: # all xaa files exist; are there also the xaa_ folder?
        xaa_files = glob.glob("x??")
        print(" -> 4 All ("+str(len(xaa_files))+") xaa files exists.")
        xaa_folder = glob.glob("x??_")
        if len(xaa_folder) == 0:
            print(" -> 5 no xaa_ folder, goint to create those and submit to cluster")
            xaa_filenames,xaa_steps = lammps_split_lammps_trajectory_in_xaa_files(filenamelammps="trj_lammps.out",filenamepos="trj_lammpsnew.out",positionsfilename = "positions.*",infilefilename = "in_file_dynamics.in",linesokformemory = 800000000,qpoints_all=qpoints_all,args=args,verbose=args.verbose)
        elif len(xaa_folder) > 0 and len(xaa_folder) == len(xaa_files):
            print(" -> 6 All ("+str(len(xaa_folder))+") xaa_ folder seem to exist -> (in case you want to calculate additional space_fft_xxx files delete xaa_ folder first) check in lastfolder if all sum_all_xxx.dat written")
            print("xaa_file:",xaa_file)
            last_folder = xaa_file[0].split("_")[0]+"_"
            print(len(qpoints_all))
            print("last_folder:",last_folder)
            check_space_ffts = glob.glob(last_folder+"/sum_all_*")
            print(len(check_space_ffts))
            if len(qpoints_all) == len(check_space_ffts):
                print(" -> 7 sum_all_new_xxx files from c++ skript seem to be done")
                if get_space_fft_from_xaa_xab_folder_parallel(N,args):
                    return check_from_space_fft(args)
            elif len(qpoints_all) < len(check_space_ffts):
                print(" -> 7 check if for the chosen qpoint all sum_all_new_xxx files done")
                for i in qpoints_all:
                    if os.path.isfile(last_folder+'/sum_all_new_'+qpointstring(i)+'.dat'):
                        pass
                    else:
                        sys.exit("it seems that "+last_folder+'/sum_all_new_'+qpointstring(i)+'.dat does not exist!')
                if get_space_fft_from_xaa_xab_folder_parallel(N,args):
                    return check_from_space_fft(args)
            else:
                print("len(qpoints_all):",len(qpoints_all),"len(check_space_ffts):",len(check_space_ffts))
                print(" -> it seems that xaa_ folder are already created but the sum_all_new is not yet written, is this job in the que?")
                sys.exit()
    else:
        print("--> 8 grep and split necessary")

    return

def check_from_space_fft(args):
    print(bcolors.FAIL + "check_from_space_fft" +bcolors.ENDC)
    if check_for_all_space_fft_and_rm_xaa_files_folders(args.supercell,args):
        check_for_all_ps_in_sql_db_and_create_if_necessary(args.supercell,args)
        print_and_save_lifetimesgood(args = args,printtoscreen=False)
        return True
    else:
        return False

def check_from_xaa_folder(args):
    print(bcolors.FAIL + "check_from_xaa_folder" +bcolors.ENDC)
    if not check_for_xaa_folder_and_make_everything_else(args.supercell,args):
        if os.path.isfile("trj_lammps.out") or os.path.isfile("trj_lammpsnew.out"):
            print(" -> 9 it seems there are no xaa files! but thre is a trj_lammps(new).out file(s) -> split trajectory and submit c++ to cluster for every xaa folder")
            xaa_filenames,xaa_steps = lammps_split_lammps_trajectory_in_xaa_files(filenamelammps="trj_lammps.out",filenamepos="trj_lammpsnew.out",positionsfilename = "positions.*",infilefilename = "in_file_dynamics.in",linesokformemory = 800000000,qpoints_all=qpoints_all,args=args,verbose=args.verbose)
            sys.exit("I am leaving here since jobs were hopefully submitted to the cluster") # since script above submits to cluster
        else:
            print(" -> 9 no xaa files and no trj_lammps(new).out file(s), returning to skript")
            return



if __name__ == '__main__':
    if args.folder == False:
        folder_all = [os.getcwd()]
    else:
        folder_all = glob.glob(os.getcwd()+"/"+args.folder+"/")
    for idx,i in enumerate(folder_all):
        print("             :",idx+1,"/",len(folder_all),i)
    print()

    ### make README
    create_READMEtxt(os.getcwd())

    for idxfolder,folder_in in enumerate(folder_all):
        os.chdir(folder_in)
        print("#"*(len(os.getcwd())+21))
        print("os.getcwd()  :",idxfolder+1,"/",len(folder_all),os.getcwd())
        print("#"*(len(os.getcwd())+21))

        ########################################
        #### get all variables, qpoints, ...
        #### pint infos to screen
        ########################################
        qpoints_all = get_correct_input_variables(args)
        print("@@@@@@@@@@@@@########@@@@@@@ init done @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        #sys.exit('hhh')
        if type(args.supercell) == bool:
            sys.exit('please specify the supercell size!')
        #print "args:",args
        #sys.exit()

        ########################################
        # check for correct format of sql.tmp database
        ########################################
        if os.path.isfile('sql.tmp'):
            import pickle
            with open('sql.tmp', 'rb') as f:
                try:
                    all = pickle.load(f, encoding='latin1')            # geht ueber alle eintraege
                except TypeError:
                    all = pickle.load(f)            # geht ueber alle eintraege
                for idn,n in enumerate(all):
                    #print 'idn',n
                    if idn == 0:
                        #print idn,n,len(n)
                        if len(n) != 9:
                            sys.exit('sql.tmp has the wrong (old?) format')
                        else:
                            print('sql.tmp has the right lenght')
                    else:
                        pass
                    break

        if args.sql:
            import pickle
            with open('sql.tmp', 'rb') as f:
                try:
                    all = pickle.load(f, encoding='latin1')            # geht ueber alle eintraege
                except TypeError:
                    all = pickle.load(f)            # geht ueber alle eintraege
                for idn,n in enumerate(all):
                    print(idn,n)
                        #if qpoint[0]==n[0] and qpoint[1]==n[1] and qpoint[2]==n[2]:
                        #    print "all idx:",idn,"||","qp:",n[0],n[1],n[2],"mdstep:",n[3],"peaks:",n[4],"ltfq:",n[5],"fq:",n[6],"lt:",n[7]
                        #if qpoint[0]==n[0] and qpoint[1]==n[1] and qpoint[2]==n[2] and cm==n[3] and n[5]!='good': # and filename==n[5]:
                        #    list_xind.append([n[3],n[4],n[5],n[6],n[7]])
                            #print 'nnn',n
                    #print 'cm',cm,len(list_xind)

            sys.exit()

        if args.exit:
            sys.exit('Exit due to specifying the -e (exit) option')

        if type(args.lifetimesaverage) != bool:
            get_lifetimessummary_from_several_folders(base="run_*/nes_of_one_mdstep:")
            sys.exit()



        ####################################################################
        # if not option -ps used, check if space_fft can be created
        ####################################################################
        #if args.make_power_spectrum == False:
        if args.space_fft: #== False:
            check_for_lammps_job()

        #sys.exit()
        if args.space_fft:
            qpstring_all = [];tmp=0;
            for idx,i in enumerate(qpoints_all):
                tmp+=equivalent_qs(i).shape[0]
                qpstring_all.append(qpointstring(i))

            get_space_fft_from_xaa_xab_folder(qpstring=qpstring_all)

        #print "ka",args.mdstepstocalc_all
        #sys.exit()
        ###############################################################
        # get args.mdstepstocalc_all from sql.db if possible
        # so that this has  not to be done in all prallel runs
        ###############################################################
        #if type(args.mdstepstocalc_all) != np.ndarray:
            #print "ka",args.mdstepstocalc_all
            #create_db(args.db) # is only created if it is not existing
            #args.mdstepstocalc_all = get_db_mdstepsall(args.db)

        #sys.exit()
        ##################################################
        # look for filenames POSITIONs to make space fft
        ##################################################
        filename_scale_by_alat_N = [ 'trj_lammpsnew.out', 'trj_lammpsnew.out_noJumps_small', 'pos', 'POSITIONs', 'out_positions_forces_out.dat', 'out_positions_forces.dat', 'out_positions.dat' ] # lammps to direct coords
        filename_scale_by_N = [ 'dum', 'posmichael' ] # michaels sim_fcc_morse (to convert to direct)

        filename = False
        scale_data = False
        for i in filename_scale_by_alat_N:
            if os.path.isfile(i) == True:
                filename = i; scale_data = float(args.alat)*args.supercell
        for i in filename_scale_by_N:
            if os.path.isfile(i) == True:
                filename = i; scale_data = args.supercell
        #print 'filename     :',filename
        #print "scale_data   :",scale_data
        #sys.exit()

        ##############################################################################
        # CHECK IF space_fft_xxxx.npy exist, if not, create it!
        ##############################################################################
        space_fft_which_go = [ 'space_fft_', 'space_fft_phase_' ]
        if args.space_fft_phase == False:
            space_fft_which_go = [ 'space_fft_']
        if args.space_fft_exp == False:
            space_fft_which_go = [ 'space_fft_phase_']
        print('space_fft_which_go:',space_fft_which_go)

        goto_get_space_fft_from_position = False
        if args.make_space_fft_from_py or args.make_space_fft_from_lammpslog or args.make_power_spectrum:
            print("@@@@@@@@@@@@@########@@@@@@@ checking existing fft_phase_ files @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
            for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
                print('idx',idx,'qpoint',qpoint)
                qpointstr = qpointstring(qpoint)
                equivalent_qs(qpoint)
                execute_timeinversion = True
                space_fft_file1 = 'space_fft_'+qpointstr+'.npy'
                space_fft_file1phase = 'space_fft_phase_'+qpointstr+'.npy'
                #@space_fft_file2 = 'sum_all_new_'+qpointstr+'.npy'
                #@space_fft_file3 = 'log_for_sum_all' # created from lammps.log
                #@space_fft_file4 = 'log.lammps'      # lammps.log with comments

                if os.path.isfile(space_fft_file1) and os.path.isfile(space_fft_file1phase):
                    print('file',space_fft_file1,'      exists already!')
                    print('file',space_fft_file1phase,'exists already!')
                else:
                    if not os.path.isfile(space_fft_file1):
                        print('file',space_fft_file1,'      DOES NOT exist!')
                    if not os.path.isfile(space_fft_file1phase):
                        print('file',space_fft_file1phase,'      DOES NOT exist!')
                    goto_get_space_fft_from_position = True

                    #sum_all = np.load(space_fft_file1)
                #@elif os.path.isfile(space_fft_file2):
                #@    print('file',space_fft_file2,'exists already!')
                #@    #print "load existing ",space_fft_file2
                #@    #sum_all = np.load(space_fft_file2)
                #@    execute_timeinversion=False;
                #@elif os.path.isfile(space_fft_file3) == True:
                #@    print("execute lammps_log_to_sum on",space_fft_file3)
                #@    sum_all = lammps_log_to_sum(filename = space_fft_file3, qpoint=qpoint)
                #@elif os.path.isfile(space_fft_file4) == True and args.make_space_fft_from_lammpslog:
                #@    print("execute lammps_log_to_sum on",space_fft_file4)
                #@    sum_all = get_space_fft_from_lammps_log(filename = space_fft_file4)
                #@elif type(filename) != bool and args.make_space_fft_from_py:
                    #sys.exit()
        ##############################################################################
        # in case fftpy needs to be done
        ##############################################################################
        if goto_get_space_fft_from_position and args.make_space_fft_from_py:
            for idxx,qpoint in enumerate(qpoints_all):
                #print('qpofqpall:',qpoint)
                if qpoint_is_longitudinal(qpoint,passhigherBZ=False,args=args):
                    print()
                    qpt = qpoint_map_longitudinal_to_transversal(qpoint,args=args)
                    print('qpt_xxx ',qpt)
                    if qpt is None:
                        qpt = np.array([0,0,0]) # in case of no clue
                    print('idxx',idxx,qpoint,'qpt (line 10688)==> qpt ==>',qpt)
                    for tq in qpt:
                        print('qpoint',qpoint,'tq',tq,'len(qpt)',len(qpt),'qpt',qpt)
            if args.verbose:
                print('filename:',filename)
                print('qpoints_all', qpoints_all)
            if type(filename) == bool:
                sys.exit("No input positins (filename) found")
            if not os.path.exists('POSITIONs'):
                os.symlink(filename,'POSITIONs')
            get_space_fft_from_positions(\
                        filename,\
                        qpoints_all=qpoints_all,\
                        mdsteps_per_chunk = 20000, \
                        scale_data=scale_data, chunkitornot = False, \
                        args = args,\
                        filenames = filename_scale_by_alat_N+filename_scale_by_N)
        #sys.exit('aaa')
        import pickle
        the_filename="sql.tmp"
        if args.make_power_spectrum:
            print("making powerspectrum")
            from multiprocessing import Manager
            import datetime

            dir0 = "ps_full"
            dir1 = "ps_smooth"
            dir2 = "ps_md_convergence"
            dir3 = "ps_fitted"
            dir4 = "ps_save"
            if not os.path.exists(dir0):
                os.makedirs(dir0)
            if not os.path.exists(dir1):
                os.makedirs(dir1)
            if not os.path.exists(dir2):
                os.makedirs(dir2)
            if not os.path.exists(dir3):
                os.makedirs(dir3)
            if not os.path.exists(dir4):
                os.makedirs(dir4)

            #sys.exit()
            ##############################################################################
            # make the powerspektrum and get linewidths/frequency
            ##############################################################################
            for idx,qpoint in enumerate(qpoints_all):  # adding one qpoint after the other is as fast as having multiple
                #sys.exit()
                for space_fft_which in space_fft_which_go:
                    manager = Manager()
                    pppx_add_db= manager.list([])
                    space_fft_to_powerspectrum_ser_par([qpoint], execute_timeinversion=True, args = args,idxx=idx,idxx_sum=len(qpoints_all),space_fft_which=space_fft_which)
                    new = [x for x in pppx_add_db]
                    if os.path.isfile(the_filename):
                        with open(the_filename, 'rb') as f:
                            old = pickle.load(f)
                        for idn,n in enumerate(new):
                            #print "nnn idx:",idn,"||","qp:",n[0],n[1],n[2],"mdstep:",n[3],"peaks:",n[4],"ltfq:",n[5],"fq:",n[6],"lt:",n[7]
                            row=False
                            for ido,o in enumerate(old):
                                if n[0]==o[0] and n[1]==o[1] and n[2]==o[2] and n[3]==o[3] and n[4]==o[4] and n[5]==o[5]:
                                    row=ido
                                    break
                            #print "row:",row
                            if type(row) != bool:  # row was found
                                old[row] = new[idn]
                            else:                  # row was not found
                                old.append(n)
                        new = old
                        with open(the_filename, 'wb') as f:
                            pickle.dump(new, f)
                    else:
                        with open(the_filename, 'wb') as f:
                            pickle.dump(new, f)
        #sys.exit()

        #print "type pppx",type(pppx_add_db)
        #np.savetxt("kkk.dat",pppx_add_db)
        if args.make_power_spectrum or args.lifetimes or args.freqs:
            print_and_save_lifetimesgood(args = args)






