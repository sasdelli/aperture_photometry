from astropy.coordinates import ICRS
from astropy import units as u
from astropy.wcs import WCS
from numpy import loadtxt, shape, mean, sort, savetxt, size, genfromtxt, copy, nan, nanmax, nanmin, isnan, sqrt, array, cos, pi
from pylab import figure
from matplotlib.pyplot import plot, savefig, xlabel, ylabel, scatter, axis, xlim, fill_between, hist
import os
import glob

def make_clean():
    temp_files = ['*.coo.*','*.mag.*','obsout','fobsout.dat','image.mag','fitparams_output','final_obsout*', '*.pdf', 'imagelist', 'shortnames']
    files=[]
    for f in temp_files:
        [files.append(ff) for ff in glob.glob(f)]
    for filename in files:
        os.remove(filename)
    return


########

make_clean()
