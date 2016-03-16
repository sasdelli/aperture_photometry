from astropy.coordinates import ICRS
from astropy import units as u
from astropy.wcs import WCS
from numpy import loadtxt, shape, mean, sort, savetxt, size, genfromtxt, copy, nan, nanmax, nanmin, isnan, sqrt, array, cos, pi, zeros, median
from pylab import figure
from matplotlib.pyplot import plot, savefig, xlabel, ylabel, scatter, axis, xlim, fill_between, hist, gca, legend
import os
import glob
from pyraf import iraf

def make_clean():
    temp_files = ['*.coo.*','*.mag.*','obsout','fobsout.dat','image.mag','fitparams_output','final_obsout*', '*.pdf', 'imagelist', 'shortnames']
    files=[]
    for f in temp_files:
        [files.append(ff) for ff in glob.glob(f)]
    for filename in files:
        os.remove(filename)
    return

def mkshortnames(string_in=''):
    file_list = glob.glob('./*'+string_in+'*.fits');file_list.sort()
    outfile = open('./shortnames','w');outfile.write('obj : '); outfile.write(" ".join([filen[2:25] for filen in file_list]) + '\n'); outfile.close()
    return


def write_file_list(string_in=''):
    file_list = glob.glob('./*'+string_in+'*.fits');file_list.sort()
    outfile = open('./imagelist','w'); outfile.write("\n".join(file_list) + '\n'); outfile.close()
    return

def run_sextractor_daofind(imagelist):
    # from http://tdc-www.harvard.edu/wcstools/sextractor/
    import subprocess
    file_list = loadtxt(imagelist,dtype='str')
    for f in file_list:
        subprocess.call("sex -c daofind.sex %s" % (f), shell=True)
        subprocess.call("mv test.cat %s.coo.1" % (f), shell=True)
    return


def fix_headers(string_in=''):
    file_list = glob.glob('./*'+string_in+'*.fits');file_list.sort()
    for name in file_list:
        fix_header(name)
    return

def find_band(date, band='G'):
    """ This finds the fits file with G band.
    """
    import pyfits
    file_list = glob.glob('*'+date+'*.fits');file_list.sort()
    for fits_file in file_list:
        hdulist = pyfits.open(fits_file)
        prihdr = hdulist[0].header
        filter = prihdr['filter1']
        if filter == band:
            file_band = fits_file
    return file_band

def fix_header(fits_file):
    """ This changes the name of the filter in the header.
    """
    import pyfits
    hdulist = pyfits.open(fits_file, mode='update')
    prihdr = hdulist[0].header
    filter = prihdr['filter1']
    changes = [[ 'Bessell-B', 'B'],
               [ 'Bessell-V', 'V'],
               [ 'SDSS-U', 'U'],
               [ 'SDSS-G', 'G'],
               [ 'SDSS-R', 'R'],
               [ 'SDSS-I', 'I'],
               [ 'SDSS-Z', 'Z']]
    for i in changes:
         if filter == i[0]:
             prihdr['filter1'] = i[1]
             hdulist.flush()
    return

def match_Stars(fits_file, save_location='./', nobsfile='obsout', standard_stars_file='standard_stars.dat', out_nobsfile='final_obsout', header_line=5, SN_coord=None, dist_treshold = 0.0004, sel_3D=True ):
    '''- fits_file is the full path of the reference fits file.
    - save_location is the directory in which the output will be saved.
    - nobsfile is the full path of the nobsfile.
    - standard_stars_file is the full path of the standard stars file.
    - phot_band is the photometric reference band.
    - index_band is the indec of the band in the standard_stars_file.
    - offset is an integer which changes the nobsfile line that is used to name the star. e.g, if g is the reference but the band order in 
    nobsfile is u then g the offset should be'1'
    - header_line is the line number of the header of the standard_stars_file
    - dist_treshold the match is discarded if the stars are farther apart than this value (degrees)
    
    The output is four images showing the matched stars plus 'final_obsout', which is the nobsfile with the star names
    in place.'''

    import pyfits
    hdulist = pyfits.open(fits_file)
    prihdr = hdulist[0].header
    phot_band = prihdr['filter1']
    index_band = phot_band.lower()
    print phot_band

    if SN_coord==None:
        SN_c = ICRS(prihdr['RA']+prihdr['DEC'], unit = (u.hourangle, u.degree))
        SN_coord = [SN_c.ra.degree, SN_c.dec.degree]

    obsout_file = nobsfile
    gri_file = standard_stars_file

    if os.path.isfile(fits_file) == True: 
        w = WCS(fits_file)
        ##
        obsout = genfromtxt(obsout_file,dtype=None, missing_values='INDEF')
        gri = genfromtxt(gri_file,dtype=None, missing_values='INDEF', skip_header=header_line-1,names=True)

        ##
        def marker_size_(in_mag):
            return 20*(nanmax(in_mag)-in_mag)/(nanmax(in_mag)-nanmin(in_mag))
        def select_phot_(in_v,phot_out,phot_list):
            out_v = copy(in_v)
            for i in range(size(out_v)): 
                if phot_list[i]!=phot_out:
                    out_v[i]=nan
            return out_v
        ##
        figure(figsize=(8,8))
        # all_pix2world is from http://docs.astropy.org/en/stable/wcs/
        lon, lat =w.all_pix2world(obsout['f4'], obsout['f5'], 0)
        
        lon = select_phot_( lon , phot_out=phot_band, phot_list=obsout['f1'])
        lat = select_phot_( lat , phot_out=phot_band, phot_list=obsout['f1'])
        scatter(lon,lat,marker='o',s=marker_size_(obsout['f6'])  )
        scatter(gri['RA'],gri['dec'],color='r',marker='o',s=marker_size_(gri[index_band])    )
        plot(SN_coord[0],SN_coord[1],'+g')
        gca().invert_yaxis()
        savefig(save_location+'out1.pdf')
        #
        # match_to_catalog_sky is from https://www.sites.google.com/site/mrpaulhancock/blog/theage-oldproblemofcross-matchingastronomicalsources
        
        indx1 = [i for i in range(size(lon)) if not isnan(lon*lat)[i]]
        
        cat1_ra = lon[indx1]
        cat1_dec = lat[indx1]
        cat2_ra =  gri['RA']
        cat2_dec = gri['dec']

        cat3_ra = [SN_coord[0]]
        cat3_dec = [SN_coord[1]]
        
        cat1 = ICRS(cat1_ra,cat1_dec,unit = (u.degree,u.degree))
        cat2 = ICRS(cat2_ra,cat2_dec,unit = (u.degree,u.degree))
        cat3 = ICRS(cat3_ra,cat3_dec,unit = (u.degree,u.degree))
        index_list,dist2d,dist3d = cat2.match_to_catalog_sky(cat1)

        index_list_SN,dist2d_3,dist3d_3 = cat3.match_to_catalog_sky(cat1)
        
        indexes_gri = []
        indexes_obsout = []
        distances=[]
        for i in range(size(gri)):
            j = indx1[index_list[i]]
            cos_dec = cos(lat[j]*pi/(180.))
            dist_ = ( sqrt(sum(array([(lon[j] - gri['RA'][i])*cos_dec  ,   lat[j] - gri['dec'][i]])**2)) )
            if dist_ < dist_treshold:
                indexes_gri.append(i)
                indexes_obsout.append(j)
            distances.append(dist_)

        if sel_3D:

            valuesA_=[]
            valuesB_=[]
            for i in indexes_gri:
                j = indx1[index_list[i]]
                valuesA_.append(obsout['f6'][j])
                valuesB_.append(gri[index_band][i])
            valuesA_=array(valuesA_);valuesB_=array(valuesB_) 
            sigma_res =  median(abs(valuesA_- valuesB_ - median(valuesA_ - valuesB_)))
            res_ =  abs(valuesA_- valuesB_ - median(valuesA_ - valuesB_))
    
            indexes_gri_new = []
            indexes_obsout = []
            distances=[]
            ii=0
            for i in indexes_gri:
                j = indx1[index_list[i]]
                cos_dec = cos(lat[j]*pi/(180.))
                dist_ = ( sqrt(sum(array([(lon[j] - gri['RA'][i])*cos_dec  ,   lat[j] - gri['dec'][i]])**2)) )
                if dist_ < dist_treshold:
                    if not sel_3D or res_[ii] < 4.*sigma_res:
                        indexes_gri_new.append(i)
                        indexes_obsout.append(j)
                distances.append(dist_)
                ii+=1
            indexes_gri = indexes_gri_new

            indexes_gri_new = []
            indexes_obsout = []
            distances=[]
            ii=0
            for i in indexes_gri:
                j = indx1[index_list[i]]
                cos_dec = cos(lat[j]*pi/(180.))
                dist_ = ( sqrt(sum(array([(lon[j] - gri['RA'][i])*cos_dec  ,   lat[j] - gri['dec'][i]])**2)) )
                if dist_ < dist_treshold:
                    if not sel_3D or res_[ii] < 3.*sigma_res:
                        indexes_gri_new.append(i)
                        indexes_obsout.append(j)
                distances.append(dist_)
                ii+=1
            indexes_gri = indexes_gri_new

        indexes_obsout_SN = nan
        j = indx1[index_list_SN[0]]
        cos_dec = cos(lat[j]*pi/(180.))
        dist_ = ( sqrt(sum(array([(lon[j] - SN_coord[0])*cos_dec  ,   lat[j] - SN_coord[1]])**2)) )
        if dist_ < dist_treshold:
            #indexes_gri.append(i)
            indexes_obsout_SN = j

        figure()
        hist(distances,1000,label='catalog stars')
        plot([dist_ ,dist_],[0,1],'g', linewidth=10, label='SN')
        plot([dist_treshold,dist_treshold],[0,3],':k', label='dist_treshold')
        xlim((0,0.002))
        xlabel('distance (degree)')
        ylabel('N stars')
        legend()
        savefig(save_location+'out2.pdf')
        #
        figure(figsize=(8,8))
        offset = list(obsout['f1']).index(phot_band)
        for i in indexes_gri:
            j = indx1[index_list[i]]
            obsout[j-offset][0]=gri[i][0]
            scatter(lon[j],lat[j],marker='o',s=marker_size_(obsout['f6'])[j]  )
            scatter(gri['RA'][i],gri['dec'][i],color='r',marker='o',s=marker_size_(gri[index_band])[i]    )
            plot(SN_coord[0],SN_coord[1],'+g')
            plot([lon[j],gri['RA'][i]],[lat[j],gri['dec'][i]])
        if not isnan(indexes_obsout_SN):
            j=indexes_obsout_SN
            obsout[j-offset][0] = 'SN'
            scatter(lon[j],lat[j],marker='o',color='g',s=marker_size_(obsout['f6'])[j]  )
        gca().invert_yaxis()
        savefig(save_location+'out3.pdf')
        #
        figure()
        for i in indexes_gri:
            j = indx1[index_list[i]]
            plot( obsout['f6'][j],gri[index_band][i],'bo')
        savefig(save_location+'out4.pdf')
        
        ## Save final_obsout            

        f=open(save_location+out_nobsfile,'w+')
        for j in range(len(obsout)):
            for k in range(len(obsout[j])):
                if k <len(obsout[j])-1:
                    f.write(str(obsout[j][k]).replace('nan','INDEF').replace('False','INDEF'))
                    for h in range(13-len(str(obsout[j][k]).replace('nan','INDEF').replace('False','INDEF'))):
                        f.write(' ')
                else:
                    f.write(str(obsout[j][k]).replace('nan','INDEF').replace('False','INDEF')+'\n')
        f.close()    

    else:
        print fits_file + 'not found'
    return


