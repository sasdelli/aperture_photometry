
====================================
Aperture photometry pipeline
====================================

## Installation

How to install the necessary codes.

* Install ureka following the instruction in http://ssb.stsci.edu/ureka/ (The version tested is 1.5.2)

* Clone the git repository on the local system:

```cmd
git clone https://github.com/sasdelli/aperture_photometry
```

## How to run the code:

Guide to use the aperture photometry pipeline.

* Create one folder for each epoch. Move the files with the format "h_e_YYYYMMDD_*.fits". If a different file_format is more convenient edit "write_dictionary.py".

* Run the ureka setup in the folder:
```cmd
ur_setup
```

* Set the XPA_METHOD for ds9, to avoid the error message "unable to verify hostname"
```cmd
export XPA_METHOD=local
```

* Running write_dictionary.py creates a clean dictionary file ("input_dict.txt"). The relevant parameters of "input_dict.txt" can be set editing the beginning of write_dictionary.py.
```cmd
python write_dictionary.py 
```

* Edit the input dictionary "input_dict.txt" (if necessary). The input parameters of the different fits files can be set independently. E.g:
```cmd
vi input_dict.txt
```

* Run the photometry code:
```cmd
python run_photometry.py
```

* If the fit of fitparams fails due to outliers, it might need to run fitparams in interactive mode. Change "interactive='no'" into "yes" within run_photometry.py.

* The final photometries are in the file final_phot_YYYYMMDD

