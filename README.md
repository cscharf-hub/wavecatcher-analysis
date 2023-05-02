# wavecatcher-analysis

This is the main analysis framework for the SiPM/PMT setups in the high energy physics group at Humboldt University of Berlin which use WaveCatcher devices for digitization.

# Requirements
CERN ROOT Release 6.24/02 or later[^1] recommended. 

It is highly recommended to install [ROOT](https://root.cern/install/#conda) with [conda](https://docs.conda.io/en/latest/miniconda.html). The easiest way to install all dependencies can be achieved by executing the included scripts as explained under [Getting started](#Getting-started).

On Windows please install first WSL and then Ubuntu from the Microsoft Store[^2]. Once done, open Ubuntu and install ROOT with conda.   

# Getting started

To get started open a linux or mac terminal and download repository:
```
git clone https://github.com/cscharf-hub/wavecatcher-analysis
```

Navigate to the repository:
```
cd wavecatcher-analysis
```

Check if you have conda and ROOT installed:
```
conda --version && root --version
``` 
If you don't, please first install conda: 
```
bash etc/scripts/install_miniconda.sh
```
and close and restart the shell once the installation has finished. Now you should see a ```base``` in front of your username in the shell. This is the default conda environment. You can create dedicated environments for the ROOT installation(s) if you like. 
ROOT can be installed with  
```
bash etc/scripts/install_root.sh
```
This will likely take a few minutes, so please be patient.

Now you have all the requirements and can compile[^3] the analysis library: 
```
make
```

This should be it. You can now execute an example macro to test if everything works:
```
root -x examples/read_exampledata.cc
```
or 
```
 root -b -x examples/timing_example.cc -q
```

# Documentation

The documentation can be found here:  
<https://wavecatcher-analysis.web.cern.ch/>

Link to repository:   
<https://github.com/cscharf-hub/wavecatcher-analysis>

# Additional remarks

You should close root with ```.q``` and restart it if you want to re-run a macro. 

If you start root from a different directory you might need to add
```
gSystem->Load("ReadRunLib.sl");
```
to your macros or to rootlogon.C (the rootlogon.C needs to be in the directory from where you start CERN ROOT).

Note that you can add ```-b``` for batch mode or ```-q``` to quit root after running your macro.

# Contact

Support and feature requests:  
christian.scharf at physik.hu-berlin.de


[^1]: Older versions of ROOT are not tested.

[^2]: To open the Microsoft store press the Windows button, type ```store``` and press enter. Now search for ```wsl``` and install it. Repeat with ```ubuntu```. If you encounter issues check [this link](https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-11-with-gui-support).

[^3]: You might need to delete the .o and .sl files by hand before you compile the code another time, depending on read-write permissions.
