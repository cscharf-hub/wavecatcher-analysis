# wavecatcher-analysis

This is the main analysis framework for the SiPM/PMT setups in the high energy physics group at HU Berlin which use WaveCatcher devices for digitization.

# Requirements
CERN ROOT Release 6.24/02 or later compiled with c++17[^1] recommended. 

It is highly recommended to install ROOT with conda:  
<https://root.cern/install/#conda>

On Windows please install WSL and Ubuntu following these intructions:   
<https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-11-with-gui-support>

# Getting started

To get started open a linux or mac terminal and download repository with:
```
git clone https://github.com/cscharf-hub/wavecatcher-analysis
```

Navigate to the repository:
```
cd wavecatcher-analysis
```

And compile[^2] the library: 
```
make
```

This should be it. You can now execute an example macro to test if everything works:
```
root -x examples/read_exampledata.cc
```

# Documentation

The documentation can be found here:  
<https://wavecatcher-analysis.web.cern.ch/>  
<https://wavecatcher-analysis.web.cern.ch/classReadRun.html>

# Additional remarks

You need to add
```
gSystem->Load("ReadRunLib.sl");
```
to your macros or to rootlogon.C (rootlogon.C needs to be in the folder from where you start CERN ROOT).

Note that you can add ```-b``` for batch mode or ```-q``` to quit root after running the analysis.

You should close root with ```.q``` and restart it if you want to re-run the code because currently not everything is deleted which might lead to crashes.

Support and feature requests:  
christian.scharf at physik.hu-berlin.de


[^1]: If your ROOT has c++11 (or 14) you need to change line 15 in makefile --std=c++17 to --std=c++11 (or 14).  
Older versions of ROOT are not tested.
[^2]: You might need to delete the .o and .sl files by hand before you compile the code another time, depending on read-write permissions.