# wavecatcher-analysis

# Getting started

To get started download repository and compile the library. To compile the library[^1] on linux or mac do: 
```
make -f makefile
```


Add
```
gSystem->Load("ReadRunLib.sl");
```
to macro or to rootlogon.C (rootlogon.C needs to be in the folder from where you start CERN ROOT).

Then execute 
```
root .x "examples/read_exampledata.cc"
```

Note that you can add ```-b``` for batch mode or ```-q``` to quit root after running the analysis.

You should close root with ```.q``` and restart it if you want to re-run the code because not everything is deleted.

A more detailed documentation can be found here:  
<https://wavecatcher-analysis.web.cern.ch/>  
<https://wavecatcher-analysis.web.cern.ch/classReadRun.html>

Support and feature requests:  
christian.scharf at physik.hu-berlin.de

# Requirements
Needs CERN ROOT Release 6.24/02 or higher compiled with c++17[^2]. 

It is recommended to install ROOT with conda:  
<https://root.cern/install/#conda>

[^1]: You might need to delete the .o and .sl files by hand when you compile the code another time, depending on read-write permissions.
[^2]: If your ROOT has c++14 you need to change line 15 in makefile --std=c++17 to --std=c++14.  
Older versions of ROOT are not supported since they miss Scale() for TGraph.


