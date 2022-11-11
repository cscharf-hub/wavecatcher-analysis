# wavecatcher-analysis

\section Getting started

To get started download repository and compile the library. To compile the library on linux or mac do: 
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

You should close root (".q") and restart it if you want to re-run the code because not everything is deleted.

A more detailed documentation can be found here:
<https://wavecatcher-analysis.web.cern.ch/>
<https://wavecatcher-analysis.web.cern.ch/classReadRun.html>

Feel free to contact <christian.scharf@physik.hu-berlin.de> for questions and feature requests.

\section Requirements
Needs CERN ROOT Release 6.24/02 or higher compiled with c++17 (otherwise need to change line 15 in makefile --std=c++17 -> --std=c++14).

It is recommended to install ROOT with conda:
<https://root.cern/install/#conda>
