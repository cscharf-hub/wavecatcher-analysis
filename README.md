# wavecatcher-analysis

This is the main analysis framework used by the experimental elementary (EE) particle physics group at Humboldt University of Berlin for SiPM/PMT setups that use WaveCatcher devices for digitization.

# Table of contents
- [Table of contents](#Table-of-contents)
- [Requirements](#Requirements)
- [Getting started](#Getting-started)
  - [Doing a custom setup](#Doing-a-custom-setup)
  - [On HU EE / DESY NAF computing infrastructure](#On-HU-EE--DESY-NAF-computing-infrastructure)
- [Documentation](#Documentation)
- [Adding custom functions and contribute](#Adding-custom-functions-and-contribute)
- [Additional remarks](#Additional-remarks)

# Requirements
We recommend using CERN ROOT Release 6.24/02 or later[^1].

It is strongly recommended to install [ROOT](https://root.cern/install/#conda) with [conda](https://docs.conda.io/en/latest/miniconda.html), as this will automatically install all dependencies and set the correct paths.
You can do so by simply executing the included scripts as explained under [Getting started](#Getting-started).

The compilation needs ```make```, which can be installed with ```sudo apt install make```. 

If you use Windows, you should first install WSL and then Ubuntu from the Microsoft Store[^2]. 
Once installed, open Ubuntu and install ROOT using conda as described below. 

It should also work on Mac if ROOT was installed using conda, but it is not tested.

# Getting started

Before downloading the repository, navigate to the directory where you want to store the analysis code using ```cd```. 
On WSL, you might prefer to save the analysis in your Windows home folder to make it more easily accessible: ```cd /mnt/c/Users/<your_user_name>/```

## Doing a custom setup
To get started, open a Linux shell and download repository:
```
git clone https://github.com/cscharf-hub/wavecatcher-analysis
```

This will create a new folder in your current directory. Navigate to the downloaded folder:
```
cd wavecatcher-analysis
```

Check if you have conda and ROOT installed:
```
conda --version && root --version
``` 
If not, first install conda:
```
bash etc/scripts/install_miniconda.sh
```
Once installation has finished, restart the shell. Now you should see a ```base``` in front of your username in the shell, which is the default conda environment. You can create dedicated environments for the ROOT installation(s) if you like. 
ROOT can be installed with  
```
bash etc/scripts/install_root.sh
```
This will likely take a few minutes, so please be patient.

Now that you have all the requirements, you can compile[^3] the analysis library: 
```
make
```

You're done! You can now execute an example macro to test if everything works:
```
root "examples/read_exampledata.cc(0)"
```
This will execute the macro ```read_exampledata.cc``` with the parameter ```0```. You can also call: 
```
 root examples/timing_example.cc -q
```
to run the timing example macro and close root once it has finished. The results will be saved in a .root file.

**It is highly advisable to compile your macros** every now and then in order to find bugs. For compilation the includes need to be set correctly and you have to ```#include <src/ReadRun.h>``` at the beginning. Now you can compile the macro by adding a ```+``` after the file name:
```
 root -b -x "examples/read_exampledata.cc+(0)" -q
```

## On the CERN infrastructure

On the CERN infrastructure everything should be set up by default, so you just need to execute:
```
git clone https://github.com/cscharf-hub/wavecatcher-analysis && cd wavecatcher-analysis && make
```


## On HU EE / DESY NAF computing infrastructure

The library can be used efficiently on our computing infrastructure, and the setup takes only a minute since ROOT is already installed there. 
This has the additional benefit that you do not need to download measurements and can do the analysis of large data on more powerful machines. 

Log into the computing cluster via ssh:
```
ssh -X <your-username>@eelg05.physik.hu-berlin.de
```
or ```@naf-atlas.desy.de``` for DESY NAF.

Clone the repository and navigate to the downloaded folder:
```
git clone https://github.com/cscharf-hub/wavecatcher-analysis && cd wavecatcher-analysis
```
Now source ROOT[^4] with either
```
source etc/scripts/root_init_ee.sh
```
at HU or 
```
source /cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh
```
on DESY NAF.

Finally, compile the library:
```
make
```

# Documentation

You can find the documentation for all classes, functions, and variables here: 

<https://wavecatcher-analysis.web.cern.ch/>   

<https://wavecatcher-analysis.web.cern.ch/classReadRun.html>

Documentation for ROOT:

<https://root.cern/doc/master/>

Direct link to histogram documentation:

<https://root.cern/doc/master/group__Histograms.html>

Interactive testing of ROOT ```Draw()``` options:

<https://root.cern/js/latest/examples.htm>

And the link to the repository:   

<https://github.com/cscharf-hub/wavecatcher-analysis>

# Adding custom functions and contribute

The easiest way to add new and/or modified functions for your analysis is to add a derived class which contains your analysis and inherits all functions and parameters from the ```ReadRun``` class. Please note that you will need to add the header file to ```Makefile``` and your new class to ```misc/LinkDef.h```. 

If you find bugs or have developed new, thoroughly tested functionality please make a merge request. 

One example how to structure a derived class can be found here:

<https://wavecatcher-analysis.web.cern.ch/classExperimental.html>

# Additional remarks

To close root, type ```.q``` and restart it if you want to re-run a macro. 

If you start root from a different directory, you might need to add
```
gSystem->Load("ReadRunLib.sl");
```
to your macros or to rootlogon.C (the rootlogon.C file needs to be in the directory from where you start root).

Note that you can add ```-b``` for batch mode or ```-q``` to quit root after running your macro.

The results of the analysis are stored in ```.root``` files. To see the content of the ```.root``` files, start a ```TBrowser```.
You can do so by typing ```rootbrowse```.[^5] Alternatively, you can double-click on the file if you have root installed locally. 
If you don't like the browser-based TBrowser, you can use the old-school TBrowser by calling ```root --web=off```. 

To do less typing you can change the default ```Browser.Name``` in ```$ROOTSYS/etc/system.rootrc``` to ```TRootBrowser``` or create an alias with ```sudo nano ~/.bash_aliases``` and add the line ```alias root='root --web=off'```.
For WSL, you could also add a Windows browser ```alias root='root --web="/mnt/c/Program\ Files\ \(x86\)/.../your-browser.exe"'```.

To update the repository to the latest version, navigate into the downloaded folder and run: ```git pull```. Then compile the repository again ```make```. 
Please note that this does not apply if there are changes in your local repository.

The analysis can be easily used with notebooks[^6]. To open an example notebook, type ```root --notebook``` and navigate to the examples folder.
Jupyter Notebook can be installed with ```pip install notebook``` and then opened with ```jupyter notebook```.

# Contact

Support and feature requests:  
christian.scharf at physik.hu-berlin.de

# Footnotes

[^1]: Older versions of root are not tested.

[^2]: To open the Microsoft store, press the Windows button, type ```store```, and press enter. Now search for ```wsl``` and install it. Repeat with ```ubuntu```. If you encounter issues check [this link](https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-11-with-gui-support).

[^3]: On some systems you might need to call ```make clean``` before compiling another time with ```make```.

[^4]: To avoid repeating this step every time you log in call ```nano ~/.rootrc``` and add line ```source /usr/local/root6/pro/bin/thisroot.sh```.

[^5]: For old ROOT versions start root with ```root``` and then type ```new TBrowser```.

[^6]: On Windows WSL, do ```jupyter notebook --generate-config```, set ```c.NotebookApp.use_redirect_file = True``` and add your default browser ```export BROWSER=<path-to-your-fav-browser>```.
