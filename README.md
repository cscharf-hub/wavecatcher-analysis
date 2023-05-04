# wavecatcher-analysis

This is the main analysis framework for the SiPM/PMT setups in the experimental elementary (EE) particle physics group at Humboldt University of Berlin which use WaveCatcher devices for digitization.

# Requirements
CERN ROOT Release 6.24/02 or later[^1] recommended.

It is highly recommended to install [ROOT](https://root.cern/install/#conda) with [conda](https://docs.conda.io/en/latest/miniconda.html). The easiest way to install all dependencies can be achieved by executing the included scripts as explained under [Getting started](#Getting-started).

Needs ```make``` for the compilation, which can be installed e. g. with ```sudo apt install make```. 

On Windows please install first WSL and then Ubuntu from the Microsoft Store[^2]. Once done, open Ubuntu and install ROOT with conda as described below. 

# Getting started

Before you download the repository, use ```cd``` to navigate to the directory where you want to store the analysis code. 
On WSL you might want to save the analysis in your Windows home folder to make it more easily accessible: ```cd /mnt/c/Users/<your_user_name>/```

## Doing a custom setup
To get started open a Linux or Mac terminal and download repository:
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
If you don't, please first install conda: 
```
bash etc/scripts/install_miniconda.sh
```
and close and restart the shell once the installation has finished. Now you should see a ```base``` in front of your username in the shell. 
This is the default conda environment. You can create dedicated environments for the ROOT installation(s) if you like. 
ROOT can be installed with  
```
bash etc/scripts/install_root.sh
```
This will likely take a few minutes, so please be patient.

Now you should have all the requirements and can compile[^3] the analysis library: 
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

## On HU EE computing infrastructure

The library can be used efficiently on our computing infrastructure and the setup only takes a minute since ROOT is already installed there. 
This has the additional benefit that you do not need to download measurements and can do the analysis of large data on more powerful machines. 
Log into the computing cluster via ssh:
```
ssh -X <your-username>@eelg05.physik.hu-berlin.de
```
Clone the repository and navigate to the downloaded folder:
```
git clone https://github.com/cscharf-hub/wavecatcher-analysis && cd wavecatcher-analysis
```
Now source ROOT[^4]:
```
source etc/scripts/root_init_ee.sh
```
And finally compile the library:
```
make
```

# Documentation

The documentation of all functions and variables can be found here:  
<https://wavecatcher-analysis.web.cern.ch/>   
<https://wavecatcher-analysis.web.cern.ch/classReadRun.html>

The documentation of ROOT can be found here:   
<https://root.cern/doc/master/>

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

The results of the analysis are stored in ```.root``` files. In order to see the content of the ```.root``` files you need to start a ```TBrowser```.
You can do so by starting ROOT with ```root``` and then typing ```new TBrowser```. You can, of course, also double click on the file if you have root installed locally. 
If you don't like the browser-based TBrowser you can use the old-school TBrowser by calling ```root --web=off```. 

To save typing, you can create an alias with ```sudo nano ~/.bash_aliases``` and add the line ```alias root='root --web=off'```. 
For WSL you could also add a Windows browser ```alias root='root --web="/mnt/c/Program\ Files\ \(x86\)/.../your-browser.exe"'```.

To update the repository to the latest version, navigate into the folder of the repository and type ```git pull``` and then compile it again ```make```. 
Please note that this does not apply if there are changes in your local repository.

The analysis can be easily used with Jupyter Notebook[^5] (see examples). 
It can be installed e. g. with ```pip install notebook``` and open with ```jupyter notebook```.

# Contact

Support and feature requests:  
christian.scharf at physik.hu-berlin.de

# Footnotes

[^1]: Older versions of ROOT are not tested.

[^2]: To open the Microsoft store press the Windows button, type ```store``` and press enter. Now search for ```wsl``` and install it. Repeat with ```ubuntu```. If you encounter issues check [this link](https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-11-with-gui-support).

[^3]: On some systems you might need to call ```make clean``` before compiling another time with ```make```.

[^4]: To avoid repeating this step every time you log in call ```nano ~/.rootrc``` and add line ```source /usr/local/root6/pro/bin/thisroot.sh```.

[^5]: On Windows WSL, do ```jupyter notebook --generate-config```, set ```c.NotebookApp.use_redirect_file = True``` and add your default browser ```export BROWSER=<path-to-your-fav-browser>```.
