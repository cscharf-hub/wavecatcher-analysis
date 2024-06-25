import ROOT
import sys
import numpy as np
ROOT.gSystem.Load('ReadRunLib.sl')

sys.path.append('./examples/Freiburg/')
from convert_data import convert_data

### For more examples see read_exampledata.py
###
### To execute, type ```python examples/Freiburg/read_Freiburg_data.py 0```.

def read_Freiburg_data(which, autoclose):
    path_raw = 'examples/Freiburg/'
    path_converted = 'examples/Freiburg/'
    path_results = 'examples/Freiburg/'

    if int(which) == 0:
        path_raw += 'run3_4w_3_20240615175403_164.gz'
        path_converted += 'run3_4w_3_20240615175403_164_converted.bin'
        path_results += 'results.root'
    else:
        print('error: path to data not specified')
        return

    # Initialize class
    mymeas = ROOT.Freiburg_DAQ(0)
    
    # Convert data
    convert_data(path_raw, path_converted)

    # Read data
    mymeas.ReadFile(path_converted, True, path_results)

    # Do some analysis
    baseline_pars = ROOT.std.vector('float')([10., 0, 50.])
    mymeas.CorrectBaselineMinSlopeRMS(baseline_pars)

    mymeas.PrintWFProjection(0, 80, -3.5, 3.5, 50)

    mymeas.PlotChannelAverages()

    # Investigate charge spectrum. should see photo electron peaks here
    intwindowminus = 10.    # lower integration window in ns rel. to max
    intwindowplus = 30.	    # upper integration window in ns rel. to max
    findmaxfrom = 80.	    # assume pulse after trigger arrives between here ...
    findmaxto = 140.		# ... and here (depends on trigger delay setting etc., for dark counts the signal is random so we look at the whole recorded time range)
    
    # Plot range
    ymin = -5
    ymax = 150
    # Save more events to root file
    ROOT.gROOT.SetBatch(True)
    for i in range(1, mymeas.nevents, int(mymeas.nevents / 10)):
           mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, i, ymin, ymax)
    ROOT.gROOT.SetBatch(False)


    # Keep plots open until script is closed manually
    if not autoclose:
        _ = input("Press enter to exit")

def main(argv):
    # Default value of which if no argument is given
    which = 0
    try:
        if len(argv) > 1:
            which = int(argv[1])
        
        autoclose = False
        if len(argv) > 2:
            autoclose = True

        print('reading measurement: ', which)
        read_Freiburg_data(which, autoclose)
    except:
        print('error, check input')

if __name__ == "__main__":
    main(sys.argv)