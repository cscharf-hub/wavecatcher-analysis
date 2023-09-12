import ROOT
import sys
import numpy as np
ROOT.gSystem.Load('ReadRunLib.sl')

### Example how to perform full analysis in Python
###
### To execute, type ```python read_exampledata.py 0```.
###
### Since the library is compiled the analysis should be fast also in Python. 
### However, at the moment of writing, there are some constraints with using Python:
### Some functions rely on C++ vectors which are similar to list and numpy array, but not exactly the same. 
### Vectors can hold only one type and the interpreter sometimes converts to the wrong types. 
### To avoid errors, always cast to ROOT.std.vector (see code below).

def read_exampledata(which, autoclose):
    path = 'examples/'
    if int(which) == 0:
        path += 'exampledata/'
    else:
        print('error: path to data not specified')
        return

    # Initialize class
    mymeas = ROOT.ReadRun(0)

    # Read data
    mymeas.ReadFile(path, True, 8, "examples/exampledata_results.root")

    # Apply baseline correction to ALL waveforms
    # Searches for the minimum of sum((y_{i+1} - y_{i})^2)+sum(y_{i+1} - y_{i})^2 over 30 ns, starting at t=0 ns until t=80 ns
    # The interpreter has problems converting python list objects into c++ vector objects
    # To do a proper cast, one has to use the ROOT.std.vector() constructor
    baseline_pars = ROOT.std.vector('float')([30., 0, 80.])
    mymeas.CorrectBaselineMinSlopeRMS(baseline_pars)
    
    # test if baseline correction worked (should be centered around 0)
    mymeas.PrintWFProjection(0, 80, -3.5, 3.5, 50);
    
    # print result baseline_correction_result
    print('baseline correction result: ', mymeas.baseline_correction_result[1][2])
    
    ##Plotting
    # Plot sums of all raw events per channel (see channel 9 has an offset)
    mymeas.PlotChannelSums()

    # Investigate charge spectrum. should see photo electron peaks here
    intwindowminus = 10.    # lower integration window in ns rel. to max
    intwindowplus = 30.	    # upper integration window in ns rel. to max
    findmaxfrom = 80.	    # assume pulse after trigger arrives between here ...
    findmaxto = 140.		# ... and here (depends on trigger delay setting etc., for dark counts the signal is random so we look at the whole recorded time range)

    # Cut out pedestal events by setting a threshold of 200 mV*ns of the integrals of the last two trigger channels
    # Note that this removes about 70% of all events for this example data!!
    # The interpreter has problems converting python list objects into c++ vector objects
    # To do a proper cast, one has to use the ROOT.std.vector() constructor
    highlow = ROOT.std.vector('bool')([False, False, False]);
    mymeas.IntegralFilter(ROOT.std.vector('float')([0., 200., 200.]), highlow, intwindowminus, intwindowplus, findmaxfrom, findmaxto)

    # Get a rough estimate of the timing of the main peaks to make sure the choice of the time window makes sense
    mymeas.PrintTimeDist(50, 170, findmaxfrom, findmaxto, 60, 1, .5)

    # Plot the average corrected waveforms per channel not taking into account events cut out by IntegralFilter()
    mymeas.PlotChannelAverages()
    
    # Plot waveforms of individual events
    example_event = 68
    # Plot range
    ymin = -5
    ymax = 150
    # Plot the waveforms for event 68 with integration window and baseline correction info
    mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, example_event, ymin, ymax)

    # Save more events to root file
    ROOT.gROOT.SetBatch(True)
    for i in range(1, mymeas.nevents, int(mymeas.nevents / 10)):
           mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, i, ymin, ymax)
    ROOT.gROOT.SetBatch(False)

    # Now select only the signal channel 9
    # The only safe way to write to existing c++ vectors seems to be vector::push_back()
    mymeas.plot_active_channels.push_back(9)
	
    ## Set starting values for the fit of a landau gauss convolution
    mymeas.PrintChargeSpectrum_pars.push_back(100); 
    mymeas.PrintChargeSpectrum_pars.push_back(9.5e3); 
    mymeas.PrintChargeSpectrum_pars.push_back(2e4); 
    mymeas.PrintChargeSpectrum_pars.push_back(700); 

    # Fit and plot the fit results for channel 9
    mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, 6e3, 1.35e4, 200, 7e3, 1.3e4, 9, 1)

    # Keep plots open until script is closed manually
    if not autoclose:
        _ = input("Press enter to exit")

def main(argv):
    # Default value of which if no argument is given, i. e. ```python read_exampledata.py```
    which = 0
    try:
        if len(argv) > 1:
            which = int(argv[1])
        
        autoclose = False
        if len(argv) > 2:
            autoclose = True

        print('reading measurement: ', which)
        read_exampledata(which, autoclose)
    except:
        print('error, check input')

if __name__ == "__main__":
    main(sys.argv)