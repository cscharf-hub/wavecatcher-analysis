#include <TROOT.h>
#include "src/ReadRun.h"

// main function, needs to have the same name as the file
void timing_example() 
{
	// path to the data
	string path("examples/exampledata/timingdata/");
	
	// initialize class
	ReadRun mymeas(0);
	
	// read data
	mymeas.ReadFile(path, true, 0, "examples/timing_example_results.root");

	// apply baseline correction 
	// CorrectBaselineMin() is optimized for SiPM measurements with a high dark count rate
	// searches for the minimum sum over 40 bins without smoothing (0.) from bin 250 (78 ns) to bin 400 (125 ns) 
	mymeas.CorrectBaselineMin(40, 0., 400, 250);

	// parameters for signal integration
	float intwindowminus = 15.;	// lower integration window in ns relative to max
	float intwindowplus = 85.;	// upper integration window in ns relative to max
	float findmaxfrom = 100.;	// assume signal from muon arrives between here ...
	float findmaxto = 140.;		// ... and here (depends on trigger delay setting)
	float rangestart = -100;	// plot range
	float rangeend = 2500;		// plot range
	int number_of_bins = 100;	// number of bins in the histogram (choose according to plot range and number of entries)

	// plot integrated signals of all channels
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, rangestart, rangeend, number_of_bins, rangestart, rangeend);

	// investigate correlation of signals between channels
	mymeas.ChargeCorrelation(intwindowminus, intwindowplus, findmaxfrom, findmaxto, rangestart, rangeend, 100, 2, 3, true);
	mymeas.ChargeCorrelation(intwindowminus, intwindowplus, findmaxfrom, findmaxto, rangestart, rangeend, 100, 2, 14, true);
	mymeas.ChargeCorrelation(intwindowminus, intwindowplus, findmaxfrom, findmaxto, rangestart, rangeend, 100, 2, 26, true);

	// filter out pedestal events by setting a threshold of 100 mV*ns of the integrals of the first two channels and 200 mv*ns of the other two
	mymeas.IntegralFilter({ 100, 100, 200, 200 }, { false, false, false, false }, intwindowminus, intwindowplus, findmaxfrom, findmaxto);

	// plot integrated signals of all channels again to see effect of IntegralFilter 
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, rangestart, rangeend, number_of_bins, rangestart, rangeend);

	// plot the average waveforms of all events that are not skipped per channel
	mymeas.PlotChannelAverages();

	// get timing for 30% CFD between t=100 ns and t=140 ns, 0 means no smoothing, true to start search at t=100 ns
	mymeas.GetTimingCFD(0.3, 100, 140, 0, true);
	
	// plot results between t=100 ns and t=140 ns and fit gauss
	mymeas.Print_GetTimingCFD(100, 140, 1, 100);

	// plot distribution of time differences between two groups of channels in range 0 ns to 8 ns with 100 bins 
	// and fit a gauss in range 1 ns to 6 ns
	mymeas.Print_GetTimingCFD_diff({ 14 }, { 26 }, 0, 8, 1, 100, 1, 6, "S");

	// apply cut for time difference between two channels (ch14 and ch26, events with time differences <1.5 ns or >5 ns will be discarded)
	mymeas.SkipEventsTimeDiffCut(14, 26, 1.5, 5);

	// plot non-skipped waveforms for certain events with integration window and timing info
	// plot range
	double ymin = -10;
	double ymax = 120;
	// activate batch mode: only write plots to root file, do not display plots immediately (would be a lot of plots)
	gROOT->SetBatch(kTRUE); 
	for (int i = 1; i < mymeas.nevents; i += static_cast<int>(mymeas.nevents / 50)) {
		if (!mymeas.skip_event[i - 1]) {
			mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, i, ymin, ymax);
		}
	}
	// deactivate batch mode
	gROOT->SetBatch(kFALSE);
}
