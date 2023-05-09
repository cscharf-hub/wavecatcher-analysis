#include <TROOT.h>
#include "src/ReadRun.h"

// call .x read_exampledata.cc(0) from root prompt or root -x "read_exampledata.cc(0)"
void read_exampledata(int which = 0) 
{
	string path("examples/");

	// select measurement data. ALL .bin files in this folder will be used!
	switch (which) { 
	case(0): {
		path += "exampledata/";
		break;
	}
	default: {
		cout << "\nerror: path to data not specified" << endl; // default
		break;
	}
	}

	// initialize class
	ReadRun mymeas(0);

	// read data
	mymeas.ReadFile(path, true, 8, "examples/exampleresults.root");

	// uncomment lines below to only plot trigger channels
	//mymeas.plot_active_channels.push_back(14);
	//mymeas.plot_active_channels.push_back(15);
	
	// apply baseline correction to ALL waveforms
	// CorrectBaseline() is the most simple method for data without dark counts (PMT data). Uses mean from 5 ns to 55 ns
	mymeas.CorrectBaselineMinSlopeRMS(100, false, 0, 250, 50);

	////plotting
	// plot sums of all raw events per channel (see channel 9 has an offset)
	mymeas.PlotChannelSums();

	// investigate charge spectrum. should see photo electron peaks here
	float intwindowminus = 10.;	// lower integration window in ns rel. to max
	float intwindowplus = 30.;	// upper integration window in ns rel. to max
	float findmaxfrom = 80.;	// assume pulse after trigger arrives between here ...
	float findmaxto = 140.;		// ... and here (depends on trigger delay setting etc., for dark counts the signal is random so we look at the whole recorded time range)

	// get a rough estimate of the timing of the main peaks to make sure the choice of the time window makes sense
	mymeas.PrintTimeDist(findmaxfrom, findmaxto, findmaxfrom - 5, findmaxto + 5, 60, 1, .5);

	// plot charge spectra (histogram of signal integrals)
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -100, 2.5e3, 200);

	// cut out pedestal events by setting a threshold of 200 mV*ns of the integrals of the last two trigger channels
	// note that this removes about 70% of all events for this example data!!
	mymeas.IntegralFilter({ 0, 200, 200 }, { false, false, false }, intwindowminus, intwindowplus, findmaxfrom, findmaxto);

	// plot the average corrected waveforms per channel not taking into account events cut out by IntegralFilter()
	mymeas.PlotChannelAverages();

	// plot waveforms of individual events
	int example_event = 68;
	//plot range
	double ymin = -5;
	double ymax = 50;
	// plot the waveforms for event 68 with integration window and baseline correction info
	mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, example_event, ymin, ymax);

	// save more events to root file
	gROOT->SetBatch(kTRUE);
	for (int i = 1; i < mymeas.nevents; i += static_cast<int>(mymeas.nevents / 10)) {
		mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, i, ymin, ymax);
	}
	gROOT->SetBatch(kFALSE);

	// now select only the signal channel 9
	mymeas.plot_active_channels = { 9 };
	
	// set starting values for the fit of a landau gauss convolution
	mymeas.PrintChargeSpectrum_pars = { 100, 9.5e3, 2e4, 700 };

	// fit and plot the fit results for channel 9
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, 6e3, 1.35e4, 200, 7e3, 1.3e4, 9, 1);
}