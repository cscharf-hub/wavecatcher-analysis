
#include <iostream>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <TCanvas.h>

using namespace std;

void timing_example(int which = 0) // main
{
	string path;
	// edit for your fs
	path = "examples/exampledata/";

	switch (which) { //specify folder to run below
	default: {
		path += "timingdata/"; 
		break;
	}
	}

	// initialize class
	ReadRun mymeas(0);
	
	// read data
	mymeas.ReadFile(path, true, 0, "timing_example_results.root");

	// smooth data with sigma=300 ps gauss to potentially improve CFD results (dangerous, always check results)
	mymeas.SmoothAll(.3);

	// apply baseline correction to ALL waveforms
	mymeas.CorrectBaselineMin(40, 0, 5., 400, 250, 0);

	// parameters for signal integration
	float intwindowminus = 15.;	// lower integration window in ns rel. to max
	float intwindowplus = 85.;	// upper integration window in ns rel. to max
	float findmaxfrom = 100.;	// assume signal from muon arrives between here ...
	float findmaxto = 140.;		// ... and here (depends on trigger delay setting)
	float rangestart = -100;	// plot range
	float rangeend = 2500;		// plot range


	// filter out pedestal events (usually with pickup RF oscillations)
	mymeas.IntegralFilter({ 100, 100, 200, 200 }, { false, false, false, false }, intwindowminus, intwindowplus, findmaxfrom, findmaxto);

	// get timing for 50% CFD between t=100 ns and t=140 ns, 0 means no smoothing, true to search forward
	mymeas.GetTimingCFD(0.5, 100, 140, 0, true);

	// apply cut for time difference between two channels (ch14 and ch26, events with time differences <1 ns or >5 ns will be cut)
	//mymeas.SkipEventsTimeDiffCut(14, 26, 1, 5, false);

	// plotting
	// plot sums of all events per channel
	mymeas.PlotChannelSums(true, false);

	// plot results between t=100 ns and t=140 ns and fit gauss
	mymeas.Print_GetTimingCFD(100, 140, 1, 100);

	// Plot distribution of time differences between two groups of channels
	mymeas.Print_GetTimingCFD_diff({ 14 }, { 26 }, 0, 8, 2, 100, 0, 8, "S");

	// plot integrated signals of all channels
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, rangestart, rangeend, 100, rangestart, rangeend, 0);
	
	// investigate correlation of signals between channels
	mymeas.ChargeCorrelation(intwindowminus, intwindowplus, findmaxfrom, findmaxto, rangestart, rangeend, 100, 2, 3, true);
	mymeas.ChargeCorrelation(intwindowminus, intwindowplus, findmaxfrom, findmaxto, rangestart, rangeend, 100, 2, 14, true);
	mymeas.ChargeCorrelation(intwindowminus, intwindowplus, findmaxfrom, findmaxto, rangestart, rangeend, 100, 2, 26, true);

	// plot non-skipped waveforms for certain events with integration window and timing info
	// plot range
	double ymin = -10;
	double ymax = 120;
	gROOT->SetBatch(kTRUE); // only write to root file, do not display plots
	for (int i = 1; i < mymeas.nevents; i += static_cast<int>(mymeas.nevents / 50)) {
		if (!mymeas.skip_event[i - 1]) {
			mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, i, ymin, ymax);
		}
	}
	gROOT->SetBatch(kFALSE);
}
