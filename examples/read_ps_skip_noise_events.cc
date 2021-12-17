#include <iostream>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <TCanvas.h>

using namespace std;

void read_ps_skip_noise_events() // main
{
	int which = 0; //select meas

	string path;

	// edit for your fs
	path = "/mnt/c/SHiP/data/";

	switch (which) { //specify folder to run below, ALL bin files in this folder will be used
	default: {
		path += "32_Sr_polyst_coated/"; // default
		break;
	}
	}

	// read data
	ReadRun mymeas(0);
	//mymeas.CorrectBaseline(0., 10.);
	mymeas.ReadFile(path, true);
	
	// only plot certain channels
	//mymeas.plot_active_channels = { 9 };

	//apply baseline correction to ALL waveforms <- NEEDED but slow when not compiled
	mymeas.CorrectBaselineMinSlopeRMS(50, 0, 3, 350, 40);

	// print events above a threshold to identify interesting events
	//mymeas.FractionEventsAboveThreshold(-7, false, false, 100, 150, true);

	// remove events with large negative values to remove events which triggered on pick-up noise
	mymeas.SkipEventsPerChannel({ -20, -6, -4 }, true);

	////plotting

	// plot sums of all events per channel
	mymeas.PlotChannelSums(true);

	// investigate charge spectrum. should see photo electron peaks here
	float intwindowminus = 50.;	// lower integration window in ns rel. to max
	float intwindowplus = 100.;	// upper integration window in ns rel. to max
	float findmaxfrom = 90.;	// assume pulse after trigger arrives between here ...
	float findmaxto = 170.;		// ... and here (depends on trigger delay setting etc., for dark counts the signal is random so we look at the whole recorded time range)

	// plot all charge spectrum of channels
	mymeas.PrintChargeSpectrum_pars = { 30, 1e3, 7e4, 2e2 };
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -100, 3000, 400, 300, 2500, 1, 1);

	// timing of maximum
	mymeas.PrintTimeDist(findmaxfrom, findmaxto, findmaxfrom - 5, findmaxto + 5, 120);

	// plot waveforms of individual events
	// plot range
	double ymin = -10;
	double ymax = 25;

	// inspect waveforms where signals are below threshold
	gROOT->SetBatch(kTRUE); // only write to root file
	for (int i = 0; i < mymeas.nevents; i++) {
		if (mymeas.skip_event[i]) mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, i + 1, ymin, ymax);
	}
	gROOT->SetBatch(kFALSE);
}