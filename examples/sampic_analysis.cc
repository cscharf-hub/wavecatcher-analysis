#include <TROOT.h>
#include <iostream>
#include "src/ReadSampic.h"

using namespace std;

// To execute, call: root sampic_analysis.cc 
void sampic_analysis() {
	// adjust the path to your file system
	string path("data/A005/");

	// initialize
	ReadSampic mymeas(0);

	// channels in 2025 testbeam data are: {2, 3, 4, 5, 8, 9, 11, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
	// invert PMT channels with negative amplitudes
	mymeas.switch_polarity_for_channels = {2, 3, 4, 5, 8, 9, 11, 13, 15};

	// read the actual data digitized with SAMPIC
	bool change_polarity = true; // invert signals for channels indicated in switch_polarity_for_channels
	mymeas.ReadFile(path, change_polarity, 0, "output.root");

	// do baseline correction using 5 ns window between 0 and 25 ns
	vector<float> window = {5, 0, 25};
	mymeas.CorrectBaselineMinSlopeRMS(window);

	// plot a signle waveform before any event building etc for checking data
	mymeas.PlotWF(1);
	
	// build events 
	mymeas.EventBuilder(10, {10}, {}, 1, false);


	float intwindowminus = 10.;	// lower integration window in ns rel. to max
	float intwindowplus = 20.;	// upper integration window in ns rel. to max
	float findmaxfrom = 5.;	// assume pulse after trigger arrives between here ...
	float findmaxto = 55.;		// ... and here (depends on trigger delay setting etc., for dark counts the signal is random so we look at the whole recorded time range)
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -100, 1000, 200);

	// cut out event with too low peak integral values
	// mymeas.IntegralFilter({ 250 }, { false }, intwindowminus, intwindowplus, findmaxfrom, findmaxto, false, false);
	
	// plot more stuff in batch mode to root file
	gROOT->SetBatch(kTRUE);
	for (int i = 1; i < mymeas.nevents; i += static_cast<int>(mymeas.nevents / 50)) {
		mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, i, -10, 150);
	}

	
	mymeas.PrintWFProjection(0, 12, -5, 5, 100);
	mymeas.PrintBaselineCorrectionResults(-1100, 300, 1400);

	mymeas.PlotChannelSums();
	mymeas.PlotChannelAverages();
	mymeas.PlotChannelAverages(true);
	
	mymeas.plot_active_channels = { 2, 3, 4, 5, 8, 9, 11, 13, 15 };
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -100, 4000, 200, 0, 0, 99, 0, false);
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -100, 4000, 200, 0, 0, 99, 0, true);
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -100, 15000, 200, 0, 0, 99, 0, false);
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -100, 15000, 200, 0, 0, 99, 0, true);
	mymeas.PlotWFHeatmaps(-10, 300, 200, "log", 0, kRainBow);

	mymeas.plot_active_channels = { 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31 };
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -100, 5000, 200, 0, 0, 99, 0, false);
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -100, 5000, 200, 0, 0, 99, 0, true);
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -100, 25000, 200, 0, 0, 99, 0, false);
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -100, 25000, 200, 0, 0, 99, 0, true);
	mymeas.PlotWFHeatmaps(-10, 400, 200, "log", 0, kRainBow);
	gROOT->SetBatch(kFALSE);
}