#include <iostream>
#include <fstream>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <TCanvas.h>

using namespace std;

void batch_read_config (int which, int read_ch = -1) // main
{
	// see batch_read() for explanations
	vector<int> Chs;

	string path = "/home/ben/";
	string folder;
	string dataname;
	string path_to_folder;
	string path_to_data;	// save results respective in data folders

	folder = "Dokumente";
	Chs = { read_ch };
	
	switch (which) {
	case(0): {dataname = "S001M001";		break; }
	case(1): {dataname = "S001M001_short";	break; }
	case(2): {dataname = "S001M003";		break; }
	case(3): {dataname = "S001M004";		break; }
	case(4): {dataname = "S001M005";		break; }
	}

	path_to_folder = path + folder + "/";	
	path_to_data = path_to_folder + dataname + "/";
	
	// initialize class
	ReadRun mymeas(0);
	mymeas.start_read_at_channel = read_ch;
	mymeas.end_read_at_channel = read_ch;
	// read data
	// Saving root file the data folder
	mymeas.ReadFile(path_to_data, true, 0, path_to_data + "/Ch_" + to_string(read_ch) + "_" + dataname + ".root");
	
	// only plot certain channels
	mymeas.plot_active_channels = Chs;

	//apply baseline correction to ALL waveforms <- NEEDED but slow when not compiled
	mymeas.CorrectBaselineMinSlopeRMS(50, false, 5, 300, 10, false);

	////plotting
	// investigate "charge" spectrum. should see photo-electron peaks here
	float intwindowminus = 20.;	// lower integration window in ns rel. to max
	float intwindowplus = 65.;	// upper integration window in ns rel. to max
	float findmaxfrom = 100.;	// assume signal from laser arrives between here ...
	float findmaxto = 160.;		// ... and here (depends on trigger delay setting)
	mymeas.PrintChargeSpectrum_pars = { 1e4, 6, .2, 10, 7, 60., -10. };
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -30, 1500, 200, -20, 1000);
	
	// plot waveforms of individual events
	//plot range
	double ymin = -5;
	double ymax = 25;
	// plot waveforms for certain events with integration window
	for (int i = 1; i < mymeas.nevents; i += static_cast<int>(mymeas.nevents / 50)) {
		mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, i, ymin, ymax);
	}
}
