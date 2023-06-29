#include "Experimental.h"


/// @brief Rebin the data to test bandwidth effects. Will combine an integer number of bins into a new, 
/// wider bin and divide the new bin content by the integer number to preserve the shape and integral. 
/// See [ROOT::TH1::Rebin()](https://root.cern.ch/doc/master/classTH1.html).
/// 
/// CAUTION: Not tested with all functions. Make sure to adjust bin numbers e. g. in 
/// CorrectBaselineMinSlopeRMS() and bin size e. g. in SmoothArray(). 
/// Know to be incompatible is PlotChannelSums() but PlotChannelAverages() works.
/// 
/// @param ngroup Integer number of bins to combine.
void Experimental::RebinAll(int ngroup) {

	SP *= static_cast<float>(ngroup);
	binNumber /= ngroup;
	cout << "\nRebinning the data to a new sampling rate of " << 1. / SP << " GS/s which corresponds to a bin size of " << SP << " ns and the data now has " << binNumber << " bins\n";
	for (int j = 0; j < nwf; j++) {
		TH1F* his = (TH1F*)rundata->At(j);
		his->Rebin(ngroup);
		his->Scale(1. / static_cast<float>(ngroup));
	}
}
/// @example timing_example_rebin.cc

/// @brief derivative of all waveforms (for measurements w/o pole-zero cancellation)
///
/// Experimental!
void Experimental::DerivativeAll() {
	// just for testing
	cout << "\nderivative of wfs";
	for (int j = 0; j < nwf; j++) {
		TH1F* his = ((TH1F*)rundata->At(j));
		double* yvals = gety(his);
		for (int i = 1; i <= his->GetNbinsX() - 1; i++) his->SetBinContent(i, yvals[i + 1] - yvals[i]);
		delete[] yvals;
		if ((j + 1) % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j + 1) / static_cast<float>(nwf) << "% -" << flush;
	}
}