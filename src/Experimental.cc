#include "Experimental.h"

/// @brief Rebin the data to test bandwidth effects. Will combine an integer number of bins into a new, 
/// wider bin and divide the new bin content by the integer number to preserve the shape and integral. 
/// 
/// CAUTION: Not tested with all functions. 
/// 
/// Please note that PlotChannelSums() is calculated while parsing the data - **before** rebinning. 
/// Use PlotChannelAverages() instead.
/// 
/// @param ngroup Integer number of bins to combine.
/// @param sigma_noise Add gaussian noise with ```sigma_noise``` to rebinned data.
/// @param seed Seed for the noise random number generator.
void Experimental::RebinAll(int ngroup, float sigma_noise, unsigned long seed) {

	SP *= static_cast<float>(ngroup);
	SP_inv = 1. / SP;
	binNumber /= ngroup;
	float norm = 1. / static_cast<float>(ngroup);
	cout	<< "\nRebinning the data to a new sampling rate of " << 1. / SP 
			<< " GS/s which corresponds to a bin size of " << SP 
			<< " ns and the data now has " << binNumber << " bins" << endl;

	#pragma omp parallel for
    for (int j = 0; j < nwf; j++) {
        vector<float> &waveform = rundata[j];
        vector<float> rebinned;
		unsigned long thread_seed = seed + j;
    	TRandom3 noise(thread_seed);

        for (size_t i = 0; i + ngroup <= waveform.size(); i += ngroup) {
            float sum = 0.;
            for (int k = 0; k < ngroup; k++) {
                sum += waveform[i + k];
            }
            float avg = sum * norm;

            if (sigma_noise != 0.) {
                avg += noise.Gaus(0., sigma_noise);
            }
            rebinned.push_back(avg);
        }
		rundata[j] = rebinned;
    }
}
/// @example timing_example_rebin.cc

/// @brief Derivative of all waveforms
///	
/// Calculates the forward differences. Reduces the number of bins by one for all waveforms. For testing purposes.
void Experimental::DerivativeAll() {
	cout << "\nForward derivative of all waveforms:" << endl;
	binNumber--; // reduce number of bins 

	#pragma omp parallel for
	for (int j = 0; j < nwf; j++) {
		vector<float> &waveform = rundata[j];
        vector<float> derivative;
		for (int i = 0; i < binNumber; i++) derivative.push_back(waveform[i + 1] - waveform[i]);
		rundata[j] = derivative;
	}
}