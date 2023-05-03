// Helper functions 
#include "ReadRun.h"

/// @brief Helper. Creates a list of .bin data files in data folder to be read in
/// @param dirname Directory
/// @param ext File extension
/// @return String of line separated file names
string ReadRun::list_files(const char* dirname, const char* ext) {

	stringstream ss;
	TSystemDirectory dir(dirname, dirname);
	TList* files = dir.GetListOfFiles();
	if (files) {
		TSystemFile* file;
		TString fname;
		TIter next(files);
		while ((file = (TSystemFile*)next())) {
			fname = file->GetName();
			if (!file->IsDirectory() && fname.EndsWith(ext)) {
				ss << fname.Data() << "\n";
				//cout << fname.Data() << "\n";
			}
		}
		TIter next2(files);
		while ((file = (TSystemFile*)next2())) {
			fname = file->GetName();
			if (!file->IsDirectory() && !fname.EndsWith(ext) && fname.Contains(ext)) {
				ss << fname.Data() << "\n";
				//cout << fname.Data() << "\n";
			}
		}
	}
	return ss.str();
}

/// @brief Helper that returns the waveform histogram for a certain channel number and a certain event number
/// @param channelnr Channel number index (not the actual channel number)
/// @param eventnr Event number
/// @param color Choose color of histogram
/// @return Waveform histogram 
TH1F* ReadRun::Getwf(int channelnr, int eventnr, int color) {
	TH1F* his;
	his = (TH1F*)rundata->At(eventnr * nchannels + channelnr);
	his->SetLineColor(color);
	his->SetMarkerColor(color);
	return his;
}

/// @brief Get array of x axis (time) for standard wavecatcher settings
/// @param shift Offset
/// @return Time array
double* ReadRun::getx(double shift) {
	double* xvals = new double[binNumber];
	for (int i = 0; i < binNumber; i++) {
		xvals[i] = static_cast<double>(SP) * static_cast<double>(i) + shift;
	}
	return xvals;
}

/// @brief Get array of y values for a certain waveform
/// @param channelnr Channel number index (not the actual channel number)
/// @param eventnr Event number
/// @return Y values of waveform
double* ReadRun::gety(int channelnr, int eventnr) {
	TH1F* his = Getwf(channelnr, eventnr);
	double* yvals = new double[binNumber];
	for (int i = 0; i < his->GetNbinsX(); i++) {
		yvals[i] = his->GetBinContent(i);
	}
	return yvals;
}

/// @brief Get array of y values for a certain waveform
/// @param his Waveform histogram
/// @return Y values of waveform
double* ReadRun::gety(TH1F* his) {
	double* yvals = new double[binNumber];
	for (int i = 0; i < his->GetNbinsX(); i++) {
		yvals[i] = his->GetBinContent(i);
	}
	return yvals;
}

/// @brief Get truncated array of y values for a certain waveform
/// @param his Waveform histogram
/// @param start_at Truncate from index
/// @param end_at Truncate to index
/// @return Truncated Y values of waveform
double* ReadRun::gety(TH1F* his, int start_at, int end_at) {
	if (start_at < 0 || start_at >= his->GetNbinsX() || end_at >= his->GetNbinsX() || end_at - start_at < 1) {
		cout << "\nError: ReadRun::gety out of range" << endl;
		return 0;
	}
	const int n_bins_new = end_at - start_at;
	double* yvals = new double[n_bins_new];
	for (int i = start_at; i < end_at; i++) {
		yvals[i - start_at] = his->GetBinContent(i);
	}
	return yvals;
}

/// @brief Translate a random number into a useful root color https://root.cern.ch/doc/master/classTColor.html
/// @param i Index of your plotting loop that is to be translated into a useful ROOT color index
/// @return ROOT color index
int ReadRun::rcolor(unsigned int i) {
	const int nclrs = 16;
	int rclrs[nclrs] = { 1, 2, 3, 4, 5, 6, 7, 13, 28, 30, 34, 38, 40, 31, 46, 49 };
	return rclrs[i - static_cast<int>(floor(i / nclrs)) * nclrs];
}

/// @brief Returns index of a certain event number (if data files are read in parallel threads)
/// @param eventnr Event number as stored in the data.
/// @return Corresponding event number in the internal data structure.
int ReadRun::GetEventIndex(int eventnr) {
	if (eventnr <= 0) eventnr = 1; // first event is 1
	if (eventnr > nevents) eventnr = nevents;
	return distance(eventnr_storage.begin(), find(eventnr_storage.begin(), eventnr_storage.end(), eventnr));
}

/// @brief  Match channel number (wavecatcher input channel) to channel index
/// @param channel_number Number of the channel as defined in the wavecatcher software
/// @return Corresponding index for this channel
int ReadRun::GetChannelIndex(int channel_number) {
	int channel_index = -1;
	for (int i = 0; i < static_cast<int>(active_channels.size()); i++) {
		if (active_channels[i] == channel_number) channel_index = i;
	}
	if (channel_index == -1) {
		cout << "\n\n\tERROR: channel " << channel_number << " does not exist in data. Will continue with first channel\n\n";
		channel_index = 0;
	}
	return channel_index;
}

/// @brief Helper to split canvas according to the number of channels to be plotted
/// @param c Canvas to be split
void ReadRun::SplitCanvas(TCanvas*& c) {
	if (plot_active_channels.empty()) c->Divide(TMath::Min(static_cast<double>(active_channels.size()), 4.), TMath::Max(TMath::Ceil(static_cast<double>(active_channels.size()) / 4.), 1.), 0, 0);
	else if (static_cast<int>(plot_active_channels.size()) > 1) c->Divide(TMath::Min(static_cast<double>(plot_active_channels.size()), 4.), TMath::Max(ceil(static_cast<double>(plot_active_channels.size()) / 4.), 1.), 0, 0);
}

/// @brief Simple linear interpolation for x
/// @param ym Y value for evaluation
/// @param x1 X1 
/// @param x2 X2
/// @param y1 Y1
/// @param y2 Y2
/// @return x value at "ym"
float ReadRun::LinearInterpolation(float ym, float x1, float x2, float y1, float y2) {
	return x1 + (ym - y1) * (x2 - x1) / (y2 - y1);
}

/// @brief Helper to perform convolution of two 1D arrays
/// 
/// Used for smoothing etc.
/// 
/// @param[in,out] result Array containing convolution result
/// @param first First array for convolution
/// @param second Second array for convolution 
/// @param size1 Size of first array (ideally size1 should be equal to size2)
/// @param size2 Size of second array
void ReadRun::Convolute(double*& result, double* first, double* second, int size1, int size2) {

	// uncomment to use sum instead of FFT
	// faster if size1<size2
	//for (int i = 0; i < size2; i++) {
	//	result[i] = 0.;
	//	for (int j = 0; j < TMath::Min(size1, i); j++) {
	//		result[i] += first[j] * second[i - j];
	//	}
	//}

	double* refirst = new double[size1];
	double* imfirst = new double[size1];
	double* resecond = new double[size1];
	double* imsecond = new double[size1];
	double* reres = new double[size1];
	double* imres = new double[size1];

	TVirtualFFT* fftfirst = TVirtualFFT::FFT(1, &size1, "R2C ES");
	fftfirst->SetPoints(first);
	fftfirst->Transform();
	fftfirst->GetPointsComplex(refirst, imfirst);
	delete fftfirst;

	TVirtualFFT* fftsecond = TVirtualFFT::FFT(1, &size1, "R2C ES");
	fftsecond->SetPoints(second);
	fftsecond->Transform();
	fftsecond->GetPointsComplex(resecond, imsecond);
	delete fftsecond;

	TComplex cofirst;
	TComplex cosecond;
	TComplex cores;

	for (int i = 0; i < size1; i++) {
		cofirst(refirst[i], imfirst[i]);
		cosecond(resecond[i], imsecond[i]);

		cores = cofirst * cosecond / static_cast<double>(size1);

		reres[i] = cores.Re();
		imres[i] = cores.Im();
	}

	//cout << "performing IFFT ... ";
	TVirtualFFT* fft_back = TVirtualFFT::FFT(1, &size1, "C2R ES");
	fft_back->SetPointsComplex(reres, imres);
	fft_back->Transform();
	fft_back->GetPoints(result);
	delete fft_back;
	delete[] imres; delete[] reres; delete[] refirst; delete[] imfirst; delete[] resecond; delete[] imsecond;
}

/// @brief Apply smoothing array of double with length nbins
/// 
/// Use with care. Method=2 is preferred. \n \n
///
/// Please note that if you want to use gaussian smoothing for data with a binning different from 0.3125 ns/bin you need to set the variable bin_size to the new bin size.
/// 
/// @param[in,out] ar Array to be smoothed.
/// @param nbins Number of bins of input.
/// @param sigma Number of bins before and after central bin for running average OR gauss sigma in ns for gauss kernel and convolution (see parameter bin_size).
/// @param method If 0 use running average (box kernel smoothing). Simple, very fast. \n 
/// If 1 use 5 sigma gaussian smoothing. This method is not central and will shift peaks. Very slow. \n
/// Else use 3 sigma gaussian kernel smoothing. Preferred method, fast.
/// @param bin_size Bin width of the array to smooth for gauss sigma. Default is .3125 for wavecatcher sampling rate. Set to 1 to change sigma unit to number of bins.
void ReadRun::SmoothArray(double*& ar, int nbins, double sigma, int method, double bin_size) {

	double* artmp = new double[nbins];
	for (int i = 0; i < nbins; i++) artmp[i] = ar[i];

	if (method == 0) {
		// calculate running average from -sigma until +sigma (sigma = number of bins)
		for (int i = 0; i < nbins; i++) {
			double mean1 = 0.;
			int nmn = 0;
			for (int k = -1 * static_cast<int>(floor(sigma)); k <= static_cast<int>(ceil(sigma)); k++) {
				if (i + k >= 0 && i + k < nbins) {
					mean1 += artmp[i + k];
					nmn++;
				}
			}
			if (nmn != 0.) {
				ar[i] = mean1 / static_cast<double>(nmn);
			}
		}
	}
	else if (method == 1) {
		// convolution with gauss clipped at +-5 sigma (very inefficient and slow)
		double* gauss = new double[nbins];

		double sum = 0.;
		double position = 0.;

		for (int i = 0; i < nbins; i++) {
			if (static_cast<double>(i) * bin_size < 5 * sigma) gauss[i] = TMath::Exp(-1. * TMath::Power((static_cast<double>(i) * bin_size - 5 * sigma), 2.) / (2. * sigma * sigma)) / (sigma * 2.506628);
			else gauss[i] = 0.;
			sum += gauss[i];
		}

		for (int i = 0; i < nbins; i++) {
			gauss[i] /= sum;
		}

		Convolute(ar, artmp, gauss, nbins, nbins);
		delete[] gauss;
	}
	else {
		// gauss kernel 3*sigma
		int nbins_3sigma = static_cast<int>(ceil(6. * sigma / bin_size));
		if (nbins_3sigma % 2 == 0) nbins_3sigma++;
		if (nbins_3sigma > 1) {
			double* gauss = new double[nbins_3sigma];
			double gauss_offset = floor(static_cast<double>(nbins_3sigma) / 2.) * bin_size;
			double denom = 2. * sigma * sigma;
			for (int i = 0; i < nbins_3sigma; i++) {
				gauss[i] = TMath::Exp(-1. * TMath::Power((static_cast<double>(i)) * bin_size - gauss_offset, 2.) / denom);
			}

			double res = 0;
			double norm = 0;
			for (int i = 0; i < nbins; i++) {
				res = 0.;
				norm = 0.;
				for (int j = max(0, nbins_3sigma / 2 - i); j < min(nbins - i + nbins_3sigma / 2, nbins_3sigma); j++) {
					res += gauss[j] * artmp[i + j - nbins_3sigma / 2];
					norm += gauss[j];
				}
				if (norm != 0.) ar[i] = res / norm;
			}
			delete[] gauss;
		}
	}
	delete[] artmp;
}

/// @brief Apply filter for array of double with length nbins
/// 
/// Experimental, can be used to highlight peaks and suppress long tails (suppresses low and high frequencies). Use Filter_test.ipynb to test parameters.
/// 
/// @param[in,out] ar Array to be filtered.
/// @param nbins Number of bins of input.
/// @param sigma1 First.
/// @param sigma2 Second.
/// @param factor Factor for negative part (<=1).
/// @param bin_size Bin width. Default is .3125. Set to 1 to get sigma in units of bins.
void ReadRun::FilterArray(double*& ar, int nbins, double sigma1, double sigma2, double factor, double bin_size) {

	double* artmp = new double[nbins];
	for (int i = 0; i < nbins; i++) artmp[i] = ar[i];

	// shifted difference of two gauss functions (~smoothed derivative)
	int nbins_2sigma = static_cast<int>(ceil((2. * sigma1 + 3. * sigma2) / bin_size));
	double* sdog = new double[nbins_2sigma];

	double denom1 = 2. * sigma1 * sigma1;
	double denom2 = 2. * sigma2 * sigma2;
	for (int i = 0; i < nbins_2sigma; i++) {
		sdog[i] = TMath::Exp(-1. * TMath::Power(static_cast<double>(i) * bin_size - 3. * sigma2, 2.) / denom1) - factor * TMath::Exp(-1. * TMath::Power(static_cast<double>(i) * bin_size - 2. * sigma2, 2.) / denom2);
	}

	double res = 0;
	double norm = 0;
	for (int i = 0; i < nbins; i++) {
		res = 0.;
		norm = 0.;
		for (int j = max(0, nbins_2sigma / 2 - i); j < min(nbins - i + nbins_2sigma / 2, nbins_2sigma); j++) {
			res += sdog[j] * artmp[i + j - nbins_2sigma / 2];
			if (sdog[j] > 0.) norm += sdog[j];
		}
		if (norm != 0.) ar[i] = res / norm;
	}
	delete[] sdog;
}