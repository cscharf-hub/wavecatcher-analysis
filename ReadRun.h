#ifndef _ReadRun
#define _ReadRun

// contact: christian.scharf@cern.ch
// some includes are probably not needed anymore

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TClonesArray.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TComplex.h>
#include <TVirtualFFT.h>
#include <TLine.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TMultiGraph.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TF1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <THStack.h>
#include <THistPainter.h>
#include <TText.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TSpectrum.h>   // peakfinder
#include <TPolyMarker.h> // peakfinder
#include <TError.h>      // root verbosity level
#include <TSystem.h>     // root verbosity level
#include <TLatex.h>      // root verbosity level

//#include <sys/resource.h>
//C, C++
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <numeric>
#include <tuple>
#include <map>

using namespace std;

class ReadRun/* : public TObject*/ {
private:

	//int maxnchannels;				// number of channels (32)
	TClonesArray* rundata;		// data array
	//TClonesArray* blh;

	double** amplValuessum;		// collects sums of all waveforms for each channel

	vector<unsigned int> eventnr_storage;	//  DORAMAS: The events will be stored here in the order they have been read

	// for baseline correction before 
	bool Using_BaselineCorrection_in_file_loop = false;
	float tCutg;
	float tCutEndg;

#pragma pack(1) // padding suppression
// struct copied from
// WaveCatcher binary -> root converter
// by manu chauveau@cenbg.in2p3.fr
// see https://owncloud.lal.in2p3.fr/public.php?service=files&t=56e4a2c53a991cb08f73d03f1ce58ba2

	struct event_data
	{
		int EventNumber;
		double EpochTime;
		unsigned int Year;
		unsigned int Month;
		unsigned int Day;
		unsigned int Hour;
		unsigned int Minute;
		unsigned int Second;
		unsigned int Millisecond;
		unsigned long long int TDCsamIndex;
		int nchannelstored;
	};

	struct channel_data_without_measurement
	{
		int channel;
		int EventIDsamIndex;
		int FirstCellToPlotsamIndex;
		short waveform[1024];
	};

	struct channel_data_with_measurement
	{
		int channel;
		int EventIDsamIndex;
		int FirstCellToPlotsamIndex;
		float	MeasuredBaseline;
		float	AmplitudeValue;
		float	ComputedCharge;
		float	RiseTimeInstant;
		float	FallTimeInstant;
		float	RawTriggerRate;
		short waveform[1024];
	};
#pragma pack() // padding suppression

public:

	// plots amplValuessum
	void PlotChannelSums(bool = true);

	// baseline correction (shifts all waveforms individually)
	void CorrectBaseline(float, float = -999);
	void CorrectBaseline_function(TH1F*, float, float, int);

	void CorrectBaselineMinSlopeRMS(int = 100, bool = false, double = 10, int = 0, int = 0, bool = false, bool = false, int = 8);

	void CorrectBaselineMin(int = 100, bool = false, double = 10, int = 0, int = 0, bool = false, int = 8);

	void FractionEventsAboveThreshold(float = 4, bool = true, bool = true, double = 0., double = 0., bool = false);

	// average all waveforms to simplify peak ID
	void SmoothAll(double = 5, bool = false);
	void DerivativeAll();

	// functions for charge spectrum
	int* GetIntWindow(TH1F*, float, float, float, float, int);
	void PrintChargeSpectrumWF(float, float, float = 0, float = 300, int = 1, float = 0., float = 0.);
	TH1F* ChargeSpectrum(int, float, float, float = 0, float = 300, float = -50, float = 600, int = 750, string = "width");
	void PrintChargeSpectrum(float, float, float = 0, float = 300, float = -50, float = 600, int = 750, float = 0., float = 0., int = 8, int = 0);
	vector<float> PrintChargeSpectrum_pars;
	void PrintChargeSpectrumPMT(float, float, float = 0, float = 300, float = -50, float = 600, int = 750);
	vector<float> PrintChargeSpectrumPMT_pars;
	void PrintChargeSpectrumPMTthreshold(float, float, float = 0, float = 300, float = -50, float = 600, int = 750, double = 4);

	// functions for time distribution
	TH1F* TimeDist(int, float = 0, float = 300, float = 0, float = 300, int = 100, int = 0);
	void PrintTimeDist(float = 0, float = 300, float = 0, float = 300, int = 100, int = 0);
	TGraph2D* MaxDist(int, float = 0, float = 300);
	void PrintMaxDist(float = 0, float = 300);

	// print FFT
	void PrintFFTWF(int = 1, float = 0., float = 0., int = 1);

	// helper functions
	stringstream list_files(const char*, const char*);	// find data files
	TH1F* Getwf(int, int, int = 1);						// channel, eventnr, color
	double* getx(double = 0.);							// x values
	double* gety(int, int);								// y values for waveform(ch, event)
	double* gety(TH1F*);								// y values for histogram

	float LinearInterpolation(float, float, float, float, float); // linear interpolation

	int GetEventIndex(int);										// get index of a triggered event (finds the correct event if files are not read sequentially)
	void SplitCanvas(TCanvas*&);								// split canvas into pads to display all active channels on one canvas
	void Convolute(double*&, double*, double*, int, int);		// convolution for filtering waveforms
	void SmoothArray(double*&, int, double = 1., bool = false);	// filtering

	ReadRun(double = 0, int = 1); // Constructor of the class with arguments to filter noise events in the cosmics setup. Default values do nothing

	void ReadFile(string, bool = false, int = 9, string = "out.root", bool = false, bool = false); // file name, bool whether or not to change sign of PMT channels (channel number>8), bool whether to save ALL waveforms to root file (only advisable for runs with small number of events)

	virtual ~ReadRun();

	string data_path;			// path to data. Can be used to save analysis results in the data folder
	//int nbinsdata;
	int nevents;				// number of triggered events
	int nchannels;
	int nwf;					// number of waveforms (nchannels*nacquisitions)

	float SP;					// ns per bin
	float pe;					// mV*ns ????
	double coef;				// ?????
	int binNumber;				// 1024 samples per waveform

	int* maxSumBin;				// For fixed integration window (triggered acquisition)

	vector<int> active_channels; // stores the numbers of the active channels
	vector<int> plot_active_channels; // stores the numbers of the active channels which should be plotted

	vector<TFitResultPtr> fit_results; // stores the fit results of all channels and all function calls in ascending order for all different PrintChargeSpectrum functions

	vector<bool> skip_event; // stores the events which should be skipped in the analysis
	double skip_event_threshold; // threshold (usually 4 mV) for PMT signal (hardcoded channel >8) to skip events where PMTs pick up radio frequency noise (NO BASELINE CORRECTION!)
	int skip_event_threshold_nch; // define how many PMT channels need to be above threshold to discard event (RF pick up should be seen by alls PMTs)
	void SkipEventsPerChannel(vector<double>, bool = false);  // in case you want to have indiviual thresholds in individual channels

	vector<vector<float>> baseline_correction_result; // store baseline values

	TFile* root_out;
	TFile* root_out_wf;

	ClassDef(ReadRun, 1)
};
#endif

class Fitf {
public:
	// as used by jan (missing after-pulses and dark counts)

	double operator() (double* x, double* p) {
		//0 - N0: Normalization (~Number of events)
		//1 - mu: for generalized poisson distribution
		//2 - lambda: Borel-branching parameter for prompt crosstalk probability 1-exp(-lambda)

		//3,4 -sigma0, sigma1
		//5 - G: gain
		//6 - B: Pedestal
		double sum = 0;
		for (int kint = 0; kint <= 50; kint++) {
			double mu = p[1];
			double lambda = p[2];

			double sigma0 = p[3];
			double sigma1 = p[4];

			double G = p[5];
			double B = p[6];

			double k = static_cast<double>(kint);
			double sigmaK = sqrt(sigma0 * sigma0 + k * sigma1 * sigma1);
			//generalized poisson envelope
			double gp = mu * TMath::Power((mu + k * lambda), k - 1) * TMath::Exp(-(mu + k * lambda)) / TMath::Factorial(kint);

			sum += p[0] * gp * (1. / sqrt(2. * TMath::Pi()) / sigmaK) * TMath::Exp(-TMath::Power(((x[0] - (k * G + B)) / sqrt(2) / sigmaK), 2));
		}
		return sum;
	};
};

class Fitf_full {
public:
	// as used by Robert Klanner
	// still missing dark counts in integration window (3.3 in paper)
	// please check for possible bugs

	double operator() (double* x, double* p) {
		//0 - N0: Normalization (~Number of events)
		//1 - mu:  mean number of photons initiating a Geiger discharg
		//2 - lambda: Borel-branching parameter for prompt crosstalk probability 1-exp(-lambda)

		//3,4 -sigma0, sigma1
		//5 - G: gain
		//6 - B: Pedestal

		//7 - alpha: after-pulsing probability
		//8 - beta: the inverse of the exponential slope of the after-pulse	PH distribution

		double sum = 0;
		for (int kint = 0; kint <= 10; kint++) {
			double mu = p[1];
			double lambda = p[2];

			double sigma0 = p[3];
			double sigma1 = p[4];

			double G = p[5];
			double B = p[6];

			double alpha = p[7];
			double beta = p[8];

			//numbe of fired cells
			double k = static_cast<double>(kint);
			//pulse width
			double sigmaK = sqrt(sigma0 * sigma0 + k * sigma1 * sigma1);
			double gausnormsigmak = 1 / (sqrt(2. * TMath::Pi()) * sigmaK);
			//generalized poisson envelope
			double gp = mu * TMath::Power((mu + k * lambda), k - 1) * TMath::Exp(-(mu + k * lambda)) / TMath::Factorial(kint);
			//gauss peak
			double gauss = gausnormsigmak * TMath::Exp(-1. * TMath::Power(x[0] - (k * G + B), 2.) / (sigmaK * sigmaK * 2.));


			if (kint == 0) {
				sum += p[0] * gp * gauss;
			}
			else {
				//after-pulse + delayed cross talk
				double bk0 = TMath::Power(1. - alpha, k);
				double bk1 = TMath::Factorial(kint) / TMath::Factorial(kint - 1) * alpha * TMath::Power(1. - alpha, k - 1);
				double pk1 = TMath::Exp(-1. * (x[0] - (k * G + B)) / beta) * gausnormsigmak / beta * sigmaK * sqrt(TMath::Pi() / 2) * (TMath::Erf((x[0] - (k * G + B)) / (sqrt(2) * sigmaK)) + 1.);


				if (kint == 1) sum += p[0] * gp * (bk0 * gauss + bk1 * pk1);
				else {
					double api2k = 0.;
					for (int ii = 2; ii < kint; ii++) {
						double iid = static_cast<double>(ii);
						double bkialpha = TMath::Factorial(kint) / (TMath::Factorial(ii) * TMath::Factorial(kint - ii)) * TMath::Power(alpha, iid) * TMath::Power((1. - alpha), (k - iid));
						double dpkidph = TMath::Power((x[0] - (k * G + B)), (iid - 1.)) / (TMath::Factorial(ii - 1) * TMath::Power(beta, iid)) * TMath::Exp(-1. * (x[0] - (k * G + B)) / beta);
						api2k += bkialpha * dpkidph;
					}
					sum += p[0] * gp * (bk0 * gauss + bk1 * pk1 + api2k);
				}
			}
		}
		return sum;
	};
};

class Fitf_biased {
public:
	// For fitting charge spectrum with biased pedestal

	double operator() (double* x, double* p) {
		//0 - N0: Normalization (~Number of events)
		//1 - mu: for generalized poisson distribution
		//2 - lambda: Borel-branching parameter for prompt crosstalk probability 1-exp(-lambda)

		//3,4 -sigma0, sigma1
		//5 - G: gain
		//6 - B: Virtual pedestal shift of pe peaks
		//7 - Pedestal scaling for biased pedestal
		//8 - Position of biased pedestal
		double sum = 0;
		for (int kint = 0; kint <= 50; kint++) {
			double mu = p[1];
			double lambda = p[2];

			double sigma0 = p[3];
			double sigma1 = p[4];

			double G = p[5];
			double B = p[6];

			double a_ped = p[7];
			double x_ped = p[8];

			double k = static_cast<double>(kint);
			double sigmaK = sqrt(sigma0 * sigma0 + k * sigma1 * sigma1);
			//generalized poisson envelope
			double gp = mu * TMath::Power((mu + k * lambda), k - 1) * TMath::Exp(-(mu + k * lambda)) / TMath::Factorial(kint);

			if (kint == 0) {
				sum += a_ped * p[0] * gp * (1. / sqrt(2. * TMath::Pi()) / sigmaK) * TMath::Exp(-TMath::Power(((x[0] - x_ped) / sqrt(2) / sigmaK), 2));
			}
			else {
				sum += p[0] * gp * (1. / sqrt(2. * TMath::Pi()) / sigmaK) * TMath::Exp(-TMath::Power(((x[0] - (k * G + B)) / sqrt(2) / sigmaK), 2));
			}
		}
		return sum;
	};
};

class Fitf_PMT {
public:
	// Gauss-Poisson
	// https://doi.org/10.1016/0168-9002(94)90183-X 

	double operator() (double* x, double* p) {
		//0 - A:		normalization to number of events in fit region
		//1 - w:		probability for type II BG
		//2 - alpha:	coefficient of exponential decrease of typ II BG
		//3 - sigma0:	sigma of pedestal
		//4 - Q0:		position of pedestal
		//5 - mu:		mean number of PE
		//6 - sigma1:	width of 1 PE peak
		//7 - Q1:		position of 1 PE peak

		double pmt_charge_spectrum = 0.;

		for (int kint = 0; kint <= 50; kint++) {
			double k = static_cast<double>(kint);

			double poiss = TMath::Power(p[5], k) * TMath::Exp(-p[5]) / TMath::Factorial(kint);

			double normgn_term = 1.;
			if (kint != 0) normgn_term /= (p[6] * TMath::Sqrt(2. * TMath::Pi() * k));
			double gn_term = (1. - p[1]) * normgn_term;
			if (kint != 0) gn_term *= TMath::Exp(-1. * TMath::Power(x[0] - k * p[7] - p[4], 2.) / (2. * k * p[6] * p[6]));
			else if (x[0] != p[4]) gn_term *= 0.; // delta function

			double Qn = p[4] + k * p[7];
			double sigman = TMath::Sqrt(p[3] * p[3] + k * p[6] * p[6]);
			double ignxe_term_exp = p[2] / 2. * TMath::Exp(-1. * p[2] * (x[0] - Qn - p[2] * sigman * sigman));

			double sgn_arg = x[0] - Qn - sigman * sigman * p[2];
			double arg_sgn = 1.;
			if (sgn_arg < 0.) arg_sgn = -1.;
			double ignxe_term_erf = TMath::Erf(fabs(p[4] - Qn - sigman * sigman * p[2]) / (sigman * TMath::Sqrt(2.))) + arg_sgn * TMath::Erf(fabs(sgn_arg) / (sigman * TMath::Sqrt(2.)));

			pmt_charge_spectrum += p[0] * poiss * (gn_term + p[1] * ignxe_term_exp * ignxe_term_erf);
		}
		if (pmt_charge_spectrum < 0.) pmt_charge_spectrum = 0.;
		return pmt_charge_spectrum;
	};
};

class Fitf_PMT_pedestal {
public:
	// Gauss-Poisson
	// https://doi.org/10.1016/0168-9002(94)90183-X 

	double operator() (double* x, double* p) {
		//0 - A:		normalization to number of events in fit region
		//1 - w:		probability for type II BG
		//2 - alpha:	coefficient of exponential decrease of typ II BG
		//3 - sigma0:	sigma of pedestal
		//4 - Q0:		position of pedestal
		//5 - mu:		mean number of PE
		//6 - sigma1:	width of 1 PE peak
		//7 - Q1:		position of 1 PE peak
		//8 - norm0:	norm of 0 PE peak

		double pmt_charge_spectrum = 0.;

		for (int kint = 0; kint <= 50; kint++) {
			double k = static_cast<double>(kint);

			double poiss = TMath::Power(p[5], k) * TMath::Exp(-p[5]) / TMath::Factorial(kint);

			double normgn_term = 1.;
			if (kint != 0) normgn_term /= (p[6] * TMath::Sqrt(2. * TMath::Pi() * k));
			else normgn_term *= p[8] / (p[3] * TMath::Sqrt(2. * TMath::Pi()));
			double gn_term = (1. - p[1]) * normgn_term;
			if (kint != 0) gn_term *= TMath::Exp(-1. * TMath::Power(x[0] - k * p[7] - p[4], 2.) / (2. * k * p[6] * p[6]));
			else  gn_term *= TMath::Exp(-1. * TMath::Power(x[0] - p[4], 2.) / (2. * p[3] * p[3]));

			double Qn = p[4] + k * p[7];
			double sigman = TMath::Sqrt(p[3] * p[3] + k * p[6] * p[6]);
			double ignxe_term_exp = p[2] / 2. * TMath::Exp(-1. * p[2] * (x[0] - Qn - p[2] * sigman * sigman));

			double sgn_arg = x[0] - Qn - sigman * sigman * p[2];
			double arg_sgn = 1.;
			if (sgn_arg < 0.) arg_sgn = -1.;
			double ignxe_term_erf = TMath::Erf(fabs(p[4] - Qn - sigman * sigman * p[2]) / (sigman * TMath::Sqrt(2.))) + arg_sgn * TMath::Erf(fabs(sgn_arg) / (sigman * TMath::Sqrt(2.)));

			pmt_charge_spectrum += p[0] * poiss * (gn_term + p[1] * ignxe_term_exp * ignxe_term_erf);
		}
		if (pmt_charge_spectrum < 0.) pmt_charge_spectrum = 0.;
		return pmt_charge_spectrum;
	};
};

class Fitf_PMT_ideal {
public:
	// Gauss-Poisson
	// https://doi.org/10.1016/0168-9002(94)90183-X 

	double operator() (double* x, double* p) {
		//0 - A_s:		norm. of PE spectrum
		//1 - mu:		mean number of PE
		//2 - sigma1:	width of 1st PE peak
		//3 - Q1:		position (gain*e) of 1st peak
		double pmt_charge_spectrum = 0.;

		for (int kint = 1; kint <= 100; kint++) {
			double k = static_cast<double>(kint);

			double poiss = TMath::Power(p[1], k) * TMath::Exp(-1. * p[1]) / TMath::Factorial(kint);

			double norm = 1. / (p[2] * TMath::Sqrt(2. * TMath::Pi() * k));

			double gauss = TMath::Exp(-1. * TMath::Power(x[0] - k * p[3], 2) / (2. * k * p[2] * p[2]));

			pmt_charge_spectrum += p[0] * poiss * norm * gauss;
		}
		return pmt_charge_spectrum;
	};
};

class Fitf_langaus {
public:
	// Landau-Gauss-convolution
	// copied from https://root.cern.ch/doc/master/langaus_8C.html 

	double operator() (double* x, double* par) {
		//Fit parameters:
		//par[0]=Width (scale) parameter of Landau density
		//par[1]=Most Probable (MP, location) parameter of Landau density
		//par[2]=Total area (integral -inf to inf, normalization constant)
		//par[3]=Width (sigma) of convoluted Gaussian function
		//
		//In the Landau distribution (represented by the CERNLIB approximation),
		//the maximum is located at x=-0.22278298 with the location parameter=0.
		//This shift is corrected within this function, so that the actual
		//maximum is identical to the MP parameter.

		// Numeric constants
		Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
		Double_t mpshift = -0.22278298;       // Landau maximum location

		// Control constants
		Double_t np = 100.0;      // number of convolution steps
		Double_t sc = 5.0;      // convolution extends to +-sc Gaussian sigmas

		// Variables
		Double_t xx;
		Double_t mpc;
		Double_t fland;
		Double_t sum = 0.0;
		Double_t xlow, xupp;
		Double_t step;
		Double_t i;


		// MP shift correction
		mpc = par[1] - mpshift * par[0];

		// Range of convolution integral
		xlow = x[0] - sc * par[3];
		xupp = x[0] + sc * par[3];

		step = (xupp - xlow) / np;

		// Convolution integral of Landau and Gaussian by sum
		for (i = 1.0; i <= np / 2; i++) {
			xx = xlow + (i - .5) * step;
			fland = TMath::Landau(xx, mpc, par[0]) / par[0];
			sum += fland * TMath::Gaus(x[0], xx, par[3]);

			xx = xupp - (i - .5) * step;
			fland = TMath::Landau(xx, mpc, par[0]) / par[0];
			sum += fland * TMath::Gaus(x[0], xx, par[3]);
		}

		return (par[2] * step * sum * invsq2pi / par[3]);
	};
};