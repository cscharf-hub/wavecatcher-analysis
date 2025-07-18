/// Main class containing the file reader and most analysis functions
#ifndef _ReadRun
#define _ReadRun

// contact: christian.scharf@physik.hu-berlin.de

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <TObject.h>
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
#include <TH1F.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <THStack.h>
#include <THistPainter.h>
#include <TText.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TSpline.h>
#include <TPaveStats.h> 
#include <TError.h>      // root verbosity level
#include <TSystem.h>
#include <TROOT.h>
#include <TLatex.h>

//C, C++
#include <stdio.h>
#include <future>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <stdexcept>
//multithreading (avoids finding omp.h manually)
#ifndef OMP_FUNC
extern "C" {
    int omp_get_max_threads();
    int omp_get_thread_num();
}
#endif

#include "utils/FitFunctions.h"
#include "utils/Helpers.h"
#include "utils/Filters.h"

using namespace std;

class ReadRun {
private:
#pragma pack(1) // byte padding suppression for WaveCatcher data format
	// structs copied from
	// WaveCatcher binary -> root converter
	// by manu chauveau@cenbg.in2p3.fr
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
	};
#pragma pack() // byte padding suppression for WaveCatcher data format

protected:
	/// @brief Index for multiple executions of the same plotting function
	int PrintChargeSpectrum_cnt;
	/// @brief Index for multiple executions of the same plotting function
	int PlotChannelAverages_cnt;
	/// @brief Index for multiple executions of the same plotting function
	int PrintWFProjection_cnt;
	/// @brief Index for multiple executions of the same plotting function
	int PlotWFHeatmaps_cnt;
	
	/// @brief Primitive check to see if data has been loaded
    virtual void checkData(bool isBaselineCorrection = false) const {
		(void)isBaselineCorrection;
        if (eventnr_storage.empty()) {
            throw runtime_error(
				"Error: No data has been loaded yet!\n \
				Please call ReadFile() before calling any functions which manipulate data.\n \
				Aborting execution."
			);
        }
    }

public:
	/// @brief Stores waveforms
	vector<vector<float>> rundata;

	/// @brief Collects sums of all waveforms for each channel
	vector<vector<float>> amplValuessum;

	/// @brief Events will be stored here in the order they have been read
	vector<unsigned int> eventnr_storage;

	// plots amplValuessum
	void PlotChannelSums(bool = false, bool = false, double = 0., double = 0., int = 2);

	void PlotChannelAverages(bool = false);

	TH2F* WFHeatmapChannel(int, float = -20, float = 200, int = 880);
	void PlotWFHeatmaps(float = -20, float = 200, int = 880, string = "", float = 0, EColorPalette = kRainBow);
	
	// baseline correction (shifts all waveforms individually in y)
	void CorrectBaseline(float, float = -999);
	void CorrectBaseline_function(vector<float>&, float, float, int);

	void CorrectBaselineMinSlopeRMS(vector<float>, double = 0, int = 2);
	void CorrectBaselineMinSlopeRMS(vector<float>, double, int, int);
	void CorrectBaselineMinSlopeRMS(int = 100, bool = false, double = 0.5, int = 0, int = 0, int = 2);

	void CorrectBaselineMin(vector<float>, double = 0, int = 2);
	void CorrectBaselineMin(vector<float>, double, int, int);
	void CorrectBaselineMin(int = 100, double = 0.5, int = 0, int = 0, int = 2);

	// functions to check baseline correction results
	TH1F* WFProjectionChannel(int, int = 0, int = 1024, float = -50, float = 50, int = 200);
	void PrintWFProjection(float = 0, float = 320, float = -50, float = 50, int = 200);

	TH1F* BaselineCorrectionResults(int, int, float = -5, float = 5, int = 200);
	void PrintBaselineCorrectionResults(float = -5, float = 5, int = 200);

	// get timing of peaks
	void GetTimingCFD(float = .3, float = 100, float = 140, double = 0., bool = true, int = 2, bool = false, bool = false);
	void SkipEventsTimeDiffCut(int, int, double, double, bool = false);

	// average all waveforms to simplify peak ID
	void SmoothAll(double, int);
	void SmoothAll(double, string = "Gaus");
	void FilterAll(double = .3, double = .9, double = .2);
	void ShiftAllToAverageCF();

	// functions for charge spectrum
	array<int, 3> GetIntWindow(TH1F*, float, float, float, float, int = 0);
	array<int, 3> GetIntWindow(const vector<float>&, float, float, float, float, int = 0);
	array<int, 3> GetIntWindow(const vector<float>&, int, int, int, int, int = 0);
	float GetPeakIntegral(TH1F*, float, float, float, float, int = 0);
	float GetPeakIntegral(const vector<float>&, float, float, float, float, int = 0);
	void PrintChargeSpectrumWF(float, float, float = 0, float = 300, int = 1, float = 0., float = 0., float = 0., float = 0.);
	TH1F* ChargeSpectrum(int, float, float, float = 0, float = 300, float = -50, float = 600, int = 750);
	void PrintChargeSpectrum(float, float, float = 0, float = 300, float = -50, float = 600, int = 750, float = 0., float = 0., int = 99, int = 0, bool = false);
	/// @brief Starting values of the fit parameters for PrintChargeSpectrum()
	vector<float> PrintChargeSpectrum_pars;


	float* ChargeList(int, float = 20, float = 80, float = 0, float = 300, bool = 1);
	void SaveChargeLists(float = 20, float = 80, float = 0, float = 300, bool = 1);
	void ChargeCorrelation(float, float, float, float, float, float, int, int, int, bool = false);

	// SiPM specific
	void PrintDCR(float = 15, float = 85, float = 0, float = 300, double = 3); // based on PMT::PrintChargeSpectrumPMTthreshold()

	// functions for time distribution
	TH1F* TimeDist(int, float = 0, float = 300, float = 0, float = 300, int = 100, int = 0, float = .3);
	void PrintTimeDist(float = 0, float = 300, float = 0, float = 300, int = 100, int = 0, float = .3);
	TGraph2D* MaxDist(int, float = 0, float = 300);
	void PrintMaxDist(float = 0, float = 300);

	TH1F* His_GetTimingCFD(int, float, float, int = -999);
	void Print_GetTimingCFD(float = 100, float = 140, int = 0, int = -999, string = "S", bool = true);
	TH1F* His_GetTimingCFD_diff(vector<int>, vector<int>, float, float, int = -999);
	void Print_GetTimingCFD_diff(vector<int>, vector<int>, float = 100, float = 140, int = 0, int = -999, float = -999, float = -999, string = "RS", bool = true);

	
	// helper functions
	TH1F* Getwf(int);							// waveform number
	virtual TH1F* Getwf(int, int, int = 1);		// channel, eventnr, color
	template<typename T>
	T* getx(double shift = 0.);					// x values
	double* gety(int);							// y values for waveform index
	double* gety(int, int);						// y values for waveform(channel, event)

	static pair<float, bool> LinearInterpolation(float, float, float, float, float, bool = false); // linear interpolation
	
	virtual int GetWaveformIndex(int, int);	// get index of waveform from event and channel
	int GetEventIndex(unsigned int);		// get index of a triggered event (finds the correct event if files are not read sequentially)
	int GetChannelIndex(int);				// get index of a certain channel
	virtual int GetCurrentChannel(int);		// get index of channel for a certain waveform
	virtual int GetCurrentEvent(int);		// get index of event for a certain waveform
	
	bool PlotChannel(int);					// check if channel should be plotted
	
	vector<bool> PolarityMap(bool, int); 	// check if channel should be inverted

	int CheckBoundsX(int);					// Check if index exists in time of waveforms
	int TimeToIndex(float);					// Convert time to the bin number of the waveform
	float IndexToTime(int);					// Convert the bin number of the waveform to the time

	/// @brief Constructor of the class
	/// @param last_bin_file Number of last .bin file to be read in.\n
	/// Set it to =>1 in order to constrain the number of .bin files to be read from the target folder.\n
	/// Intended for quick tests on a fraction of the full dataset or for batch reading if combined with min_no_of_bin_files_to_read.
	/// @param first_bin_file Number first of .bin file to be read in. \n
	/// Can be used to batch read large datasets in chunks of files.
	ReadRun(int last_bin_file = 0, int first_bin_file = 0);
	
	virtual void ReadFile(string, bool = false, int = 9, string = "out.root", bool = false, long long = -1);

	virtual ~ReadRun();

	/// @brief Path to data
	/// 
	/// Can be used to save analysis results in the data folder
	string data_path;

	/// @brief Number of last .bin file to be read in. 
	///
	/// Can be used to test analysis on a small sample of the data or for batch reading if combined with FirstBinFileToRead.
	int LastBinFileToRead;
	/// @brief Number first of .bin file to be read in. 
	///
	/// Can be used to batch read large datasets in chunks of files.
	int FirstBinFileToRead;


	/// @brief Can be used to discard the original event numbering of the data
	/// 
	/// Set to true **BEFORE** calling ReadFile() if you want to read several runs at once. 
	/// The events will be numbered in the order they are read in. 
	/// The original event numbers of the different runs will be lost.
	/// CAUTION: All .bin files of the different runs need to contain the same number of channels.
	bool discard_original_eventnr = false;

	/// @brief Do analysis only for limited range of channels to reduce memory usage
	/// 
	/// For large datasets with many channels and many events \n
	/// Only read and analyze channels from ReadRun::start_read_at_channel to ReadRun::end_read_at_channel. \n
	/// The recorded channel with the lowest wavecatcher channel number is 0 (e.g. recorded channels 3 and 4,
	///  so start would be 0 and end 1). \n 
	/// If set to -1 (default) all channels will be read in one go. \n
	/// Else channels from "start_read_at_channel" to "end_read_at_channel" will be read. \n 
	/// If "end_read_at_channel" is not defined will only read channel specified in "start_read_at_channel".
	int start_read_at_channel = -1;
	/// @brief See ReadRun::start_read_at_channel
	int end_read_at_channel = -1;

	/// @brief Number of triggered events in data
	int nevents = 0;
	/// @brief Number of active channels in data
	int nchannels;
	/// @brief Total number of waveforms read from data: number of active channels x number of events
	int nwf;
	/// @brief Maximum possible number of waveforms in data: number of active channels x number of events
	///
	/// Default is 10 million waveforms. Increase this number for very large runs.
	int maxNWF = 1e7;

	/// @brief Sampling: ns per bin of data, sampling rate 3.2 GS/s -> 0.3125 ns
	float SP = .3125;
	/// @brief 1/SP
	float SP_inv = 1 / .3125;
	/// @brief DAQ conversion factor for wavecatcher output to mV
	/// 
	/// From https://owncloud.lal.in2p3.fr/public.php?service=files&t=56e4a2c53a991cb08f73d03f1ce58ba2 
	float DAQ_factor = 250. / 4096.;
	/// @brief Number of bins (usually 1024 samples per waveform).
	int binNumber;
	/// @brief Wavecatcher hardware max. number of channels (reduce/increase if not using the 64 channel crate)
	int nChannelsWC = 64;

	/// @brief Stores bin numbers where the sum of waveforms have their maximum
	/// 
	/// Can be used for fixed integration window relative to maximum of the sum of all waveforms per channel (ReadRun::amplValuessum)
	vector<int> maxSumBin;

	/// @brief Stores the numbers of the active channels
	vector<int> active_channels;
	/// @brief Stores the numbers of the active channels which should be plotted
	///
	/// You can select the channels you want with plot_active_channels.push_back(channel_to_plot); to add them to the list. \n 
	/// If undefined all channels will be plotted.
	vector<int> plot_active_channels;
	
	/// @brief Stores the channel number where the polarity should be inverted.
	/// Example use to switch polarity for channel 0 and channel 14: ```mymeas.switch_polarity_for_channels = {0, 14};```.
	/// Needs to be called before ReadFile().
	vector<int> switch_polarity_for_channels;

	/// @brief Stores the fit results of PrintChargeSpectrum() for all channels and all function calls in ascending order 
	vector<TFitResultPtr> fit_results;

	/// @brief Stores the mean integral/lightyield from PrintChargeSpectrum() for all channels
	vector<float> mean_integral;

	/// @brief Stores the event numbers which should be skipped in the analysis
	/// 
	/// To identify events to be filtered use functions IntegralFilter(), SkipEventsPerChannel(), and SkipEventsTimeDiffCut().
	vector<bool> skip_event;
	virtual int Nevents_good(int = 0);

	void SkipEventsPerChannel(vector<float>, float = 0, float = 0, bool = false);  // in case you want to have indiviual thresholds in individual channels
	void IntegralFilter(vector<float>, vector<bool>, float, float, float = 50, float = 250, bool = false, bool = false); // Same as SkipEventsPerChannel() but filtering all events with integrals <(>) threshold
	void PrintSkippedEvents();
	void UnskipAll();
	virtual bool SkipEvent(int, int = -1);

	/// @brief Stores baseline correction results for CorrectBaseline() and related functions
	vector<vector<float>> baseline_correction_result;

	/// @brief Matrix to store timing of peaks from GetTimingCFD()
	/// 
	/// First index is the index of the waveform. \n
	/// Second index is:\n
	/// 0: CFD bin number
	/// 1: CFD time
	/// 2: Y-value of the maximum
	/// 3: Bin number of the maximum
	/// 4: Constant fraction (maximum * cf_r)
	/// 5: Start of search window
	/// 6: End of search window
	vector<vector<float>> timing_results;
	/// @brief Stores the fit results of Print_GetTimingCFD() for all channels
	vector<TFitResultPtr> timing_fit_results;

	/// @brief Stores results of analysis
	TFile* root_out;

	//other controls 

	/// @brief Set true for baseline correction during data reading.
	/// Needs to be called before ReadFile().
	bool Using_BaselineCorrection_in_file_loop = false;
	/// @brief Start of time window for baseline correction when using ReadRun::Using_BaselineCorrection_in_file_loop
	float tCutg;
	/// @brief End of time window for baseline correction when using ReadRun::Using_BaselineCorrection_in_file_loop
	float tCutEndg;

	/// @brief Calibration values to normalize charge spectrum to number of photoelectrons
	/// Chennels must be ordered as in plot_active_channels.\n
	/// The first entry must be the gain and the second the position of the pedestal.
	vector<vector<float>> PrintChargeSpectrum_cal = vector(nChannelsWC, vector<float>(2, 1));
};
#endif