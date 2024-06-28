#include "Freiburg_DAQ.h"


/// @brief Routine to read converted data by the GANDALF DAQ System at Uni Freiburg. \n
/// Call ```python examples/Freiburg/read_Freiburg_data.py 0``` to convert and analyze test data. \n
/// 
/// Note that amc_hax.py needs to be located in ```examples/Freiburg/```
/// 
/// @param path File name of (currently only a single) data file
/// @param change_polarity Defaults to true, meaning the peaks are positive.
/// @param out_file_name File name of root file to store results in.
void Freiburg_DAQ::ReadFile(string path, bool change_polarity, string out_file_name) {

	printf("+++ saving analysis results in '%s' ...\n\n", out_file_name.c_str());
	root_out = TFile::Open(out_file_name.c_str(), "recreate");

	rundata = new TClonesArray("TH1F", maxNWF); //raw data will be stored here as TH1F
	rundata->BypassStreamer(kFALSE);			//kFALSE: potentially faster read & write
	TClonesArray& testrundata = *rundata;

	ifstream file(path.c_str(), ios::binary);

	if (!file.is_open()) cerr << "Error opening file" << endl;

	int32_t rows, cols;
	file.read(reinterpret_cast<char*>(&rows), sizeof(int32_t));
	file.read(reinterpret_cast<char*>(&cols), sizeof(int32_t));

	// harcoded parameters, might want to read directly from data header in the future
	SP = 1 / 0.880; // sampling rate is 880 MSa/s (0.880 Sa/ns)
	nchannels = 8; // currently limited to 8 channels
	binNumber = static_cast<int>(cols - 2); // number of samples read from file

	float ADC_factor = 2.2 / 4.096; // convert to mV
	if (change_polarity) ADC_factor *= -1;

	amplValuessum = new double* [nchannels]; //sum of all wf for each channel
	for (int i = 0; i < nchannels; i++) {
		amplValuessum[i] = new double[binNumber]();
	}
	maxSumBin = new int[nchannels];

	int wfcounter = 0;
	int last_event_number = -1;

	for (int i = 0; i < rows; ++i) {
		vector<int32_t> row(cols);
		file.read(reinterpret_cast<char*>(row.data()), cols * sizeof(int32_t));

		int event_number = static_cast<int>(row[0]);
		int channel = static_cast<int>(row[1] / 2);

		TString name(Form("ch%02d_%05d", channel, event_number));
		TString title(Form("ch%d, event %d;t [ns];U [mV]", channel, event_number));
		auto hCh = (TH1F*)testrundata.ConstructedAt(wfcounter);
		hCh->SetName(name.Data());
		hCh->SetTitle(title.Data());
		hCh->SetBins(binNumber, -0.5 * SP, (binNumber - 0.5) * SP);

		// Fill waveforms
		for (int s = 2; s < cols; ++s) {
			float val = static_cast<float>(row[s]) * ADC_factor;

			hCh->SetBinContent(s - 1, val);
			amplValuessum[channel][s - 2] += static_cast<double>(val);
		}
		wfcounter++;

		if (event_number == 0) active_channels.push_back(static_cast<int>(channel));

		if (last_event_number != event_number) {
			last_event_number = event_number;
			skip_event.push_back(false);
			eventnr_storage.push_back(event_number);
		}
	}
	file.close();

	// get bins where the sum spectrum has its maximum for runs with fixed trigger delay and fixed 
	// integration window relative to the max of the sum spectrum (not working for DC measurement)
	for (int ch = 0; ch < nchannels; ch++) {
		double max = 0.;
		for (int i = 0; i < binNumber; i++) {
			if (amplValuessum[ch][i] > max) {
				max = amplValuessum[ch][i];
				maxSumBin[ch] = i;
			}
		}
	}

	nevents = last_event_number;
	nwf = wfcounter;
}
