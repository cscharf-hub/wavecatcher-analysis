/// \mainpage Introduction
/// This page serves as documentation of the waveform analysis framework ```wavecatcher-analysis``` for WaveCatcher setups 
/// in the experimental elementary particle physics group at the Institute of Physics at Humboldt University of Berlin. \n \n \n
/// 
/// You can find the documentation of the functions and variables at <a href="classReadRun.html"> ReadRun Class Reference</a>
/// 
/// Installation instructions can be found at <https://cscharf-hub.github.io/wavecatcher-analysis/>
/// 
/// The source code is available at <https://github.com/cscharf-hub/wavecatcher-analysis>
/// 
/// Development and maintenance: \n
/// Christian Scharf \n 
/// Contributors: \n 
/// Doramas Jimeno Sanchez \n
/// Alessia Brignoli \n
/// Ben Skodda \n 
/// Christophe Mullesch \n
/// Constantin Eckardt \n
/// Alexander Vagts \n

#include "ReadRun.h"

/// @brief Constructor of the class
ReadRun::ReadRun(int last_bin_file, int first_bin_file) {

	cout << "\ninitializing ..." << endl;

	if (last_bin_file > 0) {
		cout << "will read " << last_bin_file - first_bin_file << " .bin files from file number " 
			<< first_bin_file << " to file number " << last_bin_file << endl;
	}
	if (first_bin_file > 0) discard_original_eventnr = true;

	ROOT::EnableImplicitMT();
	TH1::AddDirectory(kFALSE);
	// init counters
	nwf = 0;
	PrintChargeSpectrum_cnt = 0;
	PlotChannelAverages_cnt = 0;
	PrintWFProjection_cnt = 0;
	PlotWFHeatmaps_cnt = 0;
	LastBinFileToRead = last_bin_file;
	FirstBinFileToRead = first_bin_file;

	root_out = new TFile();	// init results file
}

/// @brief Routine to read files created by the wavecatcher.
/// 
/// Reads the data and can already do simple data manipulation during reading: \n
///   - Change polarity of range of channels with parameters (explained below). \n 
///   - Shift all waveforms to a constant fraction such they all start at the same time (-> see ReadRun::Shift_WFs_in_file_loop ) \n 
///   - Simple baseline correction by calling CorrectBaseline() before calling ReadFile(). \n 
/// Stores the data as vectors in ```vector<vector<float>> rundata```. \n \n 
/// 
/// Reader modified from WaveCatcher binary -> root converter by manu chauveau@cenbg.in2p3.fr \n 
/// 
/// @param path Path to the data. All files in this folder containing ```.bin``` in the file name will be read in.
/// @param change_polarity Set ```true``` to change polarity (sign) of certain channels (see ```change_sign_from_to_ch_num``` below). 
/// You can also define a list of channels where the polarity should be switched with switch_polarity_for_channels.
/// @param change_sign_from_to_ch_num All channels \f$ \geq \f$ ```change_sign_from_to_ch_num``` 
/// will be inverted if ```change_polarity``` is ```true```. \n 
/// If negative number all channels \f$ \leq \f$ ```abs(change_sign_from_to_ch_num)``` will be inverted if ```change_polarity``` is ```true```.
/// @param out_file_name Name of the ```.root``` file which stores the results, e. g. ```results.root```.
/// @param debug Set ```true``` to increase the verbosity.
/// @param max_nevents_to_read Maximum number of events to be read per file. For quick testing of analysis on subset of the data.
void ReadRun::ReadFile(string path, bool change_polarity, int change_sign_from_to_ch_num, string out_file_name, bool debug, long long max_nevents_to_read) {
	if (max_nevents_to_read <= 0) max_nevents_to_read = static_cast<long long>(1e9);
	rundata.reserve(1'000'000);						// reserve space for 1M waveforms
	if (path.back() != '/') path += '/';
	data_path = path;
	
	// save results to root file
	if (out_file_name.empty()) out_file_name = "out.root";
	printf("+++ saving analysis results in '%s' ...\n\n", out_file_name.c_str());
	root_out = TFile::Open(out_file_name.c_str(), "recreate");

	// invert channels
	auto polarity_map = PolarityMap(change_polarity, change_sign_from_to_ch_num);

	// verbosity
	bool debug_header = debug;
	bool debug_data = debug;

	unsigned short output_channel;
	unsigned int output_event;
	unsigned short output_nbchannels;
	unsigned short read_channels = 0;

	//Start reading the raw data from .bin files.
	stringstream inFileList;
	inFileList << Helpers::ListFiles(path.c_str(), ".bin"); //all *.bin* files in folder path
	if (debug) cout << inFileList.str() << endl;
	string fileName;
	int file_counter = 0;
	int wfcounter = 0;
	int event_counter = 0;

	// file loop
	while (inFileList >> fileName) {
		// read only fraction/batch of the .bin files for testing or to reduce memory usage
		if (FirstBinFileToRead > 0 && FirstBinFileToRead < LastBinFileToRead && file_counter < FirstBinFileToRead) continue;
		if (LastBinFileToRead > 0 && file_counter > LastBinFileToRead) break;

		fileName = path + fileName;
		ifstream input_file(fileName.c_str(), ios::binary | ios::in);

		bool has_measurement = false;

		if (!input_file.is_open()) {
			printf("*** failed to open '%s'\n", fileName.c_str());
			continue;
		}

		if (file_counter < 10 || file_counter % 10 == 0 || debug) printf("+++ reading '%s' ...\n", fileName.c_str());

		// Header
		string header_line;
		// HEADER 1 //
		//
		// "=== DATA FILE SAVED WITH SOFTWARE VERSION: V?.??.? ==="
		//
		getline(input_file, header_line, '\n');

		if (debug_header) printf("%s\n", header_line.data());
		assert(header_line[0] == '=');

		size_t header_version_first = header_line.find_last_of('V');
		size_t header_version_last = header_line.find_first_of(' ', header_version_first);
		string software_version = header_line.substr(header_version_first, header_version_last - header_version_first);
		if (debug_header) printf("    |- data version = '%s'\n", software_version.data());
		// convert software version
		software_version.erase(0, 1);
		int v_major, v_minor, v_patch;
		istringstream software_version_iss(software_version);
		char dot_;
		software_version_iss >> v_major >> dot_ >> v_minor >> dot_ >> v_patch;

		// HEADER 2 //
		// "=== WAVECATCHER SYSTEM OF TYPE ?? WITH ?? CHANNELS AND GAIN: ??? ==="
		getline(input_file, header_line, '\n');

		if (debug_header) printf("%s\n", header_line.data());
		assert(header_line[0] == '=');

		// HEADER 3 //
		// === Rate coincidence masks ... === Posttrig in ns for SamBlock ... ===
		getline(input_file, header_line, '\n');

		if (debug_header) printf("%s\n", header_line.data());
		assert(header_line[0] == '=');

		// HEADER 4 //
		// V2.9.13: === DATA SAMPLES [1024] in Volts == NB OF CHANNELS ACQUIRED: 64 == Sampling Period: 312.5 ps == INL Correction: 1
		// V2.9.15: === DATA SAMPLES [1024] in Volts == NB OF CHANNELS ACQUIRED: 64 == Sampling Period: 312.5 ps == INL Correction: 1 == MEASUREMENTS: 0 ===
		getline(input_file, header_line, '\n');

		if (debug_header) printf("%s\n", header_line.data());
		assert(header_line[0] == '=');

		size_t nsamples_first = 1 + header_line.find_last_of('[');
		size_t nsamples_last = header_line.find_first_of(']', nsamples_first);
		string nsamples_str = header_line.substr(nsamples_first, nsamples_last - nsamples_first);

		binNumber = atoi(nsamples_str.data());
		if (debug_header) printf("    |- data sample  = %d\n", binNumber);
		if (file_counter == 0 && binNumber != 1024) {
			cout << "\nWARNING: Measurement has " << binNumber << " samples, which is non-standard. Please report any bugs!" << endl;
			cout << "If this was not intentional check the WaveCatcher settings!" << endl;
		}
		size_t waveform_bytes = static_cast<size_t>(binNumber) * sizeof(short);
		vector<short> waveform(binNumber);
		vector<float> waveform_f(binNumber);

		if (file_counter == 0) amplValuessum.resize(nChannelsWC, vector<float>(binNumber, 0.));

		size_t nchannels_first = 10 + header_line.find("ACQUIRED: ", nsamples_first);
		size_t nchannels_last = header_line.find_first_of(' ', nchannels_first);
		string nchannels_str = header_line.substr(nchannels_first, nchannels_last - nchannels_first);

		nchannels = atoi(nchannels_str.data());
		if (debug_header) printf("    |- nchannels    = %d\n", nchannels);

		size_t sp_first = 8 + header_line.find("Period:");
		size_t sp_last = header_line.find(" ps");
		float sampling_period = atof(header_line.substr(sp_first, sp_last - sp_first).data());
		if (debug_header) printf("sampling period  = %f ps\n", sampling_period);
		SP = sampling_period * 1e-3;

		// compatibility with older WC software versions
		if (v_major ==2 && v_minor == 9 && v_patch <= 13) {
			// V2.9.13 has always measurement stored (everything is set to 0 when disabled!)
			has_measurement = true; 
		}
		else {
			size_t has_measurement_first = 14 + header_line.find("MEASUREMENTS: ", nsamples_first);
			size_t has_measurement_last = header_line.find_first_of(' ', has_measurement_first);
			string has_measurement_str = header_line.substr(has_measurement_first, has_measurement_last - has_measurement_first);
			has_measurement = atoi(has_measurement_str.data());
		}

		if (debug_header) printf("    `- measurement  = %d\n", has_measurement);

		// end of header reader

		event_data an_event;

		while (input_file.read((char*)(&an_event), sizeof(an_event))) {
			//event loop
			if (debug_data) printf("%03d has %d channels\n", an_event.EventNumber, an_event.nchannelstored);

			output_event = an_event.EventNumber;
			output_nbchannels = an_event.nchannelstored;

			if (debug_data && output_event % 200 == 0) printf("EventNr: %d, nCh: %d\n", output_event, output_nbchannels);
			if (output_nbchannels > nChannelsWC) {
				cout << "ERROR:\nThe number of channels in the data is " << output_nbchannels 
					<< ", which is larger than the maximum allowed number of channels which is set to " << nChannelsWC 
					<< "\nPlease set the parameter nChannelsWC=" << output_nbchannels << endl;
			}

			// do analysis only for limited range of channels to reduce memory usage for large datasets with many channels and many events
			int start_at_ch = 0;
			if (start_read_at_channel < output_nbchannels && start_read_at_channel >= 0) start_at_ch = start_read_at_channel;
			int end_at_ch = output_nbchannels - 1;
			if (end_read_at_channel == -1 && start_read_at_channel != -1) end_read_at_channel = start_read_at_channel;
			else if (end_read_at_channel < output_nbchannels && end_read_at_channel >= 0) end_at_ch = end_read_at_channel;
			read_channels = end_at_ch - start_at_ch + 1;

			if (event_counter == 0) cout << "\nstart at ch " << start_at_ch << " end at ch " << end_at_ch << endl;

			for (int ch = 0; ch < output_nbchannels; ++ch) { // channel loop
				channel_data_with_measurement a_channel_data;
				channel_data_without_measurement a_channel_data_without_measurement;
				
				if (has_measurement) { // read with 'channel_data_with_measurement' struct
					input_file.read((char*)(&a_channel_data), sizeof(channel_data_with_measurement));
				}
				else { // read with 'channel_data_without_measurement' struct
					input_file.read((char*)(&a_channel_data_without_measurement), sizeof(channel_data_without_measurement));

					// copy the content
					a_channel_data.channel = a_channel_data_without_measurement.channel;
					a_channel_data.EventIDsamIndex = a_channel_data_without_measurement.EventIDsamIndex;
					a_channel_data.FirstCellToPlotsamIndex = a_channel_data_without_measurement.FirstCellToPlotsamIndex;
				}

				output_channel = a_channel_data.channel;
				if (debug_data) cout << "- reading channel " << output_channel << endl;
				
				// read waveform
				input_file.read((char*)waveform.data(), waveform_bytes);

				//---------------------------------------------------------------------------------------------------------------
				if (ch >= start_at_ch && ch <= end_at_ch) {
					if (event_counter == 0) active_channels.push_back(static_cast<int>(output_channel));

					float factor = DAQ_factor;
					if (polarity_map[static_cast<int>(output_channel)]) factor = -factor;

					// loop to fill waveform histograms
					for (int i = 0; i < binNumber; i++) {
						waveform_f[i] = static_cast<float>(waveform[i]) * factor;
						amplValuessum[static_cast<int>(output_channel)][i] += waveform_f[i];
					}
					rundata.push_back(waveform_f);

					// baseline correction
					if (Using_BaselineCorrection_in_file_loop) {
						CorrectBaseline_function(rundata[wfcounter], tCutg, tCutEndg, wfcounter);
					}

					wfcounter++;
				}//--------------------------------------------------------------------------------------------------------------

			} // for ch

			skip_event.push_back(false);
			if (!discard_original_eventnr) eventnr_storage.push_back(output_event); // Stores the current WaveCatcher event number
			else eventnr_storage.push_back(event_counter);
			event_counter++;
			if (event_counter >= max_nevents_to_read) {
				cout << "Stopped reading events after max_nevents_to_read=" << max_nevents_to_read << endl;
				break;
			}
		} // while an_event

		input_file.close();
		file_counter++;
	} // for file_id

	// in case there are empty channels, nchannels is the number of channels which contain data
	nchannels = read_channels;

	// get bins where the sum spectrum has its maximum for runs with fixed trigger delay and fixed 
	// integration window relative to the max of the sum spectrum (not working for DC measurement)
	for (int ch = 0; ch < nChannelsWC; ch++) {
		if (Helpers::Contains(active_channels, ch)) {
			float max_val = -9.e99;
			int i_max = 0;
			for (int i = 0; i < binNumber; i++) {
				if (amplValuessum[ch][i] > max_val) {
					max_val = amplValuessum[ch][i];
					i_max = i;
				}
			}
			maxSumBin.push_back(i_max);
		}
	}

	nevents = event_counter;
	nwf = wfcounter;

	printf("Finished reading %d files containing %d events with %d channels.\n\n", file_counter, nevents, nchannels);
}

/// @brief Destructor
ReadRun::~ReadRun() {
	//plot_active_channels.clear();
	if (root_out->IsOpen()) root_out->Close();
	cout << "\nAnalysis completed." << endl;
}

/// @brief Plot sums of all raw waveforms for each channel
/// 
/// To plot the average waveforms after baseline correction etc. use PlotChannelAverages().
/// 
/// \image html PlotChannelSums.png "This function plots the sum of all waveforms for each channel, without any corrections. Channel 9 was measured with an offset as visible here. Code in example." width=75%
/// 
/// @param smooth If true it will apply smoothing to plots. \n 
/// Do not use without very good reason as it biases the results.
/// @param normalize If true will normalize the maximum to 1.
/// @param shift Shift histogram by "shift" ns
/// @param sigma Number of bins before and after central bin for running average OR gauss sigma in ns for gauss kernel and convolution.
/// @param smooth_method 0: Use running average (box kernel smoothing). Simple, very fast. \n 
/// 2: Use 3 sigma gaussian kernel smoothing. Preferred method, fast.
void ReadRun::PlotChannelSums(bool smooth, bool normalize, double shift, double sigma, int smooth_method) {

	double* xv = getx<double>(shift);
	auto mgsums = new TMultiGraph();
	mgsums->SetTitle("channel sums; t [ns]; amplitude [mV]");
	if (normalize) mgsums->SetTitle("channel sums; t [ns]; amplitude [arb.]");

	double max_val = 0., min_val = 0.;

	for (int i = 0; i < nchannels; i++) {
		if (PlotChannel(i)) {
			double* yv = new double[binNumber];
			for (int k = 0; k < binNumber; k++) yv[k] = static_cast<double>(amplValuessum[active_channels[i]][k]);

			if (smooth) Filters::SmoothArray(yv, binNumber, sigma, smooth_method);

			TGraph* gr = new TGraph(binNumber, xv, yv);
			delete[] yv;

			double tmp_min = TMath::MinElement(gr->GetN(), gr->GetY());
			if (tmp_min < min_val) min_val = tmp_min;
			double tmp_max = TMath::MaxElement(gr->GetN(), gr->GetY());
			if (tmp_max > max_val) max_val = tmp_max;
			if (normalize) {
				for (int j = 0; j < gr->GetN(); j++) gr->SetPointY(j, gr->GetPointY(j) / max_val);
			}

			TString name(Form("channel_%02d", active_channels[i]));
			TString title(Form("Channel %d", active_channels[i]));
			gr->SetName(name.Data());
			gr->SetTitle(title.Data());
			gr->SetLineColor(Helpers::rcolor(i));
			gr->SetMarkerColor(Helpers::rcolor(i));
			mgsums->Add(gr);
		}
	}
	delete[] xv;

	auto sumc = new TCanvas("Sums", "", 600, 400);
	mgsums->Draw("AL");
	if (normalize) mgsums->GetYaxis()->SetRangeUser(-0.2, 1);
	else mgsums->GetYaxis()->SetRangeUser(min_val, max_val);
	sumc->BuildLegend(0.85, 0.70, .99, .95);
	root_out->WriteObject(mgsums, "channelsums");
	root_out->WriteObject(sumc, "channelsums_c");
}
/// @example read_exampledata.cc
/// @example read_exampledata.py

/// @brief Plot averages only of the good, corrected waveforms for each channel
/// 
/// Similar to PlotChannelSums(), but will average all non-skipped waveforms. \n
/// Can be used to inspect average waveforms after baseline correction etc. has been applied. 
/// To do so, call function after calling correction and event filter functions.
/// 
/// \image html PlotChannelAverages.png "This function plots the sum of all non-skipped waveforms for each channel, with corrections. Compare with the result for PlotChannelSums(). Code in example." width=75%
/// 
/// @param normalize If true will normalize the maximum to 1.
void ReadRun::PlotChannelAverages(bool normalize) {
	float* xv = getx<float>();
	
	auto mgav = new TMultiGraph();
	mgav->SetTitle("channel averages; t [ns]; amplitude [mV]");
	if (normalize) mgav->SetTitle("channel averages; t[ns]; amplitude[arb.]");

	float max_val = 0., min_val = 0.;

	for (int i = 0; i < nchannels; i++) {
		if (PlotChannel(i)) {

			float* yv = new float[binNumber]();

			for (int j = 0; j < nevents; j++) {
				if (!SkipEvent(j, i)) {
					for (int k = 0; k < binNumber; k++) yv[k] += rundata[GetWaveformIndex(j, i)][k];
				}
			}

			float norm = max(1.f, static_cast<float>(Nevents_good(i)));
			for (int k = 0; k < binNumber; k++) yv[k] /= norm;

			auto gr = new TGraph(binNumber, xv, yv);
			delete[] yv;

			double tmp_min = TMath::MinElement(gr->GetN(), gr->GetY());
			if (tmp_min < min_val) min_val = tmp_min;
			double tmp_max = TMath::MaxElement(gr->GetN(), gr->GetY());
			if (tmp_max > max_val) max_val = tmp_max;
			if (normalize) {
				for (int j = 0; j < gr->GetN(); j++) gr->SetPointY(j, gr->GetPointY(j) / tmp_max);
			}

			TString name(Form("channel_%02d", active_channels[i]));
			TString title(Form("Channel %d", active_channels[i]));
			gr->SetName(name.Data());
			gr->SetTitle(title.Data());
			gr->SetLineColor(Helpers::rcolor(i));
			gr->SetMarkerColor(Helpers::rcolor(i));
			mgav->Add(gr);
		}
	}
	delete[] xv;

	string cname("Averages_" + to_string(PlotChannelAverages_cnt++));
	auto avc = new TCanvas(cname.c_str(), cname.c_str(), 600, 400);
	mgav->Draw("AL");
	if (normalize) mgav->GetYaxis()->SetRangeUser(-0.2, 1);
	else mgav->GetYaxis()->SetRangeUser(min_val, max_val);
	avc->BuildLegend(0.85, 0.70, .99, .95);
	root_out->WriteObject(mgav, ("channelaverages" + to_string(PlotChannelAverages_cnt)).c_str());
	root_out->WriteObject(avc, ("channelaverages_c" + to_string(PlotChannelAverages_cnt)).c_str());
}
/// @example timing_example.cc
/// @example read_exampledata.cc
/// @example read_exampledata.py


/// @brief 2D histogram of all non-skipped waveforms for one channel
/// 
/// See PlotWFHeatmaps()
/// 
/// @param channel_index Index of the channel
/// @param ymin Min. y range
/// @param ymax Max. y range
/// @param n_bins_y Number of bins along y
/// @returns TH2F
TH2F* ReadRun::WFHeatmapChannel(int channel_index, float ymin, float ymax, int n_bins_y) {

	TString name(Form("channel__%02d", active_channels[channel_index]));
	TH2F* h2 = new TH2F(name.Data(), name.Data(), binNumber, 0, SP * static_cast<float>(binNumber), n_bins_y, ymin, ymax);
	h2->GetXaxis()->SetTitle("t [ns]");
	h2->GetYaxis()->SetTitle("I [arb.]");
	h2->GetZaxis()->SetTitle("entries");

	auto xv = getx<float>();
	#pragma omp parallel
	{
		auto h2_tmp = new TH2F("", "", binNumber, 0, SP * static_cast<float>(binNumber), n_bins_y, ymin, ymax);
		#pragma omp for
		for (int j = 0; j < nevents; j++) {
			if (!SkipEvent(j, channel_index)) {
				int wf_index = GetWaveformIndex(j, channel_index);
				for (int i = 0; i < binNumber; i++) h2_tmp->Fill(xv[i], rundata[wf_index][i]);
			}
		}
		#pragma omp critical
		h2->Add(h2_tmp);
		delete h2_tmp;
	}
	delete[] xv;
	return h2;
}

/// @brief Plot stacks of all non-skipped waveforms for all active channels
/// 
/// Classical SiPM dark count/calibration plot in the style of a phosphor oscilloscope screen. \n
/// Helpful to see the prevalence of certain signal/baseline shapes etc. \n
/// Take care not to plot too many channels at once to maintain visibility and performance (-> ReadRun::plot_active_channels).
/// 
/// \image html heatmaps.png "Waveforms of all events filled into a single histogram. Code in example." width=75%
/// 
/// @param ymin Min. y range
/// @param ymax Max. y range
/// @param n_bins_y Number of bins along y
/// @param z_opt "log": Set z logarithmic. \n
/// "COLZ0": Do not plot bins with zero entries.
/// @param z_max Max. z range. Set to zero to use auto range (default).
/// @param palette Root color palette (see https://root.cern.ch/doc/master/classTColor.html).
void ReadRun::PlotWFHeatmaps(float ymin, float ymax, int n_bins_y, string z_opt, float z_max, EColorPalette palette) {

    gStyle->SetOptStat(0);
    string name("waveforms_heatmap_" + to_string(PlotWFHeatmaps_cnt++));
    auto wfhm_c = new TCanvas(name.c_str(), name.c_str(), 600, 400);
    Helpers::SplitCanvas(wfhm_c, active_channels, plot_active_channels);

    // create histograms in parallel
	map<int, shared_future<TH2F*>> h2_future;
    for (int i = 0; i < nchannels; ++i) {
        if (PlotChannel(i)) {
            h2_future[i] = async(launch::async, [this, i, ymin, ymax, n_bins_y]() {
                return WFHeatmapChannel(i, ymin, ymax, n_bins_y);
            });
        }
    }

	// plot histograms
    int current_canvas = 0;
    for (int i = 0; i < nchannels; ++i) {
        if (PlotChannel(i)) {
            wfhm_c->cd(++current_canvas);

            gPad->SetTopMargin(.1);
            gPad->SetBottomMargin(.1);
            gPad->SetLeftMargin(.15);
            gPad->SetRightMargin(.15);
            gStyle->SetPalette(palette);

            auto h2 = h2_future[i].get();

            h2->SetContour(99);
            h2->SetStats(0);
            if (z_opt == "COLZ0") h2->Draw("COLZ0");
            else h2->Draw("CONT4Z");

            if (z_opt == "log") {
                gPad->SetLogz();
                if (z_max > 1) h2->GetZaxis()->SetRangeUser(1, z_max);
            }
            else if (z_max > 0) h2->GetZaxis()->SetRangeUser(0, z_max);
        }
    }

    wfhm_c->Update();
    root_out->WriteObject(wfhm_c, name.c_str());
}
/// @example read_exampledata.cc

/// @brief Smoothing all waveforms which are not skipped (for testing, careful when using for analysis!)
/// 
/// See Filters::SmoothArray().
/// 
/// @param sigma Number of bins before and after central bin for running average OR gauss sigma in ns.
/// @param method 0: Use running average (box kernel smoothing). Simple, very fast. \n 
/// 2: Use 3 sigma gaussian kernel smoothing. Preferred method, fast.
void ReadRun::SmoothAll(double sigma, int method) {
	cout << "Smoothing all non-skipped waveforms..." << endl;

	#pragma omp parallel for
	for (int j = 0; j < nwf; j++) {
		if (!SkipEvent(GetCurrentEvent(j), GetCurrentChannel(j))) {
			vector<double> tmp_(rundata[j].begin(), rundata[j].end());
			double* tmp_wf = tmp_.data();

			Filters::SmoothArray(tmp_wf, binNumber, sigma, method, SP);
			for (int i = 0; i < binNumber; ++i) rundata[j][i] = static_cast<float>(tmp_[i]);
		}
	}
}

/// @brief Smoothing all waveforms which are not skipped (for testing, careful when using for analysis!)
/// 
/// See Filters::SmoothArray().
/// 
/// @param sigma Number of bins before and after central bin for running average OR gauss sigma in ns.
/// @param method "Box": Use running average (box kernel smoothing). Simple, very fast. \n
/// "Gaus": Use 3 sigma gaussian kernel smoothing. Preferred method, fast. 
void ReadRun::SmoothAll(double sigma, string method) {
	cout << "\nSmoothing all non-skipped waveforms:" << endl;

	#pragma omp parallel for
	for (int j = 0; j < nwf; j++) {
		if (!SkipEvent(GetCurrentEvent(j), GetCurrentChannel(j))) {
			vector<double> tmp_(rundata[j].begin(), rundata[j].end());
			double* tmp_wf = tmp_.data();

			Filters::SmoothArray(tmp_wf, binNumber, sigma, method, SP);
			for (int i = 0; i < binNumber; ++i) rundata[j][i] = static_cast<float>(tmp_[i]);
		}
	}
}

/// @brief Filter all waveforms
/// 
/// Experimental. See Filters::ResponseFilter() or Filters::SecondOrderUnderdampedFilter(). \n
/// 
/// @param sigma1 First.
/// @param sigma2 Second.
/// @param factor Factor for negative part (0 < factor < 1). 
void ReadRun::FilterAll(double sigma1, double sigma2, double factor) {
	cout << "\nFiltering all waveforms..." << endl;

	#pragma omp parallel for
	for (int j = 0; j < nwf; j++) {
		vector<double> tmp_(rundata[j].begin(), rundata[j].end());
		double* tmp_wf = tmp_.data();

		if (factor > 0 && factor <= 1) Filters::ResponseFilter(tmp_wf, binNumber, sigma1, sigma2, factor, SP);
		else Filters::SecondOrderUnderdampedFilter(tmp_wf, binNumber, sigma1, sigma2, abs(factor), SP);
		for (int i = 0; i < binNumber; ++i) rundata[j][i] = static_cast<float>(tmp_[i]);
	}
}

/// @brief This function shifts all waveforms to the average signal starting times for each channel.
/// 
/// The signal starting times are determined with constant fraction discrimination. 
/// **Before** calling this function, please call GetTimingCFD() with suitable parameters for your data. \n
/// Additionally, PrintChargeSpectrumWF() should be called **before** calling this function since the 
/// timing reference (blue line) won't be shifted.
void ReadRun::ShiftAllToAverageCF() {
	cout << "\nShifting all waveforms to the average constant fraction time for each channel:" << endl;

	//call GetTimingCFD() in case it was not initialized
	if (static_cast<int>(timing_results.size()) == 0) GetTimingCFD();

	double* timing_mean = new double[nchannels]();

	for (int j = 0; j < nwf; j++) {
		int curr_ch = GetCurrentChannel(j);
		if (!SkipEvent(GetCurrentEvent(j), curr_ch)) timing_mean[curr_ch] += timing_results[j][0];
	}

	int* timing_mean_n = new int[nchannels];
	for (int i = 0; i < nchannels; i++) {
		double norm = max(1., static_cast<double>(Nevents_good(i)));
		timing_mean_n[i] = static_cast<int>(round(timing_mean[i] / norm));
	}
	delete[] timing_mean;

	#pragma omp parallel for
	for (int j = 0; j < nwf; j++) {
		int curr_ch = GetCurrentChannel(j);
		if (!SkipEvent(GetCurrentEvent(j), curr_ch)) {
			int shift = static_cast<int>(timing_results[j][0]) - timing_mean_n[curr_ch];
			shift = shift % binNumber;
    		if (shift < 0) shift += binNumber;
    		rotate(rundata[j].begin(), rundata[j].begin() + shift, rundata[j].end());
		}
	}
	delete[] timing_mean_n;
}

/// @brief Baseline correction constant window
/// 
/// Corrects the baseline (DC offset) of all waveforms. \n 
/// Uses the mean between t=0 and t="tCut" or between t="tCut" and t="tCutEnd" as offset. \n
/// Call method before ReadFile() if you want it to happen while reading. \n \n 
/// 
/// Most simple and fast method. \n 
/// Useful for measurements with very few background events/dark counts. \n 
/// Using a constant window for the baseline means there might be a background pulse which would lead to wrong correction. \n \n 
/// 
/// Stores results for all channels and all events in ReadRun::baseline_correction_result. \n 
/// Results will be visualized for each event in PrintChargeSpectrumWF(). 
/// 
/// @param tCut Time denoting the end or the beginning (if "tCutEnd" is set) of the integration window.
/// @param tCutEnd Time denoting the end of the integration window.
void ReadRun::CorrectBaseline(float tCut, float tCutEnd) {
	checkData(true);

	cout << "\nPerforming simple baseline correction in fixed time window."
		<< "This method is only suitable for measurements without dark counts!" << endl;
	tCutg = tCut;
	tCutEndg = tCutEnd;
	if (nwf == 0) {
		Using_BaselineCorrection_in_file_loop = true;
	}
	else {
		cout << "Baseline correction (" << nwf << " waveforms):" << endl;
		baseline_correction_result.resize(nwf, vector<float>(4));

		#pragma omp parallel for
		for (int j = 0; j < nwf; j++) {
			CorrectBaseline_function(rundata[j], tCut, tCutEnd, j);
		}
	}
}

/// @brief Helper function called by CorrectBaseline()
/// 
/// See CorrectBaseline()
/// 
void ReadRun::CorrectBaseline_function(vector<float>& waveform, float tCut, float tCutEnd, int waveform_index) {
	int iCut, iCutEnd;
	float corr = 0;

	iCut = TimeToIndex(tCut);

	if (tCutEnd <= 0) { //
		for (int i=0; i<=iCut; i++) corr += waveform[i];
		corr /= static_cast<float>(iCut);
	}
	else {
		iCutEnd = TimeToIndex(tCutEnd);
		for (int i=iCut; i<=iCutEnd; i++) corr += waveform[i]; 
		corr /= static_cast<float>(iCutEnd - iCut + 1);
	}

	// write corrected values to histograms
	if (tCut >= 0) {
		for (int i = 0; i < binNumber; i++) waveform[i] -= corr;
	}

	if (!Using_BaselineCorrection_in_file_loop) {
		baseline_correction_result[waveform_index][0] = corr;
		baseline_correction_result[waveform_index][1] = 0;
		baseline_correction_result[waveform_index][2] = tCut;
		baseline_correction_result[waveform_index][3] = tCutEnd;
	}
}

/// @brief Baseline correction method searching for non-monotonic, rather constant regions
/// 
/// Corrects the baseline (DC offset) of all waveforms. \n 
/// Determines the region of window[0] ns between window[1] ns and window[2] ns where the squared sum plus 
/// the square of the sum of the slope of the (smoothed) waveform reaches its minimum: \n \n
/// \f$\mathbf{min}\left( \sum \left(\Delta y_i \right)^2 + \left(\sum \Delta y_i \right)^2 \right) \f$ \n \n
/// 
/// Here, \f$\sum \left(\Delta y_i \right)^2 \to 0\f$ if the region is constant and 
/// \f$\left( \sum \Delta y_i \right)^2 \to 0\f$ if the region is constant or oscillating around a constant value. 
/// The second term penalizes monotonic regions with a small, but rather constant slope (e. g. tails). \n \n 
/// 
/// Slower, but versatile since it searches for the optimal baseline candidate region in a defined range. \n 
/// Will prefer constant sections of the waveform for the estimation of the baseline. \n
/// Not well suited if there is not constant baseline in the signal, which can happen if the dark count rate 
/// is so high that dark counts overlap (e. g. an array of SiPMs) or if the baseline level fluctuates. 
/// In such a case use CorrectBaselineMin() \n \n 
/// 
/// Stores results for all channels and all events in ReadRun::baseline_correction_result. \n 
/// Results will be visualized for each event in PrintChargeSpectrumWF(). \n \n 
/// 
/// @param window Vector containing {length for averaging, search start, search end} in ns. 
/// Example: {20, 10, 90} would search for the best baseline candidate from 10 ns to 90 ns, averaging over 20 ns.
/// @param sigma Number of bins before and after central bin for running average OR gauss sigma in ns for gauss
/// @param smooth_method 0: Use running average (box kernel smoothing). Simple, very fast. \n 
/// 2: Use 3 sigma gaussian kernel smoothing. Preferred method, fast.
void ReadRun::CorrectBaselineMinSlopeRMS(vector<float> window, double sigma, int smooth_method) {
	checkData(true);

	cout << "\nBaseline correction (minimum slope variation method, " << nwf << " waveforms):" << endl;
	if (window.empty()) cout << "\nWarning: Window not set in CorrectBaselineMinSlopeRMS. Will use default values." << endl;
	if (sigma != 0.) cout << "\nNotification: Using smoothing in CorrectBaselineMinSlopeRMS." << endl;

	int nbins_average = !window.empty() ? TimeToIndex(window[0]) : TimeToIndex(50.);
	int start_search_at = static_cast<int>(window.size()) > 1 ? TimeToIndex(window[1]) : 0;
	int end_search_at = static_cast<int>(window.size()) > 2 ? TimeToIndex(window[2]) : binNumber - 1;

	int end_search_loop_at = end_search_at - start_search_at - nbins_average;

	// if no valid static search window is specified, it will be dynamic from 0 ns up to 8 ns before the global maximum
	bool search_relative_to_local_max = false;
	int min_distance_from_max = TimeToIndex(8.);
	if (start_search_at < 0) {
		search_relative_to_local_max = true;
		start_search_at = 0;
		cout << "\nNotification: Using dynamic search window in CorrectBaselineMinSlopeRMS." << endl;
	}

	int nbins_search = end_search_at - start_search_at;

	baseline_correction_result.resize(nwf, vector<float>(6));

	#pragma omp parallel for
	for (int j = 0; j < nwf; j++) {
		float minchange = 1.e99;
		float sum = 0, sumsum = 0, change = 0, minsumsq = 0, sqsum = 0, minsqsum = 0, corr = 0;
		int iintwindowstart = 0, imax = 0;
		int nbins_search_ = nbins_search;
		int end_search_at_ = end_search_at;
		int end_search_loop_at_ = end_search_loop_at;

		if (search_relative_to_local_max) {
			imax = GetIntWindow(rundata[j], 0, 0, nbins_search_ + min_distance_from_max, end_search_at_)[1];
			end_search_at_ = imax - min_distance_from_max;
			nbins_search_ = end_search_at_; // starts at 0
			end_search_loop_at_ = nbins_search_ - nbins_average;
		}

		double* yvals = Helpers::gety(rundata[j], start_search_at, end_search_at_);
		// smoothing suppresses variations in slope due to noise, so the method is potentially more sensitive to excluding peaks
		if (sigma > 0) Filters::SmoothArray(yvals, nbins_search_, sigma, smooth_method);
		//calculate slope
		nbins_search_--;
		double* slope = new double[nbins_search_];
		double* slope_sq = new double[nbins_search_];
		for (int i = 0; i < nbins_search_; i++) {
			slope[i] = yvals[i + 1] - yvals[i];
			slope_sq[i] = slope[i] * slope[i];
			if (i < nbins_average) { // init
				sum += slope[i];
				sqsum += slope_sq[i];
			}
		}
		delete[] yvals;

		//find window for correction
		for (int i = 0; i < end_search_loop_at_; i++) {
			sumsum = sum * sum;
			change = sqsum + sumsum;

			if (change < minchange) {
				minchange = change;
				iintwindowstart = i + start_search_at;
				minsumsq = sumsum;
				minsqsum = sqsum;
			}

			sum -= slope[i];
			sum += slope[i + nbins_average];
			sqsum -= slope_sq[i];
			sqsum += slope_sq[i + nbins_average];
		}
		delete[] slope;
		delete[] slope_sq;

		// do correction
		for (int i=iintwindowstart; i<=iintwindowstart + nbins_average; i++) corr += rundata[j][i];
		corr /= static_cast<float>(nbins_average + 1);
		for (int i = 0; i < binNumber; i++) rundata[j][i] -= corr;

		baseline_correction_result[j][0] = corr;
		baseline_correction_result[j][1] = minchange;
		baseline_correction_result[j][2] = IndexToTime(iintwindowstart);
		baseline_correction_result[j][3] = IndexToTime(iintwindowstart + nbins_average);
		baseline_correction_result[j][4] = minsumsq;
		baseline_correction_result[j][5] = minsqsum;
	}
}
/// @example read_exampledata.cc
/// @example read_exampledata.py

/// @brief Baseline correction method searching for non-monotonic, rather constant regions
/// 
/// This is a deprecated version of CorrectBaselineMinSlopeRMS. It will be removed in future releases.
/// 
/// Corrects the baseline (DC offset) of all waveforms. \n 
/// Determines the region of "nIntegrationWindow" bins where the squared sum plus the square of the sum of 
/// the slope of the (smoothed) waveform reaches its minimum: \n \n
/// \f$\mathbf{min}\left( \sum \left(\Delta y_i \right)^2 + \left(\sum \Delta y_i \right)^2 \right) \f$ \n \n
/// 
/// Here, \f$\sum \left(\Delta y_i \right)^2 \to 0\f$ if the region is constant and 
/// \f$\left( \sum \Delta y_i \right)^2 \to 0\f$ if the region is constant or oscillating around a constant value. 
/// The second term penalizes monotonic regions with a small, but rather constant slope (e. g. tails). \n \n 
/// 
/// Slow, but versatile since it searches for the optimal baseline candidate region in a defined range. \n 
/// Will prefer constant sections of the waveform for the estimation of the baseline. \n
/// Not well suited if there is not constant baseline in the signal, which can happen if the dark count rate 
/// is so high that dark counts overlap (e. g. an array of SiPMs) or if the baseline level fluctuates. 
/// In such a case use CorrectBaselineMin() \n \n 
/// 
/// Stores results for all channels and all events in ReadRun::baseline_correction_result. \n 
/// Results will be visualized for each event in PrintChargeSpectrumWF(). \n \n 
/// 
/// @param nIntegrationWindow Number of bins used for baseline correction. The correction factor will be the 
/// signal averaged over this number bins.
/// @param smooth Deprecated! If true will apply smoothing to all waveforms. This will change all waveforms, 
/// which should rather be done with SmoothAll(). 
/// @param sigma Number of bins before and after central bin for running average OR gauss sigma in ns for gauss 
/// kernel and convolution. Set to 0 for no smoothing. Use with care!
/// @param max_bin_for_baseline Maximum bin for search window.
/// @param start_at Minimum bin for search window.
/// @param smooth_method 0: Use running average (box kernel smoothing). Simple, very fast. \n 
/// 2: Use 3 sigma gaussian kernel smoothing. Preferred method, fast.
void ReadRun::CorrectBaselineMinSlopeRMS(int nIntegrationWindow, bool smooth, double sigma, int max_bin_for_baseline, int start_at, int smooth_method) {
	cout << "WARNING: This is a deprecated version of CorrectBaselineMinSlopeRMS. "
		<< "It will be removed in future releases. Parameter bool smooth=" << smooth << " will be ignored." << endl;
	vector<float> window;
	window.push_back(IndexToTime(nIntegrationWindow));
	window.push_back(IndexToTime(start_at));
	window.push_back(IndexToTime(max_bin_for_baseline));
	CorrectBaselineMinSlopeRMS(window, sigma, smooth_method);
}
/// @example read_exampledata.cc
/// @example read_exampledata.py

/// @brief Wrapper for backwards compatibility
/// @param window 
/// @param sigma 
/// @param smooth_method 
/// @param increment Parameter has been removed
void ReadRun::CorrectBaselineMinSlopeRMS(vector<float> window, double sigma, int smooth_method, int increment) {
	(void)increment;
	cout << "WARNING: This is a deprecated version of CorrectBaselineMinSlopeRMS. "
		<< "It will be removed in future releases. Parameter increment=" << increment << " will be ignored." << endl;
	CorrectBaselineMinSlopeRMS(window, sigma, smooth_method);
}

/// @brief Baseline correction using minimum sum (\f$\propto\f$ mean) in range for correction 
/// 
/// Corrects the baseline (DC offset) of all waveforms. \n 
/// Searches for \n \n
/// \f$\mathbf{min}\left( \sum y_i \right) \f$ \n \n
/// in range {window[1], window[2]}, summing over window[0] ns. \n
/// Make sure the search range is shortly before the triggered signal is expected to arrive. \n \n 
/// 
/// Helpful for (groups of/irradiated) SiPMs with very high dark count rate (DCR) where the voltage rarely relaxes 
/// back to the constant baseline before the next dark count/signal arrives: \n
/// \f$ \Rightarrow DCR \sim 1/t_{signal} \f$ \n \n
/// 
/// Stores results for all channels and all events in ReadRun::baseline_correction_result. \n 
/// Results will be visualized for each event in PrintChargeSpectrumWF().
/// 
/// @param window Vector containing {length for averaging, search start, search end} in ns. 
/// Example: {20, 10, 90} would search for the best baseline candidate from 10 ns to 90 ns, averaging over 20 ns. 
/// @param sigma Number of bins before and after central bin for running average OR gauss sigma in ns for gauss 
/// kernel and convolution. Set to 0 for no smoothing. Use with care!
/// @param smooth_method 0: Use running average (box kernel smoothing). Simple, very fast. \n 
/// 2: Use 3 sigma gaussian kernel smoothing. Preferred method, fast.
void ReadRun::CorrectBaselineMin(vector<float> window, double sigma, int smooth_method) {
	checkData(true);

	cout << "\nBaseline correction (minimal sum method, " << nwf << " waveforms):" << endl;
	if (window.empty()) cout << "\nWarning: Window not set in CorrectBaselineMin. Will use default values." << endl;
	if (sigma != 0.) cout << "\nNotification: Using smoothing in CorrectBaselineMin." << endl;

	int nbins_average = !window.empty() ? TimeToIndex(window[0]) : TimeToIndex(10.);
	int start_search_at = static_cast<int>(window.size()) > 1 ? TimeToIndex(window[1]) : 0;
	int end_search_at = static_cast<int>(window.size()) > 2 ? TimeToIndex(window[2]) : binNumber - 1;

	int end_search_loop_at = end_search_at - start_search_at - nbins_average;

	// if no valid static search window is specified, it will be dynamic from 0 ns up to 8 ns before the global maximum
	bool search_relative_to_local_max = false;
	int min_distance_from_max = TimeToIndex(8.);
	if (start_search_at < 0 || end_search_loop_at < 0) {
		search_relative_to_local_max = true;
		start_search_at = 0;
		cout << "\nNotification: Using dynamic search window in CorrectBaselineMin." << endl;
	}

	int nbins_search = end_search_at - start_search_at;

	baseline_correction_result.resize(nwf, vector<float>(4));

	#pragma omp parallel for
	for (int j = 0; j < nwf; j++) {
		float minchange = 1e9;
		float sum = 0, corr = 0;
		int iintwindowstart = 0, imax = 0;
		int nbins_search_ = nbins_search;
		int end_search_at_ = end_search_at;
		int end_search_loop_at_ = end_search_loop_at;

		if (search_relative_to_local_max) {
			imax = GetIntWindow(rundata[j], 0, 0, nbins_search_ + min_distance_from_max, end_search_at_)[1];
			end_search_at_ = imax - min_distance_from_max;
			nbins_search_ = end_search_at_;
			end_search_loop_at_ = nbins_search_ - nbins_average;
		}

		double* yvals = Helpers::gety(rundata[j], start_search_at, end_search_at_);
		// smoothing suppresses variations in slope due to noise, so the method is potentially more sensitive to excluding peaks
		if (sigma > 0) Filters::SmoothArray(yvals, nbins_search_, sigma, smooth_method);
		//find window for correction
		for (int i = 0; i < nbins_average; i++) sum += yvals[i]; // init
		for (int i = 0; i < end_search_loop_at_; i++) {
			if (sum < minchange) {
				minchange = sum;
				iintwindowstart = i + start_search_at;
			}

			sum -= yvals[i];
			sum += yvals[i + nbins_average];
		}
		delete[] yvals;
		
		// do correction
		for (int i=iintwindowstart; i<=iintwindowstart + nbins_average; i++) corr += rundata[j][i];
		corr /= static_cast<float>(nbins_average + 1);
		for (int i = 0; i < binNumber; i++) rundata[j][i] -= corr;

		baseline_correction_result[j][0] = corr;
		baseline_correction_result[j][1] = minchange;
		baseline_correction_result[j][2] = IndexToTime(iintwindowstart);
		baseline_correction_result[j][3] = IndexToTime(iintwindowstart + nbins_average);
	}
}

/// @brief Baseline correction using minimum sum (\f$\propto\f$ mean) in range for correction 
/// 
/// This is a deprecated version of CorrectBaselineMin. It will be removed in future releases.\n
/// 
/// Corrects the baseline (DC offset) of all waveforms. \n 
/// Searches for \n \n
/// \f$\mathbf{min}\left( \sum y_i \right) \f$ \n \n
/// in range {"start_at", "max_bin_for_baseline"}, summing over "nIntegrationWindow"  bins. \n
/// Make sure the search range is shortly before the triggered signal is expected to arrive. \n \n 
/// 
/// Helpful for (groups of/irradiated) SiPMs with very high dark count rate (DCR) where the voltage rarely relaxes 
/// back to the constant baseline before the next dark count/signal arrives: \n
/// \f$ \Rightarrow DCR \sim 1/t_{signal} \f$ \n \n
/// 
/// Stores results for all channels and all events in ReadRun::baseline_correction_result. \n 
/// Results will be visualized for each event in PrintChargeSpectrumWF().
/// 
/// @param nIntegrationWindow Number of bins used for baseline correction. The correction factor will be the 
/// signal averaged over this number bins.
/// @param sigma Number of bins before and after central bin for running average OR gauss sigma in ns for gauss 
/// kernel and convolution. Set to 0 for no smoothing. Use with care!
/// @param max_bin_for_baseline Maximum bin for search window.
/// @param start_at Minimum bin for search window.
/// @param smooth_method 0: Use running average (box kernel smoothing). Simple, very fast. \n 
/// 2: Use 3 sigma gaussian kernel smoothing. Preferred method, fast.
void ReadRun::CorrectBaselineMin(int nIntegrationWindow, double sigma, int max_bin_for_baseline, int start_at, int smooth_method) {
	cout << "WARNING: This is a deprecated version of CorrectBaselineMin. It will be removed in future releases." << endl;
	vector<float> window;
	window.push_back(IndexToTime(nIntegrationWindow));
	window.push_back(IndexToTime(start_at));
	window.push_back(IndexToTime(max_bin_for_baseline));
	CorrectBaselineMin(window, sigma, smooth_method);
}
/// @example timing_example.cc

/// @brief Wrapper for backwards compatibility
/// @param window 
/// @param sigma 
/// @param smooth_method 
/// @param increment Parameter has been removed
void ReadRun::CorrectBaselineMin(vector<float> window, double sigma, int smooth_method, int increment) {
	(void)increment;
	cout << "WARNING: This is a deprecated version of CorrectBaselineMin. It will be removed in future releases." << endl;
	CorrectBaselineMin(window, sigma, smooth_method);
}

/// @brief Waveform projections for one channel
/// 
/// See PrintWFProjection() for parameters.
/// 
/// @return Histogram of the projections of non-skipped waveforms for one channel.
TH1F* ReadRun::WFProjectionChannel(int channel_index, int from_n, int to_n, float rangestart, float rangeend, int nbins) {
	from_n = CheckBoundsX(from_n);
	to_n = CheckBoundsX(to_n);

	TString name(Form("channel__%02d", active_channels[channel_index]));
	auto h1 = new TH1F(name.Data(), name.Data(), nbins, rangestart, rangeend);
	
	for (int j = 0; j < nevents; j++) {
		if (!SkipEvent(j, channel_index)) {
			for (int i = from_n; i <= to_n; i++) h1->Fill(rundata[GetWaveformIndex(j, channel_index)][i]);
		}
	}
	return h1;
}

/// @brief Plots waveform projection histograms of all channels
/// 
/// Useful to check baseline correction. Will show if the baseline correction works as intended (gaussian shape) or 
/// if the standard deviation of the baseline is too large. \n \n
/// 
/// **An asymmetry** in the distribution **might point to** many dark counts in the window or a **bad choice of parameters 
/// for the baseline correction**. In that case the baseline correction should be revisited.
/// 
/// @param from Do projection in time window (from...
/// @param to ...to) in ns. This window should reflect the search window used for baseline correction.
/// @param rangestart Plot x range start in mV.
/// @param rangeend Plot x range end in mV.
/// @param nbins Number of bins in range.
void ReadRun::PrintWFProjection(float from, float to, float rangestart, float rangeend, int nbins) {
	gStyle->SetOptFit(111);
	float default_rangestart = -10;
	float default_rangeend = 20;
	if (default_rangestart > rangestart) default_rangestart = rangestart;
	if (default_rangeend < rangeend) default_rangeend = rangeend;
	int default_nbins = static_cast<int>((default_rangeend - default_rangestart) * nbins / (rangeend - rangestart));

	string ctitle("WFProjection" + to_string(PrintWFProjection_cnt++));
	auto wf_projection_c = new TCanvas(ctitle.c_str(), ctitle.c_str(), 600, 400);
	Helpers::SplitCanvas(wf_projection_c, active_channels, plot_active_channels);
	int current_canvas = 0;

	int from_n = TimeToIndex(from);
	int to_n = TimeToIndex(to);

	for (int i = 0; i < nchannels; i++) {
		if (PlotChannel(i)) {
			wf_projection_c->cd(++current_canvas);

			auto his = WFProjectionChannel(i, from_n, to_n, default_rangestart, default_rangeend, default_nbins);
			his->GetYaxis()->SetTitle("#Entries");
			his->GetXaxis()->SetTitle("amplitude in mV");
			TString name(Form("WFProjection channel_%02d_%d", active_channels[i], PrintWFProjection_cnt));
			his->Draw();
			his->Fit("gaus", "WWM", "same");
			root_out->WriteObject(his, name.Data());
		}
	}

	Helpers::SetRangeCanvas(wf_projection_c, rangestart, rangeend);
	root_out->WriteObject(wf_projection_c, ("WFProjections" + to_string(PrintWFProjection_cnt)).c_str());
}

/// @brief Histograms of the contents of baseline_correction_result
/// 
/// See PrintBaselineCorrectionResults() for parameters.
/// 
/// @return Histogram for one channel.
TH1F* ReadRun::BaselineCorrectionResults(int channel_index, int which, float rangestart, float rangeend, int nbins) {
	if (baseline_correction_result.empty() || static_cast<int>(baseline_correction_result[0].size()) < which - 1) {
		cout << "\nError: baseline_correction_result empty. Call baseline correction first." << endl;
	}
	TString name(Form("channel__%02d", active_channels[channel_index]));
	auto h1 = new TH1F(name.Data(), name.Data(), nbins, rangestart, rangeend);

	float average_baseline = 0.;
	int counter = 0;
	for (int j = 0; j < nevents; j++) {
		if (!SkipEvent(j, channel_index)) {
			float current_baseline = baseline_correction_result[GetWaveformIndex(j, channel_index)][which];
			average_baseline += current_baseline;
			counter++;
			h1->Fill(current_baseline);
		}
	}
	cout << "Mean baseline ch" << active_channels[channel_index] << ": " << average_baseline / max(1.f, static_cast<float>(counter)) << endl;
	return h1;
}

/// @brief Print histogram of the baseline correction values for all channels
/// 
/// Currently only optimized for the correction values themselves.
/// 
/// @param rangestart Plot x range start.
/// @param rangeend Plot x range end.
/// @param nbins Number of bins in range.
void ReadRun::PrintBaselineCorrectionResults(float rangestart, float rangeend, int nbins) {
	checkData();
	gStyle->SetOptFit(111);
	float default_rangestart = -10;
	float default_rangeend = 20;
	if (default_rangestart > rangestart) default_rangestart = rangestart;
	if (default_rangeend < rangeend) default_rangeend = rangeend;
	int default_nbins = static_cast<int>((default_rangeend - default_rangestart) * nbins / (rangeend - rangestart));

	string ctitle;
	ctitle = "Correction values in mV";
	auto blc_res_c = new TCanvas(ctitle.c_str(), ctitle.c_str(), 600, 400);
	Helpers::SplitCanvas(blc_res_c, active_channels, plot_active_channels);
	int current_canvas = 0;

	for (int i = 0; i < nchannels; i++) {
		if (PlotChannel(i)) {
			blc_res_c->cd(++current_canvas);

			auto his = BaselineCorrectionResults(i, 0, default_rangestart, default_rangeend, default_nbins);
			his->GetYaxis()->SetTitle("#Entries");
			his->GetXaxis()->SetTitle(ctitle.c_str());
			his->Draw();
			his->Fit("gaus", "WWM", "same");
		}
	}
	Helpers::SetRangeCanvas(blc_res_c, rangestart, rangeend);
	root_out->WriteObject(blc_res_c, ctitle.c_str());
}

/// @brief Determine the timing of the maximum peak with constant fraction discrimination
/// 
/// Determines timing in the time interval ["start_at_t", "end_at_t"] with CFD for fraction of maximum "cf_r". \n 
/// Stores timing information for all channels and all events in ReadRun::timing_results. \n 
/// Per event results will be visualized in PrintChargeSpectrumWF(). \n 
/// Cumulative per channel results can be visualized with Print_GetTimingCFD(). \n 
/// Results can be used for time difference cuts with SkipEventsTimeDiffCut().
/// 
/// @param cf_r Fraction of maximum for CFD. \n
/// If set > 1 it will be used as a fixed threshold in mV, so it will no longer use CFD.
/// @param start_at_t Time in ns to start searching.
/// @param end_at_t Time in ns to end searching.
/// @param sigma Number of bins before and after central bin for running average OR gauss sigma in ns for gauss kernel and convolution.
/// This will bias the results! Do not use (or use very carefully, only for noisy data)! Set to 0 if you do not want to use smoothing.
/// @param find_CF_from_start If true will start search from "start_at_t" to find the first arriving photon (default setting). \n 
/// If false search backwards from the time of the maximum.  
/// @param smooth_method 0: Use running average (box kernel smoothing). Simple, very fast. \n 
/// 2: Use 3 sigma gaussian kernel smoothing. Preferred method, fast.
/// @param use_spline If false will use linear interpolation between the two bins closest to cf_r. \n
/// If true will use a 5th order spline and bisection method for interpolation. 
/// Performs a bit better in most cases compared to only linear interpolation. 
/// @param verbose Print additional warnings when CFD fails
void ReadRun::GetTimingCFD(float cf_r, float start_at_t, float end_at_t, double sigma, bool find_CF_from_start, int smooth_method, bool use_spline, bool verbose) {

	int start_at = CheckBoundsX(TimeToIndex(start_at_t));
	int end_at = CheckBoundsX(TimeToIndex(end_at_t));
	int n_range = end_at - start_at;
	
	if (cf_r <= 0) cf_r = 1;
	cout << "\nGet timing at " << (cf_r > 0 && cf_r <= 1 ? "CF=" : "threshold=");
	printf("%.2f between %.2f ns and %.2f ns (%d waveforms):\n", cf_r, start_at_t, end_at_t, nwf);

	timing_results.resize(nwf, vector<float>(8));

	double* xvals = new double[n_range]; // x values for spline interpolation
	for (int k = 0; k < n_range; k++) xvals[k] = static_cast<double>(k) + .5;
	
	#pragma omp parallel for
	for (int j = 0; j < nwf; j++) {
		double* yvals = Helpers::gety(rundata[j], start_at, end_at); // get range where to search for CFD for timing
		
		// smoothing to suppress noise, will also change timing so use with care!
		if (sigma > 0.) Filters::SmoothArray(yvals, n_range, sigma, smooth_method);

		double* max_ptr = max_element(yvals, yvals + n_range);
		float max_val = *max_ptr;
		int n_max = static_cast<int>(max_ptr - yvals);

		float cf = (cf_r <= 1) ? cf_r * max_val : cf_r;

		int i = 0;
		if (!find_CF_from_start) {
			i = n_max;
			while (i > 0 && yvals[i] > cf) i--;
		}
		else {
			i = 0;
			while (i < n_max && yvals[i] < cf) i++;
			i--;
		}

		float interpol_bin = 0.;
		pair<float, bool> lin_interpol_res = {0, false};
		if (i==0) {
			if (verbose) cout << "WARNING: CFD failed for Ch" << GetCurrentChannel(j) << ", event " << GetCurrentEvent(j) << endl;
		}
		else {
			// do interpolation for cf
			lin_interpol_res = LinearInterpolation(cf, static_cast<float>(i), static_cast<float>(i + 1), yvals[i], yvals[i + 1]);
			// go to center of bin
			interpol_bin = lin_interpol_res.first + .5;

			if (use_spline) { // use spline interpolation with tolerance epsilon*bin_size
				double epsilon = 1e-4;
				double x_low = interpol_bin - .5;
				double x_high = interpol_bin + .5;

				TSpline5* wfspl = 0;
				wfspl = new TSpline5("wf_spline", xvals, yvals, n_range, "b1e1b2e2", 0., 0., 0., 0.);

				// using bisection method: halving search window until cf is less than epsilon bins from spline value
				while (x_high - x_low > epsilon) {
					double x_mid = (x_low + x_high) / 2;
					double f_mid = wfspl->Eval(x_mid);
					if (f_mid == cf) break;

					if (f_mid > cf) x_high = x_mid;
					else x_low = x_mid;
				}
				interpol_bin = (x_low + x_high) / 2;
				delete wfspl;
			}
		}
		timing_results[j][0] = interpol_bin;											// the bin we looked for
		timing_results[j][1] = (interpol_bin + static_cast<float>(start_at)) * SP;		// cfd-time we looked for
		timing_results[j][2] = max_val;													// maximum value
		timing_results[j][3] = n_max;													// bin of maximum
		timing_results[j][4] = cf;														// constant fraction
		timing_results[j][5] = IndexToTime(start_at);									// starting time
		timing_results[j][6] = IndexToTime(end_at);								        // end time
		timing_results[j][7] = static_cast<float>(lin_interpol_res.second);				// flag will be 1 if linear interpolation worked
		delete[] yvals;
	}
	delete[] xvals;
}
/// @example timing_example.cc

/// @brief Skip events where the time difference between two channels is outside of specified range
/// 
/// Skip events where the time difference between channel "first_channel" and channel "second_channel" 
/// is less than "time_diff_min" or more than "time_diff_max" \n 
/// Needs to be called before the charge spectrum etc functions. \n 
/// Baseline correction should be called before this function.
/// 
/// @param first_channel_abs First channel.
/// @param second_channel_abs Second channel.
/// @param time_diff_min Skip events for \f$ \Delta t<t_{diff,min} \f$
/// @param time_diff_max Skip events for \f$ \Delta t>t_{diff,max} \f$
/// @param verbose Set true for extra verbosity.
void ReadRun::SkipEventsTimeDiffCut(int first_channel_abs, int second_channel_abs, double time_diff_min, double time_diff_max, bool verbose) {

	cout << "\n Removing events if the event-wise time difference between the main peaks in ch"
		<< first_channel_abs << " and ch" << second_channel_abs << " is <" << setprecision(2)
		<< time_diff_min << " ns or >" << time_diff_max << " ns" << endl;

	int counter = 0;
	int first_channel = GetChannelIndex(first_channel_abs);
	int second_channel = GetChannelIndex(second_channel_abs);

	// call GetTimingCFD() in case it was not initialized
	if (static_cast<int>(timing_results.size()) == 0) GetTimingCFD();

	// loop through events, calculate timing difference between channels and compare with cuts
	#pragma omp parallel for reduction(+:counter)
	for (int j = 0; j < nevents; j++) {
		if (!SkipEvent(j, first_channel) && !SkipEvent(j, second_channel)) {
			const int ind_first = GetWaveformIndex(j, first_channel);
			const int ind_second = GetWaveformIndex(j, second_channel);
			if (max(ind_first, ind_second) >= static_cast<int>(timing_results.size())) continue;

			float time_diff = timing_results[ind_first][1] - timing_results[ind_second][1];

			if (time_diff < time_diff_min || time_diff > time_diff_max) {
				int currevent = eventnr_storage[j];
				if (verbose) cout << "\nevent:\t" << currevent << "\tchannels:\t" << first_channel_abs << " & " << second_channel_abs << "\ttime diff:\t" << time_diff;
				skip_event[j] = true;
				counter++;
			}
		}
	}
	cout << "\t" << counter << " events will be cut out of " << nevents << endl;
}
/// @example timing_example.cc


/// @brief Skip events above/below individual thresholds per channel
/// 
/// Needs to be called before the charge spectrum etc functions. \n 
/// Baseline correction should be called before this function.
/// 
/// @param thresholds Vector should contain a threshold for each active channel saved in the data, in ascending order (ch0, ch1 ...). 
/// Negative thresholds mean events below threshold will be cut. If only one value is given this value will be used for all active channels.
/// @param rangestart Range start in ns
/// @param rangeend Range end in ns
/// @param verbose Set true for extra verbosity.
void ReadRun::SkipEventsPerChannel(vector<float> thresholds, float rangestart, float rangeend, bool verbose) { // merge with IntegralFilter()?

	if (thresholds.empty()) cout << "\nError: thresholds is empty";
	while (thresholds.size() <= active_channels.size()) thresholds.push_back(thresholds[0]);

	cout << "\n Removing events with individual amplitude threshold per channel:" << endl;
	int counter = 0;
	int n_thrshld = static_cast<int>(thresholds.size());

	int bin_start = CheckBoundsX(TimeToIndex(rangestart));
	int bin_end = CheckBoundsX(TimeToIndex(rangeend));

	#pragma omp parallel for reduction(+:counter)
	for (int j = 0; j < nwf; j++) {
		int current_event = GetCurrentEvent(j);
		int current_channel = GetCurrentChannel(j);
		if (current_event >= 0 && !SkipEvent(current_event, current_channel)) {
			if (current_channel < n_thrshld) {
				const float current_threshold = thresholds[current_channel];
				if (current_threshold == 0.) continue;

				auto max_it = max_element(rundata[j].begin() + bin_start, rundata[j].begin() + bin_end);
				auto min_it = min_element(rundata[j].begin() + bin_start, rundata[j].begin() + bin_end);

				float max_val = (max_it != rundata[j].end()) ? *max_it : 0.;
				float min_val = (min_it != rundata[j].end()) ? *min_it : 0.;

				
				if ((current_threshold > 0	&& max_val > current_threshold) || 
					(current_threshold < 0 && min_val < current_threshold)) {

					int currevent = eventnr_storage[current_event];
					if (verbose) cout << "\nevent:\t" << currevent << "\tchannel:\t" << active_channels[current_channel] << "\tthreshold\t" << current_threshold;
					skip_event[current_event] = true;
					counter++;
				}
			}
		}
	}

	cout << "\t" << counter << " events will be cut out of " << nevents << endl;
}

/// @brief Skip events with threshold on integral
/// 
/// Compare with SkipEventsPerChannel() \n
/// Needs to be called before the charge spectrum etc functions in order to have an effect on them. \n 
/// Baseline correction should be called before this function. \n
/// 
/// @param thresholds Vector should contain a threshold for each active channel saved in the data, in ascending order (ch0, ch1 ...). 
/// Negative thresholds mean events below threshold will be cut. A threshold of 0 means the channel will not be evaluated. 
/// If only one value is given this value will be used for all active channels.
/// @param g_thr  Vector should contain a bool for each active channel. True means events with integrals greater than the threshold will be cut, 
/// false means less than threshold. If only one value is given this value will be used for all active channels.
/// @param windowlow Integration time left to the maximum of the peak.
/// @param windowhi Integration time right to the maximum of the peak.
/// @param start Range for finding maximum.
/// @param end Range for finding maximum.
/// @param use_AND_condition If set true, events skipped by previously applied cuts will be un-skipped if they pass integral criterium (true && false => false, true && true => true).
/// @param verbose Set true for extra verbosity.
void ReadRun::IntegralFilter(vector<float> thresholds, vector<bool> g_thr, float windowlow, float windowhi, float start, float end, bool use_AND_condition , bool verbose) {

	if (thresholds.empty() || g_thr.empty()) cout << "\nERROR: thresholds or g_thr are empty in ReadRun::IntegralFilter().";
	while (thresholds.size() <= active_channels.size()) { thresholds.push_back(thresholds[0]); }
	while (g_thr.size() <= active_channels.size()) { g_thr.push_back(g_thr[0]); }

	cout << "\n\nRemoving events with individual integral threshold per channel:" << endl;
	int counter = 0;
	int n_thresholds = static_cast<int>(thresholds.size());

	#pragma omp parallel for reduction(+:counter)
	for (int j = 0; j < nwf; j++) {
		int currevent_counter = GetCurrentEvent(j);
		if (currevent_counter < 0) continue;
		
		if ((use_AND_condition && skip_event[currevent_counter]) || !skip_event[currevent_counter]) {
			int current_channel = GetCurrentChannel(j);
			if (current_channel >= n_thresholds) continue;
			
			float current_threshold = thresholds[current_channel];
			if (current_threshold == 0.) continue;

			auto integral = GetPeakIntegral(rundata[j], windowlow, windowhi, start, end, current_channel);

			int currevent = eventnr_storage[currevent_counter];
			// skip if above/below thresholds
			if ((g_thr[current_channel] && integral > current_threshold) || (!g_thr[current_channel] && integral < current_threshold)) {
				if (!skip_event[currevent_counter]) {
					if (verbose) cout << "Event:\t" << currevent << "\tchannel:\t" << active_channels[current_channel] << "\tthreshold\t" << thresholds[current_channel] << "\tintegral:\t" << integral << endl;
					skip_event[currevent_counter] = true;
					// go to last wf of current event
					while (GetCurrentEvent(++j) == currevent_counter); 
					j--;
					counter++;
				}
			}
			else if (use_AND_condition) { // don't skip if not
				skip_event[currevent_counter] = false;
				if (verbose) cout << "Event:\t" << currevent << "\tchannel:\t" << active_channels[current_channel] << "\thas been flagged good by integral:\t" << integral << endl;
				counter--;
			}
		}
	}
	cout << counter << " additional events will be cut out of " << nevents << " (" << static_cast<float>(100*counter/nevents) << "%)" << endl;
}
/// @example timing_example.cc
/// @example read_exampledata.cc
/// @example read_exampledata.py

/// @brief Prints a list of all skipped events into the terminal for diagnostics
void ReadRun::PrintSkippedEvents() {
	int counter = 0;
	stringstream buffer;
	for (int j = 0; j < static_cast<int>(skip_event.size()); j++) {
		if (skip_event[j]) {
			int currevent = eventnr_storage[j];
			buffer << "Event:\t" << currevent << endl;
			counter++;
		}
	}
	cout << buffer.str().c_str();
	cout << "Total number of skipped events:\t" << counter << "\tout of:\t" << nevents << endl;
}

/// @brief Sets skip_event flag to false for all events, removing any previous cuts
void ReadRun::UnskipAll() {
	const int skip_event_size = static_cast<int>(skip_event.size());
	for (int j = 0; j < skip_event_size; j++) skip_event[j] = false;
	cout << "\n\nAll event cuts were removed" << endl;
}

/// @brief Check if event should be skipped
/// @param event_index Index of the event
/// @param channel_index Index of the channel (see override for ReadSampic)
bool ReadRun::SkipEvent(int event_index, int channel_index) {
	(void)channel_index; // avoid unused parameter warning
	if (event_index >= static_cast<int>(skip_event.size()) || event_index < 0) return true;
	else return skip_event[event_index];
}

/// @brief Number of good events that are not skipped
/// @param channel Not used here, used in class ReadSampic::ReadSampic
int ReadRun::Nevents_good(int channel_index) {
	int nevents_good = 0;
	for (int i = 0; i < nevents; i++) if (!SkipEvent(i, channel_index)) nevents_good++;
	return nevents_good;
}

// functions for charge spectrum

/// @brief Determine indices for integration window for peaks
/// 
/// Default usage: Find maximum in range ("start", "end") and return bin numbers for [0] the max, 
/// [1] t_max - "windowlow", and [2] t_max + "windowhi". \n
/// If ("start" < 0 || "end" < 0) doesn't return max and integration window is fixed relative to the maximum of the sum of all waveforms: \n
/// t(max(sum_spectrum[channel])) +/- "windowhi"/"windowlow" \n 
/// If ("windowlow" == "start" && "windowhi" == "end") doesn't return max and sets fixed integration
/// window from "start" until "end" for all channels.
/// 
/// @param his Histogram to integrate.
/// @param windowlow Integration time left to the maximum of the peak.
/// @param windowhi Integration time right to the maximum of the peak.
/// @param start Range for finding maximum.
/// @param end Range for finding maximum.
/// @param channel Channel index in case the integration should be around the maximum of the sum of all waveforms
/// @return Array { max, \f$ n_{t,start} \f$ , \f$ n_{t,end} \f$ }
array<int, 3> ReadRun::GetIntWindow(TH1F* his, float windowlow, float windowhi, float start, float end, int channel) {

	int istart, iend;
	array<int, 3> foundindices = {0, 0, 0};

	if (start < 0 || end < 0) {							// fixed integration window relative to maximum of sum spectrum for each channel
		foundindices[1] = his->GetXaxis()->FindBin(his->GetXaxis()->GetBinCenter(maxSumBin[channel]) - windowlow);
		foundindices[2] = his->GetXaxis()->FindBin(his->GetXaxis()->GetBinCenter(maxSumBin[channel]) + windowhi);
	}
	else if (windowlow == start && windowhi == end) {	// fixed integration window for all channels
		foundindices[1] = his->GetXaxis()->FindBin(windowlow);
		foundindices[2] = his->GetXaxis()->FindBin(windowhi);
	}
	else {												// fixed integration window relative to maximum of each individual waveform
		istart = his->GetXaxis()->FindBin(start);
		iend = his->GetXaxis()->FindBin(end);
		foundindices[0] = istart;

		if (istart < 1 || iend > his->GetNbinsX()) {
			cout << "\nError: Start=" << istart << " or end=" << iend << " of GetIntWindow() out of range. Fix integration window." << endl;
		}

		float max_val = -9.e99;
		float curr_val = 0;
		for (int i = istart; i < iend; i++) {
			curr_val = his->GetBinContent(i);
			if (curr_val > max_val) {
				max_val = curr_val;
				foundindices[0] = i;
			}
		}

		foundindices[1] = his->GetXaxis()->FindBin(his->GetXaxis()->GetBinCenter(foundindices[0]) - windowlow);
		foundindices[2] = his->GetXaxis()->FindBin(his->GetXaxis()->GetBinCenter(foundindices[0]) + windowhi);
	}
	return foundindices;
}

/// @brief Determine indices for integration window for peaks
/// 
/// Default usage: Find maximum in range ("start", "end") and return bin numbers for [0] the max, 
/// [1] t_max - "windowlow", and [2] t_max + "windowhi". \n
/// If ("start" < 0 || "end" < 0) doesn't return max and integration window is fixed relative to the maximum of the sum of all waveforms: \n
/// t(max(sum_spectrum[channel])) +/- "windowhi"/"windowlow" \n 
/// If ("windowlow" == "start" && "windowhi" == "end") doesn't return max and sets fixed integration
/// window from "start" until "end" for all channels.
/// 
/// @param waveform Waveform to integrate.
/// @param windowlow Integration time left to the maximum of the peak.
/// @param windowhi Integration time right to the maximum of the peak.
/// @param start Range for finding maximum.
/// @param end Range for finding maximum.
/// @param channel Channel index in case the integration should be around the maximum of the sum of all waveforms
/// @return Array { max, \f$ n_{t,start} \f$ , \f$ n_{t,end} \f$ }
array<int, 3> ReadRun::GetIntWindow(const vector<float>& waveform, float windowlow, float windowhi, float start, float end, int channel) {

	array<int, 3> foundindices = {0, 0, 0};

	if (start < 0 || end < 0) {							// fixed integration window relative to maximum of sum spectrum for each channel
		foundindices[1] = maxSumBin[channel] - TimeToIndex(windowlow);
		foundindices[2] = maxSumBin[channel] - TimeToIndex(windowhi);
	}
	else if (windowlow == start && windowhi == end) {	// fixed integration window for all channels
		foundindices[1] = TimeToIndex(windowlow);
		foundindices[2] = TimeToIndex(windowlow);
	}
	else {												// fixed integration window relative to maximum of each individual waveform
		int istart = TimeToIndex(start);
		int iend = TimeToIndex(end);
		foundindices[0] = istart;

		auto max_it = max_element(waveform.begin() + istart, waveform.begin() + iend);
		foundindices[0] = static_cast<int>(distance(waveform.begin(), max_it));

		float t_max = IndexToTime(foundindices[0]) + SP * 0.5; // time of maximum at bin center
		foundindices[1] = TimeToIndex(t_max - windowlow);
		foundindices[2] = TimeToIndex(t_max + windowhi);
	}
	return foundindices;
}

/// @brief Determine indices for integration window for peaks with bin indices
/// 
/// @param waveform Waveform to integrate.
/// @param windowlow Integration time left to the maximum of the peak in bins.
/// @param windowhi Integration time right to the maximum of the peak in bins.
/// @param start Range for finding maximum in bins.
/// @param end Range for finding maximum in bins.
/// @param channel Channel index in case the integration should be around the maximum of the sum of all waveforms
/// @return Array { max, \f$ n_{t,start} \f$ , \f$ n_{t,end} \f$ }
array<int, 3> ReadRun::GetIntWindow(const vector<float>& waveform, int windowlow, int windowhi, int start, int end, int channel) {

	array<int, 3> foundindices = {0, 0, 0};

	if (start < 0 || end < 0) {							// fixed integration window relative to maximum of sum spectrum for each channel
		foundindices[1] = maxSumBin[channel] - windowlow;
		foundindices[2] = maxSumBin[channel] - windowhi;
	}
	else if (windowlow == start && windowhi == end) {	// fixed integration window for all channels
		foundindices[1] = windowlow;
		foundindices[2] = windowlow;
	}
	else {												// fixed integration window relative to maximum of each individual waveform
		int istart = CheckBoundsX(start);
		int iend = CheckBoundsX(end);
		foundindices[0] = istart;

		auto max_it = max_element(waveform.begin() + istart, waveform.begin() + iend);
		foundindices[0] = static_cast<int>(distance(waveform.begin(), max_it));

		foundindices[1] = CheckBoundsX(foundindices[0] - windowlow);
		foundindices[2] = CheckBoundsX(foundindices[0] + windowhi);
	}
	return foundindices;
}

/// @brief Calculate the integral around a peak with several options explained in GetIntWindow().
/// 
/// Calculated over an integer number of bins. \n
/// If start = end = 0 will return the amplitude of the peak. \n
/// If start < 0 || end < 0 will return the integral over the fixed integration window around 
/// the maximum of the sum spectrum for each channel. \n
/// If windowlow == start && windowhi == end will return the integral over the fixed integration window for all channels.
/// 
/// @param his Histogram with peak.
/// @param windowlow Integration time left to the maximum of the peak.
/// @param windowhi Integration time right to the maximum of the peak.
/// @param start Range for finding maximum.
/// @param end Range for finding maximum.
/// @param channel_index Channel index in case the integration should be around the maximum of the sum of all waveforms.
/// @return Integral/amplitude.
float ReadRun::GetPeakIntegral(TH1F* his, float windowlow, float windowhi, float start, float end, int channel_index) {
	auto windowind = GetIntWindow(his, windowlow, windowhi, start, end, channel_index);	// find integration window
	string integral_option(""); // For amplitude -> unit[mV].
	if (windowind[1] != windowind[2]) integral_option = "width"; // 'width' (bin width) for integral -> unit[mV x ns].
	float integral = his->Integral(windowind[1], windowind[2], integral_option.c_str());
	return integral;
}

/// @brief Calculate the integral around a peak with several options explained in GetIntWindow().
/// 
/// Calculated over an integer number of bins. \n
/// If start = end = 0 will return the amplitude of the peak. \n
/// If start < 0 || end < 0 will return the integral over the fixed integration window around 
/// the maximum of the sum spectrum for each channel. \n
/// If windowlow == start && windowhi == end will return the integral over the fixed integration window for all channels.
/// 
/// @param waveform Waveform with peak.
/// @param windowlow Integration time left to the maximum of the peak.
/// @param windowhi Integration time right to the maximum of the peak.
/// @param start Range for finding maximum.
/// @param end Range for finding maximum.
/// @param channel_index Channel index in case the integration should be around the maximum of the sum of all waveforms.
/// @return Integral/amplitude.
float ReadRun::GetPeakIntegral(const vector<float>& waveform, float windowlow, float windowhi, float start, float end, int channel_index) {
	auto windowind = GetIntWindow(waveform, windowlow, windowhi, start, end, channel_index);	// find integration window
	float integral = 0;
	for (int i = windowind[1]; i <= windowind[2]; ++i) integral += waveform[i];					// sum -> unit[mV].
	if (windowind[1] != windowind[2]) integral *= SP; 											// integral -> unit[mV x ns].
	return integral;
}

/// @brief Plot waveforms of all channels for a given event number and add the determined integration windows to the plot
/// 
/// See GetIntWindow() for explanation of parameters. \n
/// Will also add CFD timing if GetTimingCFD() was called before. 
/// 
/// image html PrintChargeSpectrumWF.png "Waveforms in all channels for a single event. Code in example."
/// 
/// @param windowlow Integrate from "windowlow" ns from max...
/// @param windowhi ...to "windowhi" ns from max.
/// @param start Find max from "start" in ns...
/// @param end ...to "end" in ns.
/// @param eventnr Event number
/// @param ymin Y axis minimum
/// @param ymax Y axis maximum
/// @param xmin X axis maximum
/// @param xmax X axis maximum
void ReadRun::PrintChargeSpectrumWF(float windowlow, float windowhi, float start, float end, int eventnr, float ymin, float ymax, float xmin, float xmax) {

	gStyle->SetOptStat(0);
	TString name(Form("waveforms_event__%05d", eventnr));
	auto intwinc = new TCanvas(name.Data(), name.Data(), 600, 400);
	Helpers::SplitCanvas(intwinc, active_channels, plot_active_channels);
	int event_index = GetEventIndex(eventnr);
	
	int current_canvas = 0;
	for (int i = 0; i < nchannels; i++) {
		if (PlotChannel(i)) {
			intwinc->cd(++current_canvas);

			auto his = Getwf(i, event_index);
			int wf_index = GetWaveformIndex(event_index, i);

			// drawing and formatting
			gPad->SetTopMargin(.01);
			int last_canvas = nchannels;
			if (!plot_active_channels.empty()) last_canvas = static_cast<int>(plot_active_channels.size());
			if (current_canvas == 1 && last_canvas < 4) gPad->SetLeftMargin(.15);
			if (current_canvas % 4 == 0 || current_canvas == last_canvas) gPad->SetRightMargin(.01);
			his->Draw("HIST");
			his->SetStats(0);

			if (wf_index !=  -1) {
				// create lines to indicate the integration window
				auto windowind = GetIntWindow(his, windowlow, windowhi, start, end, i);
				TLine* low = new TLine(his->GetXaxis()->GetBinCenter(windowind[1]), -5, his->GetXaxis()->GetBinCenter(windowind[1]), 10);
				low->SetLineColor(2);
				TLine* hi = new TLine(his->GetXaxis()->GetBinCenter(windowind[2]), -2, his->GetXaxis()->GetBinCenter(windowind[2]), 3);
				hi->SetLineColor(2);
				TLine* zero = new TLine(0, 0, 320, 0); // draw line at x=0 to check if baseline correction worked
				zero->SetLineColor(1);

				low->Draw("same");
				hi->Draw("same");
				zero->Draw("same");			

				// draw baseline and CFD parameters
				if (wf_index < static_cast<int>(baseline_correction_result.size())) {
					TLine* baselinel = new TLine(baseline_correction_result[wf_index][2], -1, baseline_correction_result[wf_index][2], 1);
					baselinel->SetLineColor(6);
					baselinel->SetLineWidth(2);
					TLine* baselineh = new TLine(baseline_correction_result[wf_index][3], -1, baseline_correction_result[wf_index][3], 1);
					baselineh->SetLineColor(6);
					baselineh->SetLineWidth(2);
					TLine* baseline = new TLine(baseline_correction_result[wf_index][2], 0, baseline_correction_result[wf_index][3], 0);
					baseline->SetLineColor(6);
					TLine* correction_value = new TLine(baseline_correction_result[wf_index][2], baseline_correction_result[wf_index][0], baseline_correction_result[wf_index][3], baseline_correction_result[wf_index][0]);
					correction_value->SetLineColor(2);

					baselinel->Draw("same");
					baselineh->Draw("same");
					baseline->Draw("same");
					correction_value->Draw("same");
				}
				
				if (wf_index < static_cast<int>(timing_results.size())) {
					TLine* timing = new TLine(timing_results[wf_index][1], -10, timing_results[wf_index][1], 100);
					timing->SetLineColor(9);
					timing->SetLineWidth(2);
					timing->Draw("same");
				}
			}

			if (ymin != 0. && ymax != 0.) his->GetYaxis()->SetRangeUser(ymin, ymax); // fix y range for better comparison 
			if (xmin != 0. && xmax != 0.) his->GetXaxis()->SetRangeUser(xmin, xmax);
		}
	}
	intwinc->Update();

	root_out->WriteObject(intwinc, name.Data());
}
/// @example timing_example.cc

/// @brief Returns array with the individual "charge"/amplitude for all events of one channel
/// 
/// Missing values evaluate to -999. See SaveChargeLists() and GetIntWindow().
/// 
/// @param channel_index Index of the channel.
/// @param windowlow Integrate from "windowlow" ns from max...
/// @param windowhi ...to "windowhi" ns from max.
/// @param start Find max from "start" in ns...
/// @param end ...to "end" in ns.
/// @param negative_vals If true will save negative values. If false will set negative values to 0.
float* ReadRun::ChargeList(int channel_index, float windowlow, float windowhi, float start, float end, bool negative_vals) {
	float* charge_list = new float[nevents]();

	#pragma omp parallel for
	for (int j = 0; j < nevents; j++) {
		int wf_index = GetWaveformIndex(j, channel_index);
		if (wf_index == -1) {
			charge_list[j] = -999.;
			continue;
		}

		charge_list[j] = GetPeakIntegral(rundata[wf_index], windowlow, windowhi, start, end, channel_index);
		if (!negative_vals && charge_list[j] < 0.) charge_list[j] = 0.;
	}
	return charge_list;
}

/// @brief Saves TGraphs to root file with the individual "charge"/amplitude for all events and all channels
/// 
/// Event with a skip_event flag will be removed. Call before filtering or call UnskipAll() before to get all events. \n
/// See GetIntWindow().
/// 
/// @param windowlow Integrate from "windowlow" ns from max...
/// @param windowhi ...to "windowhi" ns from max.
/// @param start Find max from "start" in ns...
/// @param end ...to "end" in ns.
/// @param negative_vals If true will save negative values. If false will set negative values to 0.
void ReadRun::SaveChargeLists(float windowlow, float windowhi, float start, float end, bool negative_vals) {
	float* event_list = new float[nevents];
	for (int i = 0; i < nevents; i++) event_list[i] = static_cast<float>(i);

	auto charge_list_mg = new TMultiGraph();
	if (windowlow + windowhi > 0.) charge_list_mg->SetTitle("event-wise integrals; Event number; integral [mV#timesns]");
	else charge_list_mg->SetTitle("event-wise amplitudes; Event number; amplitude [mV]");

	for (int i = 0; i < nchannels; i++) {
		if (PlotChannel(i)) {
			TString name(Form("charge_list_ch_%02d", active_channels[i]));
			float* charge_list = ChargeList(i, windowlow, windowhi, start, end, negative_vals);
			TGraph* charge_list_graph = new TGraph(nevents, event_list, charge_list);
			charge_list_graph->SetLineWidth(0);
			charge_list_graph->SetMarkerStyle(2);
			charge_list_graph->SetMarkerColor(Helpers::rcolor(i));
			charge_list_graph->SetTitle(name.Data());

			//remove skipped events
			for (int j = 0; j < nevents; j++) {
				if (skip_event[j]) charge_list_graph->RemovePoint(j);
			}

			charge_list_mg->Add(charge_list_graph);
			root_out->WriteObject(charge_list_graph, name.Data());
			delete[] charge_list;
		}
	}
	root_out->WriteObject(charge_list_mg, "all_charge_lists");
	delete[] event_list;
}

/// @brief Plot correlation of integrals/amplitudes between two channels
/// 
/// See GetIntWindow() and PrintChargeSpectrum() for parameters.
/// 
/// \image html ChargeCorrelation.png "Integrals of the signals of two channels for each event plottet against each other. A clear correlation is visible: On average, if channel 2 records a large signal also channel 3 records a large signal. Around (0,0) the correlation of the pedestals can be seen (no signal in both channels). Code in example." width=75% 
/// 
/// @param windowlow Integrate from "windowlow" ns from max...
/// @param windowhi ...to "windowhi" ns from max.
/// @param start Find max from "start" in ns...
/// @param end ...to "end" in ns. 
/// @param rangestart Plot x & y range start
/// @param rangeend Plot x & y range end
/// @param nbins Number of x & y bins of the histogram
/// @param channel1 First channel number to compare
/// @param channel2 second channel number to compare
/// @param ignore_skipped_events Set true to plot only events which passed filtering, else all events will be plotted
void ReadRun::ChargeCorrelation(float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins, int channel1, int channel2, bool ignore_skipped_events) {
	gStyle->SetOptStat(1111);
	stringstream name;
	name << "charge_correlation_ch" << channel1 << "_ch" << channel2;
	stringstream title;
	if (windowlow + windowhi > 0.) title << ";integral ch" << channel1 << " in mV#timesns;integral ch" << channel2 << " in mV#timesns;Entries";
	else title << ";amplitude ch" << channel1 << " in mV;amplitude ch" << channel2 << " in mV;Entries";

	auto charge_corr_canvas = new TCanvas(name.str().c_str(), "canvas", 600, 600);
	charge_corr_canvas->SetRightMargin(0.15);

	float* charge1 = ChargeList(GetChannelIndex(channel1), windowlow, windowhi, start, end);
	float* charge2 = ChargeList(GetChannelIndex(channel2), windowlow, windowhi, start, end);

	auto charge_corr = new TH2F(name.str().c_str(), title.str().c_str(), nbins, rangestart, rangeend, nbins, rangestart, rangeend);
	
	for (int i = 0; i < nevents; i++) {
		if (charge1[i] == -999. || charge2[i] == -999.) continue; // skip if data is missing for one of the channels
		if (!ignore_skipped_events || !skip_event[i]) charge_corr->Fill(charge1[i], charge2[i]);
	}

	charge_corr->Draw("colz");
	root_out->WriteObject(charge_corr, name.str().c_str());

	charge_corr_canvas->Update();
	charge_corr_canvas->SetGrid();
	// move stat box out of the way (causing problems since May 23?)
	//TPaveStats* stat_box = (TPaveStats*)charge_corr->FindObject("stats"); 
	//stat_box->SetX1NDC(0.6); 
	//stat_box->SetX2NDC(0.85);
	charge_corr->SetStats(0);
	charge_corr_canvas->Modified();
	name << "_c";
	root_out->WriteObject(charge_corr_canvas, name.str().c_str());
	delete[] charge1;
	delete[] charge2;
}
/// @example timing_example.cc

/// @brief Histogram of the "charge" spectrum for one channel
/// 
/// See PrintChargeSpectrum() for parameters.
/// 
/// @return Histogram for one channel.
TH1F* ReadRun::ChargeSpectrum(int channel_index, float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins) {

	TString name(Form("channel__%02d", active_channels[channel_index]));
	auto h1 = new TH1F(name.Data(), name.Data(), nbins, rangestart, rangeend);

	float pedestal = 0;
	float gain = 1;
	if (channel_index < static_cast<int>(PrintChargeSpectrum_cal.size())) {
		if (PrintChargeSpectrum_cal[channel_index][0] != 0) gain = PrintChargeSpectrum_cal[channel_index][0];
		if (PrintChargeSpectrum_cal[channel_index][1] != 1) pedestal = PrintChargeSpectrum_cal[channel_index][1];
	}

	#pragma omp parallel
	{
		auto h1_tmp = new TH1F("", "", nbins, rangestart, rangeend);
		#pragma omp for
		for (int j = 0; j < nevents; j++) {
			if (!SkipEvent(j, channel_index)) {
				float integral_value = GetPeakIntegral(rundata[GetWaveformIndex(j, channel_index)], windowlow, windowhi, start, end, channel_index);
				h1->Fill((integral_value - pedestal) / gain); 
			}
		}
		#pragma omp critical
		h1->Add(h1_tmp);
		delete h1_tmp;
	}
	return h1;
}

/// @brief Plots the "charge" spectrums of all channels
/// 
/// Integrate all pulses in range ("start", "end") from t_max - "windowlow" to t_max + "windowhi" for a given channel 
/// and return the charge histogram with x range ("rangestart", "rangeend") and the number of bins "nbins". \n 
/// It is not really charge, but either amplitude (mV) or integral (mV x ns).
/// See ChargeSpectrum() and GetIntWindow(). \n
/// Can be normalised to the number of photoelectrons by defining the calibration values in PrintChargeSpectrum_cal.
/// 
/// \image html PrintChargeSpectrum.png "Simple example of the integrated signals for a single channel. The resulting integrated signals are fitted with a Landau-Gauss convolution (-> energy deposition of a minimum ionizing particle in a thin absorber). Code in example." width=75%
/// 
/// \image html cosmics-fit-example.png "Integrated signals for a SiPM in blue and a fit with the fit function Fitf for SiPMs missing after-pulses and dark counts in red. The spectrum is not well fit by the simplified model. You can try different fit models in the example cosmics-fit.ipynb." width=50%
/// 
/// @param windowlow Integrate from "windowlow" ns from max...
/// @param windowhi ...to "windowhi" ns from max.
/// @param start Find max from "start" in ns...
/// @param end ...to "end" in ns.
/// @param rangestart Plot x range start
/// @param rangeend Plot x range end
/// @param nbins Number of bins of histogram
/// @param fitrangestart Fit range start
/// @param fitrangeend Fit range end
/// @param max_channel_nr_to_fit Fit only channels with index < "max_channel_nr_to_fit". Set to -1 to skip fitting.
/// @param which_fitf Choose fit function: \n 
/// 0 - do not fit \n
/// 1 - Fitf_langaus: landau gauss convolution for large number of photons \n 
/// 2 - Fitf_biased: if pedestal is biased because of peak finder algorithm \n 
/// 3 - Fitf_full: SiPM fit function with exponential delayed after pulsing \n 
/// 4 - Fitf_PMT_ideal: ideal %PMT fit function \n 
/// 5 - Fitf_PMT: %PMT fit function \n 
/// 6 - Fitf_PMT_pedestal: %PMT fit function with biased pedestal \n 
/// 7 - Fitf_plus_DC: default SiPM fit function + dark count spectrum (for lots of false triggers) \n 
/// else - Fitf: default SiPM fit function \n 
/// @param use_log_y Set all y axes to log scale (for dark count spectra)
void ReadRun::PrintChargeSpectrum(float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins, float fitrangestart, float fitrangeend, int max_channel_nr_to_fit, int which_fitf, bool use_log_y) {
	// print ReadRun::ChargeSpectrum for all channels optimized for SiPM signals
	checkData();
	cout << "Creating charge spectrum ..." << endl;
	
	gStyle->SetOptStat("ne");
	gStyle->SetOptFit(1111);

	if (fitrangestart == 0.)	fitrangestart = rangestart;
	if (fitrangeend == 0.)		fitrangeend = rangeend;

	string ctitle("\"charge\" spectra" + to_string(PrintChargeSpectrum_cnt++));
	auto chargec = new TCanvas(ctitle.c_str(), ctitle.c_str(), 600, 400);
	Helpers::SplitCanvas(chargec, active_channels, plot_active_channels);
	int current_canvas = 0;

	float default_rangestart = -2000;
	float default_rangeend = 30000;
	if (default_rangestart > rangestart) default_rangestart = rangestart;
	if (default_rangeend < rangeend) default_rangeend = rangeend;
	int default_nbins = static_cast<int>((default_rangeend - default_rangestart) * nbins / (rangeend - rangestart));

	// create histograms in parallel
	map<int, shared_future<TH1F*>> h1_future;
    for (int i = 0; i < nchannels; ++i) {
        if (PlotChannel(i)) {
            h1_future[i] = async(launch::async, [this, i, windowlow, windowhi, start, end, default_rangestart, default_rangeend, default_nbins]() {
                return ChargeSpectrum(i, windowlow, windowhi, start, end, default_rangestart, default_rangeend, default_nbins);
            });
        }
    }

	for (int i = 0; i < nchannels; i++) {
		if (PlotChannel(i)) {
			chargec->cd(++current_canvas);

			auto his = h1_future[i].get();
			his->GetYaxis()->SetTitle("#Entries");
			if (windowlow + windowhi > 0.) his->GetXaxis()->SetTitle("Integral in mV#timesns");
			else his->GetXaxis()->SetTitle("Amplitude in mV");
			
			if (i < static_cast<int>(PrintChargeSpectrum_cal.size()) && PrintChargeSpectrum_cal[i][0] != 1) {
				cout << "Charge spectrum for channel index " << i << " will be normalized using a gain of "
					<< PrintChargeSpectrum_cal[i][0] << " and a pedestal value of " << PrintChargeSpectrum_cal[i][1] << endl;
				his->GetXaxis()->SetTitle("Number of photoelectrons");
			}

			//store the mean integral of each channel --> used for correction factors of phi_ew analysis
			mean_integral.push_back(his->GetMean());

			//fitting
			if (i < max_channel_nr_to_fit) {
				if (which_fitf == 0) {}
				else if (which_fitf == 1) { // landau gauss convolution for large number of photons
					Fitf_langaus fitf;
					int n_par = 4;
					TF1* f = new TF1("fitf_langaus", fitf, fitrangestart, fitrangeend, n_par); f->SetLineColor(3); f->SetNpx(1000);

					f->SetParName(0, "Width");				f->SetParameter(0, 35);
					f->SetParName(1, "MPV");				f->SetParameter(1, 1000);
					f->SetParName(2, "Area");			    f->SetParameter(2, 10000);
					f->SetParName(3, "#sigma_{Gauss}");		f->SetParameter(3, 100);

					if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < min(4, static_cast<int>(PrintChargeSpectrum_pars.size())); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

					if (i < max_channel_nr_to_fit) {
						cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
						TFitResultPtr fresults = his->Fit(f, "LRS");
						fit_results.push_back(fresults);
					}
				}
				else if (which_fitf == 2) { // if pedestal is biased because of peak finder algorithm
					Fitf_biased fitf;
					int n_par = 9;
					TF1* f = new TF1("fitf_biased", fitf, fitrangestart, fitrangeend, n_par); f->SetLineColor(3); f->SetNpx(1000);

					f->SetParName(0, "N_{0}");				f->SetParameter(0, his->Integral());
					f->SetParName(1, "#mu");				f->SetParameter(1, 2.);
					f->SetParName(2, "#lambda");			f->SetParameter(2, .04);
					f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 2.1);
					f->SetParName(4, "#sigma_{1}");			f->SetParameter(4, 3.4); //f->SetParLimits(4, 1.e-9, 1.e3);
					f->SetParName(5, "Gain");				f->SetParameter(5, 30.); //f->FixParameter(5, 10.);
					f->SetParName(6, "Pedestal");			f->SetParameter(6, 2.);
					f->SetParName(7, "norm_{0}");			f->SetParameter(7, 0.7);
					f->SetParName(8, "x_{0}");				f->SetParameter(8, 5.);

					if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < min(n_par, static_cast<int>(PrintChargeSpectrum_pars.size())); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

					if (i < max_channel_nr_to_fit) {
						cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
						TFitResultPtr fresults = his->Fit(f, "LRS");
						fit_results.push_back(fresults);
					}

					// get number of excess events in the pedestal in the fit region. To get the absolute number of excess events the full pedestal needs to be inside of the fit range (fitrangestart, fitrangeend)
					//double excessEventsInPedestal = f->Integral(fitrangestart, fitrangeend)/.3125;
					//f->SetParameter(7, 1.);
					//excessEventsInPedestal -= f->Integral(fitrangestart, fitrangeend)/.3125;
					//cout << "\nNumber of excess events in the pedestal within the fit range:\t" << excessEventsInPedestal << "\n\n";
				}
				else if (which_fitf == 3) { // SiPM fit function with exponential delayed after pulsing
					Fitf_full fitf;
					int n_par = 9;
					TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, n_par); f->SetLineColor(3); f->SetNpx(1000);

					f->SetParName(0, "N_{0}");				f->SetParameter(0, his->Integral());
					f->SetParName(1, "#mu");				f->SetParameter(1, 2.);
					f->SetParName(2, "#lambda");			f->SetParameter(2, .04);
					f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 2.1);
					f->SetParName(4, "#sigma_{1}");			f->SetParameter(4, 3.4);
					f->SetParName(5, "Gain");				f->SetParameter(5, 30.); //f->FixParameter(5, 40.);
					f->SetParName(6, "Pedestal");			f->SetParameter(6, 2.);
					f->SetParName(7, "#alpha");				f->SetParameter(7, .1); //f->FixParameter(7, .2);
					f->SetParName(8, "#beta");				f->SetParameter(8, 80.); //f->FixParameter(8, 80);

					if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < min(n_par, static_cast<int>(PrintChargeSpectrum_pars.size())); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

					if (i < max_channel_nr_to_fit) {
						cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
						TFitResultPtr fresults = his->Fit(f, "LRS");
						fit_results.push_back(fresults);
					}
				}
				else if (which_fitf == 4) { // ideal PMT fit function
					Fitf_PMT_ideal fitf;
					int n_par = 4;
					TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, n_par); f->SetLineColor(3); f->SetNpx(1000);

					f->SetParName(0, "N_{0}");				f->SetParameter(0, his->Integral());
					f->SetParName(1, "#mu");				f->SetParameter(1, 1.);
					f->SetParName(2, "#sigma");				f->SetParameter(2, 5.);
					f->SetParName(3, "gain");				f->SetParameter(3, 10.);

					if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < min(n_par, static_cast<int>(PrintChargeSpectrum_pars.size())); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

					if (i < max_channel_nr_to_fit) {
						cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
						TFitResultPtr fresults = his->Fit(f, "LRS");
						fit_results.push_back(fresults);
					}
				}
				else if (which_fitf == 5) { // PMT fit function
					Fitf_PMT fitf;
					int n_par = 8;
					TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, n_par); f->SetLineColor(3); f->SetNpx(1000);

					f->SetParName(0, "N_{0}");				f->SetParameter(0, his->Integral());
					f->SetParName(1, "w");					f->SetParameter(1, .05);	f->SetParLimits(1, 1.e-99, 4.e-1); //probability for type II BG
					f->SetParName(2, "#alpha");				f->SetParameter(2, .05);	f->SetParLimits(2, 1.e-99, 5.e-2); //coefficient of exponential decrease of typ II BG
					f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 5.);		f->SetParLimits(3, 1.e-9, 1.e3);
					f->SetParName(4, "Q_{0}");				f->SetParameter(4, 0.);
					f->SetParName(5, "#mu");				f->SetParameter(5, 1.);
					f->SetParName(6, "#sigma_{1}");			f->SetParameter(6, 5.);		f->SetParLimits(6, 1.e-9, 1.e3);
					f->SetParName(7, "Q_{1}");				f->SetParameter(7, 10.);

					if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < min(n_par, static_cast<int>(PrintChargeSpectrum_pars.size())); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

					if (i < max_channel_nr_to_fit) {
						cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
						TFitResultPtr fresults = his->Fit(f, "LRS");
						fit_results.push_back(fresults);
					}
				}
				else if (which_fitf == 6) { // PMT fit function with biased pedestal
					Fitf_PMT_pedestal fitf;
					int n_par = 9;
					TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, n_par); f->SetLineColor(3); f->SetNpx(1000);

					f->SetParName(0, "A");					f->SetParameter(0, his->Integral());
					f->SetParName(1, "w");					f->SetParameter(1, .05);	f->SetParLimits(1, 1.e-9, 4.e-1);	//probability for type II BG
					f->SetParName(2, "#alpha");				f->SetParameter(2, .05);	f->SetParLimits(2, 1.e-9, 5.e-2);	//coefficient of exponential decrease of typ II BG
					f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 5.);		f->SetParLimits(3, 1.e-9, 1.e3);
					f->SetParName(4, "Q_{0}");				f->SetParameter(4, 0.);		f->SetParLimits(4, -1.e2, 1.e2);
					f->SetParName(5, "#mu");				f->SetParameter(5, 1.);		f->SetParLimits(5, 1.e-9, 1.e2);
					f->SetParName(6, "#sigma_{1}");			f->SetParameter(6, 5.);		f->SetParLimits(6, 1.e-9, 1.e3);
					f->SetParName(7, "Q_{1}");				f->SetParameter(7, 10.);	f->SetParLimits(7, 1.e-9, 1.e9);
					f->SetParName(8, "A_{0}");				f->SetParameter(8, 1.);		f->SetParLimits(8, 1.e-9, 1.e1);

					if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < min(n_par, static_cast<int>(PrintChargeSpectrum_pars.size())); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

					if (i < max_channel_nr_to_fit) {
						cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
						TFitResultPtr fresults = his->Fit(f, "LRS");
						fit_results.push_back(fresults);
					}
				}
				else if (which_fitf == 7) { // default SiPM fit function + dark count spectrum (for lots of false triggers)
					Fitf_plus_DC fitf;
					int n_par = 9;
					TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, n_par); f->SetLineColor(3); f->SetNpx(1000);

					f->SetParName(0, "A");					f->SetParameter(0, his->Integral());
					f->SetParName(1, "w");					f->SetParameter(1, .05);	f->SetParLimits(1, 1.e-9, 4.e-1);	//probability for type II BG
					f->SetParName(2, "#alpha");				f->SetParameter(2, .05);	f->SetParLimits(2, 1.e-9, 5.e-2);	//coefficient of exponential decrease of typ II BG
					f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 5.);		f->SetParLimits(3, 1.e-9, 1.e3);
					f->SetParName(4, "Q_{0}");				f->SetParameter(4, 0.);		f->SetParLimits(4, -1.e2, 1.e2);
					f->SetParName(5, "#mu");				f->SetParameter(5, 1.);		f->SetParLimits(5, 1.e-9, 1.e2);
					f->SetParName(6, "#sigma_{1}");			f->SetParameter(6, 5.);		f->SetParLimits(6, 1.e-9, 1.e3);
					f->SetParName(7, "#mu_darkcount");		f->SetParameter(7, .1);		f->SetParLimits(7, 1.e-9, 1.);
					f->SetParName(8, "N_{0}_darkcount");	f->SetParameter(8, .05);	f->SetParLimits(8, 1.e-9, .3);

					if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < min(n_par, static_cast<int>(PrintChargeSpectrum_pars.size())); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

					if (i < max_channel_nr_to_fit) {
						cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
						TFitResultPtr fresults = his->Fit(f, "LRS");
						fit_results.push_back(fresults);
					}
				}
				else { // default SiPM fit function
					Fitf fitf;
					int n_par = 7;
					TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, n_par); f->SetLineColor(3); f->SetNpx(1000);

					f->SetParName(0, "N_{0}");				f->SetParameter(0, his->Integral());
					f->SetParName(1, "#mu");				f->SetParameter(1, 2.);
					f->SetParName(2, "#lambda");			f->SetParameter(2, .04);
					f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 2.1);
					f->SetParName(4, "#sigma_{1}");			f->SetParameter(4, 3.4);
					f->SetParName(5, "Gain");				f->SetParameter(5, 30.); //f->FixParameter(5, 40.);
					f->SetParName(6, "Pedestal");			f->SetParameter(6, 2.);

					if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < min(n_par, static_cast<int>(PrintChargeSpectrum_pars.size())); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

					if (i < max_channel_nr_to_fit) {
						cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
						TFitResultPtr fresults = his->Fit(f, "LRS");
						fit_results.push_back(fresults);
					}
				}
			}
			TString name(Form("ChargeSpectrum channel_%02d_%d", active_channels[i], PrintChargeSpectrum_cnt));
			root_out->WriteObject(his, name.Data());
			his->Draw();
			if (use_log_y) gPad->SetLogy();
		}
	}

	Helpers::SetRangeCanvas(chargec, rangestart, rangeend);
	root_out->WriteObject(chargec, ("ChargeSpectra" + to_string(PrintChargeSpectrum_cnt)).c_str());
}
/// @example timing_example.cc
/// @example read_exampledata.cc
/// @example read_exampledata.py

/// @brief Calculate (SiPM) dark count rate
/// 
/// See PrintChargeSpectrum() for parameters.
/// 
/// @param threshold 1.5 photoelectron threshold
void ReadRun::PrintDCR(float windowlow, float windowhi, float rangestart, float rangeend, double threshold) {

	string unit(" mV");
	if (windowlow + windowhi > 0.) unit = " mV*ns";

	for (int i = 0; i < nchannels; i++) {
		if (PlotChannel(i)) {
			auto his = ChargeSpectrum(i, windowlow, windowhi, rangestart, rangeend, rangestart, rangeend, 500);

			stringstream lonamerate;
			lonamerate << "<0.5 pe=" << threshold << unit << " -> " << his->Integral(his->GetXaxis()->FindBin(rangestart), his->GetXaxis()->FindBin(threshold)) / his->GetEntries() / (1.e-3 * (rangeend - rangestart)) << " MHz";
			cout << "\n" << lonamerate.str().c_str() << endl;

			stringstream hinamerate;
			hinamerate << ">0.5 pe=" << threshold << unit << " -> " << his->Integral(his->GetXaxis()->FindBin(threshold) + 1, his->GetXaxis()->FindBin(rangeend)) / his->GetEntries() / (1.e-3 * (rangeend - rangestart)) << " MHz";
			cout << "\n" << hinamerate.str().c_str() << endl;
		}
	}
}

/// @brief Time distribution of maximum, CFD, or 10% - 90% rise time in a certain time window
/// 
/// See PrintTimeDist() for parameters.
/// 
TH1F* ReadRun::TimeDist(int channel_index, float from, float to, float rangestart, float rangeend, int nbins, int which, float cf_r) {
	int from_n = TimeToIndex(from);
	int to_n = TimeToIndex(to);

	TString name(Form("timedist_ch%02d", active_channels[channel_index]));
	auto h1 = new TH1F(name.Data(), name.Data(), nbins, rangestart, rangeend);
	
	for (int j = 0; j < nevents; j++) {
		if (!SkipEvent(j, channel_index)) {
			int wf_index = GetWaveformIndex(j, channel_index);
			int max_n = GetIntWindow(rundata[wf_index], 0, 0, from_n, to_n)[0];
			float max_val = rundata[wf_index][max_n];

			if (which == 0) { // time of maximum 
				h1->Fill(max_val);
			}
			else if (which == 1) { // time of 50% CFD
				do {
					max_n--;
				} while (rundata[wf_index][max_n] >= cf_r * max_val && max_n >= from_n);
				max_n++;

				h1->Fill(LinearInterpolation(cf_r * max_val, IndexToTime(max_n - 1), IndexToTime(max_n), rundata[wf_index][max_n - 1], rundata[wf_index][max_n]).first);
			}
			else { // 10%-90% rise time
				// search backwards from maximum
				int n10 = -1;
				int n90 = -1;
				do {
					max_n--;
					if (n10 == -1 && rundata[wf_index][max_n] >= .1 * max_val && rundata[wf_index][max_n - 1] <= .1 * max_val) n10 = max_n;
					if (n90 == -1 && rundata[wf_index][max_n] >= .9 * max_val && rundata[wf_index][max_n - 1] <= .9 * max_val) n90 = max_n;
				} while (rundata[wf_index][max_n] <= max_val && max_n > from_n);

				float t10 = LinearInterpolation(.1 * max_val, IndexToTime(n10 - 1), IndexToTime(n10), IndexToTime(n10 - 1), IndexToTime(n10)).first;
				float t90 = LinearInterpolation(.9 * max_val, IndexToTime(n90 - 1), IndexToTime(n90), IndexToTime(n90 - 1), IndexToTime(n90)).first;

				h1->Fill(t90 - t10);
			}
		}
	}
	if (which == 1) h1->Fit("gaus", "WWM", "same");
	return h1;
}

/// @brief Time distribution of maximum, CFD, or 10% - 90% rise time in a certain time window
/// 
/// Find peak time for a given channel in time window ["from", "to"] and return the peak time histogram 
/// with x range ["rangestart", "rangeend"] and the number of bins "nbins". \n 
/// Plots TimeDist() for all channels. \n \n  
/// 
/// For CFD it is advised to use GetTimingCFD() with Print_GetTimingCFD() instead of this function.
/// 
/// @param from Start of time interval in ns
/// @param to End of time interval in ns
/// @param rangestart Start of x range of histogram
/// @param rangeend End of x range of histogram
/// @param nbins Number of bins of histogram
/// @param which Options: \n 
/// 0 - Gives time of maximum. \n 
/// 1 - Gives constant fraction discrimination with fraction "cf_r" of maximum, searching backwards from the maximum until "from". \n 
/// else - gives the 10% - 90% rise time.
/// @param cf_r Fraction of max for CFD.
void ReadRun::PrintTimeDist(float from, float to, float rangestart, float rangeend, int nbins, int which, float cf_r) {
	gStyle->SetOptStat(1111); // 11 is title + entries

	auto time_dist_c = new TCanvas("timing of maximum", "timing of maximum", 600, 400);
	Helpers::SplitCanvas(time_dist_c, active_channels, plot_active_channels);

	int current_canvas = 0;

	for (int i = 0; i < nchannels; i++) {
		if (PlotChannel(i)) {
			time_dist_c->cd(++current_canvas);

			auto his = TimeDist(i, from, to, rangestart, rangeend, nbins, which, cf_r);
			his->GetYaxis()->SetTitle("#Entries");
			his->GetXaxis()->SetTitle("time [ns]");
			his->Draw();
			stringstream name; name << "t_{max} for " << from << "<t<" << to << " ns";
			his->SetTitle(name.str().c_str());

			TString name_save(Form("TimeDist channel_%02d", active_channels[i]));
			root_out->WriteObject(his, name_save.Data());
		}
	}

	time_dist_c->Update();
	root_out->WriteObject(time_dist_c, "TimeDist");
}
/// @example read_exampledata.cc
/// @example read_exampledata.py

/// @brief Finds maximum amplitude for a given channel in time window ["from", "to"] and creates 3d map of waveforms ordered by maxima on z axis.
///
/// Use PrintMaxDist() to plot all channels. \n
/// Use only for small datasets (<10k waveforms) as it will contain all individual points of all waveforms.
///
/// @return TGraph2D of all histograms ordered by maximum amplitude
TGraph2D* ReadRun::MaxDist(int channel_index, float from, float to) {
	// find maximum amplitude for a given channel in time window [from, to] and return 3d histogram with the number of bins nbinsy,z

	TString name(Form("maxdist_ch%02d", active_channels[channel_index]));
	TGraph2D* g3d = new TGraph2D((binNumber + 2) * nevents);
	g3d->SetTitle("waveforms; t [ns]; max. amplitude [mv]; amplitude [mV]");
	g3d->SetMarkerStyle(7);
	double* xvals = getx<double>();

	for (int j = 0; j < nevents; j++) {
		if (!SkipEvent(j, channel_index)) {
			int wf_index = GetWaveformIndex(j, channel_index);

			int bin_from = TimeToIndex(from);
			int bin_to = TimeToIndex(to);
			double max_val = rundata[wf_index][bin_from];
			for (int i = bin_from; i <= bin_to; ++i) {
				double val = rundata[wf_index][i];
				if (val > max_val) max_val = val;
			}
			
			for (int i = 0; i < binNumber; i++) g3d->SetPoint(j * binNumber + i, xvals[i], max_val, rundata[wf_index][i]);
		}
	}
	delete[] xvals;
	root_out->WriteObject(g3d, name.Data());
	return g3d;
}


/// @brief Finds maximum amplitude for a given channel in time window ["from", "to"] and creates 3d map of waveforms ordered by maxima on z axis.
/// 
/// Prints MaxDist() for all channels. \n
/// Use only for small datasets (<100k waveforms) as it will contain all individual points of all waveforms.
///
/// @param from From
/// @param to To 
void ReadRun::PrintMaxDist(float from, float to) {

	auto max_dist_c = new TCanvas("wf grouped by maximum", "wf grouped by maximum", 600, 400);
	Helpers::SplitCanvas(max_dist_c, active_channels, plot_active_channels);

	int current_canvas = 0;

	for (int i = 0; i < nchannels; i++) {
		if (PlotChannel(i)) {
			max_dist_c->cd(++current_canvas);
			auto g3d = MaxDist(i, from, to);
			g3d->Draw("AP");
		}
	}
	max_dist_c->Update();
	root_out->WriteObject(max_dist_c, "MaxDist");
}

/// @brief Plot results of GetTimingCFD()
/// 
/// See Print_GetTimingCFD() for parameters.
/// 
/// @return Timing histogram for one channel
TH1F* ReadRun::His_GetTimingCFD(int channel_index, float rangestart, float rangeend, int nbins) {

	if (nbins == -999) nbins = static_cast<int>((rangeend - rangestart) / SP);

	TString name(Form("GetTimingCFD_ch%02d", active_channels[channel_index]));
	auto his = new TH1F(name.Data(), name.Data(), nbins, rangestart, rangeend);
	for (int j = 0; j < nevents; j++) {
		if (!SkipEvent(j, channel_index)) his->Fill(timing_results[GetWaveformIndex(j, channel_index)][1]);
	}
	return his;
}

/// @brief Plot results of GetTimingCFD()
/// 
/// \image html Print_GetTimingCFD.png "Beginning of the signals for all good events determined with constant fraction discrimination. The red lines are gauss functions fitted to the distrutions. Code in example." width=75%
/// 
/// @param rangestart Start of x range for plot in ns.
/// @param rangeend End of x range for plot in ns.
/// @param do_fit If 1: fits a gaussian. \n
/// Else do not fit. \n 
/// Fit results per channel are stored in ReadRun::timing_fit_results.
/// @param nbins Number of bins for histogram. Will use 320 MHz sampling rate for binning if nbins = -999.
/// @param fitoption ROOT fit option, default is "S". See Fit(): https://root.cern/doc/master/classTH1.html
/// @param set_errors Assign errors to the bins. Will assign errors of 1 to empty bins and \f$ \sqrt(N) \f$ if they are not empty. 
/// Can improve the chi^2 fit.
void ReadRun::Print_GetTimingCFD(float rangestart, float rangeend, int do_fit, int nbins, string fitoption, bool set_errors) {

	// call GetTimingCFD() in case it was not initialized
	if (static_cast<int>(timing_results.size()) == 0) GetTimingCFD();

	gStyle->SetOptStat(1111);
	gStyle->SetOptFit(111);

	auto timing_cfd_c = new TCanvas("timing of cfd", "timing of cfd", 600, 400);
	Helpers::SplitCanvas(timing_cfd_c, active_channels, plot_active_channels);
	int current_canvas = 0;

	for (int i = 0; i < nchannels; i++) {
		if (PlotChannel(i)) {
			timing_cfd_c->cd(++current_canvas);

			auto his = His_GetTimingCFD(i, rangestart, rangeend, nbins);
			his->GetYaxis()->SetTitle("#Entries");
			his->GetXaxis()->SetTitle("time [ns]");

			if (set_errors) {
				int end = his->GetNbinsX();
				for (int k = 1; k <= end; k++) {
					if (his->GetBinContent(k) < 2) his->SetBinError(k, 1);
					else his->SetBinError(k, sqrt(his->GetBinContent(k)));
				}
			}

			his->Draw();

			if (do_fit == 1) {
				TFitResultPtr fresults = his->Fit("gaus", fitoption.c_str(), "same");
				timing_fit_results.push_back(fresults);
			}

			TString name_save(Form("Timing_cfd_channel_%02d", active_channels[i]));
			root_out->WriteObject(his, name_save.Data());
		}
	}

	timing_cfd_c->Update();
	root_out->WriteObject(timing_cfd_c, "TimingCFD");
}
/// @example timing_example.cc

/// @brief Plot timing difference between the mean timings of two channel ranges
/// 
/// See Print_GetTimingCFD_diff() for parameters.
/// 
/// @return Histogram with event-wise timing differences between two channel ranges
TH1F* ReadRun::His_GetTimingCFD_diff(vector<int> channels1, vector<int> channels2, float rangestart, float rangeend, int nbins) {

	if (nbins == -999) nbins = static_cast<int>((rangeend - rangestart) / SP);

	stringstream name;
	name << "GetTimingCFD_diff <";

	// find channel indices and assemble title
	int counter = 0;
	for (int& entry : channels2) {
		if (counter > 0) name << "&";
		auto chin2 = find(active_channels.begin(), active_channels.end(), entry);
		if (chin2 != active_channels.end()) {
			name << "ch" << entry;
			entry = chin2 - active_channels.begin();
		}
		else cout << "\n\n ERROR: channels2 = " << entry << " does not exist in data. Check parameters for Print_GetTimingCFD_diff()\n\n";
		counter++;
	}
	name << ">-<";
	counter = 0;

	for (int& entry : channels1) {
		if (counter > 0) name << "&";
		auto chin1 = find(active_channels.begin(), active_channels.end(), entry);
		if (chin1 != active_channels.end()) {
			name << "ch" << entry;
			entry = chin1 - active_channels.begin();
		}
		else cout << "\n\n ERROR: channels1 = " << entry << " does not exist in data. Check parameters for Print_GetTimingCFD_diff()\n\n";
		counter++;
	}
	name << ">";

	// fill histogram
	auto his = new TH1F(name.str().c_str(), name.str().c_str(), nbins, rangestart, rangeend);
	for (int j = 0; j < nevents; j++) {
		if (!SkipEvent(j)) {
			float mean1 = 0., mean2 = 0., cnt1 = 0., cnt2 = 0.;
			for (int i : channels1) {
				int wf_index = GetWaveformIndex(j, i);
				if (wf_index >= 0) {
					mean1 += timing_results[GetWaveformIndex(j, i)][1];
					cnt1 += 1.;
				}
			}
			for (int i : channels2) {
				int wf_index = GetWaveformIndex(j, i);
				if (wf_index >= 0) {
					mean2 += timing_results[GetWaveformIndex(j, i)][1];
					cnt2 += 1.;
				}
			}

			if (cnt1 != 0. && cnt2 !=0.) his->Fill(mean2 / cnt2 - mean1 / cnt1);
		}
	}

	return his;
}

/// @brief Plot timing difference between the mean timings of two channel ranges
/// 
/// Plots the difference between the peak times between the mean times of two ranges of channels for each event. \n
/// It calculates \f$ \Delta t = <t_{second,i}> - <t_{first,i}> \f$ . \n \n \n
/// 
/// The vectors of channels to compare are added with curly brackets:
/// > mymeas.Print_GetTimingCFD_diff({ 26, 14 }, { 19 }, 0, 20, 2, 200); \n
/// would plot \f$ \Delta t = t_{ch19} - (t_{ch26} + t_{ch14})/2 \f$ from 0 ns to 20 ns with 200 bins (100 ps bin width). 
/// Another example is given in the plot below.
/// 
/// \image html Print_GetTimingCFD_diff.png "Event-wise time differences of the start of the signals of two channels. Code in example." width=75%
/// 
/// @param channels1 Vector of first channel numbers (wavecatcher channel numbers). 
/// @param channels2 Vector of second channel numbers to compare. 
/// @param rangestart Start of x range for plot in ns.
/// @param rangeend End of x range for plot in ns.
/// @param do_fit If 1: Fit a gaussian. \n
/// If 2: Fit a gaussian-exponential convolution (Fitf_exp_gauss) to account for different arrival times of photons due to different 
/// possible light paths in the scintillator/light guide \n
/// and/or delay due to self-absorption and reemission of photons in the scintillator. \n
/// To be used for long light paths in the scintillator. See https://doi.org/10.1016/S0029-554X(79)90170-8 . \n
/// This option only works for sufficient asymmetry \f$\tau > \sigma/2\f$. 
/// Otherwise, the exponential decay time becomes too small to be fitted. \n 
/// If the asymmetry is too small (skewness<0.15) option 1 will be used by default.\n
/// If 3: Fits the sum of two gaussians where the second gauss serves as a rough background estimate. 
/// Background means events that should have been filtered out. \n
/// Else: Do not fit. \n 
/// @param nbins Number of bins for histogram.
/// @param fitrangestart Start of fitting range.
/// @param fitrangeend End of fitting range.
/// @param fitoption ROOT fitting option. Default is "RS" (chi^2). You can try to use the likelihood method with "LRS" if the data is very clean. 
/// Outliers will influence the results for the likelihood method so it is advisable to limit the fit range to exclude outliers for "LRS".
/// @param set_errors Assign errors to the bins. Will assign errors of 1 to empty bins and \f$ \sqrt(N) \f$ if they are not empty. 
/// Can improve the chi^2 fit.
void ReadRun::Print_GetTimingCFD_diff(vector<int> channels1, vector<int> channels2, float rangestart, float rangeend, int do_fit, int nbins, float fitrangestart, float fitrangeend, string fitoption, bool set_errors) {

	// call GetTimingCFD() in case it was not initialized
	if (static_cast<int>(timing_results.size()) == 0) GetTimingCFD();

	if (fitrangestart == -999) {
		fitrangestart = rangestart;
		fitrangeend = rangeend;
	}

	//gStyle->SetOptStat(1111);
	gStyle->SetOptFit(111);

	auto timing_cfd_d_c = new TCanvas("timing of cfd diff", "timing of cfd diff", 600, 400);

	auto his = His_GetTimingCFD_diff(channels1, channels2, rangestart, rangeend, nbins);
	his->GetYaxis()->SetTitle("#Entries");
	his->GetXaxis()->SetTitle("time [ns]");

	if (set_errors) {
		int end = his->GetNbinsX();
		for (int i = 1; i <= end; i++) {
			if (his->GetBinContent(i) < 2) his->SetBinError(i, 1);
			else his->SetBinError(i, sqrt(his->GetBinContent(i)));
		}
	}

	his->Draw();

	double skewness = his->GetSkewness();
	
	if (do_fit == 1 || (do_fit == 2 && abs(skewness) < .15)) {
		// gauss (default)
		TFitResultPtr fresults = his->Fit("gaus", fitoption.c_str(), "same", fitrangestart, fitrangeend);
		timing_fit_results.push_back(fresults);
		if (do_fit == 2) cout << "\nWARNING: Print_GetTimingCFD_diff\nFITTING GAUSS INSTEAD OF GAUSS x EXP CONVOLUTION BC SYMMETRY" << endl;
	}
	else if (do_fit == 2) {
		// gauss x exp convolution (effective delay from random light path and/or self-absorption and reemission)
		Fitf_exp_gauss fitf_exp_gauss;
		auto expgconv = new TF1("exp x gauss convolution", fitf_exp_gauss, fitrangestart, fitrangeend, 4);
		expgconv->SetNpx(5000);

		// this parameter describes the sigma from different light paths 
		// and/or the effective decay time constant for self-absorption and reemission
		expgconv->SetParName(0, "#tau_{eff}");		expgconv->SetParameter(0, skewness);
		if (skewness > 0) expgconv->SetParLimits(0, .15, 5.);
		else expgconv->SetParLimits(0, -5., -.15);
		//expgconv->FixParameter(0, 1.55);

		expgconv->SetParName(1, "#sigma_{gaus}");		expgconv->SetParameter(1, his->GetStdDev());
		expgconv->SetParLimits(1, 1e-1, 7.);	//expgconv->FixParameter(1, .7);

		expgconv->SetParName(2, "t_{0}");		expgconv->SetParameter(2, his->GetMean());
		expgconv->SetParLimits(2, fitrangestart, fitrangeend);	//expgconv->FixParameter(2, 6.6);

		expgconv->SetParName(3, "norm");		expgconv->SetParameter(3, his->Integral("width"));
		expgconv->SetParLimits(3, 1., 1e8);		//expgconv->FixParameter(3, 105.5);

		TFitResultPtr fresults = his->Fit(expgconv, "SR", "same");
		timing_fit_results.push_back(fresults);

		// for the phi_ew-analysis: print out the time value of the maximum of the best fit --> used to determine timing cuts
		float t_of_maximum = expgconv->GetMaximumX(-5, 5);
		cout << "Maximum of the fit is at t=" << t_of_maximum << " ns and the ";

		double max_val = expgconv->GetMaximum();
		double fwhm_x1 = expgconv->GetX(max_val / 2, fitrangestart, fitrangeend);
		double fwhm_x2 = expgconv->GetX(max_val / 2, fwhm_x1 + 1e-3, fitrangeend);
		double fwhm = fwhm_x2 - fwhm_x1;
		auto fwhm_line = new TLine(fwhm_x1, max_val/2, fwhm_x2, max_val/2);
		fwhm_line->SetLineColor(2); fwhm_line->SetLineWidth(2);
		fwhm_line->Draw("same");
		cout << "FWHM=" << fwhm << " ns" << endl;

		// TLatex l;
		// l.SetTextSize(0.025);
   		// l.DrawLatex(t_of_maximum/(rangeend - rangestart), 0.4, Form("FWHM = %.2f ns", fwhm));

		auto mean = new TLine(expgconv->GetParameter(2), 1e-2, expgconv->GetParameter(2), his->GetMaximum());
		mean->SetLineColor(1); mean->SetLineWidth(2);
		mean->Draw("same");
	}
	else if (do_fit == 3) {
		// sum of two gaussians (one as background estimate)
		auto two_gauss = new TF1("two gaussians", "gaus(0)+gaus(3)", rangestart, rangeend);
		two_gauss->SetTitle("Sum of two gauss");
		float posmax = his->GetXaxis()->GetBinCenter(his->GetMaximumBin());
		two_gauss->SetParameters(his->Integral("width"), posmax, 0.35, his->Integral("width") / 30, posmax, 2);
		two_gauss->SetParName(0, "norm_{peak}");		two_gauss->SetParName(1, "#mu_{peak}");			two_gauss->SetParName(2, "#sigma_{peak}");			two_gauss->SetParLimits(2, 1e-9, 1e2);
		two_gauss->SetParName(3, "norm_{background}");	two_gauss->SetParName(4, "#mu_{background}");	two_gauss->SetParName(5, "#sigma_{background}");	two_gauss->SetParLimits(5, 1e-9, 1e2);
		TFitResultPtr fresults = his->Fit(two_gauss, fitoption.c_str(), "same", fitrangestart, fitrangeend);
		timing_fit_results.push_back(fresults);
	}

	root_out->WriteObject(his, his->GetTitle());
	timing_cfd_d_c->Update();
	root_out->WriteObject(timing_cfd_d_c, "TimingCFD_diff");
}
/// @example timing_example.cc






/// @brief Helper that returns the waveform histogram for a certain waveform number number 
/// @param wfindex Waveform index
/// @return Waveform histogram 
TH1F* ReadRun::Getwf(int wfindex) {
	int channel = GetCurrentChannel(wfindex);
	int event_nr = GetCurrentEvent(wfindex);
	TString name(Form("ch%02d_%05d", channel, event_nr));
	TString title(Form("ch%d, event %d;t [ns];U [mV]", channel, event_nr));
	auto his = new TH1F(name.Data(), title.Data(), binNumber, 0, static_cast<float>(binNumber) * SP);
	for (int i = 1; i <= binNumber; i++) his->SetBinContent(i, rundata[wfindex][i - 1]);
	return his;
}

/// @brief Helper that returns the waveform histogram for a certain channel number and a certain event number
/// @param channelnr Channel number index (not the actual channel number)
/// @param eventnr Event number
/// @param color Choose color of histogram
/// @return Waveform histogram 
TH1F* ReadRun::Getwf(int channelnr, int eventnr, int color) {
	TString name(Form("ch%02d_%05d", channelnr, eventnr));
	TString title(Form("ch%d, event %d;t [ns];U [mV]", channelnr, eventnr));
	auto his = new TH1F(name.Data(), title.Data(), binNumber, 0, static_cast<float>(binNumber) * SP);
	for (int i = 1; i <= binNumber; i++) his->SetBinContent(i, rundata[eventnr * nchannels + channelnr][i - 1]);
	his->SetLineColor(color);
	his->SetMarkerColor(color);
	return his;
}

/// @brief Get array of x axis (time of the bin centers) for standard wavecatcher settings 
/// @tparam T double or float
/// @param shift Offset
/// @return Time array
template<typename T>
T* ReadRun::getx(double shift) {
    T* xvals = new T[binNumber];
    for (int i = 0; i < binNumber; ++i) {
        xvals[i] = static_cast<T>(SP) * (static_cast<T>(i) + 0.5) + shift;
    }
    return xvals;
}
template double* ReadRun::getx<double>(double);
template float* ReadRun::getx<float>(double);


/// @brief Get array of y values for a certain waveform
/// @param waveform_index Waveform index
/// @return Y values of waveform
double* ReadRun::gety(int waveform_index) {
	double* waveform = new double[binNumber];
	for (int i = 0; i < binNumber; i++) waveform[i] = static_cast<double>(rundata[waveform_index][i]);
	return waveform;
}

/// @brief Get array of y values for a certain waveform
/// @param channelnr Channel number index (not the actual channel number)
/// @param eventnr Event number
/// @return Y values of waveform
double* ReadRun::gety(int channelnr, int eventnr) {
	double* waveform = new double[binNumber];
	int wf_index = GetWaveformIndex(eventnr, channelnr);
	for (int i = 0; i < binNumber; i++) waveform[i] = static_cast<double>(rundata[wf_index][i]);
	return waveform;
}

/// @brief Returns index of a certain event number (if data files are read in parallel threads)
/// @param eventnr Event number as stored in the data.
/// @param channel_index Channel index as stored in the data.
/// @return Corresponding waveform number in the internal data structure.
int ReadRun::GetWaveformIndex(int eventnr, int channel_index) {
	return max(0, min(nwf - 1, eventnr * nchannels + channel_index));
}

/// @brief Returns index of a certain event number (if data files are read in parallel threads)
/// @param eventnr Event number as stored in the data.
/// @return Corresponding event number in the internal data structure.
int ReadRun::GetEventIndex(unsigned int eventnr) {
	auto event_pos = find(eventnr_storage.begin(), eventnr_storage.end(), eventnr);
	if (event_pos == eventnr_storage.end()) {
		cout << "WARNING: Event number " << eventnr << " for GetEventIndex() does not exist in data.\n"
			<< "Please check the events in the data or set discard_original_eventnr = true before calling ReadFile()." << endl;
		
		if (static_cast<int>(eventnr) < nevents) {
			eventnr = eventnr_storage[eventnr];
			event_pos = find(eventnr_storage.begin(), eventnr_storage.end(), eventnr);
			cout << "Found event number " << eventnr << " in data and will use it." << endl;
		}
	}
	return static_cast<int>(distance(eventnr_storage.begin(), event_pos));
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

/// @brief Get the current channel index for a certain waveform index
/// @param waveform_index 
/// @return Current channel index
int ReadRun::GetCurrentChannel(int waveform_index) {
	return (waveform_index - nchannels * floor(waveform_index / nchannels));
}

/// @brief Get the current event index for a certain waveform index
/// @param waveform_index 
/// @return Current event index
int ReadRun::GetCurrentEvent(int waveform_index) {
	return floor(waveform_index / nchannels);
}

/// @brief Check if a channel index should be plotted according to ReadRun::plot_active_channels
/// @param i Channel index.
bool ReadRun::PlotChannel(int i) {
	if (plot_active_channels.empty() || Helpers::Contains(plot_active_channels, active_channels[i])) return true;
	else return false;
}

/// @brief Simple linear interpolation for x
/// 
/// Will return the mean of x1 and x2 if y1=y2.
/// 
/// @param ym Y value for evaluation
/// @param x1 X1 
/// @param x2 X2
/// @param y1 Y1
/// @param y2 Y2
/// @param verbose Set true for printing errors
/// @return x value at "ym"
pair<float, bool> ReadRun::LinearInterpolation(float ym, float x1, float x2, float y1, float y2, bool verbose) {
	if (y1 == y2) return {(x1 + x2) / 2., false};
	else if ((y1 > ym && y2 > ym) || (y1 < ym && y2 < ym)) {
		if (verbose){
			cout << "\nError in LinearInterpolation: Value ym=" << ym << " out of range (" << y1 << "|" << y2 << ")." << endl;
			cout << "Will return x1. Increase window for search." << endl;
		}
		return {x1, false};
	}
	else return {x1 + (ym - y1) * (x2 - x1) / (y2 - y1), true};
}

/// @brief Channel map of polarity changes during reading. \n
/// For parameters see ReadFile(). 
vector<bool> ReadRun::PolarityMap(bool change_polarity, int change_sign_from_to_ch_num) {
	vector<bool> pol(nChannelsWC, false);
	if (!change_polarity) {
		return pol;
	}
	else if (!switch_polarity_for_channels.empty()) {
		for (int current_channel=0; current_channel<nChannelsWC; current_channel++) {
			if (Helpers::Contains(switch_polarity_for_channels, current_channel)) {
				pol[current_channel] = true;
			}
		}
		return pol;
	}
	else {
		for (int current_channel=0; current_channel<nChannelsWC; current_channel++) {
			if ((current_channel >= change_sign_from_to_ch_num) || 
				(change_sign_from_to_ch_num < 0 && current_channel <= abs(change_sign_from_to_ch_num))) {
				pol[current_channel] = true;
			}
		}
		return pol;
	}
}

/// @brief Check if index exists in time of waveforms
/// @param index Bin number
/// @return [0, binNumber)
int ReadRun::CheckBoundsX(int index) {
    return min(max(0, index), binNumber - 1);
}

/// @brief Convert time to the bin number of the waveform
/// @param time Time in ns
/// @return Bin number between [0, binNumber) corresponding to the time
int ReadRun::TimeToIndex(float time) {
   	return CheckBoundsX(static_cast<int>(round(time / SP)));
}

/// @brief Convert the bin number of the waveform to the **time of the left bin edge**
/// @param bin_index Bin number
/// @return Time in ns between [0, binNumber) * SP
float ReadRun::IndexToTime(int bin_index) {
	return static_cast<float>(CheckBoundsX(bin_index)) * SP;
}