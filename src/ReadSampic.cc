#include "ReadSampic.h"

/// @brief Routine to read files created by SAMPIC DAQ.
///
/// Reads the data and can already do simple data manipulation during reading: \n
///   - Change polarity of range of channels with parameters (explained below). \n
///   - Shift all waveforms to a constant fraction such they all start at the same time (-> see ReadRun::Shift_WFs_in_file_loop ) \n
///   - Simple baseline correction by calling CorrectBaseline() before calling ReadFile(). \n
/// Stores the data as ROOT histograms in ```TClonesArray rundata```. \n 
/// SAMPIC reference: https://doi.org/10.1109/NSSMIC.2014.7431231 \n
///
/// Please note that the absolute time is not yet implemented (for CFD timing, plotting waveforms etc.).
///
/// @param path Path to the data. All files in this folder containing ```.bin``` in the file name will be read in.
/// @param change_polarity Set ```true``` to change polarity (sign) of certain channels (see ```change_sign_from_to_ch_num``` below).
/// You can also define a list of channels where the polarity should be switched with switch_polarity_for_channels.
/// @param change_sign_from_to_ch_num All channels \f$ \geq \f$ ```change_sign_from_to_ch_num```
/// will be inverted if ```change_polarity``` is ```true```. \n
/// If negative number all channels \f$ \leq \f$ ```abs(change_sign_from_to_ch_num)``` will be inverted if ```change_polarity``` is ```true```.
/// @param out_file_name Name of the ```.root``` file which stores the results, e. g. ```results.root```.
/// @param debug Set ```true``` to increase the verbosity.
/// @param max_nfs_to_read Maximum number of events to be read per file. For quick testing of analysis on subset of the data.
void ReadSampic::ReadFile(string path, bool change_polarity, int change_sign_from_to_ch_num, string out_file_name, bool debug, long long max_nfs_to_read) {
	if (max_nfs_to_read <= 0) max_nfs_to_read = static_cast<long long>(1e9);	// default 1B waveforms
	rundata.reserve(1'000'000);												// reserve space for 1M waveforms
	hitInfo.reserve(1'000'000);
	eventsBuilt = false; 														// reset flag
	if (path.back() != '/') path += '/';
	data_path = path;
	
	// save results to root file
	if (out_file_name.empty()) out_file_name = "out.root";
	printf("+++ saving analysis results in '%s' ...\n\n", out_file_name.c_str());
	root_out = TFile::Open(out_file_name.c_str(), "recreate");

	// invert channels
	auto polarity_map = PolarityMap(change_polarity, change_sign_from_to_ch_num);

	// create list of binary files in folder
	stringstream inFileList;
	inFileList << Helpers::ListFiles(path.c_str(), ".bin");
	if (debug) cout << inFileList.str() << endl;

	string fileName;
	int file_counter = 0, wfcounter = 0;
	string data_settings("=== DATA STRUCTURE INFO  === REDUCED DATA TYPE: YES === WITHOUT WAVEFORM: NO  === TDC-LIKE FILES: NO === COMPACT BINARY DATA: YES === DATA_IN_FILE_TYPE: 1 ===");

	while (inFileList >> fileName) {
		// read only fraction/batch of the .bin files for testing or to reduce memory usage
		if (FirstBinFileToRead > 0 && FirstBinFileToRead < LastBinFileToRead && file_counter < FirstBinFileToRead) continue;
		if (LastBinFileToRead > 0 && file_counter > LastBinFileToRead) break;

		fileName = path + fileName;
		ifstream input_file(fileName.c_str());
		if (!input_file.is_open()) {
			printf("*** failed to open '%s'\n", fileName.c_str());
			continue;
		}

		if (file_counter < 10 || file_counter % 10 == 0 || debug) printf("+++ reading '%s' ...\n", fileName.c_str());

		string line;
		int header_line = 0;
		float sampling_frequency = 0;
		const float ADC_factor = 0.1; 			// conversion to mV
		int has_measurement = 1; 				// true: contains absolute time, baseline, amplitude
		unordered_set<int> active_channels_set; // store channels in data
		SP = 0.;
		while (header_line < 7) {
			getline(input_file, line);
			if (debug) cout << line.c_str() << endl;
			if (header_line == 0 && line.find("SOFTWARE VERSION: V3.5.34") == string::npos && line.find("SOFTWARE VERSION:") != string::npos) {
				cout << "WARNING: Expected software version V3.5.34 but found:\n" << line.c_str() << endl;
			}
			if (header_line == 3) {
				size_t pos = line.find("DATA_IN_FILE_TYPE:");
    			if (pos != string::npos) has_measurement = stoi(line.substr(pos + 18));
				if (line.substr(0, line.length() - 5) != data_settings.substr(0, data_settings.length() - 5)) {
					cout << "\nWARNING: Unexpected DATA STRUCTURE INFO line.\n Found " << line.c_str() << endl;
					cout << "instead of " << data_settings.c_str() << endl;
				}
				if (has_measurement == 1) cout << "Using HitStructInfoForWaveformAndMeasurements_t." << endl;
				else cout << "Using HitStructInfoForWaveformOnly_t." << endl;
			}
			if (line.find("SAMPLING FREQUENCY") != string::npos) {
				size_t pos = line.find("SAMPLING FREQUENCY ");
				size_t pos2 = line.find(" MS/s");
				if (pos != string::npos) {
					string freq_str = line.substr(pos + 19, pos2);
					sampling_frequency = stof(freq_str);
					SP = 1000/sampling_frequency;
					SP_inv = 1. / SP;
					if (file_counter == 0) cout << "Sampling frequency: " << sampling_frequency << " MS/s" << endl;
				}
			}
			header_line++;
		}
		if (SP == 0.) {
    		throw runtime_error("ERROR: File header could not be read.");
		}

		vector<short> waveform;
		vector<float> waveform_f;
		int previous_binNumber = -1;
		HitStructInfoForWaveformAndMeasurements_t hitData;
		HitStructInfoForWaveformOnly_t hitDataNoMeas;
		HitInfoReduced currentHitInfo;
		while (has_measurement ? 
				input_file.read((char *)(&hitData), sizeof(hitData)) : 
				input_file.read((char *)(&hitDataNoMeas), sizeof(hitDataNoMeas))) {
			//waveform loop
			
			if (!has_measurement) {
				currentHitInfo.HitNumber = hitDataNoMeas.HitNumber;
				currentHitInfo.Channel = static_cast<int>(hitDataNoMeas.Channel);
				currentHitInfo.FirstSampleTimeStamp = hitDataNoMeas.FirstSampleTimeStamp;
				binNumber = static_cast<int>(hitDataNoMeas.WaveformSize);
			}
			else {
				currentHitInfo.HitNumber = hitData.HitNumber;
				currentHitInfo.Channel = static_cast<int>(hitData.Channel);
				currentHitInfo.FirstSampleTimeStamp = hitData.FirstSampleTimeStamp;
				binNumber = static_cast<int>(hitData.WaveformSize);
			}

			if (debug && wfcounter < 10) {
				cout << "hit number: " << currentHitInfo.HitNumber << endl;
				cout << "channel: " << currentHitInfo.Channel << endl;
				cout << "timestamp: " << currentHitInfo.FirstSampleTimeStamp << endl;
				// other stuff that is not yet used
				// cout << "RawToTValue: " << hitData.RawToTValue << endl;
				// cout << "TOTValue: " << hitData.TOTValue << endl;
				// if (has_measurement == 1) {
					// cout << "Time: " << hitData.Time << endl;
					// cout << "Baseline: " << hitData.Baseline << endl;
					// cout << "Amplitude: " << hitData.Amplitude << endl;
				// }
				// cout << "FirstCellIndex: " << static_cast<int>(hitData.FirstCellIndex) << endl;
				cout << "samples: " << binNumber << endl;
			}

			if (wfcounter == 0 && file_counter == 0) { // init
				amplValuessum.resize(nChannelsWC, vector<float>(binNumber, 0.));
			}
			
			if (active_channels_set.insert(currentHitInfo.Channel).second) { // check if channel is new
			    active_channels.push_back(currentHitInfo.Channel);
			}

			if (binNumber != previous_binNumber) {
    			waveform.resize(binNumber);
    			waveform_f.resize(binNumber);
    			previous_binNumber = binNumber;
			}
			
			size_t waveform_bytes = static_cast<size_t>(binNumber) * sizeof(short);
			input_file.read((char *)waveform.data(), waveform_bytes);
			
			float factor = ADC_factor;
			if (polarity_map[currentHitInfo.Channel]) factor = -factor;

			for (int i = 0; i < binNumber; i++) {
				waveform_f[i] = static_cast<float>(waveform[i]) * factor;
				amplValuessum[currentHitInfo.Channel][i] += waveform_f[i];
			}
			rundata.push_back(waveform_f);
			
			currentHitInfo.Max = *max_element(waveform_f.begin(), waveform_f.end());
			currentHitInfo.Min = *min_element(waveform_f.begin(), waveform_f.end());
			hitInfo.push_back(currentHitInfo);

			wfcounter++;
			if (wfcounter > max_nfs_to_read) {
				cout << "Stopped reading waveforms after max_nfs_to_read=" << max_nfs_to_read << endl;
				break;
			}
		}

		input_file.close();
		file_counter++;
	}
	nwf = wfcounter;
	sort(active_channels.begin(), active_channels.end());
	nchannels = static_cast<int>(active_channels.size());

	printf("Finished reading %d files containing %d hits with %d channels.\n\n", file_counter, nwf, nchannels);
}


/// @brief Plot any Waveform (can be called before EventBuilder to check data)
/// @param wfNumber Number of the waveform (hit number)
/// @param ymin scale of y axis
/// @param ymax scale of y axis
void ReadSampic::PlotWF(int wfNumber, float ymin, float ymax) {
	gStyle->SetOptStat(0);
	int channel = hitInfo[wfNumber].Channel;
	TString name(Form("waveform_%05d", wfNumber));
	TString title(Form("wf%d, ch%d;t [ns];U [mV]", wfNumber, channel));
	auto intwinc = new TCanvas(name.Data(), name.Data(), 600, 400);
	auto his = new TH1F(name.Data(), title.Data(), binNumber, 0, IndexToTime(binNumber - 1));
	for (int i = 1; i <= binNumber; i++) his->SetBinContent(i, rundata[wfNumber][i - 1]);
	his->Draw("HIST");
	his->SetStats(0);
	his->GetYaxis()->SetRange(ymin, ymax);
	intwinc->Update();
	root_out->WriteObject(intwinc, name.Data());
}


/// @brief Construct events based on coincident amplitude thresholds within a time window for different channels. \n
/// It currently constructs the events from the first waveform above the set amplitude threshold. \n
/// ToDo: \n
/// - Add option to create event from most coincidences, not from first \n
/// - Add option for integral threshold? (Would require parameters for integration window -> separate function?)
/// @param coincidence_time_window Event time window in ns. All hits within time window will be added to event. \n
/// Starts at the beginning of the first channel above threshold and will include all channels 
/// which start recording before the end of the coincidence time window. This means the peak might 
/// be after the end of the coincidence time window. If you want to include only waveforms with the 
/// signal inside the time window you would need to subtract the length of a waveform (depends on SAMPIC settings) from the time window. \n
/// Must be shorter than the dead time of SAMPIC (<1000 ns).
/// @param thresholds Vector of thresholds for all channels in ascending order (lowest channel number in data is the first entry). \n
/// If you specify only a single value {val}, it will be used for all channels. The amplitude must be larger than the specified threshold.
/// @param channels Vector of channels to be used. If left empty {} all channels will be used. 
/// Provide the actual channel numbers from SAMPIC (e. g. {10, 28}).
/// @param min_coincidences Minimum number of channels above threshold. Events will be constructed if at least
/// min_coincidences channels are above their threshold.
/// @param shift_relative Shift waveforms relative to first timestamp in event \n NOT WORKING YET
void ReadSampic::EventBuilder(double coincidence_time_window, vector<float> thresholds, vector<int> channels, unsigned int min_coincidences, bool shift_relative) {
	cout << "---------Started EventBuilder()---------" << endl;
	
	coincidence_time_window = abs(coincidence_time_window);
	if (coincidence_time_window > 1000) {
		cout << "Warning: coincidence_time_window=" << coincidence_time_window << " ns must be <=1 us. Will use 1 us." << endl;
		coincidence_time_window = 1000;
	}
	// store for plotting (NOT USED AS OF NOW)
	coincidence_time = coincidence_time_window;
	
	// init 
	bool include_all_channels = false;
	if (channels.empty()) {
		channels = active_channels;
		include_all_channels = true;
	}
	else Helpers::filterChannelUserInput(channels, active_channels);
	
	min_coincidences = min(static_cast<unsigned int>(channels.size()), min_coincidences);

	float default_threshold = 0.;
	if (thresholds.size() == 1) default_threshold = thresholds[0];
	thresholds.resize(channels.size(), default_threshold);

	// reset events if called more than once
	if (eventsBuilt) {
		wf_nr_event_storage.clear();
		ch_nr_event_storage.clear();
		skip_event.clear();
		eventnr_storage.clear();
	}
	
	cout << "Coincidence time: " << coincidence_time_window << " ns" << endl;
	cout << "Minimum number of coincidences: " << min_coincidences << endl;

	// check if waveforms above thresholds
	int n_canditates = 0;
	
	#pragma omp parallel for reduction(+:n_canditates)
	for (int i=0; i<nwf; i++) {
		if (eventsBuilt) { 
			hitInfo[i].IsEventCandidate = false;
			hitInfo[i].IsEvent = false;
		}

		if (include_all_channels || Helpers::Contains(channels, hitInfo[i].Channel)) {
			int i_ch = GetChannelIndex(hitInfo[i].Channel);
			if (hitInfo[i].Max > thresholds[i_ch]) {
				hitInfo[i].IsEventCandidate = true;
				n_canditates++;
			}
		}
	}
	cout << "There are " << n_canditates << " event seed candidates out of " << nwf << " waveforms" << endl;

	// check for coincidences
	int event_counter = 0;
	for (int i=0; i<nwf; i++) {
		if (hitInfo[i].IsEventCandidate) {
			unsigned int coincidence_counter = 1;
			double coincidence_window_lo = hitInfo[i].FirstSampleTimeStamp - coincidence_time_window / 2.;
			double coincidence_window_hi = hitInfo[i].FirstSampleTimeStamp + coincidence_time_window / 2.;

			int k_lo = i;
			while (--k_lo >= 0 && hitInfo[k_lo].FirstSampleTimeStamp > coincidence_window_lo) {
				if (hitInfo[k_lo].IsEventCandidate && !hitInfo[k_lo].IsEvent) coincidence_counter++;
			}
			k_lo++;
			int k_hi = i;
			while (++k_hi < nwf && hitInfo[k_hi].FirstSampleTimeStamp < coincidence_window_hi) {
				if (hitInfo[k_hi].IsEventCandidate) coincidence_counter++;
			}
			k_hi--;

			if (coincidence_counter >= min_coincidences) {
				wf_nr_event_storage.push_back(vector<int>());
				ch_nr_event_storage.push_back(vector<int>());
				event_time_stamps.push_back(hitInfo[i].FirstSampleTimeStamp);

				for (int kk = k_lo; kk <= k_hi; kk++) {
					if (shift_relative) { // BROKEN
						int shift_bins = static_cast<int>(floor(event_time_stamps[event_counter] - hitInfo[kk].FirstSampleTimeStamp) / static_cast<double>(SP));
						cout << event_time_stamps[event_counter] - hitInfo[kk].FirstSampleTimeStamp << " shift: " << shift_bins << endl;
						// auto his = (TH1F*)rundata->At(kk);
						// Helpers::ShiftTH1(his, shift_bins);
					}
					hitInfo[kk].EventNumber = event_counter;
					hitInfo[kk].IsEvent = true;
					wf_nr_event_storage[event_counter].push_back(kk);
					ch_nr_event_storage[event_counter].push_back(hitInfo[kk].Channel);
				}
				i = k_hi; // jump to last wf in event
				eventnr_storage.push_back(event_counter);
				skip_event.push_back(false);
				event_counter++;
			}
		}
	}
	nevents = event_counter;
	eventsBuilt = true;
	cout << "Finished EventBuilder() and found " << nevents << " events in data:" << endl;
	cout << left << setw(8) << "Channel" << " | " << right << setw(10) << "Threshold" << " | " << right << setw(10) << "Events" << endl;
	cout << string(37, '-') << endl;
	for (int i = 0; i < nchannels; i++) {
		float thrshld = 0;
		auto it = find(channels.begin(), channels.end(), active_channels[i]);
		if (it != channels.end())  thrshld = thresholds[it - channels.begin()];
		
		cout << left << setw(8) << active_channels[i] << " | " << 
				right << setw(10) << thrshld << " | " << 
				right << setw(10) << Nevents_good(i) << endl;
	}
}

/// @brief Helper that returns the waveform histogram for a certain channel number and a certain event number
/// @param channelnr Channel number index (not the actual channel number)
/// @param eventnr Event number
/// @param color Choose color of histogram
/// @return Waveform histogram 
TH1F* ReadSampic::Getwf(int channelnr, int eventnr, int color) {
	checkData();
	eventnr = min(eventnr, nevents);
	channelnr = min(channelnr, static_cast<int>(active_channels.size()));

	int wfindex = -1;
	for (int i=0; i<static_cast<int>(wf_nr_event_storage[eventnr].size()); i++) {
		if (active_channels[channelnr] == hitInfo[wf_nr_event_storage[eventnr][i]].Channel) {
			wfindex = wf_nr_event_storage[eventnr][i];
			break;
		}
	}

	if (wfindex != -1) {
		int channel = hitInfo[wfindex].Channel;
		int event_nr = hitInfo[wfindex].EventNumber;
		TString name(Form("ch%02d_%05d", channel, event_nr));
		TString title(Form("ch%d, event %d;t [ns];U [mV]", channel, event_nr));
		auto his = new TH1F(name.Data(), title.Data(), binNumber, 0, SP * static_cast<float>(binNumber));
		his->SetLineColor(color);
		his->SetMarkerColor(color);
		for (int i = 1; i <= binNumber; i++) his->SetBinContent(i, rundata[wfindex][i - 1]);
		return his;
	}
	else {
		// return zeroes if there is no waveform
		static TH1F defaultHis("default", "NO DATA; t [ns]; U [mV]", binNumber, 0, SP * static_cast<float>(binNumber));
    	return &defaultHis;
	}
}


/// @brief Returns index of a certain event number (if data files are read in parallel threads)
/// @param eventnr Event number as stored in the data.
/// @param channel_index Channel index as stored in the data.
/// @return Corresponding waveform number in the internal data structure.
int ReadSampic::GetWaveformIndex(int eventnr, int channel_index) {
	int wf_index = -1;
	int ch = active_channels[channel_index];
	if (Helpers::Contains(ch_nr_event_storage[eventnr], ch)) {
		auto it = find(ch_nr_event_storage[eventnr].begin(), ch_nr_event_storage[eventnr].end(), ch);
		int index = static_cast<int>(distance(ch_nr_event_storage[eventnr].begin(), it));
		wf_index = wf_nr_event_storage[eventnr][index];
	}
	return wf_index;
}


/// @brief Get the current channel index for a certain waveform index
/// @param waveform_index 
/// @return Current channel index
int ReadSampic::GetCurrentChannel(int waveform_index) {
	return GetChannelIndex(hitInfo[waveform_index].Channel);
}

/// @brief Get the current event index for a certain waveform index
/// @param waveform_index 
/// @return Current event index. If the waveform is not assigned to an event it will return -1.
int ReadSampic::GetCurrentEvent(int waveform_index) {
	return hitInfo[waveform_index].EventNumber;
}

/// @brief Check if event has skipped flag. Will return true if event flag is false but channel does not exist in data
/// @param event_index Index of the event
/// @param channel_index Index of the channel
bool ReadSampic::SkipEvent(int event_index, int channel_index) {
	if (event_index >= static_cast<int>(skip_event.size()) || event_index < 0) { // wrong input
		return true;
	}
	else if (channel_index == -1) { // default just like parent class ReadRun
		return skip_event[event_index];
	}
	else if (!skip_event[event_index] && Helpers::Contains(ch_nr_event_storage[event_index], active_channels[channel_index])) {
		// if the channel index is passed as well, it will return false only if the channel exists in not skipped event
		return false;
	}
	else {
		return true;
	}
}