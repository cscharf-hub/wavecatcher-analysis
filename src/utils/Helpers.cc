#include "Helpers.h"

/// @brief Helper. Creates a list of .bin data files in data folder to be read in
/// @param dirname Directory
/// @param ext File extension
/// @return String of line separated file names
string Helpers::ListFiles(const char* dirname, const char* ext) {
	stringstream ss;
	TSystemDirectory dir(dirname, dirname);
	TList* files = dir.GetListOfFiles();
	
	TIter next(files);
	TObjString* objString;
	while ((objString = (TObjString*)next())) {
		const char* fileName = objString->GetString().Data();
		// skip hidden/temporary files from backup etc.
		if (fileName[0] == '.') continue;
		// Check if the filename or extension contains ".bin"
		if (strstr(fileName, ext) != nullptr) ss << fileName << "\n";
	}
	if (ss.str().empty()) throw runtime_error("Error: No .bin files found in " + static_cast<string>(dirname) + "\n");
	return ss.str();
}

/// @brief Print progress bar for a loop in steps of 10 percent
/// @param index Current loop index
/// @param length Length of loop
void Helpers::PrintProgressBar(int index, int length) {
	static int lastProgress = -1;
    int progress = (10 * (index + 1)) / length;
    if (progress != lastProgress) {
        lastProgress = progress;
        cout << "\r[" << string(progress * 5, '#')
             << string(50 - progress * 5, ' ')
             << "] " << progress * 10 << "%";
        cout.flush();
    }
    if (index + 1 == length) cout << endl << endl;
}

/// @brief Translate a random number into a useful root color https://root.cern.ch/doc/master/classTColor.html
/// @param i Index of your plotting loop that is to be translated into a useful ROOT color index
/// @return ROOT color index
int Helpers::rcolor(unsigned int i) {
	const int nclrs = 17;
	int rclrs[nclrs] = { 1, 2, 3, 4, 6, 8, 9, 13, 20, 28, 30, 34, 38, 40, 31, 46, 49 };
	return rclrs[i - static_cast<int>(floor(i / nclrs)) * nclrs];
}

/// @brief Set consistent x-axis and y-axis range for all TH1 histograms on a canvas. 
/// 
/// Will only change the axes where ```min != max```.
/// 
/// @param c Canvas Canvas
/// @param x_range_min X-axis minimum
/// @param x_range_max X-axis maximum
/// @param y_range_min Y-axis minimum
/// @param y_range_max Y-axis maximum
void Helpers::SetRangeCanvas(TCanvas*& c, double x_range_min, double x_range_max, double y_range_min, double y_range_max) {

	// Lambda function to set the axis ranges
	auto setAxisRanges = [&](TH1* his) {
		if (x_range_min != x_range_max) his->GetXaxis()->SetRangeUser(x_range_min, x_range_max);
		if (y_range_min != y_range_max) his->GetYaxis()->SetRangeUser(y_range_min, y_range_max);
		};

	// Get the list of pads on the canvas and loop over pads
	TList* pads = c->GetListOfPrimitives();
	TIter nextPad(pads);
	TObject* object;
	while ((object = nextPad())) {
		if (object->InheritsFrom(TPad::Class())) {
			TPad* pad = static_cast<TPad*>(object);
			// Set the axis ranges on the current pad
			pad->cd();
			TList* primitives = pad->GetListOfPrimitives();
			TIter nextPrimitive(primitives);
			TObject* primitive;
			while ((primitive = nextPrimitive())) {
				if (primitive->InheritsFrom(TH1::Class())) {
					setAxisRanges(static_cast<TH1*>(primitive));
				}
			}
		}
		else if (object->InheritsFrom(TH1::Class())) {
			setAxisRanges(static_cast<TH1*>(object));
		}
	}
	c->Modified();
	c->Update();
}


/// @brief Check if user input exists in data and remove channels that are not there
/// @param user_channels User input
/// @param active_channels From data
void Helpers::filterChannelUserInput(vector<int>& user_channels, const vector<int> active_channels) {
    vector<int> rmv;
    
    for (int channel : user_channels) {
        if (!Helpers::Contains(active_channels, channel)) {
            cout << "\n ------------ WARNING ------------\n" 
			<< "YOUR SELECTED CHANNEL " << channel << " DOES NOT EXIST IN DATA" << endl;
            rmv.push_back(channel);
        }
    }

    for (int channel : rmv) {
        auto it = find(user_channels.begin(), user_channels.end(), channel);
        if (it != user_channels.end()) {
            user_channels.erase(it);
        }
    }
}


/// @brief Helper to split canvas according to the number of channels to be plotted
/// @param c Canvas to be split
/// @param active_channels All channels available in the data
/// @param plot_active_channels The channels from the data that should be plotted
void Helpers::SplitCanvas(TCanvas*& c, vector<int> active_channels, vector<int> plot_active_channels) {
	// cross check if user input exists in data
	filterChannelUserInput(plot_active_channels, active_channels);

	if (plot_active_channels.empty()) {
		c->Divide(TMath::Min(static_cast<double>(active_channels.size()), 4.), TMath::Max(TMath::Ceil(static_cast<double>(active_channels.size()) / 4.), 1.), 0, 0);
	}
	else if (static_cast<int>(plot_active_channels.size()) > 1) {
		c->Divide(TMath::Min(static_cast<double>(plot_active_channels.size()), 4.), TMath::Max(ceil(static_cast<double>(plot_active_channels.size()) / 4.), 1.), 0, 0);
	}
}

/// @brief Get array of y values for a histogram
/// @tparam HistType Histogram derived from ROOT TH1
/// @param his Pointer to histogram
/// @return Array of Y values
template <typename HistType>
double* Helpers::gety(HistType* his) {
    static_assert(is_base_of<TH1, HistType>::value, "ERROR in Helpers::gety():\n Argument must be ROOT TH1 histogram.");
    double* yvals = new double[his->GetNbinsX()];
    for (int i = 0; i < his->GetNbinsX(); ++i) {
        yvals[i] = static_cast<double>(his->GetBinContent(i + 1));
    }
    return yvals;
}
template double* Helpers::gety<TH1F>(TH1F*);
template double* Helpers::gety<TH1D>(TH1D*);
template double* Helpers::gety<TH1I>(TH1I*);

/// @brief Get truncated array of y values for a certain waveform
/// @tparam HistType Histogram derived from ROOT TH1
/// @param his Waveform histogram
/// @param start_at Truncate from index
/// @param end_at Truncate to index
/// @return Truncated Y values of waveform
template <typename HistType>
double* Helpers::gety(HistType* his, int start_at, int end_at) {
	static_assert(is_base_of<TH1, HistType>::value, "ERROR in Helpers::gety():\n Argument must be ROOT TH1 histogram.");
	
	if (start_at < 0 || end_at > his->GetNbinsX() || end_at <= start_at) {
		cout << "\nError: Helpers::gety out of range" << endl;
		return nullptr;
	}
	const int n_bins_new = end_at - start_at;
	double* yvals = new double[n_bins_new];
	for (int i = start_at; i < end_at; i++) {
		yvals[i - start_at] = his->GetBinContent(i + 1);
	}
	return yvals;
}
template double* Helpers::gety<TH1F>(TH1F*, int, int);
template double* Helpers::gety<TH1D>(TH1D*, int, int);
template double* Helpers::gety<TH1I>(TH1I*, int, int);

/// @brief Get truncated array of y values for a certain waveform
/// @param waveform Waveform 
/// @param start_at Truncate from index
/// @param end_at Truncate to index
/// @return Truncated Y values of waveform
double* Helpers::gety(const vector<float>& waveform, int start_at, int end_at) {
	if (start_at < 0 || end_at > static_cast<int>(waveform.size()) || end_at <= start_at) {
		cout << "\nError: Helpers::gety out of range" << endl;
		return nullptr;
	}
	const int n_bins_new = end_at - start_at;
	double* yvals = new double[n_bins_new];
	for (int i = start_at; i < end_at; i++) {
		yvals[i - start_at] = static_cast<double>(waveform[i]);
	}
	return yvals;
}

/// @brief Shift a histogram in x \n
/// The histogram will be cycled, so bins at the end will be attached at the front and vice versa
/// @tparam HistType Histogram derived from ROOT TH1
/// @param his Histogram to be shifted
/// @param shift_bins Number of bins in x to shift by
template <typename HistType>
void Helpers::ShiftTH1(HistType*& his, int shift_bins) {
	static_assert(is_base_of<TH1, HistType>::value, "ERROR in Helpers::gety():\n Argument must be ROOT TH1 histogram.");
	int nbins = his->GetNbinsX();
	shift_bins = shift_bins % nbins;
	double* yvals = Helpers::gety(his);

	for (int i = 0; i < nbins; i++) {
		int new_bin = (i + shift_bins) % nbins;
        if (new_bin < 0) new_bin += nbins;
		his->SetBinContent(i + 1, yvals[new_bin]);
	}
	delete[] yvals;
}
template void Helpers::ShiftTH1<TH1F>(TH1F*&, int);
template void Helpers::ShiftTH1<TH1D>(TH1D*&, int);
template void Helpers::ShiftTH1<TH1I>(TH1I*&, int);