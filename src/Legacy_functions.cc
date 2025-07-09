#include "Legacy_functions.h"

/// @brief Find events with max/min above/below a certain threshold
/// 
/// Needs to be called before the charge spectrum etc functions. \n 
/// Baseline correction should be called before this function.
/// **Sequential Clone() loop makes this inefficient. Please let me know if the function is still in use and should be updated.**
/// 
/// @param threshold Threshold in mV.
/// @param max If true uses max, else uses min.
/// @param greater If true looks for events with max/min>threshold, else looks for events with max/min<threshold.
/// @param from Start search at "from" in ns.
/// @param to End search at "to" in ns.
/// @param verbose Set true for extra verbosity.
void Legacy_functions::FractionEventsAboveThreshold(float threshold, bool max, bool greater, double from, double to, bool verbose) {

	int occurences = 0;
	int occurences2ch = 0;
	int o2ch = 0;
	int current_channel = 0;
	int currevent = 0;
	int lastevent = 0;
	if (plot_active_channels.empty()) plot_active_channels = active_channels;
	vector<int> counter_abovethr(static_cast<int>(plot_active_channels.size()));
	// DORAMAS: It stores a counter of events above threshold for each channel that will be plotted

	cout << "\nSearching for fraction of events with " << (max ? "max" : "min") << (greater ? " > " : " < ") << threshold << " mV:" << endl;

	for (int j = 0; j < nwf; j++) {
		auto his = (TH1F*)(Getwf(j))->Clone(); // use Clone() to not change ranges of original histogram

		// set range where to search for amplitudes above threshold
		if (from >= 0 && to > 0) {
			his->GetXaxis()->SetRange(his->GetXaxis()->FindBin(from), his->GetXaxis()->FindBin(to));
		}

		current_channel = GetCurrentChannel(j);
		if (Helpers::Contains(plot_active_channels, active_channels[current_channel])) {
			if ((max && greater && his->GetMaximum() > threshold) || (max && !greater && his->GetMaximum() < threshold) || (!max && greater && his->GetMinimum() > threshold) || (!max && !greater && his->GetMinimum() < threshold)) {
				currevent = eventnr_storage[GetCurrentEvent(j)];

				if (verbose) cout << "\nevent:\t" << currevent << "\tchannel:\t" << active_channels[current_channel];

				// We must use 'distance' to make sure the position in 'counter_above' matches with the corresponding channel's position at 'plot_active_channels'
				counter_abovethr[distance(plot_active_channels.begin(), find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[current_channel]))] += 1;
				// This is to detect events w/ at least two channels above threshold
				if (lastevent != currevent) occurences += 1;
				if (lastevent == currevent && o2ch != occurences) {
					occurences2ch += 1;
					o2ch = occurences;
				}
				lastevent = currevent;
			}
		}
		delete his;
		Helpers::PrintProgressBar(j, nwf);
	}

	//  Loop to show the fraction of events above threshold for each channel that will be plotted
	for (int i = 0; i < static_cast<int>(plot_active_channels.size()); i++) {
		cout << "Fraction of events in channel " << plot_active_channels[i] << " above threshold: "
			<< 100. * static_cast<float>(counter_abovethr[i]) / static_cast<float>(nevents) << "%\n";
	}

	cout << "Fraction of events w/ at least 2 channels above threshold: "
		<< 100. * static_cast<float>(occurences2ch) / static_cast<float>(nevents) << "%\n"
		<< "For a total of " << nevents << " events" << endl;
}