/// Class containing SAMPIC data loader and Event Builder
#ifndef _ReadSampic
#define _ReadSampic

#include "ReadRun.h"
// #include "utils/Helpers.h"

class ReadSampic : public virtual ReadRun
{
private:
	/// @brief Coincidence time window used for EventBuilder()
	double coincidence_time;

	#pragma pack(1)
		typedef struct {
			int HitNumber;
			unsigned char Channel;
			double FirstSampleTimeStamp;
			unsigned short RawToTValue;
			float TOTValue;
			float Time;
			float Baseline;
			float Amplitude;
			unsigned char FirstCellIndex;
			unsigned char WaveformSize;
		} HitStructInfoForWaveformAndMeasurements_t;

		typedef struct {
			int HitNumber;
			unsigned char Channel;
			double FirstSampleTimeStamp;
			short RawToTValue;
			float TOTValue;
			unsigned char FirstCellIndex;
			unsigned char WaveformSize;
		} HitStructInfoForWaveformOnly_t;
	#pragma pack()

protected:
	/// @brief Check if the data is grouped in events (SAMPIC needs manual event building) 
	bool eventsBuilt = false;

	void checkData(bool isBaselineCorrection = false) const override {
		if (hitInfo.empty()) {
            throw runtime_error(
				"Error: No data has been loaded yet!\n \
				Please call ReadFile() before calling any functions which manipulate data.\n \
				Aborting execution."
			);
        }
		else if (!eventsBuilt && !isBaselineCorrection) {
			throw runtime_error(
				"Error: Sampic::EventBuilder() must be called before any analysis can be performed."
			);
		}
    }

public:
	/// @brief Initializer will call initializer of ReadRun class
	/// @param last_bin_file Number of last .bin file to be read in.\n
	/// Set it to =>1 in order to constrain the number of .bin files to be read from the target folder.\n
	/// Intended for quick tests on a fraction of the full dataset or for batch reading if combined with min_no_of_bin_files_to_read.
	/// @param first_bin_file Number first of .bin file to be read in. \n
	/// Can be used to batch read large datasets in chunks of files.
	ReadSampic(int last_bin_file = 0, int first_bin_file = 0) : ReadRun(last_bin_file, first_bin_file) {}

	void ReadFile(string, bool = true, int = 0, string = "out.root", bool = false, long long = 0) override;

	void PlotWF(int, float = 0, float = 0);

    void EventBuilder(double, vector<float> = {}, vector<int> channels = {}, unsigned int = 1, bool = false);

    TH1F* Getwf(int, int, int = 1) override;	// channel, eventnr, color
	int GetWaveformIndex(int, int) override;
	int GetCurrentChannel(int) override;
	int GetCurrentEvent(int) override;

	bool SkipEvent(int, int = -1) override;
	
	// ---------------- variables ---------------- 
	#pragma pack(1)
	/// @brief Stores additional information with the waveforms
	struct HitInfoReduced {
		int HitNumber;
		int Channel;
		double FirstSampleTimeStamp;
		int EventNumber = -1;
		float Max = 0;
		float Min = 0;
		bool IsEventCandidate = false;
		bool IsEvent = false;
	};
	#pragma pack()

	/// @brief Stores additional info about wf
	vector<HitInfoReduced> hitInfo;
	/// @brief Stores all waveform numbers assigned to events
	vector<vector<int>> wf_nr_event_storage;
	/// @brief Stores all channel numbers assigned to events
	vector<vector<int>> ch_nr_event_storage;
	/// @brief Stores all event times
	vector<double> event_time_stamps;
};
#endif