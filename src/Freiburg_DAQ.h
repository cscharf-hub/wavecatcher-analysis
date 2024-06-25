/// Class containing Freiburg_DAQ functions
#ifndef _Freiburg_DAQ
#define _Freiburg_DAQ

#include "ReadRun.h"

class Freiburg_DAQ : public virtual ReadRun {
public:
	/// @brief Initializer will call initializer of ReadRun class
	/// @param no_of_bin_files_to_read 
	Freiburg_DAQ(int no_of_bin_files_to_read) : ReadRun(no_of_bin_files_to_read) { }
	
	void ReadFile(string, bool = true, string = "out.root");
};
#endif