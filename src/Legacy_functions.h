/// Class containing legacy functions to be removed in the future. Let me know if they are still in use.
#ifndef _Legacy_functions
#define _Legacy_functions

#include "ReadRun.h"

class Legacy_functions : public virtual ReadRun {
public:
	/// @brief Initializer will call initializer of ReadRun class
	/// @param no_of_bin_files_to_read Parameter not yet implemented
	Legacy_functions(int no_of_bin_files_to_read) : ReadRun(no_of_bin_files_to_read) { }
	
	void FractionEventsAboveThreshold(float = 4, bool = true, bool = true, double = 0., double = 0., bool = false);
};
#endif