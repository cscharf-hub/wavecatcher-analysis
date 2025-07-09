/// Helper functions
#ifndef _Helpers
#define _Helpers

#include <TSystemDirectory.h>
#include <TList.h>
#include <TObjString.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TMath.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

class Helpers {
public:
	// find data files
	static string ListFiles(const char*, const char*);
	
	// progress bar
	static void PrintProgressBar(int, int);

	// use in loop, skips some poorly visible root colors (like white on white)
	static int rcolor(unsigned int);

	// set consistent ranges
	static void SetRangeCanvas(TCanvas*&, double, double, double = -999, double = -999);
	// check if user input makes sense
	static void filterChannelUserInput(vector<int>&, const vector<int>);
	// split canvas into pads to display all active channels on one canvas
	static void SplitCanvas(TCanvas*&, vector<int>, vector<int>);

	// get y values of a histogram
	template <typename HistType>
    static double* gety(HistType* his);

	// get y values of a histogram for a dedicated x range
	template <typename HistType>
	static double* gety(HistType*, int, int);
	static double* gety(const vector<float>&, int, int);

	// shift a TH1 in x
	template <typename HistType>
	static void ShiftTH1(HistType*&, int);

	/// @brief Returns true if vector vec contains value val
	/// @tparam T Type of the elements
	/// @param vec Vector to search for ->
	/// @param val Value 
	template<typename T>
	static bool Contains(const vector<T>& vec, const T& val) {
    	return find(vec.begin(), vec.end(), val) != vec.end();
	}
};
#endif