/// Helper functions
#ifndef _Helpers
#define _Helpers

#include <TSystemDirectory.h>
#include <TList.h>
#include <TSystemFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TMath.h>

#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

class Helpers {
public:
	// generate array of gauss
	static double* Gauss(int, double = 1, double = 0., double = 1., bool = false);

	// find data files
	static string ListFiles(const char*, const char*);
	
	// progress bar
	static void PrintProgressBar(int, int);

	// use in loop, skips some poorly visible root colors (like white on white)
	static int rcolor(unsigned int);

	// set consistent ranges
	static void SetRangeCanvas(TCanvas*&, double, double, double = -999, double = -999);
	// split canvas into pads to display all active channels on one canvas
	static void SplitCanvas(TCanvas*&, vector<int>, vector<int>);

	// get y values of a histogram
	static double* gety(TH1F*);
	// get y values of a histogram for a dedicated x range
	static double* gety(TH1F*, int, int);
};
#endif