#include <iostream>
using namespace std;

void batch_read() // main
{
	// example by Ben Skodda to batch run over several large runs (many events) one channel at a time, to reduce memory usage
	// call "root -b -x "batch_read.cc()" -q"

	for (int i = 0; i < 2; i++) {	//loop through selected measurements
		for (int j = 0; j < 8; j++) {	//loop through channel
			string command;
			command = "root -b 'readout.cc(";
			command += to_string(i);
			command += ",";
			command += to_string(j);
			command += ")' -q";

			cout << command;
			system(command.c_str());
		}
	}
}
