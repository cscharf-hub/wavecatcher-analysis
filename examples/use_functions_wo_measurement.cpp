#include <iostream>
#include <iomanip>

// example how to use functions of ReadRun class in a macro without loading a measurement
// an example how to use the implemented fit functions and other functions in python is given in cosmics-fit.ipynb
void use_functions_wo_measurement() {
	// initialze dummy of the ReadRun class to use its functions
	ReadRun dummy_class(0);

	// for this example an array is filled with 10/(1+(x-x0)^2):
	int n_entries = 15;
	double* test_arr = new double[n_entries];
	int x0 = static_cast<int>(n_entries / 2);
	cout << "\noriginal: ";
	for (int i = 0; i < n_entries; i++) {
		test_arr[i] = 10. / (1. + static_cast<double>((i - x0) * (i - x0)));
		cout << setprecision(2) << test_arr[i] << "\t";
	}
	
	// now the array is smoothed with a gauss kernel (method = 2) with a sigma of 1.2 using SmoothArray()
	dummy_class.SmoothArray(test_arr, n_entries, 1.2, 2, 1); // sigma = 1.2, method = 2, bin size = 1.
	
	// print content of smoothed array
	cout << "\n\nsmoothed: ";
	for (int i = 0; i < n_entries; i++) cout << setprecision(2) << test_arr[i] << "\t";
	delete[] test_arr;
	cout << endl;
}