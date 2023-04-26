
// example how to use functions of ReadRun class in a macro without loading a measurement
// an example how to use the implemented fit functions and other functions in python is given in cosmics-fit.ipynb


void use_functions_wo_measurement() {
	// for the example an array is filled with numbers and the numbers are printed:
	int n_entries = 10;
	double* test_arr = new double[n_entries];
	cout << "original numbers:" << endl;
	for (int i = 0; i < n_entries; i++) {
		test_arr[i] = i * i;
		cout << i << ": " << test_arr[i] << endl;
	}

	// now the array is smoothed with a gauss kernel (method = 2) with a sigma of 1.5 using SmoothArray()
	ReadRun dummy_class(0);
	dummy_class.SmoothArray(test_arr, n_entries, 1.5, 2, 1); // sigma = 1.5, method = 2, bin size = 1.
	
	// print content of smoothed array
	cout << "smoothed array:" << endl;
	for (int i = 0; i < n_entries; i++) cout << i << ": " << test_arr[i] << endl;
}