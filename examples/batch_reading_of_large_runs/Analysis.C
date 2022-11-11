#include <iostream>
using namespace std;

class Fitf_2 {
public:
	// sum of two spectra for event spectrum + dark count (background trigger) spectrum (missing after-pulses and dark counts)

	double operator() (double* x, double* p) {
		//0 - N0: Normalization (~Number of events)
		//1 - mu: for generalized poisson distribution
		//2 - lambda: Borel-branching parameter for prompt crosstalk probability 1-exp(-lambda)
		//3,4 -sigma0, sigma1
		//5 - G: gain
		//6 - B: Pedestal
		//7 - mu_dk: mu for dark count rate in dark events (dark events = noise/background triggers without photons in SiPMs)
		//8 - N0_dk: fraction (dark events)/(total number of events)

		double sum = 0;
		for (int kint = 0; kint <= 50; kint++) {
			double mu = p[1];
			double lambda = p[2];

			double sigma0 = p[3];
			double sigma1 = p[4];

			double G = p[5];
			double B = p[6];

			double mu_dk = p[7];
			double N0_dk = p[8];

			double k = static_cast<double>(kint);
			double sigmaK = sqrt(sigma0 * sigma0 + k * sigma1 * sigma1);
			//generalized poisson envelope
			double gp = ((1. - N0_dk) * (mu * TMath::Power((mu + k * lambda), k - 1) * TMath::Exp(-(mu + k * lambda))) + N0_dk * (mu_dk * TMath::Power((mu_dk + k * lambda), k - 1) * TMath::Exp(-(mu_dk + k * lambda)))) / TMath::Factorial(kint);

			sum += p[0] * gp * (1. / sqrt(2. * TMath::Pi()) / sigmaK) * TMath::Exp(-TMath::Power(((x[0] - (k * G + B)) / sqrt(2) / sigmaK), 2));
		}
		return sum;
	};
};


double* FITT(TH1F* hist, double range_start, double range_end, double p0, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8) {
	Fitf_2 fitf;
	TF1* f = new TF1("fitf", fitf, range_start, range_end, 9); 
	f->SetLineColor(3);

	f->SetParName(0, "N_{0}");				f->SetParameter(0, p0);
	f->SetParLimits(0, 1e4, 1e9);

	f->SetParName(1, "#mu");				f->SetParameter(1, p1);
	f->SetParLimits(1, .1, 1e2);

	f->SetParName(2, "#lambda");			f->SetParameter(2, p2);
	f->SetParLimits(2, 1e-2, 0.5);			//f->FixParameter(2, 1.);
	
	f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, p3);
	f->SetParLimits(3, 1., 25);

	f->SetParName(4, "#sigma_{1}");			f->SetParameter(4, p4);
	f->SetParLimits(4, 1., 10);

	f->SetParName(5, "Gain");				f->SetParameter(5, p5); 
	f->SetParLimits(5, 40, 70);				//f->FixParameter(5, 40.);

	f->SetParName(6, "Pedestal");			f->SetParameter(6, p6);
	f->SetParLimits(6, -10., 10.);			//f->FixParameter(6, 0.);

	f->SetParName(7, "#mu_darkcount");		f->SetParameter(7, p7);
	f->SetParLimits(7, 0., .5);				//f->FixParameter(7, .3);

	f->SetParName(8, "N_{0}_darkcount");	f->SetParameter(8, p8); 
	f->SetParLimits(8, 0., .1);				//f->FixParameter(8, .05);

	hist->Fit(f, "L");
	f->Draw("SAME");
	double* results = new double[4];

	results[0] = f->GetParameter(1);
	results[1] = f->GetParError(1);
	results[2] = f->GetParameter(5);
	results[3] = f->GetParError(5);

	return results;
}

void Analysis() {
	// TFile *outfile  = TFile::Open("Analysis_filter/results.root","RECREATE");
	int Category;
	TTree* tree = new TTree("T", "CERN 1988 staff data");
	tree->Branch("Comparison", &Category, "Category/I");
	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);

	string path = "/mnt/c/SHiP/macros/split_into_channels/";
	string folder;
	string dataname;
	string path_to_data;
	string channelname;
	string feature_tag;
	double* fit = new double[4]; //Save fit results
	double pos_fiber_arr[8][2] = { {0,0}, {1,0}, {2,1.5}, {3,1.5}, {4,12}, {5,13.5}, {6,13.5}, {7,12} };	//Positions of the Fibre {channel , position}
	double pos_source;  //Position of the source
	int color = 0;

	// THIS HAS TO BE AJUSTED MANUALLY{
	int measurement_range[2] = { 0, 0 };		//Range of measurements {fist measurement, last measurement}
	int start_new_Series = 2;
	int active_channels[] = { 0, 1, 2, 3, 4, 5, 6, 7 };	//List of channels you want to look at
	
	// fit parameters
	double fp[9] = { 7.5e5, 9., .23, 9, 7., 54., 0., 0.3, 0.05 };
	// fit range
	double range_start = -25;
	double  range_end = 1200;

	//Array to save results for later plotting  results_arr[number of measurments][number of channels][number of saved parameters]
	double results_arr[int(measurement_range[1] + 1 - measurement_range[0])][int(sizeof(active_channels) / sizeof(active_channels[0]))][6];

	auto legend_h_ch = new TLegend(0.5, 0.6, 0.9, 0.9);
	TCanvas* c2 = new TCanvas("Charge Spectra Ch2", "c2", 10, 10, 2000, 800);

	c2->Divide(4, 2);	//4X2 subplots to show 8 channels

	for (int k = 0; k <= measurement_range[1] - measurement_range[0]; k++) {		//loop over all selected measurements
		for (int which_ch = 0; which_ch < sizeof(active_channels) / sizeof(active_channels[0]); which_ch++) { //loop over all selected channel

			if (k + measurement_range[0] < start_new_Series) {
				folder = "root_files";
				color = 1;
				feature_tag = "completely wrapped";
				switch (k + measurement_range[0]) {
				case(0): {dataname = "S001M001"; pos_source = 0; break; }
				case(1): {dataname = "S001M002"; pos_source = 1.5; break; }
				}
			}
			else {
				folder = "Series_2";
				color = 2;
				feature_tag = "removed foil on the side";
				switch (k + measurement_range[0]) {
				case(2): {dataname = "S002M001"; pos_source = 0; break; }
				case(3): {dataname = "S002M002"; pos_source = 1.5; break; }
				}
			}

			//feature_tag = "Measurement"+to_string(k+measurement_range[0]); //substitute for actual useful measurement describtion
			path_to_data = path + folder + "/Ch_" + to_string(active_channels[which_ch]) + "_" + dataname + ".root";
			cout << "\n\n" << path_to_data << "\n";

			//manipulation of cannel name string
			if (which_ch < 10) channelname = "ChargeSpectrum channel_0" + to_string(active_channels[which_ch]);
			else channelname = "ChargeSpectrum channel_" + to_string(active_channels[which_ch]);


			TH1F* h_ch = new TH1F("h_ch", "0.0cm", 1000, 0, 1500);
			TFile* f3 = TFile::Open(path_to_data.c_str()); //call data from path_to_data
			h_ch = (TH1F*)f3->Get(channelname.c_str()); //call channel from this data

			h_ch->SetLineColor(color); // we need +1 here, because color 0 is white, white is not usefull for plotting ;) 
			//Set histogramm axis range
			//h_ch->GetXaxis()->SetRangeUser(0, 999);
			h_ch->GetYaxis()->SetRangeUser(0, 5500);

			c2->cd(active_channels[which_ch] + 1);	//choose subplot according to channel (+1 because subplots start with 1 not with 0)
			h_ch->Draw("SAME");	//draw histogram 

			// do fit
			fit = FITT(h_ch, range_start, range_end, fp[0], fp[1], fp[2], fp[3], fp[4], fp[5], fp[6], fp[7], fp[8]);

			if (which_ch == 0) {	//condition to plot legend only once     
				legend_h_ch->AddEntry(h_ch, feature_tag.c_str(), "l");
			}
			legend_h_ch->Draw("SAME");


			//Saving all fit results
			results_arr[k][which_ch][0] = sqrt(pow(pos_fiber_arr[active_channels[which_ch]][1] - pos_source, 2.0)); //distance between fibre and source
			results_arr[k][which_ch][1] = 0;	//error (source-fibre distance)
			results_arr[k][which_ch][2] = fit[0];	//mu
			results_arr[k][which_ch][3] = fit[1];	//error (mu)
			results_arr[k][which_ch][4] = fit[2];	//gain
			results_arr[k][which_ch][5] = fit[3];	//error (gain)
			//cout << results_arr[k][which_ch][2] << "and" << which_ch << "\n";
		}
	}


	//new plot for mu and gain
	auto* c = new TCanvas();
	auto legend_mu = new TLegend(0.3, 0.7, 0.9, 0.9);
	auto legend_gain = new TLegend(0.3, 0.7, 0.9, 0.9);
	c->Divide(2, 1);
	auto mg = new TMultiGraph();
	auto mg2 = new TMultiGraph();

	for (int k = 0; k <= measurement_range[1] - measurement_range[0]; k++) {		//loop over selected measurements
		for (unsigned int which_ch = 0; which_ch < sizeof(active_channels) / sizeof(active_channels[0]); which_ch++) {	//loop over selected channels
			if (k < start_new_Series - measurement_range[0]) { color = 1; feature_tag = "completely wrapped detector"; }
			else { color = 2; feature_tag = "removed foil on the side"; }

			c->cd(1); //Subplot for mu 	//this Line is not needed here but good to illustrate the structure
			auto* g1 = new TGraphErrors();
			g1->SetPoint(0, results_arr[k][which_ch][0], results_arr[k][which_ch][2]);	//creating a point_0 at (distance fibre-source, mu)
			g1->SetPointError(0, results_arr[k][which_ch][1], results_arr[k][which_ch][3]);// creating errorbars for point_0 (error(source-fibre sistance), errir(mu))
			g1->SetLineColor(color);
			//cout << results_arr[k][which_ch][2] << "and" << which_ch << "\n";

			//setting marker style
			g1->SetMarkerSize(1.2);
			g1->SetMarkerStyle(19 + color);
			g1->SetMarkerColor(color);
			if ((k == 0 || k == (measurement_range[1] - measurement_range[0])) && which_ch == 0) {
				legend_mu->AddEntry(g1, feature_tag.c_str(), "l");
			}
			mg->Add(g1, "P"); //add point to multigraph mg
			mg->SetTitle("Light yield for alpha; Distance to source [cm]; Detected photonnumber");

			c->cd(2); //subplot for the gain analogue to cd(1)
			auto* g2 = new TGraphErrors();
			g2->SetPoint(0, results_arr[k][which_ch][0], results_arr[k][which_ch][4]);
			g2->SetPointError(0, results_arr[k][which_ch][1], results_arr[k][which_ch][5]);
			g2->SetLineColor(color);
			g2->SetMarkerSize(1.2);
			g2->SetMarkerStyle(19 + color);
			g2->SetMarkerColor(color);
			if ((k == 0 || k == (measurement_range[1] - measurement_range[0])) && which_ch == 0) {
				legend_gain->AddEntry(g2, feature_tag.c_str(), "l");
			}
			mg2->Add(g2, "P");
			mg2->SetTitle("Gain for alpha; Distance to source [cm]; Gain");
		}
	}

	//draw multigraphs in matching subplots
	c->cd(1);
	mg->Draw("A");
	legend_mu->Draw("SAME");
	c->cd(2);
	mg2->Draw("A");
	legend_gain->Draw("SAME");
}
