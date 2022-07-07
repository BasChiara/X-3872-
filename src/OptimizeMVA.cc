#include "../include/OptimizerMVA.h"

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

int main (int argc, char** argv ){

	// INPUT params	
	if(argc < 3){
		std::cout << "USAGE : ./OptimizeMVA min_Xcut min_MRHOcut [step-size X] [step-size Mrho]"<< std::endl; 
		exit(-1);
	}
	// slicing X
	double X_cut;
	const double min_X = std::stod(argv[1]);
	double Max_X = 0.; //0.2 for large scan
	double Step_X = 0.2;
	if (argc > 3) Step_X = std::stod(argv[3]);
	if ( Step_X >= 0.05 ) Max_X = 0.2; 
	const int Nx = (int)((Max_X - min_X) / Step_X) + 1;
	
	// slicing Mrho
	double MRho_cut;
	const double min_MRho = std::stod(argv[2]);
	const double Max_MRho = 0.7001;
	double Step_MRho = 0.2;
	if (argc == 5) Step_MRho = std::stod(argv[4]);
	const int NmR = (int) ((Max_MRho - min_MRho) / Step_MRho) + 1;
	std::cout << Nx << std::endl;
	std::cout << NmR << std::endl;

	// TREE output
	TString OutFilePath = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/CutOptimization/CutOptimization";
	if( ( Step_X < 0.05 ) || ( Step_MRho < 0.05 ) ) OutFilePath.Append("Plus"); 
	OutFilePath.Append(".root");

	double Bs, Ss, Psig, Psig_err;
	TFile* OutFile = new TFile(OutFilePath, "RECREATE");
	TTree* OutTree = new TTree("CutOpt", "");
	OutTree->Branch("cut_X", &X_cut, "cut_X/D");
	OutTree->Branch("cut_MRho", &MRho_cut, "cut_MRho/D");
	OutTree->Branch("nB_Sreg", &Bs, "Bs/D");
	OutTree->Branch("nS_Sreg", &Ss, "Ss/D");
	OutTree->Branch("SignPunzi", &Psig, "SignPunzi/D");
	OutTree->Branch("SignPunzi_error", &Psig_err, "SignPunzi_error/D");
	

	std::cout << "****  START OPTIMIZATION ****\n" << std::endl;
	OptimizerMVA* opt = new OptimizerMVA();

	//opt->ApplyBDT();
	for (int kx = 0; kx < Nx; kx++){
		X_cut = min_X + kx * Step_X;
		if (X_cut > Max_X) X_cut = Max_X;
		for (int im = 0; im < NmR; im++){
			MRho_cut = min_MRho + im * Step_MRho;
			if(MRho_cut > Max_MRho) MRho_cut = Max_MRho;
			std::cout << X_cut << "\t" << MRho_cut << std::endl;
			// STEP 2 --> # BKG events extraction
			Bs = opt->BKG_NevExtraction(X_cut, MRho_cut);
			//std::cout << "#BKG events in signal region " << Bs << std::endl;
			// STEP 3 --> # SGN events extraction
			Ss = opt->SGN_NevExtraction(X_cut, MRho_cut);
			// STEP 4 --> # PUNZI SIGNIFICANCE extraction\n\n\n" << std::endl;
			Psig = opt->PunziSign(Ss, Bs, &Psig_err);
			//std::cout << "SIGNIFICANCE = " << Psig << std::endl;
			OutTree->Fill();
		}
	}
	OutFile->cd();
	OutTree->Write();
	std::cout << " \n\nWRITTEN TREE \"" << OutTree->GetName() << "\" IN FILE " << OutFilePath << std::endl;
	
	OutFile->Close();
	return 0;
}
