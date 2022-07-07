#ifndef OptimizerMVA_h
#define OptimizerMVA_h

#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
 
#include "TFile.h"
#include "TTree.h"
#include "TText.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
 
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include <TStyle.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

class OptimizerMVA{

public:

	// constructor & destructor
	OptimizerMVA();
	~OptimizerMVA();

	void SetUpInTreeSGN();
	void SetUpInTreeBKG();
	void LoadFitRanges();

	// getters
	TFile* GetFile(TString file); 
	void GetCorrelation_BDT_MR();
	int  GetSignalPreCut();
	double GetSGNfactor() {return SGNfactor;};
	double GetAvMultiplicity(const double Xcut, const double MRcut, TString dataset );
	void GetControlPlots(const double Xcut, const double MRcut);

	// analysis methods 
	int ApplyBDT();
	double BKG_NevExtraction(const double Xcut, const double MRcut);
	double SGN_NevExtraction(const double Xcut, const double MRcut);
	double PunziSign(const double S, const double B, double* PSerr);
	int    makeSGNvsBKGplot(const double Xcut, const double MRcut);
	int    makeMASSplot2D(const double Xcut, const double MRcut);


private:

	int SGNprecut;
	double SGNfactor;
// I/O file path

	TString outFile_path;

	TString inFileSGN_path, inFileBKG_path;
	TString SGNtreeName, BKGtreeName;
	TString outFileSGN_path, outFileBKG_path;
	std::string inFileFitPar_path;

// TTree
	float run, LumiBlock, event;
	float pTM_B0, SVprob, LxySign_B0, CosAlpha_B0, DR_Pi1B0, pT_Pi1, pT_Rho, D0_Rho;
	float M_Rho, M_B0, M_MuMu, M_X3872, M_K0s;
	float X;
	TFile* InFileSGN;
	TTree* TreeSGN;
	TFile* InFileBKG;
   TTree* TreeBKG;

// BDT model	
	TString BDTweightPath_;

// signal and background edges
	double SGNregL, SGNregR;
	double BKGregFarL, BKGregNearL, BKGregFarR, BKGregNearR;

		

};



#endif
