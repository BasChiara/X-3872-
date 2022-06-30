#ifndef HLTapply_h
#define HLTapply_h

//#include "./PreSelDATA2017.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include "TTree.h"  
#include "TCanvas.h"  
#include "TStyle.h"  
#include <iostream>
#include <fstream>
#include <stddef.h>
#include <vector>
#include <string>

#include <TMath.h>
#include "Math/Vector4D.h"
#include <Math/GenVector/VectorUtil.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include "TLorentzVector.h"

class HLTapply : public PreSelDATA2017{

	public:
		//CONSTRUCTOR - DECONSTRUCTOR
		HLTapply(TTree *tree=0, const TString dataset = "MC");
		virtual ~HLTapply();
		// Loop() from PreSelDATA2017
		void Loop();

		//new functions
		void GetFitParams();

		void OutTree_setup();
		int  RecoPartFillP4(const int Bidx);

		int  ApplyTriggerSelection_Muons(const int Bidx);
		int  ApplyTriggerSelection_Track(const int Bidx);
		
		bool inBKGregion(const TString ptl);
		void NK0s_per_Rho(TH1* histo);
		void NRho_per_K0s(TH1* histo);

		void PrintGenLevel(const int Bidx, const int ptlID, const bool isMCmatched_Pi1, const bool isMCmatched_Pi2, TH1* DRhisto);

	private:
		
		TString dataset_; 

		// outputs
		std::string FileFitParB0, FileFitParK0s, FileFitParX;
		TString outFileHistoPath_;
		TFile*  outFileHisto_;
		TString outFileTreePath_;
		TFile*  outFileTree_;
		TTree*  outTree_;

		// PDG particle ID
		int isMum = 13, isPip = 211, isJPsi = 443, isRho = 113, isK0s = 310, isX3872 = 9920443, isB0 = 511;
		// PDG particle mass
		const float mMuon = 0.105658, mPion = 0.1395704, mK0s = 0.497648; 
		
		// RECO tracks P4
		ROOT::Math::PtEtaPhiMVector P4_Reco_Mu1, P4_Reco_Mu2; 
		ROOT::Math::PtEtaPhiMVector P4_Reco_JPsi;
		ROOT::Math::PtEtaPhiMVector P4_Reco_Pi1, P4_Reco_Pi2;
		ROOT::Math::PtEtaPhiMVector P4_Reco_Rho;
		ROOT::Math::PtEtaPhiMVector P4_Reco_X3872;
		ROOT::Math::PtEtaPhiMVector P4_Reco_K0s;
		ROOT::Math::PtEtaPhiMVector P4_Reco_B0;
		// B0 mass sidebands
		double MB_nearLeft, MB_nearRight, MB_farLeft, MB_farRight;
		double MX_nearLeft, MX_nearRight;
		double MK0s_nearLeft, MK0s_nearRight;
		// B0 sidebands (TTree branches)
		float M_B0, M_Rho, M_X3872, M_mumu, M_K0s;	
		float pTM_B0, LxySign_B0, SVprob, CosAlpha_B0;
		float pT_Rho, pT_Pi1, DR_Pi1B0, D0_Rho;




}; //HLTapply

#endif 
