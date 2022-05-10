#include "../include/SGNvsBKGvariables.h" 

// CONSTRUCTOR-DESTRUCTOR
SGNvsBKGvariables::SGNvsBKGvariables(TTree *tree) : B0toX3872K0s(tree){ 
	
	inFilePath_ = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/B0_SGNvsBKG_anlaysis.root";

   P4_Reco_Mu1.SetM(mMuon); P4_Reco_Mu2.SetM(mMuon);
   P4_Reco_Pi1.SetM(mPion); P4_Reco_Pi2.SetM(mPion);
   P4_Reco_K0s.SetM(mK0s);
}

SGNvsBKGvariables::~SGNvsBKGvariables(){
	inFileTree_->Close();
	outFileHisto_->Close();
}





void SGNvsBKGvariables::Loop() {

   if (fChain == 0) return;
// LOAD THE SGN vs BKG TREE
	ReadTree();
   Long64_t NB0cand = inTree_->GetEntriesFast();
// LOAD Ntuples
   Long64_t Nevents = fChain->GetEntriesFast();

	std::cout << "Analysing " << NB0cand << " B0 candidates out of " << Nevents << " events\n" << std::endl;

   Long64_t Nbreak = NB0cand + 10, NPrint = NB0cand/20.; 
	Long64_t jentry;

//---- VARIABLES ----//
	int Bidx;

//---- HISTOGRAMS ----//
	int nbins = 50;
	//... Mass
	TH1F h_TOTfit_SGN_B0_M("TOTfit_SGN_B0_M", "", nbins,5., 5.5 );
	TH1F h_TOTfit_BKG_B0_M("TOTfit_BKG_B0_M", "", nbins,5., 5.5 );
	TH1F h_TOTfit_SGN_Rho_M("TOTfit_SGN_Rho_M", "", nbins,5., 5.5 );
	TH1F h_TOTfit_BKG_Rho_M("TOTfit_BKG_Rho_M", "", nbins,5., 5.5 );
	//.... pT/M
	TH1F h_TOTfit_SGN_B0_pT("TOTfit_SGN_B0_pT", "", nbins, 0, 20);
	TH1F h_TOTfit_BKG_B0_pT("TOTfit_BKG_B0_pT", "", nbins, 0, 20);
	TH1F h_TOTfit_SGN_JPsi_pT("TOTfit_SGN_JPsi_pT", "", nbins, 0, 1.);
	TH1F h_TOTfit_BKG_JPsi_pT("TOTfit_BKG_JPsi_pT", "", nbins, 0, 1.);
	TH1F h_TOTfit_SGN_Pi1_pT("TOTfit_SGN_Pi1_pT", "", nbins, 0, .3);
	TH1F h_TOTfit_BKG_Pi1_pT("TOTfit_BKG_Pi1_pT", "", nbins, 0, .3);
	TH1F h_TOTfit_SGN_Pi2_pT("TOTfit_SGN_Pi2_pT", "", nbins, 0, .2);
	TH1F h_TOTfit_BKG_Pi2_pT("TOTfit_BKG_Pi2_pT", "", nbins, 0, .2);
	TH1F h_TOTfit_SGN_Rho_pT("TOTfit_SGN_Rho_pT", "", nbins, 0, .5);
	TH1F h_TOTfit_BKG_Rho_pT("TOTfit_BKG_Rho_pT", "", nbins, 0, .5);
	TH1F h_TOTfit_SGN_X3872_pT("TOTfit_SGN_X3872_pT", "", nbins, 0, 1.);
	TH1F h_TOTfit_BKG_X3872_pT("TOTfit_BKG_X3872_pT", "", nbins, 0, 1.);
	TH1F h_TOTfit_SGN_K0s_pT("TOTfit_SGN_K0s_pT", "", nbins, 0, 0.6);
	TH1F h_TOTfit_BKG_K0s_pT("TOTfit_BKG_K0s_pT", "", nbins, 0, 0.6);

	//.... Pi-B0 DeltaR
	TH1F h_TOTfit_SGN_DR_Pi1B0_Rho("TOTfit_SGN_DR_Pi1B0_Rho", "", nbins,0, 1.);
	TH1F h_TOTfit_BKG_DR_Pi1B0_Rho("TOTfit_BKG_DR_Pi1B0_Rho", "", nbins,0, 1.);
	TH1F h_TOTfit_SGN_DR_Pi2B0_Rho("TOTfit_SGN_DR_Pi2B0_Rho", "", nbins,0, 1.);
	TH1F h_TOTfit_BKG_DR_Pi2B0_Rho("TOTfit_BKG_DR_Pi2B0_Rho", "", nbins,0, 1.);

	TH1F h_TOTfit_SGN_DR_Pi1B0_K0s("TOTfit_SGN_DR_Pi1B0_K0s", "", nbins,0, 1.);
	TH1F h_TOTfit_BKG_DR_Pi1B0_K0s("TOTfit_BKG_DR_Pi1B0_K0s", "", nbins,0, 1.);
	TH1F h_TOTfit_SGN_DR_Pi2B0_K0s("TOTfit_SGN_DR_Pi2B0_K0s", "", nbins,0, 1.);
	TH1F h_TOTfit_BKG_DR_Pi2B0_K0s("TOTfit_BKG_DR_Pi2B0_K0s", "", nbins,0, 1.);

	//.... K0s vertex
	TH1F h_TOTfit_SGN_K0s_SVp("TOTfit_SGN_K0s_SVp", "", nbins, 0., 1. );
	TH1F h_TOTfit_BKG_K0s_SVp("TOTfit_BKG_K0s_SVp", "", nbins, 0., 1. );

	//.... B0 vertex
	TH1F h_TOTfit_SGN_B0_LxySign("TOTfit_SGN_B0_LxySign", "", nbins, 0, 100);
	TH1F h_TOTfit_BKG_B0_LxySign("TOTfit_BKG_B0_LxySign", "", nbins, 0, 100);
	TH1F h_TOTfit_SGN_B0_SVchi2("TOTfit_SGN_B0_SVchi2", "", nbins, 0, 20);
	TH1F h_TOTfit_BKG_B0_SVchi2("TOTfit_BKG_B0_SVchi2", "", nbins, 0, 20);
	TH1F h_TOTfit_SGN_B0_SVp("TOTfit_SGN_B0_SVp", "", nbins, 0, 1.);
	TH1F h_TOTfit_BKG_B0_SVp("TOTfit_BKG_B0_SVp", "", nbins, 0, 1.);
	TH1F h_TOTfit_SGN_B0_cosA("TOTfit_SGN_B0_cosA", "", nbins, 0.95, 1.);
	TH1F h_TOTfit_BKG_B0_cosA("TOTfit_BKG_B0_cosA", "", nbins, 0.95, 1.);
	float rVtx;
	TH2F h_TOTfit_SGN_B0_rVSz_decayV("TOTfit_SGN_B0_rVSz_decayV", "", nbins, -10., 10., nbins, 0., 2.);
	TH2F h_TOTfit_BKG_B0_rVSz_decayV("TOTfit_BKG_B0_rVSz_decayV", "", nbins, -10., 10., nbins, 0., 2.);
	float Dphi;
	TH1F h_TOTfit_SGN_Dphi_B0RF("TOTfit_SGN_Dphi_B0RF", "", nbins, 3., 3.16);
	TH1F h_TOTfit_BKG_Dphi_B0RF("TOTfit_BKG_Dphi_B0RF", "", nbins, 3., 3.16);





   // Loop on B0 candidates 
   for (Long64_t jcand = 0; jcand < NB0cand; jcand++) {
		
		inTree_->GetEntry(jcand);		
		jentry = TriggerSel_event;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      fChain->GetEntry(jentry); 
      if (jcand + 1  == Nbreak) break;    
      if ((jcand + 1) % NPrint==0) cout << "--> " << (float)(jcand + 1)/NB0cand * 100. << " \%"<< endl;

		// ---> SIGNAL B0
		if (!(B0_SGN_idx < 0)){

			RecoPartFillP4(B0_SGN_idx);
			//... JPsi
			h_TOTfit_SGN_JPsi_pT.Fill(B0_MuMu_fitted_pt[B0_SGN_idx]/B0_finalFit_pt[B0_SGN_idx]);

			//... Pions
			h_TOTfit_SGN_Pi1_pT.Fill(P4_Reco_Pi1.Pt()/P4_Reco_B0.Pt());
			h_TOTfit_SGN_Pi2_pT.Fill(P4_Reco_Pi2.Pt()/P4_Reco_B0.Pt());

			//... Rho
			h_TOTfit_SGN_DR_Pi1B0_Rho.Fill(ROOT::Math::VectorUtil::DeltaR(P4_Reco_Pi1, P4_Reco_B0));
			h_TOTfit_SGN_DR_Pi2B0_Rho.Fill(ROOT::Math::VectorUtil::DeltaR(P4_Reco_Pi2, P4_Reco_B0));
			h_TOTfit_SGN_Rho_pT.Fill(P4_Reco_Rho.Pt()/B0_finalFit_pt[B0_SGN_idx]);
			h_TOTfit_SGN_Rho_M.Fill(B0_finalFit_Rho_mass[B0_SGN_idx]);	
			//... X(3872)
			h_TOTfit_SGN_X3872_pT.Fill(P4_Reco_X3872.Pt()/B0_finalFit_pt[B0_SGN_idx]);

			//... K0s
			h_TOTfit_SGN_DR_Pi1B0_K0s.Fill(ROOT::Math::VectorUtil::DeltaR(P4_Reco_K0trk1, P4_Reco_B0));
			h_TOTfit_SGN_DR_Pi2B0_K0s.Fill(ROOT::Math::VectorUtil::DeltaR(P4_Reco_K0trk2, P4_Reco_B0));
			h_TOTfit_SGN_K0s_pT.Fill(B0_K0s_mcFitted_pt[B0_SGN_idx]/B0_finalFit_pt[B0_SGN_idx]);
			h_TOTfit_SGN_K0s_SVp.Fill(B0_K0s_mcFitted_svprob[B0_SGN_idx]);

			//... B0
			h_TOTfit_SGN_B0_M.Fill(B0_finalFit_mass[B0_SGN_idx]);
			h_TOTfit_SGN_B0_pT.Fill(B0_finalFit_pt[B0_SGN_idx]/B0_finalFit_mass[B0_SGN_idx]);

			h_TOTfit_SGN_B0_LxySign.Fill(B0_lxySign_PV[B0_SGN_idx]);
			h_TOTfit_SGN_B0_SVchi2.Fill(B0_svchi2[B0_SGN_idx]);
			h_TOTfit_SGN_B0_SVp.Fill(B0_svprob[B0_SGN_idx]);
			h_TOTfit_SGN_B0_cosA.Fill(B0_cosAlpha_PV[B0_SGN_idx]);

			rVtx = sqrt(B0_decayVtxX[B0_SGN_idx]*B0_decayVtxX[B0_SGN_idx] + B0_decayVtxY[B0_SGN_idx]*B0_decayVtxY[B0_SGN_idx]);
			h_TOTfit_SGN_B0_rVSz_decayV.Fill( B0_decayVtxZ[B0_SGN_idx], rVtx);

			Dphi = DeltaPhi_RestFrameB0(B0_SGN_idx); 
			h_TOTfit_SGN_Dphi_B0RF.Fill(fabs(Dphi));
		}

		// ---> BACKGROUND B0
		for (UInt_t bb = 0; bb < nBKG_B0; bb++){
			Bidx = B0_BKG_idx[bb];
			RecoPartFillP4(Bidx);
			//... JPsi	
			if (!B0_BKG_isTrueJPsi[bb]) h_TOTfit_BKG_JPsi_pT.Fill(B0_MuMu_fitted_pt[Bidx]/B0_finalFit_pt[Bidx]);

			//... Rho
			if (!B0_BKG_isTrueRho[bb]){ 
				if (!B0_BKG_isTruePi1[bb]){
					h_TOTfit_BKG_DR_Pi1B0_Rho.Fill(ROOT::Math::VectorUtil::DeltaR(P4_Reco_Pi1, P4_Reco_B0));
					h_TOTfit_BKG_Pi1_pT.Fill(P4_Reco_Pi1.Pt()/P4_Reco_B0.Pt());
				}	
				if (!B0_BKG_isTruePi2[bb]){
					h_TOTfit_BKG_DR_Pi2B0_Rho.Fill(ROOT::Math::VectorUtil::DeltaR(P4_Reco_Pi2, P4_Reco_B0));
					h_TOTfit_BKG_Pi2_pT.Fill(P4_Reco_Pi2.Pt()/P4_Reco_B0.Pt());
				}
				h_TOTfit_BKG_Rho_pT.Fill(P4_Reco_Rho.Pt()/B0_finalFit_pt[Bidx]);
				h_TOTfit_BKG_Rho_M.Fill(B0_finalFit_Rho_mass[Bidx]);
			}	

			//... X(3872)
			if (!B0_BKG_isTrueJPsi[bb] || !B0_BKG_isTrueRho[bb]) h_TOTfit_BKG_X3872_pT.Fill(P4_Reco_X3872.Pt()/B0_finalFit_pt[Bidx]); 

			//.... K0s
			if (!B0_BKG_isTrueK0s[bb]){
				h_TOTfit_BKG_DR_Pi1B0_K0s.Fill(ROOT::Math::VectorUtil::DeltaR(P4_Reco_K0trk1, P4_Reco_B0));
				h_TOTfit_BKG_DR_Pi2B0_K0s.Fill(ROOT::Math::VectorUtil::DeltaR(P4_Reco_K0trk2, P4_Reco_B0));
				h_TOTfit_BKG_K0s_pT.Fill(B0_K0s_mcFitted_pt[Bidx]/B0_finalFit_pt[Bidx]);
				h_TOTfit_BKG_K0s_SVp.Fill(B0_K0s_mcFitted_svprob[Bidx]);
			}

			//... B0
			h_TOTfit_BKG_B0_M.Fill(B0_finalFit_mass[Bidx]);
			h_TOTfit_BKG_B0_pT.Fill(B0_finalFit_pt[Bidx]/B0_finalFit_mass[Bidx]);

			h_TOTfit_BKG_B0_LxySign.Fill(B0_lxySign_PV[Bidx]);
			h_TOTfit_BKG_B0_SVchi2.Fill(B0_svchi2[Bidx]);
			h_TOTfit_BKG_B0_SVp.Fill(B0_svprob[Bidx]);
			h_TOTfit_BKG_B0_cosA.Fill(B0_cosAlpha_PV[Bidx]);
			rVtx = sqrt(B0_decayVtxX[Bidx]*B0_decayVtxX[Bidx] + B0_decayVtxY[Bidx]*B0_decayVtxY[Bidx]);
			h_TOTfit_BKG_B0_rVSz_decayV.Fill( B0_decayVtxZ[Bidx], rVtx);
			Dphi = DeltaPhi_RestFrameB0(Bidx); 
			h_TOTfit_BKG_Dphi_B0RF.Fill(fabs(Dphi));

		}//on BKG B0
	}//on B0 candidates

	outFileHisto_ = new TFile("./plots/SB_variables.root", "RECREATE");
	h_TOTfit_SGN_B0_M.Write();
	h_TOTfit_BKG_B0_M.Write();

	h_TOTfit_SGN_JPsi_pT.Write();
	h_TOTfit_BKG_JPsi_pT.Write();
	h_TOTfit_SGN_Pi1_pT.Write();
	h_TOTfit_BKG_Pi1_pT.Write();
	h_TOTfit_SGN_Pi2_pT.Write();
	h_TOTfit_BKG_Pi2_pT.Write();
	h_TOTfit_SGN_Rho_pT.Write();
	h_TOTfit_BKG_Rho_pT.Write();
	h_TOTfit_SGN_X3872_pT.Write();
	h_TOTfit_BKG_X3872_pT.Write();
	h_TOTfit_SGN_K0s_pT.Write();
	h_TOTfit_BKG_K0s_pT.Write();
	h_TOTfit_SGN_B0_pT.Write();
	h_TOTfit_BKG_B0_pT.Write();

	h_TOTfit_SGN_DR_Pi1B0_Rho.Write();
	h_TOTfit_SGN_DR_Pi2B0_Rho.Write();
	h_TOTfit_BKG_DR_Pi1B0_Rho.Write();
	h_TOTfit_BKG_DR_Pi2B0_Rho.Write();
	h_TOTfit_SGN_DR_Pi1B0_K0s.Write();
	h_TOTfit_BKG_DR_Pi1B0_K0s.Write();
	h_TOTfit_SGN_DR_Pi2B0_K0s.Write();
	h_TOTfit_BKG_DR_Pi2B0_K0s.Write();

	h_TOTfit_SGN_K0s_SVp.Write();
	h_TOTfit_BKG_K0s_SVp.Write();


	h_TOTfit_SGN_B0_LxySign.Write();
	h_TOTfit_BKG_B0_LxySign.Write();
	h_TOTfit_SGN_B0_SVchi2.Write();
	h_TOTfit_BKG_B0_SVchi2.Write();
	h_TOTfit_SGN_B0_SVp.Write();
	h_TOTfit_BKG_B0_SVp.Write();
	h_TOTfit_SGN_B0_cosA.Write();
	h_TOTfit_BKG_B0_cosA.Write();
	h_TOTfit_SGN_B0_rVSz_decayV.Write();
	h_TOTfit_BKG_B0_rVSz_decayV.Write();
	h_TOTfit_SGN_Dphi_B0RF.Write();
	h_TOTfit_BKG_Dphi_B0RF.Write();

}//Loop()




int SGNvsBKGvariables::ReadTree(){

	int RETURN_VALUE = 0;

	inFileTree_ = new TFile(inFilePath_);
	if (!inFileTree_->IsOpen() ){
		std::cout << "ERROR : cannot open file "<< inFilePath_ << std::endl;
		RETURN_VALUE = 1;}
	std::cout << "Reading data from "<< inFilePath_ << "\n" << std::endl;

	inTree_ = (TTree*)inFileTree_->Get("B0_SGN_BKG");
	if ( !inTree_ ){
		std::cout<< "null pointer for TTree " << std::endl;
      exit(-1);}

	inTree_->SetBranchAddress("TriggerSel_event", &TriggerSel_event);	
	inTree_->SetBranchAddress("B0_SGN_idx", &B0_SGN_idx);
	inTree_->SetBranchAddress("nBKG_B0", &nBKG_B0);
	inTree_->SetBranchAddress("B0_BKG_idx", B0_BKG_idx);	
	inTree_->SetBranchAddress("B0_BKG_isTrueJPsi", B0_BKG_isTrueJPsi);
	inTree_->SetBranchAddress("B0_BKG_isTruePi1", B0_BKG_isTruePi1);
	inTree_->SetBranchAddress("B0_BKG_isTruePi2", B0_BKG_isTruePi2);
	inTree_->SetBranchAddress("B0_BKG_isTrueRho", B0_BKG_isTrueRho);
	inTree_->SetBranchAddress("B0_BKG_isTrueK0s", B0_BKG_isTrueK0s);

	return RETURN_VALUE; 

}//ReadTree()


int SGNvsBKGvariables::RecoPartFillP4(const int Bidx){
	int TrackQualityCheck = 1;

	//... muons P4
	if(!Muon_softId[B0_mu1_idx[Bidx]] || !Muon_softId[B0_mu2_idx[Bidx]]) TrackQualityCheck = 0; 
	P4_Reco_Mu1.SetPt(B0_finalFit_mu1_pt[Bidx]); P4_Reco_Mu1.SetEta(B0_finalFit_mu1_eta[Bidx]); P4_Reco_Mu1.SetPhi(B0_finalFit_mu1_phi[Bidx]);
	P4_Reco_Mu2.SetPt(B0_finalFit_mu2_pt[Bidx]); P4_Reco_Mu2.SetEta(B0_finalFit_mu2_eta[Bidx]); P4_Reco_Mu2.SetPhi(B0_finalFit_mu2_phi[Bidx]);
	//... JPsi P4
	P4_Reco_JPsi = P4_Reco_Mu1 + P4_Reco_Mu2; 
	//std::cout << "JPsi ... " << P4_Reco_JPsi.M() << "\t fit " << B0_finalFit_JPsi_mass[Bidx] << std::endl;
	//... Pi P4
	if(ProbeTracks_isMatchedToMuon[B0_pi1_idx[Bidx]] || ProbeTracks_isMatchedToMuon[B0_pi2_idx[Bidx]])TrackQualityCheck = 0;
	P4_Reco_Pi1.SetPt(B0_finalFit_pi1_pt[Bidx]); P4_Reco_Pi1.SetEta(B0_finalFit_pi1_eta[Bidx]); P4_Reco_Pi1.SetPhi(B0_finalFit_pi1_phi[Bidx]);
	P4_Reco_Pi2.SetPt(B0_finalFit_pi2_pt[Bidx]); P4_Reco_Pi2.SetEta(B0_finalFit_pi2_eta[Bidx]); P4_Reco_Pi2.SetPhi(B0_finalFit_pi2_phi[Bidx]);
	//... Rho P4
	P4_Reco_Rho = P4_Reco_Pi1 + P4_Reco_Pi2; //P4_Reco_Rho.SetM(B0_finalFit_Rho_mass[Bidx]);
	//... X(3872)
	P4_Reco_X3872 = P4_Reco_Rho + P4_Reco_JPsi; P4_Reco_X3872.SetM(B0_finalFit_X_mass[Bidx]);
	//std::cout << "X3872 ... " << P4_Reco_X3872.Pt() << "\t fit " << B0_finalFit_X_mass[Bidx] << std::endl;
	//... K0s track
	P4_Reco_K0trk1.SetPt(B0_K0s_nmcFitted_pi1pt[Bidx]); P4_Reco_K0trk1.SetEta(B0_K0s_nmcFitted_pi1eta[Bidx]); P4_Reco_K0trk1.SetPhi(B0_K0s_nmcFitted_pi1phi[Bidx]);	
	P4_Reco_K0trk2.SetPt(B0_K0s_nmcFitted_pi2pt[Bidx]); P4_Reco_K0trk2.SetEta(B0_K0s_nmcFitted_pi2eta[Bidx]); P4_Reco_K0trk2.SetPhi(B0_K0s_nmcFitted_pi2phi[Bidx]);	
	//... K0s
	P4_Reco_K0s.SetPt(B0_finalFit_k0s_pt[Bidx]); P4_Reco_K0s.SetEta(B0_finalFit_k0s_eta[Bidx]); P4_Reco_K0s.SetPhi(B0_finalFit_k0s_phi[Bidx]), P4_Reco_K0s.SetM(B0_K0s_mcFitted_mass[Bidx]);
	//... B0
	P4_Reco_B0.SetPt(B0_finalFit_pt[Bidx]); P4_Reco_B0.SetEta(B0_finalFit_eta[Bidx]); P4_Reco_B0.SetPhi(B0_finalFit_phi[Bidx]), P4_Reco_B0.SetM(B0_finalFit_mass[Bidx]);
	 
	return TrackQualityCheck;

}//RecoPartFillP4()


float SGNvsBKGvariables::DeltaPhi_RestFrameB0(const int Bidx){

	float Dphi;
	TLorentzVector LV_Mu1, LV_Mu2, LV_Pi1, LV_Pi2;
	LV_Mu1.SetPtEtaPhiM(B0_MuMu_prefit_mu1_pt[Bidx], B0_MuMu_prefit_mu1_eta[Bidx], B0_MuMu_prefit_mu1_phi[Bidx], mMuon);
	LV_Mu2.SetPtEtaPhiM(B0_MuMu_prefit_mu2_pt[Bidx], B0_MuMu_prefit_mu2_eta[Bidx], B0_MuMu_prefit_mu2_phi[Bidx], mMuon);
	LV_Pi1.SetPtEtaPhiM(B0_PiPi_prefit_pi1_pt[Bidx], B0_PiPi_prefit_pi1_eta[Bidx], B0_PiPi_prefit_pi1_phi[Bidx], mPion);
	LV_Pi2.SetPtEtaPhiM(B0_PiPi_prefit_pi2_pt[Bidx], B0_PiPi_prefit_pi2_eta[Bidx], B0_PiPi_prefit_pi2_phi[Bidx], mPion);
	TLorentzVector LV_Boost_X3872, LV_Boost_K0s, LV_Boost_B0;

	LV_Boost_X3872 = LV_Mu1+LV_Mu2 + LV_Pi1+LV_Pi2; 
	LV_Boost_K0s.SetPtEtaPhiM(B0_K0s_mcFitted_pt[Bidx], B0_K0s_mcFitted_eta[Bidx], B0_K0s_mcFitted_phi[Bidx], B0_K0s_nmcFitted_mass[Bidx]);
	LV_Boost_B0.SetPtEtaPhiM(P4_Reco_B0.Pt(), P4_Reco_B0.Eta(), P4_Reco_B0.Phi(), P4_Reco_B0.M());

	TVector3 BoostVectB0 = -LV_Boost_B0.BoostVector();
	LV_Boost_B0.Boost(BoostVectB0);
	//std::cout << "B0 "<< LV_Boost_B0.Px() << "\t" << LV_Boost_B0.Py() << "\t" << LV_Boost_B0.Pz() << "\t" << LV_Boost_B0.M() << std::endl;
	LV_Boost_X3872.Boost(BoostVectB0);	
	//std::cout << "X " << LV_Boost_X3872.Px() << "\t" << LV_Boost_X3872.Py() << "\t" << LV_Boost_X3872.Pz() << "\t" << LV_Boost_X3872.M() << std::endl;
	LV_Boost_K0s.Boost(BoostVectB0);
	//std::cout << "K0s "<< LV_Boost_K0s.Px() << "\t" << LV_Boost_K0s.Py() << "\t" << LV_Boost_K0s.Pz() << "\t" << LV_Boost_K0s.M() << std::endl;

	Dphi = ROOT::Math::VectorUtil::DeltaPhi(LV_Boost_X3872, LV_Boost_K0s);

	return Dphi;

}//DeltaPhi_RestFrameB0()
