#include "../include/SGNvsBKGvariables.h" 

// CONSTRUCTOR-DESTRUCTOR
SGNvsBKGvariables::SGNvsBKGvariables(TTree *tree) : B0toX3872K0s(tree){ 
	
	inFilePath_ = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/B0_SGNvsBKG_anlaysis.root";
	outFilePath_ = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/TMVAinputs.root";

   P4_Reco_Mu1.SetM(mMuon); P4_Reco_Mu2.SetM(mMuon);
   P4_Reco_Pi1.SetM(mPion); P4_Reco_Pi2.SetM(mPion);
   P4_Reco_K0s.SetM(mK0s);
}

SGNvsBKGvariables::~SGNvsBKGvariables(){
	inFileTree_->Close();
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
	nbins = 100;
	float Mlow = 5., Mhigh = 5.6;
	TH1F h_SGN_B0_M("SGN_B0_M", "", nbins+20, Mlow, Mhigh );
	TH1F h_BKG_B0_M("BKG_B0_M", "", nbins,Mlow, Mhigh );
	TH1F h_BKGr_B0_M("BKGr_B0_M", "", nbins,Mlow, Mhigh );
	TH1F h_BKGk_B0_M("BKGk_B0_M", "", nbins,Mlow, Mhigh );
	TH1F h_BKGrk_B0_M("BKGrk_B0_M", "", nbins,Mlow, Mhigh );
	
	Mlow = .4, Mhigh = 1.;
	TH1F h_SGN_Rho_M("SGN_Rho_M", "", nbins, Mlow, Mhigh );
	TH1F h_BKG_Rho_M("BKG_Rho_M", "", nbins, Mlow, Mhigh );
	TH1F h_BKGb_Rho_M("BKGb_Rho_M", "", nbins,Mlow, Mhigh );

	Mlow = 3.7, Mhigh = 4.05;	
	TH1F h_SGN_X3872_M("SGN_X3872_M", "", nbins, Mlow, Mhigh );
	Mlow = 0.38, Mhigh = 0.62;	
	TH1F h_SGN_K0s_M("SGN_K0s_M", "", nbins + 20, Mlow, Mhigh );

	TH2F h_SGN_RhovsK0s_M("SGN_RhovsK0s_M", "", nbins,  .45, .55, nbins, .4, 1.);
	TH2F h_BKGk_RhovsK0s_M("BKGk_RhovsK0s_M", "", nbins,  .45, .55 ,nbins, .4, 1.);
	TH2F h_BKGr_RhovsK0s_M("BKGr_RhovsK0s_M", "", nbins,  .45, .55 ,nbins, .4, 1.);
	TH2F h_BKGrk_RhovsK0s_M("BKGrk_RhovsK0s_M", "", nbins,  .45, .55, nbins, .4, 1.);
	

	TH2F h_SGN_B0vsX_M("SGN_B0vsX_M", "", nbins, 3.6, 4.1, nbins, 5., 5.5);
	TH2F h_BKGk_B0vsX_M("BKGk_B0vsX_M", "", nbins, 3.6, 4.1, nbins, 5., 5.5);
	TH2F h_BKG_B0vsX_M("BKG_B0vsX_M", "", nbins, 3.6, 4.1, nbins, 5., 5.5);
	
	nbins = 50;
	//.... pT/M

	TH1F h_SGN_B0_pT("SGN_B0_pT", "", nbins, 0, 20);
	TH1F h_BKG_B0_pT("BKG_B0_pT", "", 40, 2.0439124, 23.200866);
	//TH1F h_BKG_B0_pT("BKG_B0_pT", "", nbins, 0, 20);

	TH1F h_SGN_JPsi_pT("SGN_JPsi_pT", "", nbins, 0, 1.);
	TH1F h_BKG_JPsi_pT("BKG_JPsi_pT", "", nbins, 0, 1.);

	TH1F h_SGN_Pi1_pT("SGN_Pi1_pT", "", nbins, 0, .3);
	TH1F h_BKG_Pi1_pT("BKG_Pi1_pT", "", nbins, 0, .3);
	TH1F h_BKGb_Pi1_pT("BKGb_Pi1_pT", "",40, 0.016485984, 0.28147935);
	TH1F h_SGN_Pi2_pT("SGN_Pi2_pT", "", nbins, 0, .2);
	TH1F h_BKG_Pi2_pT("BKG_Pi2_pT", "", nbins, 0, .2);

	TH1F h_SGN_Rho_pT("SGN_Rho_pT", "", nbins, 0, .5);
	TH1F h_BKG_Rho_pT("BKG_Rho_pT", "", nbins, 0, .5);
	TH1F h_BKGb_Rho_pT("BKGb_Rho_pT", "", 40, 0.028098026, 0.38227733);

	TH1F h_SGN_X3872_pT("SGN_X3872_pT", "", nbins, 0.4, 1.);
	TH1F h_BKG_X3872_pT("BKG_X3872_pT", "", nbins, 0.4, 1.);

	TH1F h_SGN_K0s_pT("SGN_K0s_pT", "", nbins, 0, 0.6);
	TH1F h_BKG_K0s_pT("BKG_K0s_pT", "", nbins, 0, 0.6);
	TH1F h_BKGb_K0s_pT("BKGb_K0s_pT", "", nbins, 0, 0.6);
	
	TH1F h_SGN_LeadTrk_pT("SGN_LeadTrk_pT", "", nbins, 0., .5);
	TH1F h_BKG_LeadTrk_pT("BKG_LeadTrk_pT", "", nbins, 0., .5);

	//.... Pi-B0 DeltaR
	TH1F h_SGN_DR_Pi1B0_Rho("SGN_DR_Pi1B0_Rho", "", nbins,0, 1.);
	TH1F h_BKG_DR_Pi1B0_Rho("BKG_DR_Pi1B0_Rho", "", nbins,0, 1.);
	TH1F h_BKGb_DR_Pi1B0_Rho("BKGb_DR_Pi1B0_Rho", "", 40 ,0.0033471850, 0.74371076);
	TH1F h_SGN_DR_Pi2B0_Rho("SGN_DR_Pi2B0_Rho", "", nbins,0, 1.);
	TH1F h_BKG_DR_Pi2B0_Rho("BKG_DR_Pi2B0_Rho", "", nbins,0, 1.);

	TH1F h_SGN_DR_Pi1B0_K0s("SGN_DR_Pi1B0_K0s", "", nbins,0, 1.);
	TH1F h_BKG_DR_Pi1B0_K0s("BKG_DR_Pi1B0_K0s", "", nbins,0, 1.);
	TH1F h_SGN_DR_Pi2B0_K0s("SGN_DR_Pi2B0_K0s", "", nbins,0, 1.);
	TH1F h_BKG_DR_Pi2B0_K0s("BKG_DR_Pi2B0_K0s", "", nbins,0, 1.);

	//.... Rho vertex
	TH1F h_SGN_Rho_D0("SGN_Rho_D0", "", nbins, 0., 8.);
	TH1F h_BKG_Rho_D0("BKG_Rho_D0", "", nbins, 0., 8.);
	TH1F h_BKGb_Rho_D0("BKGb_Rho_D0", "", 40, 0.00046771113, 57.158501);
	TH1F h_SGN_Rho_D0max("SGN_Rho_D0max", "", nbins, 0., .3);
	TH1F h_BKG_Rho_D0max("BKG_Rho_D0max", "", nbins, 0., .3);
	TH1F h_BKGb_Rho_D0max("BKGb_Rho_D0max", "", nbins, 0., .3);
	TH1F h_SGN_Rho_D0min("SGN_Rho_D0min", "", nbins, 0., .01);
	TH1F h_BKG_Rho_D0min("BKG_Rho_D0min", "", nbins, 0., 0.01);
	TH1F h_BKGb_Rho_D0min("BKGb_Rho_D0min", "", nbins, 0., 0.01);

	//.... K0s vertex
	TH1F h_SGN_K0s_SVp("SGN_K0s_SVp", "", nbins, 0., 1. );
	TH1F h_BKG_K0s_SVp("BKG_K0s_SVp", "", nbins, 0., 1. );
	TH1F h_SGN_K0s_D0("SGN_K0s_D0", "", nbins, 0., 8.);
	TH1F h_BKG_K0s_D0("BKG_K0s_D0", "", nbins, 0., 8. );
	TH1F h_BKGb_K0s_D0("BKGb_K0s_D0", "", nbins, 0., 8. );
	TH2F h_SGN_K0s_D0vsPt("SGN_K0s_D0vsPt", "", nbins, 0, 0.6, nbins, 0., 8.);
	TH2F h_BKG_K0s_D0vsPt("BKG_K0s_D0vsPt", "", nbins, 0, 0.6, nbins, 0., 8.);

	//.... B0 vertex
	nbins = 40;
	TH1F h_SGN_B0_LxySign("SGN_B0_LxySign", "", nbins,0, 100); 
	TH1F h_BKG_B0_LxySign("BKG_B0_LxySign", "", nbins, 3.1187811, 464.68721);
	TH1F h_SGN_B0_SVchi2("SGN_B0_SVchi2", "", nbins, 0, 20);
	TH1F h_BKG_B0_SVchi2("BKG_B0_SVchi2", "", nbins, 0, 20);
	TH1F h_SGN_B0_SVp("SGN_B0_SVp", "", nbins, 0, 1.);
	TH1F h_BKG_B0_SVp("BKG_B0_SVp", "", nbins, 0.010026554, 1.0245656);
	TH1F h_SGN_B0_cosA("SGN_B0_cosA", "", nbins, 0.95, 1.);
	TH1F h_BKG_B0_cosA("BKG_B0_cosA", "", nbins, 0.96574142, 1.0008565);
	nbins = 50;
	float rVtx;
	TH2F h_SGN_B0_rVSz_decayV("SGN_B0_rVSz_decayV", "", nbins, -10., 10., nbins, 0., 2.);
	TH2F h_BKG_B0_rVSz_decayV("BKG_B0_rVSz_decayV", "", nbins, -10., 10., nbins, 0., 2.);
	float Dphi;
	TH1F h_SGN_Dphi_B0RF("SGN_Dphi_B0RF", "", nbins, 3., 3.16);
	TH1F h_BKG_Dphi_B0RF("BKG_Dphi_B0RF", "", nbins, 3., 3.16);


	// out TTree setup
	BookTreeTMVA();


   // Loop on B0 candidates 
   for (Long64_t jcand = 0; jcand < NB0cand; jcand++) {
		
		inTree_->GetEntry(jcand);		
		jentry = TriggerSel_event;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      fChain->GetEntry(jentry); 
      if (jcand + 1  == Nbreak) break;    
      if ((jcand + 1) % NPrint==0) cout << "--> " << Form( "%3.0f", (float)(jcand + 1)/NB0cand * 100.) << " \%"<< endl;

	// ==================== SIGNAL B0 ====================//
		if (!(B0_SGN_idx < 0)){

			RecoPartFillP4(B0_SGN_idx);

			//... JPsi
			h_SGN_JPsi_pT.Fill(B0_MuMu_fitted_pt[B0_SGN_idx]/B0_finalFit_pt[B0_SGN_idx]);

			//... Pions
			h_SGN_Pi1_pT.Fill(P4_Reco_Pi1.Pt()/P4_Reco_B0.Pt());
			h_SGN_Pi2_pT.Fill(P4_Reco_Pi2.Pt()/P4_Reco_B0.Pt());

			h_SGN_DR_Pi1B0_Rho.Fill(ROOT::Math::VectorUtil::DeltaR(P4_Reco_Pi1, P4_Reco_B0));
			h_SGN_DR_Pi2B0_Rho.Fill(ROOT::Math::VectorUtil::DeltaR(P4_Reco_Pi2, P4_Reco_B0));

			//... Rho
			h_SGN_Rho_D0.Fill(B0_PiPi_pi1_d0sig[B0_SGN_idx]);
			h_SGN_Rho_D0.Fill(B0_PiPi_pi2_d0sig[B0_SGN_idx]);
			h_SGN_Rho_D0max.Fill(B0_PiPi_pi1_maxd0PV[B0_SGN_idx]);
			h_SGN_Rho_D0min.Fill(B0_PiPi_pi1_mind0PV[B0_SGN_idx]);

			h_SGN_Rho_pT.Fill(P4_Reco_Rho.Pt()/B0_finalFit_pt[B0_SGN_idx]);
			h_SGN_Rho_M.Fill(B0_finalFit_Rho_mass[B0_SGN_idx]);	

			//... X(3872)
			h_SGN_X3872_M.Fill(P4_Reco_X3872.M());
			h_SGN_X3872_pT.Fill(P4_Reco_X3872.Pt()/B0_finalFit_pt[B0_SGN_idx]);

			//... K0s
			h_SGN_K0s_M.Fill(B0_K0s_nmcFitted_mass[B0_SGN_idx]);
			h_SGN_DR_Pi1B0_K0s.Fill(ROOT::Math::VectorUtil::DeltaR(P4_Reco_K0trk1, P4_Reco_B0));
			h_SGN_DR_Pi2B0_K0s.Fill(ROOT::Math::VectorUtil::DeltaR(P4_Reco_K0trk2, P4_Reco_B0));

			h_SGN_K0s_pT.Fill(B0_K0s_mcFitted_pt[B0_SGN_idx]/B0_finalFit_pt[B0_SGN_idx]);

			h_SGN_K0s_SVp.Fill(B0_K0s_mcFitted_svprob[B0_SGN_idx]);

			h_SGN_K0s_D0.Fill(B0_K0s_matchTrack1_D0sign[B0_SGN_idx]);
			h_SGN_K0s_D0.Fill(B0_K0s_matchTrack2_D0sign[B0_SGN_idx]);

			h_SGN_K0s_D0vsPt.Fill(B0_finalFit_k0s_pt[B0_SGN_idx]/B0_finalFit_pt[B0_SGN_idx], B0_K0s_matchTrack1_D0sign[B0_SGN_idx]);
			h_SGN_K0s_D0vsPt.Fill(B0_finalFit_k0s_pt[B0_SGN_idx]/B0_finalFit_pt[B0_SGN_idx], B0_K0s_matchTrack2_D0sign[B0_SGN_idx]);

			h_SGN_LeadTrk_pT.Fill(max(P4_Reco_K0trk1.Pt(),P4_Reco_Pi1.Pt())/B0_finalFit_pt[B0_SGN_idx]);

			//... B0
			h_SGN_B0_M.Fill(B0_finalFit_mass[B0_SGN_idx]);
			h_SGN_B0_pT.Fill(B0_finalFit_pt[B0_SGN_idx]/B0_finalFit_mass[B0_SGN_idx]);

			h_SGN_B0_LxySign.Fill(B0_lxySign_PV[B0_SGN_idx]);
			h_SGN_B0_SVchi2.Fill(B0_svchi2[B0_SGN_idx]);
			h_SGN_B0_SVp.Fill(B0_svprob[B0_SGN_idx]);
			h_SGN_B0_cosA.Fill(B0_cosAlpha_PV[B0_SGN_idx]);

			rVtx = sqrt(B0_decayVtxX[B0_SGN_idx]*B0_decayVtxX[B0_SGN_idx] + B0_decayVtxY[B0_SGN_idx]*B0_decayVtxY[B0_SGN_idx]);
			h_SGN_B0_rVSz_decayV.Fill( B0_decayVtxZ[B0_SGN_idx], rVtx);

			Dphi = DeltaPhi_RestFrameB0(B0_SGN_idx); 
			h_SGN_Dphi_B0RF.Fill(fabs(Dphi));


			//... 2D masses
			h_SGN_B0vsX_M.Fill(B0_finalFit_X_mass[B0_SGN_idx] , B0_finalFit_mass[B0_SGN_idx]);
			h_SGN_RhovsK0s_M.Fill(B0_K0s_nmcFitted_mass[B0_SGN_idx], B0_finalFit_Rho_mass[B0_SGN_idx]);



		// ================== TMVA variables ================== //

			pTM_B0_S = B0_finalFit_pt[B0_SGN_idx]/B0_finalFit_mass[B0_SGN_idx];
			SVprob_S = B0_svprob[B0_SGN_idx]; 
			SVchi2_S = B0_svchi2[B0_SGN_idx];
			LxySign_B0_S = B0_lxySign_PV[B0_SGN_idx];
			CosAlpha_B0_S = B0_cosAlpha_PV[B0_SGN_idx];

			DR_Pi1B0_S = ROOT::Math::VectorUtil::DeltaR(P4_Reco_Pi1, P4_Reco_B0);
			pT_Pi1_S = P4_Reco_Pi1.Pt()/P4_Reco_B0.Pt(); 
			pT_Rho_S = P4_Reco_Rho.Pt()/B0_finalFit_pt[B0_SGN_idx];
			D0_Rho_S = B0_PiPi_pi1_d0sig[B0_SGN_idx]; 

			M_B0_S    = B0_finalFit_mass[B0_SGN_idx];
			M_mumu_S  = B0_MuMu_fitted_mass[B0_SGN_idx];
			M_Rho_S   = B0_finalFit_Rho_mass[B0_SGN_idx];
			M_X3872_S = B0_finalFit_X_mass[B0_SGN_idx];
			M_K0s_S   = B0_K0s_nmcFitted_mass[B0_SGN_idx];

			TMVAoutTreeSGN_->Fill();
			
		}

	// ==================== BACKGROUND B0 ====================//

		for (UInt_t bb = 0; bb < nBKG_B0; bb++){
			Bidx = B0_BKG_idx[bb];
			RecoPartFillP4(Bidx);

			//... JPsi	

			//... Rho
			if (!B0_BKG_isTrueRho[bb]){ 

				h_BKG_Rho_D0.Fill(B0_PiPi_pi1_d0sig[Bidx]);
				h_BKG_Rho_D0.Fill(B0_PiPi_pi2_d0sig[Bidx]);
				h_BKG_Rho_pT.Fill(P4_Reco_Rho.Pt()/B0_finalFit_pt[Bidx]);
				h_BKG_Rho_M.Fill(B0_finalFit_Rho_mass[Bidx]);

				if (!B0_BKG_isTruePi1[bb]){
					h_BKG_DR_Pi1B0_Rho.Fill(ROOT::Math::VectorUtil::DeltaR(P4_Reco_Pi1, P4_Reco_B0));
					h_BKG_Pi1_pT.Fill(P4_Reco_Pi1.Pt()/P4_Reco_B0.Pt());
					h_BKG_Rho_D0max.Fill(B0_PiPi_pi1_maxd0PV[Bidx]);
					h_BKG_Rho_D0min.Fill(B0_PiPi_pi1_mind0PV[Bidx]);
				}	

				if (!B0_BKG_isTruePi2[bb]){
					h_BKG_DR_Pi2B0_Rho.Fill(ROOT::Math::VectorUtil::DeltaR(P4_Reco_Pi2, P4_Reco_B0));
					h_BKG_Pi2_pT.Fill(P4_Reco_Pi2.Pt()/P4_Reco_B0.Pt());
				}

			}else{

				h_SGN_Rho_D0.Fill(B0_PiPi_pi1_d0sig[Bidx]);
				h_SGN_Rho_D0.Fill(B0_PiPi_pi2_d0sig[Bidx]);
				h_SGN_Rho_pT.Fill(P4_Reco_Rho.Pt()/B0_finalFit_pt[Bidx]);
				h_SGN_Rho_M.Fill(B0_finalFit_Rho_mass[Bidx]);

			}	

			//... X(3872)
			if (!B0_BKG_isTrueJPsi[bb] || !B0_BKG_isTrueRho[bb]) h_BKG_X3872_pT.Fill(P4_Reco_X3872.Pt()/B0_finalFit_pt[Bidx]); 

			//.... K0s
			if (!B0_BKG_isTrueK0s[bb]){

				h_BKG_DR_Pi1B0_K0s.Fill(ROOT::Math::VectorUtil::DeltaR(P4_Reco_K0trk1, P4_Reco_B0));
				h_BKG_DR_Pi2B0_K0s.Fill(ROOT::Math::VectorUtil::DeltaR(P4_Reco_K0trk2, P4_Reco_B0));
				h_BKG_K0s_D0.Fill(B0_K0s_matchTrack1_D0sign[Bidx]);
				h_BKG_K0s_D0.Fill(B0_K0s_matchTrack2_D0sign[Bidx]);
				h_BKG_K0s_pT.Fill(B0_finalFit_k0s_pt[Bidx]/B0_finalFit_pt[Bidx]);
				h_BKG_K0s_SVp.Fill(B0_K0s_mcFitted_svprob[Bidx]);

				h_BKG_K0s_D0vsPt.Fill(B0_finalFit_k0s_pt[Bidx]/B0_finalFit_pt[Bidx], B0_K0s_matchTrack1_D0sign[Bidx]);
				h_BKG_K0s_D0vsPt.Fill(B0_finalFit_k0s_pt[Bidx]/B0_finalFit_pt[Bidx], B0_K0s_matchTrack2_D0sign[Bidx]);
			}else{
			
				h_SGN_DR_Pi1B0_K0s.Fill(ROOT::Math::VectorUtil::DeltaR(P4_Reco_K0trk1, P4_Reco_B0));
				h_SGN_DR_Pi2B0_K0s.Fill(ROOT::Math::VectorUtil::DeltaR(P4_Reco_K0trk2, P4_Reco_B0));
				h_SGN_K0s_D0.Fill(B0_K0s_matchTrack1_D0sign[Bidx]);
				h_SGN_K0s_D0.Fill(B0_K0s_matchTrack2_D0sign[Bidx]);
				h_SGN_K0s_pT.Fill(B0_K0s_mcFitted_pt[Bidx]/B0_finalFit_pt[Bidx]);
				h_SGN_K0s_SVp.Fill(B0_K0s_mcFitted_svprob[Bidx]);
			}

			//... B0
			//mass
			h_BKG_B0_M.Fill(B0_finalFit_mass[Bidx]);
			if(!B0_BKG_isTrueRho[bb] &&  B0_BKG_isTrueK0s[bb]) h_BKGr_B0_M.Fill(B0_finalFit_mass[Bidx]);
			if( B0_BKG_isTrueRho[bb] && !B0_BKG_isTrueK0s[bb]) h_BKGk_B0_M.Fill(B0_finalFit_mass[Bidx]);
			if(!B0_BKG_isTrueRho[bb] && !B0_BKG_isTrueK0s[bb]) h_BKGrk_B0_M.Fill(B0_finalFit_mass[Bidx]);


			h_BKG_B0_pT.Fill(B0_finalFit_pt[Bidx]/B0_finalFit_mass[Bidx]);

			h_BKG_B0_LxySign.Fill(B0_lxySign_PV[Bidx]);
			h_BKG_B0_SVchi2.Fill(B0_svchi2[Bidx]);
			h_BKG_B0_SVp.Fill(B0_svprob[Bidx]);
			h_BKG_B0_cosA.Fill(B0_cosAlpha_PV[Bidx]);

			rVtx = sqrt(B0_decayVtxX[Bidx]*B0_decayVtxX[Bidx] + B0_decayVtxY[Bidx]*B0_decayVtxY[Bidx]);
			h_BKG_B0_rVSz_decayV.Fill( B0_decayVtxZ[Bidx], rVtx);

			Dphi = DeltaPhi_RestFrameB0(Bidx); 
			h_BKG_Dphi_B0RF.Fill(fabs(Dphi));

			// ... 2D masses
			if(!B0_BKG_isTrueRho[bb] &&  B0_BKG_isTrueK0s[bb]) h_BKGr_RhovsK0s_M.Fill(B0_K0s_nmcFitted_mass[Bidx], B0_finalFit_Rho_mass[Bidx]);		
			if( B0_BKG_isTrueRho[bb] && !B0_BKG_isTrueK0s[bb]){
				h_BKGk_B0vsX_M.Fill(B0_finalFit_X_mass[Bidx] , B0_finalFit_mass[Bidx]);
				h_BKGk_RhovsK0s_M.Fill(B0_K0s_nmcFitted_mass[Bidx], B0_finalFit_Rho_mass[Bidx]);		
				h_BKG_LeadTrk_pT.Fill(max(P4_Reco_K0trk1.Pt(),P4_Reco_Pi1.Pt())/B0_finalFit_pt[Bidx]);
			}

			if(!B0_BKG_isTrueRho[bb] && !B0_BKG_isTrueK0s[bb]){
				h_BKGrk_RhovsK0s_M.Fill(B0_K0s_nmcFitted_mass[Bidx], B0_finalFit_Rho_mass[Bidx]);		
				h_BKG_LeadTrk_pT.Fill(max(P4_Reco_K0trk1.Pt(),P4_Reco_Pi1.Pt())/B0_finalFit_pt[Bidx]);
			}
			// --> DISCRIMINANT

			h_BKGb_DR_Pi1B0_Rho.Fill(ROOT::Math::VectorUtil::DeltaR(P4_Reco_Pi1, P4_Reco_B0));
			h_BKGb_Pi1_pT.Fill(P4_Reco_Pi1.Pt()/P4_Reco_B0.Pt());

			h_BKGb_Rho_D0.Fill(B0_PiPi_pi1_d0sig[Bidx]);
			h_BKGb_Rho_D0.Fill(B0_PiPi_pi2_d0sig[Bidx]);
			h_BKGb_Rho_D0max.Fill(B0_PiPi_pi1_maxd0PV[Bidx]);
			h_BKGb_Rho_D0min.Fill(B0_PiPi_pi1_mind0PV[Bidx]);
			h_BKGb_Rho_pT.Fill(P4_Reco_Rho.Pt()/B0_finalFit_pt[Bidx]);
			h_BKGb_Rho_M.Fill(B0_finalFit_Rho_mass[Bidx]);

			h_BKGb_K0s_pT.Fill(B0_K0s_mcFitted_pt[Bidx]/B0_finalFit_pt[Bidx]);
			h_BKGb_K0s_D0.Fill(B0_K0s_matchTrack1_D0sign[Bidx]);
			h_BKGb_K0s_D0.Fill(B0_K0s_matchTrack2_D0sign[Bidx]);

		// ================== TMVA variables ================== //
			
			pTM_B0_B = B0_finalFit_pt[Bidx]/B0_finalFit_mass[Bidx];
			SVprob_B = B0_svprob[Bidx]; 
			SVchi2_B = B0_svchi2[Bidx]; 
			LxySign_B0_B = B0_lxySign_PV[Bidx];
			CosAlpha_B0_B = B0_cosAlpha_PV[Bidx];

			DR_Pi1B0_B = ROOT::Math::VectorUtil::DeltaR(P4_Reco_Pi1, P4_Reco_B0);
			pT_Pi1_B = P4_Reco_Pi1.Pt()/P4_Reco_B0.Pt(); 
			pT_Rho_B = P4_Reco_Rho.Pt()/B0_finalFit_pt[Bidx];
			D0_Rho_B = B0_PiPi_pi1_d0sig[Bidx]; 

			M_B0_B    = B0_finalFit_mass[Bidx];
			M_mumu_B  = B0_MuMu_fitted_mass[Bidx];
			M_Rho_B   = B0_finalFit_Rho_mass[Bidx];
			M_X3872_B = B0_finalFit_X_mass[Bidx];
			M_K0s_B   = B0_K0s_nmcFitted_mass[Bidx];

			TMVAoutTreeBKG_->Fill();

		}//on BKG B0



	}//on EVENTS 
	

	TMVAoutFile_->cd();

	TMVAoutTreeSGN_->Write();
	TMVAoutTreeBKG_->Write();

	TMVAoutFile_->Close();

	outFileHisto_ = new TFile("./plots/SB_variables.root", "RECREATE");
	h_SGN_B0_M.Write();
	h_BKG_B0_M.Write();
	h_BKGr_B0_M.Write();
	h_BKGk_B0_M.Write();
	h_BKGrk_B0_M.Write();

	h_SGN_Rho_M.Write();
	h_BKG_Rho_M.Write();
	h_SGN_X3872_M.Write();
	h_SGN_K0s_M.Write();

	h_SGN_B0vsX_M.Write();
	h_BKGk_B0vsX_M.Write();
	h_BKG_B0vsX_M.Write();

	h_SGN_RhovsK0s_M.Write();	
	h_BKGr_RhovsK0s_M.Write();	
	h_BKGk_RhovsK0s_M.Write();	
	h_BKGrk_RhovsK0s_M.Write();	

	h_SGN_JPsi_pT.Write();
	h_BKG_JPsi_pT.Write();
	h_SGN_Pi1_pT.Write();
	h_BKG_Pi1_pT.Write();
	h_SGN_Pi2_pT.Write();
	h_BKG_Pi2_pT.Write();
	h_SGN_Rho_pT.Write();
	h_BKG_Rho_pT.Write();
	h_SGN_X3872_pT.Write();
	h_BKG_X3872_pT.Write();
	h_SGN_K0s_pT.Write();
	h_BKG_K0s_pT.Write();
	h_SGN_LeadTrk_pT.Write();
	h_BKG_LeadTrk_pT.Write();
	h_SGN_B0_pT.Write();
	h_BKG_B0_pT.Write();

	h_SGN_DR_Pi1B0_Rho.Write();
	h_SGN_DR_Pi2B0_Rho.Write();
	h_BKG_DR_Pi1B0_Rho.Write();
	h_BKG_DR_Pi2B0_Rho.Write();
	h_SGN_DR_Pi1B0_K0s.Write();
	h_BKG_DR_Pi1B0_K0s.Write();
	h_SGN_DR_Pi2B0_K0s.Write();
	h_BKG_DR_Pi2B0_K0s.Write();

	h_SGN_Rho_D0.Write();
	h_BKG_Rho_D0.Write();
	h_SGN_Rho_D0max.Write();
	h_BKG_Rho_D0max.Write();
	h_SGN_Rho_D0min.Write();
	h_BKG_Rho_D0min.Write();

	h_SGN_K0s_SVp.Write();
	h_BKG_K0s_SVp.Write();
	h_SGN_K0s_D0.Write();
	h_BKG_K0s_D0.Write();
	h_SGN_K0s_D0vsPt.Write();
	h_BKG_K0s_D0vsPt.Write();


	h_SGN_B0_LxySign.Write();
	h_BKG_B0_LxySign.Write();
	h_SGN_B0_SVchi2.Write();
	h_BKG_B0_SVchi2.Write();
	h_SGN_B0_SVp.Write();
	h_BKG_B0_SVp.Write();
	h_SGN_B0_cosA.Write();
	h_BKG_B0_cosA.Write();
	h_SGN_B0_rVSz_decayV.Write();
	h_BKG_B0_rVSz_decayV.Write();
	h_SGN_Dphi_B0RF.Write();
	h_BKG_Dphi_B0RF.Write();


	h_BKGb_DR_Pi1B0_Rho.Write();
	h_BKGb_Pi1_pT.Write();

	h_BKGb_Rho_D0.Write();
	h_BKGb_Rho_D0max.Write();
	h_BKGb_Rho_D0min.Write();
	h_BKGb_Rho_pT.Write();
	h_BKGb_Rho_M.Write();

	h_BKGb_K0s_pT.Write();
	h_BKGb_K0s_D0.Write();

	outFileHisto_->Close();


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


void SGNvsBKGvariables::BookTreeTMVA(){

	TMVAoutFile_ = new TFile(outFilePath_, "RECREATE");
	if (!TMVAoutFile_->IsOpen() ){
		std::cout << "ERROR : cannot open file "<< outFilePath_ << std::endl;
		exit(-1);}
	std::cout << "Writing output TTree in file "<< outFilePath_ << "\n" << std::endl;

	TMVAoutTreeSGN_ = new TTree("inputSIGNAL", "inputSIGNAL");
	TMVAoutTreeBKG_ = new TTree("inputBACKGROUND", "inputBACKGROUND");
	std::cout << " .... TMVA input trees setting up .... " << std::endl;
	
	TMVAoutTreeSGN_->Branch("pTM_B0", &pTM_B0_S, "pTM_B0/F");
	TMVAoutTreeSGN_->Branch("SVprob", &SVprob_S, "SVprob/F");
	TMVAoutTreeSGN_->Branch("SVchi2", &SVchi2_S, "SVchi2/F");
	TMVAoutTreeSGN_->Branch("LxySign_B0", &LxySign_B0_S, "LxySign_B0/F");
	TMVAoutTreeSGN_->Branch("CosAlpha_B0", &CosAlpha_B0_S, "CosAlpha_B0/F");
	TMVAoutTreeSGN_->Branch("DR_Pi1B0", &DR_Pi1B0_S, "DR_Pi1B0/F");
	TMVAoutTreeSGN_->Branch("pT_Pi1", &pT_Pi1_S, "pT_Pi1/F");
	TMVAoutTreeSGN_->Branch("pT_Rho", &pT_Rho_S, "pT_Rho/F");
	TMVAoutTreeSGN_->Branch("D0_Rho", &D0_Rho_S, "D0_Rho/F");
	TMVAoutTreeSGN_->Branch("M_Rho", &M_Rho_S, "M_Rho/F");
	TMVAoutTreeSGN_->Branch("M_B0", &M_B0_S, "M_B0/F");
	TMVAoutTreeSGN_->Branch("M_mumu", &M_mumu_S, "M_mumu/F");
	TMVAoutTreeSGN_->Branch("M_X3872", &M_X3872_S, "M_X3872/F");
	TMVAoutTreeSGN_->Branch("M_K0s", &M_K0s_S, "M_K0s/F");

	TMVAoutTreeBKG_->Branch("pTM_B0", &pTM_B0_B, "pTM_B0/F");
	TMVAoutTreeBKG_->Branch("SVprob", &SVprob_B, "SVprob/F");
	TMVAoutTreeBKG_->Branch("SVchi2", &SVchi2_B, "SVchi2/F");
	TMVAoutTreeBKG_->Branch("LxySign_B0", &LxySign_B0_B, "LxySign_B0/F");
	TMVAoutTreeBKG_->Branch("CosAlpha_B0", &CosAlpha_B0_B, "CosAlpha_B0/F");
	TMVAoutTreeBKG_->Branch("DR_Pi1B0", &DR_Pi1B0_B, "DR_Pi1B0/F");
	TMVAoutTreeBKG_->Branch("pT_Pi1", &pT_Pi1_B, "pT_Pi1/F");
	TMVAoutTreeBKG_->Branch("pT_Rho", &pT_Rho_B, "pT_Rho/F");
	TMVAoutTreeBKG_->Branch("D0_Rho", &D0_Rho_B, "D0_Rho/F");
	TMVAoutTreeBKG_->Branch("M_Rho", &M_Rho_B, "M_Rho/F");
	TMVAoutTreeBKG_->Branch("M_B0", &M_B0_B, "M_B0/F");
	TMVAoutTreeBKG_->Branch("M_mumu", &M_mumu_B, "M_mumu/F");
	TMVAoutTreeBKG_->Branch("M_X3872", &M_X3872_B, "M_X3872/F");
	TMVAoutTreeBKG_->Branch("M_K0s", &M_K0s_B, "M_K0s/F");


}//BookTreeTMVA()


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
	P4_Reco_K0s.SetPt(B0_finalFit_k0s_pt[Bidx]); P4_Reco_K0s.SetEta(B0_finalFit_k0s_eta[Bidx]); P4_Reco_K0s.SetPhi(B0_finalFit_k0s_phi[Bidx]); 
	P4_Reco_K0s.SetM(B0_K0s_mcFitted_mass[Bidx]);
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



