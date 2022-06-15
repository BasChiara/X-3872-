#include "../include/PrePostFit.h" 

// CONSTRUCTOR-DESTRUCTOR
PrePostFit::PrePostFit(TTree *tree) : B0toX3872K0s(tree){ 

   GenMumP4.SetM(mMuon); GenMupP4.SetM(mMuon);
   GenPimP4.SetM(mPion); GenPipP4.SetM(mPion);
   GenK0sP4.SetM(mK0s);
   
   P4_Reco_Mu1.SetM(mMuon); P4_Reco_Mu2.SetM(mMuon);
   P4_Reco_Pi1.SetM(mPion); P4_Reco_Pi2.SetM(mPion);
   P4_Reco_K0s.SetM(mK0s);
}

PrePostFit::~PrePostFit(){
   outFile_->Close();
}





void PrePostFit::Loop() {

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   Long64_t Nbreak = nentries + 10, NPrint = (int)nentries/10.; 

   // ----- VARIABLES ----- //
   int Nmatch_Mu = 0, Nmatch_JPsi = 0, Nmatch_Rho =0, Nmatch_X = 0, Nmatch_K0s = 0,  Nmatch_B0 = 0;
   bool ToCountJPsi = true, ToCountRho = true, ToCountX = true, ToCountK0s = true;
   bool isTriggerON, isFiredMu1, isFiredMu2, isFired_RhoPi1, isFired_RhoPi2, isFired_K0sPi1, isFired_K0sPi2;
   int NTriggeredB0 = 0;
   bool isRecoJPsi = false, isRecoRho = false, isRecoX = false, isRecoK0s = false, isRecoB0 = false;

   // ----- HISTOGRAMS ----- //
   int nbins = 100;
   float Mlow, Mhigh;

   // --> MUONS
   //    trigger
   TH1F h_FiredMCmatch_Mu_pT("FiredMCmatch_Mu_pT", "", nbins/2., 2,12);
   TH1F h_FiredMCmatch_Mu_eta("FiredMCmatch_Mu_eta", "", nbins/2, -3.5 , 3.5);
   TH1F h_FiredMCmatch_Mu_dr("FiredMCmatch_Mu_dr", "", nbins, 0, 0.2);
   // --> MUON PAIR (JPsi)
   Mlow = 2.9, Mhigh = 3.3;
   TH1F h_MCmatch_prefit_JPsi_M("MCmatch_prefit_JPsi_M", "", nbins, Mlow, Mhigh);
   TH1F h_MCmatch_VTXfit_JPsi_M("MCmatch_VTXfit_JPsi_M", "", nbins, Mlow, Mhigh);
   TH1F h_MCmatch_TOTfit_JPsi_M("MCmatch_TOTfit_JPsi_M", "", nbins, Mlow, Mhigh);  
   float MuMu_M;
	// --> PIONS
	TH1F h_genPi_pt("genPi_pt", "", nbins, 0., 5.);
	TH1F h_genPi_eta("genPi_eta", "", nbins, -3.5, 3.5);

   	// ... triggered
   	//KINEMATICS
   TH1F h_FiredMCmatch_MuMu_M("FiredMCmatch_MuMu_M", "", nbins, 2.8, 3.4);
   TH1F h_FiredMCmatch_MuMu_pT("FiredMCmatch_MuMu_pT", "", 20, 6, 24);
   TH1F h_FiredMCmatch_MuMu_DCA("FiredMCmatch_MuMu_DCA", "", nbins/2, 0., 0.1);
   	//VTX
   TH1F h_FiredMCmatch_MuMu_LxySign("FiredMCmatch_MuMu_LxySign", "" , 20, 0, 10);
   TH1F h_FiredMCmatch_MuMu_cosAlpha("FiredMCmatch_MuMu_cosAlpha", "", nbins/2, 0.98,1.);
   TH1F h_FiredMCmatch_MuMu_SVp("FiredMCmatch_MuMu_SVp", "", nbins/2, 0., 1.);
   // --> TRACKS
   // tracks passing the trigger
   TH1F h_FiredMCmatch_Trk_pT("FiredMCmatch_Trk_pT", "", nbins/2, 0,5);
   TH1F h_FiredMCmatch_Trk_eta("FiredMCmatch_Trk_eta", "", nbins/2. , -3.5 , 3.5);
   TH1F h_FiredMCmatch_Trk_d0Sign("FiredMCmatch_Trk_d0Sign", "", nbins/2., 0, 8.);
   TH1F h_FiredMCmatch_K0s_vs_RhoTrk_N("FiredMCmatch_K0s_vs_RhoTrk_N", "", 2, 0,2);
   // --> RHO
   Mlow = 0.2, Mhigh = 1.;
   TH1F h_MCmatch_gen_Rho_M("MCmatch_gen_Rho_M", "", nbins, Mlow, Mhigh);
   TH1F h_MCmatch_prefit_Rho_M("MCmatch_prefit_Rho_M", "", nbins, Mlow, Mhigh);
   TH1F h_MCmatch_TOTfit_Rho_M("MCmatch_TOTfit_Rho_M", "", nbins, Mlow, Mhigh);
   float PiPi_M;
   // --> X(3872)
   Mlow = 3.65, Mhigh = 4.05;
   TH1F h_MCmatch_prefit_X3872_M("MCmatch_prefit_X3872_M", "", nbins, Mlow, Mhigh);
   TH1F h_MCmatch_TOTfit_X3872_M("MCmatch_TOTfit_X3872_M", "", nbins, Mlow, Mhigh);
   float MuMuPiPi_M;
   // --> K0 SHORT
   Mlow = 0.45, Mhigh = 0.55; 
   TH1F h_MCmatch_prefit_K0s_M("MCmatch_prefit_K0s_M", "", nbins, Mlow, Mhigh);
   TH1F h_MCmatch_VTXfit_K0s_M("MCmatch_VTXfit_K0s_M", "", nbins, Mlow, Mhigh);
   TH1F h_MCmatch_TOTfit_K0s_M("MCmatch_TOTfit_K0s_M", "", nbins, Mlow, Mhigh);
   // --> B0
   Mlow = 5.1, Mhigh = 5.5;
   TH1F h_MCmatch_prefit_B0_M("MCmatch_prefit_B0_M", "", nbins, Mlow, Mhigh);
   TH1F h_MCmatch_WOMCfit_B0_M("MCmatch_WOMCfit_B0_M", "", nbins, Mlow, Mhigh);
   TH1F h_MCmatch_TOTfit_B0_M("MCmatch_TOTfit_B0_M", "", nbins, Mlow, Mhigh);
   TH1F h_Comb_WOMCfit_B0_M("Comb_WOMCfit_B0_M", "", nbins, Mlow, Mhigh);
   float MuMuPiPiK0s_M;

   // ==> TRIGGER SELECTION
   bool isOkMuMu_L0, isOkMuMu_L1, isOkMuMu_L2, isOkTrk_L0, isOkTrk_L1;
   int Nbins = 50; 
   Mlow = 5., Mhigh = 5.5;

   TH1F h_HLTrigger_B0_M("HLTrigger_B0_M", "",Nbins, Mlow, Mhigh);

   TH1F h_MuTriggerL0_B0_M("MuTriggerL0_B0_M", "",Nbins, Mlow, Mhigh);
   TH1F h_MuTriggerL1_B0_M("MuTriggerL1_B0_M", "",Nbins, Mlow, Mhigh);
   TH1F h_MuTriggerL2_B0_M("MuTriggerL2_B0_M", "",Nbins, Mlow, Mhigh);

   TH1F h_TrkTriggerL0_B0_M("TrkTriggerL0_B0_M", "", Nbins, Mlow, Mhigh);
   TH1F h_TrkTriggerL1_B0_M("TrkTriggerL1_B0_M", "", Nbins, Mlow, Mhigh);
   
   TH1F h_TotTrigger_B0_M("TotTriggerL1_B0_M", "", Nbins, Mlow, Mhigh);


   // Loop on evens
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (jentry +1  == Nbreak) break;    
      if ((jentry+1) % NPrint==0) cout << "--> #EV  " << jentry+1 << endl;

      // --> Find MC-truth
      GenPartFillP4();
		h_MCmatch_gen_Rho_M.Fill(GenRhoP4.M());
		h_genPi_pt.Fill(GenPimP4.Pt()); h_genPi_eta.Fill(GenPimP4.Eta());
		h_genPi_pt.Fill(GenPimP4.Pt()); h_genPi_eta.Fill(GenPipP4.Eta());
      // --> Match MC-truth with reco tracks
      MCmatch_Mum_Idx = -1, MCmatch_Mup_Idx = -1;
      MCmatch_Pim_Idx = -1, MCmatch_Pip_Idx = -1;
      MCmatch_K0s_Idx = -1;
      MCtruthMatching();

      // --> Search for the B0 MC-matching
      ToCountJPsi = true, ToCountRho = true, ToCountX = true, ToCountK0s = true; 
      for (UInt_t b = 0; b < nB0; b++){

	 //JPsi
	 P4_Reco_Mu1.SetPt(B0_MuMu_prefit_mu1_pt[b]); P4_Reco_Mu1.SetEta(B0_MuMu_prefit_mu1_eta[b]); P4_Reco_Mu1.SetPhi(B0_MuMu_prefit_mu1_phi[b]);
	 P4_Reco_Mu2.SetPt(B0_MuMu_prefit_mu2_pt[b]); P4_Reco_Mu2.SetEta(B0_MuMu_prefit_mu2_eta[b]); P4_Reco_Mu2.SetPhi(B0_MuMu_prefit_mu2_phi[b]);
	 MuMu_M = (P4_Reco_Mu1 + P4_Reco_Mu2).M();      
	 //Rho
	 P4_Reco_Pi1.SetPt(B0_PiPi_prefit_pi1_pt[b]); P4_Reco_Pi1.SetEta(B0_PiPi_prefit_pi1_eta[b]); P4_Reco_Pi1.SetPhi(B0_PiPi_prefit_pi1_phi[b]);
	 P4_Reco_Pi2.SetPt(B0_PiPi_prefit_pi2_pt[b]); P4_Reco_Pi2.SetEta(B0_PiPi_prefit_pi2_eta[b]); P4_Reco_Pi2.SetPhi(B0_PiPi_prefit_pi2_phi[b]);
	 PiPi_M = (P4_Reco_Pi1 + P4_Reco_Pi2).M();
	 // X(3872)
	 MuMuPiPi_M = (P4_Reco_Mu1 + P4_Reco_Mu2 + P4_Reco_Pi1 + P4_Reco_Pi2).M();
	 //K0 short
	 P4_Reco_K0s.SetPt(B0_K0s_mcFitted_pt[b]); P4_Reco_K0s.SetEta(B0_K0s_mcFitted_eta[b]); P4_Reco_K0s.SetPhi(B0_K0s_mcFitted_phi[b]);

	 //B0
	 MuMuPiPiK0s_M = (P4_Reco_Mu1 + P4_Reco_Mu2 + P4_Reco_Pi1 + P4_Reco_Pi2 + P4_Reco_K0s).M();


	 isRecoJPsi = ( (MCmatch_Mum_Idx == B0_mu1_idx[b] && MCmatch_Mup_Idx == B0_mu2_idx[b]) || (MCmatch_Mum_Idx == B0_mu2_idx[b] && MCmatch_Mup_Idx == B0_mu1_idx[b]) );
	 if (ToCountJPsi && isRecoJPsi){
	    //	  std::cout << "Reco mu1 " << B0_mu1_idx[b] << "\t mu2 " << B0_mu2_idx[b] << std::endl;
	    Nmatch_JPsi++;  ToCountJPsi = false;

	    h_MCmatch_prefit_JPsi_M.Fill(MuMu_M);
	    h_MCmatch_VTXfit_JPsi_M.Fill(B0_MuMu_fitted_mass[b]);
	    h_MCmatch_TOTfit_JPsi_M.Fill(B0_finalFit_JPsi_mass[b]); 
	 }

	 isRecoRho  = ( (MCmatch_Pim_Idx == B0_pi1_idx[b] && MCmatch_Pip_Idx == B0_pi2_idx[b]) || (MCmatch_Pim_Idx == B0_pi2_idx[b] && MCmatch_Pip_Idx == B0_pi1_idx[b]) );
	 if (ToCountRho && isRecoRho){
	    //	  std::cout << "Reco pi1 " << B0_pi1_idx[b] << "\t pi2 " << B0_pi2_idx[b] << std::endl;
	    Nmatch_Rho++;	 ToCountRho = false;

	    h_MCmatch_prefit_Rho_M.Fill(PiPi_M);
	    h_MCmatch_TOTfit_Rho_M.Fill(B0_finalFit_Rho_mass[b]);
	 }

	 isRecoX = isRecoJPsi && isRecoRho;
	 if (ToCountX && isRecoX ){
	    //	  std::cout << " X ok " << std::endl;
	    Nmatch_X++;  ToCountX = false;

	    h_MCmatch_prefit_X3872_M.Fill(MuMuPiPi_M);
	    h_MCmatch_TOTfit_X3872_M.Fill(B0_finalFit_X_mass[b]);
	 }


	 isRecoK0s  = (  MCmatch_K0s_Idx == B0_k0short_idx[b]);
	 if (isRecoK0s && ToCountK0s){
	    //	  std::cout << "Reco K0s " << MCmatch_K0s_Idx << std::endl;
	    Nmatch_K0s++;  ToCountK0s = false;

	    h_MCmatch_prefit_K0s_M.Fill(B0_K0s_prefit_mass[b]);
	    h_MCmatch_VTXfit_K0s_M.Fill(B0_K0s_nmcFitted_mass[b]);
	    h_MCmatch_TOTfit_K0s_M.Fill(B0_K0s_mcFitted_mass[b]);
	 }

	 isRecoB0 = isRecoX && isRecoK0s;
	 if(isRecoB0) {
	    Nmatch_B0++;

	    h_MCmatch_prefit_B0_M.Fill(MuMuPiPiK0s_M);
	    h_MCmatch_WOMCfit_B0_M.Fill(B0_fitted_mass_womc[b]);
	    h_MCmatch_TOTfit_B0_M.Fill(B0_finalFit_mass[b]);
	 }else h_Comb_WOMCfit_B0_M.Fill(B0_fitted_mass_womc[b]); //isRecoB0
	    
	    // --> CHECK TRIGGER when B0 is fully MC matched
	 if(isRecoB0){
	    // 	REQUIRE 2muons + 1 trk
	    isTriggerON = (bool)HLT_DoubleMu4_JpsiTrk_Displaced;
	    if (isTriggerON) NTriggeredB0++;

	    // .... muons
	    isFiredMu1 = (bool)B0_MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced[b];
	    isFiredMu2 = (bool)B0_MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced[b]; 
	    // mu1 passes the trigger?
	    if (isFiredMu1){
	       h_FiredMCmatch_Mu_pT.Fill(P4_Reco_Mu1.Pt());
	       h_FiredMCmatch_Mu_eta.Fill(P4_Reco_Mu1.Eta());
	       h_FiredMCmatch_Mu_dr.Fill(B0_MuMu_mu1_dr[b]);
	    }
	    // mu2 passes ?
	    if(isFiredMu2){ 
	       h_FiredMCmatch_Mu_pT.Fill(P4_Reco_Mu2.Pt());
	       h_FiredMCmatch_Mu_eta.Fill(P4_Reco_Mu2.Eta());
	       h_FiredMCmatch_Mu_dr.Fill(B0_MuMu_mu2_dr[b]);
	    }
	    // both pass ?
	    if(isFiredMu1 && isFiredMu2){
	       //KINEMATICS(prefit)
	       h_FiredMCmatch_MuMu_M.Fill(MuMu_M);
	       h_FiredMCmatch_MuMu_pT.Fill((P4_Reco_Mu1 + P4_Reco_Mu2).Pt());
	       h_FiredMCmatch_MuMu_DCA.Fill(B0_MuMu_DCA[b]);
	       //VTX
	       h_FiredMCmatch_MuMu_LxySign.Fill(B0_MuMu_LxySign[b]);
	       h_FiredMCmatch_MuMu_cosAlpha.Fill(B0_MuMu_cosAlpha[b]);
	       h_FiredMCmatch_MuMu_SVp.Fill(B0_MuMu_sv_prob[b]);
	    }

	    // .... tracks
	    isFired_RhoPi1 = (bool)B0_PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced[b];
	    isFired_RhoPi2 = (bool)B0_PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced[b];
	    isFired_K0sPi1 = (bool)B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced[b];
	    isFired_K0sPi2 = (bool)B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced[b];
	    // is there at least one trk?
	    if(isFired_RhoPi1 || isFired_RhoPi2 || isFired_K0sPi1 || isFired_K0sPi2){
	       // pion trk
	       if (isFired_RhoPi1 || isFired_RhoPi2){
		  h_FiredMCmatch_K0s_vs_RhoTrk_N.Fill(1.);
		  if(isFired_RhoPi1){
		     h_FiredMCmatch_Trk_pT.Fill(P4_Reco_Pi1.Pt());
		     h_FiredMCmatch_Trk_eta.Fill(P4_Reco_Pi1.Eta());
		     h_FiredMCmatch_Trk_d0Sign.Fill(B0_PiPi_pi1_d0sig[b]);
		  }
		  if(isFired_RhoPi2){
		     h_FiredMCmatch_Trk_pT.Fill(P4_Reco_Pi2.Pt());
		     h_FiredMCmatch_Trk_eta.Fill(P4_Reco_Pi2.Eta());
		     h_FiredMCmatch_Trk_d0Sign.Fill(B0_PiPi_pi2_d0sig[b]);
		  }
	       }
	       // K0s trk
	       if (isFired_K0sPi1 || isFired_K0sPi2){
		  h_FiredMCmatch_K0s_vs_RhoTrk_N.Fill(0.);
		  if(isFired_K0sPi1){
		     h_FiredMCmatch_Trk_pT.Fill(B0_K0s_matchTrack1_pt[b]);
		     h_FiredMCmatch_Trk_eta.Fill(B0_K0s_matchTrack1_eta[b]);
		     h_FiredMCmatch_Trk_d0Sign.Fill(B0_K0s_matchTrack1_D0sign[b]);
		  }
		  if(isFired_K0sPi2){
		     h_FiredMCmatch_Trk_pT.Fill(B0_K0s_matchTrack2_pt[b]);
		     h_FiredMCmatch_Trk_eta.Fill(B0_K0s_matchTrack2_eta[b]);
		     h_FiredMCmatch_Trk_d0Sign.Fill(B0_K0s_matchTrack2_D0sign[b]);
		  }
	       }

	    }// is fired a trk
	 }// isRecoB0

	 if (isRecoB0){
	    ApplyTriggerSelection(b,&isOkMuMu_L0, &isOkMuMu_L1, &isOkMuMu_L2, &isOkTrk_L0, &isOkTrk_L1);
	    if (isOkMuMu_L0) h_MuTriggerL0_B0_M.Fill(MuMuPiPiK0s_M);
	    if (isOkMuMu_L1) h_MuTriggerL1_B0_M.Fill(MuMuPiPiK0s_M);
	    if (isOkMuMu_L2) h_MuTriggerL2_B0_M.Fill(MuMuPiPiK0s_M);
	    if (isOkTrk_L0) h_TrkTriggerL0_B0_M.Fill(MuMuPiPiK0s_M);
	    if (isOkTrk_L1) h_TrkTriggerL1_B0_M.Fill(MuMuPiPiK0s_M);
	    if ( isTriggerON && isOkMuMu_L1 && isOkMuMu_L2 && isOkTrk_L1) h_TotTrigger_B0_M.Fill(MuMuPiPiK0s_M);
	 }

    }// on B0 candidates

  } // on events
  std::cout << "JPsi matching MC thruth:" << Nmatch_JPsi << std::endl;
  std::cout << "Rho matching MC thruth:" << Nmatch_Rho << std::endl;  
  std::cout << "X(3872) matching MC thruth:" << Nmatch_X << std::endl;
  std::cout << "K0s matching MC thruth:" << Nmatch_K0s << std::endl;
  std::cout << "B0 matching MC thruth:" << Nmatch_B0 << "\t triggered " << NTriggeredB0 << std::endl;

  outFile_ = new TFile("./plots/PrePostFit.root", "RECREATE");
  h_MCmatch_prefit_JPsi_M.Write();
  h_MCmatch_VTXfit_JPsi_M.Write();
  h_MCmatch_TOTfit_JPsi_M.Write();  
  h_MCmatch_gen_Rho_M.Write();
  h_genPi_pt.Write();
  h_genPi_eta.Write();
  h_MCmatch_prefit_Rho_M.Write();
  h_MCmatch_TOTfit_Rho_M.Write();
  h_MCmatch_prefit_X3872_M.Write();
  h_MCmatch_TOTfit_X3872_M.Write();
  h_MCmatch_prefit_K0s_M.Write();
  h_MCmatch_VTXfit_K0s_M.Write();
  h_MCmatch_TOTfit_K0s_M.Write();
  h_MCmatch_prefit_B0_M.Write();
  h_MCmatch_WOMCfit_B0_M.Write();
  h_MCmatch_TOTfit_B0_M.Write();
  h_Comb_WOMCfit_B0_M.Write();

  TFile* out_file_1 = new TFile("./plots/HLT_observables.root", "RECREATE");

  h_FiredMCmatch_Mu_pT.Write();
  h_FiredMCmatch_Mu_eta.Write();
  h_FiredMCmatch_Mu_dr.Write();

  h_FiredMCmatch_MuMu_M.Write();
  h_FiredMCmatch_MuMu_pT.Write();
  h_FiredMCmatch_MuMu_DCA.Write();
  h_FiredMCmatch_MuMu_LxySign.Write();
  h_FiredMCmatch_MuMu_cosAlpha.Write();
  h_FiredMCmatch_MuMu_SVp.Write();

  h_FiredMCmatch_K0s_vs_RhoTrk_N.Write();
  h_FiredMCmatch_Trk_pT.Write();
  h_FiredMCmatch_Trk_eta.Write();
  h_FiredMCmatch_Trk_d0Sign.Write();

  out_file_1->Close();

  TFile* out_file_2 = new TFile("./plots/TriggerSelection.root", "RECREATE");

  h_MuTriggerL0_B0_M.Write();
  h_MuTriggerL1_B0_M.Write();
  h_MuTriggerL2_B0_M.Write();

  h_TrkTriggerL0_B0_M.Write();
  h_TrkTriggerL1_B0_M.Write();

  h_TotTrigger_B0_M.Write();

  out_file_2->Close();



}//Loop ()



void PrePostFit::GenPartFillP4(){
   UInt_t MumIdx = -1, MupIdx = -1, PimIdx = -1 , PipIdx = -1, RhoIdx = -1, K0sIdx = -1;
   for (UInt_t g = 0; g < nGenPart; g++){
      if( (GenPart_pdgId[g] ==  isMum) &&  (GenPart_pdgId[GenPart_genPartIdxMother[g]] == isJPsi)) MumIdx = g;
      if( (GenPart_pdgId[g] == -isMum) &&  (GenPart_pdgId[GenPart_genPartIdxMother[g]] == isJPsi)) MupIdx = g;
      if( (GenPart_pdgId[g] == -isPip) &&  (GenPart_pdgId[GenPart_genPartIdxMother[g]] == isRho) 
&& (GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]] == isX3872) ) PimIdx = g;
      if( (GenPart_pdgId[g] ==  isPip) &&  (GenPart_pdgId[GenPart_genPartIdxMother[g]] == isRho)
&& (GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]] == isX3872) ) PipIdx = g;
      if( (GenPart_pdgId[g] ==  isRho) &&  (abs(GenPart_pdgId[GenPart_genPartIdxMother[g]]) == isX3872) ) RhoIdx = g;
      if( (GenPart_pdgId[g] ==  isK0s) &&  (abs(GenPart_pdgId[GenPart_genPartIdxMother[g]]) == abs(isB0))) K0sIdx = g;

   }//on generated ptl
   // mu- 				mu+
   if(MumIdx > 0){   
      GenMumP4.SetPt(GenPart_pt[MumIdx]);   
      GenMumP4.SetEta(GenPart_eta[MumIdx]);
      GenMumP4.SetPhi(GenPart_phi[MumIdx]); 
   } else std::cout << "No generated mu- is found" << std::endl;   

   if(MupIdx > 0){
      GenMupP4.SetPt(GenPart_pt[MupIdx]);
      GenMupP4.SetEta(GenPart_eta[MupIdx]);
      GenMupP4.SetPhi(GenPart_phi[MupIdx]);

   }else std::cout << "No generated mu+ is found" << std::endl;

   // pi- 				pi+   
   if(PimIdx > 0){
      GenPimP4.SetPt(GenPart_pt[PimIdx]);   
      GenPimP4.SetEta(GenPart_eta[PimIdx]); 
      GenPimP4.SetPhi(GenPart_phi[PimIdx]); 
   }else std::cout << "No generated pi- is found" << std::endl;
   if(PimIdx > 0){

      GenPipP4.SetPt(GenPart_pt[PipIdx]);
      GenPipP4.SetEta(GenPart_eta[PipIdx]);
      GenPipP4.SetPhi(GenPart_phi[PipIdx]);
   }else std::cout << "No generated pi+ is found" << std::endl;


   //Rho 
   if(RhoIdx > 0 ){
      GenRhoP4.SetPt(GenPart_pt[RhoIdx]);
      GenRhoP4.SetEta(GenPart_eta[RhoIdx]);
      GenRhoP4.SetPhi(GenPart_phi[RhoIdx]);
      GenRhoP4.SetM(GenPart_mass[RhoIdx]);
   }else std::cout << "No generated K0s is found" << std::endl;
   // K0s
   if(K0sIdx > 0 ){
      GenK0sP4.SetPt(GenPart_pt[K0sIdx]);
      GenK0sP4.SetEta(GenPart_eta[K0sIdx]);
      GenK0sP4.SetPhi(GenPart_phi[K0sIdx]);
   }else std::cout << "No generated K0s is found" << std::endl;
}// gen_Daugh_ptetaphiM()



void PrePostFit::MCtruthMatching(const bool verbose){

   const double DRmin_threshold = 0.03;
   const double DpT_threshold = 0.5;
   float DeltaPt;
 
   // ..... muons ..... //	
   ROOT::Math::PtEtaPhiMVector P4_RecoMu1, P4_RecoMu2;
   float DRminMum = 100., DR_gMum_rMu1, DR_gMum_rMu2, DRminMup = 100., DR_gMup_rMu1, DR_gMup_rMu2;
   int DRminMum_Idx = -1, DRminMup_Idx = -1; 
   // ..... pions ..... //      
   ROOT::Math::PtEtaPhiMVector P4_RecoPi1, P4_RecoPi2;
   float DRminPim = 100., DR_gPim_rPi1, DR_gPim_rPi2, DRminPip = 100., DR_gPip_rPi1, DR_gPip_rPi2;
   int DRminPim_Idx = -1, DRminPip_Idx = -1;


   // ..... K0s ..... //
   ROOT::Math::PtEtaPhiMVector P4_RecoK0s;
   float DRminK0s = 100., DR_gK0s_rK0s;
   int DRminK0s_Idx = -1;

   for (UInt_t b = 0; b < nB0; b++){
      // ..... muons ..... //
      if (Muon_softId[B0_mu1_idx[b]] && Muon_softId[B0_mu2_idx[b]]){ // softId for both muons is required

	 P4_RecoMu1.SetPt(B0_MuMu_prefit_mu1_pt[b]); P4_RecoMu1.SetEta(B0_MuMu_prefit_mu1_eta[b]); P4_RecoMu1.SetPhi(B0_MuMu_prefit_mu1_phi[b]);
	 P4_RecoMu2.SetPt(B0_MuMu_prefit_mu2_pt[b]); P4_RecoMu2.SetEta(B0_MuMu_prefit_mu2_eta[b]); P4_RecoMu2.SetPhi(B0_MuMu_prefit_mu2_phi[b]);

	 DR_gMum_rMu1 = ROOT::Math::VectorUtil::DeltaR(GenMumP4, P4_RecoMu1); // mu1
	 DR_gMup_rMu1 = ROOT::Math::VectorUtil::DeltaR(GenMupP4, P4_RecoMu1);
	 DR_gMum_rMu2 = ROOT::Math::VectorUtil::DeltaR(GenMumP4, P4_RecoMu2); // mu2
	 DR_gMup_rMu2 = ROOT::Math::VectorUtil::DeltaR(GenMupP4, P4_RecoMu2);

	 // DR(mu-, mu1) < DRmin(mu-) + m1 is nearer to mu- than mu2 + DRmin threshold
	 // mu-
	 if( (DR_gMum_rMu1 < DRminMum) && (DR_gMum_rMu1 < DR_gMum_rMu2) && (DR_gMum_rMu1 < DRmin_threshold) ){
	    DRminMum = DR_gMum_rMu1;
	    DRminMum_Idx = B0_mu1_idx[b];
	 }
	 if( (DR_gMum_rMu2 < DRminMum) && (DR_gMum_rMu2 < DR_gMum_rMu1) && (DR_gMum_rMu2 < DRmin_threshold ) ){
	    DRminMum = DR_gMum_rMu2;
	    DRminMum_Idx = B0_mu2_idx[b];
	 }
	 // mu+
	 if( (DR_gMup_rMu1 < DRminMup) && (DR_gMup_rMu1 < DR_gMup_rMu2) && (DR_gMup_rMu1 < DRmin_threshold) ){
	    DRminMup = DR_gMup_rMu1;
	    DRminMup_Idx = B0_mu1_idx[b];
	 }
	 if( (DR_gMup_rMu2 < DRminMup) && (DR_gMup_rMu2 < DR_gMup_rMu1) && (DR_gMup_rMu2 < DRmin_threshold ) ){
	    DRminMup = DR_gMup_rMu2;
	    DRminMup_Idx = B0_mu2_idx[b];
	 }

      }// ... muons //

      // ..... pions ..... //
      if (!ProbeTracks_isMatchedToMuon[B0_pi1_idx[b]] && !ProbeTracks_isMatchedToMuon[B0_pi2_idx[b]]){ // isMatchedToMuon must be false
	 P4_RecoPi1.SetPt(B0_PiPi_prefit_pi1_pt[b]); P4_RecoPi1.SetEta(B0_PiPi_prefit_pi1_eta[b]); P4_RecoPi1.SetPhi(B0_PiPi_prefit_pi1_phi[b]);
	 P4_RecoPi2.SetPt(B0_PiPi_prefit_pi2_pt[b]); P4_RecoPi2.SetEta(B0_PiPi_prefit_pi2_eta[b]); P4_RecoPi2.SetPhi(B0_PiPi_prefit_pi2_phi[b]);

	 DR_gPim_rPi1 = ROOT::Math::VectorUtil::DeltaR(GenPimP4, P4_RecoPi1);//pi1
	 DR_gPip_rPi1 = ROOT::Math::VectorUtil::DeltaR(GenPipP4, P4_RecoPi1);
	 DR_gPim_rPi2 = ROOT::Math::VectorUtil::DeltaR(GenPimP4, P4_RecoPi2);//pi2
	 DR_gPip_rPi2 = ROOT::Math::VectorUtil::DeltaR(GenPipP4, P4_RecoPi2);
	 // pi-
	 if( (DR_gPim_rPi1 < DRminPim) && (DR_gPim_rPi1 < DR_gPim_rPi2) && (DR_gPim_rPi1 < DRmin_threshold) && (DeltaPT(GenPimP4, P4_RecoPi1) < DpT_threshold)){
	    DRminPim = DR_gPim_rPi1;
	    DRminPim_Idx = B0_pi1_idx[b];
	 }
	 if( (DR_gPim_rPi2 < DRminPim) && (DR_gPim_rPi2 < DR_gPim_rPi1) && (DR_gPim_rPi2 < DRmin_threshold) && (DeltaPT(GenPimP4, P4_RecoPi2) < DpT_threshold)){
	    DRminPim = DR_gPim_rPi2;
	    DRminPim_Idx = B0_pi2_idx[b];
	 }
	 // pi+
	 if( (DR_gPip_rPi1 < DRminPip) && (DR_gPip_rPi1 < DR_gPip_rPi2) && (DR_gPip_rPi1 < DRmin_threshold) && (DeltaPT(GenPipP4, P4_RecoPi1) < DpT_threshold)){
	    DRminPip = DR_gPip_rPi1;
	    DRminPip_Idx = B0_pi1_idx[b];
	 }
	 if( (DR_gPip_rPi2 < DRminPip) && (DR_gPip_rPi2 < DR_gPip_rPi1) && (DR_gPip_rPi2 < DRmin_threshold ) && (DeltaPT(GenPipP4, P4_RecoPi2) < DpT_threshold)){
	    DRminPip = DR_gPip_rPi2;
	    DRminPip_Idx = B0_pi2_idx[b];
	 }


      }// ... pions//

      // ..... K0s ..... //
      P4_RecoK0s.SetPt(B0_K0s_mcFitted_pt[b]); P4_RecoK0s.SetEta(B0_K0s_mcFitted_eta[b]); P4_RecoK0s.SetPhi(B0_K0s_mcFitted_phi[b]);
      DR_gK0s_rK0s = ROOT::Math::VectorUtil::DeltaR(GenK0sP4, P4_RecoK0s);
      if ( (DR_gK0s_rK0s < DRminK0s) && (DR_gK0s_rK0s < DRmin_threshold)){
	 DRminK0s = DR_gK0s_rK0s;
	 DRminK0s_Idx = B0_k0short_idx[b];
      }


   }// on B0 candidates


   MCmatch_Mum_Idx = DRminMum_Idx;
   MCmatch_Mup_Idx = DRminMup_Idx;
   

   MCmatch_Pim_Idx = DRminPim_Idx;
   MCmatch_Pip_Idx = DRminPip_Idx;

   MCmatch_K0s_Idx = DRminK0s_Idx;
   if (verbose){
      std::cout << "MC matching indices " << std::endl;
      std::cout << "mu- " << MCmatch_Mum_Idx << "\t mu+ " << MCmatch_Mup_Idx << std::endl;
      std::cout << "pi- " << MCmatch_Pim_Idx << "\t pi+ " << MCmatch_Pip_Idx << std::endl;
      std::cout << "K0s " << MCmatch_K0s_Idx << std::endl;
   }

}//MCtruthMatching()

float PrePostFit::DeltaPT(ROOT::Math::PtEtaPhiMVector genV, ROOT::Math::PtEtaPhiMVector recV){

   return TMath::Abs(genV.Pt() - recV.Pt()) / genV.Pt();

}


void PrePostFit::ApplyTriggerSelection(const int Bidx,bool* isOkMuMu_L0, bool* isOkMuMu_L1, bool* isOkMuMu_L2, bool* isOkTrk_L0, bool*isOkTrk_L1){

   // TRIGGER SELECTION PARAMS
   //.... muons
   const float Min_Mu_pT = 4.,Max_Mu_eta = 2.5, Max_Mu_dr = 2.;
   const float Min_MuMu_pT = 6.9, Low_MuMu_M = 2.9,  High_MuMu_M = 3.3, Max_MuMu_DCA = 0.5;
   const float Min_MuMu_LxyS = 3, Min_MuMu_cosAlpha = 0.9, Min_MuMu_SVp = 0.1;
   //.... tracks   
   const float Min_Trk_pT = 1.2, Max_Trk_eta = 2.5, Min_Trk_D0S = 2.;
   bool isTriggerON, isFiredMu1, isFiredMu2, isFired_RhoPi1, isFired_RhoPi2, isFired_K0sPi1, isFired_K0sPi2;
   bool isSelMu1_L1, isSelMu2_L1,isSelMuMu_L0, isSelMuMu_L1, isSelMuMu_L2, isSelTrk_L0, isSelTrk_L1;

   float MuMu_M = (P4_Reco_Mu1 + P4_Reco_Mu2).M();
   float MuMuPiPiK0s_M = (P4_Reco_Mu1 + P4_Reco_Mu2 + P4_Reco_Pi1+P4_Reco_Pi2 + P4_Reco_K0s).M();
   float MuMu_pT = (P4_Reco_Mu1 + P4_Reco_Mu2).Pt();

   isTriggerON = (bool)HLT_DoubleMu4_JpsiTrk_Displaced;

   isFiredMu1 = (bool)B0_MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced[Bidx];
   isFiredMu2 = (bool)B0_MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced[Bidx]; 

   isFired_RhoPi1 = (bool)B0_PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced[Bidx];
   isFired_RhoPi2 = (bool)B0_PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced[Bidx];
   isFired_K0sPi1 = (bool)B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced[Bidx];
   isFired_K0sPi2 = (bool)B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced[Bidx];


   // Muon TRIGGER
   // LEVEL 0
   isSelMuMu_L0 = (isFiredMu1 && isFiredMu2); 
   *isOkMuMu_L0 = isSelMuMu_L0;
   // LEVEL 1
   isSelMuMu_L1 = true;
   if ( (MuMu_M < Low_MuMu_M) || (MuMu_M > High_MuMu_M) )isSelMuMu_L1 = false; // 2.9 < M(mu mu) < 3.3
   if ( MuMu_pT < Min_MuMu_pT ) isSelMuMu_L1 = false;
   if ( B0_MuMu_DCA[Bidx] > Max_MuMu_DCA) isSelMuMu_L1 = false;

   *isOkMuMu_L1 = (isSelMuMu_L0 && isSelMuMu_L1); //h_MuTriggerL1_B0_M->Fill(MuMuPiPiK0s_M);
   // LEVEL 2
   isSelMuMu_L2 = true;
   if ( B0_MuMu_LxySign[Bidx] < Min_MuMu_LxyS) isSelMuMu_L2 = false;
   if ( B0_MuMu_cosAlpha[Bidx] < Min_MuMu_cosAlpha) isSelMuMu_L2 = false;
   if ( B0_MuMu_sv_prob[Bidx] < Min_MuMu_SVp) isSelMuMu_L2 = false;

   *isOkMuMu_L2 = (isSelMuMu_L0 && isSelMuMu_L1 && isSelMuMu_L2) ; //h_MuTriggerL2_B0_M->Fill(MuMuPiPiK0s_M);

   // Track TRIGGER
   // LEVEL 0
   isSelTrk_L0 = (isFired_RhoPi1 || isFired_RhoPi2 || isFired_K0sPi1 || isFired_K0sPi2) ;
   *isOkTrk_L0 = isSelTrk_L0; 
   // LEVEL 1
   isSelTrk_L1 = false;
   if (isFired_RhoPi1){ 
      if( (P4_Reco_Pi1.Pt() > Min_Trk_pT) && (TMath::Abs(P4_Reco_Pi1.Eta()) < Max_Trk_eta) && (B0_PiPi_pi1_d0sig[Bidx] > Min_Trk_D0S) ) isSelTrk_L1 = true;
   } 
   if (isFired_RhoPi2){ 
      if( (P4_Reco_Pi2.Pt() > Min_Trk_pT) && (TMath::Abs(P4_Reco_Pi2.Eta()) < Max_Trk_eta) && (B0_PiPi_pi2_d0sig[Bidx] > Min_Trk_D0S) ) isSelTrk_L1 = true;
   }
   if (isFired_K0sPi1){
      if( (B0_K0s_matchTrack1_pt[Bidx] > Min_Trk_pT) && (TMath::Abs(B0_K0s_matchTrack1_eta[Bidx]) < Max_Trk_eta) && (B0_K0s_matchTrack1_D0sign[Bidx] > Min_Trk_D0S) ) isSelTrk_L1 = true;
   }
   if (isFired_K0sPi2){
      if( (B0_K0s_matchTrack2_pt[Bidx] > Min_Trk_pT) && (TMath::Abs(B0_K0s_matchTrack2_eta[Bidx]) < Max_Trk_eta) && (B0_K0s_matchTrack2_D0sign[Bidx] > Min_Trk_D0S) ) isSelTrk_L1 = true;
   }

   *isOkTrk_L1 = (isSelTrk_L0 && isSelTrk_L1);// h_TrkTriggerL1_B0_M->Fill(MuMuPiPiK0s_M); 

}//ApplyTriggerSelection()

