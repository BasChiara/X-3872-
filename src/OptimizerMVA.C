#include "../include/OptimizerMVA.h"

// RooFit
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooGenericPdf.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooFitResult.h"
#include "RooPlot.h"
using namespace RooFit;

OptimizerMVA::OptimizerMVA(){

	// SGN scale factor
	SGNfactor = 0.173451;

	// I/O files
	outFile_path    = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/TMVAcutApp.root"; 

	inFileSGN_path  = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/TMVAinputs.root";
	SGNtreeName     = "inputSIGNAL";
	outFileSGN_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/BDTonSGN.root";

	outFileBKG_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/BDTonBKG.root";
	inFileBKG_path  = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/data/merged/BKG_MergData17.root";
	BKGtreeName     = "B0sidebands";

	inFileFitPar_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGNfit/B0params.txt";

	LoadFitRanges();

	SetUpInTreeSGN();
	SetUpInTreeBKG();

	// BDT model
	BDTweightPath_ = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/datasetBDT/weights/TMVAClassification_BDT_nT200_S25_D4_b03_nC35.weights.xml";


}//OptimizerMVA()

OptimizerMVA::~OptimizerMVA(){ 

	InFileSGN->Close();
	InFileBKG->Close();
}

//    ==== UTILITY ====     //
TFile* OptimizerMVA::GetFile(TString file){

   TFile *File(0);
	TString Fpath("");
	if (file == "InS") Fpath = inFileSGN_path; 
	if (file == "InB") Fpath = inFileBKG_path;
	if (file == "xS")  Fpath = outFileSGN_path;
	if (file == "xB")  Fpath = outFileBKG_path;

	// Open SIGNAL file
   if (!gSystem->AccessPathName(Fpath)) {
      File = TFile::Open(Fpath); // check if file in local directory exists
   }
   if (!File) {
      std::cout << "ERROR: could not open signal-data file" << std::endl;
      exit(1);
   }

	return File; 

}//GetFile()


void OptimizerMVA::SetUpInTreeSGN(){

   InFileSGN = GetFile("InS");
   // ===== SIGNAL TREE =====
	TreeSGN = (TTree*)InFileSGN->Get(SGNtreeName);

   TreeSGN->SetBranchAddress( "pTM_B0", &pTM_B0);
   TreeSGN->SetBranchAddress( "SVprob", &SVprob);
   TreeSGN->SetBranchAddress( "LxySign_B0", &LxySign_B0);
   TreeSGN->SetBranchAddress( "CosAlpha_B0", &CosAlpha_B0);
   TreeSGN->SetBranchAddress( "DR_Pi1B0", &DR_Pi1B0);
   TreeSGN->SetBranchAddress( "pT_Pi1", &pT_Pi1);
   TreeSGN->SetBranchAddress( "pT_Rho", &pT_Rho);
	TreeSGN->SetBranchAddress( "D0_Rho", &D0_Rho);
   TreeSGN->SetBranchAddress( "M_Rho", &M_Rho);
   TreeSGN->SetBranchAddress( "M_B0", &M_B0);
   TreeSGN->SetBranchAddress( "M_mumu", &M_MuMu);
   TreeSGN->SetBranchAddress( "M_X3872", &M_X3872);
   TreeSGN->SetBranchAddress( "M_K0s", &M_K0s);

	 // friend-tree with BDT response
   TFile* outputBDT = GetFile("xS");
   TTree* Tree_xS= (TTree*)outputBDT->Get("TreeBDTx_S");
	Tree_xS->SetBranchAddress( "BDTx", &X);
	TreeSGN->AddFriend("TreeBDTx_S");


}//SetUpInTreeSGN()

void OptimizerMVA::SetUpInTreeBKG(){

	// ===== BACKGROUND TREE =====

   InFileBKG = GetFile("InB");
   TreeBKG = (TTree*)InFileBKG->Get(BKGtreeName);

   TreeBKG->SetBranchAddress( "pTM_B0", &pTM_B0);
   TreeBKG->SetBranchAddress( "SVprob", &SVprob);
   TreeBKG->SetBranchAddress( "LxySign_B0", &LxySign_B0);
   TreeBKG->SetBranchAddress( "CosAlpha_B0", &CosAlpha_B0);
   TreeBKG->SetBranchAddress( "DR_Pi1B0", &DR_Pi1B0);
   TreeBKG->SetBranchAddress( "pT_Pi1", &pT_Pi1);
   TreeBKG->SetBranchAddress( "pT_Rho", &pT_Rho);
   TreeBKG->SetBranchAddress( "D0_Rho", &D0_Rho);
   TreeBKG->SetBranchAddress( "M_Rho", &M_Rho);
   TreeBKG->SetBranchAddress( "M_B0", &M_B0);
   TreeBKG->SetBranchAddress( "M_mumu", &M_MuMu);
   TreeBKG->SetBranchAddress( "M_X3872", &M_X3872);
   TreeBKG->SetBranchAddress( "M_K0s", &M_K0s);
	
   TFile* outputBDT = GetFile("xB");
   TTree* Tree_xB= (TTree*)outputBDT->Get("TreeBDTx_B");
	Tree_xB->SetBranchAddress( "BDTx", &X);
	TreeBKG->AddFriend("TreeBDTx_B");



}// SetUpInTreeBKG()



void OptimizerMVA::LoadFitRanges(){
	
	std::string line;
	int Nline = 0;

	char ParName[30];
	double err;
	// --- B0 FIT 
	std::ifstream inFileParB0(inFileFitPar_path);	
	if(!inFileParB0.is_open()) std::cout << "ERROR cannot open " << inFileFitPar_path << std::endl;
	while(!inFileParB0.eof()){

		getline(inFileParB0, line); Nline++;
		//std::cout << Nline << "\t" << line << std::endl;
		if(line.c_str()[0] == '#') continue;
		if(Nline == 3) sscanf(line.c_str(), "%s %lf", ParName, &BKGregNearL);
		if(Nline == 4) sscanf(line.c_str(), "%s %lf", ParName, &BKGregFarL);
		if(Nline == 5) sscanf(line.c_str(), "%s %lf", ParName, &BKGregNearR);
		if(Nline == 6) sscanf(line.c_str(), "%s %lf", ParName, &BKGregFarR);
		if(Nline == 9) sscanf(line.c_str(), "%s %lf", ParName, &SGNregL);
		if(Nline ==10) sscanf(line.c_str(), "%s %lf", ParName, &SGNregR);
	}
	inFileParB0.close();


}//LoadFitRanges()

// ==== ANALYSIS TOOLS ==== //

int OptimizerMVA::ApplyBDT(){

// --> TMVA SET UP
   // This loads the library
   TMVA::Tools::Instance();
   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;
 
 
   // Create the Reader object
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
 
   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   Float_t pTM_B0, SVprob, LxySign_B0, CosAlpha_B0;  
   Float_t DR_Pi1B0, pT_Pi1, pT_Rho, D0_Rho; 
   reader->AddVariable( "pTM_B0", &pTM_B0);
   reader->AddVariable( "SVprob", &SVprob);
   reader->AddVariable( "LxySign_B0", &LxySign_B0);
   reader->AddVariable( "CosAlpha_B0", &CosAlpha_B0);
   reader->AddVariable( "DR_Pi1B0", &DR_Pi1B0);
   reader->AddVariable( "pT_Pi1", &pT_Pi1);
   reader->AddVariable( "pT_Rho", &pT_Rho);
   reader->AddVariable( "D0_Rho", &D0_Rho);
 
   // Spectator variables declared in the training have to be added to the reader, too
   Float_t M_Rho; 
   reader->AddSpectator( "M_Rho", &M_Rho);
 
   // Book the BDT method
   TString methodName = TString("BDT method"); 
	reader->BookMVA( methodName, BDTweightPath_);

// --> SIGNAL & BACKGROUND files reading
	TFile* input_signal = GetFile("InS");
	std::cout << "     TMVAClassificationApp    : Using input file (SIGNAL): " << input_signal->GetName() << std::endl;
 
   TFile *input_background = GetFile("InB");
	std::cout << "     TMVAClassificationApp    : Using input file (BACKGROUND): " << input_background->GetName() << std::endl;

// --> EVENT LOOP 
 

	float M_B0, M_MuMu, M_X3872, M_K0s;
   // ===== SIGNAL TREE =====
   std::cout << "--- Select signal sample" << std::endl;
   TTree* TreeSGN = (TTree*)input_signal->Get(SGNtreeName);

   TreeSGN->SetBranchAddress( "pTM_B0", &pTM_B0);
   TreeSGN->SetBranchAddress( "SVprob", &SVprob);
   TreeSGN->SetBranchAddress( "LxySign_B0", &LxySign_B0);
   TreeSGN->SetBranchAddress( "CosAlpha_B0", &CosAlpha_B0);
   TreeSGN->SetBranchAddress( "DR_Pi1B0", &DR_Pi1B0);
   TreeSGN->SetBranchAddress( "pT_Pi1", &pT_Pi1);
   TreeSGN->SetBranchAddress( "pT_Rho", &pT_Rho);
	TreeSGN->SetBranchAddress( "D0_Rho", &D0_Rho);
   TreeSGN->SetBranchAddress( "M_Rho", &M_Rho);
   TreeSGN->SetBranchAddress( "M_B0", &M_B0);
   TreeSGN->SetBranchAddress( "M_mumu", &M_MuMu);
   TreeSGN->SetBranchAddress( "M_X3872", &M_X3872);
   TreeSGN->SetBranchAddress( "M_K0s", &M_K0s);
	 // friend-tree with BDT response
	Float_t Xs;
	TString outFileSGN_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/BDTonSGN.root";
	TFile* output_signal = new TFile(outFileSGN_path, "RECREATE"); 
	TTree* TreeBDTx_S = new TTree("TreeBDTx_S", "friend with BDT response"); 
	TreeBDTx_S->Branch("BDTx", &Xs, "BDTx/F");


	// ===== BACKGROUND TREE =====

   std::cout << "--- Select background sample" << std::endl;
	TreeBKG = (TTree*)input_background->Get(BKGtreeName);

   TreeBKG->SetBranchAddress( "pTM_B0", &pTM_B0);
   TreeBKG->SetBranchAddress( "SVprob", &SVprob);
   TreeBKG->SetBranchAddress( "LxySign_B0", &LxySign_B0);
   TreeBKG->SetBranchAddress( "CosAlpha_B0", &CosAlpha_B0);
   TreeBKG->SetBranchAddress( "DR_Pi1B0", &DR_Pi1B0);
   TreeBKG->SetBranchAddress( "pT_Pi1", &pT_Pi1);
   TreeBKG->SetBranchAddress( "pT_Rho", &pT_Rho);
   TreeBKG->SetBranchAddress( "D0_Rho", &D0_Rho);
   TreeBKG->SetBranchAddress( "M_Rho", &M_Rho);
   TreeBKG->SetBranchAddress( "M_B0", &M_B0);
   TreeBKG->SetBranchAddress( "M_mumu", &M_MuMu);
   TreeBKG->SetBranchAddress( "M_X3872", &M_X3872);
   TreeBKG->SetBranchAddress( "M_K0s", &M_K0s);
	 // friend-tree with BDT response
	Float_t Xb;
	TFile* output_background = new TFile(outFileBKG_path, "RECREATE"); 
	TTree* TreeBDTx_B = new TTree("TreeBDTx_B", "friend with BDT response"); 
	TreeBDTx_B->Branch("BDTx", &Xb, "BDTx/F");
	


	// ---> SIGNAL evaluation 
   std::cout << "--- Processing: " << TreeSGN->GetEntries() << " signal events" << std::endl;
   TStopwatch sw;
   sw.Start();

	for (Long64_t ievt=0; ievt<TreeSGN->GetEntries();ievt++) {
 
		if (ievt%100 == 0) std::cout << "    ... Processing signal-event: " << ievt << std::endl;
		TreeSGN->GetEntry(ievt);
		Xs = reader->EvaluateMVA(methodName);

		TreeBDTx_S->Fill(); 
   }

	// ---> BACKGROUND evaluation 
   std::cout << "--- Processing: " << TreeBKG->GetEntries() << " background events" << std::endl;

	for (Long64_t ievt=0; ievt<TreeBKG->GetEntries();ievt++) {
 
		if (ievt%1000 == 0) std::cout << "    ... Processing background-event: " << ievt << std::endl;
		TreeBKG->GetEntry(ievt);
		Xb = reader->EvaluateMVA(methodName);

		TreeBDTx_B->Fill();
   }


   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();
	
	output_signal->cd();
	TreeBDTx_S->Write();

	output_background->cd();
	TreeBDTx_B->Write();

	// close files
	output_signal->Close();
	output_background->Close();
	input_signal->Close();
	input_background->Close();	



	return 0;

} //ApplyBDT()



double OptimizerMVA::BKG_NevExtraction(const double Xcut, const double MRcut){

	double Bs;

	// ==== ROOFIT SET UP ==== //
	
	RooRealVar B0m("B0m", "M(B_0)\\ [GeV]", 5.07947, 5.47934);

	B0m.setRange(   "SB1"    , BKGregFarL , BKGregNearL);       // SIDEBAND 1
	B0m.setRange(   "SB2"    , BKGregNearR, BKGregFarR);        // SIDEBAND 2
	B0m.setRange("FULLregion", BKGregFarL , BKGregFarR);        // FULLregion RANGE 
	B0m.setRange("SGNregion" , SGNregL    , SGNregR);           // SGNregion RANGE 

	RooDataSet B0m_ds("B0m_ds","B0m_ds", RooArgSet(B0m));

	// Expo + C	
	RooRealVar tau("tau", "", -15.0, -30.0 , -5.);
	RooExponential expo("expo", "", B0m, tau);
	RooProdPdf BKGdesc("BKGdesc", "BKGdesc", expo); 

	RooPolynomial Const("Const", "", B0m);
	RooProdPdf BKGbase("BKGbase", "BKGbase", Const); 

	
	// Fermi 
	RooRealVar Slope("Slope", "", 10.0, 1.0 , 30.);
	RooRealVar Flex("Flex", "", 4., 1. , 5.2);
	RooRealVar C("C", "", 0.5, 0. , 1.);
	
	RooGenericPdf Fermi("Fermi", "", "1./(1. + exp(( @0 - @1)*@2)) + @3", RooArgList(B0m, Flex, Slope, C));

	// Build PDF
	RooRealVar bkg_yield("nBKG", "", 1600, 100, 10000);
	RooRealVar f("f", "", 0., 1.);
	RooAddPdf BKGmodel("BKGmodel", "B0 sidebans shape", Fermi, bkg_yield);
	//RooAddPdf BKGmodel("BKGmodel", "B0 sidebans shape", RooArgList(expo, Const), f);


	// ==== LOOP ON EVENTS ==== // 
	int N = TreeBKG->GetEntries();
	for (int i = 0; i < N; i++){

		TreeBKG->GetEntry(i);
		if( ( X < Xcut) || ( M_Rho < MRcut) ) continue;
		B0m = M_B0;
		B0m_ds.add(RooArgSet(B0m));

	}
	B0m_ds.Print("v");
	// ==== FIT PERFORM ==== //
	
	RooFitResult *ResBKGmodelSB2 = BKGmodel.fitTo(B0m_ds, Range("SB2"), Save());
	//ResBKGmodelSB2->Print();
	RooFitResult *ResBKGmodelSB1 = BKGmodel.fitTo(B0m_ds, Range("SB1"), Save());
	//ResBKGmodelSB1->Print();
	RooFitResult *ResBKGmodel= BKGmodel.fitTo(B0m_ds, Range("SB1,SB2"), Save());
	//ResBKGmodel->Print("v");

	// ==== SAVE RESULTS ==== //
	RooPlot* B0m_fr = B0m.frame();
	B0m_ds.plotOn(B0m_fr);
	BKGmodel.plotOn(B0m_fr, Range("FULLregion"), RooFit::NormRange("SB1,SB2"));
	BKGmodel.paramOn(B0m_fr, Layout(0.60));
	B0m_fr->SetTitle(Form("X-cut = %.2f Mrho-cut = %.3f", Xcut, MRcut));

	
	TText *TxtChi2= new TText(5.2, B0m_fr->GetMaximum()*0.85, Form("Chi^2 = %.3f", B0m_fr->chiSquare()));
   TxtChi2->SetTextSize(0.035);
   TxtChi2->SetTextColor(kRed);
   B0m_fr->addObject(TxtChi2);
	//std::cout << " ---> Chi^2 = " << B0m_fr->chiSquare() << std::endl;

	TH2 *h_Corr = ResBKGmodel->correlationHist();
	
	// ==== Bs ==== //
	RooAbsReal* IntSreg = BKGmodel.createIntegral(B0m, NormSet(B0m), Range("SGNregion")); //integrate sgn region
	double Is = IntSreg->getVal();
	RooAbsReal* IntBreg = BKGmodel.createIntegral(B0m, NormSet(B0m), Range("SB1,SB2"));   //integrate bkg region
	double Ib = IntBreg->getVal();

	Bs = B0m_ds.sumEntries() * Is/Ib; 

	// ==== SAVE ON FILE ==== //
	TString outFileName("/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/CutOptimization/FitRes/FitBKG_X");
	if (Xcut < 0.) outFileName.Append(Form("m"));
	outFileName.Append(Form("%.0f_MR%.0f_Plus.root", fabs(Xcut)*100., MRcut*1000.));

	TFile* outFitFile = new TFile(outFileName, "RECREATE");
	B0m_fr->Write("FitPlot");
	ResBKGmodel->Write("FitResults");
	h_Corr->Write();

	outFitFile->Close();
	

	return Bs;

}//BKG_NevExtraction()




double OptimizerMVA::SGN_NevExtraction(const double Xcut, const double MRcut){

	int iSs = 0;
	double Ss;
	bool inSGNregion = false;
	int Nentries = TreeSGN->GetEntries();

	for(int ien = 0; ien < Nentries; ien++){
		TreeSGN->GetEntry(ien);
		inSGNregion = (M_B0 > SGNregL) && (M_B0 < SGNregR);
		if( (X < Xcut) || (M_Rho < MRcut)  || (!inSGNregion) ) continue;
		iSs ++;
	}
	Ss = iSs * SGNfactor;
	return Ss;

}//SGN_NevExtraction()
	



double OptimizerMVA::PunziSign(const double S, const double B){

	double PunziSign, Sign_denom;
	double SGNeff;
	const double a = 5.0;    // #sigmas corresp. to 5sigma significance level
	const double b = 2.0; // #sigmas corresp. to CL (90%--> 1.2816) (95% --> 1.6448) (CMStwiki --> 2.)

	SGNeff = S / (TreeSGN->GetEntries() * SGNfactor);
	
	Sign_denom = b*b + 2.*a*sqrt(B) + b*sqrt(b*b + 4.*a*sqrt(B) + 4.*B );
	PunziSign = SGNeff/Sign_denom;

	return PunziSign;

}
