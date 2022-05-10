#include "./src/B0toX3872K0s.C"
#include "./src/TriggerSelection.C"
#include "./src/SGNvsBKGvariables.C"
#include <iostream>
#include <string>
#include <TStyle.h>
#include <TCanvas.h>
#include <dirent.h>
#include <sys/stat.h>
using namespace std;


TChain* LoadTree(){

	//================ Loading files
	const string dir_name = "/eos/cms/store/group/phys_egamma/crovelli/LowPtEle/B0ToXKs_2022Apr22/BdToX3872Ks_X3872ToJpsiRho_BMuonFilter_DGamma0_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_BdToX3872Ks/220422_134814/0000/";
	const char* directory = "/eos/cms/store/group/phys_egamma/crovelli/LowPtEle/B0ToXKs_2022Apr22/BdToX3872Ks_X3872ToJpsiRho_BMuonFilter_DGamma0_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_BdToX3872Ks/220422_134814/0000";

	if (directory == NULL) {
		std::cout << "ERROR IN ACCESSING THE DIRECTORY" << std::endl;       
		return NULL;	
	}

	DIR* dir_pointer = opendir(directory);//point to the directory
	struct dirent* file = NULL; // entry in the directory
	struct stat file_stats;

	TString tree_path, tree_name = "/Events";
	int Nfiles = 0;

	//================ Creating chain                                                               

	TChain* chain =new TChain("Events");

	while((file = readdir (dir_pointer))){
		if(file == NULL){
			std::cout << "ERROR null pointer to file" << std::endl;
			return NULL;
		}

		if (strcmp(file->d_name, "xNANO_mc_") < 0) continue; // skip "." and ".."
		Nfiles ++;
		//	std::cout << file->d_name << std::endl;
		tree_path = dir_name + file->d_name + tree_name; 
		chain->Add(tree_path);
	}

	std::cout << " ... LOADING " << Nfiles << " FILES ..." << std::endl;
	cout<<" Number of events: " <<chain->GetEntries()<<std::endl;

	return chain;

}//LoadTree()



int RunTriggerSel() {

	TChain* chain = LoadTree();

	//================ Run analysis                                                                                                                         
	TriggerSelection tree( chain );

	tree.Loop();

	return 0;
}// RunTriggerSel()


int RunSgnBkgAnalysis(){ 

	TChain* chain = LoadTree();

	//================ Run analysis                                                                                                                         
	SGNvsBKGvariables tree( chain );

	tree.Loop();

	return 0;
}//RunSgnBkgAnalysis() 
