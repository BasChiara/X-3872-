#include "../include/Application.hh"

#if Application == 1
#include "../include/HLTapply.h"
#endif

//C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

//ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

using namespace std;

int main(int argc, char* argv[]) {

	// inputs from shell
	char inputFileName[500];
	//char outputFileName[500];
	char dataset[500];
	if ( argc < 2 ){
		std::cout << " missing argument: insert the file and the dataset you want to use :-)" << std::endl; 
		std::cout << " ./X3872App inputFile [dataset]" << std::endl;
		std::cout << " use ONLY the following options" << std::endl;
		std::cout << "                  inputFile [MC]" << std::endl;
		std::cout << "                  inputFile [data17B]" << std::endl;
		std::cout << "                  inputFile [data17C]" << std::endl;
		std::cout << "                  inputFile [data17D]" << std::endl;
		std::cout << "                  inputFile [data17E]" << std::endl;
		return 1;
	}
	
	strcpy(inputFileName,argv[1]);
	strcpy(dataset, argv[2]);
	
	// -------------------------
	// Loading the file from a .txt
	TChain *theChain = new TChain("Events");
	char Buffer[5000];
	char MyRootFile[10000];
	TString ChainPath("");
	std::cout << "input: " << inputFileName << std::endl;
	ifstream *inputFile = new ifstream(inputFileName);

	while( !(inputFile->eof()) ){
		inputFile->getline(Buffer,500);
		if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer)))
		{
			sscanf(Buffer,"%s",MyRootFile);
			ChainPath = TString(MyRootFile); ChainPath.Append("/Events");
			theChain->Add(TString(MyRootFile));
			std::cout << "chaining " << MyRootFile << std::endl;
		}
	}
	inputFile->close();
	delete inputFile;


#if Application == 1

	cout<<" Number of events: " << theChain->GetEntries()<<std::endl;

	HLTapply HLTapp( theChain , dataset);

	HLTapp.Loop();

#endif



	return 0;
}
