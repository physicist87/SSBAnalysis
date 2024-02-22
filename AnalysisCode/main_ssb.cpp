//////////////////////////////////////////////
//                                          //
//                                          //
//  Author: Sehwook Lee, merciful@fnal.gov  //
//                                          //
//                                          //
//////////////////////////////////////////////

#include <iostream>
#include <string>
#include <vector>
#include <ctime>

#include <TROOT.h>
#include <TUnixSystem.h>
#include <TChain.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>

#include "analysis/SSBTree.h"
#include "./interface/ssb_analysis.hpp"

using namespace std;

TROOT root ("Plots", "Program for CMS Analysis");

//argc: # of arguments, argv:array for arguments
int main(int argc, char **argv)
{
   printf("The number of options is: %i\n",argc-1);

   if(argc<2)
   {
      printf("At least, you have to set 1, 2\n");
      printf("1. Input filelist\n");
      printf("2. Output file\n");
      printf("3. GenLoop On or Off\n");
/*      printf("4. Muon Candidate Eta Cut On or Off\n");
      printf("5. Muon Candidate Eta Cut \n");
      printf("6. Charge Opposite Sign On or Off\n");
      printf("7. Dimuon Mass Cut On or Off\n");
      printf("8. Dimuon Mass Cut\n");*/
      return 1;
   }

   for (int iopt=0; iopt<argc; iopt++)
   {
      printf("Option %i = %s\n",iopt,argv[iopt]);
   }

   ///read input options

   //cout << "argc: " << argc << endl;

   char *flist = argv[1];
   printf("Input filelist = %s\n",flist);

   char *outname = argv[2];
   printf("Output file name = %s\n",outname);

   int genLoop_on = atoi(argv[3]);
   printf("Turn On or Off GenLoop : %d\n", genLoop_on);
   cout << "dkdkdk" << genLoop_on << endl;
   
//   char *logfile = argv[4];
//   printf("Output Log File Name = %s.txt\n",logfile);
/*   Bool_t trigger_pass_on =false;
   dimu_mass_cut_on = atoi(argv[3]);
   printf("Turn On or Off Trigger Pass: %d\n", trigger_pass_on);*/

/*   Double_t dimu_mass_cut ((Double_t) atof(argv[4]));
   printf("Dimuon mass cut = %f GeV\n",dimu_mass_cut);*/


   //merge files
   FILE *filelist;
   char filename[1000];
   string filelistDir, filelistName, filelistPath;

   filelistDir = "./input/";
 
   cout << endl;
   filelistName = argv[1];
   filelistPath = filelistDir + filelistName;
   filelist = fopen(filelistPath.c_str(),"r");
  
   std::vector<double> genentries_pertree;

   genentries_pertree.clear();

   std::vector<double> entries_pertree;

   entries_pertree.clear();

   while(filelist==NULL)
   {
      cout << "File not found, please try again." << endl;
      cout << "Filelist you want to use: " << filelistDir;
      cin >> filelistName;      
      filelistPath = filelistDir + filelistName;
      filelist=fopen(filelistPath.c_str(),"r");
   }

//   TChain *chgen = new TChain("demo/SSBGenTree");
   TChain *ch    = new TChain("ssbanalyzer/SSBTree"   );
   //TChain *ch    = new TChain("ssbanalyzer/SSBMiniTree"   );

   cout << endl;
   cout << "start merging file(s)" << endl;
   
   while (fscanf(filelist, "%s", filename) != EOF)
   {
      cout << "adding: " << filename << endl;
      ch->Add(filename, 0);
      entries_pertree.push_back(ch->GetEntries());
   }
   cout << "Total number of events after merging root files: " << ch->GetEntries() << endl;



   //gDirectory->Add(ch);
   //gDirectory->pwd();
   //gDirectory->ls("-l");
   //gDirectory->cd("rootree:/sync");
   //gDirectory->pwd();
   //gDirectory->GetList()->FindObject("MuID");
   //gDirectory->Print();
   //cout <<"ssibal " << gDirectory->GetPath() << endl;
   //TTree* tree = (TTree*)gDirectory->Get("sync/MuID");

   ssb_analysis *ssb = new ssb_analysis(ch);
   if ( genLoop_on == 1 ) 
   {
      
      ssb->SetInputFileName(flist);
      ssb->SetOutputFileName(outname);
      ssb->GetNtupleTotalEvent( ch->GetEntries() );
      ssb->Start( genLoop_on );
      ssb->Loop( flist );
      ssb->End();

   }
   else if ( genLoop_on == 0 )
   {
      
      ssb->SetInputFileName(flist);
      ssb->SetOutputFileName(outname);
      ssb->GetNtupleTotalEvent( ch->GetEntries() );
      ssb->Start( genLoop_on );
      ssb->Loop( flist );
      ssb->End();
      
   }
   else if ( genLoop_on == 2 )
   {
  
   }
   else { cout << "genLoop_on has wrong value..." << endl;}
   delete ssb;
   cout << "Delete ssb..." << endl << endl; 
   delete ch;
   cout << "End processing..." << endl << endl; 
       
   return 0;
}

