#ifndef _ssb_analysis_

#define _ssb_analysis_
  
#include <set>
#include <string>
#include <fstream>
#include <cassert>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <TH1.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TGraphErrors.h>
#include "TLorentzVector.h"
#include "TEnv.h"
  
#include "./../analysis/SSBTree.h"

// Textreader
#include "./../TextReader/TextReader.hpp"

#include "./../interface/ssb_eff.hpp"
#include "./../interface/ssb_cpviol.hpp"
#include "./../kinsol/TtFullLepKinSolver.hpp"
#include "./../interface/BTagCalibrationStandalone.h"
#include "./../interface/LumiReWeighting.h"

#include "./../roccor.2016.v3/RoccoR.h"

#include "./../KinSolv/analysisUtils.h"
#include "./../KinSolv/KinematicReconstruction.h"
#include "./../KinSolv/KinematicReconstructionSolution.h"

using namespace std;

class ssb_analysis : public SSBTree 
{
   public:
      //declare functions
      ssb_analysis(TTree *tree=0);
      virtual ~ssb_analysis();

      //basic frame
      virtual void Loop( char *logfile );
      void Start( int genLoopon );
      void End();
      
      //user define functions
      void SetOutputFileName(char *outname);
      void DeclareHistos();

   private:
      //put variables that you want
      char *outfile;
      TFile *fout;
      // vector for ChargeMisId

   public:

      //declare histograms

      // Cut Flow
      TH1D *h_cf_elept; 
      TH1D *h_cf_metpt; 
};
#endif

#ifdef ssb_analysis_cxx

ssb_analysis::ssb_analysis(TTree *tree)
{
   if (tree == 0)
   {
      printf("ERROR: Can't find any input tree.\n");
   }
   Init(tree);
   
}

ssb_analysis::~ssb_analysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   delete fout;
}

#endif
