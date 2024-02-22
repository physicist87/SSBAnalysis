#ifndef _ssb_cpviol_

#define _ssb_cpviol_
  
#include <set>
#include <string>
#include <fstream>
#include <cassert>
#include <map>
#include <vector>
#include <string>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TGraphErrors.h>
#include "TLorentzVector.h"
#include "TEnv.h"
  
// Textreader
#include "./../TextReader/TextReader.hpp"

using namespace std;

class SSBCPViol 
{
   public:

      SSBCPViol();
      ~SSBCPViol();

      //declare functions
      void getEff();
      double getO1Vari( TLorentzVector* top, TLorentzVector* antitop , TLorentzVector* lep_p, TLorentzVector* lep_m );
      double getO2Vari( TLorentzVector* top, TLorentzVector* antitop , TLorentzVector* bjet_m, TLorentzVector* bjet_p );
      double getO3Vari( TLorentzVector* bjet_m, TLorentzVector* bjet_p , TLorentzVector* lep_p, TLorentzVector* lep_m );
      double getO3VariJPRUp( TLorentzVector* bjet_m, TLorentzVector* bjet_p , TLorentzVector* lep_p, TLorentzVector* lep_m );
      double getO3VariJPRDown( TLorentzVector* bjet_m, TLorentzVector* bjet_p , TLorentzVector* lep_p, TLorentzVector* lep_m );
      double getO4Vari( TLorentzVector* bjet_p, TLorentzVector* bjet_m , TLorentzVector* lep_p, TLorentzVector* lep_m );
      double getObVari( TLorentzVector* bjet_m, TLorentzVector* bjet_p , TLorentzVector* lep_p, TLorentzVector* lep_m );
      double getO5Vari( TLorentzVector* bjet_m, TLorentzVector* bjet_p , TLorentzVector* lep_p, TLorentzVector* lep_m );
      double getO6Vari( TLorentzVector* bjet_m, TLorentzVector* bjet_p , TLorentzVector* lep_p, TLorentzVector* lep_m );
      double getO7Vari( TLorentzVector* top, TLorentzVector* antitop , TLorentzVector* lep_p, TLorentzVector* lep_m );
      double getO8Vari( TLorentzVector* top, TLorentzVector* antitop ,TLorentzVector* bjet_m, TLorentzVector* bjet_p , TLorentzVector* lep_p, TLorentzVector* lep_m );
      double getO9Vari( TLorentzVector* bjet_m, TLorentzVector* bjet_p , TLorentzVector* lep_p, TLorentzVector* lep_m );
      double getO10Vari( TLorentzVector* bjet_m, TLorentzVector* bjet_p , TLorentzVector* lep_p, TLorentzVector* lep_m );
      double getO11Vari( TLorentzVector* bjet_m, TLorentzVector* bjet_p , TLorentzVector* lep_p, TLorentzVector* lep_m );
      double getO12Vari( TLorentzVector* bjet_m, TLorentzVector* bjet_p , TLorentzVector* lep_p, TLorentzVector* lep_m );
      double getO13Vari( TLorentzVector* bjet_m, TLorentzVector* bjet_p , TLorentzVector* lep_p, TLorentzVector* lep_m );

      double getDetM( double a1, double a2, double a3, double a4,
                      double b1, double b2, double b3, double b4,
                      double c1, double c2, double c3, double c4,
                      double d1, double d2, double d3, double d4 );

//      double getDot( double a1, double a2, double a3, double a4, 
//                     double b1, double b2, double b3, double b4  );
      double getDot( std::vector<double> v_a1, std::vector<double> v_a2 );
      double getJPR(double j_eta,double j_pt);
      TLorentzVector* JetAngEta( TLorentzVector* jet, TString op );
      TLorentzVector* JetAngPhi( TLorentzVector* jet, TString op );
//      std::vector<double> Subtract4Vec( double a1, double a2, double a3, double a4,
//                                        double b1, double b2, double b3, double b4 );
      std::vector<double> Subtract4Vec( std::vector<double> v_a1 , std::vector<double> v_b1);
      std::vector<double> Add4Vec( std::vector<double> v_a1 , std::vector<double> v_b1);
      void test();
      // Cut step Function 

      // Text Reader for from Jaehoon. 
//      SSBEffReader->ReadFile("./configs/DimuonCorr.txt");


   private:
     
      // variable for 

      std::vector<double> v_AB;
      std::vector<double> v_PLep_m;
      std::vector<double> v_PLep_p;
      std::vector<double> v_Pb_p;
      std::vector<double> v_q_tilde;
      std::vector<double> v_p;
      std::vector<double> v_q;

      double getdot;
      // for Subtract4Vector Function 
      double subtract4Vec;

      // for Add4Vector Function
      double add4Vec;

      // for getObVari Function...
      std::vector<double> v_lep_p;
      std::vector<double> v_lep_m;
      std::vector<double> v_b_p;
      std::vector<double> v_b_m;

   
//      TextReader *SSBEffReader;

//      std::vector<double> dimu_eff;
};
#endif

