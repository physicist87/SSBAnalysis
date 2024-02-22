#ifndef _ssb_eff_

#define _ssb_eff_
  
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
#include <TGraphAsymmErrors.h>
#include "TLorentzVector.h"
#include "TEnv.h"
#include "TFile.h"  
// Textreader
#include "./../TextReader/TextReader.hpp"
// BTagging //
#include "./../interface/BTagCalibrationStandalone.h"
#include "./../interface/LumiReWeighting.h"
using namespace std;

class SSBEffCal 
{
   public:

      SSBEffCal();
      ~SSBEffCal();
      //declare functions
      void getEff();

      /////////////////
      /// Lepton SF ///
      /////////////////
      double DoubleMuon_Eff( TLorentzVector* lep1, TLorentzVector* lep2 );
      //double DoubleMuon_EffROOT( TLorentzVector* lep1, TLorentzVector* lep2, TString muidsys,  TString muisosys );
      double DoubleMuon_EffROOT( TLorentzVector* lep1, TLorentzVector* lep2, TString muidsys,  TString muisosys, TString tracksys );
      double DoubleElec_EffROOT( TLorentzVector* lep1, TLorentzVector* lep2, double ele1sueta, double ele2sueta, TString eleidsys, TString elecrecosys);
      double MuonElec_EffROOT( TLorentzVector* muon, TLorentzVector* elec, double elecsueta, TString lepidsys, TString lepisosys, TString mutrksys, TString elecrecosys );
      double TrackSF(double eta);
      double TrackSFErr(double eta,TString trksys);

      //////////////////
      /// Trigger SF ///
      //////////////////
      double TrigDiMuon_Eff( TLorentzVector* lep1, TLorentzVector* lep2, TString Sys_);
      double TrigDiElec_Eff( TLorentzVector* lep1, TLorentzVector* lep2, TString Sys_);
      double TrigMuElec_Eff( TLorentzVector* lep1, TLorentzVector* lep2, TString Sys_);
      double GetTrgEff(double lep1, double lep2 ,TString Sys_);

      ///////////////////
      /// BTagging SF ///
      ///////////////////
      double Btagging_Eff(TLorentzVector* jet, TString btag, TString btageffsys);
      double Btagging_Eff(double jpt_, double jeta_, TString btag, TString btageffsys, int jet_fl);
      double Btagging_EvenWeight(TLorentzVector* jet, double btagdisc, double btagcut );
      double Btagging_EvenWeight(std::vector<double>v_jetpt,std::vector<double>v_jeteta ,std::vector<double>v_btagdisc, double btagcut, std::vector<int> v_jetf );
      double Btagging_EvenWeightSys(std::vector<double>v_jetpt,std::vector<double>v_jeteta ,std::vector<double>v_btagdisc, double btagcut, std::vector<int> v_jetf ,TString sys_ );

   private:

      TextReader *SSBEffReader;
      TextReader *SSBConfigReader;

      TString RunPeriod;
      double eff;
      double m1pt;
      double m2pt;
      double m1eta;
      double m2eta;

      std::vector<double> dimu_eff;

      TString MuonID_File;
      TString MuonIso_File;
      TString Track_File;
      TString ElecID_File;
      TString ElecReco_File;
      TString Trig_File;
      string BTagCSV_File;

      TFile *f_muid;
      std::vector<TH2F*>v_h_muid;
      TFile *f_muiso;
      std::vector<TH2F*>v_h_muiso;
      TFile *f_track;
      TFile *f_eleid;
      TFile *f_eleiso;
      TFile *f_elereco;
      TFile *f_trg;

      TH2F* H_muid; 
      TH2F* H_muiso; 
      TGraphAsymmErrors* H_track; 
      TH2F* H_eleid;
      TH2F* H_elereco;
      TH2D* H_trig;

      /// Lumi ratio ///
      std::vector<double> v_lumi_ratio;

      /// For BTagging ///
      TFile *f_btag_loose[3];
      TFile *f_btag_medium[3];
      TFile *f_btag_tight[3];

      BTagCalibration *calib; 
      BTagCalibrationReader *reader;
      BTagCalibrationReader *readerup;
      BTagCalibrationReader *readerdown;
      BTagCalibrationReader *readerlf;
      BTagCalibrationReader *readerlfup;
      BTagCalibrationReader *readerlfdown;

      BTagCalibrationReader *reader_;
      //std::vector<BTagCalibration*> v_calib;
      std::vector<BTagCalibrationReader*> v_reader;
      std::vector<BTagCalibrationReader*> v_readerup;
      std::vector<BTagCalibrationReader*> v_readerdown;
      std::vector<BTagCalibrationReader*> v_readerlf;
      std::vector<BTagCalibrationReader*> v_readerlfup;
      std::vector<BTagCalibrationReader*> v_readerlfdown;
      std::vector<BTagCalibrationReader*> v_btagCal;

      TH2D* H_btageff;
      TH2D* H_btageff_fl[3];

      TString jetbtag;
      TString BTagEffSys;
      TString BTagSFSys;

      // Dummy TCanvas //
      TCanvas *mydummycanvas;
};
#endif

