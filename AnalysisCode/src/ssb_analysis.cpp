#define ssb_analysis_cxx

#include <iostream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "TMath.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"

#include "./../interface/ssb_analysis.hpp"
#include "./../CommonTools.hpp"

using namespace std;

void ssb_analysis::SetConfig()
{
 
   // Lep Infor ( Di Muon or Di Electron )
   lepisotype = SSBConfReader->GetText( "Iso_type" );
   lepId      = SSBConfReader->GetText( "Lep_ID"   );
   n_elep_Id  = SSBConfReader->Size( "Lep_ID_cut"  );
  
   // Muon Infor for Only MuEl channel
   Muonisotype = SSBConfReader->GetText( "MuonIso_type" );
   Muonid      = SSBConfReader->GetText( "Muon_ID"      );

   // Electron Infor Only for MuEl channel
   Elecisotype = SSBConfReader->GetText( "ElecIso_type" );
   Elecid      = SSBConfReader->GetText( "Elec_ID"      );
   n_elec_Id   = SSBConfReader->Size(    "Elec_ID_cut"  );
   // Jet Infor
   jetId      = SSBConfReader->GetText("Jet_ID");
   jetbtag    = SSBConfReader->GetText("Jet_btag");


   // Kinematic variables for Object //(not MuEl)
   lep_pt     = SSBConfReader->GetNumber( "Lep_pt"  );
   lep_eta    = SSBConfReader->GetNumber( "Lep_eta" );
   lepisocut  = SSBConfReader->GetNumber( "Iso_cut" );
   jet_pt     = SSBConfReader->GetNumber( "Jet_pt"  );
   jet_eta    = SSBConfReader->GetNumber( "Jet_eta" );
   met_cut    = SSBConfReader->GetNumber( "MET_cut" );

   // Kinematic variables for Object (MuEl)
   muon_pt     = SSBConfReader->GetNumber( "Muon_pt"     );
   muon_eta    = SSBConfReader->GetNumber( "Muon_eta"    );
   muonisocut  = SSBConfReader->GetNumber( "MuonIso_cut" );
   elec_pt     = SSBConfReader->GetNumber( "Elec_pt"     );
   elec_eta    = SSBConfReader->GetNumber( "Elec_eta"    );
   elecisocut  = SSBConfReader->GetNumber( "ElecIso_cut" );
   

   // variables for Jet Cleaning 
   muiso_type = SSBConfReader->GetText( "Muiso_for_Jetcleanning" );
   muId       = SSBConfReader->GetText( "Muid_for_Jetcleanning"  );
   eliso_type = SSBConfReader->GetText( "Eliso_for_Jetcleanning" );
   elId       = SSBConfReader->GetText( "Elid_for_Jetcleanning"  );

   // Kinematic variables for Jet Cleaning 
   mu_pt_jetcl   = SSBConfReader->GetNumber( "Mupt_for_Jetcleanning"     );
   mu_eta_jetcl  = SSBConfReader->GetNumber( "Mueta_for_Jetcleanning"    );
   mu_iso_jetcl  = SSBConfReader->GetNumber( "Muisocut_for_Jetcleanning" );
   el_pt_jetcl   = SSBConfReader->GetNumber( "Elpt_for_Jetcleanning"     );
   el_eta_jetcl  = SSBConfReader->GetNumber( "Eleta_for_Jetcleanning"    );
   el_iso_jetcl  = SSBConfReader->GetNumber( "Elisocut_for_Jetcleanning" );
   n_el_id_jetcl = SSBConfReader->Size(      "Elidcut_for_Jetcleanning"  );

   if (n_el_id_jetcl == 1)
   {
      el_id_jetcl_1 = SSBConfReader->GetNumber("Elidcut_for_Jetcleanning");
      el_id_jetcl_2 = -999;
   }
   else if (n_el_id_jetcl == 2 )
   {
      el_id_jetcl_1 = SSBConfReader->GetNumber("Elidcut_for_Jetcleanning", 1);
      el_id_jetcl_2 = SSBConfReader->GetNumber("Elidcut_for_Jetcleanning", 2);
   }
   else {cout << " Electron Jetcleaning error" <<endl; }


   /////////////////////////////
   // For KinematicSolver !!! //
   /////////////////////////////

//   double nupair[] = { 30.7137, 56.2880, 23.0744, 59.1015,24.9145 }; // for 8TeV
   double nupair[] = { 1/4.53731e+06,  5.81649e+01,  2.31911e+01,  1.19038e+04, 
                       -9.78134e+01 , -2.76894e+01, -2.01114e-02, -5.07004e+00,
                         5.76608e+01,  2.34419e+01,  1.18119e+04,  5.26711e+03,
                         1.64080e+03,  5.86524e+00, -8.07173e+01, -5.44618e+00 };
   int nupair_size = sizeof(nupair)/sizeof(nupair[0]);
   
   for (int i = 0; i < nupair_size; ++i ){ nupars_.push_back(nupair[i]); }

   topmass_begin = 50.0;
   topmass_end   = 500.0;
   topmass_step_ = 1.0;


}
void ssb_analysis::GetVariables()
{

   ///////////////////////
   /// Case of Di-Muon ///
   /////////////////////// 
   if ( TString(Decaymode).Contains( "dimuon" ) ) 
   { 

      //////////////////////
      /// Muon condition ///
      //////////////////////

      /// Muon Iso type
      if ( TString(lepisotype).Contains( "PFIsodbeta03" ) )       { v_lepton_iso = Muon_PFIsodBeta03; }
      else if ( TString(lepisotype).Contains( "PFIsodbeta04" ) )  { v_lepton_iso = Muon_PFIsodBeta04; }
      else { cout << "Muon Condition Iso Error" << endl; } 

      /// Muon ID
      if      (TString(lepId).Contains( "Loose"  ) ) { v_lepton_Id = Muon_isLoose;  }
      else if (TString(lepId).Contains( "Tight"  ) ) { v_lepton_Id = Muon_isTight;  }
      else { cout << "Muon ID Error" << endl; } 

   }
   ///////////////////////////
   /// Case of Di-Electron ///
   ///////////////////////////
   else if (TString(Decaymode).Contains( "dielec" ) )
   {
      /// Electron iso type
      if      ( TString( lepisotype ).Contains( "PFIsoRho03"   ) ) { v_lepton_iso = Elec_PFIsoRho03;   }
      else if ( TString( lepisotype ).Contains( "PFIsoRho04"   ) ) { v_lepton_iso = Elec_PFIsoRho04;   }
      else { cout << "Electron Iso type Error at Case of Di-Electron" << endl; }
   
      // Electron ID

      /// simple cut based electron ID
      if      (TString(lepId).Contains( "Loose"  ) ) { v_lepton_Id = Elec_SCB_Loose;  }
      else if (TString(lepId).Contains( "Medium" ) ) { v_lepton_Id = Elec_SCB_Medium; }
      else if (TString(lepId).Contains( "Tight"  ) ) { v_lepton_Id = Elec_SCB_Tight;  }
      else if (TString(lepId).Contains( "Veto"   ) ) { v_lepton_Id = Elec_SCB_Veto;   }
      else { cout << "Electron ID Error" << endl; }

   } 
   ////////////////////
   /// Case of MuEl ///
   ////////////////////
   else if (TString(Decaymode).Contains( "muel" ) )
   {
      /// Muon Iso type
      if      ( TString(Muonisotype).Contains( "PFIsodbeta03" ) )  { v_mlep_iso = Muon_PFIsodBeta03; }
      else if ( TString(Muonisotype).Contains( "PFIsodbeta04" ) )  { v_mlep_iso = Muon_PFIsodBeta04; }
      else { cout << "Muon Condition Iso Error at Case of MuEl" << endl; } 

      /// Muon ID
      if      ( TString(Muonid).Contains( "Loose"  ) ) { v_muon_Id = Muon_isLoose;  }
      else if ( TString(Muonid).Contains( "Tight"  ) ) { v_muon_Id = Muon_isTight;  }
      else { cout << "Muon Condition ID Error" << endl; } 

      /// Electron iso type
      if      ( TString( Elecisotype ).Contains( "PFIsoRho03" ) )   { v_elec_iso = Elec_PFIsoRho03;   }
      else if ( TString( Elecisotype ).Contains( "PFIsoRho04" ) )   { v_elec_iso = Elec_PFIsoRho04;   }
      else { cout << "Electron Iso type Error at Case of MuEl" << endl; }
   
      // Electron ID SCB
     
      if      (TString(Elecid).Contains( "Loose"  ) ) { v_electron_Id = Elec_SCB_Loose;  }
      else if (TString(Elecid).Contains( "Medium" ) ) { v_electron_Id = Elec_SCB_Medium; }
      else if (TString(Elecid).Contains( "Tight"  ) ) { v_electron_Id = Elec_SCB_Tight;  }
      else if (TString(Elecid).Contains( "Veto"   ) ) { v_electron_Id = Elec_SCB_Veto;   }
      else { cout << "Electron ID Error" << endl; } 

   }
   else { cout << "Error Decaymode !! " << endl; }

   /////////////////////////////
   /// *** Jet condition *** ///
   /////////////////////////////

   if      ( TString(jetId).Contains( "PFLoose" ) ) { v_jet_Id = Jet_PFId; jet_id = 1; }// There is no loose id UL
   else if ( TString(jetId).Contains( "PFTight" ) ) { v_jet_Id = Jet_PFId; jet_id = 1; }// There is no loose id UL
   else {cout << "Jet condition error" << endl;}

   if      ( TString(jetbtag).Contains( "CSVL"  ) )    { bdisccut = 0.244; }
   else if ( TString(jetbtag).Contains( "CSVM"  ) )    { bdisccut = 0.679; }
   else if ( TString(jetbtag).Contains( "CSVT"  ) )    { bdisccut = 0.898; }
   else if ( TString(jetbtag).Contains( "CISVL" ) )    { bdisccut = 0.423; }
   else if ( TString(jetbtag).Contains( "CISVM" ) )    { bdisccut = 0.814; }
   else if ( TString(jetbtag).Contains( "CISVT" ) )    { bdisccut = 0.941; }
   else if ( TString(jetbtag).Contains( "pfCSVV2L" ) ) { bdisccut = 0.5426; }
   else if ( TString(jetbtag).Contains( "pfCSVV2M" ) ) { bdisccut = 0.8484; }
   else if ( TString(jetbtag).Contains( "pfCSVV2T" ) ) { bdisccut = 0.9535; }
   else if ( TString(jetbtag).Contains( "deepCSVL" ) ) { bdisccut = 0.2027; }
   else if ( TString(jetbtag).Contains( "deepCSVM" ) ) { bdisccut = 0.6001; }
   else if ( TString(jetbtag).Contains( "deepCSVT" ) ) { bdisccut = 0.8819; }
   else if ( TString(jetbtag).Contains( "deepJetL" ) ) { bdisccut = 0.0508; }
   else if ( TString(jetbtag).Contains( "deepJetM" ) ) { bdisccut = 0.2598; }
   else if ( TString(jetbtag).Contains( "deepJetT" ) ) { bdisccut = 0.6502; }
   else { cout << "bscriminator error !!" << endl; }
   //////////////////////////
   // *** JetCleanning ***///
   //////////////////////////

   ///////////////////////////
   // Muon for Jet Cleaning //
   ///////////////////////////

   /// Muon iso type
   if      ( TString( muiso_type ).Contains( "PFIsodbeta03" ) ) { v_mu_iso_jcl = Muon_PFIsodBeta03; }
   else if ( TString( muiso_type ).Contains( "PFIsodbeta04" ) ) { v_mu_iso_jcl = Muon_PFIsodBeta04; }
   else { cout << "Muon Iso type error for Jet cleanning" << endl;}

   /// Muon ID
   if      ( TString( muId ).Contains( "Loose" ) ) { v_mu_Id_jcl = Muon_isLoose; }
   else if ( TString( muId ).Contains( "Tight" ) ) { v_mu_Id_jcl = Muon_isTight; }
   else { cout << "Muon ID error for Jet cleanning" << endl;}
 
   ///////////////////////////////
   // Electron for Jet Cleaning //
   ///////////////////////////////

   /// Electron iso type
   if      ( TString( eliso_type ).Contains( "PFIsoRho03"   ) ) { v_el_iso_jcl = Elec_PFIsoRho03;   }
   else if ( TString( eliso_type ).Contains( "PFIsoRho04"   ) ) { v_el_iso_jcl = Elec_PFIsoRho04;   }
   else { cout << "Electron Iso type Error for Jet Cleaning" << endl;}

   // Electron ID SCB 

   if      (TString(elId).Contains( "Loose"  ) ) { v_el_Id_jcl = Elec_SCB_Loose;  }
   else if (TString(elId).Contains( "Medium" ) ) { v_el_Id_jcl = Elec_SCB_Medium; }
   else if (TString(elId).Contains( "Tight"  ) ) { v_el_Id_jcl = Elec_SCB_Tight;  }
   else if (TString(elId).Contains( "Veto"   ) ) { v_el_Id_jcl = Elec_SCB_Veto;   }
   else { cout << "Electron ID Error" << endl; }

}

void ssb_analysis::TLVInitial()
{
      Lep1      = new TLorentzVector(); // Leading Lepton 
      Lep2      = new TLorentzVector(); // Second Leading Lepton

      TMuon     = new TLorentzVector(); // Leading Lepton 
      TElectron = new TLorentzVector(); // Second Leading Lepton

      Jet1      = new TLorentzVector(); // Leading Jet 
      Jet2      = new TLorentzVector(); // Second Leading Jet

      Jet1Up      = new TLorentzVector(); // Leading Jet 
      Jet2Up      = new TLorentzVector(); // Second Leading Jet
      Jet1Dn      = new TLorentzVector(); // Leading Jet 
      Jet2Dn      = new TLorentzVector(); // Second Leading Jet
      Jet1JERUp      = new TLorentzVector(); // Leading Jet 
      Jet2JERUp      = new TLorentzVector(); // Second Leading Jet
      Jet1JERDn      = new TLorentzVector(); // Leading Jet 
      Jet2JERDn      = new TLorentzVector(); // Second Leading Jet

      bJet1     = new TLorentzVector(); // Leading Jet 
      bJet2     = new TLorentzVector(); // Second Leading Jet
      bJet1Up     = new TLorentzVector(); // Leading Jet 
      bJet2Up     = new TLorentzVector(); // Second Leading Jet
      bJet1Dn     = new TLorentzVector(); // Leading Jet 
      bJet2Dn     = new TLorentzVector(); // Second Leading Jet

      bJet1JERUp     = new TLorentzVector(); // Leading Jet 
      bJet2JERUp     = new TLorentzVector(); // Second Leading Jet
      bJet1JERDn     = new TLorentzVector(); // Leading Jet 
      bJet2JERDn     = new TLorentzVector(); // Second Leading Jet

      Met       = new TLorentzVector(); // MET
      MetJESUp  = new TLorentzVector(); // MET
      MetJESDn  = new TLorentzVector(); // MET
      MetJERUp  = new TLorentzVector(); // MET
      MetJERDn  = new TLorentzVector(); // MET
  
      Nu1       = new TLorentzVector(); // Nutrino 1 ( For Lep1 )
      Nu2       = new TLorentzVector(); // Nutrino 2 ( For Lep2 )

      W1        = new TLorentzVector(); // W boson 1 ( For Lep1 )  
      W2        = new TLorentzVector(); // W boson 2 ( For Lep2 )

      Top       = new TLorentzVector(); // Top
      AnTop     = new TLorentzVector(); // Anti-top

      Top1      = new TLorentzVector(); // Top Leading Pt
      Top2      = new TLorentzVector(); // Top Second Leading Pt 

      bJet      = new TLorentzVector(); // bJet 
      AnbJet    = new TLorentzVector(); // Anti-bJet

      Lep       = new TLorentzVector(); // Lepton 
      AnLep     = new TLorentzVector(); // Anti-Lepton

      Nu       = new TLorentzVector(); // Lepton 
      AnNu     = new TLorentzVector(); // Anti-Lepton

      GenLep       = new TLorentzVector(); // Lepton 
      GenAnLep     = new TLorentzVector(); // Anti-Lepton

      // Lepton selection
      Mu_lep_sel  = new TLorentzVector();
      Ele_lep_sel = new TLorentzVector();

      // Veto muon and electron
      vetoMu  = new TLorentzVector();
      vetoEle = new TLorentzVector();

      // Jet
      Mu_jet_clean           = new TLorentzVector();     
      Ele_jet_clean          = new TLorentzVector();
      selected_Ele_jet_clean = new TLorentzVector();
      selected_Mu_jet_clean  = new TLorentzVector();
      TJet                   = new TLorentzVector();
}


void ssb_analysis::Loop( char *logfile )
{

   //GetTotalEvent();
   //////////
   if (fChain == 0) return;
   //////////

   //////////
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   //////////

   ///My variables
   Long64_t __tot_evt = 0;

   /// Check Total Event
   cout << "Ntuple Total Event Check !! " << NtupletotalEvent << endl;

   MCSF();
   cout << "MC_SF Check !! " << mc_sf_ << endl;

   Np_eventw_2 = 0.0;
   Nm_eventw_2 = 0.0;
   Num_cpo3_p  = 0; 
   Num_cpo3_m  = 0; 

   int count_emu = 0;
  
   ////////////////////////
   /// start event loop ///
   ////////////////////////

   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
      {
         printf("ERROR: Could not load tree!!!\n");
         break;
      }

      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
      if (jentry % 10000 == 0) printf("Event %lld\n", jentry); //%lld supports Long64_t

      __tot_evt++;
      ////////////////////////////////////////
      /// start Main Loop Function and Cut ///
      ////////////////////////////////////////

      // Make TL to use SSBTree //
      //MakeVecforTL();

      // initailizing TLorentzVector
      TLVInitial();

      evt_weight_ = 1;

      GetVariables();
      //////////////////////

      ///////////////////////////////////
      /// ** Lepton + Jet Analysis ** ///
      ///////////////////////////////////
      if ( TString(Decaymode).Contains( "dielec" ) ||
           TString(Decaymode).Contains( "dimuon" ) || 
           TString(Decaymode).Contains( "muel" )     ){

         /// Strat !///
         NumPVCount(); 
         LeptonSelector();
         LeptonOrder();
         JetSelector();
         METDefiner();
         BDsicApply();
         BJetDefiner();
         FillHisto(h_Num_PV_BeforePreSel, num_pv, evt_weight_);
         //////////////////////
         //// Event Filter ////
         //////////////////////

         if ( METFilterAPP() == false ) {continue;}
         FillHisto(h_Num_PV_AfterMetFilter, num_pv, evt_weight_);
         //if ( Filter_PV->at(0) == false ) {continue;}
         ////////////////////////////////
         /// Step 0 Trigger Selection ///
         ////////////////////////////////
         if ( Trigger() == false ) {continue;}
         FillHisto(h_Num_PV_AfterTrigger, num_pv, evt_weight_);
         if ( NumIsoLeptons() == false ){ continue;}
         FillHisto( h_EventWeight[0] , evt_weight_  );
         FillHisto( h_cf_NLeptons[0], Muon_Count, evt_weight_);
         FillHisto( h_cf_NJets[0], v_jet_idx.size(), evt_weight_ );
         FillHisto( h_cf_NPV[0]    , num_pv, evt_weight_ );
         FillHisto( h_Num_PV[0]    , num_pv, evt_weight_ );



         if ( ThirdLeptonVeto() == false ){continue;}
         if (LeptonsPtAddtional() == false ) {continue;}
         if ( DiLeptonMassCut() == false ) {continue;}
//         LeptonSFApply();

         FillHisto( h_EventWeight[1] , evt_weight_  );
         FillHisto( h_cf_NPV[1]    , num_pv        , evt_weight_ );
         FillHisto( h_cf_NLeptons[1], v_lepton_idx.size(), evt_weight_ );
         FillHisto( h_cf_Lep1pt[1] , (*Lep1).Pt()  , evt_weight_ );
         FillHisto( h_cf_Lep1eta[1], (*Lep1).Eta() , evt_weight_ );
         FillHisto( h_cf_Lep1phi[1], (*Lep1).Phi() , evt_weight_ );
         FillHisto( h_cf_Lep2pt[1] , (*Lep2).Pt()  , evt_weight_ );
         FillHisto( h_cf_Lep2eta[1], (*Lep2).Eta() , evt_weight_ );
         FillHisto( h_cf_Lep2phi[1], (*Lep2).Phi() , evt_weight_ );
         FillHisto( h_cf_dilep_inv_mass[1], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
         FillHisto( h_cf_NJets[1]  , v_jet_idx.size(), evt_weight_ );

         if(v_jet_idx.size() > 0 )
         {
            FillHisto( h_cf_Jet1pt[1] , (*Jet1).Pt()  , evt_weight_ );
            FillHisto( h_cf_Jet1eta[1], (*Jet1).Eta() , evt_weight_ );
            FillHisto( h_cf_Jet1phi[1], (*Jet1).Phi() , evt_weight_ );
            if (v_jet_idx.size() > 1)
            {
               FillHisto( h_cf_Jet2pt[1] , (*Jet2).Pt()  , evt_weight_ );
               FillHisto( h_cf_Jet2eta[1], (*Jet2).Eta() , evt_weight_ );
               FillHisto( h_cf_Jet2phi[1], (*Jet2).Phi() , evt_weight_ );
            }
         }

         FillHisto( h_cf_metpt[1] ,Met->Pt()  , evt_weight_ );
         FillHisto( h_cf_metphi[1],Met->Phi() , evt_weight_ );
         FillHisto( h_DiLepMass[1], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
         FillHisto( h_Num_PV[1], num_pv, evt_weight_ );
         FillHisto( h_Lep1pt[1] , (*Lep1).Pt()  , evt_weight_ );
         FillHisto( h_Lep1eta[1], (*Lep1).Eta() , evt_weight_ );
         FillHisto( h_Lep1phi[1], (*Lep1).Phi() , evt_weight_ );
         FillHisto( h_Lep2pt[1] , (*Lep2).Pt()  , evt_weight_ );
         FillHisto( h_Lep2eta[1], (*Lep2).Eta() , evt_weight_ );
         FillHisto( h_Lep2phi[1], (*Lep2).Phi() , evt_weight_ );
         FillHisto( h_METpt[1]   , Met->Pt()  , evt_weight_ );
         FillHisto( h_METphi[1]  , Met->Phi()  , evt_weight_ );
         FillHisto( h_Num_Jets[1]  , v_jet_idx.size(), evt_weight_ );
         FillHisto( h_Num_bJets[1], nbtagged, evt_weight_ );
         if(v_jet_idx.size() > 0 )
         {
            FillHisto( h_Jet1pt[1] , (*Jet1).Pt()  , evt_weight_ );
            FillHisto( h_Jet1eta[1], (*Jet1).Eta() , evt_weight_ );
            FillHisto( h_Jet1phi[1], (*Jet1).Phi() , evt_weight_ );
            if (v_jet_idx.size() > 1)
            {
               FillHisto( h_Jet2pt[1] , (*Jet2).Pt()  , evt_weight_ );
               FillHisto( h_Jet2eta[1], (*Jet2).Eta() , evt_weight_ );
               FillHisto( h_Jet2phi[1], (*Jet2).Phi() , evt_weight_ );
            }
         }
         ///////////////////////
         /// Z mass Veto Cut ///
         ///////////////////////
         if ( ZVetoCut() == false ) {continue;}


         FillHisto( h_EventWeight[2], evt_weight_  );
         FillHisto( h_cf_NLeptons[2], v_lepton_idx.size(), evt_weight_ );
         FillHisto( h_cf_Lep1pt[2] , (*Lep1).Pt()  , evt_weight_ );
         FillHisto( h_cf_Lep1eta[2], (*Lep1).Eta() , evt_weight_ );
         FillHisto( h_cf_Lep1phi[2], (*Lep1).Phi() , evt_weight_ );
         FillHisto( h_cf_Lep2pt[2] , (*Lep2).Pt()  , evt_weight_ );
         FillHisto( h_cf_Lep2eta[2], (*Lep2).Eta() , evt_weight_ );
         FillHisto( h_cf_Lep2phi[2], (*Lep2).Phi() , evt_weight_ );
         FillHisto( h_cf_NPV[2]    , num_pv        , evt_weight_ );
         FillHisto( h_cf_NJets[2]  , v_jet_idx.size(), evt_weight_ );
         FillHisto( h_cf_dilep_inv_mass[2], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
         FillHisto( h_cf_metpt[2] , Met->Pt() , evt_weight_ );
         FillHisto( h_cf_metphi[2], Met->Phi(), evt_weight_ );
         if(v_jet_idx.size() > 0 )
         {
            FillHisto( h_cf_Jet1pt[2] , (*Jet1).Pt()  , evt_weight_ );
            FillHisto( h_cf_Jet1eta[2], (*Jet1).Eta() , evt_weight_ );
            FillHisto( h_cf_Jet1phi[2], (*Jet1).Phi() , evt_weight_ );
            if (v_jet_idx.size() > 1)
            {
               FillHisto( h_cf_Jet2pt[2] , (*Jet2).Pt()  , evt_weight_ );
               FillHisto( h_cf_Jet2eta[2], (*Jet2).Eta() , evt_weight_ );
               FillHisto( h_cf_Jet2phi[2], (*Jet2).Phi() , evt_weight_ );
            }
         }

         FillHisto( h_Num_PV[2], num_pv, evt_weight_ );
         FillHisto( h_DiLepMass[2], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
         FillHisto( h_Lep1pt[2] , (*Lep1).Pt()  , evt_weight_ );
         FillHisto( h_Lep1eta[2], (*Lep1).Eta() , evt_weight_ );
         FillHisto( h_Lep1phi[2], (*Lep1).Phi() , evt_weight_ );
         FillHisto( h_Lep2pt[2] , (*Lep2).Pt()  , evt_weight_ );
         FillHisto( h_Lep2eta[2], (*Lep2).Eta() , evt_weight_ );
         FillHisto( h_Lep2phi[2], (*Lep2).Phi() , evt_weight_ );
         FillHisto( h_Num_Jets[2]  , v_jet_idx.size(), evt_weight_ );
         FillHisto( h_Num_bJets[2], nbtagged, evt_weight_ );
         FillHisto( h_METpt[2]   , Met->Pt()  , evt_weight_ );
         FillHisto( h_METphi[2]  , Met->Phi()  , evt_weight_ );

         if(v_jet_idx.size() > 0 )
         {
            FillHisto( h_Jet1pt[2] , (*Jet1).Pt()  , evt_weight_ );
            FillHisto( h_Jet1eta[2], (*Jet1).Eta() , evt_weight_ );
            FillHisto( h_Jet1phi[2], (*Jet1).Phi() , evt_weight_ );
            if (v_jet_idx.size() > 1)
            {
               FillHisto( h_Jet2pt[2] , (*Jet2).Pt()  , evt_weight_ );
               FillHisto( h_Jet2eta[2], (*Jet2).Eta() , evt_weight_ );
               FillHisto( h_Jet2phi[2], (*Jet2).Phi() , evt_weight_ );
            }
         }
         JetDefiner();
         /// Step 3. Num Jet >=2 ///
         if ( NumJetCut(v_jet_idx) == true ){
            FillHisto( h_EventWeight[3] , evt_weight_  );
            FillHisto( h_cf_NLeptons[3], v_lepton_idx.size(), evt_weight_ );
            FillHisto( h_cf_Lep1pt[3] , (*Lep1).Pt()  , evt_weight_ );
            FillHisto( h_cf_Lep1eta[3], (*Lep1).Eta() , evt_weight_ );
            FillHisto( h_cf_Lep1phi[3], (*Lep1).Phi() , evt_weight_ );
            FillHisto( h_cf_Lep2pt[3] , (*Lep2).Pt()  , evt_weight_ );
            FillHisto( h_cf_Lep2eta[3], (*Lep2).Eta() , evt_weight_ );
            FillHisto( h_cf_Lep2phi[3], (*Lep2).Phi() , evt_weight_ );
            FillHisto( h_cf_NPV[3]    , num_pv        , evt_weight_ );
            FillHisto( h_cf_NJets[3]  , v_jet_idx.size(), evt_weight_ );
            FillHisto( h_cf_Jet1pt[3] , (*Jet1).Pt()  , evt_weight_ );
            FillHisto( h_cf_Jet1eta[3], (*Jet1).Eta() , evt_weight_ );
            FillHisto( h_cf_Jet1phi[3], (*Jet1).Phi() , evt_weight_ );
            FillHisto( h_cf_Jet2pt[3] , (*Jet2).Pt()  , evt_weight_ );
            FillHisto( h_cf_Jet2eta[3], (*Jet2).Eta() , evt_weight_ );
            FillHisto( h_cf_Jet2phi[3], (*Jet2).Phi() , evt_weight_ );
            FillHisto( h_cf_metpt[3] ,Met->Pt()  , evt_weight_ );
            FillHisto( h_cf_metphi[3],Met->Phi() , evt_weight_ );

            FillHisto( h_cf_dilep_inv_mass[3], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
            FillHisto( h_cf_metpt[3] , Met->Pt() , evt_weight_ );
            FillHisto( h_cf_metphi[3], Met->Phi(), evt_weight_ );

            FillHisto( h_Num_PV[3], num_pv, evt_weight_ );
            FillHisto( h_DiLepMass[3], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
            FillHisto( h_Lep1pt[3] , (*Lep1).Pt()  , evt_weight_ );
            FillHisto( h_Lep1eta[3], (*Lep1).Eta() , evt_weight_ );
            FillHisto( h_Lep1phi[3], (*Lep1).Phi() , evt_weight_ );
            FillHisto( h_Lep2pt[3] , (*Lep2).Pt()  , evt_weight_ );
            FillHisto( h_Lep2eta[3], (*Lep2).Eta() , evt_weight_ );
            FillHisto( h_Lep2phi[3], (*Lep2).Phi() , evt_weight_ );
            FillHisto( h_Jet1pt[3] , (*Jet1).Pt()  , evt_weight_ );
            FillHisto( h_Jet1eta[3], (*Jet1).Eta() , evt_weight_ );
            FillHisto( h_Jet1phi[3], (*Jet1).Phi() , evt_weight_ );
            FillHisto( h_Jet2pt[3] , (*Jet2).Pt()  , evt_weight_ );
            FillHisto( h_Jet2eta[3], (*Jet2).Eta() , evt_weight_ );
            FillHisto( h_Jet2phi[3], (*Jet2).Phi() , evt_weight_ );
            FillHisto( h_Num_Jets[3], v_jet_idx.size(), evt_weight_ );
            FillHisto( h_Num_bJets[3], nbtagged, evt_weight_ );
            FillHisto( h_METpt[3]   , Met->Pt()  , evt_weight_ );
            FillHisto( h_METphi[3]  , Met->Phi()  , evt_weight_ );
   
            if ( METCut(Met) == true ) {
               FillHisto( h_EventWeight[4] , evt_weight_  );
               FillHisto( h_cf_NLeptons[4], v_lepton_idx.size(), evt_weight_ );

               FillHisto( h_cf_Lep1pt[4] , (*Lep1).Pt()  , evt_weight_ );
               FillHisto( h_cf_Lep1eta[4], (*Lep1).Eta() , evt_weight_ );
               FillHisto( h_cf_Lep1phi[4], (*Lep1).Phi() , evt_weight_ );
               FillHisto( h_cf_Lep2pt[4] , (*Lep2).Pt()  , evt_weight_ );
               FillHisto( h_cf_Lep2eta[4], (*Lep2).Eta() , evt_weight_ );
               FillHisto( h_cf_Lep2phi[4], (*Lep2).Phi() , evt_weight_ );
               FillHisto( h_cf_NPV[4]    , num_pv        , evt_weight_ );
               FillHisto( h_cf_NJets[4]  , v_jet_idx.size(), evt_weight_ );
               FillHisto( h_cf_Jet1pt[4] , (*Jet1).Pt()  , evt_weight_ );
               FillHisto( h_cf_Jet1eta[4], (*Jet1).Eta() , evt_weight_ );
               FillHisto( h_cf_Jet1phi[4], (*Jet1).Phi() , evt_weight_ );
               FillHisto( h_cf_Jet2pt[4] , (*Jet2).Pt()  , evt_weight_ );
               FillHisto( h_cf_Jet2eta[4], (*Jet2).Eta() , evt_weight_ );
               FillHisto( h_cf_Jet2phi[4], (*Jet2).Phi() , evt_weight_ );

               FillHisto( h_cf_dilep_inv_mass[4], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
               FillHisto( h_cf_metpt[4] , Met->Pt() , evt_weight_ );
               FillHisto( h_cf_metphi[4], Met->Phi(), evt_weight_ );

               FillHisto( h_Lep1pt[4]  , Lep1->Pt() , evt_weight_ );
               FillHisto( h_Lep2pt[4]  , Lep2->Pt() , evt_weight_ );
               FillHisto( h_Lep1eta[4] , Lep1->Eta(), evt_weight_ );
               FillHisto( h_Lep2eta[4] , Lep2->Eta(), evt_weight_ );
               FillHisto( h_Lep1phi[4] , Lep1->Phi(), evt_weight_ );
               FillHisto( h_Lep2phi[4] , Lep2->Phi(), evt_weight_ );

               FillHisto( h_Jet1pt[4]  , Jet1->Pt() , evt_weight_ );
               FillHisto( h_Jet2pt[4]  , Jet2->Pt() , evt_weight_ );
               FillHisto( h_Jet1eta[4] , Jet1->Eta(), evt_weight_ );
               FillHisto( h_Jet2eta[4] , Jet2->Eta(), evt_weight_ );
               FillHisto( h_Jet1phi[4] , Jet1->Phi(), evt_weight_ );
               FillHisto( h_Jet2phi[4] , Jet2->Phi(), evt_weight_ );
               FillHisto( h_METpt[4]   , Met->Pt()  , evt_weight_ );
               FillHisto( h_METphi[4]  , Met->Phi() , evt_weight_ );

               FillHisto( h_DiLepMass[4], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
               FillHisto( h_Num_PV[4], num_pv, evt_weight_ );
               FillHisto( h_Num_Jets[4], v_jet_idx.size(), evt_weight_ );
               FillHisto( h_Num_bJets[4], nbtagged, evt_weight_ );
               if ( BJetCut(v_bjet_idx) == true ){
                  FillHisto(h_bTagWeight,evt_weight_/evt_weight_beforeBtag_);
                  FillHisto( h_EventWeight[5] , evt_weight_  );
                  FillHisto( h_cf_NLeptons[5], v_lepton_idx.size(), evt_weight_ );

                  FillHisto( h_cf_Lep1pt[5] , (*Lep1).Pt()  , evt_weight_ );
                  FillHisto( h_cf_Lep1eta[5], (*Lep1).Eta() , evt_weight_ );
                  FillHisto( h_cf_Lep1phi[5], (*Lep1).Phi() , evt_weight_ );
                  FillHisto( h_cf_Lep2pt[5] , (*Lep2).Pt()  , evt_weight_ );
                  FillHisto( h_cf_Lep2eta[5], (*Lep2).Eta() , evt_weight_ );
                  FillHisto( h_cf_Lep2phi[5], (*Lep2).Phi() , evt_weight_ );
                  FillHisto( h_cf_NPV[5]    , num_pv        , evt_weight_ );
                  FillHisto( h_cf_NJets[5]  , v_jet_idx.size(), evt_weight_ );
                  FillHisto( h_cf_Jet1pt[5] , (*Jet1).Pt()  , evt_weight_ );
                  FillHisto( h_cf_Jet1eta[5], (*Jet1).Eta() , evt_weight_ );
                  FillHisto( h_cf_Jet1phi[5], (*Jet1).Phi() , evt_weight_ );
                  FillHisto( h_cf_Jet2pt[5] , (*Jet2).Pt()  , evt_weight_ );
                  FillHisto( h_cf_Jet2eta[5], (*Jet2).Eta() , evt_weight_ );
                  FillHisto( h_cf_Jet2phi[5], (*Jet2).Phi() , evt_weight_ );

                  FillHisto( h_cf_dilep_inv_mass[5], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
                  FillHisto( h_cf_metpt[5] , Met->Pt() , evt_weight_ );
                  FillHisto( h_cf_metphi[5], Met->Phi(), evt_weight_ );


                  FillHisto( h_Lep1pt[5]  , Lep1->Pt() , evt_weight_ );
                  FillHisto( h_Lep2pt[5]  , Lep2->Pt() , evt_weight_ );
                  FillHisto( h_Lep1eta[5] , Lep1->Eta(), evt_weight_ );
                  FillHisto( h_Lep2eta[5] , Lep2->Eta(), evt_weight_ );
                  FillHisto( h_Lep1phi[5] , Lep1->Phi(), evt_weight_ );
                  FillHisto( h_Lep2phi[5] , Lep2->Phi(), evt_weight_ );

                  if (TString(Decaymode).Contains("muel"))
                  {
                     FillHisto( h_Muonpt[5]  , TMuon->Pt()      , evt_weight_ );
                     FillHisto( h_Elecpt[5]  , TElectron->Pt()  , evt_weight_ );
                     FillHisto( h_Muoneta[5] , TMuon->Eta()     , evt_weight_ );
                     FillHisto( h_Eleceta[5] , TElectron->Eta() , evt_weight_ );
                     FillHisto( h_Muonphi[5] , TMuon->Phi()     , evt_weight_ );
                     FillHisto( h_Elecphi[5] , TElectron->Phi() , evt_weight_ );
                  }

                  FillHisto( h_Jet1pt[5]  , Jet1->Pt() , evt_weight_ );
                  FillHisto( h_Jet2pt[5]  , Jet2->Pt() , evt_weight_ );
                  FillHisto( h_Jet1eta[5] , Jet1->Eta(), evt_weight_ );
                  FillHisto( h_Jet2eta[5] , Jet2->Eta(), evt_weight_ );
                  FillHisto( h_Jet1phi[5] , Jet1->Phi(), evt_weight_ );
                  FillHisto( h_Jet2phi[5] , Jet2->Phi(), evt_weight_ );
                  FillHisto( h_METpt[5]   , Met->Pt()  , evt_weight_ );
                  FillHisto( h_METphi[5]  , Met->Phi() , evt_weight_ );

                  FillHisto( h_DiLepMass[5], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );

                  FillHisto( h_Num_PV[5], num_pv, evt_weight_ );
                  FillHisto( h_Num_Jets[5], v_jet_idx.size(), evt_weight_ );
                  FillHisto( h_Num_bJets[5], nbtagged, evt_weight_ );

                  //FillHisto( h_HT[5], AllJetpt, evt_weight_);


               } /// Step 5///
            } /// Step 4 ///
         }/// Step 3 ///

      }
      else if ( TString(Decaymode).Contains( "muonJet" )  ) 
      {
         cout << " -- We don't have Muon + Jet code ---" << endl;   
      }
      else {cout << "Decaymode Error : " << Decaymode << endl; }
   }//event loop
  
   printf("Total processed number of events: %lld\n", __tot_evt);


}//end Loop function

void ssb_analysis::GetNtupleTotalEvent( unsigned int totevent )// Not Used Function
{
   NtupletotalEvent = totevent;
}

void ssb_analysis::Start()
{
   fout = new TFile(Form("output/%s",outfile),"RECREATE");
//   else if      ( genLoopon == 1 ){ fout = new TFile(Form("output/%s",outfile),"UPDATE");}
   //if      ( genLoopon == 0 ){ fout = new TFile(Form("gsidcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/SSB_CPviolation/output/%s",outfile),"RECREATE");}
   //else if ( genLoopon == 1 ){ fout = new TFile(Form("gsidcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/SSB_CPviolation/output/%s",outfile),"UPDATE"  );}
//   else {cout << "genLoopon error" << endl;}
   fout->cd("");

   TDirectory *dir = gDirectory;
   dir->cd();

   DeclareHistos();
}

void ssb_analysis::DeclareHistos()
{
   for (int i =0 ; i < 10 ; i++)
   {
      h_cf_NLeptons[i]       = new TH1D(Form("_h_cf_NLeptons_%d_"      , i),Form("Num. Lepton %s"                ,cutflowName[i].Data()), 20  , 0    , 20  ); h_cf_NLeptons[i]->Sumw2();         
      h_cf_Lep1pt[i]         = new TH1D(Form("_h_cf_Lep1pt_%d_"        , i),Form("Leading Lepton pT %s"          ,cutflowName[i].Data()), 1000, 0    , 1000); h_cf_Lep1pt[i]->Sumw2();           
      h_cf_Lep1phi[i]        = new TH1D(Form("_h_cf_Lep1phi_%d_"       , i),Form("Leading Lepton phi %s"        ,cutflowName[i].Data()), 24  , -1*pi, 1*pi); h_cf_Lep1phi[i]->Sumw2();          
      h_cf_Lep1eta[i]        = new TH1D(Form("_h_cf_Lep1eta_%d_"       , i),Form("Leading Lepton eta %s"        ,cutflowName[i].Data()), 50  , -2.5 , 2.5 ); h_cf_Lep1eta[i]->Sumw2();          
      h_cf_Lep2pt[i]         = new TH1D(Form("_h_cf_Lep2pt_%d_"        , i),Form("Second Leading Lepton pT %s"   ,cutflowName[i].Data()), 1000, 0    , 1000); h_cf_Lep2pt[i]->Sumw2();           
      h_cf_Lep2phi[i]        = new TH1D(Form("_h_cf_Lep2phi_%d_"       , i),Form("Second Leading Lepton phi %s"  ,cutflowName[i].Data()), 24  , -1*pi, 1*pi); h_cf_Lep2phi[i]->Sumw2();      
      h_cf_Lep2eta[i]        = new TH1D(Form("_h_cf_Lep2eta_%d_"       , i),Form("Second Leading Lepton eta %s"  ,cutflowName[i].Data()), 50  , -2.5 , 2.5 ); h_cf_Lep2eta[i]->Sumw2();       
      h_cf_dilep_inv_mass[i] = new TH1D(Form("_h_cf_dilep_inv_mass_%d_", i),Form("Invariant Mass of Dilepton %s" ,cutflowName[i].Data()), 1000, 0    , 1000); h_cf_dilep_inv_mass[i]->Sumw2();
      h_cf_Jet1pt[i]         = new TH1D(Form("_h_cf_Jet1pt_%d_"        , i),Form("Leading Jet pT %s"             ,cutflowName[i].Data()), 1000, 0    , 1000); h_cf_Jet1pt[i]->Sumw2();       
      h_cf_Jet1phi[i]        = new TH1D(Form("_h_cf_Jet1phi_%d_"       , i),Form("Leading Jet phi %s"            ,cutflowName[i].Data()), 24  , -1*pi, 1*pi); h_cf_Jet1phi[i]->Sumw2();       
      h_cf_Jet1eta[i]        = new TH1D(Form("_h_cf_Jet1eta_%d_"       , i),Form("Leading Jet eta %s"            ,cutflowName[i].Data()), 50  , -2.5 , 2.5 ); h_cf_Jet1eta[i]->Sumw2();       
      h_cf_Jet2pt[i]         = new TH1D(Form("_h_cf_Jet2pt_%d_"        , i),Form("Second Leading Jet pT %s"      ,cutflowName[i].Data()), 1000, 0    , 1000); h_cf_Jet2pt[i]->Sumw2();        
      h_cf_Jet2phi[i]        = new TH1D(Form("_h_cf_Jet2phi_%d_"       , i),Form("Second Leading Jet phi %s"     ,cutflowName[i].Data()), 24  , -1*pi, 1*pi); h_cf_Jet2phi[i]->Sumw2();       
      h_cf_Jet2eta[i]        = new TH1D(Form("_h_cf_Jet2eta_%d_"       , i),Form("Second Leading Jet eta %s"     ,cutflowName[i].Data()), 50  , -2.5 , 2.5 ); h_cf_Jet2eta[i]->Sumw2();      
      h_cf_NJets[i]          = new TH1D(Form("_h_cf_NJets_%d_"         , i),Form("Num. of Jets %s"               ,cutflowName[i].Data()), 20  , 0    , 20  ); h_cf_NJets[i]->Sumw2();        
      h_cf_metpt[i]          = new TH1D(Form("_h_cf_metpt_%d_"         , i),Form("MET %s"                        ,cutflowName[i].Data()), 1000, 0    , 1000); h_cf_metpt[i]->Sumw2();        
      h_cf_metphi[i]         = new TH1D(Form("_h_cf_metphi_%d_"        , i),Form("MET Phi %s"                    ,cutflowName[i].Data()), 24  , -1*pi, 1*pi); h_cf_metphi[i]->Sumw2();       
      h_cf_Nbjets[i]         = new TH1D(Form("_h_cf_Nbjets_%d_"        , i),Form("Num. of b-jets %s"             ,cutflowName[i].Data()), 20  , 0    , 20  ); h_cf_Nbjets[i]->Sumw2();        
      h_cf_NPV[i]            = new TH1D(Form("_h_cf_NPV_%d_"           , i),Form("Num. of Primary Vertex %s"     ,cutflowName[i].Data()), 100 , 0    , 100 ); h_cf_NPV[i]->Sumw2();          


      h_EventWeight[i] = new TH1D(Form("h_EventWeight_%d",i), Form("Event Weight %s",cutflowName[i].Data()), 1000, -5, 5); h_EventWeight[i]->Sumw2();

      h_Lep1pt[i]  = new TH1D(Form("h_Lep1pt_%d" ,i), Form("Leading Lepton pT %s"        ,cutflowName[i].Data()), 250, 0.0, 250); h_Lep1pt[i]->Sumw2(); 
      h_Lep2pt[i]  = new TH1D(Form("h_Lep2pt_%d" ,i), Form("Second Leading Lepton pT %s" ,cutflowName[i].Data()), 250, 0.0, 250); h_Lep2pt[i]->Sumw2();
      h_Lep1eta[i] = new TH1D(Form("h_Lep1eta_%d",i), Form("Leading Lepton Eta    %s"    ,cutflowName[i].Data()), 50, -2.5, 2.5); h_Lep1eta[i]->Sumw2();
      h_Lep2eta[i] = new TH1D(Form("h_Lep2eta_%d",i), Form("Second Leading Lepton Eta %s",cutflowName[i].Data()), 50, -2.5, 2.5); h_Lep2eta[i]->Sumw2();
      h_Lep1phi[i] = new TH1D(Form("h_Lep1phi_%d",i), Form("Leading Lepton Phi %s"       ,cutflowName[i].Data()), 24, -1*pi, pi); h_Lep1phi[i]->Sumw2();
      h_Lep2phi[i] = new TH1D(Form("h_Lep2phi_%d",i), Form("Second Leading Lepton Phi %s",cutflowName[i].Data()), 24, -1*pi, pi); h_Lep2phi[i]->Sumw2();
   
      h_Muonpt[i]  = new TH1D(Form("h_Muonpt_%d",i ), Form("Muon pT %s"         ,cutflowName[i].Data()), 250, 0.0, 250); h_Muonpt[i]->Sumw2();
      h_Elecpt[i]  = new TH1D(Form("h_Elecpt_%d",i ), Form("Electron pT %s"     ,cutflowName[i].Data()), 250, 0.0, 250); h_Elecpt[i]->Sumw2();
      h_Muoneta[i] = new TH1D(Form("h_Muoneta_%d",i), Form("Muon Eta %s"        ,cutflowName[i].Data()), 50, -2.5, 2.5); h_Muoneta[i]->Sumw2();
      h_Eleceta[i] = new TH1D(Form("h_Eleceta_%d",i), Form("Electron Eta %s"    ,cutflowName[i].Data()), 50, -2.5, 2.5); h_Eleceta[i]->Sumw2();
      h_Muonphi[i] = new TH1D(Form("h_Muonphi_%d",i), Form("Muon Phi %s"        ,cutflowName[i].Data()), 24, -1*pi, pi); h_Muonphi[i]->Sumw2();
      h_Elecphi[i] = new TH1D(Form("h_Elecphi_%d",i), Form("Electron Phi %s"    ,cutflowName[i].Data()), 24, -1*pi, pi); h_Elecphi[i]->Sumw2();
   
      h_Jet1pt[i]  = new TH1D(Form("h_Jet1pt_%d", i), Form("Leading Jet pT %s"        ,cutflowName[i].Data()), 250, 0.0, 250); h_Jet1pt[i]->Sumw2();
      h_Jet2pt[i]  = new TH1D(Form("h_Jet2pt_%d", i), Form("Second Leading Jet pT %s" ,cutflowName[i].Data()), 250, 0.0, 250); h_Jet2pt[i]->Sumw2();
      h_Jet1eta[i] = new TH1D(Form("h_Jet1eta_%d",i), Form("Leading Jet Eta %s"       ,cutflowName[i].Data()), 50, -2.5, 2.5); h_Jet1eta[i]->Sumw2();
      h_Jet2eta[i] = new TH1D(Form("h_Jet2eta_%d",i), Form("Second Leading Jet Eta %s",cutflowName[i].Data()), 50, -2.5, 2.5); h_Jet2eta[i]->Sumw2();
      h_Jet1phi[i] = new TH1D(Form("h_Jet1phi_%d",i), Form("Leading Jet Phi %s"       ,cutflowName[i].Data()), 24, -1*pi, pi); h_Jet1phi[i]->Sumw2();
      h_Jet2phi[i] = new TH1D(Form("h_Jet2phi_%d",i), Form("Second Leading Jet Phi %s",cutflowName[i].Data()), 24, -1*pi, pi); h_Jet2phi[i]->Sumw2();
    
     
      h_METpt[i]  = new TH1D(Form("h_METpt_%d",i), Form("MET pT %s" ,cutflowName[i].Data()), 200, 0.0, 200); h_METpt[i]->Sumw2();
      h_METphi[i] = new TH1D(Form("h_METphi_%d",i), Form("MET Phi %s",cutflowName[i].Data()), 24, -1*pi, pi); h_METphi[i]->Sumw2();

      h_DiLepMass[i] = new TH1D(Form("h_DiLepMass_%d",i),Form("Di-Lepton Invariant Mass %s",cutflowName[i].Data()), 300, 0.0, 300); h_DiLepMass[i]->Sumw2();
      h_Num_PV[i]    = new TH1D(Form("h_Num_PV_%d",i),     Form("Num of Primary Vertex after %s",cutflowName[i].Data()), 100, 0.0, 100); h_Num_PV[i]->Sumw2();
      h_Num_Jets[i]  = new TH1D(Form("h_Num_Jets_%d",i), Form("Num. of Jets after %s",cutflowName[i].Data()), 20, 0.0, 20); h_Num_Jets[i]->Sumw2();
      h_Num_bJets[i] = new TH1D(Form("h_Num_bJets_%d",i),Form("Num. of b Jets after %s",cutflowName[i].Data()), 20, 0.0, 20); h_Num_bJets[i]->Sumw2();
      
      if(i>3){h_HT[i]      = new TH1D(Form("h_HT_%d",i), Form("Sum . pT of All Jet %s",cutflowName[i].Data()), 1000, 0.0, 1000); h_HT[i]->Sumw2();}
   }
   for (int i =0; i < 5; ++i)
   {
      h_Topo_Apla[i] = new TH1D(Form("h_Topo_Apla_%d",i), Form("Topologycal Vari. Aplanarity after %s",cutflowName[i+5].Data()), 300, -1, 2); h_Topo_Apla[i]->Sumw2();
      h_Topo_Sphe[i] = new TH1D(Form("h_Topo_Sphe_%d",i), Form("Topologycal Vari. Sphericity after %s",cutflowName[i+5].Data()), 300, -1, 2); h_Topo_Sphe[i]->Sumw2();
      h_Topo_Plan[i] = new TH1D(Form("h_Topo_Plan_%d",i), Form("Topologycal Vari. Planarity after %s",cutflowName[i+5].Data()), 300, -1, 2); h_Topo_Plan[i]->Sumw2();
      h_bTagEff[i]  = new TH1D(Form("h_bTagEff_%d",i), Form("BTagging Weight %s",cutflowName[i+5].Data()), 1000, -5, 5); h_bTagEff[i]->Sumw2();
   }
   h_PileUp      = new TH1D(Form("h_PileUp"), Form("hPileUp"), 1000, 0.0, 100); h_PileUp->Sumw2();// For Closure test // 
   h_PileUp_Up   = new TH1D(Form("h_PileUp_Up"), Form("hPileUp_Up"), 1000, 0.0, 100); h_PileUp_Up->Sumw2();// For Closure test // 
   h_PileUp_Down = new TH1D(Form("h_PileUp_Down"), Form("hPileUp_Down"), 1000, 0.0, 100); h_PileUp_Down->Sumw2();// For Closure test // 

   h_NumEl  = new TH1D(Form("h_NumEl" ), Form("hNumEl" ), 10, 0.0, 10); h_NumEl->Sumw2(); 
   h_NumMu  = new TH1D(Form("h_NumMu" ), Form("hNumMu" ), 10, 0.0, 10); h_NumMu->Sumw2(); 

   h_Num_PV_BeforePreSel = new TH1D(Form("h_Num_PV_BeforePreSel"), Form("Num. of Primary Vertex before Pre-Selection"), 100, 0.0, 100); h_Num_PV_BeforePreSel->Sumw2();// For Closure test // 
   h_Num_PV_AfterMetFilter = new TH1D(Form("h_Num_PV_AfterMetFilter"), Form("Num. of Primary Vertex After MetFilter"), 100, 0.0, 100); h_Num_PV_AfterMetFilter->Sumw2();// For Closure test // 
   h_Num_PV_AfterTrigger = new TH1D(Form("h_Num_PV_AfterTrigger"), Form("Num. of Primary Vertex After Trigger"), 100, 0.0, 100); h_Num_PV_AfterTrigger->Sumw2();// For Closure test // 

   h_Top1Mass     = new TH1D(Form("h_Top1Mass"   ), Form("Top1 Mass"   ), 1000, 0.0, 1000); h_Top1Mass->Sumw2(); 
   h_Top1pt       = new TH1D(Form("h_Top1pt"   ), Form("Top1 pt"   ), 1000, 0.0, 1000); h_Top1pt->Sumw2(); 
   h_Top1Rapidity = new TH1D(Form("h_Top1Rapidity"   ), Form("Top1 Rapidity"   ), 100, -5, 5); h_Top1Rapidity->Sumw2(); 
   h_Top1phi      = new TH1D(Form("h_Top1phi"   ), Form("Top1 phi"   ), 24, -1*pi, pi); h_Top1phi->Sumw2(); 
   h_Top1Energy   = new TH1D(Form("h_Top1Energy"   ), Form("Top1 Energy"   ), 1000, 0.0, 1000); h_Top1Energy->Sumw2(); 
   h_Top2Mass     = new TH1D(Form("h_Top2Mass" ), Form("Top2 Mass" ), 1000, 0.0, 1000); h_Top2Mass->Sumw2();
   h_Top2pt       = new TH1D(Form("h_Top2pt"   ), Form("Top2 pt"   ), 1000, 0.0, 1000); h_Top2pt->Sumw2(); 
   h_Top2Rapidity = new TH1D(Form("h_Top2Rapidity"   ), Form("Top2 Rapidity"   ), 100, -5, 5); h_Top2Rapidity->Sumw2(); 
   h_Top2phi      = new TH1D(Form("h_Top2phi"   ), Form("Top2 phi"   ), 24, -1*pi, pi); h_Top2phi->Sumw2(); 
   h_Top2Energy = new TH1D(Form("h_Top2Energy"   ), Form("Top2 Energy"   ), 1000, 0.0, 1000); h_Top2Energy->Sumw2(); 

   h_TopMass       = new TH1D(Form("h_TopMass"   ), Form("Top Mass"   ), 1000, 0.0, 1000); h_TopMass->Sumw2(); 
   h_Toppt         = new TH1D(Form("h_Toppt"   ), Form("Top pt"   ), 1000, 0.0, 1000); h_Toppt->Sumw2(); 
   h_TopRapidity   = new TH1D(Form("h_TopRapidity"   ), Form("Top Rapidity"   ), 100, -5, 5); h_TopRapidity->Sumw2(); 
   h_Topphi        = new TH1D(Form("h_Topphi"   ), Form("Top phi"   ), 24, -1*pi, pi); h_Topphi->Sumw2(); 
   h_TopEnergy     = new TH1D(Form("h_TopEnergy"   ), Form("Top Energy"   ), 1000, 0.0, 1000); h_TopEnergy->Sumw2(); 
   h_AnTopMass     = new TH1D(Form("h_AnTopMass" ), Form("AnTop Mass" ), 1000, 0.0, 1000); h_AnTopMass->Sumw2();
   h_AnToppt       = new TH1D(Form("h_AnToppt"   ), Form("AnTop pt"   ), 1000, 0.0, 1000); h_AnToppt->Sumw2(); 
   h_AnTopRapidity = new TH1D(Form("h_AnTopRapidity"   ), Form("AnTop Rapidity"   ), 100, -5, 5); h_AnTopRapidity->Sumw2(); 
   h_AnTopphi      = new TH1D(Form("h_AnTopphi"   ), Form("AnTop phi"   ), 24, -1*pi, pi); h_AnTopphi->Sumw2(); 
   h_AnTopEnergy   = new TH1D(Form("h_AnTopEnergy"   ), Form("AnTop Energy"   ), 1000, 0.0, 1000); h_AnTopEnergy->Sumw2(); 


   h_W1Mass     = new TH1D(Form("h_W1Mass"  ), Form("W1 Mass" ), 300, 0.0, 300); h_W1Mass->Sumw2();
   h_W2Mass     = new TH1D(Form("h_W2Mass"  ), Form("W2 Mass" ), 300, 0.0, 300); h_W2Mass->Sumw2();

   h_W1Mt     = new TH1D(Form("h_W1Mt"  ), Form("W1 Transverse Mass" ), 300, 0.0, 300); h_W1Mt->Sumw2();
   h_W2Mt     = new TH1D(Form("h_W2Mt"  ), Form("W2 Transverse Mass" ), 300, 0.0, 300); h_W2Mt->Sumw2();

   h_bJet1Energy = new TH1D(Form("h_bJet1Energy" ), Form("Leading bJet Energy"   ), 500, 0.0, 500); h_bJet1Energy->Sumw2(); 
   h_bJet2Energy = new TH1D(Form("h_bJet2Energy" ), Form("Second Leading bJet Energy" ), 500, 0.0, 500); h_bJet2Energy->Sumw2();

   h_bJetEnergy   = new TH1D(Form("h_bJetEnergy" ), Form("bJet Energy"   ), 1000, 0.0, 1000); h_bJetEnergy->Sumw2(); 
   h_AnbJetEnergy = new TH1D(Form("h_AnbJetEnergy" ), Form("b-barJet Energy" ), 1000, 0.0, 1000); h_AnbJetEnergy->Sumw2(); 
   h_bJetPt   = new TH1D(Form("h_bJetPt" ), Form("bJet Pt"   ), 1000, 0.0, 1000); h_bJetPt->Sumw2(); 
   h_AnbJetPt = new TH1D(Form("h_AnbJetPt" ), Form("b-barJet Pt" ), 1000, 0.0, 1000); h_AnbJetPt->Sumw2();

   h_Lep1Energy = new TH1D(Form("h_Lep1Energy" ), Form("Leading Lepton Energy"   ), 400, 0.0, 400); h_Lep1Energy->Sumw2(); 
   h_Lep2Energy = new TH1D(Form("h_Lep2Energy" ), Form("Second Leading Lepton Energy" ), 400, 0.0, 400); h_Lep2Energy->Sumw2();

   h_LepEnergy   = new TH1D(Form("h_LepEnergy" ), Form("Lepton Energy"   ), 400, 0.0, 400); h_LepEnergy->Sumw2(); 
   h_AnLepEnergy = new TH1D(Form("h_AnLepEnergy" ), Form("Anti-Lepton Energy" ), 400, 0.0, 400); h_AnLepEnergy->Sumw2();

   h_Nu1Energy = new TH1D(Form("h_Nu1Energy" ), Form("Leading Nuetrino Energy"   ), 400, 0.0, 400); h_Nu1Energy->Sumw2(); 
   h_Nu2Energy = new TH1D(Form("h_Nu2Energy" ), Form("Second Leading Nuetrino Energy" ), 400, 0.0, 400); h_Nu2Energy->Sumw2(); 

   h_NuEnergy   = new TH1D(Form("h_NuEnergy" ), Form("Nuetrino Energy"   ), 400, 0.0, 400); h_NuEnergy->Sumw2(); 
   h_AnNuEnergy = new TH1D(Form("h_AnNuEnergy" ), Form("anti-Nuetrino Energy" ), 400, 0.0, 400); h_AnNuEnergy->Sumw2();

   h_GenNuEnergy   = new TH1D(Form("h_GenNuEnergy" ), Form("Nuetrino Energy At Genrator Level"   ), 2000, 0.0, 2000); h_GenNuEnergy->Sumw2();
   h_GenAnNuEnergy = new TH1D(Form("h_GenAnNuEnergy" ), Form("anti-Nuetrino Energy At Genrator Level" ), 2000, 0.0, 2000); h_GenAnNuEnergy->Sumw2();

   h_CPO3_reco         = new TH1D(Form("h_CPO3_reco"   ), Form("CPO3_reco"   ), 200, -10, 10); h_CPO3_reco->Sumw2();
   h_CPO3_reco_JPRUp   = new TH1D(Form("h_CPO3_reco_JPRUp"     ), Form("CPO3_reco_JPRUp"     ), 200, -10, 10); h_CPO3_reco_JPRUp->Sumw2();
   h_CPO3_reco_JPRDown = new TH1D(Form("h_CPO3_reco_JPRDown"   ), Form("CPO3_reco_JPRDown"   ), 200, -10, 10); h_CPO3_reco_JPRDown->Sumw2();
   h_CPO3_reco_TopRes  = new TH1D(Form("h_CPO3_reco_TopRes"   ), Form("CPO3_reco in the Top Mass Window"   ), 200, -10, 10); h_CPO3_reco_TopRes->Sumw2();
   h_CPOb_reco         = new TH1D(Form("h_CPOb_reco"   ), Form("CPOb_reco"   ), 200, -10, 10); h_CPOb_reco->Sumw2();
   h_CPO5_reco         = new TH1D(Form("h_CPO5_reco"   ), Form("CPO5_reco"   ), 200, -10, 10); h_CPO5_reco->Sumw2();

   h_CPO3_evtweight  = new TH1D(Form("h_CPO3_evtweight"   ), Form("CPO3_evtweight"   ), 10, 0, 10); h_CPO3_evtweight->Sumw2();
   h_CPOb_evtweight  = new TH1D(Form("h_CPOb_evtweight"   ), Form("CPOb_evtweight"   ), 10, 0, 10); h_CPOb_evtweight->Sumw2();
   h_CPO5_evtweight  = new TH1D(Form("h_CPO5_evtweight"   ), Form("CPO5_evtweight"   ), 10, 0, 10); h_CPO5_evtweight->Sumw2();

   h_BjorkenX1 = new TH1D(Form("h_BjorkenX1" ), Form("Bjorken X1" ), 120, -0.1  , 1.1); h_BjorkenX1->Sumw2();  
   h_BjorkenX2 = new TH1D(Form("h_BjorkenX2" ), Form("Bjorken X2" ), 120, -0.1, 1.1); h_BjorkenX2->Sumw2(); 
   h_BjorkenX3 = new TH1D(Form("h_BjorkenX3" ), Form("Bjorken X3" ), 2000, 0.0  , 2000); h_BjorkenX3->Sumw2();

   /////////////////////////////
   /// To Check Evetn Weight ///
   /////////////////////////////

   h_LepbJetMass          = new TH1D(Form("h_LepbJetMass"), Form("Lep + bJet Invariant Mass"), 1000, 0, 1000); h_LepbJetMass->Sumw2();
   h_AnLepbJetMass        = new TH1D(Form("h_AnLepbJetMass"), Form("AnLep + bJet Invariant Mass"), 1000, 0, 1000); h_AnLepbJetMass->Sumw2();
   h2_TopMassVsLepBMass    = new TH2F(Form("TopMassVsLepBMass"), Form("Top Mass Vs (Lep+B) Mass"), 1000, 0, 1000, 1000,0,1000); h2_TopMassVsLepBMass->Sumw2();
   h2_AnTopMassVsLepBMass  = new TH2F(Form("AnTopMassVsLepBMass"), Form("Anti-Top Mass Vs (Lep+B) Mass"), 1000, 0, 1000, 1000,0,1000); h2_AnTopMassVsLepBMass->Sumw2();
   h2_AnLepBMassVsLepBMass     = new TH2F(Form("AnLepBMassVsLepBMass"), Form("(AntiLep+B) Mass Vs (Lep+B) Mass"), 1000, 0, 1000, 1000,0,1000); h2_AnLepBMassVsLepBMass->Sumw2();

   h_Diff_BBar_Pt     = new TH1D(Form("h_Diff_BBar_Pt"    ), Form("PT Difference between B and B bar "     ), 500, 0, 500); h_Diff_BBar_Pt->Sumw2();
   h_Diff_BBar_Energy = new TH1D(Form("h_Diff_BBar_Energy"), Form("Energy Difference between B and B bar"  ), 500, 0, 500); h_Diff_BBar_Energy->Sumw2();
   h_Diff_BBar_P      = new TH1D(Form("h_Diff_BBar_P"     ), Form("Momentum Difference between B and B bar"), 500, 0, 500); h_Diff_BBar_P->Sumw2();

   h_Diff_LepAnLep_Pt     = new TH1D(Form("h_Diff_LepAnLep_Pt"    ), Form("PT Difference between Lep. And AntiLep."      ), 500, 0, 500); h_Diff_LepAnLep_Pt->Sumw2();
   h_Diff_LepAnLep_Energy = new TH1D(Form("h_Diff_LepAnLep_Energy"), Form("Energy Difference between Lep. And AntiLep."  ), 500, 0, 500); h_Diff_LepAnLep_Energy->Sumw2();
   h_Diff_LepAnLep_P      = new TH1D(Form("h_Diff_LepAnLep_P"     ), Form("Momentum Difference between Lep. And AntiLep."), 500, 0, 500); h_Diff_LepAnLep_P->Sumw2();

   h_LepAnLepEngCheckO3 = new TH1D(Form("h_LepAnLepEngCheckO3"), Form("Lep1 and Lep2 pt check wheather they over the 300GeV for O3"), 100, 0, 100); h_LepAnLepEngCheckO3->Sumw2();
   h_LepAnLepEngCheckOb = new TH1D(Form("h_LepAnLepEngCheckOb"), Form("Lep1 and Lep2 pt check wheather they over the 300GeV for Ob"), 100, 0, 100); h_LepAnLepEngCheckOb->Sumw2();
   h_LepAnLepEngCheckO5 = new TH1D(Form("h_LepAnLepEngCheckO5"), Form("Lep1 and Lep2 pt check wheather they over the 300GeV for O5"), 100, 0, 100); h_LepAnLepEngCheckO5->Sumw2();

   h_PileUpCheck  = new TH1D(Form("h_PileUpCheck"),  Form("NewPileUpMethod - OriginMethod"), 4000, -2, 2); h_PileUpCheck->Sumw2();
   //////////////////////////////////
   /// For CP-Violation Variables ///
   //////////////////////////////////
   for (int i =0; i < 13; ++i)
   {
      h_Reco_CPO_[i] = new TH1D(Form("h_Reco_CPO%d",i+1 ), Form("CPO%d",i+1   ), 200, -10, 10); h_Reco_CPO_[i]->Sumw2();
      h_Reco_CPO_ReRange_[i] = new TH1D(Form("h_Reco_CPO%d_ReRange",i+1 ), Form("CPO%d",i+1   ), 40, -2, 2); h_Reco_CPO_ReRange_[i]->Sumw2();
      h_LepAnLepEngCheck_[i] = new TH1D(Form("h_LepAnLepEngCheck_%d",i+1), Form("Lep1 and Lep2 pt check wheather they over the 300GeV for CPO%d",i+1), 100, 0, 100); h_LepAnLepEngCheck_[i]->Sumw2();
   }
   ////////////////////////////////
   /// To Check O3 dilution ... ///
   ////////////////////////////////
   h_CPO3_Plus              = new TH1D(Form("h_CPO3_Plus"),              Form("CPO3 Positive"), 200, -10, 10); h_CPO3_Plus->Sumw2(); 
   h_CPO3_Minus             = new TH1D(Form("h_CPO3_Minus"),             Form("CPO3 Negative"), 200, -10, 10); h_CPO3_Minus->Sumw2(); 
   h_CPO3_JPRUp_Plus_Plus   = new TH1D(Form("h_CPO3_JPRUp_Plus_Plus"),   Form("CPO3_JPRUp Plus -> Plus"),   200, -10, 10); h_CPO3_JPRUp_Plus_Plus->Sumw2();                   
   h_CPO3_JPRUp_Plus_Minus  = new TH1D(Form("h_CPO3_JPRUp_Plus_Minus"),  Form("CPO3_JPRUp Plus -> Minus"),  200, -10, 10); h_CPO3_JPRUp_Plus_Minus->Sumw2();                   
   h_CPO3_JPRUp_Minus_Minus = new TH1D(Form("h_CPO3_JPRUp_Minus_Minus"), Form("CPO3_JPRUp Minus -> Minus"), 200, -10, 10); h_CPO3_JPRUp_Minus_Minus->Sumw2();                   
   h_CPO3_JPRUp_Minus_Plus  = new TH1D(Form("h_CPO3_JPRUp_Minus_Plus"),  Form("CPO3_JPRUp Minus -> Minus"), 200, -10, 10); h_CPO3_JPRUp_Minus_Plus->Sumw2();                   

   h_CPO3_JPRDown_Plus_Plus   = new TH1D(Form("h_CPO3_JPRDown_Plus_Plus"),   Form("CPO3_JPRDown Plus-> Plus"),   200, -10, 10); h_CPO3_JPRDown_Plus_Plus->Sumw2();
   h_CPO3_JPRDown_Plus_Minus  = new TH1D(Form("h_CPO3_JPRDown_Plus_Minus"),  Form("CPO3_JPRDown Plus-> Minus"),  200, -10, 10); h_CPO3_JPRDown_Plus_Minus->Sumw2();
   h_CPO3_JPRDown_Minus_Minus = new TH1D(Form("h_CPO3_JPRDown_Minus_Minus"), Form("CPO3_JPRDown Minus-> Minus"), 200, -10, 10); h_CPO3_JPRDown_Minus_Minus->Sumw2();
   h_CPO3_JPRDown_Minus_Plus  = new TH1D(Form("h_CPO3_JPRDown_Minus_Plus"),  Form("CPO3_JPRDown Minus-> Minus"), 200, -10, 10); h_CPO3_JPRDown_Minus_Plus->Sumw2();


   h_CPO3reco_evtweight             = new TH1D(Form("h_CPO3reco_evtweight"),             Form("Event Weight for CPO3 at Generator Level"),   10, 0, 10); h_CPO3reco_evtweight->Sumw2();
   h_CPO3reco_evtweight_EtaVariUp   = new TH1D(Form("h_CPO3reco_evtweight_EtaVariUp"),   Form("Event Weight for CPO3 at Generator Level"),   10, 0, 10); h_CPO3reco_evtweight_EtaVariUp->Sumw2();
   h_CPO3reco_evtweight_EtaVariDown = new TH1D(Form("h_CPO3reco_evtweight_EtaVariDown"), Form("Event Weight for CPO3 at Generator Level"),   10, 0, 10); h_CPO3reco_evtweight_EtaVariDown->Sumw2();
   h_CPO3reco_evtweight_PhiVariUp   = new TH1D(Form("h_CPO3reco_evtweight_PhiVariUp"),   Form("Event Weight for CPO3 at Generator Level"),   10, 0, 10); h_CPO3reco_evtweight_PhiVariUp->Sumw2();
   h_CPO3reco_evtweight_PhiVariDown = new TH1D(Form("h_CPO3reco_evtweight_PhiVariDown"), Form("Event Weight for CPO3 at Generator Level"),   10, 0, 10); h_CPO3reco_evtweight_PhiVariDown->Sumw2();

   h_CPO3_EtaVariUp   = new TH1D(Form("h_CPO3_EtaVariUp"),   Form("CPO3_EtaVariUp"),   200, -10, 10); h_CPO3_EtaVariUp->Sumw2(); 
   h_CPO3_EtaVariUp_Plus_Plus   = new TH1D(Form("h_CPO3_EtaVariUp_Plus_Plus"),   Form("CPO3_EtaVariUp Plus-> Plus"),   200, -10, 10); h_CPO3_EtaVariUp_Plus_Plus->Sumw2();
   h_CPO3_EtaVariUp_Plus_Minus  = new TH1D(Form("h_CPO3_EtaVariUp_Plus_Minus"),  Form("CPO3_EtaVariUp Plus-> Minus"),  200, -10, 10); h_CPO3_EtaVariUp_Plus_Minus->Sumw2();
   h_CPO3_EtaVariUp_Minus_Minus = new TH1D(Form("h_CPO3_EtaVariUp_Minus_Minus"), Form("CPO3_EtaVariUp Minus-> Minus"), 200, -10, 10); h_CPO3_EtaVariUp_Minus_Minus->Sumw2();
   h_CPO3_EtaVariUp_Minus_Plus  = new TH1D(Form("h_CPO3_EtaVariUp_Minus_Plus"),  Form("CPO3_EtaVariUp Minus-> Minus"), 200, -10, 10); h_CPO3_EtaVariUp_Minus_Plus->Sumw2();

   h_CPO3_EtaVariDown   = new TH1D(Form("h_CPO3_EtaVariDown"),   Form("CPO3_EtaVariDown"),   200, -10, 10); h_CPO3_EtaVariDown->Sumw2(); 
   h_CPO3_EtaVariDown_Plus_Plus   = new TH1D(Form("h_CPO3_EtaVariDown_Plus_Plus"),   Form("CPO3_EtaVariDown Plus-> Plus"),   200, -10, 10); h_CPO3_EtaVariDown_Plus_Plus->Sumw2();
   h_CPO3_EtaVariDown_Plus_Minus  = new TH1D(Form("h_CPO3_EtaVariDown_Plus_Minus"),  Form("CPO3_EtaVariDown Plus-> Minus"),  200, -10, 10); h_CPO3_EtaVariDown_Plus_Minus->Sumw2();
   h_CPO3_EtaVariDown_Minus_Minus = new TH1D(Form("h_CPO3_EtaVariDown_Minus_Minus"), Form("CPO3_EtaVariDown Minus-> Minus"), 200, -10, 10); h_CPO3_EtaVariDown_Minus_Minus->Sumw2();
   h_CPO3_EtaVariDown_Minus_Plus  = new TH1D(Form("h_CPO3_EtaVariDown_Minus_Plus"),  Form("CPO3_EtaVariDown Minus-> Minus"), 200, -10, 10); h_CPO3_EtaVariDown_Minus_Plus->Sumw2();

   h_CPO3_PhiVariUp   = new TH1D(Form("h_CPO3_PhiVariUp"),   Form("CPO3_PhiVariUp"),   200, -10, 10); h_CPO3_PhiVariUp->Sumw2(); 
   h_CPO3_PhiVariUp_Plus_Plus   = new TH1D(Form("h_CPO3_PhiVariUp_Plus_Plus"),   Form("CPO3_PhiVariUp Plus-> Plus"),   200, -10, 10); h_CPO3_PhiVariUp_Plus_Plus->Sumw2();
   h_CPO3_PhiVariUp_Plus_Minus  = new TH1D(Form("h_CPO3_PhiVariUp_Plus_Minus"),  Form("CPO3_PhiVariUp Plus-> Minus"),  200, -10, 10); h_CPO3_PhiVariUp_Plus_Minus->Sumw2();
   h_CPO3_PhiVariUp_Minus_Minus = new TH1D(Form("h_CPO3_PhiVariUp_Minus_Minus"), Form("CPO3_PhiVariUp Minus-> Minus"), 200, -10, 10); h_CPO3_PhiVariUp_Minus_Minus->Sumw2();
   h_CPO3_PhiVariUp_Minus_Plus  = new TH1D(Form("h_CPO3_PhiVariUp_Minus_Plus"),  Form("CPO3_PhiVariUp Minus-> Minus"), 200, -10, 10); h_CPO3_PhiVariUp_Minus_Plus->Sumw2();

   h_CPO3_PhiVariDown   = new TH1D(Form("h_CPO3_PhiVariDown"),   Form("CPO3_PhiVariDown"),   200, -10, 10); h_CPO3_PhiVariDown->Sumw2();
   h_CPO3_PhiVariDown_Plus_Plus   = new TH1D(Form("h_CPO3_PhiVariDown_Plus_Plus"),   Form("CPO3_PhiVariDown Plus-> Plus"),   200, -10, 10); h_CPO3_PhiVariDown_Plus_Plus->Sumw2();
   h_CPO3_PhiVariDown_Plus_Minus  = new TH1D(Form("h_CPO3_PhiVariDown_Plus_Minus"),  Form("CPO3_PhiVariDown Plus-> Minus"),  200, -10, 10); h_CPO3_PhiVariDown_Plus_Minus->Sumw2();
   h_CPO3_PhiVariDown_Minus_Minus = new TH1D(Form("h_CPO3_PhiVariDown_Minus_Minus"), Form("CPO3_PhiVariDown Minus-> Minus"), 200, -10, 10); h_CPO3_PhiVariDown_Minus_Minus->Sumw2();
   h_CPO3_PhiVariDown_Minus_Plus  = new TH1D(Form("h_CPO3_PhiVariDown_Minus_Plus"),  Form("CPO3_PhiVariDown Minus-> Minus"), 200, -10, 10); h_CPO3_PhiVariDown_Minus_Plus->Sumw2();

   h_bTagWeight   = new TH1D(Form("h_bTagWeight"),   Form("bTag Weight "),   200, -2, 2); h_bTagWeight->Sumw2();

}

void ssb_analysis::End()
{
   fout->Write();
   fout->Close();
}

void ssb_analysis::SetInputFileName( char *inname )
{
   FileName_ = inname;
   char unsco_ = '_';
   Size_t unscoIndex = FileName_.Last(unsco_);
   FileName_.Remove(unscoIndex, FileName_.Length());
   //cout << "FileName_ : " << FileName_ << endl;
   
}
void ssb_analysis::SetOutputFileName(char *outname)
{   
   outfile = outname;
}

void ssb_analysis::MCSF()
{
   if (FileName_.Contains("Data")||FileName_.Contains("Single")||FileName_.Contains("EG")){ mc_sf_ = 1.; return; }
   /// Open Xsec Tables ///
   FILE *xsecs_;
   char sampleName[1000];
   double xsec_ = -1.;
   double br_ = -1.; 
   int totalevt_  = -1.; 
   int positive_  = -1.; 
   int negative_  = -1.; 
   int posi_nega_ = -1.; 
   string xsec_dir= "./xsecAndsample/";
   string xsec_filePath = xsec_dir+ XsecTable_.Data();
   //cout << "xsec_filePath : " << xsec_filePath << endl;
   /// SampleName | TotalEvt | Positive+Negative | Xsection | Branching Fraction |
   xsecs_ = fopen(xsec_filePath.c_str(),"r");
   map<string, int> m_sam_totalevt;
   map<string, double> m_sam_xsec;
   map<string, double> m_sam_br;
   map<string, int> m_sam_positive;
   map<string, int> m_sam_negative;
   map<string, int> m_sam_posi_nega;
   if (xsecs_!=NULL) 
   { 
      //cout << "Load Xsection Table!" << endl;
      while (fscanf(xsecs_, "%s %d %d %d %d %lf %lf\n", sampleName, &totalevt_, &positive_, &negative_, &posi_nega_, &xsec_, &br_ ) != EOF)
      {
         /*cout 
         << "sampleName : " << sampleName << " totalevt_ : " << totalevt_
         << " positive_ " << positive_ << " negative_ : " << negative_ 
         << " posi_nega_ " << posi_nega_ << " xsec_ : " << xsec_ 
         << " br_ " << br_ 
         << endl;*/
         m_sam_totalevt[sampleName] = totalevt_;
         m_sam_positive[sampleName] = positive_;
         m_sam_negative[sampleName] = negative_;
         m_sam_posi_nega[sampleName] = posi_nega_;
         m_sam_xsec[sampleName] = xsec_;
         m_sam_br[sampleName] = br_;
      }
      fclose(xsecs_);
   }
   else {return;} 
   //cout << "Lumi : " << Lumi << endl;
   double lumi = Lumi/1000000;
   mc_sf_ = (m_sam_xsec[FileName_.Data()]*m_sam_br[FileName_.Data()]*lumi)/m_sam_posi_nega[FileName_.Data()];
   return;
}

// Apply MC SF To Event //
void ssb_analysis::MCSFApply()
{
   evt_weight_beforemcsf_ =1; // Initailize evt_weight_beforemcsf_ //
   evt_weight_beforemcsf_ = evt_weight_; // keep event weight // 
   if ( !TString(FileName_).Contains( "Data") ){ evt_weight_ = evt_weight_*mc_sf_; } // apply MC scale factor // 
   else {evt_weight_ = 1;}
}
void ssb_analysis::GenWeightApply()
{
   double genweight = 1.0;
   evt_weight_beforegenweight_ =1; // Initailize evt_weight_beforegenweight_ //
   evt_weight_beforegenweight_ = evt_weight_; // keep event weight //
   if ( !TString(FileName_).Contains( "Data") ){
      if (Gen_EventWeight > 0.0){genweight =1;}
      else {genweight =-1;} 
      evt_weight_ = evt_weight_*genweight;
   }
   else {evt_weight_ = 1;}
}
void ssb_analysis::PDFWeightApply()
{
   double pdfweight = 1.0;
   evt_weight_beforepdfweight_ = evt_weight_;
   if ( !TString(FileName_).Contains("Data"))
   { 
      if ( TString(FileName_).Contains("TTJets_") ){
         if ( PDFSys == 0 ) { pdfweight = 1.0; }// central //
         else if ( PDFSys > 0 )  
         {
             int index_pdf = PDFSys + 8;
             pdfweight = LHE_Weight->at(index_pdf)/LHE_Central;
             //cout << "LHE_Central : "<< LHE_Central <<" LHE_Id["<< index_pdf <<"] : " << LHE_Id->at(index_pdf) << endl;
         }
         else {cout << "Somethig Wrong in the PDF sys ... " << endl;}
      }
      else { pdfweight =1; }
      evt_weight_ = evt_weight_*pdfweight;
   }
   else { evt_weight_ = 1; }
   return;
}
void ssb_analysis::FactRenoApply()
{
   double factrenosys = 1.0;
   evt_weight_beforefactrenoweight_ =1; // Initailize evt_weight_beforefactrenoweight_ //
   evt_weight_beforefactrenoweight_ = evt_weight_; // keep event weight //
   if (LHE_Weight->size() == 0 && !TString(FileName_).Contains("Data") ) {cout << "There is no LHE_Weight!!!" <<endl;  return;}
   if (TString(FileName_).Contains("TTJets")) {
      if (FactRenoSys == 0) { factrenosys = 1; }
      else if (FactRenoSys == 1){ factrenosys = LHE_Weight->at(1)/LHE_Central; }
      else if (FactRenoSys == 2){ factrenosys = LHE_Weight->at(2)/LHE_Central; }
      else if (FactRenoSys == 3){ factrenosys = LHE_Weight->at(3)/LHE_Central; }
      else if (FactRenoSys == 4){ factrenosys = LHE_Weight->at(4)/LHE_Central; }
      else if (FactRenoSys == 5){ factrenosys = LHE_Weight->at(5)/LHE_Central; }
      else if (FactRenoSys == 6){ factrenosys = LHE_Weight->at(6)/LHE_Central; }
      else if (FactRenoSys == 7){ factrenosys = LHE_Weight->at(7)/LHE_Central; }
      else if (FactRenoSys == 8){ factrenosys = LHE_Weight->at(8)/LHE_Central; }
      else { cout << "Check out your FactRenoSys " << FactRenoSys << endl;  }
   }
   else { factrenosys =1; }
   evt_weight_ = evt_weight_*(factrenosys);
   return;
}
void ssb_analysis::FragmentApply()
{
   double fragmentsys = 1.0;
   evt_weight_beforefragmentweight_ =1; // Initailize evt_weight_beforefragmentweight_ //
   evt_weight_beforefragmentweight_ = evt_weight_; // keep event weight //
   if (TString(FileName_).Contains("TTJets")) {
      if (FragmentSys == "nominal") { fragmentsys = 1.0; }
      else if (FragmentSys == "up"){ fragmentsys = Frag_Up_Weight; }
      else if (FragmentSys == "central"){ fragmentsys = Frag_Cen_Weight; }
      else if (FragmentSys == "down"){ fragmentsys = Frag_Down_Weight; }
      else if (FragmentSys == "peterson"){ fragmentsys = Frag_Peterson_Weight; }
      else { cout << "Check out your FragmentSys " << FragmentSys << endl;  }
   }
   else { fragmentsys =1; }
   evt_weight_ = evt_weight_*(fragmentsys);
   return;
}
/// Systematic for Decay tables
void ssb_analysis::DecayTableApply()
{
   double dectabsys = 1.0;
   evt_weight_beforedectabweight_ =1; // Initailize evt_weight_beforedectabweight_ //
   evt_weight_beforedectabweight_ = evt_weight_; // keep event weight //
   if (TString(FileName_).Contains("TTJets")) {
      if (DecayTableSys == "central") { dectabsys = 1.0; }
      else if (DecayTableSys == "up"){ dectabsys = Semilep_BrUp_Weight; }
      else if (DecayTableSys == "down"){ dectabsys = Semilep_BrDown_Weight; }
      else { cout << "Check out your DecayTableSys " << DecayTableSys << endl;  }
   }
   else { dectabsys =1; }
   evt_weight_ = evt_weight_*(dectabsys);
   return;   
}

// Apply Trigger SF To Event //
void ssb_analysis::TriggerSFApply()
{
   double triggersf_ =1;
   evt_weight_beforeTrigger_ =1; // Initailize evt_weight_beforeTrigger_ //
   evt_weight_beforeTrigger_ = evt_weight_; // keep event weight // 
   if ( !TString(FileName_).Contains( "Data") )
    {  //evt_weight_ = evt_weight_*triggersf_;
/*      if ( TString(Decaymode).Contains( "dielec" ) ) { triggersf_ = 0.996; }
      else if (TString(Decaymode).Contains( "muel" )  )  { triggersf_ = 0.990;}
      else if (TString(Decaymode).Contains( "dimu" )  )  { triggersf_ = 0.996;}*/
      if ( TString(Decaymode).Contains( "dielec" ) ) { triggersf_ =SSBEffcal->TrigDiElec_Eff(Lep1,Lep2,TrigSFSys); }
      else if (TString(Decaymode).Contains( "muel" )  )  { //triggersf_ = SSBEffcal->TrigMuElec_Eff(TElectron,TMuon); 
         triggersf_ = SSBEffcal->TrigMuElec_Eff(Lep1,Lep2,TrigSFSys);
      }
      else if (TString(Decaymode).Contains( "dimu" )  )  { triggersf_ = SSBEffcal->TrigDiMuon_Eff(Lep1,Lep2,TrigSFSys);}
      else { cout << "CHECK OUT TRIGGER EFF"<< endl;}
      evt_weight_ = evt_weight_*triggersf_;
   }// apply Trigger scale factor //
   else { evt_weight_ = evt_weight_; }
   return;
}
void ssb_analysis::PileUpReWeightApply()
{
   evt_weight_beforePileup_ = 1;
   evt_weight_beforePileup_ = evt_weight_; // keep event weight // 
   double puweight_ = 1.;
   double pu_weight_central = puweight->weight( PileUp_Count_Intime );
   cout << "pu_weight_central : " << pu_weight_central << endl;
   if ( !TString(FileName_).Contains( "Data") )
   {
      puweight_ = puweight->weight( PileUp_Count_Intime );
//      if (TString(PileUpSys).Contains("central") ) { puweight_   = pu_weight_central;}
//      else {  
      cout << "PileUp sys Error ... Defalut is Weight_PileUp ... : " << PileUpSys << endl;
   }
   else {evt_weight_ = 1;}
   evt_weight_ = evt_weight_*puweight_; // apply PileUpReweight //
   return;
}
void ssb_analysis::L1PreFireApply()
{
   evt_weight_beforeL1PreFire_ = 1;
   evt_weight_beforeL1PreFire_ = evt_weight_; // keep event weight // 
   double l1prefire_ = 1.;
   double l1prefire_central = L1_PreFire_Central;
   double l1prefire_up      = L1_PreFire_Up;
   double l1prefire_dn      = L1_PreFire_Down;

   if ( !TString(FileName_).Contains( "Data") )
   {
      if (TString(L1PreFireSys).Contains("central") ) { l1prefire_   = l1prefire_central;}
      else if (TString(L1PreFireSys).Contains("up") ) { l1prefire_   = l1prefire_up; }
      else if (TString(L1PreFireSys).Contains("down") ) { l1prefire_ = l1prefire_dn; }
      else if (TString(L1PreFireSys).Contains("none") ) { l1prefire_ = 1.0; }
      else {  
      cout << "PileUp sys Error ... Defalut is Weight_PileUp ... : " << L1PreFireSys << endl;
   }
   evt_weight_ = evt_weight_*l1prefire_; } // apply PileUpReweight //
   else {evt_weight_ = 1;}
}

void ssb_analysis::NumPVCount()
{
   num_pv = 0;
   //cout << " start NumPVCount "<< endl;
   for (int i =0; i < Filter_PV->size(); ++i)
   {
      //cout <<  " Filter_PV->at(i) : " <<  Filter_PV->at(i) << endl;
      if ( Filter_PV->at(i) == true ) {num_pv++;}
   }
   //cout << "counted Num PV " << endl;
}

bool ssb_analysis::EveRun()
{
   bool everun = false;
   if (TString(FileName_).Contains( "Data" ) ) // Only Data
   {
      if ( TString(FileName_).Contains( "Jul" ) ) //Jul-Data
      {
         if ( TString(Decaymode).Contains( "dimuon" ) || TString(Decaymode).Contains( "dielec" ) || TString(Decaymode).Contains( "muel" )  )
         {
            if (Info_RunNumber >= 246908 && Info_RunNumber < 251604 ) { everun = true; }
         }
         else { cout << "YOU DON'T USE THIS FUNCTION , CHECH DECAY MODE" << endl;}
      }
      else 
      {
         if ( TString(Decaymode).Contains( "dimuon" ) || TString(Decaymode).Contains( "dielec" ) || TString(Decaymode).Contains( "muel" ) )
         {
            if (Info_RunNumber >= 251604 && Info_RunNumber <= 251883 ) { everun = true; }
         }
         else { cout << "YOU DON'T USE THIS FUNCTION , CHECH DECAY MODE" << endl;}
      }//prompt
   }
   else { everun = true;} //MC
   return everun;
}
bool ssb_analysis::ChannelIndex()
{
   bool chan_cut = true;
   int decay_idx = 0;
   if ( TString(Decaymode).Contains( "dielec" ) ) decay_idx = 22;
   if ( TString(Decaymode).Contains( "muel" ) )   decay_idx = 24;
   if ( TString(Decaymode).Contains( "dimuon" ) ) decay_idx = 26;

   if ( Channel_Idx != decay_idx && ( TString(FileName_).Contains("other") == false ) && TString(FileName_).Contains("TTJets") ) {chan_cut=false;} // only TTbar -> Di-Lepton Signal.
   if ( Channel_Idx == decay_idx && TString(FileName_).Contains("other") ) {chan_cut = false;} // only TTbar -> Di-Lepton Signal.
   return chan_cut;
}
bool ssb_analysis::METFilterAPP()
{
   bool metfilt_ = true;
   for ( int i = 0; i < METFilter_Name->size(); ++i )
   { 
      TString METFiltName = METFilter_Name->at(i);
      //cout << "METFiltName : " << METFilter_Name->at(i) << endl;
      if ( !(TString(FileName_).Contains("Data")) && 
            TString(METFiltName).Contains("Flag_eeBadScFilter")  ) {continue;}
      for( int j =0; j < v_METFilterName.size(); ++j )
      {
         if ( TString( METFiltName ).Contains(v_METFilterName.at(j)) )
         {
            if ( !(METFilter_isPass->at(i)) ){ metfilt_ = false; }
         } 
      }
   }
   for ( int i = 0; i < METFilterAdd_Name->size(); ++i )
   {
      TString METFiltName = METFilterAdd_Name->at(i);
      if ( !(TString(FileName_).Contains("Data")) && 
            TString(METFiltName).Contains("Flag_eeBadScFilter")  ) {continue;}
      for( int j =0; j < v_METFilterName.size(); ++j )
      {
         if ( TString( METFiltName ).Contains(v_METFilterName.at(j)) )
         {
            if ( !( METFilterAdd_isPass->at(i) ) ){ metfilt_ = false; }
         } 
      }
   }
   return metfilt_;
   //return true;
}

// Trigger Requirement Function
bool ssb_analysis::SelTrigger(vector<string> v_sel)
{
   TString trgName = "";
   int ptrigindex;
   bool passtrig_ = false;
   ptrigindex =0;

   for (int i =0; i < Trigger_Name->size(); i++)
   {
      //cout << "Ntuple Triggers : " <<  Trigger_Name->at(i) << endl; 
      for (int j = 0; j < v_sel.size(); j++)
      {
         trgName = v_sel[j];
   
         if ( TString( Trigger_Name->at(i) ).Contains( v_sel.at(j) ) )//TString clone 
         {
            //cout << "v_sel.at(j) ?" << v_sel.at(j) << endl;
            if ( ( Trigger_isPass->at(i)  ) && 
                !( Trigger_isError->at(i) ) && 
                 ( Trigger_isRun->at(i) )      ) 
            { 
               ptrigindex = ptrigindex+1; 
               //cout << "Trigger_isPass->at(i)   : " << Trigger_isPass->at(i)   << endl;
            }
         }
      }
   }
   //cout << "ptrigindex : "  << ptrigindex << endl;
   if ( ptrigindex > 0 ) { passtrig_ = true; }
   return passtrig_; 

}
bool ssb_analysis::Trigger()
{
   /// Variable for Trigger Function
   bool singleTrig_ = false;
   bool doubleTrig_ = false;
   bool ispassselTrig_ = false;
   bool ispassvetoTrig_ = false;
   bool trigpass = false;
   vector<string> seltrigName;
   vector<string> vetotrigName;
 
   if (!TString(FileName_).Contains( "Data" ) ) { 
       //trigpass = true; return trigpass;
       ispassvetoTrig_ = false; 
       seltrigName = trigName;
       trigpass = SelTrigger(seltrigName);
   }
   else {
      // Set veto trigger & selected trigger //
      // Channel Index //
      if ( RunPeriod.Contains("2018") )/// Only 2018, SingleEG and Double EG combined 
      {
         if (TString(Decaymode).Contains("dimuon")){ // Dimuon // 
            if ( TString(FileName_).Contains( "Single") ) {
               seltrigName = SLtrigName;
               vetotrigName = DLtrigName;
               ispassselTrig_ = SelTrigger(seltrigName);
               ispassvetoTrig_ = SelTrigger(vetotrigName);
            }
            else if( TString(FileName_).Contains( "Double") ) {
               seltrigName = DLtrigName;
               vetotrigName = SLtrigName;
               ispassselTrig_ = SelTrigger(seltrigName);
               ispassvetoTrig_ = SelTrigger(vetotrigName);
            }
            else { cout << "Check out FileName_ in Trigger ()" << endl;}
         }
         else if(TString(Decaymode).Contains("muel")) {
            if ( TString(FileName_).Contains( "MuonEG") ) {
               seltrigName = DLtrigName;
               vetotrigName = SLtrigName;
               ispassselTrig_ = SelTrigger(seltrigName);
               ispassvetoTrig_ = SelTrigger(vetotrigName);
            }
            else if( TString(FileName_).Contains( "EGamma") ) {
               seltrigName = SLtrigName;
               vetotrigName = DLtrigName;
               ispassselTrig_ = SelTrigger(seltrigName);
               ispassvetoTrig_ = SelTrigger(vetotrigName);
            }
            else { cout << "Check out FileName_ in Trigger ()" << endl;}
         }
         else if(TString(Decaymode).Contains("dielec")) { 
            //seltrigName = DLtrigName;
            //vetotrigName = SLtrigName;
            if( TString(FileName_).Contains( "EGamma") ) {
               singleTrig_ = false;
               doubleTrig_ = false;
               singleTrig_ = SelTrigger(SLtrigName);
               doubleTrig_ = SelTrigger(DLtrigName);
               if (singleTrig_&&(!doubleTrig_)) {trigpass = true;}
               else if ((!singleTrig_)&&(doubleTrig_)) {trigpass = true;}
               else {trigpass = false;}
               ispassvetoTrig_ == false;
            }
         }
         else { cout << "Check out Decaymode in 2018-Trigger()" << endl; }
      }
      else {/// 2016 (AVP, NonAPV) && 2017 
         if ( TString(FileName_).Contains( "Single") ) {
            seltrigName = SLtrigName;
            vetotrigName = DLtrigName;
            cout << "selected!! " << endl;
            ispassselTrig_ = SelTrigger(seltrigName);
            cout << "ispassselTrig_ : " << ispassselTrig_ << endl;
            cout << "veto !! " << endl;
            ispassvetoTrig_ = SelTrigger(vetotrigName);
            cout << "ispassvetoTrig_ : " << ispassvetoTrig_ << endl;
         }
         else if( TString(FileName_).Contains( "Double") || TString(FileName_).Contains( "MuonEG")) {
            seltrigName = DLtrigName;
            vetotrigName = SLtrigName;
            ispassselTrig_ = SelTrigger(seltrigName);
            ispassvetoTrig_ = SelTrigger(vetotrigName);
         }
         else { cout << "Check out FileName_ in Trigger ()" << endl;}

      }

      if (ispassvetoTrig_ == true) {trigpass = false;}
   }
   cout << "trigger : " << trigpass << endl; 
   return trigpass;
}
// Function of Muon Rocheser Correction //
TLorentzVector* ssb_analysis::ApplyRocCor(TLorentzVector* tmp,int index_)
{
   TLorentzVector* applied = tmp;// new TLorentzVector();
   return applied;
}
void ssb_analysis::CorrectedMuonCollection()
{
   //cout << "Star CorrectedMuonCollec"<< endl;
   v_CorMuons.clear();
   TLorentzVector* tmp_;
   for (int i =0 ; i < Muon->GetEntriesFast(); ++i )
   {
      tmp_ = (TLorentzVector*)Muon->At(i);
      tmp_ = ApplyRocCor(tmp_,i);
      v_CorMuons.push_back(tmp_);
   }
}
void ssb_analysis::MakeElecCollection()
{
   //cout << "Star MakeElecCollection"<< endl;
   v_Elec.clear();
   TLorentzVector* tmp_;
   for (int i =0 ; i < Elec->GetEntriesFast(); ++i )
   {
      tmp_ = (TLorentzVector*)Elec->At(i);
      tmp_ = ApplyElecSCSM(tmp_,i);
      v_Elec.push_back(tmp_);
   }

}
TLorentzVector* ssb_analysis::ApplyElecSCSM( TLorentzVector* tmp,int index_)
{
   if (FileName_.Contains("Data")) {return tmp;}
   if (EleScSmSys == "Central")    {return tmp;}
   TLorentzVector* systl_ = new TLorentzVector();
   double scsmfac_ = 1.0;
/*   if (EleScSmSys == "UpUp"){scsmfac_ = Elec_ScSmUpUp->at(index_);}
   else if (EleScSmSys == "UpDown"){scsmfac_ = Elec_ScSmUpDown->at(index_);}
   else if (EleScSmSys == "DownUp"){scsmfac_ = Elec_ScSmDownUp->at(index_);}
   else if (EleScSmSys == "DownDown"){scsmfac_ = Elec_ScSmDownDown->at(index_);}
   else {cout << "Check Out Your EleScSmSys Option !!! " << EleScSmSys << endl;}
   //cout << "scsmfac_ : " << scsmfac_ << endl;
   systl_->SetPtEtaPhiM(scsmfac_*tmp->Pt(),tmp->Eta(),tmp->Phi(),tmp->M());
   return systl_;*/
   return tmp;
}
// Function for Lepton Selection
void ssb_analysis::LeptonSelector()
{
   v_lepton_idx.clear();
   v_muon_idx.clear();
   v_electron_idx.clear();

   // dummy vector //
   v_lep_idx_temp.clear();
   v_muon_idx_temp.clear();
   v_electron_idx_temp.clear();

   /// Muon Correction ///
   CorrectedMuonCollection();
   /// Electron Class for Central & Syst.
   MakeElecCollection();
   // Case of Dimuon //
   if (TString(Decaymode).Contains("dimuon")) 
   { 
      Int_t nmu = Muon->GetEntriesFast();
      for (int i =0; i < nmu; ++i )
      {
         Mu_lep_sel = v_CorMuons.at(i);
         if ( ( v_lepton_iso->at(i) > lepisocut )   ||
              (*Mu_lep_sel).Pt()    < lep_pt        || 
              fabs( (*Mu_lep_sel).Eta() ) > lep_eta || 
              ( v_lepton_Id->at(i) == false )          ){continue;}         
         v_lep_idx_temp.push_back( i );

         if (v_lep_idx_temp.size()==1 && v_lepton_idx.size() == 0 ){  v_lepton_idx.push_back(i);}
         if (v_lep_idx_temp.size()>0)
         {
            if ( Muon_Charge->at( v_lep_idx_temp.at(0) ) != Muon_Charge->at( i ) ) { v_lepton_idx.push_back(i); }
         }
      }
   }
   // Case of Dielectron //
   else if (TString(Decaymode).Contains("dielec")) 
   { 
      Int_t nel = Elec->GetEntriesFast();

      for (int i =0; i < nel; ++i )
      {
         Ele_lep_sel = v_Elec.at(i);
         if ( (*Ele_lep_sel).Pt() != (*Ele_lep_sel).Pt() ){continue;}
         if ( (*Ele_lep_sel).Pt() < lep_pt )          {continue;}
         if ( fabs( (*Ele_lep_sel).Eta() ) > lep_eta ){continue;}   

         // Requirement of Electron Isolation with Electron Eta 
         //if (Elec_SCB_Medium->at(i) == false) {continue;}    
         //if (v_electron_Id->at(i) == false) {continue;}    
         if (v_lepton_Id->at(i) == false) {continue;}    
         if ( fabs( Elec_Supercluster_Eta->at(i) ) > 1.4442 && fabs( Elec_Supercluster_Eta->at(i) ) < 1.566 ) {continue;} 
         if ( Elec_ChargeId_GsfPx->at(i) == false ){continue;}//
         v_lep_idx_temp.push_back( i );
         if ( v_lep_idx_temp.size()==1 && v_lepton_idx.size() == 0 ){  v_lepton_idx.push_back(i); }
         if ( v_lep_idx_temp.size() > 0 )
         {
            if ( Elec_Charge->at( v_lep_idx_temp.at(0) ) != Elec_Charge->at( i ) ) { v_lepton_idx.push_back(i); }
         }

      }
   }
   // Case of MuEl //
   else if (TString(Decaymode).Contains("muel"))
   {

      Int_t nmu = Muon->GetEntriesFast();
      for (int i =0; i < nmu; ++i )
      {
         Mu_lep_sel = v_CorMuons.at(i);
         if ( ( v_mlep_iso->at(i) > muonisocut ) || (*Mu_lep_sel).Pt() < muon_pt || fabs( (*Mu_lep_sel).Eta() ) > muon_eta || v_muon_Id->at(i) == false  ){continue;}         
         //if ( ( v_mlep_iso->at(i) > muonisocut ) || (*Mu_lep_sel).Pt() < 35 || fabs( (*Mu_lep_sel).Eta() ) > muon_eta || v_muon_Id->at(i) == false  ){continue;}         
         v_muon_idx.push_back( i );
      }

      Int_t nel = Elec->GetEntriesFast();
      for (int i =0; i < nel; ++i )
      {
         Ele_lep_sel = v_Elec.at(i);

         ////
         if ( (*Ele_lep_sel).Pt() < elec_pt ){continue;}
         //if ( (*Ele_lep_sel).Pt() < 35 ){continue;}

         if ( fabs( (*Ele_lep_sel).Eta() ) > elec_eta ){continue;}

         //if (Elec_SCB_Medium->at(i) == false) {continue;}    
         if (v_electron_Id->at(i) == false) {continue;}    
         if ( Elec_Conversion->at(i) == false ){continue;}
         if ( fabs( Elec_Supercluster_Eta->at(i) ) > 1.4442 && fabs( Elec_Supercluster_Eta->at(i) ) < 1.566 ) {continue;}
         if ( Elec_ChargeId_GsfPx->at(i) == false ){continue;}        
         v_electron_idx_temp.push_back( i );

      }

      if (v_muon_idx.size() > 0 && v_electron_idx_temp.size() > 0 )
      {  
         for ( int i = 0; i < v_electron_idx_temp.size(); ++i )
         {
            if ( Muon_Charge->at( v_muon_idx.at(0) ) != Elec_Charge->at( v_electron_idx_temp.at(i) ) ) { v_electron_idx.push_back(v_electron_idx_temp.at(i) ); }
         }
      }
      
      if ( v_electron_idx.size() > 0 && v_muon_idx.size() > 0 ) 
      { 
         if ( Muon_Charge->at( v_muon_idx.at(0) ) == Elec_Charge->at( v_electron_idx.at(0) ) ) { cout << "MuEl seclection wrong !! " << endl; } 
      }

   }
   else { cout << "Lepton Selection error" << endl;}
 

   ///////////////////////////////////////////
   /// Electron Selection For Jet Cleaning /// 
   ///////////////////////////////////////////
   v_ele_idx_cand.clear();
   v_mu_idx_cand.clear();
   Int_t inele = Elec->GetEntriesFast();
   Int_t nmu   = Muon->GetEntriesFast();

   for (int ine = 0; ine <inele; ine++ )
   {  
      Ele_jet_clean = v_Elec.at(ine); 

      if (
           Elec_SCB_Veto->at(ine) == true &&
           Ele_jet_clean->Pt() > el_pt_jetcl &&
           fabs( Ele_jet_clean->Eta() ) < el_eta_jetcl &&
           Elec_Conversion->at(ine) == true
         )
      {
//         v_ele_idx_cand.push_back( ine );
         if ( fabs( Elec_Supercluster_Eta->at(ine) ) <= 1.479 )
         {
            /*if (
                 v_el_iso_jcl->at( ine ) < 0.175 &&
                 fabs( Elec_Track_GsfdXY->at(ine) ) < 0.05 &&
                 fabs( Elec_Track_GsfdZ->at(ine) ) < 0.10 
               )
            { v_ele_idx_cand.push_back( ine ); }*/

            v_ele_idx_cand.push_back( ine ); 
         } 
         else 
         {
            /*if (
                 v_el_iso_jcl->at( ine ) < 0.159 &&           
                 fabs( Elec_Track_GsfdXY->at(ine) ) < 0.10 &&
                 fabs( Elec_Track_GsfdXY->at(ine) ) < 0.20 
               )
            { v_ele_idx_cand.push_back( ine ); }*/
            v_ele_idx_cand.push_back( ine ); 
        
         }
      }

   }
     
   // Muon selectron for JetCleaning
   for ( int inmu =0; inmu < nmu; ++inmu )
   {
      Mu_jet_clean = v_CorMuons.at(inmu);
      if ( 
           (*Mu_jet_clean).Pt() > mu_pt_jetcl         && //Pt cut 
           fabs((*Mu_jet_clean).Eta()) < mu_eta_jetcl && //Eta cut
           v_mu_iso_jcl->at(inmu) < mu_iso_jetcl      && //Isolation cut...
           v_mu_Id_jcl->at(inmu) == true                 // Id cut ...
         )
      { v_mu_idx_cand.push_back( inmu ); }
   }

}
bool ssb_analysis::NumIsoLeptons()//YOU SHOULD REQUIRE THIS FUNCTION AFTER LEPTONSELETOR //
{
   bool numLeptons = true;
   if ( TString(Decaymode).Contains( "dimuon" ) ||
        TString(Decaymode).Contains( "dielec" )    ){  if ( v_lepton_idx.size() <= 1 ) {numLeptons=false;} } //Requiring  Two isolated Lepton .

   else if ( TString(Decaymode).Contains( "muel" ) ){ if (v_muon_idx.size() < 1 || v_electron_idx.size() < 1 ) {numLeptons=false; } else {}}
   else if ( TString(Decaymode).Contains( "muonJet" ) ){if( v_lepton_idx.size() != 1){numLeptons=false;} }
   else { cout << "?? something wrong " << endl; }
   return numLeptons;
}
bool ssb_analysis::LeptonsPtAddtional()//YOU SHOULD REQUIRE THIS FUNCTION AFTER NumIsoLeptons //
{
   bool lepptadd = false;
   if ( TString(Decaymode).Contains( "dimuon" ) ||
        TString(Decaymode).Contains( "dielec" ) ||
        TString(Decaymode).Contains( "muel" )      ){  if ( Lep1->Pt() > 25 && Lep2->Pt() > 20 ) {lepptadd=true;} } //Pt of Leading Lepton should be over than 25 GeV and Second Leading Lepton Pt should be over thand 20 GeV.

   else if ( TString(Decaymode).Contains( "muonJet" ) ){if( v_lepton_idx.size() == 1){lepptadd=true;} }
   else { cout << "?? something wrong " << endl; }  
   return lepptadd;
}
bool ssb_analysis::MuVeto()
{
   bool muveto;
   muveto = false;
   int vetolep = 0;
   if ( TString(Decaymode).Contains( "dimuon" ) )
   {
      for ( int ilumuon = 0; ilumuon < Muon_Count; ++ilumuon )
      {
         if ( ilumuon != v_lepton_idx.at(0) && ilumuon != v_lepton_idx.at(1) )
         {
            if ( Muon_isLoose->at( ilumuon ) == true && Muon_PFIsodBeta04->at( ilumuon ) < 0.25 )
            {   
               vetoMu = v_CorMuons.at(ilumuon);
               if ( vetoMu->Pt() > 20 && fabs(vetoMu->Eta()) < 2.4 ) { vetolep++; }
            }
         }
      }

   }
   else if ( TString(Decaymode).Contains( "dielec" ) )
   {
      // 3rd muon veto
      TLorentzVector *vetoMu = new TLorentzVector();
      for ( int ilumuon = 0; ilumuon < Muon_Count; ++ilumuon )
      {
         if ( Muon_isLoose->at( ilumuon ) == true && Muon_PFIsodBeta04->at( ilumuon ) < 0.25 )
         {   
            vetoMu = v_CorMuons.at(ilumuon);
            if ( vetoMu->Pt() > 20 && fabs(vetoMu->Eta()) < 2.4 ) { vetolep++; }
         }
      }
   }
   else if ( TString(Decaymode).Contains( "muel" ) )
   {
      for ( int ilumuon = 0; ilumuon < Muon_Count; ++ilumuon )
      {        
         if ( ilumuon != v_muon_idx.at(0) )
         {             
            if ( Muon_isLoose->at( ilumuon ) == true && Muon_PFIsodBeta04->at( ilumuon ) < 0.25 )
            {   //vetolep++;
               vetoMu = v_CorMuons.at(ilumuon);
               if ( vetoMu->Pt() > 20 && fabs(vetoMu->Eta()) < 2.4 ) { vetolep++; }
            }  
         }        
      }
   }
   else { cout << "Something Wrong At MuVeto. Check Decaymode-" << endl;}

   if ( vetolep == 0 ){ muveto = true; }

   return muveto;
}
bool ssb_analysis::ElVeto()
{
   bool elveto;
   elveto = false;
   int vetolep = 0;
   if ( TString(Decaymode).Contains( "dimuon" ) )// case of dimuon ...
   {
      for ( int iele = 0; iele < Elec_Count; ++iele )
      {
         vetoEle = v_Elec.at(iele);
         if(  vetoEle->Pt() != vetoEle->Pt()) {continue;}
         if (
             Elec_SCB_Veto->at(iele) == true &&
             vetoEle->Pt() > 20                 &&
             fabs( vetoEle->Eta() ) < 2.5       &&
             Elec_Conversion->at(iele) == true     )
         {
          vetolep++; 
         }

      }
   } // dimuon requriement...
   else if( TString(Decaymode).Contains( "dielec" ) ) 
   {
      // 3rd Electron veto ..
      for ( int iele = 0; iele < Elec_Count; ++iele )
      {
         if ( iele != v_lepton_idx.at(0) && iele != v_lepton_idx.at(1) )
         {
            vetoEle = v_Elec.at(iele);
            if (
                 Elec_SCB_Veto->at(iele) == true &&
                 vetoEle->Pt()            > 20      &&
                 fabs( vetoEle->Eta() )   < 2.5     &&
                 Elec_Conversion->at(iele) == true     )
            { 
               if ( fabs( Elec_Supercluster_Eta->at(iele) ) <= 1.479 )
               {
                  if (
                       v_el_iso_jcl->at( iele ) < 0.175 &&
                       fabs( Elec_Track_GsfdXY->at(iele) ) < 0.05 &&
                       fabs( Elec_Track_GsfdZ->at(iele) ) < 0.10 
                     ){ vetolep++; }
               
               } 
               else 
               {
                  if (
                       v_el_iso_jcl->at( iele ) < 0.159 &&           
                       fabs( Elec_Track_GsfdXY->at(iele) ) < 0.10 &&
                       fabs( Elec_Track_GsfdXY->at(iele) ) < 0.20 
                     ){ vetolep++; }
               }
            }
         }
      }

   } // case of dielectron...

   else if ( TString(Decaymode).Contains( "muel" ) ) 
   {
      for ( int iele = 0; iele < Elec_Count; ++iele )
      {
         if ( iele != v_electron_idx.at(0) )
         {
            vetoEle = v_Elec.at(iele);
            if(  vetoEle->Pt() != vetoEle->Pt()) {continue;}

            if (
                 Elec_SCB_Veto->at(iele) == true &&
                 vetoEle->Pt()            > 20   &&
                 fabs( vetoEle->Eta() )   < 2.5  &&
                 Elec_Conversion->at(iele) == true 
               )
            {
               if ( fabs( Elec_Supercluster_Eta->at(iele) ) <= 1.479 )
               {
                  if (
                       v_el_iso_jcl->at( iele ) < 0.175 &&
                       fabs( Elec_Track_GsfdXY->at(iele) ) < 0.05 &&
                       fabs( Elec_Track_GsfdZ->at(iele) ) < 0.10 
                     )
                  { vetolep++; }
               
               } 
               else 
               {
                  if (
                       v_el_iso_jcl->at( iele ) < 0.159 &&           
                       fabs( Elec_Track_GsfdXY->at(iele) ) < 0.10 &&
                       fabs( Elec_Track_GsfdXY->at(iele) ) < 0.20 
                     )
                  { vetolep++; }
               }              
            }
         }
      }
   }
   else {cout << "EleVto check !!! " << endl;}
   if ( vetolep == 0 ){ elveto = true; }

   return elveto;
}
bool ssb_analysis::ThirdLeptonVeto()
{
   bool third_veto = true;
   if ( TString(Decaymode).Contains( "dimuon" ) || TString(Decaymode).Contains( "dielec" )    )
   { 
      if ( v_lepton_idx.size() <= 1 ) {third_veto = false;;}
      //if ( v_lepton_idx.size() != 2 ) {third_veto = false;;}

      if( TString(Decaymode).Contains( "dimuon" ) ) 
      { 
         //cout << "v_lepton_idx size : " << v_lepton_idx.size() << endl;
         if ( Muon_Charge->at( v_lepton_idx.at(0) ) == Muon_Charge->at( v_lepton_idx.at(1) ) ) {third_veto = false;;}
         if ( MuVeto() == false ) {third_veto = false;;}
         if ( ElVeto() == false ) {third_veto = false;;}
      } 
      else if( TString(Decaymode).Contains( "dielec" ) ) 
      {
         if ( Elec_Charge->at( v_lepton_idx.at(0) ) == Elec_Charge->at( v_lepton_idx.at(1) ) ) {third_veto = false;;}
         if ( MuVeto() == false ) {third_veto = false;;}
         if ( ElVeto() == false ) {third_veto = false;;}
      }
   } //Requiring  Two isolated Lepton .

   else if ( TString(Decaymode).Contains( "muel" ) )
   { 
      if (v_muon_idx.size() < 1 || v_electron_idx.size() < 1 ) { third_veto = false;; } 
      if ( Muon_Charge->at( v_muon_idx.at(0) ) == Elec_Charge->at( v_electron_idx.at(0) ) ) { third_veto = false;; }
      if ( MuVeto() == false ) {third_veto = false;;}
      if ( ElVeto() == false ) {third_veto = false;;}
   }

   else { cout << "something wrong at ThirdLeptonVeto veto Function !!!!" << endl; }

   return third_veto;
}
void ssb_analysis::LeptonOrder()
{
   if ( TString(Decaymode).Contains( "dimuon" ) )
   {
      if (v_lepton_idx.size() > 1)
      {
         Lep1 = v_CorMuons.at( v_lepton_idx.at(0) );
         Lep2 = v_CorMuons.at( v_lepton_idx.at(1) );
         if (Muon_Charge->at( v_lepton_idx.at(0) ) < 0 )
         {
            Lep   = v_CorMuons.at( v_lepton_idx.at(0) );
            AnLep = v_CorMuons.at( v_lepton_idx.at(1) );
         }
         else
         {
            Lep   = v_CorMuons.at( v_lepton_idx.at(1) );
            AnLep = v_CorMuons.at( v_lepton_idx.at(0) );
         }
      }
   }
 
   else if ( TString(Decaymode).Contains( "dielec" ) )
   {
      if (v_lepton_idx.size() > 1)
      {
         Lep1 = v_Elec.at( v_lepton_idx.at(0) );
         Lep2 = v_Elec.at( v_lepton_idx.at(1) );
         if (Elec_Charge->at( v_lepton_idx.at(0) ) < 0 )
         {
            Lep   = v_Elec.at( v_lepton_idx.at(0) );
            AnLep = v_Elec.at( v_lepton_idx.at(1) );
         }
         else
         {
            Lep   = v_Elec.at( v_lepton_idx.at(1) );
            AnLep = v_Elec.at( v_lepton_idx.at(0) );
         }
      }
   }
   else if ( TString(Decaymode).Contains( "muel" ) )
   {
      if ( v_muon_idx.size() > 0 && v_electron_idx.size() > 0 )
      {
         TMuon = v_CorMuons.at( v_muon_idx.at(0) );
         TElectron = v_Elec.at( v_electron_idx.at(0) );
         if ((*TMuon).Pt() > (*TElectron).Pt())
         {
            Lep1 = v_CorMuons.at( v_muon_idx.at(0) );
            Lep2 = v_Elec.at( v_electron_idx.at(0) );
         }
         else
         {
            Lep1 = v_Elec.at( v_electron_idx.at(0) );
            Lep2 = v_CorMuons.at( v_muon_idx.at(0) );
         }
         if ( Muon_Charge->at( v_muon_idx.at(0) ) < 0 )
         {
            Lep   = v_CorMuons.at( v_muon_idx.at(0) );
            AnLep = v_Elec.at( v_electron_idx.at(0) );
         }
         else 
         {
            AnLep = v_CorMuons.at( v_muon_idx.at(0) );
            Lep   = v_Elec.at( v_electron_idx.at(0) );
         }
      }
   }
   else {cout << "Lepton TLorentzVector Error" << endl;}
}

void ssb_analysis::LeptonGetSF()
{

   lep_eff = 1.0;
   if ( TString(Decaymode).Contains( "dimuon" ) )
   {
      lep_eff = SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,LepIdSFSys,LepIsoSFSys,LepTrackSFSys); 
   }
   else if ( TString(Decaymode).Contains( "dielec" ) )
   {
      //lep_eff = SSBEffcal->DoubleElec_EffROOT(Lep1,Lep2,Elec_Supercluster_Eta->at(v_lepton_idx[0]) ,Elec_Supercluster_Eta->at(v_lepton_idx[1]) ,LepIdSFSys,LepIsoSFSys);
      lep_eff = SSBEffcal->DoubleElec_EffROOT(Lep1,Lep2,Elec_Supercluster_Eta->at(v_lepton_idx[0]) ,Elec_Supercluster_Eta->at(v_lepton_idx[1]) ,LepIdSFSys,LepRecoSFSys);// LepRecoSFSys is for Lepton //
   }
   else if ( TString(Decaymode).Contains( "muel" ) )
   {
      //lep_eff = 1.0;
      lep_eff = SSBEffcal->MuonElec_EffROOT(TMuon,TElectron,Elec_Supercluster_Eta->at(v_electron_idx[0]) ,LepIdSFSys,LepIsoSFSys,LepTrackSFSys,LepRecoSFSys);//LepRecoSFSys is for Electron
   }
   else { lep_eff = 1.0; cout << "LeptonGetSF error !!!!" << endl; }
}
void ssb_analysis::LeptonSFApply()
{

   double lep_sf_cent = 0.0;

   evt_weight_beforeLepsf_ = 1.0;
   evt_weight_beforeLepsf_ = evt_weight_; // keep event weight 
   if ( !TString(FileName_).Contains( "Data") )
   {
      LeptonGetSF(); // Getting Lepton SF //
      evt_weight_ = evt_weight_*lep_eff;
   }
   else {evt_weight_ = 1;}
}
// Jet Selection and Jet cleaning 
void ssb_analysis::JetSelector()
{
   v_jet_idx.clear();
   v_jetup_idx.clear();
   v_jetdn_idx.clear();
   v_jetresup_idx.clear();
   v_jetresdn_idx.clear();

   v_jet_TL.clear();
   v_jetup_TL.clear();
   v_jetdn_TL.clear();
   v_jetresup_TL.clear();
   v_jetresdn_TL.clear();


   Int_t injet = Jet->GetEntriesFast();

   TJet = new TLorentzVector();

   TLorentzVector *JetUp =new TLorentzVector();
   TLorentzVector *JetDown = new TLorentzVector();

   TLorentzVector *JetCen = new TLorentzVector();
   TLorentzVector *JetUpPre = new TLorentzVector();
   TLorentzVector *JetDownPre = new TLorentzVector();
   TLorentzVector* TJetResUp;// = new TLorentzVector();
   TLorentzVector* TJetResDn;// = new TLorentzVector();

   TLorentzVector *Jetpre = new TLorentzVector();

   double jetpt_  = -999; 
   double jeteta_ = -999; 
   double jetphi_ = -999; 
   double jetenergy_ = -999; 
   for (int i =0; i < injet; i++)
   {

      TJet = new TLorentzVector();
      Jetpre = new TLorentzVector();
      Jetpre = (TLorentzVector*)Jet->At(i);
      //////////////////////////////
      /// JES Systematic Method1 ///
      //////////////////////////////
      if      ( TString(JetEnSys).Contains("JetEnNorm"     ) ) { Jetpre = (TLorentzVector*)Jet->At(i); }
      else if ( TString(JetEnSys).Contains("JetEnShiftedUp") )
      {
         jetpt_     = Jetpre->Pt();
         jeteta_    = Jetpre->Eta();
         jetphi_    = Jetpre->Phi();
         jetenergy_ = Jetpre->Energy();
         /// Apply TJet
         Jetpre->SetPtEtaPhiE(jetpt_*Jet_EnShiftedUp->at(i),jeteta_,jetphi_,jetenergy_*Jet_EnShiftedUp->at(i));
      }
      else if ( TString(JetEnSys).Contains("JetEnShiftedDown") )
      {
         jetpt_     = Jetpre->Pt();
         jeteta_    = Jetpre->Eta();
         jetphi_    = Jetpre->Phi();
         jetenergy_ = Jetpre->Energy();

         Jetpre->SetPtEtaPhiE(jetpt_*Jet_EnShiftedDown->at(i),jeteta_,jetphi_,jetenergy_*Jet_EnShiftedDown->at(i));
      }
      else { cout << "Check out Jet Systematic Config ... Default is Norminal Jet ..."<< endl; Jetpre = (TLorentzVector*)Jet->At(i); }

      //////////////////////////////
      /// JER Systematic Method1 ///
      //////////////////////////////
      if      ( TString(JetResSys).Contains("JetResNorm"     ) ) { TJet = JERSmearing(Jetpre,i,"Cen"); }
      else if ( TString(JetResSys).Contains("JetResShiftedUp") )
      {
         TJet = JERSmearing(Jetpre,i,"Up");
      }
      else if ( TString(JetResSys).Contains("JetResShiftedDown") )
      {
         TJet = JERSmearing(Jetpre,i,"Down");
      }
      else { cout << "Check out Jet Systematic Config ... Default is Norminal Jet ..."<< endl; TJet = JERSmearing(Jetpre,i,"Cen"); }
    
/*      /// Jet Pt Resolution Dilution Study ///
      if (JetPtResDil != "None" &&  JetPtResDil != "none" ) {
         /// You should check your systematic options, when you study Dilution... 
         //  These part won't be correct, if you use other systematic options... 
         //  So we have to check your systematic options and dilution study options //
         TJet = ApplyJetPtRes(Jetpre,i,JetPtResDil);
      }*/
      //cout << "After Pt Res Jet Pt : " << TJet->Pt() << endl;    
      if ( JetCleaning(TJet) == true ) {
        /* if (JetPhiResDil != "None" &&  JetPhiResDil != "none" ) 
         {
         //   cout << "JetPhiResDil : " << JetPhiResDil << endl;
            TJet = ApplyJetPhiRes(TJet,i,JetPhiResDil);
         }*/
         //cout << "After Phi Res Jet Pt : " << TJet->Pt() << endl;
         v_jetdpt_res.push_back(TJet->Pt()-Jetpre->Pt());  
         v_jetdpx_res.push_back(TJet->Px()-Jetpre->Px());  
         v_jetdpy_res.push_back(TJet->Py()-Jetpre->Py());  
         ////////////////////////////////
         /// Apply Requirement of Jet ///
         ////////////////////////////////
         if ( (*TJet).Pt() > jet_pt && fabs( (*TJet).Eta() ) < jet_eta && v_jet_Id->at(i) >= jet_id )
         {
            v_jet_idx.push_back( i ); 
            v_jet_TL.push_back( TJet ); 
            if ( (*TJet).Pt() < 30. ){cout << "Something wrong " << endl;}
         }
      }

      ////////////////////////////////////////
      /// For JES Systematic with Method 2 ///
      ////////////////////////////////////////

      JetUp =new TLorentzVector();
      JetDown = new TLorentzVector();

      JetCen = new TLorentzVector();
      JetUpPre = new TLorentzVector();
      JetDownPre = new TLorentzVector();
      JetCen = (TLorentzVector*)Jet->At(i); 

      JetUpPre->SetPtEtaPhiE(JetCen->Pt()*Jet_EnShiftedUp->at(i),JetCen->Eta(),JetCen->Phi(),JetCen->Energy()*Jet_EnShiftedUp->at(i));
      JetDownPre->SetPtEtaPhiE(JetCen->Pt()*Jet_EnShiftedDown->at(i),JetCen->Eta(),JetCen->Phi(),JetCen->Energy()*Jet_EnShiftedDown->at(i));
            
      JetUp = JERSmearing(JetUpPre,i,"Cen");
      if (JetCleaning(JetUp) == true ) { 
         v_jetdpt_jesup.push_back(JetUp->Pt()-JetUpPre->Pt());  
         v_jetdpx_jesup.push_back(JetUp->Px()-JetUpPre->Px());  
         v_jetdpy_jesup.push_back(JetUp->Py()-JetUpPre->Py());  
         if ( (*JetUp).Pt() > jet_pt && fabs( (*JetUp).Eta() ) < jet_eta && v_jet_Id->at(i) >= jet_id )
         {
            v_jetup_idx.push_back( i ); 
            v_jetup_TL.push_back( JetUp ); 
            if ( (*JetUp).Pt() < 30. ){cout << "Something wrong " << endl;}
         }
      }     
      JetDown = JERSmearing(JetDownPre,i,"Cen");
      if (JetCleaning(JetDown) == true ) { 
         v_jetdpt_jesdn.push_back(JetDown->Pt()-JetDownPre->Pt());  
         v_jetdpx_jesdn.push_back(JetDown->Px()-JetDownPre->Px());  
         v_jetdpy_jesdn.push_back(JetDown->Py()-JetDownPre->Py());  
         if ( (*JetDown).Pt() > jet_pt && fabs( (*JetDown).Eta() ) < jet_eta && v_jet_Id->at(i) >= jet_id )
         {
            v_jetdn_idx.push_back( i ); 
            v_jetdn_TL.push_back( JetDown ); 
            if ( (*JetDown).Pt() < 30. ){cout << "Something wrong " << endl;}
         }
      }

      ////////////////////////////////////////
      /// For JER Systematic with Method 2 ///
      ////////////////////////////////////////

      TJetResUp = new TLorentzVector();
      TJetResDn = new TLorentzVector();

      TJetResUp = JERSmearing(Jetpre,i,"Up");
      if ( JetCleaning(TJetResUp) == true )
      {
         v_jetdpt_resup.push_back(TJetResUp->Pt()-Jetpre->Pt());  
         v_jetdpx_resup.push_back(TJetResUp->Px()-Jetpre->Px());  
         v_jetdpy_resup.push_back(TJetResUp->Py()-Jetpre->Py());  
         
         if ( (*TJetResUp).Pt() > jet_pt && fabs( (*TJetResUp).Eta() ) < jet_eta && v_jet_Id->at(i) >= jet_id )
         {
            v_jetresup_idx.push_back( i ); 
            v_jetresup_TL.push_back( TJetResUp ); 
            if ( (*TJetResUp).Pt() < 30. ){cout << "Something wrong " << endl;}
         }
      }

      TJetResDn = JERSmearing(Jetpre,i,"Down");
      if ( JetCleaning(TJetResDn) == true )
      {
         v_jetdpt_resdn.push_back(TJetResDn->Pt()-Jetpre->Pt());  
         v_jetdpx_resdn.push_back(TJetResDn->Px()-Jetpre->Px());  
         v_jetdpy_resdn.push_back(TJetResDn->Py()-Jetpre->Py());  
         
         if ( (*TJetResDn).Pt() > jet_pt && fabs( (*TJetResDn).Eta() ) < jet_eta && v_jet_Id->at(i) >= jet_id )
         {
            v_jetresdn_idx.push_back( i ); 
            v_jetresdn_TL.push_back( TJetResDn ); 
            if ( (*TJetResDn).Pt() < 30. ){cout << "Something wrong " << endl;}
         }
      }
   }

//   JetCleaning(v_jet_idx_cand, v_jet_TL_cand, v_jet_idx, v_jet_TL);
   if ( TString(FileName_).Contains( "Data" ) )
   {
      v_jetdn_idx.clear();
      v_jetup_idx.clear();
      v_jetresdn_idx.clear();
      v_jetresup_idx.clear();

      v_jetdn_TL.clear();
      v_jetup_TL.clear();
      v_jetresdn_TL.clear();
      v_jetresup_TL.clear();

      v_jetdn_idx = v_jet_idx;
      v_jetup_idx = v_jet_idx;
      v_jetresdn_idx = v_jet_idx;
      v_jetresup_idx = v_jet_idx;

      v_jetdn_TL = v_jet_TL;
      v_jetup_TL = v_jet_TL;
      v_jetresdn_TL = v_jet_TL;
      v_jetresup_TL = v_jet_TL;
   }

   return;
}
void ssb_analysis::JetCleaning(std::vector<int> v_jidx_cand, std::vector<TLorentzVector*>v_jtl_cand, std::vector<int>& v_jidx, std::vector<TLorentzVector*>&v_jtl )
{
   bool overlapjetlep = false;
   v_jidx.clear();
   v_jtl.clear();
   
   TLorentzVector* TmpJet = new TLorentzVector();
   for ( int i =0; i < v_jidx_cand.size(); ++i ) // Jet Loop (before jet cleanning )
   {
      l_jet = 0;
      overlapjetlep = false;
      TmpJet = new TLorentzVector();
      TmpJet = v_jtl_cand.at(i);
      // electron
      if ( v_ele_idx_cand.size() > 0 )
      {
         for ( int j = 0; j < v_ele_idx_cand.size(); ++j )
         {
            selected_Ele_jet_clean = new TLorentzVector();
            selected_Ele_jet_clean = v_Elec.at( v_ele_idx_cand.at(j) );
            if ( (*selected_Ele_jet_clean).DeltaR( (*TmpJet) ) < 0.4){ overlapjetlep = true; }
         }
      }

      // muon
      if ( v_mu_idx_cand.size() > 0 )
      {  
         for ( int k = 0; k < v_mu_idx_cand.size(); ++k )
         {
            selected_Mu_jet_clean = new TLorentzVector();
            selected_Mu_jet_clean = v_CorMuons.at( v_mu_idx_cand.at(k) );
            if ( (*selected_Mu_jet_clean).DeltaR( (*TmpJet) ) < 0.4){ overlapjetlep = true; }
         }
      }
      if ( overlapjetlep != true  ) 
      {
         v_jidx.push_back( v_jidx_cand.at(i) );
         v_jtl.push_back( v_jtl_cand.at(i) );
      }
   }
   return;
}
bool ssb_analysis::JetCleaning(TLorentzVector* jet_)
{
   bool overlapjetlep = false;
   // electron
   for ( int j = 0; j < v_ele_idx_cand.size(); ++j )
   {
      selected_Ele_jet_clean = new TLorentzVector();
      selected_Ele_jet_clean = v_Elec.at( v_ele_idx_cand.at(j) );
      if ( (*selected_Ele_jet_clean).DeltaR( (*jet_) ) < 0.4){ overlapjetlep = true; }
   }
   // muon
   for ( int k = 0; k < v_mu_idx_cand.size(); ++k )
   {
      selected_Mu_jet_clean = new TLorentzVector();
      selected_Mu_jet_clean = v_CorMuons.at( v_mu_idx_cand.at(k) );
      if ( (*selected_Mu_jet_clean).DeltaR( (*jet_) ) < 0.4){ overlapjetlep = true; }
   }
   return !overlapjetlep; 
}
void ssb_analysis::JetDefiner()
{
   Jet1 = new TLorentzVector();
   Jet2 = new TLorentzVector();

   if (v_jet_TL.size() >=1){
      Jet1 = v_jet_TL[0];
      if (v_jet_TL.size() > 1){ Jet2 = v_jet_TL[1]; }
   }
   return;
}
void ssb_analysis::METDefiner()
{
 
   Met = new TLorentzVector();
   /// For Data, We don't need to apply JES systematic ///  
   Met = (TLorentzVector*)MET->At(0);  

}

TLorentzVector* ssb_analysis::METSmear(std::vector<double> v_dpt, std::vector<double> v_dpx, std::vector<double> v_dpy, TLorentzVector* met_ )
{
   if ( dojer&!TString(FileName_).Contains( "Data" ) ){ 
      TLorentzVector* smearmet_ = new TLorentzVector();
      double sum_dpt =0.;
      double sum_dpx =0.;
      double sum_dpy =0.;
      // check dpt
      for (int i = 0; i< v_dpt.size(); ++i )
      {   
         sum_dpt += v_dpt[i];
         sum_dpx += v_dpx[i];
         sum_dpy += v_dpy[i];
      }
  
      smearmet_->SetPxPyPzE(met_->Px()-sum_dpx,met_->Py()-sum_dpy,0,0);
      return smearmet_;
   }
   else { return met_;}
}
// DiMuon Mass Cut : step 1
bool ssb_analysis::DiLeptonMassCut()
{
   dimu_masscut = false;

   if ( ((*Lep1)+(*Lep2)).M() > 20 ){ dimu_masscut = true; }

   return dimu_masscut;
}

// ZVeto Cut : step 2
bool ssb_analysis::ZVetoCut()
{
   zvetocut = false;

   if ( TString(Decaymode).Contains( "dimuon" ) || TString(Decaymode).Contains( "dielec" ) ) 
   {
      if ( ((*Lep1)+(*Lep2)).M() <= 76 || ((*Lep1)+(*Lep2)).M() >= 106 ){ zvetocut = true; }
   }
   else if ( TString(Decaymode).Contains( "muel" ) ){ zvetocut = true; }
   else {cout << "ZVeto Error !!" << endl;}

   return zvetocut;
}

// Num.Jet Cut : Step 3
bool ssb_analysis::NumJetCut(std::vector<int> v_jets)
{
   numjetcut = false;
   //if ( v_jets.size() == 2 ){ numjetcut = true; }
   if ( v_jets.size() >= 2 ){ numjetcut = true; }
   return numjetcut;
}

//MET Cut : step 4
bool ssb_analysis::METCut(TLorentzVector* met)
{
   metcut = false;
   if ( TString(Decaymode).Contains( "dimuon" ) || TString(Decaymode).Contains( "dielec" ) )
   {
      if (met->Pt() > 40) { metcut =true; }
   }
   else if ( TString(Decaymode).Contains( "muel" ) ){ metcut =true; }
   else {cout << "METCut Error !!" << endl;}
   return metcut;
}
void ssb_analysis::BDsicApply()
{
   v_bjet_idx.clear();
   v_bjetup_idx.clear();
   v_bjetdn_idx.clear();
   v_bjetresup_idx.clear();
   v_bjetresdn_idx.clear();
   v_bjet_TL.clear();
   v_bjetup_TL.clear();
   v_bjetdn_TL.clear();
   v_bjetresup_TL.clear();
   v_bjetresdn_TL.clear();
   nbtagged =0;
/*   for (int ije = 0; ije < v_jet_idx.size(); ije++)
   {

       if (Jet_bDisc->at( v_jet_idx.at(ije) ) > bdisccut ) { nbtagged++;v_bjet_idx.push_back( v_jet_idx.at(ije) ); v_bjet_TL.push_back( v_jet_TL.at(ije) ); }
   }
   /// For JES Up 
   for (int ije = 0; ije < v_jetup_idx.size(); ije++)
   {
       if (Jet_bDisc->at( v_jetup_idx.at(ije) ) > bdisccut ) { v_bjetup_idx.push_back( v_jetup_idx.at(ije) ); v_bjetup_TL.push_back( v_jetup_TL.at(ije) ); }
   }
   /// For JES Down 
   for (int ije = 0; ije < v_jetdn_idx.size(); ije++)
   {
       if (Jet_bDisc->at( v_jetdn_idx.at(ije) ) > bdisccut ) { v_bjetdn_idx.push_back( v_jetdn_idx.at(ije) ); v_bjetdn_TL.push_back( v_jetdn_TL.at(ije) ); }
   }
   /// For JER Up 
   for (int ije = 0; ije < v_jetresup_idx.size(); ije++)
   {
       if (Jet_bDisc->at( v_jetresup_idx.at(ije) ) > bdisccut ) { v_bjetresup_idx.push_back( v_jetresup_idx.at(ije) ); v_bjetresup_TL.push_back( v_jetresup_TL.at(ije) ); }
   }
   /// For JER Down 
   for (int ije = 0; ije < v_jetresdn_idx.size(); ije++)
   {
       if (Jet_bDisc->at( v_jetresdn_idx.at(ije) ) > bdisccut ) { v_bjetresdn_idx.push_back( v_jetresdn_idx.at(ije) ); v_bjetresdn_TL.push_back( v_jetresdn_TL.at(ije) ); }
   }
*/
   // find index for "pfDeepFlavourJetTags:probb" , "pfDeepFlavourJetTags:probbb" , "pfDeepFlavourJetTags:problepb"
   //cout << "Jet_bDisc->size() : " << Jet_bDisc->size() << endl;
   //cout << "Jet_bDisc_Name->size() : " << Jet_bDisc_Name->size() << endl;
   //cout << "Jet_bDisc_Value->size() : " << Jet_bDisc_Value->size() << endl;
   int idx_probb = -1;
   int idx_probbb = -1;
   int idx_problepb = -1;
   int idx_oldbtag = -1;// for debug//
   idx_probb = find(Jet_bDisc_Name->begin(), Jet_bDisc_Name->end(), "pfDeepCSVJetTags:probb") - Jet_bDisc_Name->begin();
   idx_probbb = find(Jet_bDisc_Name->begin(), Jet_bDisc_Name->end(), "pfDeepCSVJetTags:probbb") - Jet_bDisc_Name->begin();
   //idx_problepb = find(Jet_bDisc_Name->begin(), Jet_bDisc_Name->end(), "pfDeepCSVJetTags:problepb") - Jet_bDisc_Name->begin();
   idx_oldbtag = find(Jet_bDisc_Name->begin(), Jet_bDisc_Name->end(), "pfCombinedInclusiveSecondaryVertexV2BJetTags") - Jet_bDisc_Name->begin();
   
   //cout << "idx_probb : " << idx_probb << endl;
   //cout << "idx_probbb : " << idx_probbb << endl;
   //cout << "idx_problepb : " << idx_problepb << endl;
   //cout << "idx_oldbtag : " << idx_oldbtag << endl;
   for(int ijet = 0; ijet < v_jet_idx.size(); ++ijet){
      //cout << "Jet_bDisc : " << Jet_bDisc->at(ijet) << endl;
      //cout << "Jet_bDisc_Value: " << Jet_bDisc_Value->at(ijet*14 + idx_oldbtag) << endl;
      double btag_= Jet_bDisc_Value->at( ijet*14 + idx_probb ) + Jet_bDisc_Value->at( ijet*14 + idx_probbb ); //+ Jet_bDisc_Value->at( ijet*14 + idx_problepb );
      //cout << "btag_ " << btag_ << endl;
      if ( btag_ > bdisccut ) { nbtagged++; v_bjet_idx.push_back( v_jet_idx.at(ijet) ); v_bjet_TL.push_back( v_jet_TL.at(ijet) ); }
   }

 //itJet.bDiscriminator("pfDeepFlavourJetTags:probb") + itJet.bDiscriminator("pfDeepFlavourJetTags:probbb") + itJet.bDiscriminator("pfDeepFlavourJetTags:problepb");

}
//B Jet Cut : step 5
bool ssb_analysis::BJetCut( std::vector<int>v_bjets )
{
   bjetcut = false;
   if (v_bjets.size() >= 1) { bjetcut = true; }
   return bjetcut;
}
void ssb_analysis::BTaggigSFApply()
{
   evt_weight_beforeBtag_ = 1;
   evt_weight_beforeBtag_ = evt_weight_;
   if( TString(FileName_).Contains( "Data" ) ) {return;}
   std::vector<double> v_Jet_pT;
   std::vector<double> v_Jet_eta;
   std::vector<double> v_Jet_bDisc;
   std::vector<int> v_Jet_flav;
   v_Jet_pT.clear();
   v_Jet_eta.clear();
   v_Jet_bDisc.clear();
   v_Jet_flav.clear();
   TLorentzVector* jet_ = new TLorentzVector();
   for ( int ijet = 0; ijet < v_jet_idx.size(); ++ijet )
   {
      int idx_ = v_jet_idx[ijet];
      jet_ = (TLorentzVector*)v_jet_TL[ijet];
      v_Jet_pT.push_back(jet_->Pt()); 
      v_Jet_eta.push_back(jet_->Eta()); 
      v_Jet_bDisc.push_back(Jet_bDisc->at(idx_)); 
      v_Jet_flav.push_back(abs(Jet_HadronFlavour->at(idx_))); 
   } 

   //if( !TString(FileName_).Contains( "Data" ) ) {evt_weight_ = SSBEffcal->Btagging_EvenWeight( v_Jet_pT, v_Jet_eta, v_Jet_bDisc, bdisccut,v_Jet_flav )*evt_weight_;}
   if( !TString(FileName_).Contains( "Data" ) ) {evt_weight_ = SSBEffcal->Btagging_EvenWeightv2( v_Jet_pT, v_Jet_eta, v_Jet_bDisc, bdisccut,v_Jet_flav , "central")*evt_weight_;}
   else { evt_weight_ = 1; }

}

void ssb_analysis::BJetDefiner(std::vector<int> v_jets, std::vector<int> v_bjets)
{
   //////////////////
   // Define BJets //
   //////////////////
   bcandi = -9;
   if( v_bjet_TL.size() > 1)
   {
      bJet1 = v_bjet_TL.at(0); 
      bJet2 = v_bjet_TL.at(1);
   }

   else if (  v_bjet_idx.size() > 0 )
   {
      bJet1 = v_bjet_TL.at(0);
      for ( int k = 0; k < v_jet_idx.size(); ++k)
      {
         if ( v_jet_idx.at(k) != v_bjet_idx.at(0) ){ bJet2 =v_jet_TL[k]; break; }
      }
   }
 
   // Bjet TLorentzVector
   else if ( v_bjet_idx.size() == 0 )
   {
      if ( v_jet_idx.size() > 1 )
      {
         bJet1 = v_jet_TL[0];
         bJet2 = v_jet_TL[1];
      }
   }
   else { cout << "you should check out nbjet cut!! " << endl; }
}
void ssb_analysis::BJetDefiner()
{
   //////////////////
   // Define BJets //
   //////////////////
   /// For Central Jet ///
   if( v_bjet_TL.size() > 1)
   {
      bJet1 = v_bjet_TL.at(0); 
      bJet2 = v_bjet_TL.at(1);
   }
   else if ( v_bjet_TL.size() > 0 )
   {  
      bJet1 = v_bjet_TL.at(0);
      for ( int i = 0; i < v_jet_idx.size(); ++i )
      {
         if ( v_jet_idx[i] != v_bjet_idx[0] ){ bJet2 = v_jet_TL[i]; break; }
      }
   }
   else if ( v_bjet_TL.size() == 0 )
   {
      bJet1 = Jet1;
      bJet2 = Jet2;
   }
   else { cout << "you should check out nbjet cut!! Error about v_bjet_TL !! " << endl; }
   /// For JES ShiftedUp Jet ///
   if( v_bjetup_TL.size() > 1)
   {
      bJet1Up = v_bjetup_TL.at(0); 
      bJet2Up = v_bjetup_TL.at(1);
   }
   else if ( v_bjetup_TL.size() > 0 )
   {  
      bJet1Up = v_bjetup_TL.at(0);
      for ( int i = 0; i < v_jetup_idx.size(); ++i )
      {
         if ( v_jetup_idx[i] != v_bjetup_idx[0] ){ bJet2Up = v_jetup_TL[i]; break; }
      }
   }
   else if ( v_bjetup_TL.size() == 0 )
   {
      bJet1Up = Jet1Up;
      bJet2Up = Jet2Up;
   }
   else { cout << "you should check out nbjet cut!! Error about v_bjetup_TL !! " << endl; }

}
/////////////////////
//Double batagging //
/////////////////////
bool ssb_analysis::DoubleBtag( std::vector<int>v_bjets )
{
   doublebtag = false;
   if ( v_bjets.size() >=2 ) {doublebtag = true;}
   return doublebtag;
}
//Kinematic Solver to reconstruct top
void ssb_analysis::KinSol()
{

   xconstraint = 0; yconstraint  = 0;
 
   ksolweight_ = -1.0;
   /// initialize array
   //cout << "At KinSol -- " << "Lep->Px() : " << Lep->Px() << " AnLep->Px() : " << AnLep->Px() << " bJet1->Px() : " << bJet1->Px() << " bJet2->Px() : " << bJet2->Px() << " Met->Px() : " << Met->Px()<< endl;
   xconstraint = Lep->Px() + AnLep->Px() + bJet1->Px() + bJet2->Px() + Met->Px();
   yconstraint = Lep->Py() + AnLep->Py() + bJet1->Py() + bJet2->Py() + Met->Py();

   ssbflsolver->SetConstraints( xconstraint , yconstraint );  

   TtFullLepKinSolver::NeutrinoSolution nuSol1  = ssbflsolver->getNuSolution( (*Lep),   (*AnLep), (*bJet1), (*bJet2) );
   TtFullLepKinSolver::NeutrinoSolution nuSol2  = ssbflsolver->getNuSolution( (*Lep),   (*AnLep), (*bJet2), (*bJet1) );


   if ( nuSol1.weight >nuSol2.weight )
   {
      ksolweight_ = nuSol1.weight;
      (*AnNu)     = nuSol1.neutrino;
      (*Nu)       = nuSol1.neutrinoBar;
      (*bJet)     = (*bJet2);
      (*AnbJet)   = (*bJet1);
   }
   else 
   {  
      ksolweight_ = nuSol2.weight;
      (*AnNu)     = nuSol2.neutrino; 
      (*Nu)       = nuSol2.neutrinoBar;
      (*bJet)     = (*bJet1);
      (*AnbJet)   = (*bJet2);
   }
   //(*Top)      = (*Lep)   + (*AnbJet) + (*AnNu);
   //(*AnTop)    = (*AnLep) + (*bJet)   + (*Nu);
   (*AnTop)     = (*Lep)   + (*AnbJet) + (*AnNu);
   (*Top)       = (*AnLep) + (*bJet)   + (*Nu);
   (*W1)        = (*Lep)   + (*AnNu);
   (*W2)        = (*AnLep) + (*Nu);


}
void ssb_analysis::KinSolSys( TLorentzVector* lep, TLorentzVector* alep, TLorentzVector* bj1, TLorentzVector* bj2, TLorentzVector* met,TString Sys_ )
{

   ksolweight_ = -1.0;
   /// initialize array

   xconstraint = lep->Px() + alep->Px() + bj1->Px() + bj2->Px() + met->Px();
   yconstraint = lep->Py() + alep->Py() + bj1->Py() + bj2->Py() + met->Py();

   ssbflsolver->SetConstraints( xconstraint , yconstraint );  

   TtFullLepKinSolver::NeutrinoSolution nuSol1  = ssbflsolver->getNuSolution( (*lep),   (*alep), (*bj1), (*bj2) );
   TtFullLepKinSolver::NeutrinoSolution nuSol2  = ssbflsolver->getNuSolution( (*lep),   (*alep), (*bj2), (*bj1) );

   if ( nuSol1.weight >nuSol2.weight )
   {
      ksolweight_ = nuSol1.weight;
      (*AnNu)     = nuSol1.neutrino;
      (*Nu)       = nuSol1.neutrinoBar;
      (*bJet)     = (*bj2);
      (*AnbJet)   = (*bj1);
   }
   else 
   {  
      ksolweight_ = nuSol2.weight;
      (*AnNu)     = nuSol2.neutrino; 
      (*Nu)       = nuSol2.neutrinoBar;
      (*bJet)     = (*bj1);
      (*AnbJet)   = (*bj2);
   }
   (*AnTop)     = (*lep)  + (*AnbJet) + (*AnNu);
   (*Top)       = (*alep) + (*bJet)   + (*Nu);
   (*W1)        = (*lep)  + (*AnNu);
   (*W2)        = (*alep) + (*Nu);
}
void ssb_analysis::Bjorken( TLorentzVector *lep, TLorentzVector *alep,TLorentzVector *bj,TLorentzVector *abj,TLorentzVector *nu,TLorentzVector *anu )
{
   x1_bj = -999;
   x2_bj = -999;
   x3_bj = -999;
   x1_bj = ( lep->Pz() + alep->Pz() + bj->Pz() + abj->Pz() + nu->Pz() + anu->Pz() 
           + lep->E()  + alep->E()  + bj->E()  + abj->E()  + nu->E()  + anu->E()  )/(2*4*1000);
   x2_bj = ( lep->E()  + alep->E()  + bj->E()  + abj->E()  + nu->E()  + anu->E() 
           - lep->Pz() - alep->Pz() - bj->Pz() - abj->Pz() - nu->Pz() - anu->Pz() )/(2*4*1000);
   x3_bj = ( x1_bj + x2_bj )*(4*1000);
   
   
}

void ssb_analysis::LepAnLepMisCharge()
{
   laecO3 = 100; laecOb = 100; laecO5 = 100;
   if ( TString(Decaymode).Contains( "dimuon" ) )
   {
      if ( Lep->Pt() < 300 )
      {
         if ( AnLep->Pt() < 300 ){ laecO3 = 1.5; laecOb = 1.5; laecO5 = 1.5; }
         else                    { laecO3 = 3.5; laecOb = 3.5; laecO5 = 3.5; }
      }
      else 
      {
         if ( AnLep->Pt() < 300 ){ laecO3 = 5.5; laecOb = 5.5; laecO5 = 5.5; }
         else                    { laecO3 = 7.5; laecOb = 7.5; laecO5 = 7.5; }
      }
   }
   else if ( TString(Decaymode).Contains( "dielec" ) )
   {
      if( Lep->Pt() >20 && Lep->Pt() < 50 )
      {
         if( fabs( Lep->Eta() ) < 1.479  ) 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  ){
               if( fabs( AnLep->Eta() ) < 1.479 ) { laecO3 = 1.5; laecOb = 1.5; laecO5 = 1.5; }
               else                               { laecO3 = 3.5; laecOb = 3.5; laecO5 = 3.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { laecO3 = 5.5; laecOb = 5.5; laecO5 = 5.5; }
               else                                { laecO3 = 7.5; laecOb = 7.5; laecO5 = 7.5; }
            } //
         } // Lep eta < 1.479
         else // Lep 1.479  < eta < 2.5 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  )
            {
               if( fabs( AnLep->Eta() ) < 1.479  ) { laecO3 = 11.5; laecOb = 11.5; laecO5 = 11.5; }
               else                                { laecO3 = 13.5; laecOb = 13.5; laecO5 = 13.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { laecO3 = 15.5; laecOb = 15.5; laecO5 = 15.5; }
               else                                { laecO3 = 17.5; laecOb = 17.5; laecO5 = 17.5; }

            } //
         }
      }
      else // Lep >= 50GeV 
      {
         if( fabs( Lep->Eta() ) < 1.479  ) 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  ){
               if( fabs( AnLep->Eta() ) < 1.479 ) { laecO3 = 21.5; laecOb = 21.5; laecO5 = 21.5; }
               else                               { laecO3 = 23.5; laecOb = 23.5; laecO5 = 23.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { laecO3 = 25.5; laecOb = 25.5; laecO5 = 25.5; }
               else                                { laecO3 = 27.5; laecOb = 27.5; laecO5 = 27.5; }
            } //
         } // Lep eta < 1.479
         else // Lep 1.479  < eta < 2.5 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  )
            {
               if( fabs( AnLep->Eta() ) < 1.479  ) { laecO3 = 31.5; laecOb = 31.5; laecO5 = 31.5; }
               else                                { laecO3 = 33.5; laecOb = 33.5; laecO5 = 33.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { laecO3 = 35.5; laecOb = 35.5; laecO5 = 35.5; }
               else                                { laecO3 = 37.5; laecOb = 37.5; laecO5 = 37.5; }

            } //
         }
      }
   }
   else if ( TString(Decaymode).Contains( "muel" ) )
   {
      if (TMuon->Pt() < 300 )
      {
         if( TElectron->Pt() > 20 && TElectron->Pt() < 50 )
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { laecO3 = 1.5; laecOb = 1.5; laecO5 = 1.5; } // Barrel
            else                                    { laecO3 = 3.5; laecOb = 3.5; laecO5 = 3.5; } // EndCap
         }
         else // Electron Pt > 50 GeV 
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { laecO3 = 5.5; laecOb = 5.5; laecO5 = 5.5; } // Barrel
            else                                    { laecO3 = 7.5; laecOb = 7.5; laecO5 = 7.5; } // Endcap
         }
      } 
      else 
      {
         if( TElectron->Pt() > 20 && TElectron->Pt() < 50 )
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { laecO3 = 11.5; laecOb = 11.5; laecO5 = 11.5; } // Barrel
            else                                    { laecO3 = 13.5; laecOb = 13.5; laecO5 = 13.5; } // EndCap
         }
         else // Electron Pt > 50 GeV 
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { laecO3 = 15.5; laecOb = 15.5; laecO5 = 15.5; } // Barrel
            else                                    { laecO3 = 17.5; laecOb = 17.5; laecO5 = 17.5; } // Endcap
         }

      }

   }

}
void ssb_analysis::NewLepAnLepMisCharge()
{
   int disc_leptons = 0;
   int value_mis = 0;
   v_MisCharge.clear();
   if ( v_recocp_O.size() != 13) {return;} /// After getting CPVarialbes, You can execute this fucntion ...///
   if ( TString(Decaymode).Contains( "dimuon" ) )
   {
      if ( Lep->Pt() < 300 )
      {
         if ( AnLep->Pt() < 300 ){ disc_leptons = 1.5;}
         else                    { disc_leptons = 3.5;}
      }
      else 
      {
         if ( AnLep->Pt() < 300 ){ disc_leptons = 5.5;}
         else                    { disc_leptons = 7.5;}
      }
   }
   else if ( TString(Decaymode).Contains( "dielec" ) )
   {
      if( Lep->Pt() >20 && Lep->Pt() < 50 )
      {
         if( fabs( Lep->Eta() ) < 1.479  ) 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  ){
               if( fabs( AnLep->Eta() ) < 1.479 ) { disc_leptons = 1.5; }
               else                               { disc_leptons = 3.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { disc_leptons = 5.5; }
               else                                { disc_leptons = 7.5; }
            } //
         } // Lep eta < 1.479
         else // Lep 1.479  < eta < 2.5 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  )
            {
               if( fabs( AnLep->Eta() ) < 1.479  ) { disc_leptons = 11.5; }
               else                                { disc_leptons = 13.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { disc_leptons = 15.5; }
               else                                { disc_leptons = 17.5; }

            } //
         }
      }
      else // Lep >= 50GeV 
      {
         if( fabs( Lep->Eta() ) < 1.479  ) 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  ){
               if( fabs( AnLep->Eta() ) < 1.479 ) { disc_leptons = 21.5; }
               else                               { disc_leptons = 23.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { disc_leptons = 25.5; }
               else                                { disc_leptons = 27.5; }
            } //
         } // Lep eta < 1.479
         else // Lep 1.479  < eta < 2.5 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  )
            {
               if( fabs( AnLep->Eta() ) < 1.479  ) { disc_leptons = 31.5; }
               else                                { disc_leptons = 33.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { disc_leptons = 35.5; }
               else                                { disc_leptons = 37.5; }

            } //
         }
      }
   }
   else if ( TString(Decaymode).Contains( "muel" ) )
   {
      if (TMuon->Pt() < 300 )
      {
         if( TElectron->Pt() > 20 && TElectron->Pt() < 50 )
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { disc_leptons = 1.5; } // Barrel
            else                                    { disc_leptons = 3.5; } // EndCap
         }
         else // Electron Pt > 50 GeV 
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { disc_leptons = 5.5; } // Barrel
            else                                    { disc_leptons = 7.5; } // Endcap
         }
      } 
      else 
      {
         if( TElectron->Pt() > 20 && TElectron->Pt() < 50 )
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { disc_leptons = 11.5; } // Barrel
            else                                    { disc_leptons = 13.5; } // EndCap
         }
         else // Electron Pt > 50 GeV 
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { disc_leptons = 15.5; } // Barrel
            else                                    { disc_leptons = 17.5; } // Endcap
         }

      }

   }
   /// Check Up Positive and Negative !!!
   for ( int i = 0; i < v_recocp_O.size(); ++i )
   {
      if (v_recocp_O[i] > 0.0000) {
         value_mis = disc_leptons;
      }
      else { value_mis = 50+disc_leptons; }
      v_MisCharge.push_back(value_mis);
   }


}
double ssb_analysis::JetPtResol(TLorentzVector *jet)
{
   double jpr_;
   double jpt_;
   double jeta_;
   jpt_ = jet->Pt();
   jeta_ = jet->Eta();
   jpr_ = ssbcpviol->getJPR(jeta_,jpt_);
   return jpr_;
}

double ssb_analysis::CalTopPtRewight(TString Opt_)
{
   double topptreweight =1.0;
   if (Opt_.Contains("None")){return topptreweight;}
   else if (Opt_.Contains("Apply"))
   {
      if ( TString(FileName_).Contains( "TTJets" ) )
      {
         TLorentzVector* genTop   = new TLorentzVector();
         TLorentzVector* genAnTop = new TLorentzVector();
         if (GenTop->GetEntriesFast() >0 && GenAnTop->GetEntriesFast() >0 ){
            genTop   = (TLorentzVector*)GenTop->At(0);
            genAnTop = (TLorentzVector*)GenAnTop->At(0);
         }
         else {
         // For Herwig Case //
            genTop = FindGenPar(6);
            genAnTop = FindGenPar(-6);
         }
         double a = 0.0615;
         double b = -0.0005;
         double toppt = genTop->Pt();
         double antoppt = genAnTop->Pt();
         double topsf = exp(a+b*toppt);
         double antopsf = exp(a+b*antoppt);
         topptreweight = sqrt(topsf*antopsf);
      }
      else {topptreweight = 1.0;}
   }
   else {return 1.0;}
   return topptreweight;
}
void ssb_analysis::TopPtReweightApply()
{
   double evt_weight_beforetopptweight_ = evt_weight_;
   double topptreweight =1.0;
   if ( TString(FileName_).Contains( "TTJets" ) )
   {
      topptreweight = CalTopPtRewight(TopPtSys);
   }
   else {topptreweight =1;}
   evt_weight_ = evt_weight_*topptreweight;

   return;
}
TLorentzVector* ssb_analysis::JERSmearing(TLorentzVector *jet, int idx_, TString op_)
{
   double jerfrac_ = 1.;
   if ( op_.Contains( "Cen") == true ) { jerfrac_ = Jet_EnergyResolution_SF->at(idx_); }
   else if ( op_ == "Up") { jerfrac_ = Jet_EnergyResolution_SFUp->at(idx_); }
   else if ( op_ == "Down") { jerfrac_ = Jet_EnergyResolution_SFDown->at(idx_); }
   else {cout << "Check out your JERSmeraing option !! " << endl;jerfrac_ =1.0; }

   TLorentzVector* jetsmerd = new TLorentzVector();
   if ( dojer&!TString(FileName_).Contains( "Data" ) )
   {
      jetsmerd->SetPtEtaPhiE(jet->Pt()*jerfrac_,
                             jet->Eta(),jet->Phi(),
                            jet->Energy()*jerfrac_);
   }
   else {
      jetsmerd = jet;
   }
   //cout << "before smeared jet pt : " << jet->Pt()  << " Eta : " << jet->Eta() << " phi " << jet->Phi() << " Energy : " << jet->Energy() << endl;
   //cout << "after smeared jetsmerd pt : " << jetsmerd->Pt()  << " Eta : " << jetsmerd->Eta() << " phi " << jetsmerd->Phi() << " Energy : " << jetsmerd->Energy() << endl;
   return  jetsmerd;
}
TLorentzVector* ssb_analysis::ApplyJetPtRes(TLorentzVector *jet, int idx_, TString op_)
{
   double jerfrac_ = 1.;
   if ( op_.Contains( "None") == true ) { return jet; }
   else if ( op_ == "Up") { jerfrac_ = 1+Jet_EnergyResolution_MC->at(idx_); }
   else if ( op_ == "Down") { jerfrac_ = 1-Jet_EnergyResolution_MC->at(idx_); }
   else {cout << "Check out your JERSmeraing option !! " << endl;jerfrac_ =1.0; }
   //cout << "jerfrac_ : " << jerfrac_ << endl;
   TLorentzVector* variedjet = new TLorentzVector();
   if ( !TString(FileName_).Contains( "Data" ) )
   {
      variedjet->SetPtEtaPhiE(jet->Pt()*jerfrac_,
                             jet->Eta(),jet->Phi(),
                            jet->Energy()*jerfrac_);
   }
   else {
      variedjet = jet;
   }
   //cout << "before smeared jet pt : " << jet->Pt()  << " Eta : " << jet->Eta() << " phi " << jet->Phi() << " Energy : " << jet->Energy() << endl;
   //cout << "after smeared variedjet pt : " << variedjet->Pt()  << " Eta : " << variedjet->Eta() << " phi " << variedjet->Phi() << " Energy : " << variedjet->Energy() << endl;
   return  variedjet;
}
TLorentzVector* ssb_analysis::ApplyJetPhiRes(TLorentzVector *jet, int idx_, TString op_)
{
   double jerfrac_ = 1.;
   if ( op_.Contains( "None") == true ) { return jet; }
   else if ( op_ == "Up") { jerfrac_ = 1+Jet_PhiResolution_MC->at(idx_); }
   else if ( op_ == "Down") { jerfrac_ = 1-Jet_PhiResolution_MC->at(idx_); }
   else {cout << "Check out your JERSmeraing option !! " << endl;jerfrac_ =1.0; }
   TLorentzVector* variedjet = new TLorentzVector();
   if ( !TString(FileName_).Contains( "Data" ) )
   {
      variedjet->SetPtEtaPhiE(jet->Pt(),
                             jet->Eta(),jet->Phi()*jerfrac_,
                            jet->Energy());
   }
   else {
      variedjet = jet;
   }
   //cout << "before smeared jet pt : " << jet->Pt()  << " Eta : " << jet->Eta() << " phi " << jet->Phi() << " Energy : " << jet->Energy() << endl;
   //cout << "after smeared variedjet pt : " << variedjet->Pt()  << " Eta : " << variedjet->Eta() << " phi " << variedjet->Phi() << " Energy : " << variedjet->Energy() << endl;
   return  variedjet;
}

bool ssb_analysis::FindSameObj(TLorentzVector* ref, TLorentzVector* tar)
{
   bool sameobj = false;
   if (ref->DeltaR(*tar) == 0 && ref->Pt() == tar->Pt() && ref->Eta() == tar->Eta()){sameobj = true;}
   return sameobj;
}

/*void ssb_analysis::ApplyJetPtRes(TString op_)
{
   double bjerfrac1_ = 1.;
   double bjerfrac2_ = 1.;
   //cout << "Check out v_bjet size "<< v_bjet_idx.size() << endl;
   if (v_bjet_idx.size() <2){cout << "Not 2 bjets... " << endl; return;}
   int idxb1 = v_bjet_idx[0]; 
   int idxb2 = v_bjet_idx[1];

   double bj1pt = bJet1->Pt();
   double bj1eta = bJet1->Eta();
   double bj1phi = bJet1->Phi();
   double bj1energy = bJet1->Energy();
   double bj1px = bJet1->Px();
   double bj1py = bJet1->Py();

   double bj2pt = bJet2->Pt();
   double bj2eta = bJet2->Eta();
   double bj2phi = bJet2->Phi();
   double bj2energy = bJet2->Energy();
   double bj2px = bJet2->Px();
   double bj2py = bJet2->Py();
   std::vector<double> v_diff_pt;
   std::vector<double> v_diff_px;
   std::vector<double> v_diff_py;
   if ( op_.Contains( "None") == true ) { return; }
   else if ( op_ == "Up") 
   {
      bjerfrac1_ = 1+Jet_EnergyResolution_MC->at(idxb1); 
      bjerfrac2_ = 1+Jet_EnergyResolution_MC->at(idxb2); 
   }
   else if ( op_ == "Down") {
      bjerfrac1_ = 1-Jet_EnergyResolution_MC->at(idxb1); 
      bjerfrac2_ = 1-Jet_EnergyResolution_MC->at(idxb2); 
   }
   else {cout << "Check out your JERSmeraing option !! " << endl;bjerfrac1_ =1.0; bjerfrac2_ =1.0; }
   ////////////////////////////
   /// Change bJet1 & bJet2 ///
   ////////////////////////////
   bJet1->SetPtEtaPhiE(bj1pt*bjerfrac1_,bj1eta,bj1phi,bj1energy*bjerfrac1_);
   bJet2->SetPtEtaPhiE(bj2pt*bjerfrac2_,bj2eta,bj2phi,bj2energy*bjerfrac2_);

   cout 
   << "before apply pt res bJet1 pt : " << bj1pt << " eta : "  << bj1eta << " phi : " << bj1phi << " energy : " << bj1energy  
   << endl;
   cout 
   << "after apply pt res bJet1 pt : " << bJet1->Pt() << " eta : "  << bJet1->Eta() << " phi : " << bJet1->Phi() << " energy : " << bJet1->Energy()
   << endl;
   v_diff_pt.push_back(bJet1->Pt() - bj1pt);
   v_diff_px.push_back(bJet1->Px() - bj1px);
   v_diff_py.push_back(bJet1->Py() - bj1py);
   /////////////////////
   /// Calculate MET /// 
   /////////////////////
   Met = METSmear(v_diff_pt,v_diff_px,v_diff_py,Met);
}*/

void ssb_analysis::ClearVectors()
{
   v_lepton_Id->clear();
   v_lepton_iso->clear();
   v_jet_Id->clear();
   v_mlep_iso->clear();
   v_elec_iso->clear();
   v_muon_Id->clear();  // for MuEl channel
   v_electron_Id->clear(); // for ElEl and MuEl channels

   /// User variables for JetCleanning cut
   v_mu_iso_jcl->clear();
   v_mu_Id_jcl->clear();
   v_el_iso_jcl->clear();
   v_el_Id_jcl->clear();
}

void ssb_analysis::SetGenLepAnLep()
{
   GenLep   = new TLorentzVector(); 
   GenAnLep = new TLorentzVector();
   for (int i = 0; i < GenPar_pdgId->size(); ++i)
   {
      if ( GenPar_pdgId->at(i) == 2212 ) {continue;}
      TLorentzVector* par = (TLorentzVector*)GenPar->At(i);
      if ( TString(Decaymode).Contains( "dielec" ) ) {
         if ( GenPar_pdgId->at(i) == 11 ) { GenLep = par; }
         if ( GenPar_pdgId->at(i) == -11 ) { GenAnLep = par; }
      }
      else if ( TString(Decaymode).Contains( "dimuon" ) ) {
         if ( GenPar_pdgId->at(i) == 13 ) { GenLep = par; }
         if ( GenPar_pdgId->at(i) == -13 ) { GenAnLep = par; }
      }
      else if ( TString(Decaymode).Contains( "muel" ) ) {
         if ( GenPar_pdgId->at(i) == 11 ) { GenLep = par; }
         if ( GenPar_pdgId->at(i) == -11 ) { GenAnLep = par; }
         if ( GenPar_pdgId->at(i) == 13 ) { GenLep = par; }
         if ( GenPar_pdgId->at(i) == -13 ) { GenAnLep = par; }
      }
      else { cout << "Check Your Decay Mode !!!" << endl;return;}
   }
//   cout << "GenLep pt : " << GenLep->Pt() << endl;
//   cout << "GenAnLep pt : " << GenAnLep->Pt() << endl;

}
TLorentzVector* ssb_analysis::FindGenPar(int pdgid_)
{
    TLorentzVector* tarpar_ = new TLorentzVector();
    bool find_par = false;
    for (int i = 0; i < GenPar_pdgId->size(); ++i)
    {
       if ( GenPar_pdgId->at(i) == pdgid_ ) { tarpar_ = (TLorentzVector*)GenPar->At(i); find_par = true;} 
    }
    if (find_par == false)cout << "There is no target partice !! pdgId : " << pdgid_  << endl;
    return tarpar_;
}

void ssb_analysis::BandBbarJetDiff()
{
   /// This Function is for dilution factor study of difference of response between b and bbar jet ///
   /// Fraction of DeltaR is 0.001.
   if ( TString(FileName_).Contains( "Data" ) ){return;}
   if ( bAndbBarDil == false ){return;}
//   cout << "bAndbBarDil : " << bAndbBarDil  << endl;
   double frac = 0.0011;
   //double bjpx_,bjpy_,bjpz_,bje_; 
   //double anbjpx_,anbjpy_,anbjpz_,anbje_;
   double bjpt_,bjeta_,bjphi_,bje_; 
   double anbjpt_,anbjeta_,anbjphi_,anbje_;
 
   double fac_bj   = 1-0.5*0.0011;
   double fac_anbj = 1+0.5*0.0011;
   // Get BJet 4 - momentum
   /*bjpx_ = bJet->Px();
   bjpy_ = bJet->Py();
   bjpz_ = bJet->Pz();*/
   bjpt_ = bJet->Pt();
   bjeta_ = bJet->Eta();
   bjphi_ = bJet->Phi();
   bje_ = bJet->Energy();

   // Get AnBJet 4 - momentum
   /*anbjpx_ = AnbJet->Px();
   anbjpy_ = AnbJet->Py();
   anbjpz_ = AnbJet->Pz();*/
   anbjpt_ = AnbJet->Pt();
   anbjeta_ = AnbJet->Eta();
   anbjphi_ = AnbJet->Phi();
   anbje_ = AnbJet->Energy();

   /*cout << " ----------------- before bAndbbarJet Diff---------------------- " << endl;
   cout << "bJet px : " << bJet->Px() << " py : " << bJet->Py() << " pz : " <<  bJet->Pz() << " energy : " << bJet->Energy() 
   << "AnbJet px : " << AnbJet->Px() << " py : " << AnbJet->Py() << " pz : " <<  AnbJet->Pz() << " energy : " << AnbJet->Energy() << endl;*/

   // Apply Diff. of Response
   //bJet->SetPxPyPzE(bjpx_*fac_bj, bjpy_*fac_bj , bjpz_*fac_bj, bje_*fac_bj);
   bJet->SetPtEtaPhiE(bjpt_*fac_bj, bjeta_, bjphi_, bje_*fac_bj);
   //AnbJet->SetPxPyPzE(anbjpx_*fac_anbj,anbjpy_*fac_anbj,anbjpz_*fac_anbj,anbje_*fac_anbj);
   AnbJet->SetPtEtaPhiE(anbjpt_*fac_anbj,anbjeta_,anbjphi_,anbje_*fac_anbj);

   // For Debug 
   /*cout << " ----------------- after bAndbbarJet Diff---------------------- " << endl;
   cout << " fac_bj ? : " << fac_bj << " fac_anbj : " << fac_anbj << endl;
   cout << "bJet px : " << bJet->Px() << " py : " << bJet->Py() << " pz : " <<  bJet->Pz() << " energy : " << bJet->Energy() 
   << "AnbJet px : " << AnbJet->Px() << " py : " << AnbJet->Py() << " pz : " <<  AnbJet->Pz() << " energy : " << AnbJet->Energy() << endl;*/


   // Apply Top 
   (*Top)   = (*AnLep) + (*bJet)   + (*Nu);
   (*AnTop) = (*Lep)   + (*AnbJet) + (*AnNu);

}
void ssb_analysis::SetUpKINObs()
{
   isKinSol=false;
   v_leptons_VLV.clear();
   v_jets_VLV.clear();
   v_bjets_VLV.clear();
   v_lepidx_KIN.clear();
   v_anlepidx_KIN.clear();
   v_jetidx_KIN.clear();
   v_bjetidx_KIN.clear();
   v_btagging_KIN.clear();
   /// lepton ///
   v_leptons_VLV.push_back(common::TLVtoLV(*Lep));
   v_lepidx_KIN.push_back(0);
   v_leptons_VLV.push_back(common::TLVtoLV(*AnLep));
   v_anlepidx_KIN.push_back(1);
   /*cout << "jet info. for kinematic solver"  << endl;
   cout << v_jet_idx.size() << " " << v_jet_TL.size() << endl;*/
   const KinematicReconstruction* kinematicReconstruction(0); 
   kinematicReconstruction = new KinematicReconstruction(1, true);

   const LV met_LV = common::TLVtoLV(*Met);

   for (int ijet = 0; ijet < v_jet_idx.size(); ++ijet)
   {
      int idx_jet = v_jet_idx[ijet];
      if(Jet_bDisc->at(idx_jet)  > bdisccut ){
         v_bjets_VLV.push_back(common::TLVtoLV(*v_jet_TL[ijet]));
         //v_bjetidx_KIN.push_back(v_jet_idx[ijet]);
         v_bjetidx_KIN.push_back(ijet);
      }
      v_jets_VLV.push_back(common::TLVtoLV(*v_jet_TL[ijet]));
      v_jetidx_KIN.push_back(ijet);
      v_btagging_KIN.push_back(Jet_bDisc->at(idx_jet));
   }
   KinematicReconstructionSolutions kinematicReconstructionSolutions = kinematicReconstruction->solutions(v_lepidx_KIN,v_anlepidx_KIN, v_jetidx_KIN,  v_bjetidx_KIN,  v_leptons_VLV, v_jets_VLV, v_btagging_KIN, met_LV);
   //cout << "Num Sol : " << kinematicReconstructionSolutions.numberOfSolutions() << endl;
   //cout << "MET ? " << met_LV.pt() << endl;
   if (kinematicReconstructionSolutions.numberOfSolutions())
   {
      isKinSol= true;
      LV top1 = kinematicReconstructionSolutions.solution().top();
      LV top2 = kinematicReconstructionSolutions.solution().antiTop();
      LV bjet1 = kinematicReconstructionSolutions.solution().bjet();
      LV bjet2 = kinematicReconstructionSolutions.solution().antiBjet();
      LV neutrino1 = kinematicReconstructionSolutions.solution().neutrino();
      LV neutrino2 = kinematicReconstructionSolutions.solution().antiNeutrino();
      //kinematicReconstructionSolutions.solution().print();
      //Top = new TLorentzVector(common::LVtoTLV(top1));   
      Top       = new TLorentzVector(common::LVtoTLV(top1));
      AnTop     = new TLorentzVector(common::LVtoTLV(top2));
      bJet      = new TLorentzVector(common::LVtoTLV(bjet1));
      AnbJet    = new TLorentzVector(common::LVtoTLV(bjet2));
      Nu        = new TLorentzVector(common::LVtoTLV(neutrino1));
      AnNu      = new TLorentzVector(common::LVtoTLV(neutrino2));

      (*W1)        = (*Lep)  + (*AnNu);
      (*W2)        = (*AnLep)  + (*Nu);

   }
   //delete 
   delete kinematicReconstruction;
}
//void ssb_analysis::SetUpKINObsSyst()
void ssb_analysis::SetUpKINObsSyst(vector<int> v_jetsys_idx, vector<TLorentzVector*> v_jetsys_TL, TLorentzVector* metsys)
{
   isKinSol=false;
   v_leptons_VLV.clear();
   v_jets_VLV.clear();
   v_bjets_VLV.clear();
   v_lepidx_KIN.clear();
   v_anlepidx_KIN.clear();
   v_jetidx_KIN.clear();
   v_bjetidx_KIN.clear();
   v_btagging_KIN.clear();
   /// lepton ///
   v_leptons_VLV.push_back(common::TLVtoLV(*Lep));
   v_lepidx_KIN.push_back(0);
   v_leptons_VLV.push_back(common::TLVtoLV(*AnLep));
   v_anlepidx_KIN.push_back(1);
   //cout << "jet info. for kinematic solver"  << endl;
   //cout << v_jet_idx.size() << " " << v_jetsys_TL.size() << endl;
   const KinematicReconstruction* kinematicReconstruction(0); 
   kinematicReconstruction = new KinematicReconstruction(1, true);

   const LV met_LV = common::TLVtoLV(*metsys);

   for (int ijet = 0; ijet < v_jetsys_idx.size(); ++ijet)
   {
      int idx_jet = v_jetsys_idx[ijet];
      if(Jet_bDisc->at(idx_jet)  > bdisccut ){
         v_bjets_VLV.push_back(common::TLVtoLV(*v_jetsys_TL[ijet]));
         //v_bjetidx_KIN.push_back(v_jetsys_idx[ijet]);
         v_bjetidx_KIN.push_back(ijet);
      }
      v_jets_VLV.push_back(common::TLVtoLV(*v_jetsys_TL[ijet]));
      v_jetidx_KIN.push_back(ijet);
      v_btagging_KIN.push_back(Jet_bDisc->at(idx_jet));
   }
   KinematicReconstructionSolutions kinematicReconstructionSolutions = kinematicReconstruction->solutions(v_lepidx_KIN,v_anlepidx_KIN, v_jetidx_KIN,  v_bjetidx_KIN,  v_leptons_VLV, v_jets_VLV, v_btagging_KIN, met_LV);
   if (kinematicReconstructionSolutions.numberOfSolutions())
   {
      isKinSol= true;
      LV top1 = kinematicReconstructionSolutions.solution().top();
      LV top2 = kinematicReconstructionSolutions.solution().antiTop();
      LV bjet1 = kinematicReconstructionSolutions.solution().bjet();
      LV bjet2 = kinematicReconstructionSolutions.solution().antiBjet();
      LV neutrino1 = kinematicReconstructionSolutions.solution().neutrino();
      LV neutrino2 = kinematicReconstructionSolutions.solution().antiNeutrino();
      //Top = new TLorentzVector(common::LVtoTLV(top1));   
      Top       = new TLorentzVector(common::LVtoTLV(top1));
      AnTop     = new TLorentzVector(common::LVtoTLV(top2));
      bJet      = new TLorentzVector(common::LVtoTLV(bjet1));
      AnbJet    = new TLorentzVector(common::LVtoTLV(bjet2));
      Nu        = new TLorentzVector(common::LVtoTLV(neutrino1));
      AnNu      = new TLorentzVector(common::LVtoTLV(neutrino2));

      (*W1)        = (*Lep)  + (*AnNu);
      (*W2)        = (*AnLep)  + (*Nu);
      //cout << "Top pt " << Top->Pt() << endl;
   }
   delete kinematicReconstruction;
}
