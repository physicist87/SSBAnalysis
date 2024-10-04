//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 16 11:30:18 2024 by ROOT version 6.14/09
// from TTree SSBTree/Tree for Physics Analyses at CMS
// found on file: dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/Run2FULL/2016PreVFP/MC/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_20240205_143158/240205_053742/0000/SSBTree_108.root
//////////////////////////////////////////////////////////

#ifndef SSBTree_h
#define SSBTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "TClonesArray.h"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

using namespace std;

class SSBTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Info_EventNumber;
   Int_t           Info_Luminosity;
   Int_t           Info_RunNumber;
   //Char_t          Info_isData;
   Bool_t          Info_isData;
   Int_t           Channel_Idx;
   Int_t           Channel_Idx_Final;
   Int_t           Channel_Lepton_Count;
   Int_t           Channel_Lepton_Count_Final;
   Int_t           Channel_Jets;
   Int_t           Channel_Jets_Abs;
   Double_t        L1_PreFire_Central;
   Double_t        L1_PreFire_Up;
   Double_t        L1_PreFire_Down;
   vector<string>  *METFilter_Name;
   vector<bool>    *METFilter_isError;
   vector<bool>    *METFilter_isPass;
   vector<bool>    *METFilter_isRun;
   vector<string>  *METFilterAdd_Name;
   vector<bool>    *METFilterAdd_isPass;
   TClonesArray    *Elec;
   TClonesArray    *RawElec;
   vector<int>     *Elec_Charge;
   vector<bool>    *Elec_ChargeId_GsfCtf;
   vector<bool>    *Elec_ChargeId_GsfCtfPx;
   vector<bool>    *Elec_ChargeId_GsfPx;
   vector<double>  *Elec_Charge_CtfTr;
   vector<double>  *Elec_Charge_GsfTr;
   vector<double>  *Elec_Charge_Px;
   vector<bool>    *Elec_Conversion;
   Int_t           Elec_Count;
   vector<int>     *Elec_Id_Loose;
   vector<int>     *Elec_Id_RobustHighEnergy;
   vector<int>     *Elec_Id_RobustLoose;
   vector<int>     *Elec_Id_RobustTight;
   vector<int>     *Elec_Id_Tight;
   vector<int>     *Elec_Inner_Hit;
   vector<bool>    *Elec_MVA_Loose;
   vector<bool>    *Elec_MVA_Medium;
   vector<bool>    *Elec_MVA_Tight;
   vector<float>   *Elec_MVA_Values;
   vector<bool>    *Elec_MVA_NonIso_Loose;
   vector<bool>    *Elec_MVA_NonIso_Medium;
   vector<bool>    *Elec_MVA_NonIso_Tight;
   vector<float>   *Elec_MVA_Iso_V;
   vector<float>   *Elec_MVA_nonIso_V;
   vector<int>     *Elec_MVA_Categories;
   vector<float>   *Elec_MVA_HZZ_Values;
   vector<int>     *Elec_MVA_HZZ_Categories;
   vector<double>  *Elec_PFIsoRho03;
   vector<double>  *Elec_PFIsoRho04;
   vector<bool>    *Elec_PFIsoValid;
   vector<double>  *Elec_PFIsodBeta03;
   vector<double>  *Elec_PFIsodBeta04;
   vector<bool>    *Elec_SCB_Loose;
   vector<bool>    *Elec_SCB_Medium;
   vector<bool>    *Elec_SCB_Tight;
   vector<bool>    *Elec_SCB_Veto;
   vector<float>   *Elec_SCB_dEtaIn;
   vector<float>   *Elec_SCB_dPhiIn;
   vector<float>   *Elec_SCB_hOverE;
   vector<bool>    *Elec_SCB_HEEP;
   vector<float>   *Elec_SCB_ooEmooP;
   vector<float>   *Elec_SCB_sigmaIetaIeta;
   vector<float>   *Elec_Scale_StatUp;
   vector<float>   *Elec_Scale_StatDown;
   vector<float>   *Elec_Scale_SystUp;
   vector<float>   *Elec_Scale_SystDown;
   vector<float>   *Elec_GainUp;
   vector<float>   *Elec_GainDown;
   vector<float>   *Elec_RhoUp;
   vector<float>   *Elec_RhoDown;
   vector<float>   *Elec_PhiUp;
   vector<float>   *Elec_PhiDown;
   vector<double>  *Elec_Supercluster_Eta;
   vector<double>  *Elec_Track_CtfdXY;
   vector<double>  *Elec_Track_CtfdZ;
   vector<double>  *Elec_Track_GsfdXY;
   vector<double>  *Elec_Track_GsfdZ;
   vector<int>     *Elec_pdgId;
   vector<double>  *Elec_relIso03;
   vector<double>  *Elec_relIso04;
   vector<double>  *Elec_Rnd;
   Char_t          Filter_Greedy_Muon;
   Char_t          Filter_Inconsistent_MuonPt;
   Char_t          Filter_PFReco_Muon;
   vector<bool>    *Filter_PV;
   TClonesArray    *GenJet;
   Int_t           GenJet_Count;
   vector<float>   *GenJet_ECalEnergy;
   vector<float>   *GenJet_HCalEnergy;
   TClonesArray    *GenMET;
   Int_t           GenMET_Count;
   TClonesArray    *GenPar;
   Int_t           GenPar_Count;
   vector<int>     *GenPar_Dau1_Idx;
   vector<int>     *GenPar_Dau2_Idx;
   vector<int>     *GenPar_Dau_Counter;
   vector<int>     *GenPar_Idx;
   vector<int>     *GenPar_Mom1_Idx;
   vector<int>     *GenPar_Mom2_Idx;
   vector<int>     *GenPar_Mom_Counter;
   vector<int>     *GenPar_Status;
   vector<int>     *GenPar_pdgId;
   Double_t        Gen_EventWeight;
   TClonesArray    *GenTop;
   TClonesArray    *GenAnTop;
   TClonesArray    *GenBJet;
   Int_t           GenBJet_Count;
   TClonesArray    *GenBHad;
   Int_t           GenBHad_Count;
   vector<int>     *GenBHad_pdgId;
   vector<int>     *GenBHad_Flavour;
   vector<int>     *GenBHad_FromTopWeakDecay;
   TClonesArray    *Jet;
   TClonesArray    *RawJet;
   vector<int>     *Jet_Charge;
   Int_t           Jet_Count;
   vector<double>  *Jet_EnShiftedDown;
   vector<double>  *Jet_EnShiftedUp;
   vector<int>     *Jet_PartonFlavour;
   vector<int>     *Jet_HadronFlavour;
   vector<int>     *Jet_PFId;
   vector<int>     *Jet_PFIdVeto;
   vector<float>   *Jet_PhiResolution_MC;
   vector<float>   *Jet_PhiResolution_DATA;
   vector<double>  *Jet_EnergyResolution_MC;
   vector<double>  *Jet_EnergyResolution_DATA;
   vector<double>  *Jet_EnergyResolution_SF;
   vector<double>  *Jet_EnergyResolution_SFDown;
   vector<double>  *Jet_EnergyResolution_SFUp;
   vector<int>     *Jet_PileUpId;
   vector<float>   *Jet_PileUpMVA;
   vector<float>   *Jet_bDisc;
   vector<string>  *Jet_bDisc_Name;
   vector<float>   *Jet_bDisc_Value;
   vector<bool>    *Jet_isJet;
   Double_t        LHE_Central;
   vector<int>     *LHE_Id;
   vector<double>  *LHE_Weight;
   Double_t        Frag_Cen_Weight;
   Double_t        Frag_Up_Weight;
   Double_t        Frag_Down_Weight;
   Double_t        Frag_Peterson_Weight;
   Double_t        Semilep_BrUp_Weight;
   Double_t        Semilep_BrDown_Weight;
   TClonesArray    *MET;
   Double_t        MET_Significance;
   TClonesArray    *Muon;
   TClonesArray    *GenMuon;
   vector<int>     *Muon_Charge;
   Int_t           Muon_Count;
   vector<double>  *Muon_PFIsodBeta03;
   vector<double>  *Muon_PFIsodBeta04;
   vector<bool>    *Muon_isHighPt;
   vector<bool>    *Muon_isHighTrkPt;
   vector<bool>    *Muon_isLoose;
   vector<bool>    *Muon_isMedium;
   vector<bool>    *Muon_isMediumPrompt;
   vector<bool>    *Muon_isSoft;
   vector<bool>    *Muon_isTight;
   vector<bool>    *Muon_isPFIsoVeryLoose;
   vector<bool>    *Muon_isPFIsoLoose;
   vector<bool>    *Muon_isPFIsoMedium;
   vector<bool>    *Muon_isPFIsoTight;
   vector<bool>    *Muon_isPFIsoVeryTight;
   vector<bool>    *Muon_isPFIsoVeryVeryTight;
   vector<int>     *Muon_pdgId;
   vector<double>  *Muon_rand1;
   vector<double>  *Muon_rand2;
   vector<double>  *Muon_relIso03;
   vector<double>  *Muon_relIso04;
   vector<int>     *Muon_trackerLayers;
   vector<double>  *Muon_tuneP_Pt;
   vector<double>  *Muon_tuneP_Eta;
   vector<double>  *Muon_tuneP_Phi;
   vector<int>     *Muon_tuneP_Charge;
   vector<double>  *PDFWeight_BjorkenX1;
   vector<double>  *PDFWeight_BjorkenX2;
   vector<double>  *PDFWeight_Cent;
   vector<double>  *PDFWeight_Cent_Down;
   vector<double>  *PDFWeight_Cent_Up;
   vector<int>     *PDFWeight_Id1;
   vector<int>     *PDFWeight_Id2;
   vector<double>  *PDFWeight_PDF1;
   vector<double>  *PDFWeight_PDF2;
   vector<double>  *PDFWeight_Q;
   vector<double>  *PDFWeight_Var1_Down;
   vector<double>  *PDFWeight_Var1_Up;
   vector<double>  *PDFWeight_Var2_Down;
   vector<double>  *PDFWeight_Var2_Up;
   TClonesArray    *Photon;
   Int_t           Photon_Count;
   vector<bool>    *Photon_SCB_Loose;
   vector<bool>    *Photon_SCB_Medium;
   vector<bool>    *Photon_SCB_Tight;
   vector<float>   *Photon_ChaHadIso;
   vector<float>   *Photon_NeuHadIso;
   vector<float>   *Photon_PhoIso;
   vector<float>   *Photon_WorstChargedIso;
   vector<float>   *Photon_ChaHadEffArea;
   vector<float>   *Photon_NeuHadEffArea;
   vector<float>   *Photon_PhoHadEffArea;
   vector<double>  *Photon_R9;
   vector<double>  *Photon_HoverE;
   vector<double>  *Photon_SuperCluster_Eta;
   vector<double>  *Photon_SuperCluster_Phi;
   vector<double>  *Photon_SuperCluster_EtaWidth;
   vector<double>  *Photon_SuperCluster_PhiWidth;
   vector<double>  *Photon_SuperCluster_Energy;
   vector<double>  *Photon_SuperCluster_RawEnergy;
   vector<bool>    *Photon_ElectronVeto;
   vector<double>  *Photon_Full5x5_SigmaIetaIeta;
   vector<bool>    *Photon_MVANonTrig_Tight;
   Int_t           PV_Count;
   Int_t           PileUp_Count_Interaction;
   Float_t         PileUp_Count_Intime;
   vector<double>  *Rho;
   vector<string>  *Trigger_Name;
   vector<unsigned int> *Trigger_PreScale;
   vector<bool>    *Trigger_isError;
   vector<bool>    *Trigger_isPass;
   vector<bool>    *Trigger_isRun;
   vector<double>  *Vertex_SumPtSquare;
   vector<double>  *Vertex_X;
   vector<double>  *Vertex_X_Error;
   vector<double>  *Vertex_Y;
   vector<double>  *Vertex_Y_Error;
   vector<double>  *Vertex_Z;
   vector<double>  *Vertex_Z_Error;

   // List of branches
   TBranch        *b_Info_EventNumber;   //!
   TBranch        *b_Info_Luminosity;   //!
   TBranch        *b_Info_RunNumber;   //!
   TBranch        *b_Info_isData;   //!
   TBranch        *b_Channel_Idx;   //!
   TBranch        *b_Channel_Idx_Final;   //!
   TBranch        *b_Channel_Lepton_Count;   //!
   TBranch        *b_Channel_Lepton_Count_Final;   //!
   TBranch        *b_Channel_Jets;   //!
   TBranch        *b_Channel_Jets_Abs;   //!
   TBranch        *b_L1_PreFire_Central;   //!
   TBranch        *b_L1_PreFire_Up;   //!
   TBranch        *b_L1_PreFire_Down;   //!
   TBranch        *b_METFilter_Name;   //!
   TBranch        *b_METFilter_isError;   //!
   TBranch        *b_METFilter_isPass;   //!
   TBranch        *b_METFilter_isRun;   //!
   TBranch        *b_METFilterAdd_Name;   //!
   TBranch        *b_METFilterAdd_isPass;   //!
   TBranch        *b_Elec;   //!
   TBranch        *b_RawElec;   //!
   TBranch        *b_Elec_Charge;   //!
   TBranch        *b_Elec_ChargeId_GsfCtf;   //!
   TBranch        *b_Elec_ChargeId_GsfCtfPx;   //!
   TBranch        *b_Elec_ChargeId_GsfPx;   //!
   TBranch        *b_Elec_Charge_CtfTr;   //!
   TBranch        *b_Elec_Charge_GsfTr;   //!
   TBranch        *b_Elec_Charge_Px;   //!
   TBranch        *b_Elec_Conversion;   //!
   TBranch        *b_Elec_Count;   //!
   TBranch        *b_Elec_Id_Loose;   //!
   TBranch        *b_Elec_Id_RobustHighEnergy;   //!
   TBranch        *b_Elec_Id_RobustLoose;   //!
   TBranch        *b_Elec_Id_RobustTight;   //!
   TBranch        *b_Elec_Id_Tight;   //!
   TBranch        *b_Elec_Inner_Hit;   //!
   TBranch        *b_Elec_MVA_Loose;   //!
   TBranch        *b_Elec_MVA_Medium;   //!
   TBranch        *b_Elec_MVA_Tight;   //!
   TBranch        *b_Elec_MVA_Values;   //!
   TBranch        *b_Elec_MVA_NonIso_Loose;   //!
   TBranch        *b_Elec_MVA_NonIso_Medium;   //!
   TBranch        *b_Elec_MVA_NonIso_Tight;   //!
   TBranch        *b_Elec_MVA_Iso_V;   //!
   TBranch        *b_Elec_MVA_nonIso_V;   //!
   TBranch        *b_Elec_MVA_Categories;   //!
   TBranch        *b_Elec_MVA_HZZ_Values;   //!
   TBranch        *b_Elec_MVA_HZZ_Categories;   //!
   TBranch        *b_Elec_PFIsoRho03;   //!
   TBranch        *b_Elec_PFIsoRho04;   //!
   TBranch        *b_Elec_PFIsoValid;   //!
   TBranch        *b_Elec_PFIsodBeta03;   //!
   TBranch        *b_Elec_PFIsodBeta04;   //!
   TBranch        *b_Elec_SCB_Loose;   //!
   TBranch        *b_Elec_SCB_Medium;   //!
   TBranch        *b_Elec_SCB_Tight;   //!
   TBranch        *b_Elec_SCB_Veto;   //!
   TBranch        *b_Elec_SCB_dEtaIn;   //!
   TBranch        *b_Elec_SCB_dPhiIn;   //!
   TBranch        *b_Elec_SCB_hOverE;   //!
   TBranch        *b_Elec_SCB_HEEP;   //!
   TBranch        *b_Elec_SCB_ooEmooP;   //!
   TBranch        *b_Elec_SCB_sigmaIetaIeta;   //!
   TBranch        *b_Elec_Scale_StatUp;   //!
   TBranch        *b_Elec_Scale_StatDown;   //!
   TBranch        *b_Elec_Scale_SystUp;   //!
   TBranch        *b_Elec_Scale_SystDown;   //!
   TBranch        *b_Elec_GainUp;   //!
   TBranch        *b_Elec_GainDown;   //!
   TBranch        *b_Elec_RhoUp;   //!
   TBranch        *b_Elec_RhoDown;   //!
   TBranch        *b_Elec_PhiUp;   //!
   TBranch        *b_Elec_PhiDown;   //!
   TBranch        *b_Elec_Supercluster_Eta;   //!
   TBranch        *b_Elec_Track_CtfdXY;   //!
   TBranch        *b_Elec_Track_CtfdZ;   //!
   TBranch        *b_Elec_Track_GsfdXY;   //!
   TBranch        *b_Elec_Track_GsfdZ;   //!
   TBranch        *b_Elec_pdgId;   //!
   TBranch        *b_Elec_relIso03;   //!
   TBranch        *b_Elec_relIso04;   //!
   TBranch        *b_Elec_Rnd;   //!
   TBranch        *b_Filter_Greedy_Muon;   //!
   TBranch        *b_Filter_Inconsistent_MuonPt;   //!
   TBranch        *b_Filter_PFReco_Muon;   //!
   TBranch        *b_Filter_PV;   //!
   TBranch        *b_GenJet;   //!
   TBranch        *b_GenJet_Count;   //!
   TBranch        *b_GenJet_ECalEnergy;   //!
   TBranch        *b_GenJet_HCalEnergy;   //!
   TBranch        *b_GenMET;   //!
   TBranch        *b_GenMET_Count;   //!
   TBranch        *b_GenPar;   //!
   TBranch        *b_GenPar_Count;   //!
   TBranch        *b_GenPar_Dau1_Idx;   //!
   TBranch        *b_GenPar_Dau2_Idx;   //!
   TBranch        *b_GenPar_Dau_Counter;   //!
   TBranch        *b_GenPar_Idx;   //!
   TBranch        *b_GenPar_Mom1_Idx;   //!
   TBranch        *b_GenPar_Mom2_Idx;   //!
   TBranch        *b_GenPar_Mom_Counter;   //!
   TBranch        *b_GenPar_Status;   //!
   TBranch        *b_GenPar_pdgId;   //!
   TBranch        *b_Gen_EventWeight;   //!
   TBranch        *b_GenTop;   //!
   TBranch        *b_GenAnTop;   //!
   TBranch        *b_GenBJet;   //!
   TBranch        *b_GenBJet_Count;   //!
   TBranch        *b_GenBHad;   //!
   TBranch        *b_GenBHad_Count;   //!
   TBranch        *b_GenBHad_pdgId;   //!
   TBranch        *b_GenBHad_Flavour;   //!
   TBranch        *b_GenBHad_FromTopWeakDecay;   //!
   TBranch        *b_Jet;   //!
   TBranch        *b_RawJet;   //!
   TBranch        *b_Jet_Charge;   //!
   TBranch        *b_Jet_Count;   //!
   TBranch        *b_Jet_EnShiftedDown;   //!
   TBranch        *b_Jet_EnShiftedUp;   //!
   TBranch        *b_Jet_PartonFlavour;   //!
   TBranch        *b_Jet_HadronFlavour;   //!
   TBranch        *b_Jet_PFId;   //!
   TBranch        *b_Jet_PFIdVeto;   //!
   TBranch        *b_Jet_PhiResolution_MC;   //!
   TBranch        *b_Jet_PhiResolution_DATA;   //!
   TBranch        *b_Jet_EnergyResolution_MC;   //!
   TBranch        *b_Jet_EnergyResolution_DATA;   //!
   TBranch        *b_Jet_EnergyResolution_SF;   //!
   TBranch        *b_Jet_EnergyResolution_SFDown;   //!
   TBranch        *b_Jet_EnergyResolution_SFUp;   //!
   TBranch        *b_Jet_PileUpId;   //!
   TBranch        *b_Jet_PileUpMVA;   //!
   TBranch        *b_Jet_bDisc;   //!
   TBranch        *b_Jet_bDisc_Name;   //!
   TBranch        *b_Jet_bDisc_Value;   //!
   TBranch        *b_Jet_isJet;   //!
   TBranch        *b_LHE_Central;   //!
   TBranch        *b_LHE_Id;   //!
   TBranch        *b_LHE_Weight;   //!
   TBranch        *b_Frag_Cen_Weight;   //!
   TBranch        *b_Frag_Up_Weight;   //!
   TBranch        *b_Frag_Down_Weight;   //!
   TBranch        *b_Frag_Peterson_Weight;   //!
   TBranch        *b_Semilep_BrUp_Weight;   //!
   TBranch        *b_Semilep_BrDown_Weight;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_MET_Significance;   //!
   TBranch        *b_Muon;   //!
   TBranch        *b_GenMuon;   //!
   TBranch        *b_Muon_Charge;   //!
   TBranch        *b_Muon_Count;   //!
   TBranch        *b_Muon_PFIsodBeta03;   //!
   TBranch        *b_Muon_PFIsodBeta04;   //!
   TBranch        *b_Muon_isHighPt;   //!
   TBranch        *b_Muon_isHighTrkPt;   //!
   TBranch        *b_Muon_isLoose;   //!
   TBranch        *b_Muon_isMedium;   //!
   TBranch        *b_Muon_isMediumPrompt;   //!
   TBranch        *b_Muon_isSoft;   //!
   TBranch        *b_Muon_isTight;   //!
   TBranch        *b_Muon_isPFIsoVeryLoose;   //!
   TBranch        *b_Muon_isPFIsoLoose;   //!
   TBranch        *b_Muon_isPFIsoMedium;   //!
   TBranch        *b_Muon_isPFIsoTight;   //!
   TBranch        *b_Muon_isPFIsoVeryTight;   //!
   TBranch        *b_Muon_isPFIsoVeryVeryTight;   //!
   TBranch        *b_Muon_pdgId;   //!
   TBranch        *b_Muon_rand1;   //!
   TBranch        *b_Muon_rand2;   //!
   TBranch        *b_Muon_relIso03;   //!
   TBranch        *b_Muon_relIso04;   //!
   TBranch        *b_Muon_trackerLayers;   //!
   TBranch        *b_Muon_tuneP_Pt;   //!
   TBranch        *b_Muon_tuneP_Eta;   //!
   TBranch        *b_Muon_tuneP_Phi;   //!
   TBranch        *b_Muon_tuneP_Charge;   //!
   TBranch        *b_PDFWeight_BjorkenX1;   //!
   TBranch        *b_PDFWeight_BjorkenX2;   //!
   TBranch        *b_PDFWeight_Cent;   //!
   TBranch        *b_PDFWeight_Cent_Down;   //!
   TBranch        *b_PDFWeight_Cent_Up;   //!
   TBranch        *b_PDFWeight_Id1;   //!
   TBranch        *b_PDFWeight_Id2;   //!
   TBranch        *b_PDFWeight_PDF1;   //!
   TBranch        *b_PDFWeight_PDF2;   //!
   TBranch        *b_PDFWeight_Q;   //!
   TBranch        *b_PDFWeight_Var1_Down;   //!
   TBranch        *b_PDFWeight_Var1_Up;   //!
   TBranch        *b_PDFWeight_Var2_Down;   //!
   TBranch        *b_PDFWeight_Var2_Up;   //!
   TBranch        *b_Photon;   //!
   TBranch        *b_Photon_Count;   //!
   TBranch        *b_Photon_SCB_Loose;   //!
   TBranch        *b_Photon_SCB_Medium;   //!
   TBranch        *b_Photon_SCB_Tight;   //!
   TBranch        *b_Photon_ChaHadIso;   //!
   TBranch        *b_Photon_NeuHadIso;   //!
   TBranch        *b_Photon_PhoIso;   //!
   TBranch        *b_Photon_WorstChargedIso;   //!
   TBranch        *b_Photon_ChaHadEffArea;   //!
   TBranch        *b_Photon_NeuHadEffArea;   //!
   TBranch        *b_Photon_PhoHadEffArea;   //!
   TBranch        *b_Photon_R9;   //!
   TBranch        *b_Photon_HoverE;   //!
   TBranch        *b_Photon_SuperCluster_Eta;   //!
   TBranch        *b_Photon_SuperCluster_Phi;   //!
   TBranch        *b_Photon_SuperCluster_EtaWidth;   //!
   TBranch        *b_Photon_SuperCluster_PhiWidth;   //!
   TBranch        *b_Photon_SuperCluster_Energy;   //!
   TBranch        *b_Photon_SuperCluster_RawEnergy;   //!
   TBranch        *b_Photon_ElectronVeto;   //!
   TBranch        *b_Photon_Full5x5_SigmaIetaIeta;   //!
   TBranch        *b_Photon_MVANonTrig_Tight;   //!
   TBranch        *b_PV_Count;   //!
   TBranch        *b_PileUp_Count_Interaction;   //!
   TBranch        *b_PileUp_Count_Intime;   //!
   TBranch        *b_Rho;   //!
   TBranch        *b_Trigger_Name;   //!
   TBranch        *b_Trigger_PreScale;   //!
   TBranch        *b_Trigger_isError;   //!
   TBranch        *b_Trigger_isPass;   //!
   TBranch        *b_Trigger_isRun;   //!
   TBranch        *b_Vertex_SumPtSquare;   //!
   TBranch        *b_Vertex_X;   //!
   TBranch        *b_Vertex_X_Error;   //!
   TBranch        *b_Vertex_Y;   //!
   TBranch        *b_Vertex_Y_Error;   //!
   TBranch        *b_Vertex_Z;   //!
   TBranch        *b_Vertex_Z_Error;   //!

   SSBTree(TTree *tree=0);
   virtual ~SSBTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef SSBTree_cxx
SSBTree::SSBTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {/*
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/Run2FULL/2016PreVFP/MC/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_20240205_143158/240205_053742/0000/SSBTree_108.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/Run2FULL/2016PreVFP/MC/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_20240205_143158/240205_053742/0000/SSBTree_108.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/Run2FULL/2016PreVFP/MC/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_20240205_143158/240205_053742/0000/SSBTree_108.root:/ssbanalyzer");
      dir->GetObject("SSBTree",tree);

   */}
   Init(tree);
}

SSBTree::~SSBTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SSBTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SSBTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void SSBTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   METFilter_Name = 0;
   METFilter_isError = 0;
   METFilter_isPass = 0;
   METFilter_isRun = 0;
   METFilterAdd_Name = 0;
   METFilterAdd_isPass = 0;
   Elec = 0;
   RawElec = 0;
   Elec_Charge = 0;
   Elec_ChargeId_GsfCtf = 0;
   Elec_ChargeId_GsfCtfPx = 0;
   Elec_ChargeId_GsfPx = 0;
   Elec_Charge_CtfTr = 0;
   Elec_Charge_GsfTr = 0;
   Elec_Charge_Px = 0;
   Elec_Conversion = 0;
   Elec_Id_Loose = 0;
   Elec_Id_RobustHighEnergy = 0;
   Elec_Id_RobustLoose = 0;
   Elec_Id_RobustTight = 0;
   Elec_Id_Tight = 0;
   Elec_Inner_Hit = 0;
   Elec_MVA_Loose = 0;
   Elec_MVA_Medium = 0;
   Elec_MVA_Tight = 0;
   Elec_MVA_Values = 0;
   Elec_MVA_NonIso_Loose = 0;
   Elec_MVA_NonIso_Medium = 0;
   Elec_MVA_NonIso_Tight = 0;
   Elec_MVA_Iso_V = 0;
   Elec_MVA_nonIso_V = 0;
   Elec_MVA_Categories = 0;
   Elec_MVA_HZZ_Values = 0;
   Elec_MVA_HZZ_Categories = 0;
   Elec_PFIsoRho03 = 0;
   Elec_PFIsoRho04 = 0;
   Elec_PFIsoValid = 0;
   Elec_PFIsodBeta03 = 0;
   Elec_PFIsodBeta04 = 0;
   Elec_SCB_Loose = 0;
   Elec_SCB_Medium = 0;
   Elec_SCB_Tight = 0;
   Elec_SCB_Veto = 0;
   Elec_SCB_dEtaIn = 0;
   Elec_SCB_dPhiIn = 0;
   Elec_SCB_hOverE = 0;
   Elec_SCB_HEEP = 0;
   Elec_SCB_ooEmooP = 0;
   Elec_SCB_sigmaIetaIeta = 0;
   Elec_Scale_StatUp = 0;
   Elec_Scale_StatDown = 0;
   Elec_Scale_SystUp = 0;
   Elec_Scale_SystDown = 0;
   Elec_GainUp = 0;
   Elec_GainDown = 0;
   Elec_RhoUp = 0;
   Elec_RhoDown = 0;
   Elec_PhiUp = 0;
   Elec_PhiDown = 0;
   Elec_Supercluster_Eta = 0;
   Elec_Track_CtfdXY = 0;
   Elec_Track_CtfdZ = 0;
   Elec_Track_GsfdXY = 0;
   Elec_Track_GsfdZ = 0;
   Elec_pdgId = 0;
   Elec_relIso03 = 0;
   Elec_relIso04 = 0;
   Elec_Rnd = 0;
   Filter_PV = 0;
   GenJet = 0;
   GenJet_ECalEnergy = 0;
   GenJet_HCalEnergy = 0;
   GenMET = 0;
   GenPar = 0;
   GenPar_Dau1_Idx = 0;
   GenPar_Dau2_Idx = 0;
   GenPar_Dau_Counter = 0;
   GenPar_Idx = 0;
   GenPar_Mom1_Idx = 0;
   GenPar_Mom2_Idx = 0;
   GenPar_Mom_Counter = 0;
   GenPar_Status = 0;
   GenPar_pdgId = 0;
   GenTop = 0;
   GenAnTop = 0;
   GenBJet = 0;
   GenBHad = 0;
   GenBHad_pdgId = 0;
   GenBHad_Flavour = 0;
   GenBHad_FromTopWeakDecay = 0;
   Jet = 0;
   RawJet = 0;
   Jet_Charge = 0;
   Jet_EnShiftedDown = 0;
   Jet_EnShiftedUp = 0;
   Jet_PartonFlavour = 0;
   Jet_HadronFlavour = 0;
   Jet_PFId = 0;
   Jet_PFIdVeto = 0;
   Jet_PhiResolution_MC = 0;
   Jet_PhiResolution_DATA = 0;
   Jet_EnergyResolution_MC = 0;
   Jet_EnergyResolution_DATA = 0;
   Jet_EnergyResolution_SF = 0;
   Jet_EnergyResolution_SFDown = 0;
   Jet_EnergyResolution_SFUp = 0;
   Jet_PileUpId = 0;
   Jet_PileUpMVA = 0;
   Jet_bDisc = 0;
   Jet_bDisc_Name = 0;
   Jet_bDisc_Value = 0;
   Jet_isJet = 0;
   LHE_Id = 0;
   LHE_Weight = 0;
   MET = 0;
   Muon = 0;
   GenMuon = 0;
   Muon_Charge = 0;
   Muon_PFIsodBeta03 = 0;
   Muon_PFIsodBeta04 = 0;
   Muon_isHighPt = 0;
   Muon_isHighTrkPt = 0;
   Muon_isLoose = 0;
   Muon_isMedium = 0;
   Muon_isMediumPrompt = 0;
   Muon_isSoft = 0;
   Muon_isTight = 0;
   Muon_isPFIsoVeryLoose = 0;
   Muon_isPFIsoLoose = 0;
   Muon_isPFIsoMedium = 0;
   Muon_isPFIsoTight = 0;
   Muon_isPFIsoVeryTight = 0;
   Muon_isPFIsoVeryVeryTight = 0;
   Muon_pdgId = 0;
   Muon_rand1 = 0;
   Muon_rand2 = 0;
   Muon_relIso03 = 0;
   Muon_relIso04 = 0;
   Muon_trackerLayers = 0;
   Muon_tuneP_Pt = 0;
   Muon_tuneP_Eta = 0;
   Muon_tuneP_Phi = 0;
   Muon_tuneP_Charge = 0;
   PDFWeight_BjorkenX1 = 0;
   PDFWeight_BjorkenX2 = 0;
   PDFWeight_Cent = 0;
   PDFWeight_Cent_Down = 0;
   PDFWeight_Cent_Up = 0;
   PDFWeight_Id1 = 0;
   PDFWeight_Id2 = 0;
   PDFWeight_PDF1 = 0;
   PDFWeight_PDF2 = 0;
   PDFWeight_Q = 0;
   PDFWeight_Var1_Down = 0;
   PDFWeight_Var1_Up = 0;
   PDFWeight_Var2_Down = 0;
   PDFWeight_Var2_Up = 0;
   Photon = 0;
   Photon_SCB_Loose = 0;
   Photon_SCB_Medium = 0;
   Photon_SCB_Tight = 0;
   Photon_ChaHadIso = 0;
   Photon_NeuHadIso = 0;
   Photon_PhoIso = 0;
   Photon_WorstChargedIso = 0;
   Photon_ChaHadEffArea = 0;
   Photon_NeuHadEffArea = 0;
   Photon_PhoHadEffArea = 0;
   Photon_R9 = 0;
   Photon_HoverE = 0;
   Photon_SuperCluster_Eta = 0;
   Photon_SuperCluster_Phi = 0;
   Photon_SuperCluster_EtaWidth = 0;
   Photon_SuperCluster_PhiWidth = 0;
   Photon_SuperCluster_Energy = 0;
   Photon_SuperCluster_RawEnergy = 0;
   Photon_ElectronVeto = 0;
   Photon_Full5x5_SigmaIetaIeta = 0;
   Photon_MVANonTrig_Tight = 0;
   Rho = 0;
   Trigger_Name = 0;
   Trigger_PreScale = 0;
   Trigger_isError = 0;
   Trigger_isPass = 0;
   Trigger_isRun = 0;
   Vertex_SumPtSquare = 0;
   Vertex_X = 0;
   Vertex_X_Error = 0;
   Vertex_Y = 0;
   Vertex_Y_Error = 0;
   Vertex_Z = 0;
   Vertex_Z_Error = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Info_EventNumber", &Info_EventNumber, &b_Info_EventNumber);
   fChain->SetBranchAddress("Info_Luminosity", &Info_Luminosity, &b_Info_Luminosity);
   fChain->SetBranchAddress("Info_RunNumber", &Info_RunNumber, &b_Info_RunNumber);
   fChain->SetBranchAddress("Info_isData", &Info_isData, &b_Info_isData);
   fChain->SetBranchAddress("Channel_Idx", &Channel_Idx, &b_Channel_Idx);
   fChain->SetBranchAddress("Channel_Idx_Final", &Channel_Idx_Final, &b_Channel_Idx_Final);
   fChain->SetBranchAddress("Channel_Lepton_Count", &Channel_Lepton_Count, &b_Channel_Lepton_Count);
   fChain->SetBranchAddress("Channel_Lepton_Count_Final", &Channel_Lepton_Count_Final, &b_Channel_Lepton_Count_Final);
   fChain->SetBranchAddress("Channel_Jets", &Channel_Jets, &b_Channel_Jets);
   fChain->SetBranchAddress("Channel_Jets_Abs", &Channel_Jets_Abs, &b_Channel_Jets_Abs);
   fChain->SetBranchAddress("L1_PreFire_Central", &L1_PreFire_Central, &b_L1_PreFire_Central);
   fChain->SetBranchAddress("L1_PreFire_Up", &L1_PreFire_Up, &b_L1_PreFire_Up);
   fChain->SetBranchAddress("L1_PreFire_Down", &L1_PreFire_Down, &b_L1_PreFire_Down);
   fChain->SetBranchAddress("METFilter_Name", &METFilter_Name, &b_METFilter_Name);
   fChain->SetBranchAddress("METFilter_isError", &METFilter_isError, &b_METFilter_isError);
   fChain->SetBranchAddress("METFilter_isPass", &METFilter_isPass, &b_METFilter_isPass);
   fChain->SetBranchAddress("METFilter_isRun", &METFilter_isRun, &b_METFilter_isRun);
   fChain->SetBranchAddress("METFilterAdd_Name", &METFilterAdd_Name, &b_METFilterAdd_Name);
   fChain->SetBranchAddress("METFilterAdd_isPass", &METFilterAdd_isPass, &b_METFilterAdd_isPass);
   fChain->SetBranchAddress("Elec", &Elec, &b_Elec);
   fChain->SetBranchAddress("RawElec", &RawElec, &b_RawElec);
   fChain->SetBranchAddress("Elec_Charge", &Elec_Charge, &b_Elec_Charge);
   fChain->SetBranchAddress("Elec_ChargeId_GsfCtf", &Elec_ChargeId_GsfCtf, &b_Elec_ChargeId_GsfCtf);
   fChain->SetBranchAddress("Elec_ChargeId_GsfCtfPx", &Elec_ChargeId_GsfCtfPx, &b_Elec_ChargeId_GsfCtfPx);
   fChain->SetBranchAddress("Elec_ChargeId_GsfPx", &Elec_ChargeId_GsfPx, &b_Elec_ChargeId_GsfPx);
   fChain->SetBranchAddress("Elec_Charge_CtfTr", &Elec_Charge_CtfTr, &b_Elec_Charge_CtfTr);
   fChain->SetBranchAddress("Elec_Charge_GsfTr", &Elec_Charge_GsfTr, &b_Elec_Charge_GsfTr);
   fChain->SetBranchAddress("Elec_Charge_Px", &Elec_Charge_Px, &b_Elec_Charge_Px);
   fChain->SetBranchAddress("Elec_Conversion", &Elec_Conversion, &b_Elec_Conversion);
   fChain->SetBranchAddress("Elec_Count", &Elec_Count, &b_Elec_Count);
   fChain->SetBranchAddress("Elec_Id_Loose", &Elec_Id_Loose, &b_Elec_Id_Loose);
   fChain->SetBranchAddress("Elec_Id_RobustHighEnergy", &Elec_Id_RobustHighEnergy, &b_Elec_Id_RobustHighEnergy);
   fChain->SetBranchAddress("Elec_Id_RobustLoose", &Elec_Id_RobustLoose, &b_Elec_Id_RobustLoose);
   fChain->SetBranchAddress("Elec_Id_RobustTight", &Elec_Id_RobustTight, &b_Elec_Id_RobustTight);
   fChain->SetBranchAddress("Elec_Id_Tight", &Elec_Id_Tight, &b_Elec_Id_Tight);
   fChain->SetBranchAddress("Elec_Inner_Hit", &Elec_Inner_Hit, &b_Elec_Inner_Hit);
   fChain->SetBranchAddress("Elec_MVA_Loose", &Elec_MVA_Loose, &b_Elec_MVA_Loose);
   fChain->SetBranchAddress("Elec_MVA_Medium", &Elec_MVA_Medium, &b_Elec_MVA_Medium);
   fChain->SetBranchAddress("Elec_MVA_Tight", &Elec_MVA_Tight, &b_Elec_MVA_Tight);
   fChain->SetBranchAddress("Elec_MVA_Values", &Elec_MVA_Values, &b_Elec_MVA_Values);
   fChain->SetBranchAddress("Elec_MVA_NonIso_Loose", &Elec_MVA_NonIso_Loose, &b_Elec_MVA_NonIso_Loose);
   fChain->SetBranchAddress("Elec_MVA_NonIso_Medium", &Elec_MVA_NonIso_Medium, &b_Elec_MVA_NonIso_Medium);
   fChain->SetBranchAddress("Elec_MVA_NonIso_Tight", &Elec_MVA_NonIso_Tight, &b_Elec_MVA_NonIso_Tight);
   fChain->SetBranchAddress("Elec_MVA_Iso_V", &Elec_MVA_Iso_V, &b_Elec_MVA_Iso_V);
   fChain->SetBranchAddress("Elec_MVA_nonIso_V", &Elec_MVA_nonIso_V, &b_Elec_MVA_nonIso_V);
   fChain->SetBranchAddress("Elec_MVA_Categories", &Elec_MVA_Categories, &b_Elec_MVA_Categories);
   fChain->SetBranchAddress("Elec_MVA_HZZ_Values", &Elec_MVA_HZZ_Values, &b_Elec_MVA_HZZ_Values);
   fChain->SetBranchAddress("Elec_MVA_HZZ_Categories", &Elec_MVA_HZZ_Categories, &b_Elec_MVA_HZZ_Categories);
   fChain->SetBranchAddress("Elec_PFIsoRho03", &Elec_PFIsoRho03, &b_Elec_PFIsoRho03);
   fChain->SetBranchAddress("Elec_PFIsoRho04", &Elec_PFIsoRho04, &b_Elec_PFIsoRho04);
   fChain->SetBranchAddress("Elec_PFIsoValid", &Elec_PFIsoValid, &b_Elec_PFIsoValid);
   fChain->SetBranchAddress("Elec_PFIsodBeta03", &Elec_PFIsodBeta03, &b_Elec_PFIsodBeta03);
   fChain->SetBranchAddress("Elec_PFIsodBeta04", &Elec_PFIsodBeta04, &b_Elec_PFIsodBeta04);
   fChain->SetBranchAddress("Elec_SCB_Loose", &Elec_SCB_Loose, &b_Elec_SCB_Loose);
   fChain->SetBranchAddress("Elec_SCB_Medium", &Elec_SCB_Medium, &b_Elec_SCB_Medium);
   fChain->SetBranchAddress("Elec_SCB_Tight", &Elec_SCB_Tight, &b_Elec_SCB_Tight);
   fChain->SetBranchAddress("Elec_SCB_Veto", &Elec_SCB_Veto, &b_Elec_SCB_Veto);
   fChain->SetBranchAddress("Elec_SCB_dEtaIn", &Elec_SCB_dEtaIn, &b_Elec_SCB_dEtaIn);
   fChain->SetBranchAddress("Elec_SCB_dPhiIn", &Elec_SCB_dPhiIn, &b_Elec_SCB_dPhiIn);
   fChain->SetBranchAddress("Elec_SCB_hOverE", &Elec_SCB_hOverE, &b_Elec_SCB_hOverE);
   fChain->SetBranchAddress("Elec_SCB_HEEP", &Elec_SCB_HEEP, &b_Elec_SCB_HEEP);
   fChain->SetBranchAddress("Elec_SCB_ooEmooP", &Elec_SCB_ooEmooP, &b_Elec_SCB_ooEmooP);
   fChain->SetBranchAddress("Elec_SCB_sigmaIetaIeta", &Elec_SCB_sigmaIetaIeta, &b_Elec_SCB_sigmaIetaIeta);
   fChain->SetBranchAddress("Elec_Scale_StatUp", &Elec_Scale_StatUp, &b_Elec_Scale_StatUp);
   fChain->SetBranchAddress("Elec_Scale_StatDown", &Elec_Scale_StatDown, &b_Elec_Scale_StatDown);
   fChain->SetBranchAddress("Elec_Scale_SystUp", &Elec_Scale_SystUp, &b_Elec_Scale_SystUp);
   fChain->SetBranchAddress("Elec_Scale_SystDown", &Elec_Scale_SystDown, &b_Elec_Scale_SystDown);
   fChain->SetBranchAddress("Elec_GainUp", &Elec_GainUp, &b_Elec_GainUp);
   fChain->SetBranchAddress("Elec_GainDown", &Elec_GainDown, &b_Elec_GainDown);
   fChain->SetBranchAddress("Elec_RhoUp", &Elec_RhoUp, &b_Elec_RhoUp);
   fChain->SetBranchAddress("Elec_RhoDown", &Elec_RhoDown, &b_Elec_RhoDown);
   fChain->SetBranchAddress("Elec_PhiUp", &Elec_PhiUp, &b_Elec_PhiUp);
   fChain->SetBranchAddress("Elec_PhiDown", &Elec_PhiDown, &b_Elec_PhiDown);
   fChain->SetBranchAddress("Elec_Supercluster_Eta", &Elec_Supercluster_Eta, &b_Elec_Supercluster_Eta);
   fChain->SetBranchAddress("Elec_Track_CtfdXY", &Elec_Track_CtfdXY, &b_Elec_Track_CtfdXY);
   fChain->SetBranchAddress("Elec_Track_CtfdZ", &Elec_Track_CtfdZ, &b_Elec_Track_CtfdZ);
   fChain->SetBranchAddress("Elec_Track_GsfdXY", &Elec_Track_GsfdXY, &b_Elec_Track_GsfdXY);
   fChain->SetBranchAddress("Elec_Track_GsfdZ", &Elec_Track_GsfdZ, &b_Elec_Track_GsfdZ);
   fChain->SetBranchAddress("Elec_pdgId", &Elec_pdgId, &b_Elec_pdgId);
   fChain->SetBranchAddress("Elec_relIso03", &Elec_relIso03, &b_Elec_relIso03);
   fChain->SetBranchAddress("Elec_relIso04", &Elec_relIso04, &b_Elec_relIso04);
   fChain->SetBranchAddress("Elec_Rnd", &Elec_Rnd, &b_Elec_Rnd);
   fChain->SetBranchAddress("Filter_Greedy_Muon", &Filter_Greedy_Muon, &b_Filter_Greedy_Muon);
   fChain->SetBranchAddress("Filter_Inconsistent_MuonPt", &Filter_Inconsistent_MuonPt, &b_Filter_Inconsistent_MuonPt);
   fChain->SetBranchAddress("Filter_PFReco_Muon", &Filter_PFReco_Muon, &b_Filter_PFReco_Muon);
   fChain->SetBranchAddress("Filter_PV", &Filter_PV, &b_Filter_PV);
   fChain->SetBranchAddress("GenJet", &GenJet, &b_GenJet);
   fChain->SetBranchAddress("GenJet_Count", &GenJet_Count, &b_GenJet_Count);
   fChain->SetBranchAddress("GenJet_ECalEnergy", &GenJet_ECalEnergy, &b_GenJet_ECalEnergy);
   fChain->SetBranchAddress("GenJet_HCalEnergy", &GenJet_HCalEnergy, &b_GenJet_HCalEnergy);
   fChain->SetBranchAddress("GenMET", &GenMET, &b_GenMET);
   fChain->SetBranchAddress("GenMET_Count", &GenMET_Count, &b_GenMET_Count);
   fChain->SetBranchAddress("GenPar", &GenPar, &b_GenPar);
   fChain->SetBranchAddress("GenPar_Count", &GenPar_Count, &b_GenPar_Count);
   fChain->SetBranchAddress("GenPar_Dau1_Idx", &GenPar_Dau1_Idx, &b_GenPar_Dau1_Idx);
   fChain->SetBranchAddress("GenPar_Dau2_Idx", &GenPar_Dau2_Idx, &b_GenPar_Dau2_Idx);
   fChain->SetBranchAddress("GenPar_Dau_Counter", &GenPar_Dau_Counter, &b_GenPar_Dau_Counter);
   fChain->SetBranchAddress("GenPar_Idx", &GenPar_Idx, &b_GenPar_Idx);
   fChain->SetBranchAddress("GenPar_Mom1_Idx", &GenPar_Mom1_Idx, &b_GenPar_Mom1_Idx);
   fChain->SetBranchAddress("GenPar_Mom2_Idx", &GenPar_Mom2_Idx, &b_GenPar_Mom2_Idx);
   fChain->SetBranchAddress("GenPar_Mom_Counter", &GenPar_Mom_Counter, &b_GenPar_Mom_Counter);
   fChain->SetBranchAddress("GenPar_Status", &GenPar_Status, &b_GenPar_Status);
   fChain->SetBranchAddress("GenPar_pdgId", &GenPar_pdgId, &b_GenPar_pdgId);
   fChain->SetBranchAddress("Gen_EventWeight", &Gen_EventWeight, &b_Gen_EventWeight);
   fChain->SetBranchAddress("GenTop", &GenTop, &b_GenTop);
   fChain->SetBranchAddress("GenAnTop", &GenAnTop, &b_GenAnTop);
   fChain->SetBranchAddress("GenBJet", &GenBJet, &b_GenBJet);
   fChain->SetBranchAddress("GenBJet_Count", &GenBJet_Count, &b_GenBJet_Count);
   fChain->SetBranchAddress("GenBHad", &GenBHad, &b_GenBHad);
   fChain->SetBranchAddress("GenBHad_Count", &GenBHad_Count, &b_GenBHad_Count);
   fChain->SetBranchAddress("GenBHad_pdgId", &GenBHad_pdgId, &b_GenBHad_pdgId);
   fChain->SetBranchAddress("GenBHad_Flavour", &GenBHad_Flavour, &b_GenBHad_Flavour);
   fChain->SetBranchAddress("GenBHad_FromTopWeakDecay", &GenBHad_FromTopWeakDecay, &b_GenBHad_FromTopWeakDecay);
   fChain->SetBranchAddress("Jet", &Jet, &b_Jet);
   fChain->SetBranchAddress("RawJet", &RawJet, &b_RawJet);
   fChain->SetBranchAddress("Jet_Charge", &Jet_Charge, &b_Jet_Charge);
   fChain->SetBranchAddress("Jet_Count", &Jet_Count, &b_Jet_Count);
   fChain->SetBranchAddress("Jet_EnShiftedDown", &Jet_EnShiftedDown, &b_Jet_EnShiftedDown);
   fChain->SetBranchAddress("Jet_EnShiftedUp", &Jet_EnShiftedUp, &b_Jet_EnShiftedUp);
   fChain->SetBranchAddress("Jet_PartonFlavour", &Jet_PartonFlavour, &b_Jet_PartonFlavour);
   fChain->SetBranchAddress("Jet_HadronFlavour", &Jet_HadronFlavour, &b_Jet_HadronFlavour);
   fChain->SetBranchAddress("Jet_PFId", &Jet_PFId, &b_Jet_PFId);
   fChain->SetBranchAddress("Jet_PFIdVeto", &Jet_PFIdVeto, &b_Jet_PFIdVeto);
   fChain->SetBranchAddress("Jet_PhiResolution_MC", &Jet_PhiResolution_MC, &b_Jet_PhiResolution_MC);
   fChain->SetBranchAddress("Jet_PhiResolution_DATA", &Jet_PhiResolution_DATA, &b_Jet_PhiResolution_DATA);
   fChain->SetBranchAddress("Jet_EnergyResolution_MC", &Jet_EnergyResolution_MC, &b_Jet_EnergyResolution_MC);
   fChain->SetBranchAddress("Jet_EnergyResolution_DATA", &Jet_EnergyResolution_DATA, &b_Jet_EnergyResolution_DATA);
   fChain->SetBranchAddress("Jet_EnergyResolution_SF", &Jet_EnergyResolution_SF, &b_Jet_EnergyResolution_SF);
   fChain->SetBranchAddress("Jet_EnergyResolution_SFDown", &Jet_EnergyResolution_SFDown, &b_Jet_EnergyResolution_SFDown);
   fChain->SetBranchAddress("Jet_EnergyResolution_SFUp", &Jet_EnergyResolution_SFUp, &b_Jet_EnergyResolution_SFUp);
   fChain->SetBranchAddress("Jet_PileUpId", &Jet_PileUpId, &b_Jet_PileUpId);
   fChain->SetBranchAddress("Jet_PileUpMVA", &Jet_PileUpMVA, &b_Jet_PileUpMVA);
   fChain->SetBranchAddress("Jet_bDisc", &Jet_bDisc, &b_Jet_bDisc);
   fChain->SetBranchAddress("Jet_bDisc_Name", &Jet_bDisc_Name, &b_Jet_bDisc_Name);
   fChain->SetBranchAddress("Jet_bDisc_Value", &Jet_bDisc_Value, &b_Jet_bDisc_Value);
   fChain->SetBranchAddress("Jet_isJet", &Jet_isJet, &b_Jet_isJet);
   fChain->SetBranchAddress("LHE_Central", &LHE_Central, &b_LHE_Central);
   fChain->SetBranchAddress("LHE_Id", &LHE_Id, &b_LHE_Id);
   fChain->SetBranchAddress("LHE_Weight", &LHE_Weight, &b_LHE_Weight);
   fChain->SetBranchAddress("Frag_Cen_Weight", &Frag_Cen_Weight, &b_Frag_Cen_Weight);
   fChain->SetBranchAddress("Frag_Up_Weight", &Frag_Up_Weight, &b_Frag_Up_Weight);
   fChain->SetBranchAddress("Frag_Down_Weight", &Frag_Down_Weight, &b_Frag_Down_Weight);
   fChain->SetBranchAddress("Frag_Peterson_Weight", &Frag_Peterson_Weight, &b_Frag_Peterson_Weight);
   fChain->SetBranchAddress("Semilep_BrUp_Weight", &Semilep_BrUp_Weight, &b_Semilep_BrUp_Weight);
   fChain->SetBranchAddress("Semilep_BrDown_Weight", &Semilep_BrDown_Weight, &b_Semilep_BrDown_Weight);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("MET_Significance", &MET_Significance, &b_MET_Significance);
   fChain->SetBranchAddress("Muon", &Muon, &b_Muon);
   fChain->SetBranchAddress("GenMuon", &GenMuon, &b_GenMuon);
   fChain->SetBranchAddress("Muon_Charge", &Muon_Charge, &b_Muon_Charge);
   fChain->SetBranchAddress("Muon_Count", &Muon_Count, &b_Muon_Count);
   fChain->SetBranchAddress("Muon_PFIsodBeta03", &Muon_PFIsodBeta03, &b_Muon_PFIsodBeta03);
   fChain->SetBranchAddress("Muon_PFIsodBeta04", &Muon_PFIsodBeta04, &b_Muon_PFIsodBeta04);
   fChain->SetBranchAddress("Muon_isHighPt", &Muon_isHighPt, &b_Muon_isHighPt);
   fChain->SetBranchAddress("Muon_isHighTrkPt", &Muon_isHighTrkPt, &b_Muon_isHighTrkPt);
   fChain->SetBranchAddress("Muon_isLoose", &Muon_isLoose, &b_Muon_isLoose);
   fChain->SetBranchAddress("Muon_isMedium", &Muon_isMedium, &b_Muon_isMedium);
   fChain->SetBranchAddress("Muon_isMediumPrompt", &Muon_isMediumPrompt, &b_Muon_isMediumPrompt);
   fChain->SetBranchAddress("Muon_isSoft", &Muon_isSoft, &b_Muon_isSoft);
   fChain->SetBranchAddress("Muon_isTight", &Muon_isTight, &b_Muon_isTight);
   fChain->SetBranchAddress("Muon_isPFIsoVeryLoose", &Muon_isPFIsoVeryLoose, &b_Muon_isPFIsoVeryLoose);
   fChain->SetBranchAddress("Muon_isPFIsoLoose", &Muon_isPFIsoLoose, &b_Muon_isPFIsoLoose);
   fChain->SetBranchAddress("Muon_isPFIsoMedium", &Muon_isPFIsoMedium, &b_Muon_isPFIsoMedium);
   fChain->SetBranchAddress("Muon_isPFIsoTight", &Muon_isPFIsoTight, &b_Muon_isPFIsoTight);
   fChain->SetBranchAddress("Muon_isPFIsoVeryTight", &Muon_isPFIsoVeryTight, &b_Muon_isPFIsoVeryTight);
   fChain->SetBranchAddress("Muon_isPFIsoVeryVeryTight", &Muon_isPFIsoVeryVeryTight, &b_Muon_isPFIsoVeryVeryTight);
   fChain->SetBranchAddress("Muon_pdgId", &Muon_pdgId, &b_Muon_pdgId);
   fChain->SetBranchAddress("Muon_rand1", &Muon_rand1, &b_Muon_rand1);
   fChain->SetBranchAddress("Muon_rand2", &Muon_rand2, &b_Muon_rand2);
   fChain->SetBranchAddress("Muon_relIso03", &Muon_relIso03, &b_Muon_relIso03);
   fChain->SetBranchAddress("Muon_relIso04", &Muon_relIso04, &b_Muon_relIso04);
   fChain->SetBranchAddress("Muon_trackerLayers", &Muon_trackerLayers, &b_Muon_trackerLayers);
   fChain->SetBranchAddress("Muon_tuneP_Pt", &Muon_tuneP_Pt, &b_Muon_tuneP_Pt);
   fChain->SetBranchAddress("Muon_tuneP_Eta", &Muon_tuneP_Eta, &b_Muon_tuneP_Eta);
   fChain->SetBranchAddress("Muon_tuneP_Phi", &Muon_tuneP_Phi, &b_Muon_tuneP_Phi);
   fChain->SetBranchAddress("Muon_tuneP_Charge", &Muon_tuneP_Charge, &b_Muon_tuneP_Charge);
   fChain->SetBranchAddress("PDFWeight_BjorkenX1", &PDFWeight_BjorkenX1, &b_PDFWeight_BjorkenX1);
   fChain->SetBranchAddress("PDFWeight_BjorkenX2", &PDFWeight_BjorkenX2, &b_PDFWeight_BjorkenX2);
   fChain->SetBranchAddress("PDFWeight_Cent", &PDFWeight_Cent, &b_PDFWeight_Cent);
   fChain->SetBranchAddress("PDFWeight_Cent_Down", &PDFWeight_Cent_Down, &b_PDFWeight_Cent_Down);
   fChain->SetBranchAddress("PDFWeight_Cent_Up", &PDFWeight_Cent_Up, &b_PDFWeight_Cent_Up);
   fChain->SetBranchAddress("PDFWeight_Id1", &PDFWeight_Id1, &b_PDFWeight_Id1);
   fChain->SetBranchAddress("PDFWeight_Id2", &PDFWeight_Id2, &b_PDFWeight_Id2);
   fChain->SetBranchAddress("PDFWeight_PDF1", &PDFWeight_PDF1, &b_PDFWeight_PDF1);
   fChain->SetBranchAddress("PDFWeight_PDF2", &PDFWeight_PDF2, &b_PDFWeight_PDF2);
   fChain->SetBranchAddress("PDFWeight_Q", &PDFWeight_Q, &b_PDFWeight_Q);
   fChain->SetBranchAddress("PDFWeight_Var1_Down", &PDFWeight_Var1_Down, &b_PDFWeight_Var1_Down);
   fChain->SetBranchAddress("PDFWeight_Var1_Up", &PDFWeight_Var1_Up, &b_PDFWeight_Var1_Up);
   fChain->SetBranchAddress("PDFWeight_Var2_Down", &PDFWeight_Var2_Down, &b_PDFWeight_Var2_Down);
   fChain->SetBranchAddress("PDFWeight_Var2_Up", &PDFWeight_Var2_Up, &b_PDFWeight_Var2_Up);
   fChain->SetBranchAddress("Photon", &Photon, &b_Photon);
   fChain->SetBranchAddress("Photon_Count", &Photon_Count, &b_Photon_Count);
   fChain->SetBranchAddress("Photon_SCB_Loose", &Photon_SCB_Loose, &b_Photon_SCB_Loose);
   fChain->SetBranchAddress("Photon_SCB_Medium", &Photon_SCB_Medium, &b_Photon_SCB_Medium);
   fChain->SetBranchAddress("Photon_SCB_Tight", &Photon_SCB_Tight, &b_Photon_SCB_Tight);
   fChain->SetBranchAddress("Photon_ChaHadIso", &Photon_ChaHadIso, &b_Photon_ChaHadIso);
   fChain->SetBranchAddress("Photon_NeuHadIso", &Photon_NeuHadIso, &b_Photon_NeuHadIso);
   fChain->SetBranchAddress("Photon_PhoIso", &Photon_PhoIso, &b_Photon_PhoIso);
   fChain->SetBranchAddress("Photon_WorstChargedIso", &Photon_WorstChargedIso, &b_Photon_WorstChargedIso);
   fChain->SetBranchAddress("Photon_ChaHadEffArea", &Photon_ChaHadEffArea, &b_Photon_ChaHadEffArea);
   fChain->SetBranchAddress("Photon_NeuHadEffArea", &Photon_NeuHadEffArea, &b_Photon_NeuHadEffArea);
   fChain->SetBranchAddress("Photon_PhoHadEffArea", &Photon_PhoHadEffArea, &b_Photon_PhoHadEffArea);
   fChain->SetBranchAddress("Photon_R9", &Photon_R9, &b_Photon_R9);
   fChain->SetBranchAddress("Photon_HoverE", &Photon_HoverE, &b_Photon_HoverE);
   fChain->SetBranchAddress("Photon_SuperCluster_Eta", &Photon_SuperCluster_Eta, &b_Photon_SuperCluster_Eta);
   fChain->SetBranchAddress("Photon_SuperCluster_Phi", &Photon_SuperCluster_Phi, &b_Photon_SuperCluster_Phi);
   fChain->SetBranchAddress("Photon_SuperCluster_EtaWidth", &Photon_SuperCluster_EtaWidth, &b_Photon_SuperCluster_EtaWidth);
   fChain->SetBranchAddress("Photon_SuperCluster_PhiWidth", &Photon_SuperCluster_PhiWidth, &b_Photon_SuperCluster_PhiWidth);
   fChain->SetBranchAddress("Photon_SuperCluster_Energy", &Photon_SuperCluster_Energy, &b_Photon_SuperCluster_Energy);
   fChain->SetBranchAddress("Photon_SuperCluster_RawEnergy", &Photon_SuperCluster_RawEnergy, &b_Photon_SuperCluster_RawEnergy);
   fChain->SetBranchAddress("Photon_ElectronVeto", &Photon_ElectronVeto, &b_Photon_ElectronVeto);
   fChain->SetBranchAddress("Photon_Full5x5_SigmaIetaIeta", &Photon_Full5x5_SigmaIetaIeta, &b_Photon_Full5x5_SigmaIetaIeta);
   fChain->SetBranchAddress("Photon_MVANonTrig_Tight", &Photon_MVANonTrig_Tight, &b_Photon_MVANonTrig_Tight);
   fChain->SetBranchAddress("PV_Count", &PV_Count, &b_PV_Count);
   fChain->SetBranchAddress("PileUp_Count_Interaction", &PileUp_Count_Interaction, &b_PileUp_Count_Interaction);
   fChain->SetBranchAddress("PileUp_Count_Intime", &PileUp_Count_Intime, &b_PileUp_Count_Intime);
   fChain->SetBranchAddress("Rho", &Rho, &b_Rho);
   fChain->SetBranchAddress("Trigger_Name", &Trigger_Name, &b_Trigger_Name);
   fChain->SetBranchAddress("Trigger_PreScale", &Trigger_PreScale, &b_Trigger_PreScale);
   fChain->SetBranchAddress("Trigger_isError", &Trigger_isError, &b_Trigger_isError);
   fChain->SetBranchAddress("Trigger_isPass", &Trigger_isPass, &b_Trigger_isPass);
   fChain->SetBranchAddress("Trigger_isRun", &Trigger_isRun, &b_Trigger_isRun);
   fChain->SetBranchAddress("Vertex_SumPtSquare", &Vertex_SumPtSquare, &b_Vertex_SumPtSquare);
   fChain->SetBranchAddress("Vertex_X", &Vertex_X, &b_Vertex_X);
   fChain->SetBranchAddress("Vertex_X_Error", &Vertex_X_Error, &b_Vertex_X_Error);
   fChain->SetBranchAddress("Vertex_Y", &Vertex_Y, &b_Vertex_Y);
   fChain->SetBranchAddress("Vertex_Y_Error", &Vertex_Y_Error, &b_Vertex_Y_Error);
   fChain->SetBranchAddress("Vertex_Z", &Vertex_Z, &b_Vertex_Z);
   fChain->SetBranchAddress("Vertex_Z_Error", &Vertex_Z_Error, &b_Vertex_Z_Error);
   Notify();
}

Bool_t SSBTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SSBTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SSBTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SSBTree_cxx
