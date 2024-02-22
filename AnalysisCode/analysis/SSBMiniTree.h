//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Sep 16 20:12:13 2017 by ROOT version 5.34/32
// from TTree SSBMiniTree/CMSAnalysis
// found on file: dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/Run2016/MiniTree/MC/Version2/TTJets_Signal/TTJets_Signal_1.root
//////////////////////////////////////////////////////////

#ifndef SSBMiniTree_h
#define SSBMiniTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TClonesArray.h>
#include <vector>
#include <vector>
#include <vector>
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.
using namespace std;
class SSBMiniTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Info_EventNumber;
   Int_t           Info_RunNumber;
   Int_t           Info_Luminosity;
   Int_t           Channel_Idx;
   Double_t        Gen_EventWeight;
   TClonesArray    *GenTop;
   TClonesArray    *GenAnTop;
   vector<string>  *Trigger_Name;
   vector<bool>    *Trigger_isError;
   vector<bool>    *Trigger_isPass;
   vector<bool>    *Trigger_isRun;
   vector<bool>    *Filter_PV;
   vector<double>  *Vertex_SumPtSquare;
   vector<double>  *Vertex_X;
   vector<double>  *Vertex_X_Error;
   vector<double>  *Vertex_Y;
   vector<double>  *Vertex_Y_Error;
   vector<double>  *Vertex_Z;
   vector<double>  *Vertex_Z_Error;
   vector<string>  *METFilter_Name;
   vector<bool>    *METFilter_isError;
   vector<bool>    *METFilter_isPass;
   vector<bool>    *METFilter_isRun;
   vector<string>  *METFilterAdd_Name;
   vector<bool>    *METFilterAdd_isPass;
   TClonesArray    *Jet;
   Int_t           Jet_Count;
   vector<double>  *Jet_EnShiftedDown;
   vector<double>  *Jet_EnShiftedUp;
   vector<int>     *Jet_PartonFlavour;
   vector<int>     *Jet_HadronFlavour;
   vector<int>     *Jet_PFId;
   vector<float>   *Jet_bDisc;
   vector<double>  *Jet_EnergyResolution_SF;
   vector<double>  *Jet_EnergyResolution_SFDown;
   vector<double>  *Jet_EnergyResolution_SFUp;
   TClonesArray    *Elec;
   Int_t           Elec_Count;
   vector<double>  *Elec_PFIsoRho03;
   vector<double>  *Elec_PFIsoRho04;
   vector<bool>    *Elec_SCB_Loose;
   vector<bool>    *Elec_SCB_Medium;
   vector<bool>    *Elec_SCB_Tight;
   vector<bool>    *Elec_SCB_Veto;
   vector<double>  *Elec_Supercluster_Eta;
   vector<double>  *Elec_Track_GsfdXY;
   vector<double>  *Elec_Track_GsfdZ;
   vector<int>     *Elec_Charge;
   vector<bool>    *Elec_Conversion;
   vector<bool>    *Elec_ChargeId_GsfPx;
   vector<bool>    *Elec_ChargeId_GsfCtfPx;
   TClonesArray    *Muon;
   vector<int>     *Muon_Charge;
   Int_t           Muon_Count;
   vector<double>  *Muon_PFIsodBeta03;
   vector<double>  *Muon_PFIsodBeta04;
   vector<bool>    *Muon_isLoose;
   vector<bool>    *Muon_isTight;
   TClonesArray    *MET;
   vector<double>  *MET_JetEnShiftedUp_PT;
   vector<double>  *MET_JetEnShiftedUp_Phi;
   vector<double>  *MET_JetEnShiftedDown_PT;
   vector<double>  *MET_JetEnShiftedDown_Phi;
   TClonesArray    *METMUEGCleanCor;
   vector<double>  *METMUEGCleanCor_JetEnShiftedUp_PT;
   vector<double>  *METMUEGCleanCor_JetEnShiftedUp_Phi;
   vector<double>  *METMUEGCleanCor_JetEnShiftedDown_PT;
   vector<double>  *METMUEGCleanCor_JetEnShiftedDown_Phi;
   TClonesArray    *METMUCleanCor;
   vector<double>  *METMUCleanCor_JetEnShiftedUp_PT;
   vector<double>  *METMUCleanCor_JetEnShiftedUp_Phi;
   vector<double>  *METMUCleanCor_JetEnShiftedDown_PT;
   vector<double>  *METMUCleanCor_JetEnShiftedDown_Phi;
   Float_t         PileUp_Count_Intime;

   // List of branches
   TBranch        *b_Info_EventNumber;   //!
   TBranch        *b_Info_RunNumber;   //!
   TBranch        *b_Info_Luminosity;   //!
   TBranch        *b_Channel_Idx;   //!
   TBranch        *b_Gen_EventWeight;   //!
   TBranch        *b_GenTop;   //!
   TBranch        *b_GenAnTop;   //!
   TBranch        *b_Trigger_Name;   //!
   TBranch        *b_Trigger_isError;   //!
   TBranch        *b_Trigger_isPass;   //!
   TBranch        *b_Trigger_isRun;   //!
   TBranch        *b_Filter_PV;   //!
   TBranch        *b_Vertex_SumPtSquare;   //!
   TBranch        *b_Vertex_X;   //!
   TBranch        *b_Vertex_X_Error;   //!
   TBranch        *b_Vertex_Y;   //!
   TBranch        *b_Vertex_Y_Error;   //!
   TBranch        *b_Vertex_Z;   //!
   TBranch        *b_Vertex_Z_Error;   //!
   TBranch        *b_METFilter_Name;   //!
   TBranch        *b_METFilter_isError;   //!
   TBranch        *b_METFilter_isPass;   //!
   TBranch        *b_METFilter_isRun;   //!
   TBranch        *b_METFilterAdd_Name;   //!
   TBranch        *b_METFilterAdd_isPass;   //!
   TBranch        *b_Jet;   //!
   TBranch        *b_Jet_Count;   //!
   TBranch        *b_Jet_EnShiftedDown;   //!
   TBranch        *b_Jet_EnShiftedUp;   //!
   TBranch        *b_Jet_PartonFlavour;   //!
   TBranch        *b_Jet_HadronFlavour;   //!
   TBranch        *b_Jet_PFId;   //!
   TBranch        *b_Jet_bDisc;   //!
   TBranch        *b_Jet_EnergyResolution_SF;   //!
   TBranch        *b_Jet_EnergyResolution_SFDown;   //!
   TBranch        *b_Jet_EnergyResolution_SFUp;   //!
   TBranch        *b_Elec;   //!
   TBranch        *b_Elec_Count;   //!
   TBranch        *b_Elec_PFIsoRho03;   //!
   TBranch        *b_Elec_PFIsoRho04;   //!
   TBranch        *b_Elec_SCB_Loose;   //!
   TBranch        *b_Elec_SCB_Medium;   //!
   TBranch        *b_Elec_SCB_Tight;   //!
   TBranch        *b_Elec_SCB_Veto;   //!
   TBranch        *b_Elec_Supercluster_Eta;   //!
   TBranch        *b_Elec_Track_GsfdXY;   //!
   TBranch        *b_Elec_Track_GsfdZ;   //!
   TBranch        *b_Elec_Charge;   //!
   TBranch        *b_Elec_Conversion;   //!
   TBranch        *b_Elec_ChargeId_GsfPx;   //!
   TBranch        *b_Elec_ChargeId_GsfCtfPx;   //!
   TBranch        *b_Muon;   //!
   TBranch        *b_Muon_Charge;   //!
   TBranch        *b_Muon_Count;   //!
   TBranch        *b_Muon_PFIsodBeta03;   //!
   TBranch        *b_Muon_PFIsodBeta04;   //!
   TBranch        *b_Muon_isLoose;   //!
   TBranch        *b_Muon_isTight;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_MET_JetEnShiftedUp_PT;   //!
   TBranch        *b_MET_JetEnShiftedUp_Phi;   //!
   TBranch        *b_MET_JetEnShiftedDown_PT;   //!
   TBranch        *b_MET_JetEnShiftedDown_Phi;   //!
   TBranch        *b_METMUEGCleanCor;   //!
   TBranch        *b_METMUEGCleanCor_JetEnShiftedUp_PT;   //!
   TBranch        *b_METMUEGCleanCor_JetEnShiftedUp_Phi;   //!
   TBranch        *b_METMUEGCleanCor_JetEnShiftedDown_PT;   //!
   TBranch        *b_METMUEGCleanCor_JetEnShiftedDown_Phi;   //!
   TBranch        *b_METMUCleanCor;   //!
   TBranch        *b_METMUCleanCor_JetEnShiftedUp_PT;   //!
   TBranch        *b_METMUCleanCor_JetEnShiftedUp_Phi;   //!
   TBranch        *b_METMUCleanCor_JetEnShiftedDown_PT;   //!
   TBranch        *b_METMUCleanCor_JetEnShiftedDown_Phi;   //!
   TBranch        *b_PileUp_Count_Intime;   //!

   SSBMiniTree(TTree *tree=0);
   virtual ~SSBMiniTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef SSBMiniTree_cxx
SSBMiniTree::SSBMiniTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
/*      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/Run2016/MiniTree/MC/Version2/TTJets_Signal/TTJets_Signal_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/Run2016/MiniTree/MC/Version2/TTJets_Signal/TTJets_Signal_1.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/Run2016/MiniTree/MC/Version2/TTJets_Signal/TTJets_Signal_1.root:/ssbanalyzer");
      dir->GetObject("SSBMiniTree",tree);
*/
   }
   Init(tree);
}

SSBMiniTree::~SSBMiniTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SSBMiniTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SSBMiniTree::LoadTree(Long64_t entry)
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

void SSBMiniTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   GenTop = 0;
   GenAnTop = 0;
   Trigger_Name = 0;
   Trigger_isError = 0;
   Trigger_isPass = 0;
   Trigger_isRun = 0;
   Filter_PV = 0;
   Vertex_SumPtSquare = 0;
   Vertex_X = 0;
   Vertex_X_Error = 0;
   Vertex_Y = 0;
   Vertex_Y_Error = 0;
   Vertex_Z = 0;
   Vertex_Z_Error = 0;
   METFilter_Name = 0;
   METFilter_isError = 0;
   METFilter_isPass = 0;
   METFilter_isRun = 0;
   METFilterAdd_Name = 0;
   METFilterAdd_isPass = 0;
   Jet = 0;
   Jet_EnShiftedDown = 0;
   Jet_EnShiftedUp = 0;
   Jet_PartonFlavour = 0;
   Jet_HadronFlavour = 0;
   Jet_PFId = 0;
   Jet_bDisc = 0;
   Jet_EnergyResolution_SF = 0;
   Jet_EnergyResolution_SFDown = 0;
   Jet_EnergyResolution_SFUp = 0;
   Elec = 0;
   Elec_PFIsoRho03 = 0;
   Elec_PFIsoRho04 = 0;
   Elec_SCB_Loose = 0;
   Elec_SCB_Medium = 0;
   Elec_SCB_Tight = 0;
   Elec_SCB_Veto = 0;
   Elec_Supercluster_Eta = 0;
   Elec_Track_GsfdXY = 0;
   Elec_Track_GsfdZ = 0;
   Elec_Charge = 0;
   Elec_Conversion = 0;
   Elec_ChargeId_GsfPx = 0;
   Elec_ChargeId_GsfCtfPx = 0;
   Muon = 0;
   Muon_Charge = 0;
   Muon_PFIsodBeta03 = 0;
   Muon_PFIsodBeta04 = 0;
   Muon_isLoose = 0;
   Muon_isTight = 0;
   MET = 0;
   MET_JetEnShiftedUp_PT = 0;
   MET_JetEnShiftedUp_Phi = 0;
   MET_JetEnShiftedDown_PT = 0;
   MET_JetEnShiftedDown_Phi = 0;
   METMUEGCleanCor = 0;
   METMUEGCleanCor_JetEnShiftedUp_PT = 0;
   METMUEGCleanCor_JetEnShiftedUp_Phi = 0;
   METMUEGCleanCor_JetEnShiftedDown_PT = 0;
   METMUEGCleanCor_JetEnShiftedDown_Phi = 0;
   METMUCleanCor = 0;
   METMUCleanCor_JetEnShiftedUp_PT = 0;
   METMUCleanCor_JetEnShiftedUp_Phi = 0;
   METMUCleanCor_JetEnShiftedDown_PT = 0;
   METMUCleanCor_JetEnShiftedDown_Phi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Info_EventNumber", &Info_EventNumber, &b_Info_EventNumber);
   fChain->SetBranchAddress("Info_RunNumber", &Info_RunNumber, &b_Info_RunNumber);
   fChain->SetBranchAddress("Info_Luminosity", &Info_Luminosity, &b_Info_Luminosity);
   fChain->SetBranchAddress("Channel_Idx", &Channel_Idx, &b_Channel_Idx);
   fChain->SetBranchAddress("Gen_EventWeight", &Gen_EventWeight, &b_Gen_EventWeight);
   fChain->SetBranchAddress("GenTop", &GenTop, &b_GenTop);
   fChain->SetBranchAddress("GenAnTop", &GenAnTop, &b_GenAnTop);
   fChain->SetBranchAddress("Trigger_Name", &Trigger_Name, &b_Trigger_Name);
   fChain->SetBranchAddress("Trigger_isError", &Trigger_isError, &b_Trigger_isError);
   fChain->SetBranchAddress("Trigger_isPass", &Trigger_isPass, &b_Trigger_isPass);
   fChain->SetBranchAddress("Trigger_isRun", &Trigger_isRun, &b_Trigger_isRun);
   fChain->SetBranchAddress("Filter_PV", &Filter_PV, &b_Filter_PV);
   fChain->SetBranchAddress("Vertex_SumPtSquare", &Vertex_SumPtSquare, &b_Vertex_SumPtSquare);
   fChain->SetBranchAddress("Vertex_X", &Vertex_X, &b_Vertex_X);
   fChain->SetBranchAddress("Vertex_X_Error", &Vertex_X_Error, &b_Vertex_X_Error);
   fChain->SetBranchAddress("Vertex_Y", &Vertex_Y, &b_Vertex_Y);
   fChain->SetBranchAddress("Vertex_Y_Error", &Vertex_Y_Error, &b_Vertex_Y_Error);
   fChain->SetBranchAddress("Vertex_Z", &Vertex_Z, &b_Vertex_Z);
   fChain->SetBranchAddress("Vertex_Z_Error", &Vertex_Z_Error, &b_Vertex_Z_Error);
   fChain->SetBranchAddress("METFilter_Name", &METFilter_Name, &b_METFilter_Name);
   fChain->SetBranchAddress("METFilter_isError", &METFilter_isError, &b_METFilter_isError);
   fChain->SetBranchAddress("METFilter_isPass", &METFilter_isPass, &b_METFilter_isPass);
   fChain->SetBranchAddress("METFilter_isRun", &METFilter_isRun, &b_METFilter_isRun);
   fChain->SetBranchAddress("METFilterAdd_Name", &METFilterAdd_Name, &b_METFilterAdd_Name);
   fChain->SetBranchAddress("METFilterAdd_isPass", &METFilterAdd_isPass, &b_METFilterAdd_isPass);
   fChain->SetBranchAddress("Jet", &Jet, &b_Jet);
   fChain->SetBranchAddress("Jet_Count", &Jet_Count, &b_Jet_Count);
   fChain->SetBranchAddress("Jet_EnShiftedDown", &Jet_EnShiftedDown, &b_Jet_EnShiftedDown);
   fChain->SetBranchAddress("Jet_EnShiftedUp", &Jet_EnShiftedUp, &b_Jet_EnShiftedUp);
   fChain->SetBranchAddress("Jet_PartonFlavour", &Jet_PartonFlavour, &b_Jet_PartonFlavour);
   fChain->SetBranchAddress("Jet_HadronFlavour", &Jet_HadronFlavour, &b_Jet_HadronFlavour);
   fChain->SetBranchAddress("Jet_PFId", &Jet_PFId, &b_Jet_PFId);
   fChain->SetBranchAddress("Jet_bDisc", &Jet_bDisc, &b_Jet_bDisc);
   fChain->SetBranchAddress("Jet_EnergyResolution_SF", &Jet_EnergyResolution_SF, &b_Jet_EnergyResolution_SF);
   fChain->SetBranchAddress("Jet_EnergyResolution_SFDown", &Jet_EnergyResolution_SFDown, &b_Jet_EnergyResolution_SFDown);
   fChain->SetBranchAddress("Jet_EnergyResolution_SFUp", &Jet_EnergyResolution_SFUp, &b_Jet_EnergyResolution_SFUp);
   fChain->SetBranchAddress("Elec", &Elec, &b_Elec);
   fChain->SetBranchAddress("Elec_Count", &Elec_Count, &b_Elec_Count);
   fChain->SetBranchAddress("Elec_PFIsoRho03", &Elec_PFIsoRho03, &b_Elec_PFIsoRho03);
   fChain->SetBranchAddress("Elec_PFIsoRho04", &Elec_PFIsoRho04, &b_Elec_PFIsoRho04);
   fChain->SetBranchAddress("Elec_SCB_Loose", &Elec_SCB_Loose, &b_Elec_SCB_Loose);
   fChain->SetBranchAddress("Elec_SCB_Medium", &Elec_SCB_Medium, &b_Elec_SCB_Medium);
   fChain->SetBranchAddress("Elec_SCB_Tight", &Elec_SCB_Tight, &b_Elec_SCB_Tight);
   fChain->SetBranchAddress("Elec_SCB_Veto", &Elec_SCB_Veto, &b_Elec_SCB_Veto);
   fChain->SetBranchAddress("Elec_Supercluster_Eta", &Elec_Supercluster_Eta, &b_Elec_Supercluster_Eta);
   fChain->SetBranchAddress("Elec_Track_GsfdXY", &Elec_Track_GsfdXY, &b_Elec_Track_GsfdXY);
   fChain->SetBranchAddress("Elec_Track_GsfdZ", &Elec_Track_GsfdZ, &b_Elec_Track_GsfdZ);
   fChain->SetBranchAddress("Elec_Charge", &Elec_Charge, &b_Elec_Charge);
   fChain->SetBranchAddress("Elec_Conversion", &Elec_Conversion, &b_Elec_Conversion);
   fChain->SetBranchAddress("Elec_ChargeId_GsfPx", &Elec_ChargeId_GsfPx, &b_Elec_ChargeId_GsfPx);
   fChain->SetBranchAddress("Elec_ChargeId_GsfCtfPx", &Elec_ChargeId_GsfCtfPx, &b_Elec_ChargeId_GsfCtfPx);
   fChain->SetBranchAddress("Muon", &Muon, &b_Muon);
   fChain->SetBranchAddress("Muon_Charge", &Muon_Charge, &b_Muon_Charge);
   fChain->SetBranchAddress("Muon_Count", &Muon_Count, &b_Muon_Count);
   fChain->SetBranchAddress("Muon_PFIsodBeta03", &Muon_PFIsodBeta03, &b_Muon_PFIsodBeta03);
   fChain->SetBranchAddress("Muon_PFIsodBeta04", &Muon_PFIsodBeta04, &b_Muon_PFIsodBeta04);
   fChain->SetBranchAddress("Muon_isLoose", &Muon_isLoose, &b_Muon_isLoose);
   fChain->SetBranchAddress("Muon_isTight", &Muon_isTight, &b_Muon_isTight);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("MET_JetEnShiftedUp_PT", &MET_JetEnShiftedUp_PT, &b_MET_JetEnShiftedUp_PT);
   fChain->SetBranchAddress("MET_JetEnShiftedUp_Phi", &MET_JetEnShiftedUp_Phi, &b_MET_JetEnShiftedUp_Phi);
   fChain->SetBranchAddress("MET_JetEnShiftedDown_PT", &MET_JetEnShiftedDown_PT, &b_MET_JetEnShiftedDown_PT);
   fChain->SetBranchAddress("MET_JetEnShiftedDown_Phi", &MET_JetEnShiftedDown_Phi, &b_MET_JetEnShiftedDown_Phi);
   fChain->SetBranchAddress("METMUEGCleanCor", &METMUEGCleanCor, &b_METMUEGCleanCor);
   fChain->SetBranchAddress("METMUEGCleanCor_JetEnShiftedUp_PT", &METMUEGCleanCor_JetEnShiftedUp_PT, &b_METMUEGCleanCor_JetEnShiftedUp_PT);
   fChain->SetBranchAddress("METMUEGCleanCor_JetEnShiftedUp_Phi", &METMUEGCleanCor_JetEnShiftedUp_Phi, &b_METMUEGCleanCor_JetEnShiftedUp_Phi);
   fChain->SetBranchAddress("METMUEGCleanCor_JetEnShiftedDown_PT", &METMUEGCleanCor_JetEnShiftedDown_PT, &b_METMUEGCleanCor_JetEnShiftedDown_PT);
   fChain->SetBranchAddress("METMUEGCleanCor_JetEnShiftedDown_Phi", &METMUEGCleanCor_JetEnShiftedDown_Phi, &b_METMUEGCleanCor_JetEnShiftedDown_Phi);
   fChain->SetBranchAddress("METMUCleanCor", &METMUCleanCor, &b_METMUCleanCor);
   fChain->SetBranchAddress("METMUCleanCor_JetEnShiftedUp_PT", &METMUCleanCor_JetEnShiftedUp_PT, &b_METMUCleanCor_JetEnShiftedUp_PT);
   fChain->SetBranchAddress("METMUCleanCor_JetEnShiftedUp_Phi", &METMUCleanCor_JetEnShiftedUp_Phi, &b_METMUCleanCor_JetEnShiftedUp_Phi);
   fChain->SetBranchAddress("METMUCleanCor_JetEnShiftedDown_PT", &METMUCleanCor_JetEnShiftedDown_PT, &b_METMUCleanCor_JetEnShiftedDown_PT);
   fChain->SetBranchAddress("METMUCleanCor_JetEnShiftedDown_Phi", &METMUCleanCor_JetEnShiftedDown_Phi, &b_METMUCleanCor_JetEnShiftedDown_Phi);
   fChain->SetBranchAddress("PileUp_Count_Intime", &PileUp_Count_Intime, &b_PileUp_Count_Intime);
   Notify();
}

Bool_t SSBMiniTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SSBMiniTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SSBMiniTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SSBMiniTree_cxx
