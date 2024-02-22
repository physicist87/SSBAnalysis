void Draw(){
    gStyle->SetPaintTextFormat("4.3f");
    for (int i = 0; i < 3; ++i)
    {
       TString Channel = "";
       TString channel = "";
       if (i == 0){
          Channel = "MuMu";
          channel = "mumu";
       }
       if (i == 1){
          Channel = "ElEl";
          channel = "ee";
       }
       if (i == 2){
          Channel = "MuEl";
          channel = "emu";
       }
       TString fileName = ;
       //TFile* f1 = TFile::Open("./triggerSummary_mumu.root");
       TFile* f1 = TFile::Open(Form("./triggerSummary_%s.root",channel.Data()));
       TH2F* tmp = (TH2F*)f1->Get("scalefactor_eta2d_with_syst");
       tmp->SetTitle("");
       TCanvas* c1 = new TCanvas("c1");
       c1->cd();
       tmp->Draw("colzTexte");
       c1->SaveAs(Form("TriggerSF_%s.pdf",Channel.Data()));
    }
}
