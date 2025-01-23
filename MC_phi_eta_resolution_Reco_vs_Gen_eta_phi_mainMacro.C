#include <vector>
#include <sstream>
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TKey.h>
#include <iostream>
#include "JetCorrector.h"

void MC_phi_eta_resoultion(int infilInd=1) {

  using namespace std;
  gROOT->Reset();

  
  //TFile *HiFor_file = new TFile("/eos/cms/store/group/phys_heavyions/uacharya/simu/merged_JRA_all.root");
  TFile *HiFor_file = new TFile(Form("/eos/cms/store/group/phys_heavyions/uacharya/simu/recon_JRA_files/JRA_%d.root",infilInd));

  TDirectory *dir = (TDirectory*)HiFor_file->Get(Form("ak4pf"));
  if (!dir) {
    cout << "Directory ak4PFJetAnalyzer not found!" << endl;
    return;
  }

  TTree *tree;
  dir->GetObject("t", tree);
  if (!tree) {
    cout << "Tree t not found!" << endl;
    return;
  }
 
  //return;
  
  
  // vector<string> JEC_input_file;
  // JEC_input_file.push_back("/afs/cern.ch/user/u/uacharya/eos/cms_analaysis/JEC_analysis/JRA_step3harvest_L2Relative_AK4PF.txt");
  // JetCorrector JECCor(JEC_input_file);


  vector<string> JEC_input_file;
  //JEC_input_file.push_back("/afs/cern.ch/user/u/uacharya/eos/cms_analaysis/JEC_analysis/JRA_step3harvest_L2Relative_AK4PF.txt");
  JEC_input_file.push_back("/afs/cern.ch/user/u/uacharya/MC_Truth_JEC/CMSSW_14_1_0/src/MC-truth-JEC/condor/Files/ParallelMCL1_L2Relative_AK4PF.txt");
  JetCorrector JECCor(JEC_input_file);

  const int mx = 100000;

				
  //Declaration of leaves types
  vector<int>     npus;
  vector<float>   tnpus;
  vector<float>   zpositions;
  vector<int>     bxns;
  vector<float>   sumpt_lowpt;
  vector<float>   sumpt_highpt;
  vector<int>     ntrks_lowpt;
  vector<int>     ntrks_highpt;
  vector<float>   rhos;
  Float_t         rho;
  Float_t         pthat;
  Float_t         beta;
  Float_t         betaStar;
  Float_t         weight;
  Float_t         refpvz;
  Float_t         pudensity;
  Float_t         gpudensity;
  Long64_t        npv;
  Long64_t        run;
  Long64_t        lumi;
  Long64_t        evt;
  UChar_t         nref;
  vector<unsigned char> refrank;
  vector<int>     refpdgid;
  vector<float>   refe;
  vector<float>*   refpt=nullptr;
  vector<float>*   refeta=nullptr;
  vector<float>*   refphi=nullptr;
  vector<float>   refy;
  vector<float>   refdrjt;
  vector<float>   refarea;
  vector<float>   jte;
  vector<float>*  jtpt=nullptr;
  vector<float>*  jteta=nullptr;
  vector<float>*  jtphi=nullptr;
  vector<float>*  jty=nullptr;
  vector<float>   jtjec;
  vector<float>   jtarea;
  vector<float>*   jtchf=nullptr;
  vector<float>*   jtnhf=nullptr;
  vector<float>*   jtnef=nullptr;
  vector<float>*   jtcef=nullptr;
  vector<float>*   jtmuf=nullptr;
  vector<float>   jthfhf;
  vector<float>   jthfef;
  vector<int>     refnMult;
  vector<int>     refchMult;
  vector<int>     jtnMult;
  vector<int>*    jtchMult=nullptr;
  vector<float>   refdzvtx;
  // // Set branch addresses.
  // // tree->SetBranchAddress("npus",npus);
  // // tree->SetBranchAddress("tnpus",tnpus);
  // // tree->SetBranchAddress("zpositions",zpositions);
  // // tree->SetBranchAddress("bxns",bxns);
  // // tree->SetBranchAddress("sumpt_lowpt",sumpt_lowpt);
  // // tree->SetBranchAddress("sumpt_highpt",sumpt_highpt);
  // // tree->SetBranchAddress("ntrks_lowpt",ntrks_lowpt);
  // // tree->SetBranchAddress("ntrks_highpt",ntrks_highpt);
  // // tree->SetBranchAddress("rhos",rhos);
  // // tree->SetBranchAddress("rho",rho);
  // // tree->SetBranchAddress("pthat",pthat);
  // // tree->SetBranchAddress("beta",beta);
  // // tree->SetBranchAddress("betaStar",betaStar);
  // // tree->SetBranchAddress("weight",weight);
  // // tree->SetBranchAddress("refpvz",refpvz);
  // // tree->SetBranchAddress("pudensity",pudensity);
  // // tree->SetBranchAddress("gpudensity",gpudensity);
  // // tree->SetBranchAddress("npv",npv);
  // // tree->SetBranchAddress("run",run);
  // // tree->SetBranchAddress("lumi",lumi);
  // // tree->SetBranchAddress("evt",evt);
  tree->SetBranchAddress("nref",&nref);
  // // tree->SetBranchAddress("refrank",refrank);
  // // tree->SetBranchAddress("refpdgid",refpdgid);
  // // tree->SetBranchAddress("refe",refe);
  tree->SetBranchAddress("refpt",&refpt);
  tree->SetBranchAddress("refeta",&refeta);
  tree->SetBranchAddress("refphi",&refphi);
  // // tree->SetBranchAddress("refy",refy);
  // // tree->SetBranchAddress("refdrjt",refdrjt);
  // // tree->SetBranchAddress("refarea",refarea);
  // // tree->SetBranchAddress("jte",jte);
  tree->SetBranchAddress("jtpt",&jtpt);
  tree->SetBranchAddress("jteta",&jteta);
  tree->SetBranchAddress("jtphi",&jtphi);
  tree->SetBranchAddress("jty",&jty);
  // // tree->SetBranchAddress("jtjec",jtjec);
  // // tree->SetBranchAddress("jtarea",jtarea);
  tree->SetBranchAddress("jtchf",&jtchf);
  tree->SetBranchAddress("jtnhf",&jtnhf);
  tree->SetBranchAddress("jtnef",&jtnef);
  tree->SetBranchAddress("jtcef",&jtcef);
  tree->SetBranchAddress("jtmuf",&jtmuf);
  // // tree->SetBranchAddress("jthfhf",jthfhf);
  // // tree->SetBranchAddress("jthfef",jthfef);
  // // tree->SetBranchAddress("refnMult",refnMult);
  // // tree->SetBranchAddress("refchMult",refchMult);
  // // tree->SetBranchAddress("jtnMult",jtnMult);
  tree->SetBranchAddress("jtchMult",&jtchMult);
  // // tree->SetBranchAddress("refdzvtx",refdzvtx);


  				   
 
  Long64_t nentries = tree->GetEntries();
  const int nEta = 50;
  float etaBin[] = {-5.191,  -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.051, -1.93, -1.791, -1.653, -1.566, -1.479, -1.392, -1.305, -1.174, -1.044, -0.913, -0.783, -0.652, -0.522, -0.391, -0.261, -0.130, 0.0, 0.130, 0.261, 0.391, 0.522, 0.652, 0.783, 0.913, 1.044, 1.174, 1.305, 1.392, 1.479, 1.566, 1.653, 1.791, 1.93, 2.051, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 5.191};

  const int nptBin = 5;
  double avgPtBin[]={50.0,90.0,120.0,170.0,250.0,1000.0};
  TH1F* hRecoPt[nEta][nptBin];
  TH1F* hGenPt[nEta][nptBin];
  TH1F* hGenEta_RecoEta[nEta][nptBin];
  TH1F* hGenPhi_RecoPhi[nEta][nptBin];

  for (int i = 0; i < nEta; i++) {
    for (int ipt = 0; ipt <nptBin; ipt++) {
      hRecoPt[i][ipt] = new TH1F(Form("hRecoPt_etaBin_%d_%d", i,ipt), Form("pT_Reco_in_Eta_Bin_%.1f_to_%.1f", etaBin[i], etaBin[i+1]), 1000,0,1000);//1000,avgPtBin);
      hRecoPt[i][ipt] ->Sumw2();
      hGenPt[i][ipt] = new TH1F(Form("hGenpT_etaBin_%d_%d",i,ipt), Form("pT_Gen_in_Eta_Bin_%.1f_to_%.1f", etaBin[i], etaBin[i+1]), 1000,0,1000);//1000,avgPtBin);//1000,0,1000);
      hGenPt[i][ipt] ->Sumw2();
      hGenPhi_RecoPhi[i][ipt] = new TH1F(Form("hGenPhi_RecoPhi_etaBin_%d_%d", i,ipt), Form("RecoPhi_GenPhi_difference_in_Eta_Bin_%.1f_to_%.1f", etaBin[i], etaBin[i+1]), 400, -2, 2);
      hGenPhi_RecoPhi[i][ipt] ->Sumw2();
      hGenEta_RecoEta[i][ipt] = new TH1F(Form("hGenEta_RecoEta_etaBin_%d_%d", i,ipt), Form("RecoEta_GenEta_difference_in_Eta_Bin_%.1f_to_%.1f", etaBin[i], etaBin[i+1]), 400, -2, 2);
      hGenEta_RecoEta[i][ipt] ->Sumw2();
    }
  }
  // const int nEta=4;
  float nRatio13;
  float nRatio1325;
  float nRatio2530;
  float nRatio3050;
 

  // Event Loop
  bool alternate=true;
  for (Long64_t i = 0; i < nentries; i++) {
    //for (Long64_t i = 0; i < 1000; i++) {
    //for (Long64_t i = 0; i < 500000; i++) {
    //for (Long64_t i = 500000; i < nentries; i++) {
    tree->GetEntry(i);

    float ptResponse;
    float maxRawPt_Gen= -9999;
    float secondmaxRawPt_Gen= -9999;
    float maxRawPt_Reco= -9999;
    float secondmaxRawPt_Reco= -9999;
  
    
    for (int j = 0; j < nref; j++) {
      
      float pt_gen = (*refpt)[j];
      float eta_gen = (*refeta)[j];
      float phi_gen = (*refphi)[j];
      float pt_raw = (*jtpt)[j];
      float eta_raw = (*jteta)[j];
      float phi_raw = (*jtphi)[j];
      //if(pt_gen>maxRawPt_Gen)
      if(1)
	{
	  maxRawPt_Gen=pt_gen;
	}
      
      JECCor.SetJetPT(pt_raw);
      JECCor.SetJetEta(eta_raw);
      JECCor.SetJetPhi(phi_raw);
      //cout<<i<<'\t'<<j<<'\t'<<"raw: " << pt <<'\t'<< eta << '\t'<< phi <<endl;
      float pt = JECCor.GetCorrectedPT();
      float eta =eta_raw;
      float phi = phi_raw;
      //if(pt_raw>maxRawPt_Reco)
      if(1)
	{
	  maxRawPt_Reco=pt;
 	}

      for (int ib = 0; ib < nEta; ib++) {
	for(int ipt=0;ipt<nptBin; ipt++){
	  if(eta >= etaBin[ib] && eta < etaBin[ib+1]){
	    if(maxRawPt_Gen >= avgPtBin[ipt] && maxRawPt_Gen < avgPtBin[ipt+1]){
	      hGenPhi_RecoPhi[ib][ipt]->Fill(phi-phi_gen);
	      hGenPhi_RecoPhi[ib][ipt]->SetTitle("MC_gen_reco_#phi_differnece;#phi_{Reco}-#phi_{Gen};");
	      hGenEta_RecoEta[ib][ipt]->Fill(eta-eta_gen);
	      hGenEta_RecoEta[ib][ipt]->SetTitle("MC_gen_reco_#eta_differnece;#eta_{Reco}-#eta_{Gen};");

	      hRecoPt[ib][ipt]->Fill(maxRawPt_Reco);
	      hRecoPt[ib][ipt]->SetTitle("MC_Reco_pT;p_{T}^{Reco};");

	      hGenPt[ib][ipt]->Fill(maxRawPt_Gen);
	      hGenPt[ib][ipt]->SetTitle("MC_gen_pT;p_{T}^{Gen};");

	    }
	  }
	}
      }   
    }
  }

 
#if 1

  // Save Output
  TFile *outputFile = new TFile(Form("/afs/cern.ch/user/u/uacharya/eos/cms_analaysis/JEC_analysis/simuMC_output2/JRA_output_%d.root",infilInd), "RECREATE");

  for (int ib = 0; ib < nEta; ib++) {
    for(int ipt=0;ipt<nptBin; ipt++){

      hGenPhi_RecoPhi[ib][ipt]->Write();
      hGenEta_RecoEta[ib][ipt]->Write();
      // hRecoPt[ib][ipt]->Write();
      // hGenPt[ib][ipt]->Write();
    }
  }
  outputFile->Close();
#endif

  
}
