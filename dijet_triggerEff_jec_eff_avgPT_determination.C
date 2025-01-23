//void createHistos_fromminAOD(int infilInd=1,int inDir=2){
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

void dijet_triggerEff_jec_eff_avgPT_determination(int infilInd=1){


  //////////////////////////////////////////////////////////
  //   This file has been automatically generated 
  //     (Fri Oct  4 18:23:21 2024 by ROOT version6.32.04)
  //   from TTree t/ Jet Analysis Tree
  //   found on file: root://eoshome-u.cern.ch//eos/user/u/uacharya/phys_heavyions/2024ppref_hiforest/Forest/JetMET0/Run2024H-PromptReco-v1/MINIAOD/JetMET0/pp_ref_lowpileup_20241003/241003_184418/0000/HiForestMiniAOD_1.root
  //////////////////////////////////////////////////////////

  using namespace std;
  //Reset ROOT and connect tree file
  gROOT->Reset();
  //JetMet0_V2//
  //TFile *HiFor_file = new TFile(Form("/eos/cms/store/group/phys_heavyions/uacharya/JetMet0_v2/HiForestMiniAOD_%d.root",infilInd));


  //JetMet1_V2//

  //TFile *HiFor_file = new TFile(Form("/eos/cms/store/group/phys_heavyions/uacharya/JetMet1_v2/HiForestMiniAOD_%d.root",infilInd));
  TFile *HiFor_file = new TFile(Form("/afs/cern.ch/user/u/uacharya/eos/phys_heavyions/2024ppref_hiforest/ZeroBias_Sample/SpecialZeroBias5/Run2024I-PromptReco-v2/MINIAOD/SpecialZeroBias5/pp_ref_lowpileup_SpecialZerobias5_V2/241120_180213/0000/HiForestMiniAOD_%d.root",infilInd));

   
  // TDirectory * dir = (TDirectory*)HiFor_file->Get(Form("ak%dPFJetAnalyzer",inDir));
  TDirectory * dir = (TDirectory*)HiFor_file->Get(Form("ak4PFJetAnalyzer"));
  if (!dir) {
    cout << "Directory ak4PFJetAnalyzer not found!" << endl;
    return;
  }
   
  TTree *tree;
  dir->GetObject("t",tree);

  if (!tree) {
    cout << "Tree t not found!" << endl;
  } else {
    cout << "Tree successfully loaded!" << endl;
  }
  TDirectory * dirTrigg = (TDirectory*)HiFor_file->Get(Form("hltanalysis"));
  if (!dirTrigg) {
    cout << "Directory hltanalysis not found!" << endl;
    return;
  }
    
  TTree *treeTrigg;
  dirTrigg->GetObject("HltTree",treeTrigg);

  if (!treeTrigg) {
    cout << "Tree HltTree not found!" << endl;
  } else {
    cout << "Tree successfully loaded!" << endl;
  }


  vector<string> JEC_input_file;
  JEC_input_file.push_back("/afs/cern.ch/user/u/uacharya/eos/cms_analaysis/JEC_analysis/JRA_step3harvest_L2Relative_AK4PF.txt");
  //JetCorrector::JetCorrector;
  JetCorrector JECCor(JEC_input_file);
  
  //tree->Print();
  const int mx=100000;
  //Declaration of leaves types
  Int_t           HLT_ZeroBias_v12;
  Int_t           HLT_PFJet110_v11;
  Int_t           HLT_PFJet40_ZeroBiasCopy_v1;
  Int_t           HLT_PFJet80_L1SingleJet60_v1;
  Int_t           HLT_SpecialZeroBias_v5;
  Int_t           L1_ZeroBias_copy;
  Int_t           run;
  Int_t           evt;
  Int_t           lumi;
  Int_t           nref;
  Int_t           ncalo;
  Float_t         rawpt[mx];
  Float_t         jtpt[mx];
  Float_t         jteta[mx];
  Float_t         jty[mx];
  Float_t         jtphi[mx];
  Float_t         jtpu[432];
  Float_t         jtm[432];
  Float_t         jtarea[432];
  //Int_t           ncalo;
  Float_t         calopt[88];
  Float_t         caloeta[88];
  Float_t         calophi[88];
  Float_t         jtPfCHF[mx];
  Float_t         jtPfNHF[mx];
  Float_t         jtPfCEF[mx];
  Float_t         jtPfNEF[mx];
  Float_t         jtPfMUF[mx];
  Int_t           jtPfCHM[mx];
  Int_t           jtPfNHM[mx];
  Int_t           jtPfCEM[mx];
  Int_t           jtPfNEM[mx];
  Int_t           jtPfMUM[mx];
  Float_t         jttau1[432];
  Float_t         jttau2[432];
  Float_t         jttau3[432];
  Float_t         discr_deepCSV[432];
  Float_t         discr_pfJP[432];
   
  // Set branch addresses.
  treeTrigg->SetBranchAddress("HLT_ZeroBias_v12",&HLT_ZeroBias_v12);
  treeTrigg->SetBranchAddress("HLT_PFJet110_v11",&HLT_PFJet110_v11);
  treeTrigg->SetBranchAddress("HLT_PFJet40_ZeroBiasCopy_v1",&HLT_PFJet40_ZeroBiasCopy_v1);
  treeTrigg->SetBranchAddress("HLT_PFJet80_L1SingleJet60_v1",&HLT_PFJet80_L1SingleJet60_v1);
  treeTrigg->SetBranchAddress("HLT_SpecialZeroBias_v5",&HLT_SpecialZeroBias_v5);
  treeTrigg->SetBranchAddress("L1_ZeroBias_copy",&L1_ZeroBias_copy);

  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("evt",&evt);
  tree->SetBranchAddress("lumi",&lumi);
  tree->SetBranchAddress("nref",&nref);
  tree->SetBranchAddress("ncalo",&ncalo);
  tree->SetBranchAddress("rawpt",rawpt);
  tree->SetBranchAddress("jtpt",jtpt);
  tree->SetBranchAddress("jteta",jteta);
  tree->SetBranchAddress("jty",jty);
  tree->SetBranchAddress("jtphi",jtphi);
  tree->SetBranchAddress("jtpu",jtpu);
  tree->SetBranchAddress("jtm",jtm);
  tree->SetBranchAddress("jtarea",jtarea);
  //tree->SetBranchAddress("ncalo",&ncalo);
  tree->SetBranchAddress("calopt",calopt);
  tree->SetBranchAddress("caloeta",caloeta);
  tree->SetBranchAddress("calophi",calophi);
  tree->SetBranchAddress("jtPfCHF",jtPfCHF);
  tree->SetBranchAddress("jtPfNHF",jtPfNHF);
  tree->SetBranchAddress("jtPfCEF",jtPfCEF);
  tree->SetBranchAddress("jtPfNEF",jtPfNEF);
  tree->SetBranchAddress("jtPfMUF",jtPfMUF);
  tree->SetBranchAddress("jtPfCHM",jtPfCHM);
  tree->SetBranchAddress("jtPfNHM",jtPfNHM);
  tree->SetBranchAddress("jtPfCEM",jtPfCEM);
  tree->SetBranchAddress("jtPfNEM",jtPfNEM);
  tree->SetBranchAddress("jtPfMUM",jtPfMUM);
  tree->SetBranchAddress("jttau1",jttau1);
  tree->SetBranchAddress("jttau2",jttau2);
  tree->SetBranchAddress("jttau3",jttau3);
#if 1

  /*
  //     This is the loop skeleton
  //       To read only selected branches, Insert statements like:
  // tree->SetBranchStatus("*",0);  // disable all branches
  // TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname
  */
  //Long64_t nentries = tree->GetEntries();
  Long64_t nbytes = 0;
  Long64_t nentries = treeTrigg->GetEntries();

  TH1F *hrawPt_unCor = new TH1F("hrawPt_unCor","hrawPt_unCor;Raw_Leading_P_{T};",1000,0,500);
  TH1F *hrawPt40Zero= new TH1F("hrawPt40Zero","hrawPt40Zero;Corrected_Leading_P_{T};",1000,0,500);
  hrawPt40Zero->Sumw2();
  TH1F *hrawPt_trig40Zero = new TH1F("hrawPt_trig40Zero","hrawPt_trig40Zero;Corrected_Leading_P_{T};",1000,0,500);
  hrawPt_trig40Zero->Sumw2();
  TH1F *hrawPt110= new TH1F("hrawPt110","hrawPt110;Corrected_Leading_P_{T};",1000,0,500);
  hrawPt110->Sumw2();
  TH1F *hrawPt_trig110 = new TH1F("hrawPt_trig110","hrawPt_trig110;Corrected_Leading_P_{T};",1000,0,500);
  hrawPt_trig110->Sumw2();
  TH1F *hrawPt4080= new TH1F("hrawPt4080","hrawPt4080;Corrected_Leading_P_{T};",1000,0,500);
  hrawPt4080->Sumw2();
  TH1F *hrawPt_trig4080 = new TH1F("hrawPt_trig4080","hrawPt_trig4080;Corrected_Leading_P_{T};",1000,0,500);
  hrawPt_trig4080->Sumw2();
  Float_t jtPf_CHF,jtPf_MUF,jtPf_NHF,jtPf_NEF;
  Int_t jtPf_CHM;
  //for (Long64_t i=0; i<500;i++) {
  for (Long64_t i=0; i<nentries;i++) {
    nbytes += tree->GetEntry(i);
    treeTrigg ->GetEntry(i);
    // cout << "Processing file " << " entry " << i + 1 <<" "<<"run: "<< run<<endl;
    //if(run==386749||run==386753){
    float maxRawPt= -9999;
    float secondmaxRawPt= -9999;
    float maxRawPt_unCor= -9999;
    float maxRawPt_trig;
    Float_t rawPt,jtPt,jtEta,jtPhi,jtY;
    Float_t Corrected_rawPT;
    Float_t Correction_rawPT;
    float avgPt;
    bool trig110Condition=HLT_PFJet110_v11==1;
    bool trig40Condition=HLT_PFJet40_ZeroBiasCopy_v1==1;
    bool trig80Condition=HLT_PFJet80_L1SingleJet60_v1==1;
    // Float_t lead_pt=rawpt[0];
    for(int j=0;j<nref;j++)
      {
	jtEta=jteta[j];
	jtPhi=jtphi[j];
	jtPf_CHF=jtPfCHF[j];
	jtPf_MUF=jtPfMUF[j];
	jtPf_NHF=jtPfNHF[j];
	jtPf_CHM=jtPfCHM[j];
	jtPf_NEF=jtPfNEF[j];
	 
	if((TMath::Abs(jtEta)<2) && jtPf_MUF<0.001 && ( jtPf_CHF > 0.01 && jtPf_CHF < 0.99) && jtPf_CHM >0 && jtPf_NHF < 0.9 && jtPf_NEF <0.9)
	  //if( jtPf_MUF<0.001 && ( jtPf_CHF > 0.01 && jtPf_CHF < 0.99) && jtPf_CHM >0 && jtPf_NHF < 0.9 && jtPf_NEF <0.9)
	  //if(1)
	  {
	    rawPt=rawpt[j];
	    JECCor.SetJetPT(rawPt);
	    JECCor.SetJetEta(jtEta);
	    JECCor.SetJetPhi(jtPhi);
	    Corrected_rawPT=JECCor.GetCorrectedPT();
	    Correction_rawPT=JECCor.GetCorrection();
	    //This is only for the cross check to see if the jetcorrection is working or not	     
	    if(1)
	      //if(HLT_PFJet40_ZeroBiasCopy_v1==1 || HLT_PFJet80_L1SingleJet60_v1==1)
	      {
		if(rawPt>maxRawPt_unCor)
		  {
		    maxRawPt_unCor=rawPt;
		    //if(HLT_PFJet40_ZeroBiasCopy_v1==1 || HLT_PFJet80_L1SingleJet60_v1==1){cout<<maxRawPt_unCor<<endl;}

		  }
	      }
	    /////////////////////////////////////////////////////
	    /////////////////////////////////////////////
	    if(1)
	      //if(HLT_PFJet40_ZeroBiasCopy_v1==1 || HLT_PFJet80_L1SingleJet60_v1==1)
	      {
		if(Corrected_rawPT>maxRawPt)
		  {
		    secondmaxRawPt=maxRawPt;
		    maxRawPt=Corrected_rawPT;
		  }
		else if(Corrected_rawPT>secondmaxRawPt)
		  {
		    secondmaxRawPt=Corrected_rawPT;
		  }
		avgPt=0.5*(maxRawPt+secondmaxRawPt);

	      }

	    
	    
	  }
      }

    if(L1_ZeroBias_copy==1)//||HLT_SpecialZeroBias_v5)
     
      {
	//cout<<" HLT_deno: "<<maxRawPt<<'\t'<<secondmaxRawPt<<'\t'<<avgPt<<endl;
	hrawPt40Zero->Fill(avgPt);
	hrawPt40Zero->SetTitle(";p_{T}^{Avg}; Events with L1_ZeroBias_copy_trigger");

      }

   
    if((HLT_PFJet40_ZeroBiasCopy_v1==1) && (L1_ZeroBias_copy==1))//||HLT_SpecialZeroBias_v5))
      {
	//cout<<" HLT_nume: "<<maxRawPt<<'\t'<<secondmaxRawPt<<'\t'<<avgPt<<endl;
	hrawPt_trig40Zero->Fill(avgPt);
	hrawPt_trig40Zero->SetTitle(";p_{T}^{Avg}; Events with HLT_PFJet40_ZeroBiasCopy_v1_Trigger_and L1_ZeroBias_copy_trigger");
      }

    
    
    if(HLT_PFJet40_ZeroBiasCopy_v1==1 || HLT_PFJet80_L1SingleJet60_v1==1)
      {
	// avgPt=0.5*(maxRawPt+secondmaxRawPt);
	cout<<" HLT_deno: "<<maxRawPt<<'\t'<<secondmaxRawPt<<'\t'<<avgPt<<endl;
	// if(avgPt<=0)continue;
	hrawPt110->Fill(avgPt);
	hrawPt110->SetTitle(";p_{T}^{Avg}; Events with HLT_PFJet40_0r_80_Trigger");

      }
    //cout<<i<<" "<<"jtEta: "<<jtEta<<" " <<maxRawPt_unCor<<" "<<"maxRawPt: "<< maxRawPt<<'\t'<<secondmaxRawPt<<endl;

    //hrawPt_unCor->Fill(maxRawPt_unCor);
   
    //if(1)
    //if(HLT_PFJet110_v11==1 &&(HLT_PFJet40_ZeroBiasCopy_v1==1 || HLT_PFJet80_L1SingleJet60_v1==1))
    if(HLT_PFJet110_v11==1 &&(HLT_PFJet40_ZeroBiasCopy_v1==1 || HLT_PFJet80_L1SingleJet60_v1==1))
      {
	// avgPt=0.5*(maxRawPt+secondmaxRawPt);
	// if(avgPt<=0)continue;
	cout<<" HLT_nume: "<<maxRawPt<<'\t'<<secondmaxRawPt<<'\t'<<avgPt<<endl;
	hrawPt_trig110->Fill(avgPt);
	hrawPt_trig110->SetTitle(";p_{T}^{Avg}; Events with HLT_PFJet40_0r_80_Trigger_and PFJet110_Trigger");
	//cout<<"hlttrig :"<<maxRawPt<<endl;
      }


    if(HLT_PFJet40_ZeroBiasCopy_v1==1)
      {
	cout<<" HLT_deno: "<<maxRawPt<<'\t'<<secondmaxRawPt<<'\t'<<avgPt<<endl;
	hrawPt4080->Fill(avgPt);
	hrawPt4080->SetTitle(";p_{T}^{Avg}; Events with HLT_PFJet40_0r_80_Trigger");

      }

   
    if(HLT_PFJet40_ZeroBiasCopy_v1==1 && HLT_PFJet80_L1SingleJet60_v1==1)
      {
	cout<<" HLT_nume: "<<maxRawPt<<'\t'<<secondmaxRawPt<<'\t'<<avgPt<<endl;
	hrawPt_trig4080->Fill(avgPt);
	hrawPt_trig4080->SetTitle(";p_{T}^{Avg}; Events with HLT_PFJet40_and_80_Trigger");
      }

  
  }
  TH1F *hrawPt_ratio40Zero = (TH1F*)hrawPt_trig40Zero->Clone("hrawPt_ratio40Zero");
  hrawPt_ratio40Zero->SetTitle("Ratio of Trigger Passed to All Events; Corrected_leading P_{T}; Ratio");
  hrawPt_ratio40Zero->Divide(hrawPt40Zero);
  hrawPt_ratio40Zero->SetTitle(";p_{T}^{avg}; #frac{Events_with_L1_ZeroBias_and_PFJet40_trig}{Events_with_L1_ZeroBias_Trig}");

  TH1F *hrawPt_ratio110 = (TH1F*)hrawPt_trig110->Clone("hrawPt_ratio110");
  hrawPt_ratio110->SetTitle("Ratio of Trigger Passed to All Events; Corrected_leading P_{T}; Ratio");
  hrawPt_ratio110->Divide(hrawPt110);
  hrawPt_ratio110->SetTitle(";p_{T}^{avg}; #frac{Events with HLT_PFJet40_0r_80_and_110_Trig}{Events with HLT_PFJet40_or80_Trig}");
  
  TH1F *hrawPt_ratio4080 = (TH1F*)hrawPt_trig4080->Clone("hrawPt_ratio4080");
  hrawPt_ratio4080->SetTitle("Ratio of Trigger Passed to All Events; Corrected_leading P_{T}; Ratio");
  hrawPt_ratio4080->Divide(hrawPt4080);
  hrawPt_ratio4080->SetTitle(";p_{T}^{avg}; #frac{Events with HLT_PFJet40_and_80_and_110_Trig}{Events with HLT_PFJet40}");

  TFile *outputFile= new TFile(Form("/afs/cern.ch/user/u/uacharya/eos/cms_analaysis/JEC_analysis/zerobias_miniAOD/HiForestMiniAOD_zerobias_output_%d.root",infilInd),"RECREATE");

  //TFile *outputFile= new TFile(Form("/afs/cern.ch/user/u/uacharya/eos/cms_analaysis/JEC_analysis/histo_eff_JEC/HiForestMiniAOD_JM1_v2_output_%d_4.root",infilInd),"RECREATE");
  // TFile *outputFile= new TFile(Form("/afs/cern.ch/user/u/uacharya/eos/cms_analaysis/JEC_analysis/HiForestMiniAOD_JM0_v2_output_%d_4.root",infilInd),"RECREATE");
  //hrawPt_unCor->Write();
  hrawPt40Zero->Write();
  hrawPt_trig40Zero->Write();
  hrawPt_ratio40Zero->Write();
  hrawPt110->Write();
  hrawPt_trig110->Write();
  hrawPt_ratio110->Write();
  hrawPt4080->Write();
  hrawPt_trig4080->Write();
  hrawPt_ratio4080->Write();
  
  //outputFile->Write();
  //gSystem->Sleep(100);
  outputFile->Close();
  

#endif
}
  
