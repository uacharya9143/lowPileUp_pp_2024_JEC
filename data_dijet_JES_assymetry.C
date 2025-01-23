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
//#include "JetCorrector.h"
#include "/afs/cern.ch/user/u/uacharya/eos/cms_analaysis/JEC_analysis/JetCorrector.h"

void data_dijet_imbalance_asymmetry_L2Residual_JM0(int infilInd=1) {

  using namespace std;
  gROOT->Reset();

  // File and Histogram Setup
 
  //JetMet0_V2//
  TFile *HiFor_file = new TFile(Form("/eos/cms/store/group/phys_heavyions/uacharya/JetMet0_v2/HiForestMiniAOD_%d.root", infilInd));

  //JetMet1_V2//
  //TFile *HiFor_file = new TFile(Form("/eos/cms/store/group/phys_heavyions/uacharya/JetMet1_v2/HiForestMiniAOD_%d.root",infilInd));

  
  TDirectory *dir = (TDirectory*)HiFor_file->Get(Form("ak4PFJetAnalyzer"));
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
  //JEC_input_file.push_back("/afs/cern.ch/user/u/uacharya/eos/cms_analaysis/JEC_analysis/JRA_step3harvest_L2Relative_AK4PF.txt");
  JEC_input_file.push_back("/afs/cern.ch/user/u/uacharya/MC_Truth_JEC/CMSSW_14_1_0/src/MC-truth-JEC/condor/Files/ParallelMCL1_L2Relative_AK4PF.txt");
  ///afs/cern.ch/user/u/uacharya/MC_Truth_JEC/CMSSW_14_1_0/src/MC-truth-JEC/condor/Files/ParallelMCL1_L2Relative_AK4PF.txt");
  JetCorrector JECCor(JEC_input_file);

  const int mx = 100000;
  Int_t nref;
  Float_t rawpt[mx], jteta[mx], jtphi[mx];
  Int_t           HLT_PFJet110_v11;
  Int_t           HLT_PFJet40_ZeroBiasCopy_v1;
  Int_t           HLT_PFJet80_L1SingleJet60_v1;
  Int_t           run;
  Int_t           evt;
  Int_t           lumi;
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
  treeTrigg->SetBranchAddress("HLT_PFJet110_v11",&HLT_PFJet110_v11);
  treeTrigg->SetBranchAddress("HLT_PFJet40_ZeroBiasCopy_v1",&HLT_PFJet40_ZeroBiasCopy_v1);
  treeTrigg->SetBranchAddress("HLT_PFJet80_L1SingleJet60_v1",&HLT_PFJet80_L1SingleJet60_v1);
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("evt",&evt);
  tree->SetBranchAddress("lumi",&lumi);
  tree->SetBranchAddress("nref", &nref);
  tree->SetBranchAddress("rawpt", rawpt);
  tree->SetBranchAddress("jteta", jteta);
  tree->SetBranchAddress("jtphi", jtphi);
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
  
  // int nEta = 22;
  // float etaBin[] = {-5.0, -3.5, -3.0, -2.7, -2.4, -2.0, -1.7, -1.3, -1.0, -0.7, -0.4, 0.0, 0.4, 0.7, 1.0, 1.3, 1.7, 2.0, 2.4, 2.7, 3.0, 3.5, 5.0};
  // const int nptBin = 5;
  // double avgPtBin[]={50.0,90.0,120.0,160.0,200.0,500.0};

  // int nEta = 36;
 //  float etaBin[] = {-5.191, -3.839, -3.489, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -1.93, -1.653, -1.479, -1.305, -1.044, -0.783, -0.522, -0.261, 0.0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
 // const int nptBin = 5;
 //  double avgPtBin[]={50.0,90.0,120.0,170.0,250.0,1000.0};

  // int nEta = 60;

  // float etaBin[] = {-5.191, -4.515, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.051, -1.93, -1.791, -1.653, -1.566, -1.479, -1.392, -1.305, -1.174, -1.044, -0.913, -0.783, -0.652, -0.522, -0.391, -0.261, -0.130, 0.0, 0.130, 0.261, 0.391, 0.522, 0.652, 0.783, 0.913, 1.044, 1.174, 1.305, 1.392, 1.479, 1.566, 1.653, 1.791, 1.93, 2.051, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.515, 5.191};

  const int nEta = 50;
  float etaBin[] = {-5.191,  -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.051, -1.93, -1.791, -1.653, -1.566, -1.479, -1.392, -1.305, -1.174, -1.044, -0.913, -0.783, -0.652, -0.522, -0.391, -0.261, -0.130, 0.0, 0.130, 0.261, 0.391, 0.522, 0.652, 0.783, 0.913, 1.044, 1.174, 1.305, 1.392, 1.479, 1.566, 1.653, 1.791, 1.93, 2.051, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 5.191};

  const int nptBin = 5;
  double avgPtBin[]={50.0,90.0,120.0,170.0,250.0,1000.0};

  int eventCount[]={0};
  
  TH1F *hPt_eta = new TH1F("hPt_eta", "Leading Jet Pt by Eta Bin", nEta, etaBin);
  TH1F *hPt_eta13 = new TH1F("hPt_eta13", "Second Jet Pt by Eta Bin", nEta, etaBin);
  TH1F *hPt_eta_third = new TH1F("hPt_eta_third", "Third Jet Pt by Eta Bin", nEta, etaBin);
  TH1F *hPt_eta_ratio;
  TH1F *heventCount = new TH1F("heventCount", "Event Count by Eta Bin", nEta, etaBin);
  TH1F *hPt_avg = new TH1F("hPt_avg", "p_{T} avg by Eta Bin", nEta, etaBin);
  TH1F *h_alpha = new TH1F("h_alpha", "#alpha values by Eta Bin", nEta, etaBin);
  TH1F *heventCount_maxPt = new TH1F("heventCount_maxPt", "Event Count for Leading Jet (maxPt) by Eta Bin", nEta, etaBin);
  TH1F *heventCount_secondMaxPt = new TH1F("heventCount_secondMaxPt", "Event Count for Second Jet (secondMaxPt) by Eta Bin", nEta, etaBin);

  Long64_t nentries = tree->GetEntries();
  Float_t jtPf_CHF,jtPf_MUF,jtPf_NHF,jtPf_NEF, jtPf_CEF;
  Int_t jtPf_CHM;

 
  //TH1F* hAssymBins[nEta];
  TH1F* hAssymBins[nptBin][nEta];
  TH1F* hAssymBins40[nptBin][nEta];
  TH1F* hAssymBins80[nptBin][nEta];
  TH1F* hAssymBins110[nptBin][nEta];
  TH1F* hPt_Eta[nEta];
  TH1F* hPt_Eta13[nEta];

  
  for (int ipt = 0; ipt <nptBin; ipt++) {
    for (int i = 0; i < nEta; i++) {
      hAssymBins[ipt][i] = new TH1F(Form("hAssym_etaBin_%d_%d", ipt,i), Form("Asymmetry_in_Eta_Bin_%.3f_to_%.3f", etaBin[i], etaBin[i+1]), 200, -1, 1);
      hAssymBins[ipt][i] ->Sumw2();
      
      hAssymBins40[ipt][i] = new TH1F(Form("hAssym_etaBin40_%d_%d", ipt,i), Form("Asymmetry in Eta Bin %.3f to %.3f", etaBin[i], etaBin[i+1]), 200, -1, 1);
      hAssymBins40[ipt][i] ->Sumw2();
      hAssymBins80[ipt][i] = new TH1F(Form("hAssym_etaBin80_%d_%d", ipt,i), Form("Asymmetry in Eta Bin %.3f to %.3f", etaBin[i], etaBin[i+1]), 200, -1, 1);
      hAssymBins80[ipt][i] ->Sumw2();
      hAssymBins110[ipt][i] = new TH1F(Form("hAssym_etaBin110_%d_%d", ipt,i), Form("Asymmetry in Eta Bin %.3f to %.3f", etaBin[i], etaBin[i+1]), 200, -1, 1);
      hAssymBins110[ipt][i] ->Sumw2();
    }
  }
  
  for (int i = 0; i < nEta; i++) {
    hPt_Eta[i] = new TH1F(Form("hPt_Eta_%d", i), Form("MaxPt Eta Bin %.1f to %.1f", etaBin[i], etaBin[i+1]), 1000, 0, 1000);
    hPt_Eta13[i] = new TH1F(Form("hPt_Eta13_%d", i), Form("SecondMaxPt in Eta Bin %.1f to %.1f", etaBin[i], etaBin[i+1]), 1000, 0, 1000);
  }
  TH1F *hResponse = new TH1F("hResponse", "Response by Eta Bin", nEta, etaBin);
  float asymmetryMeans[] = {0};
  // Event Loop
  TH1F *hrawPt_unCor = new TH1F("hrawPt_unCor","hrawPt_unCor;Raw_Leading_P_{T};",1000,0,1000);
  TH1F *hrawPt = new TH1F("hrawPt","hrawPt;Corrected_Leading_P_{T};",1000,0,1000);
  
  bool alternate1= true;
  bool alternate2= true;
  bool alternate3= true;
   
  for (Long64_t i = 0; i < nentries; i++) {
    //for (Long64_t i = 0; i < 1000; i++) {
    tree->GetEntry(i);
    treeTrigg ->GetEntry(i);
    float maxPt = -9999, maxEta = -9999, maxPhi = -9999;
    float secondMaxPt = -9999, secondMaxEta = -9999, secondMaxPhi = -9999;
    float thirdMaxPt = -9999, thirdMaxEta = -9999, thirdMaxPhi = -9999;
    float assym;
    float pt_cor[nref];
    float eta_cor[nref];
    float phi_cor[nref];
    float avg_pT;
    float alpha_Val;

    bool trig110Condition=(HLT_PFJet110_v11==1);
    bool trig40Condition=(HLT_PFJet40_ZeroBiasCopy_v1==1);
    bool trig80Condition=(HLT_PFJet80_L1SingleJet60_v1==1);
 

    // Leading Jet Selection: Identify leading jet with highest pt within the eta bins
    for (int j = 0; j < nref; j++) {

      jtPf_CHF=jtPfCHF[j];
      jtPf_MUF=jtPfMUF[j];
      jtPf_NHF=jtPfNHF[j];
      jtPf_CHM=jtPfCHM[j];
      jtPf_NEF=jtPfNEF[j];
      jtPf_CEF=jtPfCEF[j];
      //if( jtPf_MUF<0.001 && jtPf_CHF > 0.01 && jtPf_CHF < 0.99  && jtPf_NHF < 0.9 && jtPf_NEF <0.9) {  
      if( jtPf_MUF<0.80 && jtPf_CHF > 0.01  && jtPf_NHF < 0.99 && jtPf_NEF <0.9 && jtPf_CEF<0.80 && jtPf_CHM>0) {
	
	JECCor.SetJetPT(rawpt[j]);
	JECCor.SetJetEta(jteta[j]);
	JECCor.SetJetPhi(jtphi[j]);
	//float pt1 = rawpt[j];
	float pt = JECCor.GetCorrectedPT();
	pt_cor[j]=pt;
	float eta = jteta[j];
	eta_cor[j]=eta;
	float phi = jtphi[j];
	phi_cor[j]=phi;
	//cout<<i<<'\t'<<j<<'\t'<<rawpt[j]<<'\t'<<pt<<endl;
	hrawPt_unCor->Fill(rawpt[j]);
	hrawPt->Fill(pt_cor[j]);

	if (pt > maxPt) {
	  thirdMaxPt = secondMaxPt;
	  thirdMaxEta = secondMaxEta;
	  thirdMaxPhi = secondMaxPhi;
	  //cout<<i<<'\t'<<j<<'\t'<<"third1: " << thirdMaxPt <<'\t'<< thirdMaxEta << '\t'<< thirdMaxPhi <<endl;

	  secondMaxPt = maxPt;
	  secondMaxEta = maxEta;
	  secondMaxPhi = maxPhi;
	  //cout<<i<<'\t'<<j<<'\t'<<"second1: " << secondMaxPt <<'\t'<< secondMaxEta << '\t'<< secondMaxPhi <<endl;

	  maxPt = pt;
	  maxEta = eta;
	  maxPhi = phi;
	  //cout<<i<<'\t'<<j<<'\t'<<"first: " << pt <<'\t'<< eta << '\t'<< phi <<endl;

	} else if (pt > secondMaxPt) {
	  // Update third jet before changing secondMaxPt
	  thirdMaxPt = secondMaxPt;
	  thirdMaxEta = secondMaxEta;
	  thirdMaxPhi = secondMaxPhi;

	  secondMaxPt = pt;
	  secondMaxEta = eta;
	  secondMaxPhi = phi;

	} else if (pt > thirdMaxPt) {
	  // Only update third jet
	  thirdMaxPt = pt;
	  thirdMaxEta = eta;
	  thirdMaxPhi = phi;

	}
      }
        
    }
 
  
    //if(HLT_PFJet110_v11==1 &&(HLT_PFJet40_ZeroBiasCopy_v1==1 || HLT_PFJet80_L1SingleJet60_v1==1))
    float alphaValue=0.35;

#if 1
    if(trig40Condition ){
      for (int ib = 0; ib < nEta; ib++) {
	
	if( fabs(maxEta)<=1.3 &&fabs(secondMaxEta)>1.3){
	  float delPhi = TMath::Abs(maxPhi - secondMaxPhi);
	  bool phiCondition = (delPhi >2.7);
	  if(!phiCondition) continue;
	  avg_pT=0.5*(maxPt+secondMaxPt);
	  alpha_Val=thirdMaxPt/avg_pT;
	  bool alphaCondition1=(alpha_Val<alphaValue);
	  if((avg_pT>=avgPtBin[0] &&avg_pT<avgPtBin[1]) && alphaCondition1){
	    if(secondMaxEta >= etaBin[ib] && secondMaxEta < etaBin[ib+1]) {
	      assym=(secondMaxPt-maxPt)/(secondMaxPt+maxPt);
	      if(fabs(assym)>0.6) continue;	      
	      hAssymBins[0][ib]->Fill(assym);
	    }
	  }
	}
	else if(fabs(secondMaxEta)<=1.3 &&fabs(maxEta)>1.3){
	  float delPhi = TMath::Abs(maxPhi - secondMaxPhi);
	  bool phiCondition = (delPhi >2.7);
	  if(!phiCondition) continue;
	  avg_pT=0.5*(maxPt+secondMaxPt);
	  alpha_Val=thirdMaxPt/avg_pT;
	  bool alphaCondition2=(alpha_Val<alphaValue);
	  if((avg_pT>=avgPtBin[0] &&avg_pT<avgPtBin[1]) && alphaCondition2){
	    if( maxEta>= etaBin[ib] &&  maxEta< etaBin[ib+1]) {
	      //if(assym>0.55) {cout<<i<<'\t'<<ib<<'\t'<<assym<<endl;}

	      assym=(maxPt-secondMaxPt)/(maxPt+secondMaxPt);
	      if(fabs(assym)>0.6) continue;
	      hAssymBins[0][ib]->Fill(assym);
	    }
	  }
	}
	else if((fabs(secondMaxEta)<=1.3 &&fabs(maxEta)<=1.3)){
	  float delPhi = TMath::Abs(maxPhi - secondMaxPhi);
	  bool phiCondition = (delPhi >2.7);
	  if(!phiCondition) continue;
	  avg_pT=0.5*(maxPt+secondMaxPt);
	  alpha_Val=thirdMaxPt/avg_pT;
	  bool alphaCondition3=(alpha_Val<alphaValue);
	  if(avg_pT>=avgPtBin[0] &&avg_pT<avgPtBin[1] && alphaCondition3){
	    if( maxEta >= etaBin[ib] &&  maxEta < etaBin[ib+1] && secondMaxEta >= etaBin[ib] &&  secondMaxEta < etaBin[ib+1])
	      {
		assym= alternate1?(maxPt-secondMaxPt)/(maxPt+secondMaxPt):(secondMaxPt-maxPt)/(maxPt+secondMaxPt);
		if(fabs(assym)>0.6) continue;
		hAssymBins[0][ib]->Fill(assym);
		alternate1=!alternate1;
		break;
	      }
	    else if( maxEta>= etaBin[ib] &&  maxEta< etaBin[ib+1])
	      {
		assym=(maxPt-secondMaxPt)/(maxPt+secondMaxPt);
		if(fabs(assym)>0.6) continue;
		hAssymBins[0][ib]->Fill(assym);
	      }
	    else if(secondMaxEta>= etaBin[ib] &&  secondMaxEta< etaBin[ib+1]) {
	      assym=(secondMaxPt-maxPt)/(maxPt+secondMaxPt);
	      if(fabs(assym)>0.6) continue;
	      hAssymBins[0][ib]->Fill(assym);
	    }	    
	  }
	}   
      }
    }    
    
#endif

#if 1
    if(trig80Condition){
      for (int ib = 0; ib < nEta; ib++) {
	
	if( fabs(maxEta)<=1.3 &&fabs(secondMaxEta)>1.3){
	  float delPhi = TMath::Abs(maxPhi - secondMaxPhi);
	  bool phiCondition = (delPhi >2.7);
	  if(!phiCondition) continue;
	  avg_pT=0.5*(maxPt+secondMaxPt);
	  alpha_Val=thirdMaxPt/avg_pT;
	  bool alphaCondition81=(alpha_Val<alphaValue);
	  if((avg_pT>=avgPtBin[1] &&avg_pT<avgPtBin[2]) && alphaCondition81){
	    if(secondMaxEta >= etaBin[ib] && secondMaxEta < etaBin[ib+1]) {
	      assym=(secondMaxPt-maxPt)/(secondMaxPt+maxPt);
	      if(fabs(assym)>0.6) continue;
	      hAssymBins[1][ib]->Fill(assym);
	    }
	  }
	}
	else if(fabs(secondMaxEta)<=1.3 &&fabs(maxEta)>1.3){
	  float delPhi = TMath::Abs(maxPhi - secondMaxPhi);
	  bool phiCondition = (delPhi >2.7);
	  if(!phiCondition) continue;
	  avg_pT=0.5*(maxPt+secondMaxPt);
	  alpha_Val=thirdMaxPt/avg_pT;
	  bool alphaCondition82=(alpha_Val<alphaValue);
	  if((avg_pT>=avgPtBin[1] &&avg_pT<avgPtBin[2]) && alphaCondition82){
	    if( maxEta>= etaBin[ib] &&  maxEta< etaBin[ib+1]) {
	      assym=(maxPt-secondMaxPt)/(maxPt+secondMaxPt);
	      if(fabs(assym)>0.6) continue;
	      hAssymBins[1][ib]->Fill(assym);
	    }
	  }
	}
	else if((fabs(secondMaxEta)<=1.3 &&fabs(maxEta)<=1.3)){
	  float delPhi = TMath::Abs(maxPhi - secondMaxPhi);
	  bool phiCondition = (delPhi >2.7);
	  if(!phiCondition) continue;
	  avg_pT=0.5*(maxPt+secondMaxPt);
	  alpha_Val=thirdMaxPt/avg_pT;
	  bool alphaCondition82=(alpha_Val<alphaValue);
	  if((avg_pT>=avgPtBin[1] &&avg_pT<avgPtBin[2]) && alphaCondition82){
	    if(maxEta >= etaBin[ib] &&  maxEta < etaBin[ib+1] && secondMaxEta >= etaBin[ib] &&  secondMaxEta < etaBin[ib+1]){
	      assym= alternate2?(maxPt-secondMaxPt)/(maxPt+secondMaxPt):(secondMaxPt-maxPt)/(maxPt+secondMaxPt);
	      if(fabs(assym)>0.6) continue;
	      hAssymBins[1][ib]->Fill(assym);
	      alternate2=!alternate2;
	      break;
	    }
	    else if( maxEta>= etaBin[ib] &&  maxEta< etaBin[ib+1]) {
	      assym=(maxPt-secondMaxPt)/(maxPt+secondMaxPt);
	      if(fabs(assym)>0.6) continue;
	      hAssymBins[1][ib]->Fill(assym);
	    }
	    else if(secondMaxEta>= etaBin[ib] &&  secondMaxEta< etaBin[ib+1]) {
	      assym=(secondMaxPt-maxPt)/(maxPt+secondMaxPt);
	      if(fabs(assym)>0.6) continue;
	      hAssymBins[1][ib]->Fill(assym);
	      //cout<<"case2: "<<etaBin[ib]<<'\t'<<etaBin[ib+1]<<'\t'<<secondMaxEta<<'\t'<<secondMaxPt<<'\t'<<maxEta<<'\t'<<maxPt<<'\t'<<avg_pT<<'\t'<<assym<<endl;
	    }
	  }
	}
      }	
    }    
    
#endif

#if 1
    if(trig110Condition ){
      for(int ipt=2;ipt<nptBin; ipt++){
	for (int ib = 0; ib < nEta; ib++) {
	  if( fabs(maxEta)<=1.3 &&fabs(secondMaxEta)>1.3){
	    float delPhi = TMath::Abs(maxPhi - secondMaxPhi);
	    bool phiCondition = (delPhi >2.7);
	    if(!phiCondition) continue;
	    avg_pT=0.5*(maxPt+secondMaxPt);
	    alpha_Val=thirdMaxPt/avg_pT;
	    bool alphaCondition1101=(alpha_Val<alphaValue);
	    if(!alphaCondition1101) continue;
	    if((avg_pT >= avgPtBin[ipt] && avg_pT < avgPtBin[ipt+1])){
	      if(secondMaxEta >= etaBin[ib] && secondMaxEta < etaBin[ib+1]) {
		assym=(secondMaxPt-maxPt)/(secondMaxPt+maxPt);
		if(fabs(assym)>0.6) continue;
		hAssymBins[ipt][ib]->Fill(assym);
		//cout<<secondMaxEta<<'\t'<<secondMaxPt<<'\t'<<maxEta<<'\t'<<maxPt<<'\t'<<avg_pT<<endl;
	      }	  
	    }
	  }
	  else if(fabs(secondMaxEta)<=1.3 &&fabs(maxEta)>1.3){
	    float delPhi = TMath::Abs(maxPhi - secondMaxPhi);
	    bool phiCondition = (delPhi >2.7);
	    if(!phiCondition) continue;
	    avg_pT=0.5*(maxPt+secondMaxPt);
	    alpha_Val=thirdMaxPt/avg_pT;
	    bool alphaCondition1102=(alpha_Val<alphaValue);
	    if(!alphaCondition1102) continue;
	    if(avg_pT>=avgPtBin[ipt] && avg_pT<avgPtBin[ipt+1]){
	      if( maxEta>= etaBin[ib] &&  maxEta< etaBin[ib+1]){
		assym=(maxPt-secondMaxPt)/(maxPt+secondMaxPt);
		if(fabs(assym)>0.6) continue;
		hAssymBins[ipt][ib]->Fill(assym);
		//cout<<secondMaxEta<<'\t'<<secondMaxPt<<'\t'<<maxEta<<'\t'<<maxPt<<'\t'<<avg_pT<<endl;
	      }
	    }
	  }
	  else if((fabs(secondMaxEta)<=1.3 &&fabs(maxEta)<=1.3)){
	    float delPhi = TMath::Abs(maxPhi - secondMaxPhi);
	    bool phiCondition = (delPhi >2.7);
	    if(!phiCondition) continue;
	    avg_pT=0.5*(maxPt+secondMaxPt);
	    alpha_Val=thirdMaxPt/avg_pT;
	    bool alphaCondition1103=(alpha_Val<alphaValue);
	    if(!alphaCondition1103) continue;
	    if(avg_pT>=avgPtBin[ipt] &&avg_pT<avgPtBin[ipt+1]){
	      if(maxEta >= etaBin[ib] &&  maxEta < etaBin[ib+1] && secondMaxEta >= etaBin[ib] &&  secondMaxEta < etaBin[ib+1]){
		assym= alternate3?(maxPt-secondMaxPt)/(maxPt+secondMaxPt):(secondMaxPt-maxPt)/(maxPt+secondMaxPt);
		if(fabs(assym)>0.6) continue;
		hAssymBins[ipt][ib]->Fill(assym);
		//cout<<"case1: "<<i<<'\t'<<etaBin[ib]<<'\t'<<etaBin[ib+1]<<'\t'<<secondMaxEta<<'\t'<<secondMaxPt<<'\t'<<maxEta<<'\t'<<maxPt<<'\t'<<avg_pT<<'\t'<<assym<<endl;
		alternate3=!alternate3;
		break;
	      }

	      else if( maxEta>= etaBin[ib] &&  maxEta< etaBin[ib+1]) {
		assym=(maxPt-secondMaxPt)/(maxPt+secondMaxPt);
		if(fabs(assym)>0.6) continue;
		hAssymBins[ipt][ib]->Fill(assym);
		//cout<<"case2: "<<i<<'\t'<<etaBin[ib]<<'\t'<<etaBin[ib+1]<<'\t'<<secondMaxEta<<'\t'<<secondMaxPt<<'\t'<<maxEta<<'\t'<<maxPt<<'\t'<<avg_pT<<'\t'<<assym<<endl;
	      }
	      else if(secondMaxEta>= etaBin[ib] &&  secondMaxEta< etaBin[ib+1]) {
		assym=(secondMaxPt-maxPt)/(maxPt+secondMaxPt);
		if(fabs(assym)>0.6) continue;
		hAssymBins[ipt][ib]->Fill(assym);
		//cout<<"case3: "<<i<<'\t'<<etaBin[ib]<<'\t'<<etaBin[ib+1]<<'\t'<<secondMaxEta<<'\t'<<secondMaxPt<<'\t'<<maxEta<<'\t'<<maxPt<<'\t'<<avg_pT<<'\t'<<assym<<endl;
	      }
	    }
	  }	    	
	}
      }
    }
      
#endif

  }
  
  TFile *outputFile = new TFile(Form("/afs/cern.ch/user/u/uacharya/eos/cms_analaysis/JEC_analysis/miniAOD_output2/HiForestMiniAOD_JM0_output_%d.root", infilInd), "RECREATE");
  //hrawPt_unCor->Write();
  //hrawPt->Write();
  for (int ipt = 0; ipt < nptBin; ipt++) {
   for (int ib = 0; ib < nEta; ib++) {
      hAssymBins[ipt][ib]->Write();
      //hAssymBins40[ipt][ib]->Write();
      //hAssymBins80[ipt][ib]->Write();
      //hAssymBins110[ipt][ib]->Write();
    }
  }
  // h_alpha->Write();
  outputFile->Close();
}
