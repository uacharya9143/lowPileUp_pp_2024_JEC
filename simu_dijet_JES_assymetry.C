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

void simu_dijet_imbalance_asymmetry_L2Residual(int infilInd=1) {

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
  //JEC_input_file.push_back("/afs/cern.ch/user/u/uacharya/MC_Truth_JEC/CMSSW_14_1_0/src/MC-truth-JEC/condor/Files/ParallelMCL1_L2Relative_AK4PF.txt");
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
  vector<float>   refpt;
  vector<float>   refeta;
  vector<float>   refphi;
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
  // // tree->SetBranchAddress("refpt",refpt);
  // // tree->SetBranchAddress("refeta",refeta);
  // // tree->SetBranchAddress("refphi",refphi);
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

 
  //int nEta = 60;
  //float etaBin[] = {-5.0, -3.5, -3.0, -2.7, -2.4, -2.0, -1.7, -1.3, -1.0, -0.7, -0.4, 0.0, 0.4, 0.7, 1.0, 1.3, 1.7, 2.0, 2.4, 2.7, 3.0, 3.5, 5.0};
  //float etaBin[] = {-5.191, -4.515, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.051, -1.93, -1.791, -1.653, -1.566, -1.479, -1.392, -1.305, -1.174, -1.044, -0.913, -0.783, -0.652, -0.522, -0.391, -0.261, -0.130, 0.0, 0.130, 0.261, 0.391, 0.522, 0.652, 0.783, 0.913, 1.044, 1.174, 1.305, 1.392, 1.479, 1.566, 1.653, 1.791, 1.93, 2.051, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.515, 5.191};


   const int nEta = 50;
  //float etaBin[] = {-5.0, -3.5, -3.0, -2.7, -2.4, -2.0, -1.7, -1.3, -1.0, -0.7, -0.4, 0.0, 0.4, 0.7, 1.0, 1.3, 1.7, 2.0, 2.4, 2.7, 3.0, 3.5, 5.0};
  float etaBin[] = {-5.191,  -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.051, -1.93, -1.791, -1.653, -1.566, -1.479, -1.392, -1.305, -1.174, -1.044, -0.913, -0.783, -0.652, -0.522, -0.391, -0.261, -0.130, 0.0, 0.130, 0.261, 0.391, 0.522, 0.652, 0.783, 0.913, 1.044, 1.174, 1.305, 1.392, 1.479, 1.566, 1.653, 1.791, 1.93, 2.051, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 5.191};
 
 const int nptBin = 5;
  double avgPtBin[]={50.0,90.0,120.0,170.0,250.0,1000.0};
  
  TH1F *hPt_eta = new TH1F("hPt_eta", "Leading Jet Pt by Eta Bin", nEta, etaBin);
  TH1F *hPt_eta13 = new TH1F("hPt_eta13", "Second Jet Pt by Eta Bin", nEta, etaBin);
  TH1F *hPt_eta_third = new TH1F("hPt_eta_third", "Third Jet Pt by Eta Bin", nEta, etaBin);
  TH1F *hPt_eta_ratio;
  Long64_t nentries = tree->GetEntries();
  Float_t jtPf_CHF,jtPf_MUF,jtPf_NHF,jtPf_NEF, jtPf_CEF;
  Int_t jtPf_CHM;

  TH1F* hAssymBins[nptBin][nEta];
 TH1F *hrawPt_unCor = new TH1F("hrawPt_unCor","hrawPt_unCor;Raw_Leading_P_{T};",1000,0,1000);
  TH1F *hrawPt = new TH1F("hrawPt","hrawPt;Corrected_Leading_P_{T};",1000,0,1000);
 
  for (int ipt = 0; ipt <nptBin; ipt++) {
    for (int i = 0; i < nEta; i++) {
      hAssymBins[ipt][i] = new TH1F(Form("hAssym_etaBin_%d_%d", ipt,i), Form("Asymmetry in Eta Bin %.3f to %.3f", etaBin[i], etaBin[i+1]), 200, -1, 1);
      hAssymBins[ipt][i] ->Sumw2();
    }
  }
  // Event Loop
  bool alternate=true;
  for (Long64_t i = 0; i < nentries; i++) {
  //for (Long64_t i = 0; i < 1000; i++) {
  //for (Long64_t i = 0; i < 500000; i++) {
  //for (Long64_t i = 500000; i < nentries; i++) {
    tree->GetEntry(i);

    float maxPt = -9999, maxEta = -9999, maxPhi = -9999;
    float secondMaxPt = -9999, secondMaxEta = -9999, secondMaxPhi = -9999;
    float thirdMaxPt = -9999, thirdMaxEta = -9999, thirdMaxPhi = -9999;
    float assym;
    float avg_pT;
    float alpha_Val;
     float pt_cor[nref];
    float eta_cor[nref];
    float phi_cor[nref];
  
    for (int j = 0; j < nref; j++) {
      
      jtPf_CHF=(*jtchf)[j];
      jtPf_MUF=(*jtmuf)[j];
      jtPf_NHF=(*jtnhf)[j];
      //jtPf_CHM=jtPfCHM[j];
      jtPf_NEF=(*jtnef)[j];
      jtPf_CEF=(*jtcef)[j];
      jtPf_CHM=(*jtchMult)[j];
      
      // jtPf_CHF=jtchf[j];
      // jtPf_MUF=jtmuf[j];
      // jtPf_NHF=jtnhf[j];
      // //jtPf_CHM=jtPfCHM[j];
      // jtPf_NEF=jtnef[j];
 
      //if( jtPf_MUF<0.001 && jtPf_CHF > 0.01 && jtPf_CHF < 0.99   && jtPf_NHF < 0.99 && jtPf_NEF <0.9) {  
      if( jtPf_MUF<0.80 && jtPf_CHF > 0.01  && jtPf_NHF < 0.99 && jtPf_NEF <0.9 && jtPf_CEF<0.80 && jtPf_CHM>0) {  

	float pt_raw = (*jtpt)[j];
	float eta_raw = (*jteta)[j];
	float phi_raw = (*jtphi)[j];

	JECCor.SetJetPT(pt_raw);
	JECCor.SetJetEta(eta_raw);
	JECCor.SetJetPhi(phi_raw);
	//cout<<i<<'\t'<<j<<'\t'<<"raw: " << pt <<'\t'<< eta << '\t'<< phi <<endl;
	float pt = JECCor.GetCorrectedPT();
	float eta =eta_raw;
	float phi = phi_raw;
	hrawPt_unCor->Fill(pt_raw);
	hrawPt->Fill(pt);


	if (pt > maxPt) {
	  // Update second and third jets before changing maxPt
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
	  //cout<<i<<'\t'<<j<<'\t'<<"third2: " << thirdMaxPt <<'\t'<< thirdMaxEta << '\t'<< thirdMaxPhi <<endl;

	  secondMaxPt = pt;
	  secondMaxEta = eta;
	  secondMaxPhi = phi;
	  //cout<<i<<'\t'<<j<<'\t'<<"second2: " << secondMaxPt <<'\t'<< secondMaxEta << '\t'<< secondMaxPhi <<endl;

	} else if (pt > thirdMaxPt) {
	  // Only update third jet
	  thirdMaxPt = pt;
	  thirdMaxEta = eta;
	  thirdMaxPhi = phi;
	  //cout<<i<<'\t'<<j<<'\t'<<"third3: " << thirdMaxPt <<'\t'<< thirdMaxEta << '\t'<< thirdMaxPhi <<endl;

	}

      }
    }
  

#if 1
    if(1){
      float alphaValue=0.35;
      
      for(int ipt=0;ipt<nptBin; ipt++){
	for (int ib = 0; ib < nEta; ib++) {
	  if( fabs(maxEta)<1.3 &&fabs(secondMaxEta)>1.3){
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
		//cout<<"case40: "<<i<<'\t'<<etaBin[ib]<<'\t'<<etaBin[ib+1]<<'\t'<<secondMaxEta<<'\t'<<secondMaxPt<<'\t'<<maxEta<<'\t'<<maxPt<<'\t'<<avg_pT<<'\t'<<assym<<endl;
	      }	  
	    }
	  }
	  else if(fabs(secondMaxEta)<1.3 &&fabs(maxEta)>1.3){
	    float delPhi = TMath::Abs(maxPhi - secondMaxPhi);
	    bool phiCondition = (delPhi >2.7);
	    if(!phiCondition) continue;
	    avg_pT=0.5*(maxPt+secondMaxPt);
	    alpha_Val=thirdMaxPt/avg_pT;
	    bool alphaCondition1102=(alpha_Val<alphaValue);
	    if(!alphaCondition1102) continue;
	    //if((avg_pT>=avgPtBin[0] &&avg_pT<avgPtBin[1])){
	    if(avg_pT>=avgPtBin[ipt] && avg_pT<avgPtBin[ipt+1]){
	      if( maxEta>= etaBin[ib] &&  maxEta< etaBin[ib+1]){
		assym=(maxPt-secondMaxPt)/(maxPt+secondMaxPt);
		if(fabs(assym)>0.6) continue;	      
		hAssymBins[ipt][ib]->Fill(assym);
		//cout<<"case80: "<<i<<'\t'<<etaBin[ib]<<'\t'<<etaBin[ib+1]<<'\t'<<secondMaxEta<<'\t'<<secondMaxPt<<'\t'<<maxEta<<'\t'<<maxPt<<'\t'<<avg_pT<<'\t'<<assym<<endl;
	      }
	    }
	  }
	  else if((fabs(secondMaxEta)<1.3 &&fabs(maxEta)<1.3)){
	    float delPhi = TMath::Abs(maxPhi - secondMaxPhi);
	    bool phiCondition = (delPhi >2.7);
	    if(!phiCondition) continue;
	    avg_pT=0.5*(maxPt+secondMaxPt);
	    alpha_Val=thirdMaxPt/avg_pT;
	    bool alphaCondition1103=(alpha_Val<alphaValue);
	    if(!alphaCondition1103) continue;
	    if(avg_pT>=avgPtBin[ipt] &&avg_pT<avgPtBin[ipt+1]){
	      if(maxEta >= etaBin[ib] &&  maxEta < etaBin[ib+1] && secondMaxEta >= etaBin[ib] &&  secondMaxEta < etaBin[ib+1]){
		assym= alternate ? (maxPt-secondMaxPt)/(maxPt+secondMaxPt) : (secondMaxPt-maxPt)/(maxPt+secondMaxPt);
		//assym=(maxPt-secondMaxPt)/(maxPt+secondMaxPt);
		if(fabs(assym)>0.6) continue;	      
		hAssymBins[ipt][ib]->Fill(assym);
		//filled = true;
		//cout<<"case1: "<<i<<'\t'<<etaBin[ib]<<'\t'<<etaBin[ib+1]<<'\t'<<secondMaxEta<<'\t'<<secondMaxPt<<'\t'<<maxEta<<'\t'<<maxPt<<'\t'<<avg_pT<<'\t'<<assym<<endl;
		alternate = !alternate;
		break;
	      }
	      else if( maxEta>= etaBin[ib] &&  maxEta< etaBin[ib+1]) {
		assym=(maxPt-secondMaxPt)/(maxPt+secondMaxPt);
		if(fabs(assym)>0.6) continue;	      
		hAssymBins[ipt][ib]->Fill(assym);
		//cout<<"case2: "<<i<<'\t'<<etaBin[ib]<<'\t'<<etaBin[ib+1]<<'\t'<<secondMaxEta<<'\t'<<secondMaxPt<<'\t'<<maxEta<<'\t'<<maxPt<<'\t'<<avg_pT<<'\t'<<assym<<endl;
		//break;
	      }
	      else if(secondMaxEta>= etaBin[ib] &&  secondMaxEta< etaBin[ib+1]) {
		assym=(secondMaxPt-maxPt)/(maxPt+secondMaxPt);
		if(fabs(assym)>0.6) continue;	      
		hAssymBins[ipt][ib]->Fill(assym);
		//cout<<"case3: "<<i<<'\t'<<etaBin[ib]<<'\t'<<etaBin[ib+1]<<'\t'<<secondMaxEta<<'\t'<<secondMaxPt<<'\t'<<maxEta<<'\t'<<maxPt<<'\t'<<avg_pT<<'\t'<<assym<<endl;
		//break;
	      }
	    }
	  }	    	
	}
      }
    }
      
#endif
  }
 
#if 1

  // Save Output
  TFile *outputFile = new TFile(Form("/afs/cern.ch/user/u/uacharya/eos/cms_analaysis/JEC_analysis/simuMC_output/JRA_output_%d.root",infilInd), "RECREATE");
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
#endif

  
}
