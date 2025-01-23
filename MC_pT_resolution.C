#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TKey.h>
#include <iostream>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <fstream> // Required for file operations
Double_t fgaus1(Double_t *x, Double_t *par){

  Double_t xx = x[0];

  Double_t fitval = 
    par[0]*par[1]*(1.0/(sqrt(2.*3.14159)*par[3]))*TMath::Gaus(xx,par[2],par[3])+par[4];
  // par[0]*par[1]*TMath::Gaus(xx,par[2],par[3]);
  return fitval;
}
Double_t fitfunctjec(Double_t *x, Double_t *par){

  Double_t xx = x[0];

  Double_t fitval = sqrt(pow(par[0],2)+pow(par[1],2)/xx+pow((par[2]/xx),2)+pow((par[3]/xx),3));

  return fitval;
}
void MC_pT_resolution(){
#if 1
  ofstream outFile_MC("jer_sigma_rootfiles_skiplowestPtBin/MC_pTDependent_resolution_lowPU_pp13TeV_with_newStatistics.txt");
  if (!outFile_MC.is_open()) {
    cerr << "Error: Could not open file for writing." << endl;
    return;
  }
  outFile_MC << "{1 JetEta 1 JetPt sqrt(pow(p0,2)+pow(p1,2)/x+pow((p2/x),2)+pow((p3/x),3)) Resolution}\n";
#endif

gROOT->Reset();
  // const int npt= 5;
  gStyle->SetOptStat(0);
  
  float alphaBin[]={0.0,0.15,0.20,0.25,0.30,0.35};
  const int nEta = 50;
  float etaBin[] = {-5.191,  -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.051, -1.93, -1.791, -1.653, -1.566, -1.479, -1.392, -1.305, -1.174, -1.044, -0.913, -0.783, -0.652, -0.522, -0.391, -0.261, -0.130, 0.0, 0.130, 0.261, 0.391, 0.522, 0.652, 0.783, 0.913, 1.044, 1.174, 1.305, 1.392, 1.479, 1.566, 1.653, 1.791, 1.93, 2.051, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 5.191};
 const int npt = 5;
 double avgPtBin[]={50.0,90.0,120.0,170.0,250.0,1000.0};
 TF1 *ftsp1[nEta][npt];
 TF1 *ftsp2[nEta];

  

 TFile *fIn= new TFile("pTdependent/pT_reco_pt_Gen_Ratio_JER.root");
 TH1F *hpt_gen_reco[nEta][npt];
 TH1F *hpt_gen_reco_fit[nEta][npt];
 TH1F * hSigmaPlot[nEta];
 //TCanvas *cA[nEta];
 TCanvas *cAf[nEta];
 float sigmaVal[nEta][npt];
 float dsigmaVal[nEta][npt];
 TCanvas *cAs=new TCanvas("cAs","cAs",10,10,1200,1200);
 cAs->Divide(5,10);
 for(int ita=0;ita<nEta;ita++){
   //cA[ita]= new TCanvas(Form("cA%d",ita),Form("cA%d",ita),10,10,1200,800);
   //cA[ita]->Divide(3,2);
   cAf[ita]= new TCanvas(Form("cAf%d",ita),Form("cAf%d",ita),10,10,1200,800);
   cAf[ita]->Divide(3,2);
   hSigmaPlot[ita]=new TH1F(Form("hSigmaPlot_%d", ita), Form("hSigma_plot_ptReco_ptGen_%.1f_to_%.1f", etaBin[ita], etaBin[ita+1]), npt,avgPtBin );
   
   for(int ipt=0;ipt<npt; ipt++){
     // cA[ita]->cd(ipt+1);
     hpt_gen_reco[ita][ipt]=(TH1F*)fIn->Get(Form("hGenpT_RecoPt_etaBin_%d_%d",ita,ipt));
     hpt_gen_reco[ita][ipt]->SetTitle(Form("GEN_RECO_PT_Ratio_%.1f < p_{T} < %.1f;p_{T}^{Reco}/p_{T}^{Gen};",avgPtBin[ipt],avgPtBin[ipt+1]));
     //hpt_gen_reco[ita][ipt]->Draw();
     
     cAf[ita]->cd(ipt+1);
     gPad->SetGridy();
     gPad->SetGridx();
     hpt_gen_reco_fit[ita][ipt]=(TH1F *)hpt_gen_reco[ita][ipt]->Clone();
     hpt_gen_reco_fit[ita][ipt]->SetTitle(Form("GEN_RECO_PT_Ratio_fit%.1f < p_{T} < %.1f;p_{T}^{Reco}/p_{T}^{Gen};",avgPtBin[ipt],avgPtBin[ipt+1]));
     
     ftsp1[ita][ipt]=new TF1(Form("ftsp1%d%d",ita,ipt),fgaus1,0,2,5);
     ftsp1[ita][ipt]->SetLineColor(2);
     ftsp1[ita][ipt]->SetParameter(0,hpt_gen_reco_fit[ita][ipt]->GetBinWidth(1));
     ftsp1[ita][ipt]->SetParameter(1, hpt_gen_reco_fit[ita][ipt]->GetMaximum());
     ftsp1[ita][ipt]->SetParameter(2, hpt_gen_reco_fit[ita][ipt]->GetMean());
     ftsp1[ita][ipt]->SetParLimits(2, 0.5,1.5);  // Mean near zero
     ftsp1[ita][ipt]->SetParameter(3, hpt_gen_reco_fit[ita][ipt]->GetRMS());
     ftsp1[ita][ipt]->SetParameter(4, 0.05);
     hpt_gen_reco_fit[ita][ipt]->Fit(ftsp1[ita][ipt],"Q");
     cAf[ita]->SaveAs(Form("pTdependent/pT_resolution_perEtaBin_fit_%1f_to_%1f.pdf",etaBin[ita],etaBin[ita+1]),"RECREATE");


     
     cout<<ita<<'\t'<<ipt<<'\t'<<ftsp1[ita][ipt]->GetParameter(3)<<'\t'<<hpt_gen_reco_fit[ita][ipt]->GetStdDev()<<endl;
     // hpt_gen_reco_fit[ita][ipt]->Fit("gaus","Q");
     sigmaVal[ita][ipt]=ftsp1[ita][ipt]->GetParameter(3);
     dsigmaVal[ita][ipt]=ftsp1[ita][ipt]->GetParError(3);
     hSigmaPlot[ita]->SetBinContent(ipt+1,sigmaVal[ita][ipt]);
     hSigmaPlot[ita]->SetBinError(ipt+1,dsigmaVal[ita][ipt]);
     hSigmaPlot[ita]->SetMarkerStyle(20);
     hSigmaPlot[ita]->SetMarkerColor(4);
     hSigmaPlot[ita]->SetLineColor(4);
     hSigmaPlot[ita]->SetTitle(Form(";p_{T}^{Gen};#sigma_{p_{T}}"));
     hSigmaPlot[ita]->GetYaxis()->SetRangeUser(0.03,0.2);
     if(ita==1||ita==2||ita==3||ita==4||ita==46||ita==47||ita==48)
       {
	 hSigmaPlot[ita]->GetYaxis()->SetRangeUser(0.03,0.3);

       }
   }
   cAs->cd(ita+1);
   hSigmaPlot[ita]->Draw();
   ftsp2[ita]=new TF1(Form("ftsp2%d",ita),fitfunctjec,50,1000,4);
   ftsp2[ita]->SetLineColor(2);
   ftsp2[ita]->SetParameter(0,1.50);
   ftsp2[ita]->SetParameter(1,0.50);
   ftsp2[ita]->SetParameter(2,-2.50);
   ftsp2[ita]->SetParameter(3,1.50);
   hSigmaPlot[ita]->Fit(ftsp2[ita],"Q");
   outFile_MC << etaBin[ita] <<'\t'
	     << etaBin[ita+1] <<'\t'
	     << 6 <<'\t'
	     << 50 <<'\t'
	     << 1000 <<'\t'
	     << ftsp2[ita]->GetParameter(0) <<'\t'
	     << ftsp2[ita]->GetParameter(1)<<'\t'
	     << ftsp2[ita]->GetParameter(2)<<'\t'
	     << ftsp2[ita]->GetParameter(3) <<'\n';

 }
   cAs->SaveAs(Form("pTdependent/pT_resolution_perEtaBin_vs_pT.pdf"),"RECREATE");
 outFile_MC.close();
   TCanvas *cnn= new TCanvas();
   hSigmaPlot[48]->Draw();





}
