#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TKey.h>
#include <iostream>
#include "PlotHelper4.h"

Double_t fgaus1(Double_t *x, Double_t *par){

  Double_t xx = x[0];

  Double_t fitval = 
    par[0]*par[1]*(1.0/(sqrt(2.*3.14159)*par[3]))*TMath::Gaus(xx,par[2],par[3]);//+par[4];
  // par[0]*par[1]*TMath::Gaus(xx,par[2],par[3]);
  return fitval;
}
Double_t fgaus2(Double_t *x, Double_t *par){

  Double_t xx = x[0];

  Double_t fitval = 
    par[0]*par[1]*(1.0/(sqrt(2.*3.14159)*par[3]))*TMath::Gaus(xx,par[2],par[3]);//+par[4];
  //par[0]*par[1]*TMath::Gaus(xx,par[2],par[3]);
  return fitval;
}

void dijet_resolution_assymetry_and_alphaVal(const char *alphaCond="015", double alphaCd=0.15){


  gROOT->Reset();

  //####################################################################################################
  //###################### WIth AsymCut of |Assym|<0.55########################################
  //####################################################################################################
#if 1
  TFile *fIn1= new TFile(Form("data_res_files/data_JER_StdMthd_%s.root",alphaCond));

  TFile *fIn2= new TFile(Form("simu_res_files/simu_JER_StdMthd_%s.root",alphaCond));
#endif

  //####################################################################################################

  
  /* const int nEta = 22; */
  /* float etaBin[] = {-5.0, -3.5, -3.0, -2.7, -2.4, -2.0, -1.7, -1.3, -1.0, -0.7, -0.4, 0.0, 0.4, 0.7, 1.0, 1.3, 1.7, 2.0, 2.4, 2.7, 3.0, 3.5, 5.0}; */

  /* const int npt = 5; */
  /* double avgPtBin[]={50.0,90.0,120.0,160.0,200.0,500.0}; */
  
  const int nEta = 50;
  float etaBin[] = {-5.191,  -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.051, -1.93, -1.791, -1.653, -1.566, -1.479, -1.392, -1.305, -1.174, -1.044, -0.913, -0.783, -0.652, -0.522, -0.391, -0.261, -0.130, 0.0, 0.130, 0.261, 0.391, 0.522, 0.652, 0.783, 0.913, 1.044, 1.174, 1.305, 1.392, 1.479, 1.566, 1.653, 1.791, 1.93, 2.051, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 5.191};

  const int npt = 5;
  double avgPtBin[]={50.0,90.0,120.0,170.0,250.0,1000.0};

  //TH1F *hAssym[30];
  TH1F *hAssym1[npt][50];
  TH1F *hAss_nent1[npt];
  float meanAssym1[npt][50];
  float dmeanAssym1[npt][50];
  float response1[npt][50];
  float dresponse1[npt][50];
  TH1F *hsigma1[npt];
  float nEntries_data[npt][50];
  TH1F *hAssym2[npt][50];
  float meanAssym2[npt][50];
  float dmeanAssym2[npt][50];

  float stdDev1[npt][50];
  float dstdDev1[npt][50];
  float resolution1[npt][50];
  
  float response2[npt][50];
  float dresponse2[npt][50];
  TH1F *hsigma2[npt];
  
  float stdDev2[npt][50];
  float dstdDev2[npt][50];
  float resolution2[npt][50];
  
  float nEntries_sim[npt][50];
  TH1F *hAss_nent2[npt];

  TH1F * hDataSimDiff[npt];
  TH1F * hDataSimDiff_nEntries[npt];
  //TH2F *hResponse=new TH2F("hResponse","hResponse",50,-5,5,20,0,2);
  TCanvas *cA[npt];//= new TCanvas("cA","cA",10,10,750,750);
  TCanvas *cAs[npt];//= new TCanvas("cA","cA",10,10,750,750);

  //PdfFileHelper PdfFile1("Assymetry_Plots.pdf");
  //PdfFileHelper PdfFile2("Response.pdf");
  //PdfFileHelper PdfFile3("Response_Ratio_Plots.pdf");

  TF1 *ftsp1[npt][50];
  TF1 *ftsp2[npt][50];

  
  for(int ipt=0;ipt<npt; ipt++){
    
    cA[ipt]= new TCanvas(Form("cA%d",ipt),Form("cA%d",ipt),15,10,800,800);
    cA[ipt]->Divide(5,10);
    cAs[ipt]= new TCanvas(Form("cAs%d",ipt),Form("cAs%d",ipt),15,10,800,800);
    cAs[ipt]->Divide(5,10);
    hsigma1[ipt]=new TH1F(Form("hsigma_data_%d",ipt),Form("hsigma_data_%d",ipt),nEta,etaBin);
    hsigma1[ipt]->Sumw2();
    hsigma1[ipt] ->SetTitle(Form("Data @ #alpha < %.2f;#eta;#sigma_{A}",alphaCd));

    hsigma2[ipt]=new TH1F(Form("hsigma_simu_%d",ipt),Form("hsigma_simu_%d",ipt),nEta,etaBin);
    hsigma2[ipt]->Sumw2();
    hsigma2[ipt] ->SetTitle(Form("MC_Simulation  @ #alpha < %.2f;#eta;#sigma_{A}",alphaCd));


    hAss_nent1[ipt]=new TH1F(Form("hAss_nent1%d",ipt),Form("hAss_nent1%d",ipt),nEta,etaBin);
    hAss_nent1[ipt]->Sumw2();
    hAss_nent1[ipt] ->SetTitle(Form("Data @ #alpha < %.2f;#eta;No_of_events",alphaCd));

    hAss_nent2[ipt]=new TH1F(Form("hAss_nent2%d",ipt),Form("hAss_nent2%d",ipt),nEta,etaBin);
    hAss_nent2[ipt]->Sumw2();
    hAss_nent2[ipt] ->SetTitle(Form("MC_Simulation  @ #alpha < %.2f;#eta;No_of_events",alphaCd));

    for(int ib=0;ib<nEta;ib++){
      //hAssym[ipt][ib]=(TH1F *) fIn->Get(Form("hAssym_etaBin110_%d_%d",ipt,ib));
       if(ib==0 ||ib==1||ib==2||ib==3||ib==46||ib==47||ib==48||ib==49) continue;
      gStyle->SetOptStat(0);
      hAssym1[ipt][ib]=(TH1F *) fIn1->Get(Form("hAssym_etaBin_%d_%d",ipt,ib));
      nEntries_data[ipt][ib]=hAssym1[ipt][ib]->GetEntries();
      hAssym1[ipt][ib]->Rebin(5);

      cA[ipt]->cd(ib+1);
      cA[ipt]->SaveAs(Form("jerHistos/Data_Asymmetry_@_\u03B1_%s_%d.pdf",alphaCond,ipt),"RECREATE");
      //cA[ipt]->SaveAs(Form("jerHistos_skiplowestPtBin/Data_Asymmetry_@_\u03B1_%s_%d.pdf",alphaCond,ipt),"RECREATE");
      //cA[ipt]->SaveAs(Form("jer_histos2/Data_Asymmetry_@_\u03B1_%s_%d_noAssymCut.pdf",alphaCond,ipt),"RECREATE");
      gPad->SetGridy();
      gPad->SetGridx();
 
      ftsp1[ipt][ib]=new TF1(Form("ftsp1%d%d",ipt,ib),fgaus1,-0.6,0.6,4);
      ftsp1[ipt][ib]->SetLineColor(2);
       ftsp1[ipt][ib]->SetParameter(0,hAssym1[ipt][ib]->GetBinWidth(1));
      // ftsp1[ipt][ib]->FixParameter(0,hAssym1[ipt][ib]->GetBinWidth(1));
      // ftsp1[ipt][ib]->SetParameter(1,1000);
      //ftsp1[ipt][ib]->SetParameter(2,0.01);
      //ftsp1[ipt][ib]->SetParameter(3,0.15);
     ftsp1[ipt][ib]->SetParameter(1, hAssym1[ipt][ib]->GetMaximum());
      ftsp1[ipt][ib]->SetParameter(2, hAssym1[ipt][ib]->GetMean());
      ftsp1[ipt][ib]->SetParLimits(2, -0.1, 0.1);  // Mean near zero
      ftsp1[ipt][ib]->SetParameter(3, hAssym1[ipt][ib]->GetRMS());
      //ftsp1[ipt][ib]->SetParLimits(3, 0.0, 0.5);  // Reasonable sigma range
      //ftsp1[ipt][ib]->SetParameter(4, -0.5);
     hAssym1[ipt][ib]->Fit(ftsp1[ipt][ib],"QRL");
      hAssym1[ipt][ib]->Draw();
    int fitStatus1 = hAssym1[ipt][ib]->Fit(ftsp1[ipt][ib], "QS");
      if (fitStatus1 != 0) {
	cout << "Data_Fit failed for ipt: " << ipt << ", ib: " << ib << endl;
}
     
      // meanAssym1[ipt][ib]=hAssym1[ipt][ib]->GetMean();     
      // dmeanAssym1[ipt][ib]=hAssym1[ipt][ib]->GetMeanError();

      meanAssym1[ipt][ib]=ftsp1[ipt][ib]->GetParameter(2);     
      dmeanAssym1[ipt][ib]=ftsp1[ipt][ib]->GetParError(2);

      if(1)
	{
	  //stdDev1[ipt][ib]=sqrt(2.0)*hAssym1[ipt][ib]->GetStdDev();
	  //dstdDev1[ipt][ib]=stdDev1[ipt][ib]/sqrt(2 * hAssym1[ipt][ib]->GetEntries());
	  stdDev1[ipt][ib]=sqrt(2.0)*ftsp1[ipt][ib]->GetParameter(3);
	  dstdDev1[ipt][ib]=sqrt(2.0)*ftsp1[ipt][ib]->GetParError(3);
	  cout<<ipt<<" "<<ib<<" "<<" data: " <<"fit: "<< sqrt(2.0)*ftsp1[ipt][ib]->GetParameter(3)<< " "<<"hist: " << sqrt(2.0)*hAssym1[ipt][ib]->GetStdDev()<<endl;
	  
	}
     
      
      hsigma1[ipt]->SetBinContent(ib+1,stdDev1[ipt][ib]);
      hsigma1[ipt]->SetBinError(ib+1,dstdDev1[ipt][ib]);
      hsigma1[ipt]->GetYaxis()->SetRangeUser(0.0,0.5);

      hsigma1[ipt]->SetMarkerStyle(8);
      hsigma1[ipt]->SetMarkerColor(ipt+1);
      hsigma1[ipt]->SetLineColor(ipt+1);
      if(ipt==4)
	{
	  hsigma1[ipt]->SetMarkerColor(6);
	  hsigma1[ipt]->SetLineColor(6);
 
	}
      cAs[ipt]->cd(ib+1);
       cAs[ipt]->SaveAs(Form("jerHistos/MC_Asymmetry_@_\u03B1_%s_%d.pdf",alphaCond,ipt),"RECREATE");
      //cAs[ipt]->SaveAs(Form("jerHistos_skiplowestPtBin/MC_Asymmetry_@_\u03B1_%s_%d.pdf",alphaCond,ipt),"RECREATE");
       
      gPad->SetGridy();
      gPad->SetGridx();
    
      hAssym2[ipt][ib]=(TH1F *) fIn2->Get(Form("hAssym_etaBin_%d_%d",ipt,ib));
      hAssym2[ipt][ib]->Rebin(5);

      ftsp2[ipt][ib]=new TF1(Form("ftsp2%d%d",ipt,ib),fgaus2,-0.6,0.6,4);
      ftsp2[ipt][ib]->SetLineColor(2);
      ftsp2[ipt][ib]->SetParameter(0,hAssym2[ipt][ib]->GetBinWidth(1));
      ftsp2[ipt][ib]->FixParameter(0,hAssym2[ipt][ib]->GetBinWidth(1));
      /* ftsp2[ipt][ib]->SetParameter(1,1000); */
      /* ftsp2[ipt][ib]->SetParameter(2,0.01); */
      /* ftsp2[ipt][ib]->SetParameter(3,0.15); */
      ftsp2[ipt][ib]->SetParameter(1, hAssym2[ipt][ib]->GetMaximum());
      ftsp2[ipt][ib]->SetParameter(2, hAssym2[ipt][ib]->GetMean());
      ftsp1[ipt][ib]->SetParLimits(2, -0.1, 0.1);  // Mean near zero
      ftsp2[ipt][ib]->SetParameter(3, hAssym2[ipt][ib]->GetRMS());
      //ftsp2[ipt][ib]->SetParLimits(3, 0.0, 0.35);  // Reasonable sigma range
      hAssym2[ipt][ib]->Fit(ftsp2[ipt][ib],"QRL");
      //ftsp2[ipt][ib]->SetParameter(4, -0.5);
    hAssym2[ipt][ib]->Draw();
      int fitStatus2 = hAssym2[ipt][ib]->Fit(ftsp2[ipt][ib], "QS");
      if (fitStatus2 != 0) {
	cout << " Simu_Fit failed for ipt: " << ipt << ", ib: " << ib << endl;
}

      if(1)
   	{
	  //stdDev2[ipt][ib]=sqrt(2.0)*hAssym2[ipt][ib]->GetStdDev();
	  //dstdDev2[ipt][ib]=stdDev2[ipt][ib]/sqrt(2 * hAssym2[ipt][ib]->GetEntries());
	  stdDev2[ipt][ib]=sqrt(2.0)*ftsp2[ipt][ib]->GetParameter(3);
	  dstdDev2[ipt][ib]=sqrt(2.0)*ftsp2[ipt][ib]->GetParError(3);
	  cout<<ipt<<" "<<ib<<" "<<" sim: " <<"fit: "<< sqrt(2.0)*ftsp2[ipt][ib]->GetParameter(3)<< " "<<"hist: " << sqrt(2.0)*hAssym2[ipt][ib]->GetStdDev()<<endl;

	}
      

     
      hsigma2[ipt]->SetBinContent(ib+1,stdDev2[ipt][ib]);
      hsigma2[ipt]->SetBinError(ib+1,dstdDev2[ipt][ib]);
      hsigma2[ipt]->GetYaxis()->SetRangeUser(0.0,0.5);

      hsigma2[ipt]->SetMarkerStyle(8);
      hsigma2[ipt]->SetMarkerColor(ipt+1);
      hsigma2[ipt]->SetLineColor(ipt+1);
      if(ipt==4)
	{
	  hsigma2[ipt]->SetMarkerColor(6);
	  hsigma2[ipt]->SetLineColor(6);
 
	}
    }
  }
  TLegend *leg= new TLegend(0.4,0.65,0.6,0.9);
  leg->AddEntry(hsigma2[0]," 50 < p_{T}_{avg} < 90 ","pl");
  leg->AddEntry(hsigma2[1]," 90 < p_{T}_{avg} < 120 ","pl");
  leg->AddEntry(hsigma2[2]," 120 < p_{T}_{avg} < 170 ","pl");
  leg->AddEntry(hsigma2[3]," 170 < p_{T}_{avg} < 250 ","pl");
  leg->AddEntry(hsigma2[4]," 250 < p_{T}_{avg} < 1000 ","pl");
  
  TCanvas *cR2= new TCanvas("cR2","cR2",15,10,800,800);
  cR2->Divide(1,2);
  cR2->cd(1);
  gPad->SetGridy();
  gPad->SetGridx();
  hsigma1[0]->Draw();
  hsigma1[1]->Draw();//"same");
  hsigma1[2]->Draw("same");
  hsigma1[3]->Draw("same");
  hsigma1[4]->Draw("same");
  //leg->Draw();

  cR2->cd(2);
  gPad->SetGridy();
  gPad->SetGridx();
  hsigma2[0]->Draw();
  hsigma2[1]->Draw();//"same");
  hsigma2[2]->Draw("same");
  hsigma2[3]->Draw("same");
  hsigma2[4]->Draw("same");
  leg->Draw();

  cR2->SaveAs(Form("jerHistos/JER_sigma_hist_with_fit_@_\u03B1_%s.pdf",alphaCond),"RECREATE");
  //cR2->SaveAs(Form("jerHistos_skiplowestPtBin/JER_sigma_hist_with_fit_@_\u03B1_%s.pdf",alphaCond),"RECREATE");
  //cR2->SaveAs(Form("jerHistos/JER_sigma_hist_without_fit_@_\u03B1_%s.pdf",alphaCond),"RECREATE");
 
#if 1
   TFile *outputFile = new TFile(Form("jer_sigma_rootfiles/dijet_jer_standard_sigma_plot_alphaCond_newStat_noReBin_fit%s.root", alphaCond), "RECREATE");
  //TFile *outputFile = new TFile(Form("jer_sigma_rootfiles/dijet_jer_standard_sigma_plot_alphaCond_newStat_nofit%s.root", alphaCond), "RECREATE");
			 
  for (int ipt = 0; ipt < npt; ipt++) {
    //if(ipt==0)continue;
    hsigma1[ipt]->Write();
    hsigma2[ipt]->Write();
  }
  outputFile->Close();
#endif

}
