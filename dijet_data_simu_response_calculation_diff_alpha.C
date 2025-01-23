#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TKey.h>
#include <iostream>
#include "PlotHelper4.h"

void dijet_data_simu_response_calculation_diff_alpha(const char *alphaCond="015", double alphaCd=0.15){


  gROOT->Reset();


 
  /* /\* TFile *fIn1= new TFile(Form("data/dijet_trigger_response_vs_eta_alpha%s_2.root",alphaCond)); *\/ */

  /* /\* TFile *fIn2= new TFile(Form("simu/dijet_balance_after_JEC_simu_alphaCond%s_final2.root",alphaCond)); *\/ */

  /* TFile *fIn1= new TFile(Form("data_with_new_statistics/data_assym_l2Residuals_alphaCond_%s.root",alphaCond)); */
  /* TFile *fIn2= new TFile(Form("simu_with_new_statistics/simu_assym_l2Residuals_alphaCond_%s.root",alphaCond)); */


  TFile *fIn1= new TFile(Form("data_with_new_statistics/data_L2Residual_Asymmerty_%s.root",alphaCond));
  TFile *fIn2= new TFile(Form("simu_with_new_statistics/simu_L2Residual_Asymmerty_%s.root",alphaCond));
  
//TFile *fIn= new TFile("jec_dijet_balance_without_JEC.root");
   //const int nEta = 30;
   //float etaBin[] = {-5.0, -4.7, -4.4, -4.0, -3.7, -3.3, -3.0, -2.7, -2.4, -2.0, -1.7, -1.3, -1.0, -0.7, -0.4, 0.0, 0.4, 0.7, 1.0, 1.3, 1.7, 2.0, 2.4, 2.7, 3.0, 3.4, 3.7, 4.0, 4.4, 4.7, 5.0};

  /*   const int nEta = 22; */
  /* float etaBin[] = {-5.0, -3.5, -3.0, -2.7, -2.4, -2.0, -1.7, -1.3, -1.0, -0.7, -0.4, 0.0, 0.4, 0.7, 1.0, 1.3, 1.7, 2.0, 2.4, 2.7, 3.0, 3.5, 5.0}; */

  /* const int npt = 5; */
  /* double avgPtBin[]={50.0,90.0,120.0,160.0,200.0,500.0}; */


 int nEta = 50;
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
  TH1F *hResponse1[npt];
  float nEntries_data[npt][50];
  TH1F *hAssym2[npt][50];
  float meanAssym2[npt][50];
  float dmeanAssym2[npt][50];

  float dstdDev1[npt][50];
  float response2[npt][50];
  float dresponse2[npt][50];
  TH1F *hResponse2[npt];
  float dstdDev2[npt][50];
  float nEntries_sim[npt][50];
  TH1F *hAss_nent2[npt];

  TH1F * hDataSimDiff[npt];
  TH1F * hDataSimDiff_nEntries[npt];
  //TH2F *hResponse=new TH2F("hResponse","hResponse",50,-5,5,20,0,2);
  TCanvas *cA[npt];//= new TCanvas("cA","cA",10,10,750,750);
  TCanvas *cAs[npt];//= new TCanvas("cA","cA",10,10,750,750);
  //PdfFileHelper PdfFile1("Assymetry_Plots.pdf");
  //PdfFileHelper PdfFile2("Response.pdf");
  // PdfFileHelper PdfFile3("Response_Ratio_Plots.pdf");
  for(int ipt=0;ipt<npt; ipt++){
    
    cA[ipt]= new TCanvas(Form("cA%d",ipt),Form("cA%d",ipt),15,10,800,800);
    cA[ipt]->Divide(10,5);
    cAs[ipt]= new TCanvas(Form("cAs%d",ipt),Form("cAs%d",ipt),15,10,800,800);
    cAs[ipt]->Divide(10,5);
    hResponse1[ipt]=new TH1F(Form("hResponse_dat_%d",ipt),Form("hResponse_dat_%d",ipt),nEta,etaBin);
    hResponse1[ipt]->Sumw2();
    hResponse1[ipt] ->SetTitle(Form("Data @ #alpha < %.2f;#eta;#frac{1+<A>}{1-<A>}",alphaCd));

    hResponse2[ipt]=new TH1F(Form("hResponse_sim_%d",ipt),Form("hResponse_sim_%d",ipt),nEta,etaBin);
    hResponse2[ipt]->Sumw2();
    hResponse2[ipt] ->SetTitle(Form("MC_Simulation  @ #alpha < %.2f;#eta;#frac{1+<A>}{1-<A>}",alphaCd));


    hAss_nent1[ipt]=new TH1F(Form("hAss_nent1%d",ipt),Form("hAss_nent1%d",ipt),nEta,etaBin);
    hAss_nent1[ipt]->Sumw2();
    hAss_nent1[ipt] ->SetTitle(Form("Data @ #alpha < %.2f;#eta;No_of_events",alphaCd));

    hAss_nent2[ipt]=new TH1F(Form("hAss_nent2%d",ipt),Form("hAss_nent2%d",ipt),nEta,etaBin);
    hAss_nent2[ipt]->Sumw2();
    hAss_nent2[ipt] ->SetTitle(Form("MC_Simulation  @ #alpha < %.2f;#eta;No_of_events",alphaCd));

   for(int ib=0;ib<nEta;ib++){
     //hAssym[ipt][ib]=(TH1F *) fIn->Get(Form("hAssym_etaBin110_%d_%d",ipt,ib));
     gStyle->SetOptStat(0);
     hAssym1[ipt][ib]=(TH1F *) fIn1->Get(Form("hAssym_etaBin_%d_%d",ipt,ib));
     nEntries_data[ipt][ib]=hAssym1[ipt][ib]->GetEntries();
     meanAssym1[ipt][ib]=hAssym1[ipt][ib]->GetMean();
     //dmeanAssym1[ipt][ib]=hAssym1[ipt][ib]->GetRMSError();
     dstdDev1[ipt][ib]=hAssym1[ipt][ib]->GetStdDev();
     dmeanAssym1[ipt][ib]=hAssym1[ipt][ib]->GetMeanError();
     // if((1-meanAssym[ipt][ib])>0)
     if(1)
       {response1[ipt][ib]=(1+meanAssym1[ipt][ib])/(1-meanAssym1[ipt][ib]);}
      //if(dmeanAssym1[ipt][ib]>0){
      if((1-meanAssym1[ipt][ib])>0){
	dresponse1[ipt][ib]=(2*dmeanAssym1[ipt][ib])/pow((1-meanAssym1[ipt][ib]),2);
	}
      cA[ipt]->cd(ib+1);
      // cA[ipt]->SaveAs(Form("histPlots_noAsymCut/Data_Asymmetry_@_\u03B1_%s_%d.pdf",alphaCond,ipt),"RECREATE");
      gPad->SetGridy();
      gPad->SetGridx();
      hAssym1[ipt][ib]->Draw();
      hAssym1[ipt][ib]->SetMarkerStyle(7);
      hAssym1[ipt][ib]->SetMarkerColor(4);
      hAssym1[ipt][ib]->SetLineColor(4);
     
      //hAssym[ipt][ib]->GetXaxis()->SetRangeUser(0,1);
      
      hResponse1[ipt]->SetBinContent(ib+1,response1[ipt][ib]);
      hResponse1[ipt]->SetBinError(ib+1,dresponse1[ipt][ib]);

      hAss_nent1[ipt]->SetBinContent(ib+1,nEntries_data[ipt][ib]);
      //hAss_nent1[ipt]->Scale(1.0/hAss_nent1[ipt]->GetMaximum());
      hAss_nent1[ipt]->SetMarkerStyle(8);
      hAss_nent1[ipt]->SetMarkerColor(ipt+1);
      hAss_nent1[ipt]->SetLineColor(ipt+1);
      if(ipt==4)
	{
	  hAss_nent1[ipt]->SetMarkerColor(6);
	  hAss_nent1[ipt]->SetLineColor(6);
 
	}
      //hResponse1[ipt]->SetBinError(ib+1,dmeanAssym1[ipt][ib]);
      hResponse1[ipt]->SetMarkerStyle(8);
      hResponse1[ipt]->SetMarkerColor(ipt+1);
      hResponse1[ipt]->SetLineColor(ipt+1);
      if(ipt==4)
	{
	  hResponse1[ipt]->SetMarkerColor(6);
	  hResponse1[ipt]->SetLineColor(6);
 
	}
      hResponse1[ipt]->GetYaxis()->SetRangeUser(0.5,1.3);

      hAssym2[ipt][ib]=(TH1F *) fIn2->Get(Form("hAssym_etaBin_%d_%d",ipt,ib));
      nEntries_sim[ipt][ib]=hAssym2[ipt][ib]->GetEntries();
      meanAssym2[ipt][ib]=hAssym2[ipt][ib]->GetMean();
      //dmeanAssym2[ipt][ib]=hAssym2[ipt][ib]->GetRMSError();
      dmeanAssym2[ipt][ib]=hAssym2[ipt][ib]->GetMeanError();
      dstdDev2[ipt][ib]=hAssym1[ipt][ib]->GetStdDev();
    
      if(1)
	{response2[ipt][ib]=(1+meanAssym2[ipt][ib])/(1-meanAssym2[ipt][ib]);}

      /* if(dmeanAssym2[ipt][ib]>0){ */
      /* 	  dresponse2[ipt][ib]=response2[ipt][ib]*dmeanAssym2[ipt][ib]/meanAssym2[ipt][ib]; */
      /* 	} */
      if((1-meanAssym2[ipt][ib])>0){
	dresponse2[ipt][ib]=(2*dmeanAssym2[ipt][ib])/pow((1-meanAssym2[ipt][ib]),2);
	}
      cAs[ipt]->cd(ib+1);
      //cAs[ipt]->SaveAs(Form("histPlots_noAsymCut/MC_Asymmetry_@_\u03B1_%s_%d.pdf",alphaCond,ipt),"RECREATE");
     gPad->SetGridy();
      gPad->SetGridx();
      hAssym2[ipt][ib]->Draw();
      hAssym2[ipt][ib]->SetMarkerStyle(7);
      hAssym2[ipt][ib]->SetMarkerColor(2);
      hAssym2[ipt][ib]->SetLineColor(2);
       //hAssym[ipt][ib]->GetXaxis()->SetRangeUser(0,1);
      
      hResponse2[ipt]->SetBinContent(ib+1,response2[ipt][ib]);
      hResponse2[ipt]->SetBinError(ib+1,dresponse2[ipt][ib]);
      hAss_nent2[ipt]->SetBinContent(ib+1,nEntries_sim[ipt][ib]);
      // hAss_nent2[ipt]->Scale(1.0/hAss_nent2[ipt]->GetMaximum());
      //hAss_nent2[ipt]->Scale(1.0/nEntries_data[4][ib]);
      //hAss_nent2[ipt]->Scale(1.0/hAss_nent2[4]);
      hAss_nent2[ipt]->SetMarkerStyle(8);
      hAss_nent2[ipt]->SetMarkerColor(ipt+1);
      hAss_nent2[ipt]->SetLineColor(ipt+1);
      if(ipt==4)
	{
	  hAss_nent2[ipt]->SetMarkerColor(6);
	  hAss_nent2[ipt]->SetLineColor(6);
 
	}
         //hResponse2[ipt]->SetBinError(ib+1,dmeanAssym2[ipt][ib]);
      hResponse2[ipt]->SetMarkerStyle(8);
      hResponse2[ipt]->SetMarkerColor(ipt+1);
      hResponse2[ipt]->SetLineColor(ipt+1);
      if(ipt==4)
	{
	  hResponse2[ipt]->SetMarkerColor(6);
	  hResponse2[ipt]->SetLineColor(6);
 
	}
      hResponse2[ipt]->GetYaxis()->SetRangeUser(0.5,1.3);
      hDataSimDiff[ipt]=(TH1F *)hResponse2[ipt]->Clone("hDataSimDiff");
      hDataSimDiff[ipt]->SetName(Form("hDataSimDiff_ptBins_%d",ipt));
      hDataSimDiff[ipt]->Divide(hResponse1[ipt]);
      hDataSimDiff[ipt]->SetTitle(Form("MC_Data_ratio @ #alpha < %.2f;#eta;#frac{MC}{Data}",alphaCd));
      hDataSimDiff[ipt]->GetYaxis()->SetRangeUser(0.7,1.35);

      hDataSimDiff_nEntries[ipt]=(TH1F *)hAss_nent2[ipt]->Clone("hDataSimDiff_nEntries");
      hDataSimDiff_nEntries[ipt]->SetName(Form("hDataSimDiff_nEntries_ptBins_%d",ipt));
      hDataSimDiff_nEntries[ipt]->Divide(hAss_nent1[ipt]);
      hDataSimDiff_nEntries[ipt]->SetTitle(Form("MC_Data_ratio @ #alpha < %.2f;#eta;#frac{MC}{Data}",alphaCd));
      hDataSimDiff_nEntries[ipt]->GetYaxis()->SetRangeUser(0.7,1.35);
  
    }

   
  }
  
  TCanvas *cR2= new TCanvas("cR2","cR2",15,10,800,800);
  cR2->Divide(1,2);
  cR2->cd(1);
  gPad->SetGridy();
  gPad->SetGridx();
  hResponse1[0]->Draw();
  hResponse1[1]->Draw("same");
  hResponse1[2]->Draw("same");
  hResponse1[3]->Draw("same");
  hResponse1[4]->Draw("same");
  
  cR2->cd(2);
  gPad->SetGridy();
  gPad->SetGridx();
  hResponse2[0]->Draw();
  hResponse2[1]->Draw("same");
  hResponse2[2]->Draw("same");
  hResponse2[3]->Draw("same");
  hResponse2[4]->Draw("same");

  
  TLegend *leg= new TLegend(0.4,0.12,0.7,0.45);
  leg->AddEntry(hResponse2[0]," 50 < p_{T}_{avg} < 90 ","p");
  leg->AddEntry(hResponse2[1]," 90 < p_{T}_{avg} < 120 ","p");
  leg->AddEntry(hResponse2[2]," 120 < p_{T}_{avg} < 170 ","p");
  leg->AddEntry(hResponse2[3]," 170 < p_{T}_{avg} < 250 ","p");
  leg->AddEntry(hResponse2[4]," 250 < p_{T}_{avg} < 1000 ","p");
  leg->Draw();
  cR2->SaveAs(Form("data_simu_diff_rootfiles_alphaCond/Response_@_\u03B1_%s.pdf",alphaCond),"RECREATE");

  TCanvas *cR3= new TCanvas("cR3","cR3",15,10,750,750);
  hDataSimDiff[0]->Draw();
  hDataSimDiff[1]->Draw("same");
  hDataSimDiff[2]->Draw("same");
  hDataSimDiff[3]->Draw("same");
  hDataSimDiff[4]->Draw("same");
  gPad->SetGridy();
  gPad->SetGridx();
  gPad->SetLeftMargin(0.12); // Moves the plot to the right

  TLegend *leg2= new TLegend(0.35,0.16,0.6,0.38);
  leg2->AddEntry(hDataSimDiff[0]," 50 < p_{T}_{avg} < 90 ","p");
  leg2->AddEntry(hDataSimDiff[1]," 90 < p_{T}_{avg} < 120 ","p");
  leg2->AddEntry(hDataSimDiff[2]," 120 < p_{T}_{avg} < 170 ","p");
  leg2->AddEntry(hDataSimDiff[3]," 170 < p_{T}_{avg} < 250 ","p");
  leg2->AddEntry(hDataSimDiff[4]," 250 < p_{T}_{avg} < 1000 ","p");
  leg2->Draw();
  cR3->SaveAs(Form("data_simu_diff_rootfiles_alphaCond/MC_Data_Ratio_@_\u03B1_%s_Vs_\u03B7.pdf",alphaCond),"RECREATE");
  //leg2->Draw();

  TCanvas *cR4= new TCanvas("cR4","cR4",15,10,800,800);
  cR4->Divide(1,2);
  cR4->cd(1);
  gPad->SetGridy();
  gPad->SetGridx();
  hAss_nent1[2]->Draw("P");
  hAss_nent1[1]->Draw("sameP");
  hAss_nent1[4]->Draw("sameP");
  hAss_nent1[3]->Draw("sameP");
  hAss_nent1[0]->Draw("sameP");
  
  cR4->cd(2);
  gPad->SetGridy();
  gPad->SetGridx();
  hAss_nent2[4]->Draw("P");
  hAss_nent2[1]->Draw("sameP");
  hAss_nent2[2]->Draw("sameP");
  hAss_nent2[3]->Draw("sameP");
  hAss_nent2[0]->Draw("sameP");
  TLegend *leg4= new TLegend(0.75,0.75,0.87,0.88);
  leg4->AddEntry(hAss_nent2[0]," 50 < p_{T}_{avg} < 90 ","pl");
  leg4->AddEntry(hAss_nent2[1]," 90 < p_{T}_{avg} < 120 ","pl");
  leg4->AddEntry(hAss_nent2[2]," 120 < p_{T}_{avg} < 170 ","pl");
  leg4->AddEntry(hAss_nent2[3]," 170 < p_{T}_{avg} < 250 ","pl");
  leg4->AddEntry(hAss_nent2[4]," 250 < p_{T}_{avg} < 1000 ","pl");
  leg4->Draw();
  cR4->SaveAs(Form("data_simu_diff_rootfiles_alphaCond/NEvents_in_Asymmetry_plot_@_\u03B1_%s_Vs_\u03B7.pdf",alphaCond),"RECREATE");
  
#if 1
  TFile *outputFile = new TFile(Form("data_simu_diff_rootfiles_alphaCond/dijet_balance_l2Residual_alphaCond_%s.root", alphaCond), "RECREATE");
			 
  for (int ipt = 0; ipt < npt; ipt++) {
      hDataSimDiff[ipt]->Write();
      hResponse1[ipt]->Write();
      hResponse2[ipt]->Write();
      hAss_nent1[ipt]->Write();
      hAss_nent2[ipt]->Write();
  }
  outputFile->Close();
#endif

}
