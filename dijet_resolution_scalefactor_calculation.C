#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TKey.h>
#include <iostream>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <fstream> // Required for file operations

Double_t fitfunct(Double_t *x, Double_t *par){

  Double_t xx = x[0];

  Double_t fitval = 1.0/(par[0]+par[1]*TMath::Log10(0.01*xx)+(par[2]*(10.0/xx)));

  return fitval;
}

Double_t fitfunctjec(Double_t *x, Double_t *par){

  Double_t xx = x[0];

  Double_t fitval = sqrt(pow(par[0],2)+pow(par[1],2)/xx+pow((par[2]/xx),2)+pow((par[3]/xx),3));

  return fitval;
}

void dijet_resolution_scalefactor_calculation(){
  gROOT->Reset();
  // const int npt= 5;
  gStyle->SetOptStat(0);

#if 0
  ofstream outFile_MC("jer_sigma_rootfiles_skiplowestPtBin/MC_sigma_resolution_lowPU_pp13TeV_with_newStatistics.txt");
  if (!outFile_MC.is_open()) {
    cerr << "Error: Could not open file for writing." << endl;
    return;
  }
  outFile_MC << "{1 JetEta 1 JetPt sqrt(pow(p0,2)+pow(p1,2)/x+pow((p2/x),2)) MC_resolution_Correction}\n";
#endif
#if 1
  
ofstream outFile_RES("jer_sigma_rootfiles_skiplowestPtBin/JER_resolution_lowPU_pp13TeV_with_newStatistics_SF.txt");
  if (!outFile_RES.is_open()) {
    cerr << "Error: Could not open file for writing." << endl;
    return;
  }
  outFile_RES << "{1 JetEta 0 None ScaleFactor}\n";
#endif
  
  const int nfile= 5;
  TFile *fIn[nfile];
  float alphaBin[]={0.0,0.15,0.20,0.25,0.30,0.35};
  const int nEta = 50;
  float etaBin[] = {-5.191,  -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.051, -1.93, -1.791, -1.653, -1.566, -1.479, -1.392, -1.305, -1.174, -1.044, -0.913, -0.783, -0.652, -0.522, -0.391, -0.261, -0.130, 0.0, 0.130, 0.261, 0.391, 0.522, 0.652, 0.783, 0.913, 1.044, 1.174, 1.305, 1.392, 1.479, 1.566, 1.653, 1.791, 1.93, 2.051, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 5.191};
float etaBin2[] = {-5.20,  -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.051, -1.93, -1.791, -1.653, -1.566, -1.479, -1.392, -1.305, -1.174, -1.044, -0.913, -0.783, -0.652, -0.522, -0.391, -0.261, -0.130, 0.0, 0.130, 0.261, 0.391, 0.522, 0.652, 0.783, 0.913, 1.044, 1.174, 1.305, 1.392, 1.479, 1.566, 1.653, 1.791, 1.93, 2.051, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 5.20};

  const int npt = 4;
  double avgPtBin[]={90.0,120.0,170.0,250.0,1000.0};
 
  //double avgPtBin[]={50.0,90.0,120.0,160.0,200.0,500.0};
  const int npt_4b = 4;
  double avgPtBin2[]={90.0,120.0,170.0,250.0,500.0};
  double avgPtBin3[]={90.0,120.0,170.0,250.0,500.0};



#if 1

  //sigma from the gaussian fit//
  fIn[0]= new TFile("jer_sigma_rootfiles_skiplowestPtBin/dijet_jer_standard_sigma_plot_alphaCond_newStat_fit015.root");
  fIn[1]= new TFile("jer_sigma_rootfiles_skiplowestPtBin/dijet_jer_standard_sigma_plot_alphaCond_newStat_fit020.root");
  fIn[2]= new TFile("jer_sigma_rootfiles_skiplowestPtBin/dijet_jer_standard_sigma_plot_alphaCond_newStat_fit025.root");
  fIn[3]= new TFile("jer_sigma_rootfiles_skiplowestPtBin/dijet_jer_standard_sigma_plot_alphaCond_newStat_fit030.root");
  fIn[4]= new TFile("jer_sigma_rootfiles_skiplowestPtBin/dijet_jer_standard_sigma_plot_alphaCond_newStat_fit035.root");
#endif
 
  
  TH1F *hDatPt[nfile][npt];
  TH1F *hSimPt[nfile][npt];
  /* const int nEta = 22; */
  /* float etaBin[] = {-5.0, -3.5, -3.0, -2.7, -2.4, -2.0, -1.7, -1.3, -1.0, -0.7, -0.4, 0.0, 0.4, 0.7, 1.0, 1.3, 1.7, 2.0, 2.4, 2.7, 3.0, 3.5, 5.0}; */
  /* float etaBin2[] = {-5.20, -3.5, -3.0, -2.7, -2.4, -2.0, -1.7, -1.3, -1.0, -0.7, -0.4, 0.0, 0.4, 0.7, 1.0, 1.3, 1.7, 2.0, 2.4, 2.7, 3.0, 3.5, 5.20}; */

  
  
  float dat_per_alphaVal_perpt_perEtaval[nfile][npt][nEta];
  float ddat_per_alphaVal_perpt_perEtaval[nfile][npt][nEta];

  float Sim_per_alphaVal_perpt_perEtaval[nfile][npt][nEta];
  float dSim_per_alphaVal_perpt_perEtaval[nfile][npt][nEta];

  TH1F *hdat_sim_vs_alpha[npt][nEta];
  float nAlphaValue_dat[npt][nEta][nfile];
  float dnAlphaValue_dat[npt][nEta][nfile];
  float nAlphaValue_sim[npt][nEta][nfile];
  float dnAlphaValue_sim[npt][nEta][nfile];

  for( int iA=0; iA<nfile;iA++){
    for( int ipt=0; ipt<npt;ipt++){
      hDatPt[iA][ipt]=(TH1F *) fIn[iA]->Get(Form("hsigma_data_%d",ipt+1));
      hSimPt[iA][ipt]=(TH1F *) fIn[iA]->Get(Form("hsigma_simu_%d",ipt+1));
      for(int ita=0;ita<hDatPt[iA][ipt]->GetNbinsX();ita++){
	//if(ita==0 ||ita==1||ita==20||ita==21) continue;
	dat_per_alphaVal_perpt_perEtaval[iA][ipt][ita]=hDatPt[iA][ipt]->GetBinContent(ita+1);
	ddat_per_alphaVal_perpt_perEtaval[iA][ipt][ita]=hDatPt[iA][ipt]->GetBinError(ita+1);
 	Sim_per_alphaVal_perpt_perEtaval[iA][ipt][ita]=hSimPt[iA][ipt]->GetBinContent(ita+1);
	dSim_per_alphaVal_perpt_perEtaval[iA][ipt][ita]=hSimPt[iA][ipt]->GetBinError(ita+1);
      }


    }
    

  }

  for( int ipt=0; ipt<npt;ipt++){
    for(int ita=0;ita<nEta;ita++){
      for( int iA=0; iA<nfile;iA++){
	//if(ita==0 ||ita==1||ita==20||ita==21) continue;

	nAlphaValue_dat[ipt][ita][iA]=dat_per_alphaVal_perpt_perEtaval[iA][ipt][ita];
	dnAlphaValue_dat[ipt][ita][iA]=ddat_per_alphaVal_perpt_perEtaval[iA][ipt][ita];
	nAlphaValue_sim[ipt][ita][iA]=Sim_per_alphaVal_perpt_perEtaval[iA][ipt][ita];
	dnAlphaValue_sim[ipt][ita][iA]=dSim_per_alphaVal_perpt_perEtaval[iA][ipt][ita];
      }
    }
  }

  float dat_alpha0[npt][nEta];
  float ddat_alpha0[npt][nEta];
  float sim_alpha0[npt][nEta];
  float dsim_alpha0[npt][nEta];

  //float dat_sim_ratio_alpha0[npt][nEta];
  //float ddat_sim_ratio_alpha0[npt][nEta];

  // Fit the graph with a linear function (polynomial of degree 1)
  TF1* fitFunc1[npt][nEta];
  TF1* fitFunc2[npt][nEta];
   
  TCanvas *cAd[npt];
  TCanvas *cAs[npt];
  // TCanvas *cA[npt];

  //#if 1
 
  TGraphErrors* graph_vs_alpha_data[npt][nEta];
  TGraphErrors* graph_vs_alpha_sim[npt][nEta];
  //const int totBins=6;
  //float Xval[]={0.075,0.10,0.15,0.20,0.25,0.30};
  const int totBins=5;
  float Xval[]={0.15,0.20,0.25,0.30,0.35};
  float dataSim_ratio_for_same_eta_and_pt_but_diff_alpha[npt][nEta][totBins];
  float ddataSim_ratio_for_same_eta_and_pt_but_diff_alpha[npt][nEta][totBins];
  TH1F *hdata_alpha0_per_pt[npt];
  TH1F *hsim_alpha0_per_pt[npt];
  TH1F *hdata_simRatio_alpha0_Vs_Pt_perEta[nEta];
  for( int ipt=0; ipt<npt;ipt++){
    cAd[ipt]= new TCanvas(Form("cAd%d",ipt),Form("cAd%d",ipt),10,10,1200,800);
    cAd[ipt]->Divide(5,10);
    cAs[ipt]= new TCanvas(Form("cAs%d",ipt),Form("cAs%d",ipt),10,10,1200,800);
    cAs[ipt]->Divide(5,10);
    //cA[ipt]= new TCanvas(Form("cA%d",ipt),Form("cA%d",ipt),10,10,1200,800);
    //   cA[ipt]->Divide(5,10);

    hdata_alpha0_per_pt[ipt]=(TH1F*)hDatPt[0][ipt]->Clone();
    hdata_alpha0_per_pt[ipt]->Reset();
    hdata_alpha0_per_pt[ipt]->SetTitle(Form("hdata_alpha0_per_pt_for_zero_alpha_%.1f < avg_p_{T} < %.1f;#eta;#sigma_{A}_{#alpha #rightarrow 0}",avgPtBin[ipt],avgPtBin[ipt+1]));

    hsim_alpha0_per_pt[ipt]=(TH1F*)hSimPt[0][ipt]->Clone();
    hsim_alpha0_per_pt[ipt]->Reset();
    hsim_alpha0_per_pt[ipt]->SetTitle(Form("hsim_alpha0_per_pt_for_zero_alpha_%.1f < avg_p_{T} < %.1f;#eta;#sigma_{A}_{#alpha #rightarrow 0}",avgPtBin[ipt],avgPtBin[ipt+1]));


    
    for(int ita=0;ita<nEta;ita++){
      // if(ita==0 ||ita==1||ita==20||ita==21) continue;
      if(ita==0 ||ita==1||ita==2||ita==3||ita==46||ita==47||ita==48||ita==49) continue;

      cAd[ipt]->cd(ita+1);
      // if(ita==0&&ita==21)continue;
      graph_vs_alpha_data[ipt][ita] = new TGraphErrors(totBins, Xval, nAlphaValue_dat[ipt][ita],0,dnAlphaValue_dat[ipt][ita]);
      graph_vs_alpha_data[ipt][ita]->SetTitle(Form("Ratio_vs_#alpha@(pt %d, eta %d);#alpha;#sigma_{A}", ipt, ita));

      graph_vs_alpha_sim[ipt][ita] = new TGraphErrors(totBins, Xval, nAlphaValue_sim[ipt][ita],0,dnAlphaValue_sim[ipt][ita]);
      graph_vs_alpha_sim[ipt][ita]->SetTitle(Form("Ratio_vs_#alpha@(pt %d, eta %d);#alpha;#sigma_{A}", ipt, ita));
      /* if(ita==0 || ita ==18 ||ita ==19 ||ita==21) */
      /* 	{ */
      /* 	  graph_vs_alpha[ipt][ita]->GetYaxis()->SetRangeUser(0.75,1.5); */
      /* 	} */
      /* else */
      graph_vs_alpha_data[ipt][ita]->GetYaxis()->SetRangeUser(0.0,0.35);
      
      graph_vs_alpha_data[ipt][ita]->SetMarkerStyle(8);
      graph_vs_alpha_data[ipt][ita]->SetMarkerColor(2);
      graph_vs_alpha_data[ipt][ita]->SetLineColor(2);
      graph_vs_alpha_data[ipt][ita]->Draw("AP");
      fitFunc1[ipt][ita]= new TF1(Form("fitFunc1_%d_%d",ipt,ita), "pol1", 0.05, 0.35);
      fitFunc1[ipt][ita]->SetLineColor(2);
      graph_vs_alpha_data[ipt][ita]->Fit(fitFunc1[ipt][ita],"QR");
      cAs[ipt]->cd(ita+1);

      /* if(ipt==3 || ipt==4) */
      /* 	{ */
	
      /* 	  if(ita==0 ||ita==21) */
      /* 	    { */
      /* 	      graph_vs_alpha_sim[ipt][ita]->GetYaxis()->SetRangeUser(0.0,1.5); */
      /* 	    } */
      /* 	  else graph_vs_alpha_sim[ipt][ita]->GetYaxis()->SetRangeUser(0.0,0.25); */
      /* 	} */
      /* else */
      graph_vs_alpha_sim[ipt][ita]->GetYaxis()->SetRangeUser(0.0,0.25);
      graph_vs_alpha_sim[ipt][ita]->SetMarkerStyle(23);
      graph_vs_alpha_sim[ipt][ita]->SetMarkerColor(4);
      graph_vs_alpha_sim[ipt][ita]->SetLineColor(4);

      fitFunc2[ipt][ita]= new TF1(Form("fitFunc2_%d_%d",ipt,ita), "pol1", 0.05, 0.35);
      fitFunc2[ipt][ita]->SetLineColor(4);
      graph_vs_alpha_sim[ipt][ita]->Fit(fitFunc2[ipt][ita],"QR");
      graph_vs_alpha_sim[ipt][ita]->Draw("AP");
      /* cA[ipt]->cd(ita+1); */
      /* graph_vs_alpha_data[ipt][ita]->Draw("AP"); */
      /* graph_vs_alpha_sim[ipt][ita]->Draw("SameAP"); */
  

      // Store the fit result for data
      TFitResultPtr fitResult1 = graph_vs_alpha_data[ipt][ita]->Fit(fitFunc1[ipt][ita], "QS"); // "QS" for quiet mode, storing fit result

      double x1_value = 0;
      //getvalue for alpha==0 rom the pol1 fit.//
      dat_alpha0[ipt][ita]=fitFunc1[ipt][ita]->GetParameter(0);
      sim_alpha0[ipt][ita]=fitFunc2[ipt][ita]->GetParameter(0);
       
      
      // Get the intercept and its error
      double p01 = fitFunc1[ipt][ita]->GetParameter(0);      // Intercept
      double p01_error = fitFunc1[ipt][ita]->GetParError(0); // Uncertainty of intercept

      // Get the slope and its error
      double p11 = fitFunc1[ipt][ita]->GetParameter(1);      // Slope
      double p11_error = fitFunc1[ipt][ita]->GetParError(1); // Uncertainty of slope

      // Get the covariance matrix
      double covariance1 = fitResult1->CovMatrix(0, 1); // Covariance between p0 and p1

      // Calculate the error at x = 0
      // double error_at_x01 = sqrt(p01_error * p01_error + x1_value * x1_value * p11_error * p11_error + 2 * x1_value * covariance1);
      ddat_alpha0[ipt][ita] = sqrt(p01_error * p01_error + x1_value * x1_value * p11_error * p11_error + 2 * x1_value * covariance1);
      //#########################################################################################################
      //#########################################################################################################

      // Store the fit result for simu//
      TFitResultPtr fitResult2 = graph_vs_alpha_sim[ipt][ita]->Fit(fitFunc2[ipt][ita], "QS"); // "QS" for quiet mode, storing fit result

      double x2_value = 0;
     
      // Get the intercept and its error
      double p02 = fitFunc2[ipt][ita]->GetParameter(0);      // Intercept
      double p02_error = fitFunc2[ipt][ita]->GetParError(0); // Uncertainty of intercept

      // Get the slope and its error
      double p12 = fitFunc2[ipt][ita]->GetParameter(1);      // Slope
      double p12_error = fitFunc2[ipt][ita]->GetParError(1); // Uncertainty of slope

      // Get the covariance matrix
      double covariance2 = fitResult2->CovMatrix(0, 1); // Covariance between p0 and p1

      // Calculate the error at x = 0
      // double error_at_x01 = sqrt(p01_error * p01_error + x1_value * x1_value * p11_error * p11_error + 2 * x1_value * covariance1);
      dsim_alpha0[ipt][ita] = sqrt(p02_error * p02_error + x2_value * x2_value * p12_error * p12_error + 2 * x2_value * covariance2);
         
      //#########################################################################################################
      //#########################################################################################################   //FIll the histogram with the ratio of the MC and Data when alpha==0//
#if 1   
      hdata_alpha0_per_pt[ipt]->SetBinContent(ita+1, dat_alpha0[ipt][ita]);
      hdata_alpha0_per_pt[ipt]->SetBinError(ita+1, ddat_alpha0[ipt][ita]);
      hsim_alpha0_per_pt[ipt]->SetBinContent(ita+1, sim_alpha0[ipt][ita]);
      hsim_alpha0_per_pt[ipt]->SetBinError(ita+1, dsim_alpha0[ipt][ita]);

      // hdata_simRatio_alpha0[ipt]->SetBinError(ita+1, error_at_x0);//ddat_sim_ratio_alpha0[ipt][ita]);


	
      hdata_alpha0_per_pt[ipt]->SetMarkerStyle(8);
      hdata_alpha0_per_pt[ipt]->SetMarkerColor(ipt+1);
      hdata_alpha0_per_pt[ipt]->SetLineColor(ipt+1);
      if(ipt==4)
	{
	  hdata_alpha0_per_pt[ipt]->SetMarkerColor(6);
	  hdata_alpha0_per_pt[ipt]->SetLineColor(6);
	     
	}
      hdata_alpha0_per_pt[ipt]->GetYaxis()->SetRangeUser(0.0,0.5);

      hsim_alpha0_per_pt[ipt]->SetMarkerStyle(8);
      hsim_alpha0_per_pt[ipt]->SetMarkerColor(ipt+1);
      hsim_alpha0_per_pt[ipt]->SetLineColor(ipt+1);
      if(ipt==4)
	{
	  hsim_alpha0_per_pt[ipt]->SetMarkerColor(6);
	  hsim_alpha0_per_pt[ipt]->SetLineColor(6);
	     
	}
      hsim_alpha0_per_pt[ipt]->GetYaxis()->SetRangeUser(0.0,0.5);

    
      //cA[ipt]->SaveAs(Form("histPlots/MC_Data_Ratio_vs_\u03B1_%.2f_to_%.2f.pdf",avgPtBin[ipt],avgPtBin[ipt+1]),"RECREATE");
#endif
    }
      
  }
  TLegend *leg= new TLegend(0.35,0.6,0.75,0.85);//0.7,0.9,0.89,0.99);
  //leg->AddEntry(hsim_alpha0_per_pt[0]," 50.0 < p_{T}_{avg} < 90.0 ","pl");
  leg->AddEntry(hsim_alpha0_per_pt[0]," 90.0 < p_{T}_{avg} < 120.0 ","pl");
  leg->AddEntry(hsim_alpha0_per_pt[1]," 120.0 < p_{T}_{avg} < 170.0 ","pl");
  leg->AddEntry(hsim_alpha0_per_pt[2]," 170.0 < p_{T}_{avg} < 250.0 ","pl");
  leg->AddEntry(hsim_alpha0_per_pt[3]," 250.0 < p_{T}_{avg} < 1000.0 ","pl");
  TCanvas *cRR=new TCanvas();
  cRR->Divide(1,2);
  gPad->SetGridx();
  gPad->SetGridy();
  cRR->cd(1);
  hdata_alpha0_per_pt[0]->Draw();//"P");
  hdata_alpha0_per_pt[1]->Draw("same");//P");
  hdata_alpha0_per_pt[2]->Draw("same");//P");
  hdata_alpha0_per_pt[3]->Draw("same");//P");
  //hdata_alpha0_per_pt[4]->Draw("same");//P");
  leg->Draw();
  cRR->cd(2);
  hsim_alpha0_per_pt[0]->Draw();//"P");
  hsim_alpha0_per_pt[1]->Draw("same");//P");
  hsim_alpha0_per_pt[2]->Draw("same");//P");
  hsim_alpha0_per_pt[3]->Draw("same");//P");
  //hsim_alpha0_per_pt[4]->Draw("same");//P");
  
  /* cRR->SaveAs(Form("histPlots/MC_Data_Ratio_@_\u03B1_0_VS_\u03B7.pdf"),"RECREATE"); */

  /* leg->Draw(); */


  /*  // Plot graph of Ma/Data@ alpha=0 vs Pt for each ita bins// */
  float dat_alpha0_perEta_Vs_pT[nEta][npt];
  float ddat_alpha0_perEta_Vs_pT[nEta][npt];
  float sim_alpha0_perEta_Vs_pT[nEta][npt];
  float dsim_alpha0_perEta_Vs_pT[nEta][npt];
  TH1F *hdata_alpha0_Vs_Pt_perEta[nEta];
  TH1F *hsim_alpha0_Vs_Pt_perEta[nEta];

  float dat_alpha0_perEta_Vs_pT_skip_Pt0[nEta][npt_4b];
  float ddat_alpha0_perEta_Vs_pT_skip_Pt0[nEta][npt_4b];
  float sim_alpha0_perEta_Vs_pT_skip_Pt0[nEta][npt_4b];
  float dsim_alpha0_perEta_Vs_pT_skip_Pt0[nEta][npt_4b];
  TH1F *hdata_alpha0_Vs_Pt_perEta_skip_Pt0[nEta];
  TH1F *hsim_alpha0_Vs_Pt_perEta_skip_Pt0[nEta];

    
  TF1 *fit01[nEta];
  TCanvas *cPP= new TCanvas(Form("cPP"),Form("cPP"),10,10,1200,800);
  cPP->Divide(5,10);

  // Eliminating first ptBin (50-90) due to statistical issue//
  
  for(int ita=0;ita<nEta;ita++){
    if(ita==0 ||ita==1||ita==2||ita==3||ita==46||ita==47||ita==48||ita==49) continue;

    hdata_alpha0_Vs_Pt_perEta[ita]= new TH1F(Form("hdata_alpha0_Vs_Pt_perEta_%.1f_%.1f",etaBin[ita],etaBin[ita+1]),Form("hdata_alpha0_Vs_Pt_perEta_%.1f_%.1f",etaBin[ita],etaBin[ita+1]),npt,avgPtBin2);
    hsim_alpha0_Vs_Pt_perEta[ita]= new TH1F(Form("hsim_alpha0_Vs_Pt_perEta_%.1f_%.1f",etaBin[ita],etaBin[ita+1]),Form("hsim_alpha0_Vs_Pt_perEta_%.1f_%.1f",etaBin[ita],etaBin[ita+1]),npt,avgPtBin2);
    for( int ipt=0; ipt<npt;ipt++){
      //cout<<ita<<'\t'<<ipt<<'\t'<<dat_alpha0_perEta_Vs_pT[ita][ipt]<<'\t'<<ddat_alpha0_perEta_Vs_pT[ita][ipt]<<endl;
      dat_alpha0_perEta_Vs_pT[ita][ipt]=dat_alpha0[ipt][ita];
      ddat_alpha0_perEta_Vs_pT[ita][ipt]=ddat_alpha0[ipt][ita];
      sim_alpha0_perEta_Vs_pT[ita][ipt]=sim_alpha0[ipt][ita];
      dsim_alpha0_perEta_Vs_pT[ita][ipt]=dsim_alpha0[ipt][ita];
 
      hdata_alpha0_Vs_Pt_perEta[ita]->SetBinContent(ipt+1,dat_alpha0_perEta_Vs_pT[ita][ipt]);
      hdata_alpha0_Vs_Pt_perEta[ita]->SetBinError(ipt+1,ddat_alpha0_perEta_Vs_pT[ita][ipt]);
      hdata_alpha0_Vs_Pt_perEta[ita]->SetMarkerStyle(8);
      hdata_alpha0_Vs_Pt_perEta[ita]->SetMarkerColor(2);
      hdata_alpha0_Vs_Pt_perEta[ita]->SetLineColor(2);
      hdata_alpha0_Vs_Pt_perEta[ita]->SetTitle(";p_{T}^{avg};#sigma_{#alpha #rightarrow 0}");
      hdata_alpha0_Vs_Pt_perEta[ita]->GetYaxis()->SetRangeUser(0.0,0.35);

      hsim_alpha0_Vs_Pt_perEta[ita]->SetBinContent(ipt+1,sim_alpha0_perEta_Vs_pT[ita][ipt]);
      hsim_alpha0_Vs_Pt_perEta[ita]->SetBinError(ipt+1,dsim_alpha0_perEta_Vs_pT[ita][ipt]);
      hsim_alpha0_Vs_Pt_perEta[ita]->SetMarkerStyle(8);
      hsim_alpha0_Vs_Pt_perEta[ita]->SetMarkerColor(4);
      hsim_alpha0_Vs_Pt_perEta[ita]->SetLineColor(4);
      hsim_alpha0_Vs_Pt_perEta[ita]->SetTitle(";p_{T}^{avg};#sigma_{#alpha #rightarrow 0}");
      hsim_alpha0_Vs_Pt_perEta[ita]->GetYaxis()->SetRangeUser(0.0,0.35);



    }
    cPP->cd(ita+1);
    hdata_alpha0_Vs_Pt_perEta[ita]->Draw();
    hsim_alpha0_Vs_Pt_perEta[ita]->Draw("same");
    TLegend *leg1= new TLegend(0.55,0.6,0.75,0.85);//0.7,0.9,0.89,0.99);
    leg1->AddEntry(hdata_alpha0_Vs_Pt_perEta[ita]," Data ","pl");
    leg1->AddEntry(hsim_alpha0_Vs_Pt_perEta[ita]," MC ","pl");
    leg1->Draw();

	 
  }

  TCanvas *cPP11= new TCanvas();
  hdata_alpha0_Vs_Pt_perEta[15]->Draw();
  hsim_alpha0_Vs_Pt_perEta[15]->Draw("same");
  ///////////////////////////////////////////////////////////////////////////////////
  TCanvas *cRA= new TCanvas(Form("cRA"),Form("cRA"),10,10,1200,800);
  cRA->Divide(5,10);
  float dat_sim_ratio_alpha0[nEta];
  float ddat_sim_ratio_alpha0[nEta];
  TH1F * hAbs_data_simRatio_alpha0_allbin = new TH1F(Form("hAbs_data_simRatio_alpha_allbin"), Form("hAbs_data_simRatio_alpha_allbin"), nEta, etaBin);//0, 5.0);

  for(int ita=0;ita<nEta;ita++){
    //if(ita==0 ||ita==1||ita==20||ita==21) continue;
    if(ita==0 ||ita==1||ita==2||ita==3||ita==46||ita==47||ita==48||ita==49) continue;
    cRA->cd(ita+1);
      
    fit01[ita]=new TF1(Form("fit01%d",ita),fitfunctjec,50,500,4);
    fit01[ita]->SetLineColor(2);
    fit01[ita]->SetParameter(0,1.50);
    fit01[ita]->SetParameter(1,0.50);
    fit01[ita]->SetParameter(2,-2.50);
    fit01[ita]->SetParameter(3,1.50);
    hdata_simRatio_alpha0_Vs_Pt_perEta[ita]= (TH1F*)hdata_alpha0_Vs_Pt_perEta[ita]->Clone("hdata_simRatio_alpha0_Vs_Pt_perEta");//Form("hdata_simRatio_alpha0_Vs_Pt_perEta_%d",ita));
    hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->SetName(Form("hdata_simRatio_alpha0_Vs_Pt_perEta_%d",ita));

    hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->Divide(hsim_alpha0_Vs_Pt_perEta[ita]);//Taking Datato MC ratio//
    
    hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->SetMarkerStyle(8);
    hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->SetMarkerColor(4);
    hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->SetLineColor(4);
    hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->Draw();
    hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->GetYaxis()->SetRangeUser(0.2,2.0);
    
    /* hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->Fit(fit01[ita],"QR"); */
    /* dat_sim_ratio_alpha0[ita]=fit01[ita]->GetParameter(0); */
    /* ddat_sim_ratio_alpha0[ita]=fit01[ita]->GetParError(0); */

    //hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->Fit("pol0");
    // Perform the fit
    // Perform the fit and get a TFitResultPtr
    TFitResultPtr fitResult = hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->Fit("pol0", "S");

    // Check if the fit was successful
    if (fitResult->IsValid()) {
      // Get the TF1 object from the histogram
      TF1 *fitFunction = hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->GetFunction("pol0");
      if (fitFunction) {
        // Get the value of p0 (constant term)
        double p0 = fitFunction->GetParameter(0);
        std::cout << "p0 = " << p0 << std::endl;

        // Get the uncertainty of p0
        double p0Error = fitFunction->GetParError(0);
	dat_sim_ratio_alpha0[ita]=fitFunction->GetParameter(0);//p0;//fit01[ita]->GetParameter(0);
	ddat_sim_ratio_alpha0[ita]=fitFunction->GetParError(0);//p0Error;//fit01[ita]->GetParError(0);

        std::cout << "Error on p0 = " << p0Error << std::endl;
      } else {
        std::cerr << "Fit function 'pol0' not found!" << std::endl;
      }
    } else {
      std::cerr << "Fit failed!" << std::endl;
    }
    
    hAbs_data_simRatio_alpha0_allbin->SetBinContent(ita+1,dat_sim_ratio_alpha0[ita]);
    hAbs_data_simRatio_alpha0_allbin->SetBinError(ita+1,ddat_sim_ratio_alpha0[ita]);
    //cRA->SaveAs(Form("histPlots/MC_Data_Ratio_@_\u03B1_0_VS_pT_per_\u03B7.pdf"),"RECREATE");
    outFile_RES << etaBin2[ita] <<'\t'
	     << etaBin2[ita+1] <<'\t'
	     << 3 <<'\t'
	     << dat_sim_ratio_alpha0[ita] <<'\t'
	     << dat_sim_ratio_alpha0[ita]-ddat_sim_ratio_alpha0[ita]<<'\t'
	     << dat_sim_ratio_alpha0[ita]+ddat_sim_ratio_alpha0[ita] <<'\n';
    
  }
  outFile_RES.close();
  TCanvas *cRat=new TCanvas();
  hdata_simRatio_alpha0_Vs_Pt_perEta[19]->Draw();
  
  //outFile.close();
  TCanvas *cRr2=new TCanvas();
  hAbs_data_simRatio_alpha0_allbin->SetMarkerStyle(8);
  hAbs_data_simRatio_alpha0_allbin->SetMarkerColor(2);
  hAbs_data_simRatio_alpha0_allbin->SetLineColor(2);
  hAbs_data_simRatio_alpha0_allbin->GetYaxis()->SetRangeUser(0.0,3.0);
  //hAbs_data_simRatio_alpha0_allbin->SetTitle(";#eta;#sigma_{#alpha #rightarrow 0}");
  hAbs_data_simRatio_alpha0_allbin->SetTitle(";#eta;Scale Factor");
  hAbs_data_simRatio_alpha0_allbin->Draw();
  gPad->SetGridx();
  gPad->SetGridy();

  //float AvetaBin[] = {0.0, 0.4, 0.7, 1.0, 1.3, 1.7, 2.0, 2.4, 2.7, 3.0, 3.5, 5.0};
  float AvetaBin[] = {0.0, 0.130, 0.261, 0.391, 0.522, 0.652, 0.783, 0.913, 1.044, 1.174, 1.305, 1.392, 1.479, 1.566, 1.653, 1.791, 1.93, 2.051, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 5.191};
  float nval_ratio_etabin[nEta];
  float dnval_ratio_etabin[nEta];
  TH1F * hAbs_data_simRatio_alpha0 = new TH1F(Form("hAbs_data_simRatio_alpha"), Form("hAbs_data_simRatio_alpha"), nEta / 2, AvetaBin);//0, 5.0);
  //hAbs_data_simRatio_alpha0->SetTitle(";|#eta|;#sigma_{#alpha #rightarrow 0}");
  hAbs_data_simRatio_alpha0->SetTitle(";|#eta|;Scale Factor");
  for (int ib = 0; ib < nEta / 2; ib++) {
   
    int posBin = nEta / 2 + ib;  // Positive eta bin
    int negBin = nEta / 2 - ib - 1;  // Negative eta bin
    if(ib==21||ib==22|| ib==23||ib==24)continue;

    nval_ratio_etabin[posBin]=dat_sim_ratio_alpha0[posBin];
    nval_ratio_etabin[negBin]=dat_sim_ratio_alpha0[negBin];
    dnval_ratio_etabin[posBin]=ddat_sim_ratio_alpha0[posBin];
    dnval_ratio_etabin[negBin]=ddat_sim_ratio_alpha0[negBin];
    float avg_ratio=( nval_ratio_etabin[posBin]+ nval_ratio_etabin[negBin])/2.0;
    float davg_ratio=avg_ratio*sqrt(pow(dnval_ratio_etabin[posBin]/nval_ratio_etabin[posBin],2)+pow(dnval_ratio_etabin[negBin]/nval_ratio_etabin[negBin],2));
    cout<<nval_ratio_etabin[negBin]<<'\t'<<nval_ratio_etabin[posBin]<<'\t'<<avg_ratio<<endl;
    hAbs_data_simRatio_alpha0->SetBinContent(ib+1,avg_ratio);
    hAbs_data_simRatio_alpha0->SetBinError(ib+1,davg_ratio);
    hAbs_data_simRatio_alpha0->SetMarkerStyle(8);
    hAbs_data_simRatio_alpha0->SetMarkerColor(2);
    hAbs_data_simRatio_alpha0->SetLineColor(2);
    hAbs_data_simRatio_alpha0->GetYaxis()->SetRangeUser(0.0,3.0);

  }
  TCanvas *cR2=new TCanvas();
  hAbs_data_simRatio_alpha0->Draw();
  gPad->SetGridx();
  gPad->SetGridy();







}
