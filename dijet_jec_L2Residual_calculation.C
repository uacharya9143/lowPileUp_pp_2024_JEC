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
void dijet_jec_L2Residual_calculation(){
  gROOT->Reset();
  const int npt= 5;
  gStyle->SetOptStat(0);

  ofstream outFile("data_simu_diff_rootfiles_alphaCond/L2Residual_lowPU_pp13TeV_with_newStatistics.txt");
  if (!outFile.is_open()) {
    cerr << "Error: Could not open file for writing." << endl;
    return;
  }
  outFile << "{1 JetEta 1 JetPt 1./([p0]+[p1]*log10(0.01*x)+[p2]/(x/10.0)) Correction L2Residual}\n";
  
   const int nfile= 5;
   TFile *fIn[nfile];
   //float alphaBin[]={0.0,0.10,0.15,0.20,0.25,0.30};
  //fIn[0]= new TFile("alphaCond/dijet_balance_alphaCond_0075.root");
  /* fIn[0]= new TFile("alphaCond/dijet_balance_alphaCond_010_newMethod.root"); */
  /* fIn[1]= new TFile("alphaCond/dijet_balance_alphaCond_015_newMethod.root"); */
  /* fIn[2]= new TFile("alphaCond/dijet_balance_alphaCond_020_newMethod.root"); */
  /* fIn[3]= new TFile("alphaCond/dijet_balance_alphaCond_025_newMethod.root"); */
  /* fIn[4]= new TFile("alphaCond/dijet_balance_alphaCond_030_newMethod.root"); */
   float alphaBin[]={0.0,0.15,0.20,0.25,0.30,0.35};
  fIn[0]= new TFile("data_simu_diff_rootfiles_alphaCond/dijet_balance_l2Residual_alphaCond_015.root");
  fIn[1]= new TFile("data_simu_diff_rootfiles_alphaCond/dijet_balance_l2Residual_alphaCond_020.root");
  fIn[2]= new TFile("data_simu_diff_rootfiles_alphaCond/dijet_balance_l2Residual_alphaCond_025.root");
  fIn[3]= new TFile("data_simu_diff_rootfiles_alphaCond/dijet_balance_l2Residual_alphaCond_030.root");
  fIn[4]= new TFile("data_simu_diff_rootfiles_alphaCond/dijet_balance_l2Residual_alphaCond_035.root");
 
  
  TH1F *hDatSimPt[nfile][npt];
  /* const int nEta = 22; */
  /* float etaBin[] = {-5.0, -3.5, -3.0, -2.7, -2.4, -2.0, -1.7, -1.3, -1.0, -0.7, -0.4, 0.0, 0.4, 0.7, 1.0, 1.3, 1.7, 2.0, 2.4, 2.7, 3.0, 3.5, 5.0}; */
  /* float etaBin2[] = {-5.20, -3.5, -3.0, -2.7, -2.4, -2.0, -1.7, -1.3, -1.0, -0.7, -0.4, 0.0, 0.4, 0.7, 1.0, 1.3, 1.7, 2.0, 2.4, 2.7, 3.0, 3.5, 5.20}; */

  const int nEta = 50;
  float etaBin[] = {-5.191,  -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.051, -1.93, -1.791, -1.653, -1.566, -1.479, -1.392, -1.305, -1.174, -1.044, -0.913, -0.783, -0.652, -0.522, -0.391, -0.261, -0.130, 0.0, 0.130, 0.261, 0.391, 0.522, 0.652, 0.783, 0.913, 1.044, 1.174, 1.305, 1.392, 1.479, 1.566, 1.653, 1.791, 1.93, 2.051, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 5.191};
  float etaBin2[] = {-5.20,  -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.051, -1.93, -1.791, -1.653, -1.566, -1.479, -1.392, -1.305, -1.174, -1.044, -0.913, -0.783, -0.652, -0.522, -0.391, -0.261, -0.130, 0.0, 0.130, 0.261, 0.391, 0.522, 0.652, 0.783, 0.913, 1.044, 1.174, 1.305, 1.392, 1.479, 1.566, 1.653, 1.791, 1.93, 2.051, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 5.20};

  //float nAlphaBin=7;
  float datSim_ratio_per_alphaVal_perpt_perEtaval[nfile][npt][nEta];
  float ddatSim_ratio_per_alphaVal_perpt_perEtaval[nfile][npt][nEta];
  //float dataSim_ratio_for_same_eta_and_pt_but_diff_alpha[nfile];
  TH1F *hdat_sim_vs_alpha[npt][nEta];
  float nAlphaValue[npt][nEta][nfile];
  float dnAlphaValue[npt][nEta][nfile];
  float nAlphaValue_ratio_per_alpha030[npt][nEta][nfile];
  float dnAlphaValue_ratio_per_alpha030[npt][nEta][nfile];

  for( int iA=0; iA<nfile;iA++){
    for( int ipt=0; ipt<npt;ipt++){
      hDatSimPt[iA][ipt]=(TH1F *) fIn[iA]->Get(Form("hDataSimDiff_ptBins_%d",ipt));
      for(int ita=0;ita<hDatSimPt[iA][ipt]->GetNbinsX();ita++){
	datSim_ratio_per_alphaVal_perpt_perEtaval[iA][ipt][ita]=hDatSimPt[iA][ipt]->GetBinContent(ita+1);
	ddatSim_ratio_per_alphaVal_perpt_perEtaval[iA][ipt][ita]=hDatSimPt[iA][ipt]->GetBinError(ita+1);

      }


    }
    

  }
  //Normalization with alpha=0.3
  for( int ipt=0; ipt<npt;ipt++){
    for(int ita=0;ita<nEta;ita++){
      for( int iA=0; iA<nfile;iA++){
	nAlphaValue[ipt][ita][iA]=datSim_ratio_per_alphaVal_perpt_perEtaval[iA][ipt][ita];
	dnAlphaValue[ipt][ita][iA]=ddatSim_ratio_per_alphaVal_perpt_perEtaval[iA][ipt][ita];
	
	nAlphaValue_ratio_per_alpha030[ipt][ita][iA]=datSim_ratio_per_alphaVal_perpt_perEtaval[iA][ipt][ita]/datSim_ratio_per_alphaVal_perpt_perEtaval[3][ipt][ita];
	dnAlphaValue_ratio_per_alpha030[ipt][ita][iA]=nAlphaValue_ratio_per_alpha030[ipt][ita][iA]*sqrt(pow(ddatSim_ratio_per_alphaVal_perpt_perEtaval[iA][ipt][ita]/datSim_ratio_per_alphaVal_perpt_perEtaval[iA][ipt][ita],2)+pow(ddatSim_ratio_per_alphaVal_perpt_perEtaval[3][ipt][ita]/datSim_ratio_per_alphaVal_perpt_perEtaval[3][ipt][ita],2));

      }
    }
  }

  float dat_sim_ratio_alpha0[npt][nEta];
  float ddat_sim_ratio_alpha0[npt][nEta];

  // Fit the graph with a linear function (polynomial of degree 1)
  TF1* fitFunc[npt][nEta];
    
  TCanvas *cA[npt];
  TCanvas *cAM[npt];
 
#if 1
  /* double avgPtBin[]={50.0,90.0,120.0,160.0,200.0,500.0}; */
  /* double avgPtBin2[]={50.0,90.0,120.0,160.0,200.0,300.0}; */
 double avgPtBin[]={50.0,90.0,120.0,170.0,250.0,1000.0};
  double avgPtBin2[]={50.0,90.0,120.0,170.0,250.0,500.0};

  TGraphErrors* graph_vs_alpha[npt][nEta];
  TGraphErrors* graph_vs_alpha_for_alpha30[npt][nEta];
  const int totBins=5;
  float Xval[]={0.15,0.20,0.25,0.30,0.35};
  float dataSim_ratio_for_same_eta_and_pt_but_diff_alpha[npt][nEta][totBins];
  float ddataSim_ratio_for_same_eta_and_pt_but_diff_alpha[npt][nEta][totBins];
  TH1F *hdata_simRatio_alpha0[npt];
  TH1F *hdata_simRatio_alpha0_newMethod[npt];
  for( int ipt=0; ipt<npt;ipt++){
    cA[ipt]= new TCanvas(Form("cA%d",ipt),Form("cA%d",ipt),10,10,1200,800);
    cA[ipt]->Divide(10,5);
   cAM[ipt]= new TCanvas(Form("cAM%d",ipt),Form("cAM%d",ipt),10,10,1200,800);
    cAM[ipt]->Divide(10,5);

    hdata_simRatio_alpha0[ipt]=(TH1F*)hDatSimPt[0][ipt]->Clone();
    hdata_simRatio_alpha0[ipt]->Reset();
    hdata_simRatio_alpha0[ipt]->SetTitle(Form("hdata_simRatio_for_zero_alpha_oldMethod;#eta;<MC/Data>_{#alpha #rightarrow 0}"));
    hdata_simRatio_alpha0[ipt]->SetName(Form("hdata_simRatio_alpha0_oldMethod_%d",ipt));

    hdata_simRatio_alpha0_newMethod[ipt]=(TH1F*)hdata_simRatio_alpha0[ipt]->Clone();
    hdata_simRatio_alpha0_newMethod[ipt]->Reset();
    hdata_simRatio_alpha0_newMethod[ipt]->SetTitle(Form("hdata_simRatio_for_zero_alpha_newMethod;#eta;<MC/Data>_{#alpha #rightarrow 0}"));
    hdata_simRatio_alpha0_newMethod[ipt]->SetName(Form("hdata_simRatio_alpha0_newMethod_%d",ipt));

    for(int ita=0;ita<nEta;ita++){
      cA[ipt]->cd(ita+1);
      // if(ita==0&&ita==21)continue;
      graph_vs_alpha[ipt][ita] = new TGraphErrors(totBins, Xval, nAlphaValue[ipt][ita],0,dnAlphaValue[ipt][ita]);
      graph_vs_alpha[ipt][ita]->SetTitle(Form("Ratio_vs_#alpha@(pt %d, eta %d);#alpha;MC/Data", ipt, ita));
      //if(ita==0 || ita ==18 ||ita ==19 ||ita==21)
      if(1)
	{
	  graph_vs_alpha[ipt][ita]->GetYaxis()->SetRangeUser(0.75,1.5);
	}
      else graph_vs_alpha[ipt][ita]->GetYaxis()->SetRangeUser(0.9,1.2);
      
      graph_vs_alpha[ipt][ita]->SetMarkerStyle(8);
      graph_vs_alpha[ipt][ita]->SetMarkerColor(4);
      graph_vs_alpha[ipt][ita]->SetLineColor(4);
      graph_vs_alpha[ipt][ita]->Draw("AP");

      fitFunc[ipt][ita]= new TF1(Form("fitFunc_%d_%d",ipt,ita), "pol1", 0.05, 0.35);
      graph_vs_alpha[ipt][ita]->Fit(fitFunc[ipt][ita],"QR");

      cAM[ipt]->cd(ita+1);

      graph_vs_alpha_for_alpha30[ipt][ita] = new TGraphErrors(totBins, Xval, nAlphaValue_ratio_per_alpha030[ipt][ita],0,dnAlphaValue_ratio_per_alpha030[ipt][ita]);
      graph_vs_alpha_for_alpha30[ipt][ita]->SetTitle(Form("Ratio_vs_#alpha@(pt %d, eta %d);#alpha;MC/Data", ipt, ita));
      graph_vs_alpha_for_alpha30[ipt][ita]->SetMarkerStyle(8);
      if(ipt==4)
	{
	  graph_vs_alpha_for_alpha30[ipt][ita]->SetMarkerColor(6);
	  graph_vs_alpha_for_alpha30[ipt][ita]->SetLineColor(6);
	}
      else
	{
	  graph_vs_alpha_for_alpha30[ipt][ita]->SetMarkerColor(ipt+1);
	  graph_vs_alpha_for_alpha30[ipt][ita]->SetLineColor(ipt+1);
	}
      graph_vs_alpha_for_alpha30[ipt][ita]->Draw("AP");
      // Store the fit result
      TFitResultPtr fitResult = graph_vs_alpha[ipt][ita]->Fit(fitFunc[ipt][ita], "QS"); // "QS" for quiet mode, storing fit result

      double x_value = 0;

      //dat_sim_ratio_alpha0[ipt][ita]=fitFunc[ipt][ita]->Eval(0);//Evaluate the x for alpha->0

      dat_sim_ratio_alpha0[ipt][ita]=fitFunc[ipt][ita]->GetParameter(0);//Evaluate the x for alpha->0


      //ddat_sim_ratio_alpha0[ipt][ita]= fitFunc[ipt][ita]->GetParError(0);
      //cout<<ipt<<'\t'<<ita<<'\t'<<dat_sim_ratio_alpha0[ipt][ita]<<'\t'<<ddat_sim_ratio_alpha0[ipt][ita]<<endl;

      // Get the intercept and its error
      double p0 = fitFunc[ipt][ita]->GetParameter(0);      // Intercept
      double p0_error = fitFunc[ipt][ita]->GetParError(0); // Uncertainty of intercept

      // Get the slope and its error
      double p1 = fitFunc[ipt][ita]->GetParameter(1);      // Slope
      double p1_error = fitFunc[ipt][ita]->GetParError(1); // Uncertainty of slope

      // Get the covariance matrix
      //double covariance = fitFunc[ipt][ita]->GetCovarianceMatrixElement(0, 1); // Covariance between p0 and p1
      double covariance = fitResult->CovMatrix(0, 1); // Covariance between p0 and p1

      // Calculate the error at x = 0
      double error_at_x0 = sqrt(p0_error * p0_error + x_value * x_value * p1_error * p1_error + 2 * x_value * covariance);
      ddat_sim_ratio_alpha0[ipt][ita] = sqrt(p0_error * p0_error + x_value * x_value * p1_error * p1_error + 2 * x_value * covariance);
      

    }
    //cA[ipt]->SaveAs(Form("histPlots/MC_Data_Ratio_vs_\u03B1_%.2f_to_%.2f.pdf",avgPtBin[ipt],avgPtBin[ipt+1]),"RECREATE");

      
  }
  
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////

  TMultiGraph* multiGraph_vs_alpha_for_alpha30[nEta];
  TMultiGraph* multiGraph_vs_alpha_for_alpha30_[nEta];
  TF1* fitFunc_nm[nEta];
  float kFSRVal[nEta];
  float c_factor[nEta][npt];
  float dc_factor[nEta][npt];

  TCanvas *cM= new TCanvas(Form("cM"),Form("cM"),10,10,1200,800);
  cM->Divide(6,4);
  TLegend *lg= new TLegend(0.5,0.7,0.88,0.9);//0.7,0.9,0.89,0.99);
  lg->AddEntry(hdata_simRatio_alpha0[0]," 50.0 < p_{T}_{avg} < 90.0 ","pl");
  lg->AddEntry(hdata_simRatio_alpha0[1]," 90.0 < p_{T}_{avg} < 120.0 ","pl");
  lg->AddEntry(hdata_simRatio_alpha0[2]," 120.0 < p_{T}_{avg} < 170.0 ","pl");
  lg->AddEntry(hdata_simRatio_alpha0[3]," 170.0 < p_{T}_{avg} < 250.0 ","pl");
  lg->AddEntry(hdata_simRatio_alpha0[4]," 250.0 < p_{T}_{avg} < 1000.0 ","pl");
  gPad->SetGridx();
  gPad->SetGridy();
  for(int ita=0;ita<nEta;ita++){
    cM->cd(ita+1);
    multiGraph_vs_alpha_for_alpha30[ita]=new TMultiGraph();
    multiGraph_vs_alpha_for_alpha30[ita]->SetTitle(Form("[%.3f < #eta < %.3f];#alpha;MC/Data", etaBin[ita],etaBin[ita+1]));

    for (int i = 0; i < 5; i++) { // Assuming 5 graphs per eta bin
      if (graph_vs_alpha_for_alpha30[i][ita] && graph_vs_alpha_for_alpha30[i][ita]->GetN() > 0) {
        multiGraph_vs_alpha_for_alpha30[ita]->Add(graph_vs_alpha_for_alpha30[i][ita], "AP");
      }
    }
    /* multiGraph_vs_alpha_for_alpha30[ita]->Add(graph_vs_alpha_for_alpha30[0][ita], "AP"); */
    /* multiGraph_vs_alpha_for_alpha30[ita]->Add(graph_vs_alpha_for_alpha30[1][ita], "AP"); */
    /* multiGraph_vs_alpha_for_alpha30[ita]->Add(graph_vs_alpha_for_alpha30[2][ita], "AP"); */
    /* multiGraph_vs_alpha_for_alpha30[ita]->Add(graph_vs_alpha_for_alpha30[3][ita], "AP"); */
    /* multiGraph_vs_alpha_for_alpha30[ita]->Add(graph_vs_alpha_for_alpha30[4][ita], "AP"); */
    multiGraph_vs_alpha_for_alpha30[ita]->Draw("A");
    fitFunc_nm[ita]= new TF1(Form("fitFunc_nm_%d",ita), "pol1", 0.05, 0.35);

    multiGraph_vs_alpha_for_alpha30[ita]->Fit(fitFunc_nm[ita],"FQ");
    //multiGraph_vs_alpha_for_alpha30[ita]->Fit("pol1","FQ");
    cout<<ita<<'\t'<<fitFunc_nm[ita]->GetParameter(0)<<'\t'<<fitFunc_nm[ita]->GetParameter(1)<<endl;
    //if(ita==0 || ita ==1 ||ita ==20 ||ita==21)
     if(1)
     {
	multiGraph_vs_alpha_for_alpha30[ita]->GetYaxis()->SetRangeUser(0.85,1.2);
      }
    else
      {
	multiGraph_vs_alpha_for_alpha30[ita]->GetYaxis()->SetRangeUser(0.95,1.05);
      }
    kFSRVal[ita]=fitFunc_nm[ita]->GetParameter(1);
    //cM->SaveAs(Form("histPlots/multigraph_MC_Data_Ratio_vs_\u03B1_normalized_with_MCData_@_\u03B1_03.pdf"),"RECREATE");

    lg->Draw();
    
  }
  TCanvas *cT= new TCanvas();
  multiGraph_vs_alpha_for_alpha30[17]->Draw("A");
  multiGraph_vs_alpha_for_alpha30[17]->GetYaxis()->SetRangeUser(0.95,1.05);

  for( int ipt=0; ipt<npt;ipt++){
    for(int ita=0;ita<nEta;ita++){
      
      c_factor[ita][ipt]=datSim_ratio_per_alphaVal_perpt_perEtaval[3][ipt][ita]-(kFSRVal[ita]*datSim_ratio_per_alphaVal_perpt_perEtaval[3][ipt][ita]*0.3);
      dc_factor[ita][ipt]=c_factor[ita][ipt]*sqrt(pow(ddatSim_ratio_per_alphaVal_perpt_perEtaval[3][ipt][ita]/datSim_ratio_per_alphaVal_perpt_perEtaval[3][ipt][ita],2));
      //cout<<ita<<'\t'<<kFSRVal[ita]<<'\t'<<"c_factor: "<<c_factor[ita][ipt]<<'\t'<<"dc_factor: "<<dc_factor[ita][ipt]<<endl;
    }
  }


  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////

  TH1F * hAbs_data_simRatio_alpha0[npt];
  TH1F * hAbs_data_simRatio_alpha0_newMethod[npt];
  TH1F * hAbs_old_Vs_newMethod[npt];
  // float AvetaBin[] = {0.0, 0.4, 0.7, 1.0, 1.3, 1.7, 2.0, 2.4, 2.7, 3.0, 3.5, 5.0};
  float AvetaBin[] = {0.0, 0.130, 0.261, 0.391, 0.522, 0.652, 0.783, 0.913, 1.044, 1.174, 1.305, 1.392, 1.479, 1.566, 1.653, 1.791, 1.93, 2.051, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 5.20};
  float nval_ratio_etabin[npt][nEta];
  float dnval_ratio_etabin[npt][nEta];

  float dat_sim_ratio_alpha0_perEta_Vs_pT_newMethod[npt][nEta];
 float  ddat_sim_ratio_alpha0_perEta_Vs_pT_newMethod[npt][nEta];
  float nval_ratio_etabin_newMethod[npt][nEta];
  float dnval_ratio_etabin_newMethod[npt][nEta];
  for( int ipt=0; ipt<npt;ipt++){
    hAbs_data_simRatio_alpha0[ipt] = new TH1F(Form("hAbs_data_simRatio_alpha0_oldMethod%d", ipt), Form("hAbs_data_simRatio_alpha0_oldMethod%d", ipt), nEta / 2, AvetaBin);//0, 5.0);
    hAbs_data_simRatio_alpha0[ipt]->SetTitle("hAbs_data_simRatio_alpha0_oldMethod;|#eta|;<MC/Data>_{#alpha #rightarrow 0}");

    hAbs_data_simRatio_alpha0_newMethod[ipt] = new TH1F(Form("hAbs_data_simRatio_alpha0_newMethod%d", ipt), Form("hAbs_data_simRatio_alpha0_newMethod%d", ipt), nEta / 2, AvetaBin);//0, 5.0);
    hAbs_data_simRatio_alpha0_newMethod[ipt]->SetTitle("hAbs_data_simRatio_alpha0_newMethod;|#eta|;<MC/Data>_{#alpha #rightarrow 0}");

    for (int ita = 0; ita < nEta ; ita++) {

      hdata_simRatio_alpha0[ipt]->SetBinContent(ita+1, dat_sim_ratio_alpha0[ipt][ita]);
      hdata_simRatio_alpha0[ipt]->SetBinError(ita+1, ddat_sim_ratio_alpha0[ipt][ita]);

      // hdata_simRatio_alpha0[ipt]->SetBinError(ita+1, error_at_x0);//ddat_sim_ratio_alpha0[ipt][ita]);


	
      hdata_simRatio_alpha0[ipt]->SetMarkerStyle(8);
      hdata_simRatio_alpha0[ipt]->SetMarkerColor(ipt+1);
      hdata_simRatio_alpha0[ipt]->SetLineColor(ipt+1);
      if(ipt==4)
	{
	  hdata_simRatio_alpha0[ipt]->SetMarkerColor(6);
	  hdata_simRatio_alpha0[ipt]->SetLineColor(6);
	     
	}
      hdata_simRatio_alpha0[ipt]->GetYaxis()->SetRangeUser(0.5,1.5);
      ///////////////////////////////////////////////////
      dat_sim_ratio_alpha0_perEta_Vs_pT_newMethod[ipt][ita]=c_factor[ita][ipt];//values for alpha->0 from extrapolation
      ddat_sim_ratio_alpha0_perEta_Vs_pT_newMethod[ipt][ita]=dc_factor[ita][ipt];

      hdata_simRatio_alpha0_newMethod[ipt]->SetBinContent(ita+1,dat_sim_ratio_alpha0_perEta_Vs_pT_newMethod[ipt][ita] );
      hdata_simRatio_alpha0_newMethod[ipt]->SetBinError(ita+1, ddat_sim_ratio_alpha0_perEta_Vs_pT_newMethod[ipt][ita]);
      hdata_simRatio_alpha0_newMethod[ipt]->SetMarkerStyle(8);
      hdata_simRatio_alpha0_newMethod[ipt]->SetMarkerColor(ipt+1);
      hdata_simRatio_alpha0_newMethod[ipt]->SetLineColor(ipt+1);
      if(ipt==4)
	{
	  hdata_simRatio_alpha0_newMethod[ipt]->SetMarkerColor(6);
	  hdata_simRatio_alpha0_newMethod[ipt]->SetLineColor(6);
	     
	}
      hdata_simRatio_alpha0_newMethod[ipt]->GetYaxis()->SetRangeUser(0.5,1.5);
      
    }
    
    
    for (int ib = 0; ib < nEta / 2; ib++) {
      int posBin = nEta / 2 + ib;  // Positive eta bin
      int negBin = nEta / 2 - ib - 1;  // Negative eta bin
      //if(ib==11)continue;

      nval_ratio_etabin[ipt][posBin]=dat_sim_ratio_alpha0[ipt][posBin];
      nval_ratio_etabin[ipt][negBin]=dat_sim_ratio_alpha0[ipt][negBin];
      dnval_ratio_etabin[ipt][posBin]=ddat_sim_ratio_alpha0[ipt][posBin];
      dnval_ratio_etabin[ipt][negBin]=ddat_sim_ratio_alpha0[ipt][negBin];

      nval_ratio_etabin_newMethod[ipt][posBin]=dat_sim_ratio_alpha0_perEta_Vs_pT_newMethod[ipt][posBin];
      nval_ratio_etabin_newMethod[ipt][negBin]=dat_sim_ratio_alpha0_perEta_Vs_pT_newMethod[ipt][negBin];
      dnval_ratio_etabin_newMethod[ipt][posBin]=ddat_sim_ratio_alpha0_perEta_Vs_pT_newMethod[ipt][posBin];
      dnval_ratio_etabin_newMethod[ipt][negBin]=ddat_sim_ratio_alpha0_perEta_Vs_pT_newMethod[ipt][negBin];
     
      
      float avg_ratio=(nval_ratio_etabin[ipt][posBin]+ nval_ratio_etabin[ipt][negBin])/2.0;
      float davg_ratio=avg_ratio*sqrt(pow(dnval_ratio_etabin[ipt][posBin]/nval_ratio_etabin[ipt][posBin],2)+pow(dnval_ratio_etabin[ipt][negBin]/nval_ratio_etabin[ipt][negBin],2));

      float avg_ratio_newMethod=(nval_ratio_etabin_newMethod[ipt][posBin]+nval_ratio_etabin_newMethod[ipt][negBin])/2.0;
      float davg_ratio_newMethod=avg_ratio_newMethod*sqrt(pow(dnval_ratio_etabin_newMethod[ipt][posBin]/nval_ratio_etabin_newMethod[ipt][posBin],2)+pow(dnval_ratio_etabin_newMethod[ipt][negBin]/nval_ratio_etabin_newMethod[ipt][negBin],2));

	
      hAbs_data_simRatio_alpha0[ipt]->SetBinContent(ib+1,avg_ratio);
      hAbs_data_simRatio_alpha0[ipt]->SetBinError(ib+1,davg_ratio);

      
      hAbs_data_simRatio_alpha0[ipt]->SetMarkerStyle(8);
      hAbs_data_simRatio_alpha0[ipt]->SetMarkerColor(ipt+1);
      hAbs_data_simRatio_alpha0[ipt]->SetLineColor(ipt+1);
      
      if(ipt==4)
	{
	  hAbs_data_simRatio_alpha0[ipt]->SetMarkerColor(6);
	  hAbs_data_simRatio_alpha0[ipt]->SetLineColor(6);
	     
	}
      hAbs_data_simRatio_alpha0[ipt]->GetYaxis()->SetRangeUser(0.5,1.3);
      /////////////////////////////////////////////
      hAbs_data_simRatio_alpha0_newMethod[ipt]->SetBinContent(ib+1,avg_ratio_newMethod);
      hAbs_data_simRatio_alpha0_newMethod[ipt]->SetBinError(ib+1,davg_ratio_newMethod);
      hAbs_data_simRatio_alpha0_newMethod[ipt]->SetMarkerStyle(8);
      hAbs_data_simRatio_alpha0_newMethod[ipt]->SetMarkerColor(ipt+1);
      hAbs_data_simRatio_alpha0_newMethod[ipt]->SetLineColor(ipt+1);
      if(ipt==4)
	{
	  hAbs_data_simRatio_alpha0_newMethod[ipt] ->SetMarkerColor(6);
	  hAbs_data_simRatio_alpha0_newMethod[ipt] ->SetLineColor(6);

	}
      hAbs_data_simRatio_alpha0_newMethod[ipt]->GetYaxis()->SetRangeUser(0.5,1.5);


      hAbs_old_Vs_newMethod[ipt]=(TH1F *)hAbs_data_simRatio_alpha0_newMethod[ipt]->Clone(Form("hAbs_old_Vs_newMethod_%d",ipt));
     hAbs_old_Vs_newMethod[ipt]->SetTitle("with_kFSR_correction_Vs_without_kFSR_correction;|#eta|;with_kFSR/without_kFSR_{#alpha #rightarrow 0}");
     hAbs_old_Vs_newMethod[ipt]->Divide(hAbs_data_simRatio_alpha0[ipt]);
      hAbs_old_Vs_newMethod[ipt]->SetMarkerStyle(8);
      hAbs_old_Vs_newMethod[ipt]->SetMarkerColor(ipt+1);
      hAbs_old_Vs_newMethod[ipt]->SetLineColor(ipt+1);
      if(ipt==4)
	{
	  hAbs_old_Vs_newMethod[ipt] ->SetMarkerColor(6);
	  hAbs_old_Vs_newMethod[ipt] ->SetLineColor(6);

	}
      hAbs_old_Vs_newMethod[ipt]->GetYaxis()->SetRangeUser(0.5,1.3);

      
    }


    
  }


   TLegend *leg= new TLegend(0.35,0.1,0.8,0.35);//0.7,0.9,0.89,0.99);
  leg->AddEntry(hdata_simRatio_alpha0[0]," 50.0 < p_{T}_{avg} < 90.0 ","pl");
  leg->AddEntry(hdata_simRatio_alpha0[1]," 90.0 < p_{T}_{avg} < 120.0 ","pl");
  leg->AddEntry(hdata_simRatio_alpha0[2]," 120.0 < p_{T}_{avg} < 170.0 ","pl");
  leg->AddEntry(hdata_simRatio_alpha0[3]," 170.0 < p_{T}_{avg} < 250.0 ","pl");
  leg->AddEntry(hdata_simRatio_alpha0[4]," 250.0 < p_{T}_{avg} < 1000.0 ","pl");
  TCanvas *cRR=new TCanvas();
  /* cRR->Divide(2,1); */
 /*  cRR->cd(1); */
 /*  hdata_simRatio_alpha0[0]->Draw();//"P"); */
 /*  hdata_simRatio_alpha0[1]->Draw("same");//P"); */
 /*  hdata_simRatio_alpha0[2]->Draw("same");//P"); */
 /*  hdata_simRatio_alpha0[3]->Draw("same");//P"); */
 /*  hdata_simRatio_alpha0[4]->Draw("same");//P"); */
 /*  gPad->SetGridx(); */
 /*  gPad->SetGridy(); */
 /*  leg->Draw(); */
 /* //cRR->SaveAs(Form("histPlots/MC_Data_Ratio_@_\u03B1_0_VS_\u03B7.pdf"),"RECREATE"); */
 /*  cRR->cd(2); */
  
  hdata_simRatio_alpha0_newMethod[0]->Draw();//"P");
  hdata_simRatio_alpha0_newMethod[1]->Draw("same");//P");
  hdata_simRatio_alpha0_newMethod[2]->Draw("same");//P");
  hdata_simRatio_alpha0_newMethod[3]->Draw("same");//P");
  hdata_simRatio_alpha0_newMethod[4]->Draw("same");//P");
  gPad->SetGridx();
  gPad->SetGridy();
  
  leg->Draw();
   cRR->SaveAs(Form("data_simu_diff_rootfiles_alphaCond/with_kFSR_MC_Data_Ratio_@_\u03B1_0_VS_\u03B7.pdf"),"RECREATE");

  
  TCanvas *cR2=new TCanvas();
  
 TLegend *leg2= new TLegend(0.35,0.1,0.8,0.35);//0.7,0.9,0.89,0.99);
  leg2->AddEntry(hdata_simRatio_alpha0[0]," 50.0 < p_{T}_{avg} < 90.0 ","pl");
  leg2->AddEntry(hdata_simRatio_alpha0[1]," 90.0 < p_{T}_{avg} < 120.0 ","pl");
  leg2->AddEntry(hdata_simRatio_alpha0[2]," 120.0 < p_{T}_{avg} < 170.0 ","pl");
  leg2->AddEntry(hdata_simRatio_alpha0[3]," 170.0 < p_{T}_{avg} < 250.0 ","pl");
  leg2->AddEntry(hdata_simRatio_alpha0[4]," 250.0 < p_{T}_{avg} < 1000.0 ","pl");
  /* cR2->Divide(2,1); */
  /* cR2->cd(1); */
  /* gPad->SetGridx(); */
  /* gPad->SetGridy(); */
  
  /* hAbs_data_simRatio_alpha0[0]->Draw(""); */
  /* hAbs_data_simRatio_alpha0[1]->Draw("same"); */
  /* hAbs_data_simRatio_alpha0[2]->Draw("same"); */
  /* hAbs_data_simRatio_alpha0[3]->Draw("same"); */
  /* hAbs_data_simRatio_alpha0[4]->Draw("same"); */
  /* leg2->Draw(); */
  
  /* cR2->cd(2); */
  hAbs_data_simRatio_alpha0_newMethod[0]->Draw("");
  hAbs_data_simRatio_alpha0_newMethod[1]->Draw("same");
  hAbs_data_simRatio_alpha0_newMethod[2]->Draw("same");
  hAbs_data_simRatio_alpha0_newMethod[3]->Draw("same");
  hAbs_data_simRatio_alpha0_newMethod[4]->Draw("same");
  gPad->SetGridx();
  gPad->SetGridy();
  leg2->Draw();
  cR2->SaveAs(Form("data_simu_diff_rootfiles_alphaCond/absoulte_with_kFSR_MC_Data_Ratio_@_\u03B1_0_VS_Absoulte_\u03B7.pdf"),"RECREATE");

  TCanvas *cRA=new TCanvas();
  gPad->SetGridx();
  gPad->SetGridy();
  hAbs_old_Vs_newMethod[0]->Draw();
  hAbs_old_Vs_newMethod[1]->Draw("same");
  hAbs_old_Vs_newMethod[2]->Draw("same");
  hAbs_old_Vs_newMethod[3]->Draw("same");
  hAbs_old_Vs_newMethod[4]->Draw("same");
  leg2->Draw();
  cRA->SaveAs(Form("data_simu_diff_rootfiles_alphaCond/Ratio_absoulte_with_without_kFSR_MC_Data_Ratio_@_\u03B1_0_VS_absloute_\u03B7.pdf"),"RECREATE");

  //cR2->SaveAs(Form("histPlots/MC_Data_Ratio_@_\u03B1_0_VS_Absoulte_\u03B7.pdf"),"RECREATE");

  ////////////////////////////////////////////////////////////////



  

  
  // Plot graph of Ma/Data@ alpha=0 vs Pt for each ita bins//
  float dat_sim_ratio_alpha0_perEta_Vs_pT[nEta][npt];
  float ddat_sim_ratio_alpha0_perEta_Vs_pT[nEta][npt];
  //TGraphErrors* graph_vs_[npt][nEta];
  TH1F *hdata_simRatio_alpha0_Vs_Pt_perEta[nEta];
  TF1 *fit01[nEta];
  TCanvas *cPP= new TCanvas(Form("cPP"),Form("cPP"),10,10,1200,800);
  cPP->Divide(10,5);
  
    for(int ita=0;ita<nEta;ita++){
      hdata_simRatio_alpha0_Vs_Pt_perEta[ita]= new TH1F(Form("hdata_simRatio_alpha0_Vs_Pt_perEta_%.3f_%.3f",etaBin[ita],etaBin[ita+1]),Form("hdata_simRatio_alpha0_Vs_Pt_perEta_%.1f_%.1f",etaBin[ita],etaBin[ita+1]),npt,avgPtBin2);
        for( int ipt=0; ipt<npt;ipt++){
	  dat_sim_ratio_alpha0_perEta_Vs_pT[ita][ipt]=dat_sim_ratio_alpha0[ipt][ita];//values for alpha->0 from extrapolation
	  ddat_sim_ratio_alpha0_perEta_Vs_pT[ita][ipt]=ddat_sim_ratio_alpha0[ipt][ita];
	  hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->SetBinContent(ipt+1,dat_sim_ratio_alpha0_perEta_Vs_pT[ita][ipt]);
	  hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->SetBinError(ipt+1,ddat_sim_ratio_alpha0_perEta_Vs_pT[ita][ipt]);
	  hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->SetMarkerStyle(8);
	  hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->SetMarkerColor(4);
	  hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->SetLineColor(4);
	  
	    
	  hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->SetTitle(Form("[%.3f < #eta < %.3f];p_{T}^{avg};<MC/Data>_{#alpha #rightarrow 0}", etaBin[ita],etaBin[ita+1]));
	}
    }
    //Without KFSR Method//
    for(int ita=0;ita<nEta;ita++){
      cPP->cd(ita+1);
      
      fit01[ita]=new TF1(Form("fit01%d",ita),fitfunct,50,1000,3);
      fit01[ita]->SetLineColor(2);
      fit01[ita]->SetParameter(0,1.50);
      fit01[ita]->SetParameter(1,1.50);
      fit01[ita]->SetParameter(2,1.50);
      hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->Fit(fit01[ita],"QR");
      cPP->SaveAs(Form("data_simu_diff_rootfiles_alphaCond/MC_Data_Ratio_@_\u03B1_0_VS_pT_per_\u03B7.pdf"),"RECREATE");

      //if(ita==0|| ita==21){
      if(1){
	hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->GetYaxis()->SetRangeUser(0.75,1.5);
      }
      else
	{
	  hdata_simRatio_alpha0_Vs_Pt_perEta[ita]->GetYaxis()->SetRangeUser(0.75,1.25);
	}
     
    }

	
    // outFile.close();
	
    //###############################################################################//
    //###############################################################################//
    //###############################################################################//
    //###############################################################################//


 // Plot graph of Ma/Data@ alpha=0 vs Pt for each ita bins//
  float dat_sim_ratio_alpha0_perEta_Vs_pT_cfact[nEta][npt];
  float ddat_sim_ratio_alpha0_perEta_Vs_pT_cfact[nEta][npt];
  //TGraphErrors* graph_vs_[npt][nEta];
  TH1F *hdata_simRatio_alpha0_Vs_Pt_perEta_cfact[nEta];
  TF1 *fit01_cfact[nEta];
  TCanvas *cPPf= new TCanvas(Form("cPPf"),Form("cPPf"),10,10,1200,800);
  cPPf->Divide(10,5);
  
    for(int ita=0;ita<nEta;ita++){
      hdata_simRatio_alpha0_Vs_Pt_perEta_cfact[ita]= new TH1F(Form("hdata_simRatio_alpha0_Vs_Pt_perEta_cfact_%.3f-%.3f",etaBin[ita],etaBin[ita+1]),Form("hdata_simRatio_alpha0_Vs_Pt_perEta_cfact_%.1f_%.1f",etaBin[ita],etaBin[ita+1]),npt,avgPtBin2);
        for( int ipt=0; ipt<npt;ipt++){
	  dat_sim_ratio_alpha0_perEta_Vs_pT_cfact[ita][ipt]=c_factor[ita][ipt];//values for alpha->0 from extrapolation
	  ddat_sim_ratio_alpha0_perEta_Vs_pT_cfact[ita][ipt]=dc_factor[ita][ipt];
	  hdata_simRatio_alpha0_Vs_Pt_perEta_cfact[ita]->SetBinContent(ipt+1,dat_sim_ratio_alpha0_perEta_Vs_pT_cfact[ita][ipt]);
	  hdata_simRatio_alpha0_Vs_Pt_perEta_cfact[ita]->SetBinError(ipt+1,ddat_sim_ratio_alpha0_perEta_Vs_pT_cfact[ita][ipt]);
	  hdata_simRatio_alpha0_Vs_Pt_perEta_cfact[ita]->SetMarkerStyle(8);
	  hdata_simRatio_alpha0_Vs_Pt_perEta_cfact[ita]->SetMarkerColor(1);
	  hdata_simRatio_alpha0_Vs_Pt_perEta_cfact[ita]->SetLineColor(1);
	  
	    
	  hdata_simRatio_alpha0_Vs_Pt_perEta_cfact[ita]->SetTitle(Form("[%.3f < #eta < %.3f];p_{T}^{avg};<MC/Data>_{#alpha #rightarrow 0}", etaBin[ita],etaBin[ita+1]));
	}
    }
     //With KFSR Method//
   for(int ita=0;ita<nEta;ita++){
      cPPf->cd(ita+1);
      
      fit01_cfact[ita]=new TF1(Form("fit01_cfact%d",ita),fitfunct,50,1000,3);
      fit01_cfact[ita]->SetLineColor(2);
      fit01_cfact[ita]->SetParameter(0,1.50);
      fit01_cfact[ita]->SetParameter(1,1.50);
      fit01_cfact[ita]->SetParameter(2,1.50);
      hdata_simRatio_alpha0_Vs_Pt_perEta_cfact[ita]->Fit(fit01_cfact[ita],"QR");
      

      //if(ita==0|| ita==21){
      if(1){
	hdata_simRatio_alpha0_Vs_Pt_perEta_cfact[ita]->GetYaxis()->SetRangeUser(0.75,1.5);
      }
      else
	{
	  hdata_simRatio_alpha0_Vs_Pt_perEta_cfact[ita]->GetYaxis()->SetRangeUser(0.75,1.5);
	}
  
      cPPf->SaveAs(Form("data_simu_diff_rootfiles_alphaCond/MC_Data_Ratio_@_\u03B1_0_VS_pT_per_\u03B7_cFact_newMethod.pdf"),"RECREATE");
//cout<< etaBin2[ita] <<'\t'<< etaBin2[ita+1] <<'\t'<< 5 <<'\t'<< 50 <<'\t'<< 500 <<'\t'<< fit01[ita]->GetParameter(0) <<'\t'<< fit01[ita]->GetParameter(1) <<'\t'<< fit01[ita]->GetParameter(2) <<endl;
      outFile<< etaBin2[ita] <<'\t'
	     << etaBin2[ita+1] <<'\t'
	     << 5 <<'\t'<< 50 <<'\t'
	     << 1000 <<'\t'
	     << fit01_cfact[ita]->GetParameter(0) <<'\t'
	     << fit01_cfact[ita]->GetParameter(1) <<'\t'
	     << fit01_cfact[ita]->GetParameter(2) <<'\n';
    }
   TCanvas *cPPf1= new TCanvas();
   hdata_simRatio_alpha0_Vs_Pt_perEta_cfact[17]->Draw();
 
	
    outFile.close();
	

    
    
    
#endif




}



