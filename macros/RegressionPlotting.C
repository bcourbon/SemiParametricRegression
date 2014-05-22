#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveLabel.h"
#include "TLegend.h"

#include <iostream>
#include <fstream>

void RegressionPlotting(bool doEBEE=true, bool doEB=false, bool doEE=false){

  TString name;
  if (doEBEE) name="resultsEBEE.root";
  else if (doEB) name="resultsEB.root";
  else if (doEE) name="resultsEE.root";

  TFile *fResults = TFile::Open(name);



  //Plot eraw/etrue and ecor/etrue to quantify regression performance
  
  TH1 *hecor = (TH1*)fResults->Get("hCor_true");
  TH1 *heraw = (TH1*)fResults->Get("hRaw_true");
  hecor->SetLineColor(kRed);
  heraw->SetLineColor(kBlue);
  hecor->GetXaxis()->SetRangeUser(0.8,1.2);  

  TCanvas *cresponse = new TCanvas;
  hecor->Draw("HIST");
  heraw->Draw("HISTSAME");
  hecor->SetLineWidth(2);
  heraw->SetLineWidth(2);
  TLegend *leg = new TLegend(0.12,0.78,0.32,0.88);
  leg->AddEntry(heraw,"Eraw/Etrue","l");
  leg->AddEntry(hecor,"Ecor/Etrue","l");
  leg->Draw();
  hecor->SetTitle("E/Etrue");
  hecor->GetXaxis()->SetTitle("E/Etrue");
  gStyle->SetOptStat(0);
  TPaveLabel *OldPerf = new TPaveLabel(0.63,0.88,0.88,0.78, Form("Peak = %.3f , RMS = %.3f ",heraw->GetXaxis()->GetBinCenter(heraw->GetMaximumBin()),heraw->GetRMS()),"NDC");
  TPaveLabel *NewPerf = new TPaveLabel(0.63,0.76,0.88,0.66, Form("Peak = %.3f , RMS = %.3f ",hecor->GetXaxis()->GetBinCenter(hecor->GetMaximumBin()),hecor->GetRMS()),"NDC");
  OldPerf->SetBorderSize(1);
  OldPerf->SetTextSize(0.25);
  OldPerf->SetTextColor(kBlue);
  NewPerf->SetBorderSize(1);
  NewPerf->SetTextSize(0.25);
  NewPerf->SetTextColor(kRed);
  OldPerf->Draw();
  NewPerf->Draw();
  cresponse->SaveAs("plots/Performance.png");
  


  //Plot ecor/eraw
  TH1 *hecor_raw = (TH1*)fResults->Get("hCor_raw");
  TCanvas *canCor=new TCanvas("c2","",600,600);
  hecor_raw->Draw();
  hecor_raw->GetXaxis()->SetTitle("Ecor/Eraw");
  hecor_raw->SetTitle("Ecor/Eraw");  
  canCor->SaveAs("plots/Correction.png");



  // Plot Correction vs input variables

  TCanvas *can=new TCanvas("c3","",1300,700);
  can->Divide(3,2);

  can->cd(1);
  TProfile *profCorEta=(TProfile*)fResults->Get("profEta");
  profCorEta->Draw();
  profCorEta->GetXaxis()->SetTitle("Eta");
  profCorEta->GetYaxis()->SetTitle("Ecor/Eraw");
  profCorEta->SetTitle("Profile of Ecor/Eraw vs Eta");
  profCorEta->SetMarkerStyle(3);

  can->cd(2);
  TProfile *profCorphi= (TProfile*)fResults->Get("profPhi");
  profCorphi->Draw();
  profCorphi->GetXaxis()->SetTitle("Phi");
  profCorphi->GetYaxis()->SetTitle("Ecor/Eraw");
  profCorphi->SetTitle("Profile of Ecor/Eraw vs Phi");
  profCorphi->SetMarkerStyle(3);

  can->cd(3);
  TProfile *profCoretawidth=(TProfile*)fResults->Get("profEtaWidth");
  profCoretawidth->Draw();
  profCoretawidth->GetXaxis()->SetTitle("Eta Width");
  profCoretawidth->GetYaxis()->SetTitle("Ecor/Eraw");
  profCoretawidth->SetTitle("Profile of Ecor/Eraw vs Eta Width");
  profCoretawidth->SetMarkerStyle(3);

  can->cd(4);
  TProfile *profCorphiwidth=(TProfile*)fResults->Get("profPhiWidth");
  profCorphiwidth->Draw();
  profCorphiwidth->GetXaxis()->SetTitle("Phi Width");
  profCorphiwidth->GetYaxis()->SetTitle("Ecor/Eraw");
  profCorphiwidth->SetTitle("Profile of Ecor/Eraw vs Phi Width");
  profCorphiwidth->SetMarkerStyle(3);

  can->cd(5);
  TProfile *profCorr9= (TProfile*)fResults->Get("profR9");
  profCorr9->Draw();
  profCorr9->GetXaxis()->SetTitle("R9");
  profCorr9->GetYaxis()->SetTitle("Ecor/Eraw");
  profCorr9->SetTitle("Profile of Ecor/Eraw vs R9");
  profCorr9->SetMarkerStyle(3);

  can->cd(6);
  TProfile *profCornVtx=(TProfile*)fResults->Get("profNvtx");
  profCornVtx->Draw();
  profCornVtx->GetXaxis()->SetTitle("nVtx");
  profCornVtx->GetYaxis()->SetTitle("Ecor/Eraw");
  profCornVtx->SetTitle("Profile of Ecor/Eraw vs nVtx");
  profCornVtx->SetMarkerStyle(3);



  can->SaveAs("plots/Correction_Inputs.png"); 


  //Plot Bias and Resolution vs Pt

  TCanvas *canPt=new TCanvas("c4","",1100,500);

  canPt->Divide(2,0);

  canPt->cd(1);
  TProfile *profCorPt=(TProfile*)fResults->Get("profPtCor");
  profCorPt->Draw();
  profCorPt->GetXaxis()->SetTitle("Pt");
  profCorPt->GetYaxis()->SetTitle("Ecor/Etrue");
  profCorPt->SetTitle("Profile of Ecor/Etrue vs Pt");
  profCorPt->SetMarkerStyle(3);

  canPt->cd(2);
  TProfile *profResPt=(TProfile*)fResults->Get("profPtRes");
  profResPt->Draw();
  profResPt->GetXaxis()->SetTitle("Pt");
  profResPt->GetYaxis()->SetTitle("|Ecor-Etrue|/Eraw");
  profResPt->SetTitle("Profile of SigmaE/E vs Pt");
  profResPt->SetMarkerStyle(3);

  canPt->SaveAs("plots/Correction_Pt.png"); 


}

