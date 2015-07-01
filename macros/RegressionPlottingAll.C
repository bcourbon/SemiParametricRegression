#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TText.h"
#include "RooDoubleCBFast.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TColor.h"

#include <iostream>
#include <fstream>

using namespace RooFit ;

//Perform a fit with a double CB and keep mean and mean error values

RooDoubleCBFast* FitWithDCB(RooDataHist *Hist, RooRealVar *x, double &mean, double &meanerr){

  RooRealVar *mu=new RooRealVar("mu","mu",1,0.5,1.5);
  RooRealVar *sig=new RooRealVar("sig","sig",0.05,0.001,1);
  RooRealVar *alpha1=new RooRealVar("alpha1","alpha1",2.,0,10);
  RooRealVar *alpha2=new RooRealVar("alpha2","alpha2",1.,0,10);
  RooRealVar *n1=new RooRealVar("n1","n2",2,0,20);
  RooRealVar *n2=new RooRealVar("n2","n2",2,0,20);
  RooDoubleCBFast *DoubleCB=new RooDoubleCBFast("doubleCB","doubleCB",*x,*mu,*sig,*alpha1,*n1,*alpha2,*n2);
  DoubleCB->fitTo(*Hist);
  mean=mu->getValV();
  meanerr=mu->getError();
  return DoubleCB;

}

//effsigma function from Chris
double effSigma(TH1 * hist)
{

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
    return 0.;
  }
  
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }

  Int_t ierr=0;
  Int_t ismin=999;
  
  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }   
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;
  
  return widmin;
  
}



void RegressionPlottingAll(bool doEBEE=true, bool doEB=false, bool doEE=false){

  TString name1;
  TString name2;
  TString name3;
  string dir;

  if (doEBEE){ 
  name1="results_noPU_new.root";
  name2="results_bx50_new.root";
  name3="results_bx25_new.root";
  dir="plots";
  }
  if (doEB){ 
 name1="resultsEB_noPU_new.root";
  name2="resultsEB_bx50_new.root";
  name3="resultsEB_bx25_new.root";
  dir="plotsEB";
  }
  if (doEE){ 
 name1="resultsEE_noPU_new.root";
  name2="resultsEE_bx50_new.root";
  name3="resultsEE_bx25_new.root";
  dir="plotsEE";
  }

  TFile *fResults1 = TFile::Open(name1);
  TFile *fResults2 = TFile::Open(name2);
  TFile *fResults3 = TFile::Open(name3);

  //Plot eraw/etrue and ecor/etrue to quantify regression performance
  
  //Get Histograms
  TH1 *hecor1 = (TH1*)fResults1->Get("hCor_true");
  TH1 *heraw1 = (TH1*)fResults1->Get("hRaw_true");
  TH1 *hecor2 = (TH1*)fResults2->Get("hCor_true");
  TH1 *heraw2 = (TH1*)fResults2->Get("hRaw_true");
  TH1 *hecor3 = (TH1*)fResults3->Get("hCor_true");
  TH1 *heraw3 = (TH1*)fResults3->Get("hRaw_true");
  hecor1->Scale(1./hecor1->Integral());
  heraw1->Scale(1./heraw1->Integral());
  hecor2->Scale(1./hecor2->Integral());
  heraw2->Scale(1./heraw2->Integral());
  hecor3->Scale(1./hecor3->Integral());
  heraw3->Scale(1./heraw3->Integral());
  //Compute peak of the distributions by fitting histograms with double Crystal Ball
  double PeakRaw1,PeakRawErr1,PeakCor1,PeakCorErr1;
  double PeakRaw2,PeakRawErr2,PeakCor2,PeakCorErr2;
  double PeakRaw3,PeakRawErr3,PeakCor3,PeakCorErr3;
  RooRealVar *bias=new RooRealVar("Ecor/Etrue","Ecor/Etrue",0.8,1.2);
  const RooArgList *var=new RooArgList(*bias,"");
  RooPlot *plot = bias->frame(Range(0.9,1.07));
  plot->SetTitle("");
  plot->GetYaxis()->SetTitle("");
  RooDataHist *Hraw1=new RooDataHist("","",*var,heraw1);
  RooDoubleCBFast *DoubleCBraw1=FitWithDCB(Hraw1,bias,PeakRaw1,PeakRawErr1); 
  Hraw1->plotOn(plot,MarkerSize(0.8),MarkerColor(kBlue+2),MarkerStyle(24));
  DoubleCBraw1->plotOn(plot,LineColor(kBlue+2),LineStyle(2),LineWidth(2),Name("hraw1"));
  RooDataHist *Hcor1=new RooDataHist("","",*var,hecor1);
  RooDoubleCBFast *DoubleCBcor1=FitWithDCB(Hcor1,bias,PeakCor1,PeakCorErr1); 
  Hcor1->plotOn(plot,MarkerSize(0.6),MarkerColor(kBlue+2));
  DoubleCBcor1->plotOn(plot,LineColor(kBlue+2),LineWidth(2),Name("hcor1"));
  RooDataHist *Hraw2=new RooDataHist("","",*var,heraw2);
  RooDoubleCBFast *DoubleCBraw2=FitWithDCB(Hraw2,bias,PeakRaw2,PeakRawErr2); 
  Hraw2->plotOn(plot,MarkerSize(0.8),MarkerColor(kRed+1),MarkerStyle(24));
  DoubleCBraw2->plotOn(plot,LineColor(kRed+1),LineStyle(2),LineWidth(2),Name("hraw2"));
  RooDataHist *Hcor2=new RooDataHist("","",*var,hecor2);
  RooDoubleCBFast *DoubleCBcor2=FitWithDCB(Hcor2,bias,PeakCor2,PeakCorErr2); 
  Hcor2->plotOn(plot,MarkerSize(0.6),MarkerColor(kRed+1));
  DoubleCBcor2->plotOn(plot,LineColor(kRed+1),LineWidth(2),Name("hcor2"));
  RooDataHist *Hraw3=new RooDataHist("","",*var,heraw3);
  RooDoubleCBFast *DoubleCBraw3=FitWithDCB(Hraw3,bias,PeakRaw3,PeakRawErr3); 
  Hraw3->plotOn(plot,MarkerSize(0.8),MarkerColor(kOrange+1),MarkerStyle(24));
  DoubleCBraw3->plotOn(plot,LineColor(kOrange+1),LineStyle(2),LineWidth(2),Name("hraw3"));
  RooDataHist *Hcor3=new RooDataHist("","",*var,hecor3);
  RooDoubleCBFast *DoubleCBcor3=FitWithDCB(Hcor3,bias,PeakCor3,PeakCorErr3); 
  Hcor3->plotOn(plot,MarkerSize(0.6),MarkerColor(kOrange+1));
  DoubleCBcor3->plotOn(plot,LineColor(kOrange+1),LineWidth(2),Name("hcor3"));
  //Compute Effective Sigma of the distributions
  /*double SigEffRaw1=effSigma(heraw1);
  double SigEffCor1=effSigma(hecor1);
  double SigEffRaw2=effSigma(heraw2);
  double SigEffCor2=effSigma(hecor2);
  double SigEffRaw3=effSigma(heraw3);
  double SigEffCor3=effSigma(hecor3);
  double SigEffErr=0.001;*/
  //Legend and StatBox
  TLegend *leg = new TLegend(0.11,0.39,0.4,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(plot->findObject("hraw1"),"Eraw/Etrue, no PU","l");
  leg->AddEntry(plot->findObject("hcor1"),"Ecor/Etrue, no PU","l");
  leg->AddEntry(plot->findObject("hraw2"),"Eraw/Etrue, bx 50","l");
  leg->AddEntry(plot->findObject("hcor2"),"Ecor/Etrue, bx 50","l");
  leg->AddEntry(plot->findObject("hraw3"),"Eraw/Etrue, bx 25","l");
  leg->AddEntry(plot->findObject("hcor3"),"Ecor/Etrue, bx 25","l");
  //leg->SetTextFont(62);
  //leg->SetTextSize(0.06);
  plot->addObject(leg);
  /*TPaveText *stats = new TPaveText(0.58,0.39,0.89,0.89,"NDC");
  stats->SetBorderSize(0);
  stats->SetFillColor(0);
  stats->AddText(Form("Peak Raw = %.4f +/- %.0e",PeakRaw,PeakRawErr));
  stats->AddText(Form("Peak Cor = %.4f +/- %.0e",PeakCor,PeakCorErr));
  stats->AddText(Form("Sigma Eff Raw = %.3f +/- %.3f",SigEffRaw,SigEffErr));
  stats->AddText(Form("Sigma Eff Cor = %.3f +/- %.3f",SigEffCor,SigEffErr));
  stats->SetTextSize(0.035);
  plot->addObject(stats);*/
  //Plot
  TCanvas *cresponse=new TCanvas;
  plot->Draw();
  cresponse->SaveAs(Form("plotsAll/%s/Performance.pdf",dir.c_str()));
  cresponse->SaveAs(Form("plotsAll/%s/Performance.png",dir.c_str()));
  

  //Plot ecor/eraw
/*  TH1 *hecor_raw = (TH1*)fResults->Get("hCor_raw");
  TCanvas *canCor=new TCanvas;
  hecor_raw->Draw("HIST");
  hecor_raw->SetLineWidth(2);
  hecor_raw->GetXaxis()->SetTitle("Ecor/Eraw");
  hecor_raw->SetTitle("Ecor/Eraw");  
  canCor->SaveAs(Form("plotsAll/%s/Correction.png",dir.c_str()));
*/

  //Plot |ecor-etrue|/eraw and |eraw-etrue|/eraw
  TH1 *hres1 = (TH1*)fResults1->Get("hRes");
  TH1 *hres2 = (TH1*)fResults2->Get("hRes");
  TH1 *hres3 = (TH1*)fResults3->Get("hRes");
  hres1->Scale(1./hres1->Integral());
  hres2->Scale(1./hres2->Integral());
  hres3->Scale(1./hres3->Integral());
  TCanvas *canR=new TCanvas;
  canR->SetLogy();
  hres1->Draw("HIST,SAME");
  hres1->SetLineWidth(2);
  hres1->SetLineColor(kBlue+2);
  hres1->GetXaxis()->SetRangeUser(0,0.05);
  hres2->Draw("HIST,SAME");
  hres2->SetLineWidth(2);
  hres2->SetLineColor(kRed+1);
  hres2->GetXaxis()->SetRangeUser(0,0.05);
  hres3->Draw("HIST,SAME");
  hres3->SetLineWidth(2);
  hres3->SetLineColor(kOrange+1);
  hres3->GetXaxis()->SetRangeUser(0,0.05);
  hres1->GetXaxis()->SetTitle("#sigma(E)/E");
  hres1->SetTitle("");
  hres1->GetYaxis()->SetTitle("");
  hres1->SetStats(0);
TLegend *leg2 = new TLegend(0.59,0.11,0.89,0.41);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->AddEntry(hres1,"#sigma(E)/E, no PU","l");
  leg2->AddEntry(hres2,"#sigma(E)/E, bx 50","l");
  leg2->AddEntry(hres3,"#sigma(E)/E, bx 25","l");  
  leg2->Draw(); 
 canR->SaveAs(Form("plotsAll/%s/Resolution.pdf",dir.c_str()));
 canR->SaveAs(Form("plotsAll/%s/Resolution.png",dir.c_str()));


  // Plot Correction vs input variables

/*
  TCanvas *canCoreta=new TCanvas;
  gStyle->SetOptStat(0);
  TProfile *profCorEta=(TProfile*)fResults->Get("profEta");
  TProfile *profRawEta=(TProfile*)fResults->Get("profEta0");
  profRawEta->Draw();
  profCorEta->Draw("same");
  profRawEta->GetXaxis()->SetTitle("Eta");
  profRawEta->GetYaxis()->SetTitle("E/Eraw");
  profRawEta->SetTitle("Profile of E/Eraw vs Eta");
  profCorEta->SetMarkerStyle(8);
  profCorEta->SetLineColor(kRed);
  profCorEta->SetMarkerColor(kRed);
  profRawEta->SetMarkerStyle(4);
  profRawEta->SetLineColor(kBlue);
  profRawEta->SetMarkerColor(kBlue);
  TLegend *leg2 = new TLegend(0.68,0.12,0.88,0.27);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->AddEntry(profRawEta,"Eraw/Etrue","lp");
  leg2->AddEntry(profCorEta,"Ecor/Etrue","lp");
  leg2->SetTextFont(62);
  //leg2->SetTextSize(0.06);
  leg2->Draw();
  canCoreta->SaveAs(Form("plotsAll/%s/CorEta.png",dir.c_str()));


  TCanvas *canCorPhi=new TCanvas;
  gStyle->SetOptStat(0);
  TProfile *profCorPhi=(TProfile*)fResults->Get("profPhi");
  TProfile *profRawPhi=(TProfile*)fResults->Get("profPhi0");
  profRawPhi->Draw();
  profCorPhi->Draw("same");
  profRawPhi->GetXaxis()->SetTitle("Phi");
  profRawPhi->GetYaxis()->SetTitle("E/Eraw");
  profRawPhi->SetTitle("Profile of E/Eraw vs Phi");
  profCorPhi->SetMarkerStyle(8);
  profCorPhi->SetLineColor(kRed);
  profCorPhi->SetMarkerColor(kRed);
  profRawPhi->SetMarkerStyle(4);
  profRawPhi->SetLineColor(kBlue);
  profRawPhi->SetMarkerColor(kBlue);
  leg2->Draw();
  canCorPhi->SaveAs(Form("plotsAll/%s/CorPhi.png",dir.c_str()));

  TCanvas *canCorEtaWidth=new TCanvas;
  gStyle->SetOptStat(0);
  TProfile *profCorEtaWidth=(TProfile*)fResults->Get("profEtaWidth");
  TProfile *profRawEtaWidth=(TProfile*)fResults->Get("profEtaWidth0");
  profRawEtaWidth->Draw();
  profCorEtaWidth->Draw("same");
  profRawEtaWidth->GetXaxis()->SetTitle("EtaWidth");
  profRawEtaWidth->GetYaxis()->SetTitle("E/Eraw");
  profRawEtaWidth->SetTitle("Profile of E/Eraw vs EtaWidth");
  profCorEtaWidth->SetMarkerStyle(8);
  profCorEtaWidth->SetLineColor(kRed);
  profCorEtaWidth->SetMarkerColor(kRed);
  profRawEtaWidth->SetMarkerStyle(4);
  profRawEtaWidth->SetLineColor(kBlue);
  profRawEtaWidth->SetMarkerColor(kBlue);
  leg2->Draw();
  canCorEtaWidth->SaveAs(Form("plotsAll/%s/CorEtaWidth.png",dir.c_str()));

  TCanvas *canCorPhiWidth=new TCanvas;
  gStyle->SetOptStat(0);
  TProfile *profCorPhiWidth=(TProfile*)fResults->Get("profPhiWidth");
  TProfile *profRawPhiWidth=(TProfile*)fResults->Get("profPhiWidth0");
  profRawPhiWidth->Draw();
  profCorPhiWidth->Draw("same");
  profRawPhiWidth->GetXaxis()->SetTitle("PhiWidth");
  profRawPhiWidth->GetYaxis()->SetTitle("E/Eraw");
  profRawPhiWidth->SetTitle("Profile of E/Eraw vs PhiWidth");
  profCorPhiWidth->SetMarkerStyle(8);
  profCorPhiWidth->SetLineColor(kRed);
  profCorPhiWidth->SetMarkerColor(kRed);
  profRawPhiWidth->SetMarkerStyle(4);
  profRawPhiWidth->SetLineColor(kBlue);
  profRawPhiWidth->SetMarkerColor(kBlue);
  leg2->Draw();
  canCorPhiWidth->SaveAs(Form("plotsAll/%s/CorPhiWidth.png",dir.c_str()));

  TCanvas *canCorR9=new TCanvas;
  gStyle->SetOptStat(0);
  TProfile *profCorR9=(TProfile*)fResults->Get("profR9");
  TProfile *profRawR9=(TProfile*)fResults->Get("profR90");
  profRawR9->Draw();
  profCorR9->Draw("same");
  profRawR9->GetXaxis()->SetTitle("R9");
  profRawR9->GetYaxis()->SetTitle("E/Eraw");
  profRawR9->SetTitle("Profile of E/Eraw vs R9");
  profCorR9->SetMarkerStyle(8);
  profCorR9->SetLineColor(kRed);
  profCorR9->SetMarkerColor(kRed);
  profRawR9->SetMarkerStyle(4);
  profRawR9->SetLineColor(kBlue);
  profRawR9->SetMarkerColor(kBlue);
  leg2->Draw();
  canCorR9->SaveAs(Form("plotsAll/%s/CorR9.png",dir.c_str()));

  TCanvas *canCorNvtx=new TCanvas;
  gStyle->SetOptStat(0);
  TProfile *profCorNvtx=(TProfile*)fResults->Get("profNvtx");
  TProfile *profRawNvtx=(TProfile*)fResults->Get("profNvtx0");
  profRawNvtx->Draw();
  profCorNvtx->Draw("same");
  profRawNvtx->GetXaxis()->SetTitle("Nvtx");
  profRawNvtx->GetYaxis()->SetTitle("E/Eraw");
  profRawNvtx->SetTitle("Profile of E/Eraw vs Nvtx");
  profCorNvtx->SetMarkerStyle(8);
  profCorNvtx->SetLineColor(kRed);
  profCorNvtx->SetMarkerColor(kRed);
  profRawNvtx->SetMarkerStyle(4);
  profRawNvtx->SetLineColor(kBlue);
  profRawNvtx->SetMarkerColor(kBlue);
  leg2->Draw();
  canCorNvtx->SaveAs(Form("plotsAll/%s/CorNvtx.png",dir.c_str()));

  TCanvas *canCorPt=new TCanvas;
  gStyle->SetOptStat(0);
  TProfile *profCorPt=(TProfile*)fResults->Get("profPt");
  TProfile *profRawPt=(TProfile*)fResults->Get("profPt0");
  profRawPt->Draw();
  profCorPt->Draw("same");
  profRawPt->GetXaxis()->SetTitle("Pt");
  profRawPt->GetYaxis()->SetTitle("E/Eraw");
  profRawPt->SetTitle("Profile of E/Eraw vs Pt");
  profCorPt->SetMarkerStyle(8);
  profCorPt->SetLineColor(kRed);
  profCorPt->SetMarkerColor(kRed);
  profRawPt->SetMarkerStyle(4);
  profRawPt->SetLineColor(kBlue);
  profRawPt->SetMarkerColor(kBlue);
  leg2->Draw();
  canCorPt->SaveAs(Form("plotsAll/%s/CorPt.png",dir.c_str()));


*/

  //Plot Bias and Resolution vs Pt

 /* TCanvas *canResPt=new TCanvas;
  TProfile *profResPt=(TProfile*)fResults->Get("profPtRes");
  profResPt->Draw();
  profResPt->GetXaxis()->SetTitle("Pt");
  profResPt->GetYaxis()->SetTitle("|Ecor-Etrue|/Eraw");
  profResPt->SetTitle("Profile of SigmaE/E vs Pt");
  profResPt->SetMarkerStyle(3);
  profResPt->SetLineColor(kRed);
  canResPt->SaveAs(Form("plotsAll/%s/ResPt.png",dir.c_str()));
*/

//Plot Bias and resolution in bins of Pt and eta


  double errpt[9];
  double pt[9];
  double sigeff11[9];
  double errsigeff11[9];
  double sigeff12[9];
  double errsigeff12[9];
  double sigeff13[9];
  double errsigeff13[9];
  double sigeff14[9];
  double errsigeff14[9];
  double mean11[9];
  double errmean11[9];
  double mean12[9];
  double errmean12[9];
  double mean13[9];
  double errmean13[9];
  double mean14[9];
  double errmean14[9];
  double sigeff21[9];
  double errsigeff21[9];
  double sigeff22[9];
  double errsigeff22[9];
  double sigeff23[9];
  double errsigeff23[9];
  double sigeff24[9];
  double errsigeff24[9];
  double mean21[9];
  double errmean21[9];
  double mean22[9];
  double errmean22[9];
  double mean23[9];
  double errmean23[9];
  double mean24[9];
  double errmean24[9];
  double sigeff31[9];
  double errsigeff31[9];
  double sigeff32[9];
  double errsigeff32[9];
  double sigeff33[9];
  double errsigeff33[9];
  double sigeff34[9];
  double errsigeff34[9];
  double mean31[9];
  double errmean31[9];
  double mean32[9];
  double errmean32[9];
  double mean33[9];
  double errmean33[9];
  double mean34[9];
  double errmean34[9];

/*  TCanvas *can10 = new TCanvas("can7","can7",1200,700);
  can10->Divide(5,4);
  TCanvas *can11 = new TCanvas("can8","can8",1200,700);
  can11->Divide(5,4);
  TCanvas *can12 = new TCanvas("can9","can9",1200,700);
  can12->Divide(5,4);
  TCanvas *can13 = new TCanvas("can10","can10",1200,700);
  can13->Divide(5,4);
*/
  for (int i=1;i<10;i++){
  	
	pt[i-1]=(i+1)*20.;
        errpt[i-1]=20;

	TH1 *HistPt11=(TH1*)fResults1->Get(Form("hBias_Pt%i_eta1",(i+1)*20));
  	//RooDataHist *HistPtclone11=new RooDataHist("","",*var,HistPt11);
  	//RooDoubleCBFast *CB11=FitWithDCB(HistPtclone11,bias,mean11[i-1],errmean11[i-1]);  
  	sigeff11[i-1]=effSigma(HistPt11);
	errsigeff11[i-1]=0.4/199;

	TH1 *HistPt12=(TH1*)fResults1->Get(Form("hBias_Pt%i_eta2",(i+1)*20));
  	//RooDataHist *HistPtclone12=new RooDataHist("","",*var,HistPt12);
  	//RooDoubleCBFast *CB12=FitWithDCB(HistPtclone12,bias,mean12[i-1],errmean12[i-1]);  
  	sigeff12[i-1]=effSigma(HistPt12);
	errsigeff12[i-1]=0.4/99;

	TH1 *HistPt13=(TH1*)fResults1->Get(Form("hBias_Pt%i_eta3",(i+1)*20));
  	//RooDataHist *HistPtclone13=new RooDataHist("","",*var,HistPt13);
  	//RooDoubleCBFast *CB13=FitWithDCB(HistPtclone13,bias,mean13[i-1],errmean13[i-1]);  
  	sigeff13[i-1]=effSigma(HistPt13);
	errsigeff13[i-1]=0.4/99;
 

	TH1 *HistPt14=(TH1*)fResults1->Get(Form("hBias_Pt%i_eta4",(i+1)*20));
  	//RooDataHist *HistPtclone14=new RooDataHist("","",*var,HistPt14);
	//RooDoubleCBFast *CB14=FitWithDCB(HistPtclone14,bias,mean14[i-1],errmean14[i-1]);  
  	sigeff14[i-1]=effSigma(HistPt14);
	errsigeff14[i-1]=0.4/99;	

	TH1 *HistPt21=(TH1*)fResults2->Get(Form("hBias_Pt%i_eta1",(i+1)*20));
  	//RooDataHist *HistPtclone21=new RooDataHist("","",*var,HistPt21);
  	//RooDoubleCBFast *CB21=FitWithDCB(HistPtclone21,bias,mean21[i-1],errmean21[i-1]);  
  	sigeff21[i-1]=effSigma(HistPt21);
	errsigeff21[i-1]=0.4/199;

	TH1 *HistPt22=(TH1*)fResults2->Get(Form("hBias_Pt%i_eta2",(i+1)*20));
  	//RooDataHist *HistPtclone22=new RooDataHist("","",*var,HistPt22);
  	//RooDoubleCBFast *CB22=FitWithDCB(HistPtclone22,bias,mean22[i-1],errmean22[i-1]);  
  	sigeff22[i-1]=effSigma(HistPt22);
	errsigeff22[i-1]=0.4/99;

	TH1 *HistPt23=(TH1*)fResults2->Get(Form("hBias_Pt%i_eta3",(i+1)*20));
  	//RooDataHist *HistPtclone23=new RooDataHist("","",*var,HistPt23);
  	//RooDoubleCBFast *CB23=FitWithDCB(HistPtclone23,bias,mean23[i-1],errmean23[i-1]);  
  	sigeff23[i-1]=effSigma(HistPt23);
	errsigeff23[i-1]=0.4/99;
 

	TH1 *HistPt24=(TH1*)fResults2->Get(Form("hBias_Pt%i_eta4",(i+1)*20));
  	//RooDataHist *HistPtclone24=new RooDataHist("","",*var,HistPt24);
	//RooDoubleCBFast *CB24=FitWithDCB(HistPtclone24,bias,mean24[i-1],errmean24[i-1]);  
  	sigeff24[i-1]=effSigma(HistPt24);
	errsigeff24[i-1]=0.4/99;

	TH1 *HistPt31=(TH1*)fResults3->Get(Form("hBias_Pt%i_eta1",(i+1)*20));
  	//RooDataHist *HistPtclone31=new RooDataHist("","",*var,HistPt31);
  	//RooDoubleCBFast *CB31=FitWithDCB(HistPtclone31,bias,mean31[i-1],errmean31[i-1]);  
  	sigeff31[i-1]=effSigma(HistPt31);
	errsigeff31[i-1]=0.4/199;

	TH1 *HistPt32=(TH1*)fResults3->Get(Form("hBias_Pt%i_eta2",(i+1)*20));
  	//RooDataHist *HistPtclone32=new RooDataHist("","",*var,HistPt32);
  	//RooDoubleCBFast *CB32=FitWithDCB(HistPtclone32,bias,mean32[i-1],errmean32[i-1]);  
  	sigeff32[i-1]=effSigma(HistPt32);
	errsigeff32[i-1]=0.4/99;

	TH1 *HistPt33=(TH1*)fResults3->Get(Form("hBias_Pt%i_eta3",(i+1)*20));
  	//RooDataHist *HistPtclone33=new RooDataHist("","",*var,HistPt33);
  	//RooDoubleCBFast *CB33=FitWithDCB(HistPtclone33,bias,mean33[i-1],errmean33[i-1]);  
  	sigeff33[i-1]=effSigma(HistPt33);
	errsigeff33[i-1]=0.4/99;
 

	TH1 *HistPt34=(TH1*)fResults3->Get(Form("hBias_Pt%i_eta4",(i+1)*20));
  	//RooDataHist *HistPtclone34=new RooDataHist("","",*var,HistPt34);
	//RooDoubleCBFast *CB34=FitWithDCB(HistPtclone34,bias,mean34[i-1],errmean34[i-1]);  
  	sigeff34[i-1]=effSigma(HistPt34);
	errsigeff34[i-1]=0.4/99;
  }



  TCanvas *can15=new TCanvas;
  can15->Divide(2,2);
  can15->cd(1);
  TGraphErrors *graphRes_Pt11 = new TGraphErrors(9,pt,sigeff11,errpt,errsigeff11);
  graphRes_Pt11->SetMarkerStyle(8);
  graphRes_Pt11->SetMarkerColor(kBlue+2);
  graphRes_Pt11->SetLineColor(kBlue+2);
  TGraphErrors *graphRes_Pt21 = new TGraphErrors(9,pt,sigeff21,errpt,errsigeff21);
  graphRes_Pt21->SetMarkerStyle(8);
  graphRes_Pt21->SetMarkerColor(kRed+1);
  graphRes_Pt21->SetLineColor(kRed+1);
  TGraphErrors *graphRes_Pt31 = new TGraphErrors(9,pt,sigeff31,errpt,errsigeff31);
  graphRes_Pt31->SetMarkerStyle(8);
  graphRes_Pt31->SetMarkerColor(kOrange+1);
  graphRes_Pt31->SetLineColor(kOrange+1);
  TMultiGraph *graphRes_Pt1=new TMultiGraph();
  graphRes_Pt1->Add(graphRes_Pt11);
  graphRes_Pt1->Add(graphRes_Pt21);
  graphRes_Pt1->Add(graphRes_Pt31);
  graphRes_Pt1->Draw("AP");
  //graphRes_Pt1->GetXaxis()->SetTitle("pt");
  //graphRes_Pt1->GetYaxis()->SetTitle("Ecor/Etrue effective width");
  graphRes_Pt1->SetTitle("Ecor/Etrue effective width vs Pt for |#eta|<1;pt;Ecor/Etrue effective width"); 
  TLegend *leg3 = new TLegend(0.69,0.69,0.89,0.89);
  leg3->SetBorderSize(0);
  leg3->SetFillColor(0);
  leg3->AddEntry(graphRes_Pt11, "no PU","lp");
  leg3->AddEntry(graphRes_Pt21, "bx 50","lp");
  leg3->AddEntry(graphRes_Pt31, "bx 25","lp");
  leg3->Draw(); 
  
  can15->cd(2);
  TGraphErrors *graphRes_Pt12 = new TGraphErrors(9,pt,sigeff12,errpt,errsigeff12);
  graphRes_Pt12->SetMarkerStyle(8);
  graphRes_Pt12->SetMarkerColor(kBlue+2);
  graphRes_Pt12->SetLineColor(kBlue+2);
  TGraphErrors *graphRes_Pt22 = new TGraphErrors(9,pt,sigeff22,errpt,errsigeff22);
  graphRes_Pt22->SetMarkerStyle(8);
  graphRes_Pt22->SetMarkerColor(kRed+1);
  graphRes_Pt22->SetLineColor(kRed+1);
  TGraphErrors *graphRes_Pt32 = new TGraphErrors(9,pt,sigeff32,errpt,errsigeff32);
  graphRes_Pt32->SetMarkerStyle(8);
  graphRes_Pt32->SetMarkerColor(kOrange+1);
  graphRes_Pt32->SetLineColor(kOrange+1);
  TMultiGraph *graphRes_Pt2=new TMultiGraph();
  graphRes_Pt2->Add(graphRes_Pt12);
  graphRes_Pt2->Add(graphRes_Pt22);
  graphRes_Pt2->Add(graphRes_Pt32);
  graphRes_Pt2->Draw("AP");
  //graphRes_Pt2->GetXaxis()->SetTitle("pt");
  //graphRes_Pt2->GetYaxis()->SetTitle("Ecor/Etrue effective width");
  graphRes_Pt2->SetTitle("Ecor/Etrue effective width vs Pt for 1<|#eta|<1.5;pt;Ecor/Etrue effective width"); 
  leg3->Draw();

  can15->cd(3);
  TGraphErrors *graphRes_Pt13 = new TGraphErrors(9,pt,sigeff13,errpt,errsigeff13);
  graphRes_Pt13->SetMarkerStyle(8);
  graphRes_Pt13->SetMarkerColor(kBlue+2);
  graphRes_Pt13->SetLineColor(kBlue+2);
  TGraphErrors *graphRes_Pt23 = new TGraphErrors(9,pt,sigeff23,errpt,errsigeff23);
  graphRes_Pt23->SetMarkerStyle(8);
  graphRes_Pt23->SetMarkerColor(kRed+1);
  graphRes_Pt23->SetLineColor(kRed+1);
  TGraphErrors *graphRes_Pt33 = new TGraphErrors(9,pt,sigeff33,errpt,errsigeff33);
  graphRes_Pt33->SetMarkerStyle(8);
  graphRes_Pt33->SetMarkerColor(kOrange+1);
  graphRes_Pt33->SetLineColor(kOrange+1);
  TMultiGraph *graphRes_Pt3=new TMultiGraph();
  graphRes_Pt3->Add(graphRes_Pt13);
  graphRes_Pt3->Add(graphRes_Pt23);
  graphRes_Pt3->Add(graphRes_Pt33);
  graphRes_Pt3->Draw("AP");
  //graphRes_Pt3->GetXaxis()->SetTitle("pt");
  //graphRes_Pt3->GetYaxis()->SetTitle("Ecor/Etrue effective width");
  graphRes_Pt3->SetTitle("Ecor/Etrue effective width vs Pt for 1.5<|#eta|<2;pt;Ecor/Etrue effective width"); 
  leg3->Draw(); 

  can15->cd(4);
  TGraphErrors *graphRes_Pt14 = new TGraphErrors(9,pt,sigeff14,errpt,errsigeff14);
  graphRes_Pt14->SetMarkerStyle(8);
  graphRes_Pt14->SetMarkerColor(kBlue+2);
  graphRes_Pt14->SetLineColor(kBlue+2);
  TGraphErrors *graphRes_Pt24 = new TGraphErrors(9,pt,sigeff24,errpt,errsigeff24);
  graphRes_Pt24->SetMarkerStyle(8);
  graphRes_Pt24->SetMarkerColor(kRed+1);
  graphRes_Pt24->SetLineColor(kRed+1);
  TGraphErrors *graphRes_Pt34 = new TGraphErrors(9,pt,sigeff34,errpt,errsigeff34);
  graphRes_Pt34->SetMarkerStyle(8);
  graphRes_Pt34->SetMarkerColor(kOrange+1);
  graphRes_Pt34->SetLineColor(kOrange+1);
  TMultiGraph *graphRes_Pt4=new TMultiGraph();
  graphRes_Pt4->Add(graphRes_Pt14);
  graphRes_Pt4->Add(graphRes_Pt24);
  graphRes_Pt4->Add(graphRes_Pt34);
  graphRes_Pt4->Draw("AP");
  //graphRes_Pt4->GetXaxis()->SetTitle("pt");
  //graphRes_Pt4->GetYaxis()->SetTitle("Ecor/Etrue effective width");
  graphRes_Pt4->SetTitle("Ecor/Etrue effective width vs Pt for 1.5<|#eta|;pt;Ecor/Etrue effective width"); 
  leg3->Draw();

  can15->SaveAs(Form("plotsAll/%s/Res_Pt_Eta.pdf",dir.c_str()));
  can15->SaveAs(Form("plotsAll/%s/Res_Pt_Eta.png",dir.c_str()));

/*
  TCanvas *can14=new TCanvas;
  can14->Divide(2,2);
  can14->cd(1);
  TGraphErrors *graphCor_Pt1 = new TGraphErrors(9,pt,mean1,errpt,errmean1);
  graphCor_Pt1->Draw("AL*");
  graphCor_Pt1->SetMarkerStyle(8);
  graphCor_Pt1->GetXaxis()->SetTitle("pt");
  graphCor_Pt1->GetYaxis()->SetTitle("Ecor/Etrue mean");
  graphCor_Pt1->SetTitle("Ecor/Etrue peak value vs Pt for 0<|eta|<1"); 
  can14->cd(2);
  TGraphErrors *graphCor_Pt2 = new TGraphErrors(9,pt,mean2,errpt,errmean2);
  graphCor_Pt2->Draw("AL*");
  graphCor_Pt2->SetMarkerStyle(8);
  graphCor_Pt2->GetXaxis()->SetTitle("pt");
  graphCor_Pt2->GetYaxis()->SetTitle("Ecor/Etrue mean");
  graphCor_Pt2->SetTitle("Ecor/Etrue peak value vs Pt for 1<|eta|<1.48"); 
  can14->cd(3);
  TGraphErrors *graphCor_Pt3 = new TGraphErrors(9,pt,mean3,errpt,errmean3);
  graphCor_Pt3->Draw("AL*");
  graphCor_Pt3->SetMarkerStyle(8);
  graphCor_Pt3->GetXaxis()->SetTitle("pt");
  graphCor_Pt3->GetYaxis()->SetTitle("Ecor/Etrue mean");
  graphCor_Pt3->SetTitle("Ecor/Etrue peak value vs Pt for 1.48<|eta|<2"); 
  can14->cd(4);
  TGraphErrors *graphCor_Pt4 = new TGraphErrors(9,pt,mean4,errpt,errmean4);
  graphCor_Pt4->Draw("AL*");
  graphCor_Pt4->SetMarkerStyle(8);
  graphCor_Pt4->GetXaxis()->SetTitle("pt");
  graphCor_Pt4->GetYaxis()->SetTitle("Ecor/Etrue mean");
  graphCor_Pt4->SetTitle("Ecor/Etrue peak value vs Pt for 2<|eta|<2.5"); 
  can14->SaveAs(Form("plotsAll/%s/Bias_Pt_Eta.png",dir.c_str()));
*/
}

