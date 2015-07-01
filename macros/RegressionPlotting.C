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
#include "TText.h"
#include "RooDoubleCBFast.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "TGraphErrors.h"
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



void RegressionPlotting(bool doEBEE=true, bool doEB=false, bool doEE=false){

  TString name;
  string dir;

  if (doEBEE){ 
  name="results_bx50.root";
  dir="plots";
  }
  if (doEB){ 
  name="resultsEB_bx50.root";
  dir="plotsEB";
  }
  if (doEE){ 
  name="resultsEE_bx50.root";
  dir="plotsEE";
  }

  TFile *fResults = TFile::Open(name);

  //Plot eraw/etrue and ecor/etrue to quantify regression performance
  
  //Get Histograms
  TH1 *hecor = (TH1*)fResults->Get("hCor_true");
  TH1 *heraw = (TH1*)fResults->Get("hRaw_true");
  //Compute peak of the distributions by fitting histograms with double Crystal Ball
  double PeakRaw,PeakRawErr,PeakCor,PeakCorErr;
  RooRealVar *bias=new RooRealVar("Ecor/Etrue","Ecor/Etrue",0.8,1.2);
  const RooArgList *var=new RooArgList(*bias,"");
  RooPlot *plot = bias->frame();
  //plot->SetTitle("E/");
  plot->GetYaxis()->SetTitle("");
  RooDataHist *Hraw=new RooDataHist("","",*var,heraw);
  RooDoubleCBFast *DoubleCBraw=FitWithDCB(Hraw,bias,PeakRaw,PeakRawErr); 
  Hraw->plotOn(plot,MarkerSize(0.8));
  DoubleCBraw->plotOn(plot,LineColor(kBlue),Name("hraw"));
  RooDataHist *Hcor=new RooDataHist("","",*var,hecor);
  RooDoubleCBFast *DoubleCBcor=FitWithDCB(Hcor,bias,PeakCor,PeakCorErr); 
  Hcor->plotOn(plot,MarkerSize(0.8));
  DoubleCBcor->plotOn(plot,LineColor(kRed),Name("hcor"));
  //Compute Effective Sigma of the distributions
  double SigEffRaw=effSigma(heraw);
  double SigEffCor=effSigma(hecor);
  double SigEffErr=0.001;
  //Legend and StatBox
  TLegend *leg = new TLegend(0.11,0.6,0.4,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(plot->findObject("hraw"),"Eraw/Etrue","l");
  leg->AddEntry(plot->findObject("hcor"),"Ecor/Etrue","l");
  leg->SetTextFont(62);
  leg->SetTextSize(0.06);
  plot->addObject(leg);
  TPaveText *stats = new TPaveText(0.58,0.5,0.89,0.89,"NDC");
  stats->SetBorderSize(0);
  stats->SetFillColor(0);
  stats->AddText(Form("Peak Raw = %.4f +/- %.0e",PeakRaw,PeakRawErr));
  stats->AddText(Form("Peak Cor = %.4f +/- %.0e",PeakCor,PeakCorErr));
  stats->AddText(Form("Sigma Eff Raw = %.3f +/- %.3f",SigEffRaw,SigEffErr));
  stats->AddText(Form("Sigma Eff Cor = %.3f +/- %.3f",SigEffCor,SigEffErr));
  stats->SetTextSize(0.035);
  plot->addObject(stats);
  //Plot
  TCanvas *cresponse=new TCanvas;
  plot->Draw();
  cresponse->SaveAs(Form("plots_bx50/%s/Performance.pdf",dir.c_str()));
  cresponse->SaveAs(Form("plots_bx50/%s/Performance.png",dir.c_str()));
  

  //Plot ecor/eraw
/*  TH1 *hecor_raw = (TH1*)fResults->Get("hCor_raw");
  TCanvas *canCor=new TCanvas;
  hecor_raw->Draw("HIST");
  hecor_raw->SetLineWidth(2);
  hecor_raw->GetXaxis()->SetTitle("Ecor/Eraw");
  hecor_raw->SetTitle("Ecor/Eraw");  
  canCor->SaveAs(Form("plots_bx50/%s/Correction.png",dir.c_str()));
*/

  //Plot |ecor-etrue|/eraw and |eraw-etrue|/eraw
  TH1 *hres = (TH1*)fResults->Get("hRes");
  TCanvas *canR=new TCanvas;
  hres->Draw("HIST");
  hres->SetLineWidth(2);
  hres->SetLineColor(kRed);
  hres->GetXaxis()->SetTitle("SigmaE/E");
  hres->GetYaxis()->SetTitle("");
  //hres->SetTitle("Per-Event Resolution");  
  canR->SaveAs(Form("plots_bx50/%s/Resolution.pdf",dir.c_str()));
  canR->SaveAs(Form("plots_bx50/%s/Resolution.png",dir.c_str()));


  // Plot Correction vs input variables


  TCanvas *canCoreta=new TCanvas;
  gStyle->SetOptStat(0);
  TProfile *profCorEta=(TProfile*)fResults->Get("profEta");
  TProfile *profRawEta=(TProfile*)fResults->Get("profEta0");
  profRawEta->SetMinimum(0.6);
  profRawEta->SetMaximum(1.1);
  profRawEta->Draw();
  profCorEta->Draw("same");
  profRawEta->GetXaxis()->SetTitle("Eta");
  profRawEta->GetYaxis()->SetTitle("E/Eraw");
  //profRawEta->SetTitle("Profile of E/Eraw vs Eta");
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
  canCoreta->SaveAs(Form("plots_bx50/%s/CorEta.pdf",dir.c_str()));
  canCoreta->SaveAs(Form("plots_bx50/%s/CorEta.png",dir.c_str()));


  TCanvas *canCorPhi=new TCanvas;
  gStyle->SetOptStat(0);
  TProfile *profCorPhi=(TProfile*)fResults->Get("profPhi");
  TProfile *profRawPhi=(TProfile*)fResults->Get("profPhi0");
  profRawPhi->SetMinimum(0.6);
  profRawPhi->SetMaximum(1.1);
  profRawPhi->Draw();
  profCorPhi->Draw("same");
  profRawPhi->GetXaxis()->SetTitle("Phi");
  profRawPhi->GetYaxis()->SetTitle("E/Eraw");
  //profRawPhi->SetTitle("Profile of E/Eraw vs Phi");
  profCorPhi->SetMarkerStyle(8);
  profCorPhi->SetLineColor(kRed);
  profCorPhi->SetMarkerColor(kRed);
  profRawPhi->SetMarkerStyle(4);
  profRawPhi->SetLineColor(kBlue);
  profRawPhi->SetMarkerColor(kBlue);
  leg2->Draw();
  canCorPhi->SaveAs(Form("plots_bx50/%s/CorPhi.pdf",dir.c_str()));
  canCorPhi->SaveAs(Form("plots_bx50/%s/CorPhi.png",dir.c_str()));

  TCanvas *canCorEtaWidth=new TCanvas;
  gStyle->SetOptStat(0);
  TProfile *profCorEtaWidth=(TProfile*)fResults->Get("profEtaWidth");
  TProfile *profRawEtaWidth=(TProfile*)fResults->Get("profEtaWidth0");
  profRawEtaWidth->Draw();
  profCorEtaWidth->Draw("same");
  profRawEtaWidth->GetXaxis()->SetTitle("EtaWidth");
  profRawEtaWidth->GetYaxis()->SetTitle("E/Eraw");
  //profRawEtaWidth->SetTitle("Profile of E/Eraw vs EtaWidth");
  profCorEtaWidth->SetMarkerStyle(8);
  profCorEtaWidth->SetLineColor(kRed);
  profCorEtaWidth->SetMarkerColor(kRed);
  profRawEtaWidth->SetMarkerStyle(4);
  profRawEtaWidth->SetLineColor(kBlue);
  profRawEtaWidth->SetMarkerColor(kBlue);
  leg2->Draw();
  canCorEtaWidth->SaveAs(Form("plots_bx50/%s/CorEtaWidth.pdf",dir.c_str()));
  canCorEtaWidth->SaveAs(Form("plots_bx50/%s/CorEtaWidth.png",dir.c_str()));

  TCanvas *canCorPhiWidth=new TCanvas;
  gStyle->SetOptStat(0);
  TProfile *profCorPhiWidth=(TProfile*)fResults->Get("profPhiWidth");
  TProfile *profRawPhiWidth=(TProfile*)fResults->Get("profPhiWidth0");
  profRawPhiWidth->SetMinimum(0.4);
  profRawPhiWidth->SetMaximum(1.1);
  profRawPhiWidth->Draw();
  profCorPhiWidth->Draw("same");
  profRawPhiWidth->GetXaxis()->SetTitle("PhiWidth");
  profRawPhiWidth->GetYaxis()->SetTitle("E/Eraw");
  //profRawPhiWidth->SetTitle("Profile of E/Eraw vs PhiWidth");
  profCorPhiWidth->SetMarkerStyle(8);
  profCorPhiWidth->SetLineColor(kRed);
  profCorPhiWidth->SetMarkerColor(kRed);
  profRawPhiWidth->SetMarkerStyle(4);
  profRawPhiWidth->SetLineColor(kBlue);
  profRawPhiWidth->SetMarkerColor(kBlue);
  leg2->Draw();
  canCorPhiWidth->SaveAs(Form("plots_bx50/%s/CorPhiWidth.pdf",dir.c_str()));
  canCorPhiWidth->SaveAs(Form("plots_bx50/%s/CorPhiWidth.png",dir.c_str()));

  TCanvas *canCorR9=new TCanvas;
  gStyle->SetOptStat(0);
  TProfile *profCorR9=(TProfile*)fResults->Get("profR9");
  TProfile *profRawR9=(TProfile*)fResults->Get("profR90");
  profRawR9->SetMinimum(0.2);
  profRawR9->SetMaximum(1.1);
  profRawR9->Draw();
  profCorR9->Draw("same");
  profRawR9->GetXaxis()->SetTitle("R9");
  profRawR9->GetYaxis()->SetTitle("E/Eraw");
  //profRawR9->SetTitle("Profile of E/Eraw vs R9");
  profCorR9->SetMarkerStyle(8);
  profCorR9->SetLineColor(kRed);
  profCorR9->SetMarkerColor(kRed);
  profRawR9->SetMarkerStyle(4);
  profRawR9->SetLineColor(kBlue);
  profRawR9->SetMarkerColor(kBlue);
  leg2->Draw();
  canCorR9->SaveAs(Form("plots_bx50/%s/CorR9.pdf",dir.c_str()));
  canCorR9->SaveAs(Form("plots_bx50/%s/CorR9.png",dir.c_str()));

  TCanvas *canCorNvtx=new TCanvas;
  gStyle->SetOptStat(0);
  TProfile *profCorNvtx=(TProfile*)fResults->Get("profNvtx");
  TProfile *profRawNvtx=(TProfile*)fResults->Get("profNvtx0");
  profRawNvtx->SetMinimum(0.5);
  profRawNvtx->SetMaximum(1.2);
  profRawNvtx->Draw();
  profCorNvtx->Draw("same");
  profRawNvtx->GetXaxis()->SetTitle("Nvtx");
  profRawNvtx->GetYaxis()->SetTitle("E/Eraw");
  //profRawNvtx->SetTitle("Profile of E/Eraw vs Nvtx");
  profCorNvtx->SetMarkerStyle(8);
  profCorNvtx->SetLineColor(kRed);
  profCorNvtx->SetMarkerColor(kRed);
  profRawNvtx->SetMarkerStyle(4);
  profRawNvtx->SetLineColor(kBlue);
  profRawNvtx->SetMarkerColor(kBlue);
  leg2->Draw();
  canCorNvtx->SaveAs(Form("plots_bx50/%s/CorNvtx.pdf",dir.c_str()));
  canCorNvtx->SaveAs(Form("plots_bx50/%s/CorNvtx.png",dir.c_str()));

  TCanvas *canCorPt=new TCanvas;
  gStyle->SetOptStat(0);
  TProfile *profCorPt=(TProfile*)fResults->Get("profPt");
  TProfile *profRawPt=(TProfile*)fResults->Get("profPt0");
  profRawPt->SetMinimum(0.5);
  profRawPt->SetMaximum(1.2);
  profRawPt->Draw();
  profCorPt->Draw("same");
  profRawPt->GetXaxis()->SetTitle("Pt");
  profRawPt->GetYaxis()->SetTitle("E/Eraw");
  //profRawPt->SetTitle("Profile of E/Eraw vs Pt");
  profCorPt->SetMarkerStyle(8);
  profCorPt->SetLineColor(kRed);
  profCorPt->SetMarkerColor(kRed);
  profRawPt->SetMarkerStyle(4);
  profRawPt->SetLineColor(kBlue);
  profRawPt->SetMarkerColor(kBlue);
  leg2->Draw();
  canCorPt->SaveAs(Form("plots_bx50/%s/CorPt.pdf",dir.c_str()));
  canCorPt->SaveAs(Form("plots_bx50/%s/CorPt.png",dir.c_str()));




  //Plot Bias and Resolution vs Pt

 /* TCanvas *canResPt=new TCanvas;
  TProfile *profResPt=(TProfile*)fResults->Get("profPtRes");
  profResPt->Draw();
  profResPt->GetXaxis()->SetTitle("Pt");
  profResPt->GetYaxis()->SetTitle("|Ecor-Etrue|/Eraw");
  profResPt->SetTitle("Profile of SigmaE/E vs Pt");
  profResPt->SetMarkerStyle(3);
  profResPt->SetLineColor(kRed);
  canResPt->SaveAs(Form("plots_bx50/%s/ResPt.png",dir.c_str()));
*/

//Plot Bias and resolution in bins of Pt and eta


  double errpt[9];
  double pt[9];
  double sigeff1[9];
  double errsigeff1[9];
  double sigeff2[9];
  double errsigeff2[9];
  double sigeff3[9];
  double errsigeff3[9];
  double sigeff4[9];
  double errsigeff4[9];
  double mean1[9];
  double errmean1[9];
  double mean2[9];
  double errmean2[9];
  double mean3[9];
  double errmean3[9];
  double mean4[9];
  double errmean4[9];

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

	TH1 *HistPt1=(TH1*)fResults->Get(Form("hBias_Pt%i_eta1",(i+1)*20));
  	RooDataHist *HistPtclone1=new RooDataHist("","",*var,HistPt1);
  	RooDoubleCBFast *CB1=FitWithDCB(HistPtclone1,bias,mean1[i-1],errmean1[i-1]);  
  	sigeff1[i-1]=effSigma(HistPt1);
	errsigeff1[i-1]=0.4/199;
  	//can10->cd(i+1);
	//RooPlot *plot1=bias->frame(Title(Form("%i<Pt<%i and 0<|eta|<1",i*10+5,(i+1)*10+5)));
        //HistPtclone1->plotOn(plot1);
	//CB1->plotOn(plot1);
	//plot1->Draw();

	TH1 *HistPt2=(TH1*)fResults->Get(Form("hBias_Pt%i_eta2",(i+1)*20));
  	RooDataHist *HistPtclone2=new RooDataHist("","",*var,HistPt2);
  	RooDoubleCBFast *CB2=FitWithDCB(HistPtclone2,bias,mean2[i-1],errmean2[i-1]);  
  	sigeff2[i-1]=effSigma(HistPt2);
	errsigeff2[i-1]=0.4/99;
  	//can11->cd(i+1);
	//RooPlot *plot2=bias->frame(Title(Form("%i<Pt<%i and 1<|eta|<1.5",i*10+5,(i+1)*10+5)));
        //HistPtclone2->plotOn(plot2);
	//CB2->plotOn(plot2);
	//plot2->Draw();

	TH1 *HistPt3=(TH1*)fResults->Get(Form("hBias_Pt%i_eta3",(i+1)*20));
  	RooDataHist *HistPtclone3=new RooDataHist("","",*var,HistPt3);
  	RooDoubleCBFast *CB3=FitWithDCB(HistPtclone3,bias,mean3[i-1],errmean3[i-1]);  
  	sigeff3[i-1]=effSigma(HistPt3);
	errsigeff3[i-1]=0.4/99;
  	//can12->cd(i+1);
	//RooPlot *plot3=bias->frame(Title(Form("%i<Pt<%i and 1.5<|eta|<2",i*10+5,(i+1)*10+5)));
        //HistPtclone3->plotOn(plot3);
	//CB3->plotOn(plot3);
	//plot3->Draw();

	TH1 *HistPt4=(TH1*)fResults->Get(Form("hBias_Pt%i_eta4",(i+1)*20));
  	RooDataHist *HistPtclone4=new RooDataHist("","",*var,HistPt4);
	RooDoubleCBFast *CB4=FitWithDCB(HistPtclone4,bias,mean4[i-1],errmean4[i-1]);  
  	sigeff4[i-1]=effSigma(HistPt4);
	errsigeff4[i-1]=0.4/99;
  	//can13->cd(i+1);
	//RooPlot *plot4=bias->frame(Title(Form("%i<Pt<%i and 2<|eta|<2.5",i*10+5,(i+1)*10+5)));
        //HistPtclone4->plotOn(plot4);
	//CB4->plotOn(plot4);
	//plot4->Draw();  	

  }

/*  can10->SaveAs(Form("plots_bx50/%s/PtDistribution1.png",dir.c_str()));
  can11->SaveAs(Form("plots_bx50/%s/PtDistribution2.png",dir.c_str()));
  can12->SaveAs(Form("plots_bx50/%s/PtDistribution3.png",dir.c_str()));
  can13->SaveAs(Form("plots_bx50/%s/PtDistribution4.png",dir.c_str()));
*/


  TCanvas *can15=new TCanvas;
  can15->Divide(2,2);
  can15->cd(1);
  TGraphErrors *graphRes_Pt1 = new TGraphErrors(9,pt,sigeff1,errpt,errsigeff1);
  //graphRes_Pt1->SetMinimum(0.);
  //graphRes_Pt1->SetMaximum(0.05);
  graphRes_Pt1->Draw("AL*");
  graphRes_Pt1->SetMarkerStyle(8);
  graphRes_Pt1->GetXaxis()->SetTitle("pt");
  graphRes_Pt1->GetYaxis()->SetTitle("Ecor/Etrue effective width");
  graphRes_Pt1->SetTitle("Ecor/Etrue effective width vs Pt for 0<|eta|<1"); 
  can15->cd(2);
  TGraphErrors *graphRes_Pt2 = new TGraphErrors(9,pt,sigeff2,errpt,errsigeff2);
  //graphRes_Pt2->SetMinimum(0.);
  //graphRes_Pt2->SetMaximum(0.05);
  graphRes_Pt2->Draw("AL*");
  graphRes_Pt2->SetMarkerStyle(8);
  graphRes_Pt2->GetXaxis()->SetTitle("pt");
  graphRes_Pt2->GetYaxis()->SetTitle("Ecor/Etrue effective width");
  graphRes_Pt2->SetTitle("Ecor/Etrue effective width vs Pt for 1<|eta|<1.48"); 
  can15->cd(3);
  TGraphErrors *graphRes_Pt3 = new TGraphErrors(9,pt,sigeff3,errpt,errsigeff3);
  //graphRes_Pt3->SetMinimum(0.);
  //graphRes_Pt3->SetMaximum(0.05);
  graphRes_Pt3->Draw("AL*");
  graphRes_Pt3->SetMarkerStyle(8);
  graphRes_Pt3->GetXaxis()->SetTitle("pt");
  graphRes_Pt3->GetYaxis()->SetTitle("Ecor/Etrue effective width");
  graphRes_Pt3->SetTitle("Ecor/Etrue effective width vs Pt for 1.48<|eta|<2"); 
  can15->cd(4);
  TGraphErrors *graphRes_Pt4 = new TGraphErrors(9,pt,sigeff4,errpt,errsigeff4);
  //graphRes_Pt4->SetMinimum(0.);
  //graphRes_Pt4->SetMaximum(0.05);
  graphRes_Pt4->Draw("AL*");
  graphRes_Pt4->SetMarkerStyle(8);
  graphRes_Pt4->GetXaxis()->SetTitle("pt");
  graphRes_Pt4->GetYaxis()->SetTitle("Ecor/Etrue effective width");
  graphRes_Pt4->SetTitle("Ecor/Etrue effective width vs Pt for 2<|eta|<2.5"); 
  can15->SaveAs(Form("plots_bx50/%s/Resolution_Pt_Eta.pdf",dir.c_str()));
  can15->SaveAs(Form("plots_bx50/%s/Resolution_Pt_Eta.png",dir.c_str()));

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
  can14->SaveAs(Form("plots_bx50/%s/Bias_Pt_Eta.pdf",dir.c_str()));
  can14->SaveAs(Form("plots_bx50/%s/Bias_Pt_Eta.png",dir.c_str()));

}

