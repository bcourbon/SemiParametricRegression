#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooHybridBDTAutoPdf.h"
#include "RooFormulaVar.h"
#include "RooProdPdf.h"
#include "RooUniform.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TProfile.h"
#include "RooAddPdf.h"
#include "RooNDKeysPdf.h"
#include "RooExtendPdf.h"
#include "RooMinimizer.h"
#include "TStyle.h"
#include "TFile.h"
#include "TNtuple.h"
#include "HybridGBRForest.h"
#include "RooProduct.h"
#include "RooChebychev.h"
#include "RooBernstein.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
//#include "HZZ2L2QRooPdfs.h"
#include "RooDoubleCBFast.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooCBShape.h"
#include "RooWorkspace.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPaveLabel.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TLegend.h"
#include "RooRandom.h"
#include "RooAddition.h"
#include "TSystem.h"
#include "RooLinearVar.h"


using namespace RooFit;
 

void RegressionTesting(bool dobarrel=true) {
  
  TString fname;

  if (dobarrel) 
    fname = "wereg_ph_eb_bx50_new.root";
  else if (!dobarrel) 
    fname = "wereg_ph_ee_bx50_new.root";
  
  
  TFile *fws = TFile::Open(fname); 
  RooWorkspace *ws = (RooWorkspace*)fws->Get("wereg");


  ws->Print();


  //read variables from workspace
  RooGBRTargetFlex *meantgt = static_cast<RooGBRTargetFlex*>(ws->arg("sigmeant"));  
  
  RooRealVar *tgtvar = ws->var("tgtvar");
  RooRealVar *Pt = new RooRealVar("genPt","genPt",1.);
  Pt->setRange(0.,200.);
  Pt->setBins(100);  

  RooArgList vars;
  vars.add(meantgt->FuncVars());
  vars.add(*tgtvar);
  vars.add(*Pt);

 
  //read testing dataset from TTree
  RooRealVar weightvar("weightvar","",1.);

  TTree *dtree;
  
    //TFile *fdin = TFile::Open("Trees/RegressionPhoton_ntuple_bx50_new.root");
    TFile *fdin = TFile::Open("Trees/RegressionPhoton_ntuple_139.root");
    TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("gedPhotonTree");
    dtree = (TTree*)ddir->Get("RegressionTree");       
  
  //selection cuts for testing
  TCut selcut;
  if (dobarrel) 
    selcut = "genPt>20. && scIsEB"; 
  else
    selcut = "genPt>20. && !scIsEB"; 
  
  //TCut selweight = "xsecweight(procidx)*puweight(numPU,procidx)";
  TCut prescale10 = "(eventNumber%10==0)";
  TCut prescale10alt = "(eventNumber%10==1)";
  TCut prescale25 = "(eventNumber%25==0)";
  TCut prescale100 = "(eventNumber%100==0)";  
  TCut prescale1000 = "(eventNumber%1000==0)";  
  TCut evenevents = "(eventNumber%2==0)";
  TCut oddevents = "(eventNumber%2==1)";
  TCut prescale100alt = "(eventNumber%100==1)";
  TCut prescale1000alt = "(eventNumber%1000==1)";
  TCut prescale50alt = "(eventNumber%50==1)";
  
    weightvar.SetTitle(selcut);
  

  //make testing dataset
  RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",dtree,vars,weightvar);   

  //weightvar.SetTitle(oddevents*selcut);
  weightvar.SetTitle(oddevents*selcut);

  //make reduced testing dataset for integration over conditional variables
  //RooDataSet *hdatasmall = RooTreeConvert::CreateDataSet("hdatasmall",dtree,vars,weightvar);     
    
  //retrieve full pdf from workspace
  RooAbsPdf *sigpdf = ws->pdf("sigpdf");
  
  //input variables
  RooRealVar *scetavar = ws->var("var_1");
  scetavar->setBins(100);
  scetavar->setRange(-3.,3.);

  RooRealVar *scphivar = ws->var("var_2");
  scphivar->setBins(100);
  scphivar->setRange(-3.2,3.2);

  RooRealVar *scr9var = ws->var("var_3");
  scr9var->setBins(100);
  scr9var->setRange(0,1.2);

  RooRealVar *scphiwidthvar = ws->var("var_5");
  scphiwidthvar->setBins(100);
  scphiwidthvar->setRange(0,0.1);

  RooRealVar *scetawidthvar = ws->var("var_4");
  scetawidthvar->setBins(100);
  scetawidthvar->setRange(0,0.03);

  RooRealVar *nVtxvar = ws->var("var_9");
  nVtxvar->setBins(50);
  nVtxvar->setRange(0,50);


  //regressed output functions
  RooAbsReal *sigmeanlim = ws->function("sigmeanlim");
  RooAbsReal *sigwidthlim = ws->function("sigwidthlim");
  RooAbsReal *signlim = ws->function("signlim");
  RooAbsReal *sign2lim = ws->function("sign2lim");

  //formula for corrected energy/true energy ( 1.0/(etrue/eraw) * regression mean)
  RooFormulaVar ecor("ecor","","1./(@0)*@1",RooArgList(*tgtvar,*sigmeanlim));
  RooRealVar *ecorvar = (RooRealVar*)hdata->addColumn(ecor);
  ecorvar->setRange(0.,2.);
  ecorvar->setBins(800);
  
  //formula for raw energy/true energy (1.0/(etrue/eraw))
  RooFormulaVar raw("raw","","1./@0",RooArgList(*tgtvar));
  RooRealVar *rawvar = (RooRealVar*)hdata->addColumn(raw);
  rawvar->setRange(0.,2.);
  rawvar->setBins(800);

 //formula for raw resolution (|eraw/etrue-1|)
  RooFormulaVar rawres("rawres","","sqrt(pow((@0-1),2))",RooArgList(*rawvar));
  RooRealVar *rawresvar = (RooRealVar*)hdata->addColumn(rawres);
  rawresvar->setRange(0.,0.05);
  rawresvar->setBins(2500);

  //clone data and add regression outputs for plotting
  RooDataSet *hdataclone = new RooDataSet(*hdata,"hdataclone");
  RooRealVar *meanvar = (RooRealVar*)hdataclone->addColumn(*sigmeanlim);
  //meanvar->setRange(0.8,1.2);
  //meanvar->setBins(800);
  RooRealVar *widthvar = (RooRealVar*)hdataclone->addColumn(*sigwidthlim);
  //widthvar->setRange(0.,0.05);
  //widthvar->setBins(800);

  //formula for corrected resolution (|eraw/etrue-1|)
  RooFormulaVar corres("corres","","1./(@0)*@1",RooArgList(*meanvar,*widthvar));
  RooRealVar *corresvar = (RooRealVar*)hdataclone->addColumn(corres);
  //corresvar->setRange(0.,0.1);
  //corresvar->setBins(500);

  TString name;
  if (dobarrel) 
    name = "resultsEB_bx50_new.root";
  else if (!dobarrel) 
    name = "resultsEE_bx50_new.root";

  TFile *file=new TFile(name, "RECREATE");

  //create histograms for eraw/etrue and ecor/etrue to quantify regression performance
  TH1 *heraw = hdata->createHistogram("hraw",*rawvar,Binning(800,0.,2.));
  TH1 *hecor = hdata->createHistogram("hecor",*ecorvar);  
  hecor->Write("hCor_true");
  heraw->Write("hRaw_true"); 
  
  //create histograms for ecor/eraw
 // TH1 *hecor_raw = hdataclone->createHistogram("hraw",*meanvar,Binning(500,0.8,1.2));
 // hecor_raw->Write("hCor_raw");

  //create histograms for |ecor-etrue|/eraw and |eraw-etrue|/eraw
  TH1 *hres = hdataclone->createHistogram("hres",*corresvar,Binning(500,0,0.1));
  hres->Write("hRes");

  cout<<"+++++++++++++++++++++++++++++++++++++++++="<<endl;
  hdataclone->Print();

  TH1 *hecor_res1=hdataclone->createHistogram("hecor_res1",*ecorvar,Binning(1000,0.5,1.5),"var_1<0.5","");
  hecor_res1->Write("hecor_Res1");



 // TH1 *hrawres = hdataclone->createHistogram("hresraw",*rawresvar,Binning(500,0,0.1));
 // hrawres->Write("hRawRes");

  // Performance vs input variables
  
 /* TH2 *hCor_eta=hdataclone->createHistogram(*scetavar,*meanvar);
  TProfile *profCorEta= new TProfile();
  profCorEta=hCor_eta->ProfileX("",1,-1,"s");
  profCorEta->SetMinimum(0.9);
  profCorEta->SetMaximum(1.1);
  profCorEta->Write("profEta");

  TH2 *hCor_phi=hdataclone->createHistogram(*scphivar,*meanvar);
  TProfile *profCorphi= new TProfile();
  profCorphi=hCor_phi->ProfileX("",1,-1,"s");
  profCorphi->SetMinimum(0.9);
  profCorphi->SetMaximum(1.1);
  profCorphi->Write("profPhi");

  TH2 *hCor_r9=hdataclone->createHistogram(*scr9var,*meanvar);
  TProfile *profCorr9= new TProfile();
  profCorr9=hCor_r9->ProfileX("",1,-1,"s");
  profCorr9->SetMinimum(0.9);
  profCorr9->SetMaximum(1.2);
  profCorr9->Write("profR9");

  TH2 *hCor_phiwidth=hdataclone->createHistogram(*scphiwidthvar,*meanvar);
  TProfile *profCorphiwidth= new TProfile();
  profCorphiwidth=hCor_phiwidth->ProfileX("",1,-1,"s");
  profCorphiwidth->SetMinimum(0.9);
  profCorphiwidth->SetMaximum(1.2);
  profCorphiwidth->Write("profPhiWidth");

  TH2 *hCor_etawidth=hdataclone->createHistogram(*scetawidthvar,*meanvar);
  TProfile *profCoretawidth= new TProfile();
  profCoretawidth=hCor_etawidth->ProfileX("",1,-1,"s");
  profCoretawidth->SetMinimum(0.9);
  profCoretawidth->SetMaximum(1.1);
  profCoretawidth->Write("profEtaWidth");

  TH2 *hCor_nVtx=hdataclone->createHistogram(*nVtxvar,*meanvar);
  TProfile *profCornVtx= new TProfile();
  profCornVtx=hCor_nVtx->ProfileX("",1,-1,"s");
  profCornVtx->SetMinimum(0.9);
  profCornVtx->SetMaximum(1.1);
  profCornVtx->Write("profNvtx");
*/

  TH2 *hCor_eta=hdataclone->createHistogram(*scetavar,*ecorvar);
  TProfile *profCorEta= new TProfile();
  profCorEta=hCor_eta->ProfileX("",1,-1,"s");
  //profCorEta->SetMinimum(0.9);
  //profCorEta->SetMaximum(1.1);
  profCorEta->Write("profEta");

  TH2 *hCor_phi=hdataclone->createHistogram(*scphivar,*ecorvar);
  TProfile *profCorphi= new TProfile();
  profCorphi=hCor_phi->ProfileX("",1,-1,"s");
  //profCorphi->SetMinimum(0.9);
  //profCorphi->SetMaximum(1.1);
  profCorphi->Write("profPhi");

  TH2 *hCor_r9=hdataclone->createHistogram(*scr9var,*ecorvar);
  TProfile *profCorr9= new TProfile();
  profCorr9=hCor_r9->ProfileX("",1,-1,"s");
  //profCorr9->SetMinimum(0.9);
  //profCorr9->SetMaximum(1.2);
  profCorr9->Write("profR9");

  TH2 *hCor_phiwidth=hdataclone->createHistogram(*scphiwidthvar,*ecorvar);
  TProfile *profCorphiwidth= new TProfile();
  profCorphiwidth=hCor_phiwidth->ProfileX("",1,-1,"s");
  //profCorphiwidth->SetMinimum(0.9);
  //profCorphiwidth->SetMaximum(1.2);
  profCorphiwidth->Write("profPhiWidth");

  TH2 *hCor_etawidth=hdataclone->createHistogram(*scetawidthvar,*ecorvar);
  TProfile *profCoretawidth= new TProfile();
  profCoretawidth=hCor_etawidth->ProfileX("",1,-1,"s");
  //profCoretawidth->SetMinimum(0.9);
  //profCoretawidth->SetMaximum(1.1);
  profCoretawidth->Write("profEtaWidth");

  TH2 *hCor_nVtx=hdataclone->createHistogram(*nVtxvar,*ecorvar);
  TProfile *profCornVtx= new TProfile();
  profCornVtx=hCor_nVtx->ProfileX("",1,-1,"s");
  //profCornVtx->SetMinimum(0.9);
  //profCornVtx->SetMaximum(1.1);
  profCornVtx->Write("profNvtx");

  TH2 *hCor_Pt=hdataclone->createHistogram(*Pt,*ecorvar);
  TProfile *profCorPt= new TProfile();
  profCorPt=hCor_Pt->ProfileX("",1,-1,"s");
  //profCorPt->SetMinimum(0.9);
  //profCorPt->SetMaximum(1.1);
  profCorPt->Write("profPt");

  TH2 *hRaw_eta=hdataclone->createHistogram(*scetavar,*rawvar);
  TProfile *profRawEta= new TProfile();
  profRawEta=hRaw_eta->ProfileX("",1,-1,"s");
  //profRawEta->SetMinimum(0.9);
  //profRawEta->SetMaximum(1.1);
  profRawEta->Write("profEta0");

  TH2 *hRaw_phi=hdataclone->createHistogram(*scphivar,*rawvar);
  TProfile *profRawphi= new TProfile();
  profRawphi=hRaw_phi->ProfileX("",1,-1,"s");
  //profRawphi->SetMinimum(0.9);
  //profRawphi->SetMaximum(1.1);
  profRawphi->Write("profPhi0");

  TH2 *hRaw_r9=hdataclone->createHistogram(*scr9var,*rawvar);
  TProfile *profRawr9= new TProfile();
  profRawr9=hRaw_r9->ProfileX("",1,-1,"s");
  //profRawr9->SetMinimum(0.9);
  //profRawr9->SetMaximum(1.2);
  profRawr9->Write("profR90");

  TH2 *hRaw_phiwidth=hdataclone->createHistogram(*scphiwidthvar,*rawvar);
  TProfile *profRawphiwidth= new TProfile();
  profRawphiwidth=hRaw_phiwidth->ProfileX("",1,-1,"s");
  //profRawphiwidth->SetMinimum(0.9);
  //profRawphiwidth->SetMaximum(1.2);
  profRawphiwidth->Write("profPhiWidth0");

  TH2 *hRaw_etawidth=hdataclone->createHistogram(*scetawidthvar,*rawvar);
  TProfile *profRawetawidth= new TProfile();
  profRawetawidth=hRaw_etawidth->ProfileX("",1,-1,"s");
  //profRawetawidth->SetMinimum(0.9);
  //profRawetawidth->SetMaximum(1.1);
  profRawetawidth->Write("profEtaWidth0");

  TH2 *hRaw_nVtx=hdataclone->createHistogram(*nVtxvar,*rawvar);
  TProfile *profRawnVtx= new TProfile();
  profRawnVtx=hRaw_nVtx->ProfileX("",1,-1,"s");
  //profRawnVtx->SetMinimum(0.9);
  //profRawnVtx->SetMaximum(1.1);
  profRawnVtx->Write("profNvtx0");

  TH2 *hRaw_Pt=hdataclone->createHistogram(*Pt,*rawvar);
  TProfile *profRawPt= new TProfile();
  profRawPt=hRaw_Pt->ProfileX("",1,-1,"s");
  //profRawPt->SetMinimum(0.9);
  //profRawPt->SetMaximum(1.1);
  profRawPt->Write("profPt0");

  //Plot SigmaE/E vs Pt
/*
  TH2 *hRaw_Pt=hdataclone->createHistogram(*Pt,*ecorvar);
  TProfile *profCorPt= new TProfile();
  profCorPt=hCor_Pt->ProfileX("",1,-1,"s");
  profCorPt->SetMinimum(0.9);
  profCorPt->SetMaximum(1.1);
  profCorPt->Write("profPtCor");
 
  TH2 *hRes_Pt=hdataclone->createHistogram(*Pt,*widthvar);
  TProfile *profResPt= new TProfile();
  profResPt=hRes_Pt->ProfileX("",1,-1,"s");
  profResPt->SetMinimum(0.);
  profResPt->SetMaximum(0.03);
  profResPt->Write("profPtRes");  

  TH2 *hRawRes_Pt=hdataclone->createHistogram(*Pt,*rawresvar);
  TProfile *profRawResPt= new TProfile();
  profRawResPt=hRawRes_Pt->ProfileX("",1,-1,"s");
  profRawResPt->SetMinimum(0.);
  profRawResPt->SetMaximum(0.03);
  profRawResPt->Write("profPtRawRes");
     
  TH2 *hRes_Pt1=hdataclone->createHistogram(*Pt,*widthvar,100,100,"var_1<1 && var_1>-1","");
  TProfile *profResPt1= new TProfile();Weights
  profResPt1=hRes_Pt1->ProfileX("",1,-1,"s");
  profResPt1->SetMinimum(0.);
  profResPt1->SetMaximum(0.03);
  profResPt1->Write("profPtRes1");

  TH2 *hRes_Pt2=hdataclone->createHistogram(*Pt,*widthvar,100,100,"(var_1<1.5 && var_1>1) || (var_1>-1.5 && var_1<-1) ","");
  TProfile *profResPt2= new TProfile();
  profResPt2=hRes_Pt2->ProfileX("",1,-1,"s");
  profResPt2->SetMinimum(0.);
  profResPt2->SetMaximum(0.03);
  profResPt2->Write("profPtRes2");
 
  TH2 *hRes_Pt3=hdataclone->createHistogram(*Pt,*widthvar,100,100,"(var_1<2 && var_1>1.5) || (var_1>-2 && var_1<-1.5) ","");
  TProfile *profResPt3= new TProfile();
  profResPt3=hRes_Pt3->ProfileX("",1,-1,"s");
  profResPt3->SetMinimum(0.);
  profResPt3->SetMaximum(0.03);
  profResPt3->Write("profPtRes3");

  TH2 *hRes_Pt4=hdataclone->createHistogram(*Pt,*widthvar,100,100,"(var_1<2.5 && var_1>2) || (var_1>-2.5 && var_1<-2) ","");
  TProfile *profResPt4= new TProfile();
  profResPt4=hRes_Pt4->ProfileX("",1,-1,"s");
  profResPt4->SetMinimum(0.);
  profResPt4->SetMaximum(0.03);
  profResPt4->Write("profPtRes4");

*/

//Plot Ecor/Etrue in Pt and Eta bins



  for (int i=1;i<10;i++){
	RooDataSet *dataPt1=(RooDataSet*) hdata->reduce(Form("genPt>%i && genPt<%i && var_1<1 && var_1>-1",i*20+5,(i+1)*20+5));
	RooDataSet *dataPt2=(RooDataSet*) hdata->reduce(Form("genPt>%i && genPt<%i && ((var_1<1.5 && var_1>1) || (var_1>-1.5 && var_1<-1))",i*20+5,(i+1)*20+5));
	RooDataSet *dataPt3=(RooDataSet*) hdata->reduce(Form("genPt>%i && genPt<%i && ((var_1<2 && var_1>1.5) || (var_1>-2 && var_1<-1.5))",i*20+5,(i+1)*20+5));
	RooDataSet *dataPt4=(RooDataSet*) hdata->reduce(Form("genPt>%i && genPt<%i && ((var_1<2.5 && var_1>2) || (var_1>-2.5 && var_1<-2))",i*20+5,(i+1)*20+5));
	TH1* histPt1= dataPt1->createHistogram(Form("hBias_Pt%i_eta1",(i+1)*10),*ecorvar,Binning(199,0.8,1.2));
	TH1* histPt2= dataPt2->createHistogram(Form("hBias_Pt%i_eta2",(i+1)*10),*ecorvar,Binning(99,0.8,1.2));
	TH1* histPt3= dataPt3->createHistogram(Form("hBias_Pt%i_eta3",(i+1)*10),*ecorvar,Binning(99,0.8,1.2));
	TH1* histPt4= dataPt4->createHistogram(Form("hBias_Pt%i_eta4",(i+1)*10),*ecorvar,Binning(99,0.8,1.2));
	histPt1->Write(Form("hBias_Pt%i_eta1",(i+1)*20));
	histPt2->Write(Form("hBias_Pt%i_eta2",(i+1)*20));
	histPt3->Write(Form("hBias_Pt%i_eta3",(i+1)*20));
	histPt4->Write(Form("hBias_Pt%i_eta4",(i+1)*20));
	//delete dataPt1,histPt1,dataPt2,histPt2,dataPt3,histPt3,dataPt4,histPt4;
}


   
  file->Close();

}
