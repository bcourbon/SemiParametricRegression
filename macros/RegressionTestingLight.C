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
 

void RegressionTestingLight(bool dobarrel=true) {
  
  TString fname;

  if (dobarrel) 
    fname = "wereg_ph_eb_bx50_new.root";
  else if (!dobarrel) 
    fname = "wereg_ph_ee_bx50_new.root";
  
  TFile *fws = TFile::Open(fname); 
  RooWorkspace *ws = (RooWorkspace*)fws->Get("wereg");

  //read variables from workspace
  RooGBRTargetFlex *meantgt = static_cast<RooGBRTargetFlex*>(ws->arg("sigmeant"));  
  
  RooRealVar *tgtvar = ws->var("tgtvar"); //Etrue/Eraw
  RooRealVar *erawvar = ws->var("var_0"); //Eraw

  RooArgList vars;
  vars.add(meantgt->FuncVars());
  vars.add(*tgtvar);
 
  //read testing dataset from TTree
  RooRealVar weightvar("weightvar","",1.);

  TTree *dtree;
  
    TFile *fdin = TFile::Open("Trees/RegressionPhoton_ntuple_bx50_new.root");
    TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("gedPhotonTree");
    dtree = (TTree*)ddir->Get("RegressionTree");       
  
  //selection cuts for testing
  TCut selcut;
  if (dobarrel) 
    selcut = "genPt>20. && scIsEB"; 
  else
    selcut = "genPt>20. && !scIsEB"; 
  
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

  weightvar.SetTitle(oddevents*selcut);    
    
  //retrieve full pdf from workspace
  RooAbsPdf *sigpdf = ws->pdf("sigpdf");
  
  //regressed output functions
  RooAbsReal *sigmeanlim = ws->function("sigmeanlim"); //Ecor/Eraw
  RooAbsReal *sigwidthlim = ws->function("sigwidthlim"); //sigma(E)/Eraw
  RooAbsReal *signlim = ws->function("signlim");
  RooAbsReal *sign2lim = ws->function("sign2lim");

  RooRealVar *meanvar = (RooRealVar*)hdata->addColumn(*sigmeanlim);
  RooRealVar *widthvar = (RooRealVar*)hdata->addColumn(*sigwidthlim);

  //Ecor/Etrue ( =1.0/(etrue/eraw) * regression mean)
  RooFormulaVar cor("cor","","1./(@0)*@1",RooArgList(*tgtvar,*meanvar));
  RooRealVar *corvar = (RooRealVar*)hdata->addColumn(cor);
  TH1 *hcor = hdata->createHistogram("hcor",*corvar);  
  
  //Eraw/Etrue ( =1.0/(etrue/eraw))
  RooFormulaVar raw("raw","","1./@0",RooArgList(*tgtvar));
  RooRealVar *rawvar = (RooRealVar*)hdata->addColumn(raw);
  TH1 *hraw = hdata->createHistogram("hraw",*rawvar);  

  //sigma(E)/Ecor (=1.0/(Ecor/Etrue) * regression width) 
  RooFormulaVar res("res","","1./(@0)*@1",RooArgList(*meanvar,*widthvar));
  RooRealVar *resvar = (RooRealVar*)hdataclone->addColumn(res);
  TH1 *hres = hdata->createHistogram("hres",*resvar);  

  //Ecor = Eraw * regression mean)
  RooFormulaVar ecor("cor","","@0*@1",RooArgList(*eraw,*meanvar));
  RooRealVar *ecorvar = (RooRealVar*)hdata->addColumn(ecor);  
  TH1 *hecor = hdata->createHistogram("hecor",*ecorvar);  

  TH1 *heraw = hdata->createHistogram("heraw",*erawvar);  

  //Sigma(E) = Eraw * regression width)
  RooFormulaVar eres("eres","","@0*@1",RooArgList(*eraw,*widthvar));
  RooRealVar *ecorvar = (RooRealVar*)hdata->addColumn(ecor);
  TH1 *heres = hdata->createHistogram("heres",*eresvar);  

  TString filename;
  if (dobarrel) 
    filename = "resultsEB_bx25_test.root";
  else if (!dobarrel) 
    filename = "resultsEE_bx25_test.root";

  TFile *file=new TFile(name, "RECREATE");

  hecor->Write("Ecor");
  heraw->Write("Eraw");
  heres->Write("Res");
  hcor->Write("Ecor/Etrue");
  hraw->Write("Eraw/Etrue");
  hres->Write("Res/Ecor");

  file->Close();




}
