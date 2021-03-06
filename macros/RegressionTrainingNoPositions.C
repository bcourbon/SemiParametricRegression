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
#include "RooAddPdf.h"
#include "RooNDKeysPdf.h"
#include "RooExtendPdf.h"
#include "RooMinimizer.h"
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
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TFile.h"
#include "TStyle.h"
 
using namespace RooFit;
  
double getweight(TFile *file, double xsec) {
 
  TDirectory *dir = (TDirectory*)file->FindObjectAny("AnaFwkMod");
  TH1D *hallevts = (TH1D*)dir->Get("hDAllEvents");
  
  return xsec/hallevts->GetSumOfWeights();
  
}

float xsecweights[50];
float xsecweight(int procidx=0) {
  return xsecweights[procidx];
}

void initweights(TChain *chain, float *xsecs, float lumi) {
 
  TObjArray *files = chain->GetListOfFiles();
  for (int i=0; i<files->GetEntries(); ++i) {    
    TFile *file = TFile::Open(files->At(i)->GetTitle(),"READ");
    
    xsecweights[i] = getweight(file,lumi*xsecs[i]);
    
    file->Close();    
  } 
  
  chain->SetAlias("procidx","This->GetTreeNumber()");
  
}

void RegressionTrainingNoLocal(bool dobarrel=true) {
   
  
  //build vectors with list of input variables
  std::vector<std::string> *varsf = new std::vector<std::string>;
  varsf->push_back("scRawEnergy");
  //varsf->push_back("scEta");
  //varsf->push_back("scPhi");
  varsf->push_back("scSeedR9");  
  varsf->push_back("scEtaWidth");
  varsf->push_back("scPhiWidth");  
  varsf->push_back("N_ECALClusters");
  varsf->push_back("hoe");
  varsf->push_back("rho");
  varsf->push_back("nVtx");  
 
  varsf->push_back("scSeedEta-scEta");
  varsf->push_back("atan2(sin(scSeedPhi-scPhi),cos(scSeedPhi-scPhi))");

  varsf->push_back("scSeedRawEnergy/scRawEnergy");
  
  varsf->push_back("scSeedE3x3/scSeedE5x5");
  varsf->push_back("scSeedSigmaIetaIeta");   
  varsf->push_back("scSeedSigmaIphiIphi");   
  varsf->push_back("scSeedSigmaIetaIphi");
  varsf->push_back("scSeedEmax/scSeedE5x5");
  varsf->push_back("scSeedE2nd/scSeedE5x5");
  varsf->push_back("scSeedETop/scSeedE5x5");
  varsf->push_back("scSeedEBottom/scSeedE5x5");
  varsf->push_back("scSeedELeft/scSeedE5x5");
  varsf->push_back("scSeedERight/scSeedE5x5");
  varsf->push_back("scSeedE2x5max/scSeedE5x5");
  varsf->push_back("scSeedE2x5Left/scSeedE5x5");
  varsf->push_back("scSeedE2x5Right/scSeedE5x5");
  varsf->push_back("scSeedE2x5Top/scSeedE5x5");
  varsf->push_back("scSeedE2x5Bottom/scSeedE5x5");
  
  std::vector<std::string> *varseb = new std::vector<std::string>(*varsf);
  std::vector<std::string> *varsee = new std::vector<std::string>(*varsf);
  
  varseb->push_back("scSeedE5x5/scSeedRawEnergy");
  
  varseb->push_back("scSeedCryIeta_glob");
  varseb->push_back("scSeedCryIphi_glob");
  varseb->push_back("(scSeedCryIeta_glob-1*abs(scSeedCryIeta_glob)/scSeedCryIeta_glob)%5");
  varseb->push_back("(scSeedCryIphi_glob-1)%2");       
  varseb->push_back("(abs(scSeedCryIeta_glob)<=25)*((scSeedCryIeta_glob-1*abs(scSeedCryIeta_glob)/scSeedCryIeta_glob)%25) + (abs(scSeedCryIeta_glob)>25)*((scSeedCryIeta_glob-26*abs(scSeedCryIeta_glob)/scSeedCryIeta_glob)%20)");
  varseb->push_back("(scSeedCryIphi_glob-1)%20"); 
  //varseb->push_back("scSeedCryPhi");
  //varseb->push_back("scSeedCryEta");

  varseb->push_back("scSeedCryIeta_glob");
  varseb->push_back("scSeedCryIphi_glob");
  varseb->push_back("scSeedCryIeta_glob-ixClusterSeed_log");
  varseb->push_back("scSeedCryIphi_glob-iyClusterSeed_log");


  varsee->push_back("scPreshowerEnergy/scRawEnergy");
  varsee->push_back("scPreshowerEnergyPlane1/scRawEnergy");
  varsee->push_back("scPreshowerEnergyPlane2/scRawEnergy");
    
  //varsee->push_back("scSeedCryX");
  //varsee->push_back("scSeedCryY");

  varsee->push_back("scSeedCryIx_glob");
  varsee->push_back("scSeedCryIy_glob");
  varsee->push_back("scSeedCryIx_glob-ixClusterSeed_log");
  varsee->push_back("scSeedCryIy_glob-iyClusterSeed_log");


  //select appropriate input list for barrel or endcap
  std::vector<std::string> *varslist;
  if (dobarrel) varslist = varseb;
  else varslist = varsee;
  
  //create RooRealVars for each input variable
  RooArgList vars;
  for (unsigned int ivar=0; ivar<varslist->size(); ++ivar) {
    RooRealVar *var = new RooRealVar(TString::Format("var_%i",ivar),varslist->at(ivar).c_str(),0.);
    vars.addOwned(*var);
  }
  
  //make list of input variable RooRealVars
  RooArgList condvars(vars);
  
  //create RooRealVar for target
  RooRealVar *tgtvar = new RooRealVar("tgtvar","genEnergy/scRawEnergy",1.);
  if (!dobarrel) tgtvar->SetTitle("genEnergy/(scRawEnergy + scPreshowerEnergy)");  
  
  //add target to full list
  vars.addOwned(*tgtvar);
    
  //RooRealVar for event weight 
  RooRealVar weightvar("weightvar","",1.);

  //Initialize TChains with event weights if needed
  TString treeloc;


  TChain *tree;

     tree = new TChain("gedPhotonTree/RegressionTree");
     tree->Add("Trees/RegressionPhoton_ntuple_bx25_v6.root");
     //tree->Add("Trees/test_Ben.root");
 /*   
  float xsecs[50];
  xsecs[0] = 0.001835*81930.0;
    xsecs[1] = 0.05387*8884.0;    
    initweights(tree,xsecs,1.);  
    
    double weightscale = xsecweights[1];
    xsecweights[0] /= weightscale;
    xsecweights[1] /= weightscale;
*/
  
  //training selection cut
  TCut selcut;
  if (dobarrel) {
    selcut = "genPt>11. && scIsEB"/* && scRawEnergy/genEnergy>0.25"*/; 
  }
  else {
    selcut = "genPt>11. && !scIsEB"/* && scRawEnergy/genEnergy>0.25"*/;     
  }
  
  
  //TCut selweight = "xsecweight(procidx)";
    TCut prescale10 = "(eventNumber%10==0)";
  TCut prescale20 = "(eventNumber%20==0)";
  TCut prescale25 = "(eventNumber%25==0)";
  TCut prescale50 = "(eventNumber%50==0)";
  TCut prescale100 = "(eventNumber%100==0)";  
  TCut prescale1000 = "(eventNumber%1000==0)";  
  TCut evenevents = "(eventNumber%2==0)";
  TCut oddevents = "(eventNumber%2==1)";  


  //weightvar title used for per-event weights and selection cuts

    //weightvar.SetTitle(evenevents/**selweight*/*selcut);
    weightvar.SetTitle(evenevents/**selweight*/*selcut);

  //create RooDataSet from TChain
  RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",tree,vars,weightvar);   
  

  //RooRealVars corresponding to regressed parameters (sigma, mean, left tail parameter, right tail parameter)
  RooRealVar sigwidthtvar("sigwidthtvar","",0.01);
  sigwidthtvar.setConstant(false);
  
  RooRealVar sigmeantvar("sigmeantvar","",1.);
  sigmeantvar.setConstant(false); 
  
  RooRealVar signvar("signvar","",3.);
  signvar.setConstant(false);       
  
  RooRealVar sign2var("sign2var","",3.);
  sign2var.setConstant(false);     

  //define non-parametric functions for each regressed parameter
  RooGBRFunctionFlex *sigwidthtfunc = new RooGBRFunctionFlex("sigwidthtfunc","");
  RooGBRFunctionFlex *sigmeantfunc = new RooGBRFunctionFlex("sigmeantfunc","");
  RooGBRFunctionFlex *signfunc = new RooGBRFunctionFlex("signfunc","");
  RooGBRFunctionFlex *sign2func = new RooGBRFunctionFlex("sign2func","");

  //define mapping of input variables to non-parametric functions (in this case trivial since all 4 functions depend on the same inputs, but this is not a requirement)
  RooGBRTargetFlex *sigwidtht = new RooGBRTargetFlex("sigwidtht","",*sigwidthtfunc,sigwidthtvar,condvars);  
  RooGBRTargetFlex *sigmeant = new RooGBRTargetFlex("sigmeant","",*sigmeantfunc,sigmeantvar,condvars);  
  RooGBRTargetFlex *signt = new RooGBRTargetFlex("signt","",*signfunc,signvar,condvars);  
  RooGBRTargetFlex *sign2t = new RooGBRTargetFlex("sign2t","",*sign2func,sign2var,condvars);  

  //define list of mapped functions to regress
  RooArgList tgts;
  tgts.add(*sigwidtht);
  tgts.add(*sigmeant);
  tgts.add(*signt);
  tgts.add(*sign2t);  
  
  //define transformations corresponding to parameter bounds for non-parametric outputs  
  RooRealConstraint sigwidthlim("sigwidthlim","",*sigwidtht,0.0002,0.5);
  RooRealConstraint sigmeanlim("sigmeanlim","",*sigmeant,0.2,2.0);
  RooRealConstraint signlim("signlim","",*signt,1.01,50.); 
  RooRealConstraint sign2lim("sign2lim","",*sign2t,1.01,50.); 

  //define pdf, which depends on transformed outputs (and is intended to be treated as a conditional pdf over the
  //regression inputs in this case)
  //The actual pdf below is a double crystal ball, with crossover points alpha_1 and alpha_2 set constant, but all other
  //parameters regressed
  RooDoubleCBFast sigpdf("sigpdf","",*tgtvar,sigmeanlim,sigwidthlim,RooConst(2.),signlim,RooConst(1.),sign2lim);
  
  //dummy variable
  RooConstVar etermconst("etermconst","",0.);  
   
  //dummy variable
  RooRealVar r("r","",1.);
  r.setConstant();

  //define list of pdfs
  std::vector<RooAbsReal*> vpdf;
  vpdf.push_back(&sigpdf);  

  //define list of training datasets
  std::vector<RooAbsData*> vdata;
  vdata.push_back(hdata);     
  
  //define minimum event weight per tree node
  double minweight = 200;//200
  std::vector<double> minweights;
  minweights.push_back(minweight);
  

  //run training
  if (1) {
    RooHybridBDTAutoPdf bdtpdfdiff("bdtpdfdiff","",tgts,etermconst,r,vdata,vpdf);
    bdtpdfdiff.SetMinCutSignificance(5.);
    //bdtpdfdiff.SetPrescaleInit(100);
    bdtpdfdiff.SetShrinkage(0.1);
    bdtpdfdiff.SetMinWeights(minweights);
    bdtpdfdiff.SetMaxNodes(750);
    bdtpdfdiff.TrainForest(1000);   
  }
     
  //create workspace and output to file
  RooWorkspace *wereg = new RooWorkspace("wereg");
  wereg->import(sigpdf);

    
  if (dobarrel)
    wereg->writeToFile("wereg_ph_eb_bx25_v6_nopositions.root");    
  else if (!dobarrel)
    wereg->writeToFile("wereg_ph_ee_bx25_v6_nopositions.root");    
  
  
  return;
  
  
}
