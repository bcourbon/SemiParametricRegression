#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TMath.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TVector3.h"
#include "TProfile.h"
#include "TChain.h"
#include "TROOT.h"

#include <iostream>
#include <sstream>
#include <vector>

void DumpVariablesPerScenario(int sc){

gROOT->ProcessLine("#include <vector>");

   TTree *tree;
   TFile *f; 
   if (sc==1) f = TFile::Open("Trees/RegressionPhoton_ntuple_noPU_new.root");
   if (sc==2) f = TFile::Open("Trees/RegressionPhoton_ntuple_bx50_new.root");
   if (sc==3) f = TFile::Open("Trees/RegressionPhoton_ntuple_bx25_new.root");
   TDirectory *d = (TDirectory*)f->FindObjectAny("gedPhotonTree");
   tree = (TTree*)d->Get("RegressionTree"); 

   float pt;
   float eta;
   float phi;
   float r9;
   float sieie;
   float etaWidth;
   float phiWidth;
   int isMatched;
   float Etrue;
   float Eraw;
   float rho;
   int nvtx;
   float hoe;

   TBranch        *b_pt;  
   TBranch        *b_eta;  
   TBranch        *b_phi;  
   TBranch        *b_r9; 
   TBranch        *b_sieie;
   TBranch        *b_EtaWidth;
   TBranch        *b_PhiWidth;
   TBranch        *b_isMatched;
   TBranch        *b_Etrue;
   TBranch        *b_Eraw;
   TBranch        *b_rho;
   TBranch        *b_nvtx;
   TBranch        *b_hoe;

   pt=0;
   eta=0;
   phi=0;
   r9=0;
   sieie=0;
   etaWidth=0;
   phiWidth=0;
   isMatched=0;
   Etrue=0;
   Eraw=0;
   nvtx=0;
   rho=0;
   hoe=0;

   tree->SetBranchAddress("genPt", &pt, &b_pt);
   tree->SetBranchAddress("genEta", &eta, &b_eta);  
   tree->SetBranchAddress("genPhi", &phi, &b_phi); 
   tree->SetBranchAddress("scSeedR9", &r9, &b_r9); 
   tree->SetBranchAddress("scSeedSigmaIetaIeta", &sieie, &b_sieie); 
   tree->SetBranchAddress("scPhiWidth", &phiWidth, &b_PhiWidth); 
   tree->SetBranchAddress("scEtaWidth", &etaWidth, &b_EtaWidth); 
   tree->SetBranchAddress("isMatched", &isMatched, &b_isMatched);
   tree->SetBranchAddress("genEnergy", &Etrue, &b_Etrue);  
   tree->SetBranchAddress("scRawEnergy", &Eraw, &b_Eraw);
   tree->SetBranchAddress("rho", &rho, &b_rho);    
   tree->SetBranchAddress("nVtx", &nvtx, &b_nvtx);    
   tree->SetBranchAddress("hoe", &hoe, &b_hoe);    

   TH1F *h_efficiency_Eta_N = new TH1F("h_efficiency_Eta_N","h_efficiency_Eta_N",300,-3,3);
   TH1F *h_efficiency_Eta_D = new TH1F("h_efficiency_Eta_D","h_efficiency_Eta_D",300,-3,3);
   TH1F *h_efficiency_Eta = new TH1F("h_efficiency_Eta","",300,-3,3);

   TH1F *h_pt = new TH1F("gen Pt","",300,0,300);
   TH1F *h_eta = new TH1F("Eta","",300,-3,3);
   TH1F *h_phi = new TH1F("Phi","",314,-3.14,3.14);
   TH1F *h_r9_EB = new TH1F("R9 (EB)","",440,0,1.1);
   TH1F *h_r9_EE = new TH1F("R9 (EE)","",440,0,1.1);
   TH1F *h_sieie_EB = new TH1F("SigIeIe (EB)","",200,0,0.02);
   TH1F *h_sieie_EE = new TH1F("SigIeIe (EE)","",200,0,0.04);
   TH1F *h_PhiWidth_EB = new TH1F("Phi Width (EB)","",200,0,0.1);
   TH1F *h_PhiWidth_EE = new TH1F("Phi Width (EE)","",200,0,0.1);
   TH1F *h_EtaWidth_EB = new TH1F("Eta Width (EB)","",500,0,0.05);
   TH1F *h_EtaWidth_EE = new TH1F("Eta Width (EE)","",500,0,0.05);
   TH1F *h_Energy_EB = new TH1F("Eraw/Etrue (EB)","",400,0,2);
   TH1F *h_Energy_EE = new TH1F("Eraw/Etrue (EE)","",400,0,2);
   TH1F *h_hoe_EB = new TH1F("H/E (EB)","",500,0,0.5);
   TH1F *h_hoe_EE = new TH1F("H/E (EE)","",500,0,0.5);
   TH1F *h_rho = new TH1F("rho","",100,0,100);
   TH1F *h_nvtx = new TH1F("nvtx","",100,0,100);

for (Long64_t jentry=0; jentry<tree->GetEntries();jentry++) {
   if (jentry%100000 == 0) cout << "Event number "  << jentry << endl ;
   Long64_t ientry =  tree->LoadTree(jentry);
   //Read the branches
   b_pt->GetEntry(ientry); 
   b_eta->GetEntry(ientry); 
   b_phi->GetEntry(ientry); 
   b_r9->GetEntry(ientry); 
   b_sieie->GetEntry(ientry); 
   b_PhiWidth->GetEntry(ientry);
   b_EtaWidth->GetEntry(ientry);  
   b_isMatched->GetEntry(ientry);
   b_Eraw->GetEntry(ientry);
   b_Etrue->GetEntry(ientry);
   b_hoe->GetEntry(ientry); 
   b_rho->GetEntry(ientry);
   b_nvtx->GetEntry(ientry);     


   if(pt>10){

   	h_efficiency_Eta_D->Fill(eta);

   	if(isMatched==1){
		h_efficiency_Eta_N->Fill(eta);
		h_pt->Fill(pt);
		h_eta->Fill(eta);
		h_phi->Fill(phi);
		h_rho->Fill(rho);
		h_nvtx->Fill(nvtx);
		if(eta<1.48){
			h_r9_EB->Fill(r9);
			h_sieie_EB->Fill(sieie);
			h_PhiWidth_EB->Fill(phiWidth);
			h_EtaWidth_EB->Fill(etaWidth);
			h_hoe_EB->Fill(hoe);
			h_Energy_EB->Fill(Eraw/Etrue);
		}
		if(eta>1.52){
			h_r9_EE->Fill(r9);
			h_sieie_EE->Fill(sieie);
			h_PhiWidth_EE->Fill(phiWidth);
			h_EtaWidth_EE->Fill(etaWidth);
			h_hoe_EE->Fill(hoe);
			h_Energy_EE->Fill(Eraw/Etrue);
		}
	}


   }

}

   h_efficiency_Eta->Add(h_efficiency_Eta_N);
   h_efficiency_Eta->Divide(h_efficiency_Eta_D);

   for (int i=0;i<300;i++){
   if(h_efficiency_Eta->GetBinContent(i)!=0) h_efficiency_Eta->SetBinError(i,sqrt(h_efficiency_Eta->GetBinContent(i)*(1-h_efficiency_Eta->GetBinContent(i))/  h_efficiency_Eta_N->GetBinContent(i)));
   }

   TFile *file = new TFile(Form("PhotonVariables%i.root",sc),"RECREATE"); 

  	
   h_efficiency_Eta->Write("h_Efficiency_Eta");
   h_pt->Write("h_Pt");
   h_eta->Write("h_Eta");
   h_phi->Write("h_Phi");
   h_rho->Write("h_Rho");
   h_nvtx->Write("h_Nvtx");
   h_r9_EB->Write("h_R9_EB");
   h_sieie_EB->Write("h_Sieie_EB");
   h_PhiWidth_EB->Write("h_PhiWidth_EB");
   h_EtaWidth_EB->Write("h_EtaWidth_EB");
   h_hoe_EB->Write("h_Hoe_EB");
   h_Energy_EB->Write("h_Energy_EB");
   h_r9_EE->Write("h_R9_EE");
   h_sieie_EE->Write("h_Sieie_EE");
   h_PhiWidth_EE->Write("h_PhiWidth_EE");
   h_EtaWidth_EE->Write("h_EtaWidth_EE");
   h_hoe_EE->Write("h_Hoe_EE");
   h_Energy_EE->Write("h_Energy_EE");

   file->Close();

}

void DumpVariables(){

   DumpVariablesPerScenario(1);
   DumpVariablesPerScenario(2);
   DumpVariablesPerScenario(3);

   TFile *f1 = TFile::Open("PhotonVariables1.root");
   TFile *f2 = TFile::Open("PhotonVariables2.root");
   TFile *f3 = TFile::Open("PhotonVariables3.root");

   //Plotting Reconstruction Efficency vs Eta

   TCanvas *cEffEta=new TCanvas;
   TH1 *hEffEta1=(TH1*)f1->Get("h_Efficiency_Eta");
   TH1 *hEffEta2=(TH1*)f2->Get("h_Efficiency_Eta");
   TH1 *hEffEta3=(TH1*)f3->Get("h_Efficiency_Eta");
   hEffEta1->GetYaxis()->SetRangeUser(0.6,1.1);
   hEffEta1->Draw("E");
   hEffEta2->Draw("E,same");
   hEffEta3->Draw("E,same");
   hEffEta1->SetLineColor(kBlue+2);
   hEffEta2->SetLineColor(kRed+1);
   hEffEta3->SetLineColor(kOrange+1);
   hEffEta1->SetLineWidth(2);
   hEffEta2->SetLineWidth(2);
   hEffEta3->SetLineWidth(2);
   hEffEta1->GetXaxis()->SetTitle("#eta");
   hEffEta1->GetYaxis()->SetTitle("reconstruction efficiency");
   hEffEta1->SetStats(kFALSE);
   TLegend *leg = new TLegend(0.73,0.15,0.88,0.35);
   leg->AddEntry(hEffEta1,"no PU","l");
   leg->AddEntry(hEffEta2,"bx 50","l");
   leg->AddEntry(hEffEta3,"bx 25","l");
   leg->Draw();
   cEffEta->SaveAs("plotsVariables/Efficiency_Eta.png");
   cEffEta->SaveAs("plotsVariables/Efficiency_Eta.pdf");
   cEffEta->SaveAs("plotsVariables/Efficiency_Eta.root");

   TCanvas *cPt=new TCanvas;
   TH1 *hPt1=(TH1*)f1->Get("h_Pt");
   TH1 *hPt2=(TH1*)f2->Get("h_Pt");
   TH1 *hPt3=(TH1*)f3->Get("h_Pt");
   hPt1->Scale(1./hPt1->Integral());
   hPt2->Scale(1./hPt2->Integral());
   hPt3->Scale(1./hPt3->Integral());
   hPt1->Draw();
   hPt2->Draw("same");
   hPt3->Draw("same");
   hPt1->SetLineColor(kBlue+2);
   hPt2->SetLineColor(kRed+1);
   hPt3->SetLineColor(kOrange+1);
   hPt1->SetLineWidth(2);
   hPt2->SetLineWidth(2);
   hPt3->SetLineWidth(2);
   hPt1->GetXaxis()->SetTitle("p_T");
   hPt1->SetStats(kFALSE);
   leg->Draw();
   cPt->SaveAs("plotsVariables/Pt.png");
   cPt->SaveAs("plotsVariables/Pt.pdf");
   cPt->SaveAs("plotsVariables/Pt.root");

   TCanvas *cEta=new TCanvas;
   TH1 *hEta1=(TH1*)f1->Get("h_Eta");
   TH1 *hEta2=(TH1*)f2->Get("h_Eta");
   TH1 *hEta3=(TH1*)f3->Get("h_Eta");
   hEta1->Scale(1./hEta1->Integral());
   hEta2->Scale(1./hEta2->Integral());
   hEta3->Scale(1./hEta3->Integral());
   hEta1->Draw();
   hEta2->Draw("same");
   hEta3->Draw("same");
   hEta1->SetLineColor(kBlue+2);
   hEta2->SetLineColor(kRed+1);
   hEta3->SetLineColor(kOrange+1);
   hEta1->SetLineWidth(2);
   hEta2->SetLineWidth(2);
   hEta3->SetLineWidth(2);
   hEta1->GetXaxis()->SetTitle("#eta");
   hEta1->SetStats(kFALSE);
   leg->Draw();
   cEta->SaveAs("plotsVariables/Eta.png");
   cEta->SaveAs("plotsVariables/Eta.pdf");
   cEta->SaveAs("plotsVariables/Eta.root");

   TCanvas *cPhi=new TCanvas;
   TH1 *hPhi1=(TH1*)f1->Get("h_Phi");
   TH1 *hPhi2=(TH1*)f2->Get("h_Phi");
   TH1 *hPhi3=(TH1*)f3->Get("h_Phi");
   hPhi1->Scale(1./hPhi1->Integral());
   hPhi2->Scale(1./hPhi2->Integral());
   hPhi3->Scale(1./hPhi3->Integral());
   hPhi1->Draw();
   hPhi2->Draw("same");
   hPhi3->Draw("same");
   hPhi1->SetLineColor(kBlue+2);
   hPhi2->SetLineColor(kRed+1);
   hPhi3->SetLineColor(kOrange+1);
   hPhi1->SetLineWidth(2);
   hPhi2->SetLineWidth(2);
   hPhi3->SetLineWidth(2);
   hPhi1->GetXaxis()->SetTitle("#phi");
   hPhi1->SetStats(kFALSE);
   leg->Draw();
   cPhi->SaveAs("plotsVariables/Phi.png");
   cPhi->SaveAs("plotsVariables/Phi.pdf");
   cPhi->SaveAs("plotsVariables/Phi.root");

   TCanvas *cRho=new TCanvas;
   TH1 *hRho1=(TH1*)f1->Get("h_Rho");
   TH1 *hRho2=(TH1*)f2->Get("h_Rho");
   TH1 *hRho3=(TH1*)f3->Get("h_Rho");
   hRho1->Scale(1./hRho1->Integral());
   hRho2->Scale(1./hRho2->Integral());
   hRho3->Scale(1./hRho3->Integral());
   hRho1->Draw();
   hRho2->Draw("same");
   hRho3->Draw("same");
   hRho1->SetLineColor(kBlue+2);
   hRho2->SetLineColor(kRed+1);
   hRho3->SetLineColor(kOrange+1);
   hRho1->SetLineWidth(2);
   hRho2->SetLineWidth(2);
   hRho3->SetLineWidth(2);
   hRho1->GetXaxis()->SetTitle("#rho");
   hRho1->SetStats(kFALSE);
   leg->Draw();
   cRho->SaveAs("plotsVariables/Rho.png");
   cRho->SaveAs("plotsVariables/Rho.pdf");
   cRho->SaveAs("plotsVariables/Rho.root");

   TCanvas *cNvtx=new TCanvas;
   cNvtx->SetLogy();
   TH1 *hNvtx1=(TH1*)f1->Get("h_Nvtx");
   TH1 *hNvtx2=(TH1*)f2->Get("h_Nvtx");
   TH1 *hNvtx3=(TH1*)f3->Get("h_Nvtx");
   hNvtx1->Scale(1./hNvtx1->Integral());
   hNvtx2->Scale(1./hNvtx2->Integral());
   hNvtx3->Scale(1./hNvtx3->Integral());
   hNvtx1->GetXaxis()->SetRangeUser(0,60);
   hNvtx1->Draw();
   hNvtx2->Draw("same");
   hNvtx3->Draw("same");

   hNvtx1->SetLineColor(kBlue+2);
   hNvtx2->SetLineColor(kRed+1);
   hNvtx3->SetLineColor(kOrange+1);
   hNvtx1->SetLineWidth(2);
   hNvtx2->SetLineWidth(2);
   hNvtx3->SetLineWidth(2);
   hNvtx1->GetXaxis()->SetTitle("N_Vtx");
   hNvtx1->SetStats(kFALSE);
   leg->Draw();
   cNvtx->SaveAs("plotsVariables/Nvtx.png");
   cNvtx->SaveAs("plotsVariables/Nvtx.pdf");
   cNvtx->SaveAs("plotsVariables/Nvtx.root");

   TCanvas *cR9_EB=new TCanvas;
   TH1 *hR9_EB1=(TH1*)f1->Get("h_R9_EB");
   TH1 *hR9_EB2=(TH1*)f2->Get("h_R9_EB");
   TH1 *hR9_EB3=(TH1*)f3->Get("h_R9_EB");
   hR9_EB1->Scale(1./hR9_EB1->Integral());
   hR9_EB2->Scale(1./hR9_EB2->Integral());
   hR9_EB3->Scale(1./hR9_EB3->Integral());
   hR9_EB1->Draw();
   hR9_EB2->Draw("same");
   hR9_EB3->Draw("same");
   hR9_EB1->SetLineColor(kBlue+2);
   hR9_EB2->SetLineColor(kRed+1);
   hR9_EB3->SetLineColor(kOrange+1);
   hR9_EB1->SetLineWidth(2);
   hR9_EB2->SetLineWidth(2);
   hR9_EB3->SetLineWidth(2);
   hR9_EB1->GetXaxis()->SetTitle("r9");
   hR9_EB1->SetStats(kFALSE);
   TLegend *leg2 = new TLegend(0.13,0.15,0.28,0.35);
   leg2->AddEntry(hEffEta1,"no PU","l");
   leg2->AddEntry(hEffEta2,"bx 50","l");
   leg2->AddEntry(hEffEta3,"bx 25","l");
   leg2->Draw();
   cR9_EB->SaveAs("plotsVariables/R9_EB.png");
   cR9_EB->SaveAs("plotsVariables/R9_EB.pdf");
   cR9_EB->SaveAs("plotsVariables/R9_EB.root");

   TCanvas *cSieie_EB=new TCanvas;
   TH1 *hSieie_EB1=(TH1*)f1->Get("h_Sieie_EB");
   TH1 *hSieie_EB2=(TH1*)f2->Get("h_Sieie_EB");
   TH1 *hSieie_EB3=(TH1*)f3->Get("h_Sieie_EB");
   hSieie_EB1->Scale(1./hSieie_EB1->Integral());
   hSieie_EB2->Scale(1./hSieie_EB2->Integral());
   hSieie_EB3->Scale(1./hSieie_EB3->Integral());
   hSieie_EB1->Draw();
   hSieie_EB2->Draw("same");
   hSieie_EB3->Draw("same");
   hSieie_EB1->SetLineColor(kBlue+2);
   hSieie_EB2->SetLineColor(kRed+1);
   hSieie_EB3->SetLineColor(kOrange+1);
   hSieie_EB1->SetLineWidth(2);
   hSieie_EB2->SetLineWidth(2);
   hSieie_EB3->SetLineWidth(2);
   hSieie_EB1->GetXaxis()->SetTitle("#sigma_{i#eta i#eta}");
   hSieie_EB1->SetStats(kFALSE);
   leg->Draw();
   cSieie_EB->SaveAs("plotsVariables/Sieie_EB.png");
   cSieie_EB->SaveAs("plotsVariables/Sieie_EB.pdf");
   cSieie_EB->SaveAs("plotsVariables/Sieie_EB.root");

   TCanvas *cPhiWidth_EB=new TCanvas;
   TH1 *hPhiWidth_EB1=(TH1*)f1->Get("h_PhiWidth_EB");
   TH1 *hPhiWidth_EB2=(TH1*)f2->Get("h_PhiWidth_EB");
   TH1 *hPhiWidth_EB3=(TH1*)f3->Get("h_PhiWidth_EB");
   hPhiWidth_EB1->Scale(1./hPhiWidth_EB1->Integral());
   hPhiWidth_EB2->Scale(1./hPhiWidth_EB2->Integral());
   hPhiWidth_EB3->Scale(1./hPhiWidth_EB3->Integral());
   hPhiWidth_EB1->Draw();
   hPhiWidth_EB2->Draw("same");
   hPhiWidth_EB3->Draw("same");
   hPhiWidth_EB1->SetLineColor(kBlue+2);
   hPhiWidth_EB2->SetLineColor(kRed+1);
   hPhiWidth_EB3->SetLineColor(kOrange+1);
   hPhiWidth_EB1->SetLineWidth(2);
   hPhiWidth_EB2->SetLineWidth(2);
   hPhiWidth_EB3->SetLineWidth(2);
   hPhiWidth_EB1->GetXaxis()->SetTitle("#phi width");
   hPhiWidth_EB1->SetStats(kFALSE);
   leg->Draw();
   cPhiWidth_EB->SaveAs("plotsVariables/PhiWidth_EB.png");
   cPhiWidth_EB->SaveAs("plotsVariables/PhiWidth_EB.pdf");
   cPhiWidth_EB->SaveAs("plotsVariables/PhiWidth_EB.root");

   TCanvas *cEtaWidth_EB=new TCanvas;
   TH1 *hEtaWidth_EB1=(TH1*)f1->Get("h_EtaWidth_EB");
   TH1 *hEtaWidth_EB2=(TH1*)f2->Get("h_EtaWidth_EB");
   TH1 *hEtaWidth_EB3=(TH1*)f3->Get("h_EtaWidth_EB");
   hEtaWidth_EB1->Scale(1./hEtaWidth_EB1->Integral());
   hEtaWidth_EB2->Scale(1./hEtaWidth_EB2->Integral());
   hEtaWidth_EB3->Scale(1./hEtaWidth_EB3->Integral());
   hEtaWidth_EB1->Draw();
   hEtaWidth_EB2->Draw("same");
   hEtaWidth_EB3->Draw("same");
   hEtaWidth_EB1->SetLineColor(kBlue+2);
   hEtaWidth_EB2->SetLineColor(kRed+1);
   hEtaWidth_EB3->SetLineColor(kOrange+1);
   hEtaWidth_EB1->SetLineWidth(2);
   hEtaWidth_EB2->SetLineWidth(2);
   hEtaWidth_EB3->SetLineWidth(2);
   hEtaWidth_EB1->GetXaxis()->SetTitle("#eta width");
   hEtaWidth_EB1->SetStats(kFALSE);
   leg->Draw();
   cEtaWidth_EB->SaveAs("plotsVariables/EtaWidth_EB.png");
   cEtaWidth_EB->SaveAs("plotsVariables/EtaWidth_EB.pdf");
   cEtaWidth_EB->SaveAs("plotsVariables/EtaWidth_EB.root");

   TCanvas *cHoe_EB=new TCanvas;
   cHoe_EB->SetLogy();
   TH1 *hHoe_EB1=(TH1*)f1->Get("h_Hoe_EB");
   TH1 *hHoe_EB2=(TH1*)f2->Get("h_Hoe_EB");
   TH1 *hHoe_EB3=(TH1*)f3->Get("h_Hoe_EB");
   hHoe_EB1->Scale(1./hHoe_EB1->Integral());
   hHoe_EB2->Scale(1./hHoe_EB2->Integral());
   hHoe_EB3->Scale(1./hHoe_EB3->Integral());
   hHoe_EB1->Draw();
   hHoe_EB2->Draw("same");
   hHoe_EB3->Draw("same");
   hHoe_EB1->SetLineColor(kBlue+2);
   hHoe_EB2->SetLineColor(kRed+1);
   hHoe_EB3->SetLineColor(kOrange+1);
   hHoe_EB1->SetLineWidth(2);
   hHoe_EB2->SetLineWidth(2);
   hHoe_EB3->SetLineWidth(2);
   hHoe_EB1->GetXaxis()->SetTitle("H/E");
   hHoe_EB1->SetStats(kFALSE);
   leg->Draw();
   cHoe_EB->SaveAs("plotsVariables/Hoe_EB.png");
   cHoe_EB->SaveAs("plotsVariables/Hoe_EB.pdf");
   cHoe_EB->SaveAs("plotsVariables/Hoe_EB.root");

   TCanvas *cEnergy_EB=new TCanvas;
   TH1 *hEnergy_EB1=(TH1*)f1->Get("h_Energy_EB");
   TH1 *hEnergy_EB2=(TH1*)f2->Get("h_Energy_EB");
   TH1 *hEnergy_EB3=(TH1*)f3->Get("h_Energy_EB");
   hEnergy_EB1->Scale(1./hEnergy_EB1->Integral());
   hEnergy_EB2->Scale(1./hEnergy_EB2->Integral());
   hEnergy_EB3->Scale(1./hEnergy_EB3->Integral());
   hEnergy_EB1->Draw();
   hEnergy_EB2->Draw("same");
   hEnergy_EB3->Draw("same");
   hEnergy_EB1->SetLineColor(kBlue+2);
   hEnergy_EB2->SetLineColor(kRed+1);
   hEnergy_EB3->SetLineColor(kOrange+1);
   hEnergy_EB1->SetLineWidth(2);
   hEnergy_EB2->SetLineWidth(2);
   hEnergy_EB3->SetLineWidth(2);
   hEnergy_EB1->GetXaxis()->SetTitle("E_raw / E_true");
   hEnergy_EB1->SetStats(kFALSE);
   leg->Draw();
   cEnergy_EB->SaveAs("plotsVariables/Energy_EB.png");
   cEnergy_EB->SaveAs("plotsVariables/Energy_EB.pdf");
   cEnergy_EB->SaveAs("plotsVariables/Energy_EB.root");

   TCanvas *cEnergy_EB_log=new TCanvas;
   cEnergy_EB_log->SetLogy();
   TH1 *hEnergy_EB_log1=(TH1*)f1->Get("h_Energy_EB");
   TH1 *hEnergy_EB_log2=(TH1*)f2->Get("h_Energy_EB");
   TH1 *hEnergy_EB_log3=(TH1*)f3->Get("h_Energy_EB");
   hEnergy_EB_log1->Scale(1./hEnergy_EB_log1->Integral());
   hEnergy_EB_log2->Scale(1./hEnergy_EB_log2->Integral());
   hEnergy_EB_log3->Scale(1./hEnergy_EB_log3->Integral());
   hEnergy_EB_log1->Draw();
   hEnergy_EB_log2->Draw("same");
   hEnergy_EB_log3->Draw("same");
   hEnergy_EB_log1->SetLineColor(kBlue+2);
   hEnergy_EB_log2->SetLineColor(kRed+1);
   hEnergy_EB_log3->SetLineColor(kOrange+1);
   hEnergy_EB_log1->SetLineWidth(2);
   hEnergy_EB_log2->SetLineWidth(2);
   hEnergy_EB_log3->SetLineWidth(2);
   hEnergy_EB_log1->GetXaxis()->SetTitle("E_raw / E_true");
   hEnergy_EB_log1->SetStats(kFALSE);
   leg->Draw();
   cEnergy_EB_log->SaveAs("plotsVariables/Energy_EB_log.png");
   cEnergy_EB_log->SaveAs("plotsVariables/Energy_EB_log.pdf");
   cEnergy_EB_log->SaveAs("plotsVariables/Energy_EB_log.root");

   TCanvas *cR9_EE=new TCanvas;
   TH1 *hR9_EE1=(TH1*)f1->Get("h_R9_EE");
   TH1 *hR9_EE2=(TH1*)f2->Get("h_R9_EE");
   TH1 *hR9_EE3=(TH1*)f3->Get("h_R9_EE");
   hR9_EE1->Scale(1./hR9_EE1->Integral());
   hR9_EE2->Scale(1./hR9_EE2->Integral());
   hR9_EE3->Scale(1./hR9_EE3->Integral());
   hR9_EE1->Draw();
   hR9_EE2->Draw("same");
   hR9_EE3->Draw("same");
   hR9_EE1->SetLineColor(kBlue+2);
   hR9_EE2->SetLineColor(kRed+1);
   hR9_EE3->SetLineColor(kOrange+1);
   hR9_EE1->SetLineWidth(2);
   hR9_EE2->SetLineWidth(2);
   hR9_EE3->SetLineWidth(2);
   hR9_EE1->GetXaxis()->SetTitle("r9");
   hR9_EE1->SetStats(kFALSE);
   leg2->Draw();
   cR9_EE->SaveAs("plotsVariables/R9_EE.png");
   cR9_EE->SaveAs("plotsVariables/R9_EE.pdf");
   cR9_EE->SaveAs("plotsVariables/R9_EE.root");

   TCanvas *cSieie_EE=new TCanvas;
   TH1 *hSieie_EE1=(TH1*)f1->Get("h_Sieie_EE");
   TH1 *hSieie_EE2=(TH1*)f2->Get("h_Sieie_EE");
   TH1 *hSieie_EE3=(TH1*)f3->Get("h_Sieie_EE");
   hSieie_EE1->Scale(1./hSieie_EE1->Integral());
   hSieie_EE2->Scale(1./hSieie_EE2->Integral());
   hSieie_EE3->Scale(1./hSieie_EE3->Integral());
   hSieie_EE1->Draw();
   hSieie_EE2->Draw("same");
   hSieie_EE3->Draw("same");
   hSieie_EE1->SetLineColor(kBlue+2);
   hSieie_EE2->SetLineColor(kRed+1);
   hSieie_EE3->SetLineColor(kOrange+1);
   hSieie_EE1->SetLineWidth(2);
   hSieie_EE2->SetLineWidth(2);
   hSieie_EE3->SetLineWidth(2);
   hSieie_EE1->GetXaxis()->SetTitle("#sigma_{i#eta i#eta}");
   hSieie_EE1->SetStats(kFALSE);
   leg->Draw();
   cSieie_EE->SaveAs("plotsVariables/Sieie_EE.png");
   cSieie_EE->SaveAs("plotsVariables/Sieie_EE.pdf");
   cSieie_EE->SaveAs("plotsVariables/Sieie_EE.root");

   TCanvas *cPhiWidth_EE=new TCanvas;
   TH1 *hPhiWidth_EE1=(TH1*)f1->Get("h_PhiWidth_EE");
   TH1 *hPhiWidth_EE2=(TH1*)f2->Get("h_PhiWidth_EE");
   TH1 *hPhiWidth_EE3=(TH1*)f3->Get("h_PhiWidth_EE");
   hPhiWidth_EE1->Scale(1./hPhiWidth_EE1->Integral());
   hPhiWidth_EE2->Scale(1./hPhiWidth_EE2->Integral());
   hPhiWidth_EE3->Scale(1./hPhiWidth_EE3->Integral());
   hPhiWidth_EE1->Draw();
   hPhiWidth_EE2->Draw("same");
   hPhiWidth_EE3->Draw("same");
   hPhiWidth_EE1->SetLineColor(kBlue+2);
   hPhiWidth_EE2->SetLineColor(kRed+1);
   hPhiWidth_EE3->SetLineColor(kOrange+1);
   hPhiWidth_EE1->SetLineWidth(2);
   hPhiWidth_EE2->SetLineWidth(2);
   hPhiWidth_EE3->SetLineWidth(2);
   hPhiWidth_EE1->GetXaxis()->SetTitle("#phi width");
   hPhiWidth_EE1->SetStats(kFALSE);
   leg->Draw();
   cPhiWidth_EE->SaveAs("plotsVariables/PhiWidth_EE.png");
   cPhiWidth_EE->SaveAs("plotsVariables/PhiWidth_EE.pdf");
   cPhiWidth_EE->SaveAs("plotsVariables/PhiWidth_EE.root");

   TCanvas *cEtaWidth_EE=new TCanvas;
   TH1 *hEtaWidth_EE1=(TH1*)f1->Get("h_EtaWidth_EE");
   TH1 *hEtaWidth_EE2=(TH1*)f2->Get("h_EtaWidth_EE");
   TH1 *hEtaWidth_EE3=(TH1*)f3->Get("h_EtaWidth_EE");
   hEtaWidth_EE1->Scale(1./hEtaWidth_EE1->Integral());
   hEtaWidth_EE2->Scale(1./hEtaWidth_EE2->Integral());
   hEtaWidth_EE3->Scale(1./hEtaWidth_EE3->Integral());
   hEtaWidth_EE1->Draw();
   hEtaWidth_EE2->Draw("same");
   hEtaWidth_EE3->Draw("same");
   hEtaWidth_EE1->SetLineColor(kBlue+2);
   hEtaWidth_EE2->SetLineColor(kRed+1);
   hEtaWidth_EE3->SetLineColor(kOrange+1);
   hEtaWidth_EE1->SetLineWidth(2);
   hEtaWidth_EE2->SetLineWidth(2);
   hEtaWidth_EE3->SetLineWidth(2);
   hEtaWidth_EE1->GetXaxis()->SetTitle("#eta width");
   hEtaWidth_EE1->SetStats(kFALSE);
   leg->Draw();
   cEtaWidth_EE->SaveAs("plotsVariables/EtaWidth_EE.png");
   cEtaWidth_EE->SaveAs("plotsVariables/EtaWidth_EE.pdf");
   cEtaWidth_EE->SaveAs("plotsVariables/EtaWidth_EE.root");

   TCanvas *cHoe_EE=new TCanvas;
   cHoe_EE->SetLogy();
   TH1 *hHoe_EE1=(TH1*)f1->Get("h_Hoe_EE");
   TH1 *hHoe_EE2=(TH1*)f2->Get("h_Hoe_EE");
   TH1 *hHoe_EE3=(TH1*)f3->Get("h_Hoe_EE");
   hHoe_EE1->Scale(1./hHoe_EE1->Integral());
   hHoe_EE2->Scale(1./hHoe_EE2->Integral());
   hHoe_EE3->Scale(1./hHoe_EE3->Integral());
   hHoe_EE1->Draw();
   hHoe_EE2->Draw("same");
   hHoe_EE3->Draw("same");
   hHoe_EE1->SetLineColor(kBlue+2);
   hHoe_EE2->SetLineColor(kRed+1);
   hHoe_EE3->SetLineColor(kOrange+1);
   hHoe_EE1->SetLineWidth(2);
   hHoe_EE2->SetLineWidth(2);
   hHoe_EE3->SetLineWidth(2);
   hHoe_EE1->GetXaxis()->SetTitle("H/E");
   hHoe_EE1->SetStats(kFALSE);
   leg->Draw();
   cHoe_EE->SaveAs("plotsVariables/Hoe_EE.png");
   cHoe_EE->SaveAs("plotsVariables/Hoe_EE.pdf");
   cHoe_EE->SaveAs("plotsVariables/Hoe_EE.root");

   TCanvas *cEnergy_EE=new TCanvas;
   TH1 *hEnergy_EE1=(TH1*)f1->Get("h_Energy_EE");
   TH1 *hEnergy_EE2=(TH1*)f2->Get("h_Energy_EE");
   TH1 *hEnergy_EE3=(TH1*)f3->Get("h_Energy_EE");
   hEnergy_EE1->Scale(1./hEnergy_EE1->Integral());
   hEnergy_EE2->Scale(1./hEnergy_EE2->Integral());
   hEnergy_EE3->Scale(1./hEnergy_EE3->Integral());
   hEnergy_EE1->Draw();
   hEnergy_EE2->Draw("same");
   hEnergy_EE3->Draw("same");
   hEnergy_EE1->SetLineColor(kBlue+2);
   hEnergy_EE2->SetLineColor(kRed+1);
   hEnergy_EE3->SetLineColor(kOrange+1);
   hEnergy_EE1->SetLineWidth(2);
   hEnergy_EE2->SetLineWidth(2);
   hEnergy_EE3->SetLineWidth(2);
   hEnergy_EE1->GetXaxis()->SetTitle("E_raw / E_true");
   hEnergy_EE1->SetStats(kFALSE);
   leg->Draw();
   cEnergy_EE->SaveAs("plotsVariables/Energy_EE.png");
   cEnergy_EE->SaveAs("plotsVariables/Energy_EE.pdf");
   cEnergy_EE->SaveAs("plotsVariables/Energy_EE.root");

   TCanvas *cEnergy_EE_log=new TCanvas;
   cEnergy_EE_log->SetLogy();
   TH1 *hEnergy_EE_log1=(TH1*)f1->Get("h_Energy_EE");
   TH1 *hEnergy_EE_log2=(TH1*)f2->Get("h_Energy_EE");
   TH1 *hEnergy_EE_log3=(TH1*)f3->Get("h_Energy_EE");
   hEnergy_EE_log1->Scale(1./hEnergy_EE_log1->Integral());
   hEnergy_EE_log2->Scale(1./hEnergy_EE_log2->Integral());
   hEnergy_EE_log3->Scale(1./hEnergy_EE_log3->Integral());
   hEnergy_EE_log1->Draw();
   hEnergy_EE_log2->Draw("same");
   hEnergy_EE_log3->Draw("same");
   hEnergy_EE_log1->SetLineColor(kBlue+2);
   hEnergy_EE_log2->SetLineColor(kRed+1);
   hEnergy_EE_log3->SetLineColor(kOrange+1);
   hEnergy_EE_log1->SetLineWidth(2);
   hEnergy_EE_log2->SetLineWidth(2);
   hEnergy_EE_log3->SetLineWidth(2);
   hEnergy_EE_log1->GetXaxis()->SetTitle("E_raw / E_true");
   hEnergy_EE_log1->SetStats(kFALSE);
   leg->Draw();
   cEnergy_EE_log->SaveAs("plotsVariables/Energy_EE_log.png");
   cEnergy_EE_log->SaveAs("plotsVariables/Energy_EE_log.pdf");
   cEnergy_EE_log->SaveAs("plotsVariables/Energy_EE_log.root");

}



