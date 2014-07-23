#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TBranch.h"
#include "math.h"
#include <cstdlib>

void RegVars(int, bool, char*);

int main(int argc, char *argv[]){
  int arg1=300, arg2=1;
  char *altFile; sprintf(altFile, "");
  if(argc>1) arg1=atoi(argv[1]);
  if(argc>2) arg2=atoi(argv[2]);
  if(argc>3) altFile=argv[3];
  if(argc==1) printf("RegVars [mass] [splitFile] [altFile=\"\"]\n");
  else RegVars(arg1,arg2,altFile);
  return 0;
}

void RegVars(int mass, bool splitFile, char *altFile){
  char fileNameBase[80];
  if(altFile[0]=='\0') sprintf(fileNameBase,"files/tree_Signal.root");
  else sprintf(fileNameBase,"files/%s.root",altFile);
  TFile *inFile = new TFile(TString::Format("%s",fileNameBase),"READ");
  TTree *inTree = (TTree*)inFile->Get(TString::Format("Radion_m%i_8TeV",mass));
  //TTree *inTree = (TTree*)inFile->Get(TString::Format("MSSM_H_m%i_8TeV",mass));
  //TTree *inTree = (TTree*)inFile->Get(TString::Format("ggHH_8TeV"));
  char inTreeName[80];
  sprintf(inTreeName,"%s",inTree->GetName());
  TFile *outFile = new TFile (TString::Format("files/%s_regInput.root",inTreeName),"RECREATE");
  TTree *outTree = new TTree("EventTree","Regression data");

  const int maxP = 15;
  int event;
  float rho25, MET, METphi;
  float jet_pt[maxP], jet_eta[maxP], jet_phi[maxP], jet_e[maxP], jet_secVtxPt[maxP], jet_secVtxM[maxP], jet_secVtx3dL[maxP], jet_secVtx3deL[maxP], jet_softLeptPt[maxP], jet_softLeptPtRel[maxP], jet_softLeptDR[maxP], jet_emfrac[maxP], jet_hadfrac[maxP], jet_chadfrac[maxP], jet_nhadfrac[maxP], jet_phofrac[maxP], jet_mufrac[maxP], jet_elefrac[maxP], jet_JECUnc[maxP], jet_leadTrackPt[maxP], jet_genPt[maxP];
  int jet_bgenMatched[maxP], jet_flavour[maxP], jet_softLeptIdLooseMu[maxP], jet_softLeptIdEle95[maxP], jet_nNeutrals[maxP], jet_nCharged[maxP];
  float bjet_pt, bjet_phi, bjet_eta, bjet_mt, bjet_e, bjet_secVtxPt, bjet_secVtxM, bjet_secVtx3dL, bjet_secVtx3deL, bjet_emfrac, bjet_hadfrac, bjet_chadfrac, bjet_nhadfrac, bjet_phofrac, bjet_mufrac, bjet_elefrac, bjet_JECUnc, bjet_leadTrackPt, bjet_softLeptPt, bjet_softLeptPtRel, bjet_softLeptDR, bjet_genPt;
  int bjet_bgenMatched, bjet_flavour, bjet_softLeptIdLooseMu, bjet_softLeptIdEle95, bjet_nNeutrals, bjet_nCharged;

  inTree->SetBranchAddress("event", &event);
  inTree->SetBranchAddress("rho", &rho25);
  inTree->SetBranchAddress("met_corr_pfmet", &MET);
  inTree->SetBranchAddress("met_corr_phi_pfmet", &METphi);

  for(int i=0; i<maxP; ++i){
    inTree->SetBranchAddress(TString::Format("j%i_pt",i+1),&(jet_pt[i]));
    inTree->SetBranchAddress(TString::Format("j%i_phi",i+1),&jet_phi[i]);
    inTree->SetBranchAddress(TString::Format("j%i_eta",i+1),&jet_eta[i]);
    inTree->SetBranchAddress(TString::Format("j%i_e",i+1),&jet_e[i]);
    inTree->SetBranchAddress(TString::Format("j%i_bgenMatched",i+1),&jet_bgenMatched[i]);
    inTree->SetBranchAddress(TString::Format("j%i_secVtxPt",i+1),&jet_secVtxPt[i]);
    inTree->SetBranchAddress(TString::Format("j%i_secVtxM",i+1),&jet_secVtxM[i]);
    inTree->SetBranchAddress(TString::Format("j%i_secVtx3dL",i+1),&jet_secVtx3dL[i]);
    inTree->SetBranchAddress(TString::Format("j%i_secVtx3deL",i+1),&jet_secVtx3deL[i]);
    inTree->SetBranchAddress(TString::Format("j%i_emfrac",i+1),&jet_emfrac[i]);
    inTree->SetBranchAddress(TString::Format("j%i_hadfrac",i+1),&jet_hadfrac[i]);
    inTree->SetBranchAddress(TString::Format("j%i_nNeutrals",i+1),&jet_nNeutrals[i]);
    inTree->SetBranchAddress(TString::Format("j%i_nCharged",i+1),&jet_nCharged[i]);
    inTree->SetBranchAddress(TString::Format("j%i_chadfrac",i+1),&jet_chadfrac[i]);
    inTree->SetBranchAddress(TString::Format("j%i_nhadfrac",i+1),&jet_nhadfrac[i]);
    inTree->SetBranchAddress(TString::Format("j%i_phofrac",i+1),&jet_phofrac[i]);
    inTree->SetBranchAddress(TString::Format("j%i_mufrac",i+1),&jet_mufrac[i]);
    inTree->SetBranchAddress(TString::Format("j%i_elefrac",i+1),&jet_elefrac[i]);
    inTree->SetBranchAddress(TString::Format("j%i_JECUnc",i+1),&jet_JECUnc[i]);
    inTree->SetBranchAddress(TString::Format("j%i_leadTrackPt",i+1),&jet_leadTrackPt[i]);
    inTree->SetBranchAddress(TString::Format("j%i_softLeptPt",i+1),&jet_softLeptPt[i]);
    inTree->SetBranchAddress(TString::Format("j%i_softLeptPtRel",i+1),&jet_softLeptPtRel[i]);
    inTree->SetBranchAddress(TString::Format("j%i_softLeptDR",i+1),&jet_softLeptDR[i]);
    inTree->SetBranchAddress(TString::Format("j%i_softLeptIdLooseMu",i+1),&jet_softLeptIdLooseMu[i]);
    inTree->SetBranchAddress(TString::Format("j%i_softLeptIdEle95",i+1),&jet_softLeptIdEle95[i]);
    inTree->SetBranchAddress(TString::Format("j%i_genPt",i+1),&jet_genPt[i]);
    inTree->SetBranchAddress(TString::Format("j%i_flavour",i+1),&jet_flavour[i]);
  }

  outTree->Branch("event", &event, "event/I");
  outTree->Branch("rho", &rho25, "rho/F");
  outTree->Branch("MET", &MET, "MET/F");
  outTree->Branch("METphi", &METphi, "METphi/F");
  outTree->Branch("bjet_pt",		    &bjet_pt,		    "bjet_pt/F");		    
  outTree->Branch("bjet_phi",		    &bjet_phi,		    "bjet_phi/F");		    
  outTree->Branch("bjet_eta",		    &bjet_eta,		    "bjet_eta/F");		    
  outTree->Branch("bjet_e",		    &bjet_e,		    "bjet_e/F");		    
  outTree->Branch("bjet_mt",		    &bjet_mt,		    "bjet_mt/F");		    
  outTree->Branch("bjet_bgenMatched",	    &bjet_bgenMatched,	    "bjet_bgenMatched/I");	    
  outTree->Branch("bjet_secVtxPt",	    &bjet_secVtxPt,	    "bjet_secVtxPt/F");	    
  outTree->Branch("bjet_secVtxM",	    &bjet_secVtxM,	    "bjet_secVtxM/F");	    
  outTree->Branch("bjet_secVtx3dL",	    &bjet_secVtx3dL,	    "bjet_secVtx3dL/F");	    
  outTree->Branch("bjet_secVtx3deL",	    &bjet_secVtx3deL,	    "bjet_secVtx3deL/F");	    
  outTree->Branch("bjet_emfrac",	    &bjet_emfrac,	    "bjet_emfrac/F");	    
  outTree->Branch("bjet_hadfrac",	    &bjet_hadfrac,	    "bjet_hadfrac/F");	    
  outTree->Branch("bjet_nNeutrals",	    &bjet_nNeutrals,	    "bjet_nNeutrals/I");	    
  outTree->Branch("bjet_nCharged",	    &bjet_nCharged,	    "bjet_nCharged/I");	    
  outTree->Branch("bjet_chadfrac",	    &bjet_chadfrac,	    "bjet_chadfrac/F");	    
  outTree->Branch("bjet_nhadfrac",	    &bjet_nhadfrac,	    "bjet_nhadfrac/F");	    
  outTree->Branch("bjet_phofrac",	    &bjet_phofrac,	    "bjet_phofrac/F");	    
  outTree->Branch("bjet_mufrac",	    &bjet_mufrac,	    "bjet_mufrac/F");	    
  outTree->Branch("bjet_elefrac",	    &bjet_elefrac,	    "bjet_elefrac/F");	    
  outTree->Branch("bjet_JECUnc",	    &bjet_JECUnc,	    "bjet_JECUnc/F");	    
  outTree->Branch("bjet_leadTrackPt",	    &bjet_leadTrackPt,	    "bjet_leadTrackPt/F");	    
  outTree->Branch("bjet_softLeptPt",	    &bjet_softLeptPt,	    "bjet_softLeptPt/F");	    
  outTree->Branch("bjet_softLeptPtRel",	    &bjet_softLeptPtRel,    "bjet_softLeptPtRel/F");	    	    
  outTree->Branch("bjet_softLeptDR",	    &bjet_softLeptDR,	    "bjet_softLeptDR/F");	    
  outTree->Branch("bjet_softLeptIdLooseMu", &bjet_softLeptIdLooseMu,"bjet_softLeptIdLooseMu/I");  
  outTree->Branch("bjet_softLeptIdEle95",   &bjet_softLeptIdEle95,  "bjet_softLeptIdEle95/I");     
  outTree->Branch("bjet_genPt",		    &bjet_genPt,	    "bjet_genPt/F");		    	    
  outTree->Branch("bjet_flavour",           &bjet_flavour,          "bjet_flavour/I");              
  
  TH1D *bPt = new TH1D("bPt","bPt",250,0,1000);
  unsigned nEntries = inTree->GetEntries();
  for(unsigned i=0; i<nEntries; ++i){
    inTree->GetEntry(i);//printf("%i %i %.0f\n",i,i%15,jet_bgenMatched[i%15]);
    for(unsigned j=0; j<maxP; ++j){
      if((jet_flavour[j]==5 || jet_flavour[j]==-5) && jet_bgenMatched[j]==1){
	bjet_pt 		    = jet_pt[j];		   
	bjet_phi		    = jet_phi[j];		   
	bjet_eta		    = jet_eta[j];		   
	bjet_e	    	            = jet_e[j];		   
	bjet_mt	    	            = pow(jet_e[j],2)-pow(jet_pt[j]*(1+sinh(jet_eta[j])),2) > 0 ? sqrt(pow(jet_e[j],2)-pow(jet_pt[j]*(1+sinh(jet_eta[j])),2)) : sqrt(pow(jet_pt[j]*(1+sinh(jet_eta[j])),2)-pow(jet_e[j],2)) ;
	bjet_bgenMatched	    = jet_bgenMatched[j];	   
	bjet_secVtxPt	   	    = jet_secVtxPt[j];	   
	bjet_secVtxM	   	    = jet_secVtxM[j];	   
	bjet_secVtx3dL	            = jet_secVtx3dL[j];	   
	bjet_secVtx3deL	            = jet_secVtx3deL[j];	   
	bjet_emfrac	   	    = jet_emfrac[j];	   
	bjet_hadfrac	   	    = jet_hadfrac[j];	   
	bjet_nNeutrals	            = jet_nNeutrals[j];	   
	bjet_nCharged	   	    = jet_nCharged[j];	   
	bjet_chadfrac	   	    = jet_chadfrac[j];	   
	bjet_nhadfrac	   	    = jet_nhadfrac[j];	   
	bjet_phofrac	   	    = jet_phofrac[j];	   
	bjet_mufrac	   	    = jet_mufrac[j];	   
	bjet_elefrac	   	    = jet_elefrac[j];	   
	bjet_JECUnc	   	    = jet_JECUnc[j];	   
	bjet_leadTrackPt	    = jet_leadTrackPt[j];	   
	bjet_softLeptPt	            = jet_softLeptPt[j];	   
	bjet_softLeptPtRel   	    = jet_softLeptPtRel[j];   
	bjet_softLeptDR	            = jet_softLeptDR[j];	  
	bjet_softLeptIdLooseMu	    = jet_softLeptIdLooseMu[j];
	bjet_softLeptIdEle95 	    = jet_softLeptIdEle95[j]; 
	bjet_genPt	   	    = jet_genPt[j];	   
	bjet_flavour         	    = jet_flavour[j];         
	outTree->Fill();
      }
    }
  }

  char thecut[100]; sprintf(thecut,"bjet_pt>20. && abs(bjet_eta)<2.5 && bjet_genPt>0.");
  printf("%i of %i jets examined are associated to a b quark.\n",outTree->GetEntries(),maxP*inTree->GetEntries());
  printf("%i of %i jets pass the quality cuts.\n",outTree->GetEntries(thecut),outTree->GetEntries());

  if(splitFile){
    outTree->Project("bPt","bjet_pt","");
    float bjet_ptWeight;
    outTree->SetBranchAddress("bjet_pt", &bjet_pt);
    TBranch *weight = outTree->Branch("bjet_ptWeight", &bjet_ptWeight, "bjet_ptWeight/F");
    nEntries = outTree->GetEntries(); 
    for(unsigned i=0; i<nEntries; ++i){
      outTree->GetEntry(i);
      bjet_ptWeight = bPt->GetBinContent(bPt->FindBin(bjet_pt))>0 ? 1.0/bPt->GetBinContent(bPt->FindBin(bjet_pt)) : 1.0;
      weight->Fill();
    }
    outFile->cd();
    bPt->Write("bPt");
  }
  
  outTree->Write();
  outFile->Close();
  inFile->Close();

  if(splitFile){
    TFile *inFile2 = new TFile(TString::Format("files/%s_regInput.root",inTreeName),"READ");
    TTree *newTree = (TTree*)inFile2->Get("EventTree");
    TFile *outFileEven = new TFile(TString::Format("files/%s_regInput_evenEvents.root",inTreeName),"RECREATE");
    TTree *newTreeEven = newTree->CopyTree("event % 2 == 0");
    newTreeEven->Write();
    outFileEven->Close();
    TFile *outFileOdd = new TFile(TString::Format("files/%s_regInput_oddEvents.root",inTreeName),"RECREATE");
    TTree *newTreeOdd = newTree->CopyTree("event % 2 != 0");
    newTreeOdd->Write();
    outFileOdd->Close();
    inFile2->Close();
  }

}
