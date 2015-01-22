#include "TROOT.h"
#include "TSystem.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "THStack.h"
#include "TTree.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TFrame.h"
#include "TStyle.h"
#include "TCut.h"
#include "TError.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TGaxis.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <utility>
#include <string>
#include <sstream>
#include "/afs/cern.ch/work/h/hebda/HggHbb/CMSSW_5_2_5/src/ggAnalysis/regression/src/tdrstyle.C"

// RooFit headers
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooVoigtian.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"

using namespace std;
using namespace RooFit;

int doAltSample=0;
bool doHistApprox=true;
void RegVal(char*, int);
void plotParameters(TCanvas *c, int Rmass, RooArgList *param, RooArgList *paramCorr);//for voit
void plotParameters(TCanvas *c, int Rmass, RooArgList *param, RooArgList *paramCorr, int fittype);//0 for gaus, 1 for voit, 2 for gausscb
void plotParameters(TCanvas *c, int Rmass, double *arrayUncorr, int nEventsUncorr, double *arrayCorr, int nEventsCorr, RooRealVar *sig1, RooRealVar *sig2 );
void sortArray(double*, const int);
double sigma_eff(double*, int, double&);
double sigma_eff(TH1*);

int main(int argc, char *argv[])
{
  if(argc==1){
    printf("RegVal [treeName] [Rmass]\n");
    return 1;
  }

  char name[80]; sprintf(name, "%s", argv[1]);

  RegVal(name, atoi(argv[2]));

  return 0;

}

void RegVal(char *treeName, int Rmass)
{
  gErrorIgnoreLevel = kError;
  gROOT->SetStyle("Plain");
  gROOT->ProcessLine(".L src/tdrstyle.C");
  setTDRStyle();
  int mass=125;

  //TCut thecut = "j1_pt>20. && j2_pt>20. && abs(j1_eta)<2.5 && abs(j2_eta)<2.5"; already applied
  TCut thecut = "1";

  TFile *file = new TFile("files/tree_Signal_regression.root");
  TTree *tree = (TTree*)file->Get(treeName);

  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  TGaxis::SetMaxDigits(3);

  // Observables
  double lowVal=40, highVal=200;
  //double lowVal=60., highVal=180.;
  int nbins=(highVal-lowVal)/3;

  //restrict the range
  double fitRangeLow=60, fitRangeHigh=180;
  TCut rangeCut = TString::Format("Mjj>%.0f && Mjj<%.0f",fitRangeLow,fitRangeHigh).Data();
  TCut rangeCutCorr = TString::Format("Mjjcorr>%.0f && Mjjcorr<%.0f",fitRangeLow,fitRangeHigh).Data();

  RooRealVar Hmass("Mjj", "M_{jj}", lowVal, highVal, "GeV");
  Hmass.setBins(nbins);
  RooRealVar HmassCorr("Mjjcorr", "M_{jj}", lowVal, highVal, "GeV");
  HmassCorr.setBins(nbins);
  RooDataSet dataset("radion", "radion", tree, RooArgList(Hmass,HmassCorr), thecut);

  double mean_jj = dataset.mean(Hmass);
  double mean_jjCorr = dataset.mean(HmassCorr);
  double rms_jj = dataset.rmsVar(Hmass)->getVal();
  double rms_jjCorr = dataset.rmsVar(HmassCorr)->getVal();
  double mean_jjErr = dataset.meanVar(Hmass)->getError();//rms_jj/sqrt(dataset.numEntries());
  double mean_jjCorrErr = dataset.meanVar(HmassCorr)->getError();//rms_jjCorr/sqrt(dataset.numEntries());
  double rms_jjErr = dataset.rmsVar(Hmass)->getError();
  double rms_jjCorrErr = dataset.rmsVar(HmassCorr)->getError();
  double oldMean, oldMeanErr, newMean, newMeanErr, resChangeV, resChangeErrV, resChangeG, resChangeErrG, resChangeS, resChangeErrS, resChangeGCB, resChangeErrGCB;

  bool GAUSS = true;
  bool GAUSSCB = true;
  bool VOIGT = true;

  TH1D* HmassUncorrHist = new TH1D("h1","h1",nbins,lowVal,highVal);
  TH1D* HmassCorrHist = new TH1D("h2","h2",nbins,lowVal,highVal);
  if(doHistApprox){
    tree->Draw("Mjj>>h1",thecut);
    tree->Draw("Mjjcorr>>h2", thecut);
  }

  tree->Draw("Mjj",thecut+rangeCut);
  double* HmassUncorrArray = tree->GetV1();
  int sizeOfHmassUncorrArray = tree->GetEntries(thecut+rangeCut);
  tree->Draw("Mjjcorr",thecut+rangeCutCorr);
  double* HmassCorrArray = tree->GetV1();
  int sizeOfHmassCorrArray = tree->GetEntries(thecut+rangeCutCorr);
  sortArray(HmassUncorrArray,sizeOfHmassUncorrArray);
  sortArray(HmassCorrArray,sizeOfHmassCorrArray);

  if(VOIGT)
    {
      RooRealVar mu_voigt("mu_voigt", "mean (uncorr.)", mean_jj, mean_jj-rms_jj, mean_jj+rms_jj, "GeV");
      RooRealVar width_voigt("width_voigt", "width (uncorr.)", rms_jj, 0.1*rms_jj, 5.*rms_jj, "GeV");
      RooRealVar sigma_voigt("sigma_voigt", "sigma (uncorr.)", 5., 0.0001, 20., "GeV");
      RooVoigtian voigt("voigt", "voigt", Hmass, mu_voigt, width_voigt, sigma_voigt);

      RooRealVar mu_voigtCorr("mu_voigtCorr", "mean (corr.)", mean_jjCorr, mean_jjCorr-rms_jjCorr, mean_jjCorr+rms_jjCorr, "GeV");
      RooRealVar width_voigtCorr("width_voigtCorr", "width (corr.)", rms_jjCorr, 0.1*rms_jjCorr, 5.*rms_jjCorr, "GeV");
      RooRealVar sigma_voigtCorr("sigma_voigtCorr", "sigma (corr.)", 5., 0.0001, 20., "GeV");
      RooVoigtian voigtCorr("voigtCorr", "voigtCorr", HmassCorr, mu_voigtCorr, width_voigtCorr, sigma_voigtCorr);

      RooPlot * jj_frame = Hmass.frame(lowVal,highVal);
      RooPlot * jjCorr_frame = HmassCorr.frame(lowVal,highVal);
      dataset.plotOn(jj_frame, LineColor(kRed+2), MarkerColor(kRed), MarkerStyle(8), MarkerSize(1.0), Binning(nbins,lowVal,highVal));
      dataset.plotOn(jjCorr_frame, LineColor(kBlue+2), MarkerColor(kBlue), MarkerStyle(8), MarkerSize(1.0), Binning(nbins,lowVal,highVal));
      RooFitResult *f = voigt.fitTo(dataset, Save(), Range(fitRangeLow,fitRangeHigh));
      voigt.plotOn(jj_frame, LineColor(kRed+2), LineWidth(3));
      RooFitResult * fCorr = voigtCorr.fitTo(dataset, Save(), Range(fitRangeLow,fitRangeHigh));
      voigtCorr.plotOn(jjCorr_frame, LineColor(kBlue+2), LineWidth(3));
      RooArgList p = f->floatParsFinal();
      RooArgList pCorr = fCorr->floatParsFinal();
      jjCorr_frame->GetYaxis()->SetTitleOffset(1.25);
      jjCorr_frame->SetTitle("");
      jjCorr_frame->Draw();
      jj_frame->Draw("same");
      double fg = 2. * sigma_voigt.getVal() * sqrt(2. * log(2));
      double fl = 2. * width_voigt.getVal();
      double fv = 0.5346 * fl + sqrt(0.2166 * pow(fl, 2.) + pow(fg, 2.));
      double fgCorr = 2. * sigma_voigtCorr.getVal() * sqrt(2. * log(2));
      double flCorr = 2. * width_voigtCorr.getVal();
      double fvCorr = 0.5346 * flCorr + sqrt(0.2166 * pow(flCorr, 2.) + pow(fgCorr, 2.));
      double res= fv / mu_voigt.getVal() * 100. / (2. * sqrt(2. * log(2.)));
      double resErr= res * sqrt(pow(mu_voigt.getError()/mu_voigt.getVal(),2) + pow(sigma_voigt.getError()/sigma_voigt.getVal(),2) + pow(width_voigt.getError()/width_voigt.getVal(),2));
      double resCorr= fvCorr / mu_voigtCorr.getVal() * 100. / (2. * sqrt(2. * log(2.)));
      double resCorrErr= res * sqrt(pow(mu_voigtCorr.getError()/mu_voigtCorr.getVal(),2) + pow(sigma_voigtCorr.getError()/sigma_voigtCorr.getVal(),2) + pow(width_voigtCorr.getError()/width_voigtCorr.getVal(),2));
      resChangeV = (res-resCorr)/res * 100.;
      resChangeErrV = fabs(resChangeV) * sqrt(pow(resCorrErr/resCorr,2) + pow(resErr/res,2));

      plotParameters( c1, Rmass, &p, &pCorr);
      if(Rmass>0) c1->Print(Form("plots/resImp_voigt_M-%i.pdf",Rmass));
      else c1->Print(Form("plots/resImp_voigt_SM.pdf"));
      c1->Clear();

    }

  if(GAUSSCB)
    {
      RooRealVar mu_gauss("mu_gauss", "mean (uncorr.)", mean_jj, mean_jj-rms_jj, mean_jj+rms_jj, "GeV");
      RooRealVar sigma_gauss("sigma_gauss", "sigma_gauss (uncorr.)", rms_jj, 0.01*rms_jj, 1000.*rms_jj, "GeV");
      RooGaussian gauss("gauss", "gauss", Hmass, mu_gauss, sigma_gauss);
      RooRealVar sigma_cb("sigma_cb", "sigma_cb (uncorr.)", rms_jj, 0.01*rms_jj, 1000.*rms_jj, "GeV");
      RooRealVar alpha("alpha","alpha (uncorr.)", 1.0, 0.5, 3);
      RooRealVar n("n", "n (uncorr.)", 4.0, 0.5, 10);
      RooCBShape cb("cb","cb",Hmass, mu_gauss,sigma_cb,alpha,n);
      RooRealVar frac("frac","frac (uncorr.)",0.2, 0, 0.35);
      RooAddPdf gausscb("gausscb","gausscb",gauss,cb,frac);

      RooRealVar mu_gaussCorr("mu_gaussCorr", "mean (corr.)", mean_jjCorr, mean_jjCorr-rms_jjCorr, mean_jjCorr+rms_jjCorr, "GeV");
      RooRealVar sigma_gaussCorr("sigma_gaussCorr", "sigma_gauss (corr.)", rms_jjCorr, 0.01*rms_jjCorr, 1000.*rms_jjCorr, "GeV");
      RooGaussian gaussCorr("gaussCorr", "gaussCorr", HmassCorr, mu_gaussCorr, sigma_gaussCorr);
      RooRealVar sigma_cbCorr("sigma_cbCorr", "sigma_cb (corr.)", rms_jjCorr, 0.01*rms_jjCorr, 1000.*rms_jjCorr, "GeV");
      RooRealVar alphaCorr("alphaCorr","alpha (corr.)", 1.0, 0.5, 3);
      RooRealVar nCorr("nCorr", "n (corr.)", 4.0, 0.5, 10);
      RooCBShape cbCorr("cbCorr","cbCorr",HmassCorr, mu_gaussCorr,sigma_cbCorr,alphaCorr,nCorr);
      RooRealVar fracCorr("fracCorr","frac (corr.)",0.2, 0, 0.35);
      RooAddPdf gausscbCorr("gausscbCorr","gausscbCorr",gaussCorr,cbCorr,fracCorr);

      RooPlot * jj_frame = Hmass.frame(lowVal,highVal);
      RooPlot * jjCorr_frame = HmassCorr.frame(lowVal,highVal);
      dataset.plotOn(jj_frame, LineColor(kRed+2), MarkerColor(kRed), MarkerStyle(8), MarkerSize(1.0), Binning(nbins,lowVal,highVal));
      dataset.plotOn(jjCorr_frame, LineColor(kBlue+2), MarkerColor(kBlue), MarkerStyle(8), MarkerSize(1.0), Binning(nbins,lowVal,highVal));
      RooFitResult *f = gausscb.fitTo(dataset, Save(), Range(fitRangeLow,fitRangeHigh));
      gausscb.plotOn(jj_frame, LineColor(kRed+2), LineWidth(3));
      RooFitResult * fCorr = gausscbCorr.fitTo(dataset, Save(), Range(fitRangeLow,fitRangeHigh));
      gausscbCorr.plotOn(jjCorr_frame, LineColor(kBlue+2), LineWidth(3));
      RooArgList p = f->floatParsFinal();
      RooArgList pCorr = fCorr->floatParsFinal();
      jjCorr_frame->GetYaxis()->SetTitleOffset(1.25);
      jjCorr_frame->SetTitle("");
      jjCorr_frame->Draw();
      jj_frame->Draw("same");
      double res = sigma_cb.getVal() / mu_gauss.getVal() * 100.;
      double resErr= res * sqrt(pow(mu_gauss.getError()/mu_gauss.getVal(),2) + pow(sigma_cb.getError()/sigma_cb.getVal(),2));
      double resCorr = sigma_cbCorr.getVal() / mu_gaussCorr.getVal() * 100.;
      double resCorrErr= res * sqrt(pow(mu_gaussCorr.getError()/mu_gaussCorr.getVal(),2) + pow(sigma_cbCorr.getError()/sigma_cbCorr.getVal(),2));
      
      oldMean = mu_gauss.getVal();
      oldMeanErr = mu_gauss.getError();
      newMean = mu_gaussCorr.getVal();
      newMeanErr = mu_gaussCorr.getError();
      resChangeGCB = (res-resCorr)/res * 100.;
      resChangeErrGCB = fabs(resChangeGCB) * sqrt(pow(resCorrErr/resCorr,2) + pow(resErr/res,2));
      plotParameters( c1, Rmass, &p, &pCorr, 2);
      if(Rmass>0) c1->Print(Form("plots/resImp_gauscb_M-%i.pdf",Rmass));
      else c1->Print(Form("plots/resImp_gauscb_SM.pdf"));
      c1->Clear();
    }

   if(GAUSS)
    {
      RooRealVar mu_gauss("mu_gauss", "mean (uncorr.)", mean_jj, mean_jj-rms_jj, mean_jj+rms_jj, "GeV");
      RooRealVar sigma_gauss("sigma_gauss", "sigma (uncorr.)", rms_jj, .01*rms_jj, 5.*rms_jj, "GeV");
      RooGaussian gauss("gauss", "gauss", Hmass, mu_gauss, sigma_gauss);

      RooRealVar mu_gaussCorr("mu_gaussCorr", "mean (corr.)", mean_jjCorr, mean_jjCorr-rms_jjCorr, mean_jjCorr+rms_jjCorr, "GeV");
      RooRealVar sigma_gaussCorr("sigma_gaussCorr", "sigma (corr.)", rms_jjCorr, .01*rms_jjCorr, 5.*rms_jjCorr, "GeV");
      RooGaussian gaussCorr("gaussCorr", "gaussCorr", HmassCorr, mu_gaussCorr, sigma_gaussCorr);

      RooPlot * jj_frame = Hmass.frame(lowVal,highVal);
      RooPlot * jjCorr_frame = HmassCorr.frame(lowVal,highVal);
      dataset.plotOn(jj_frame, LineColor(kRed+2), MarkerColor(kRed), MarkerStyle(8), MarkerSize(1.0), Binning(nbins,lowVal,highVal));
      dataset.plotOn(jjCorr_frame, LineColor(kBlue+2), MarkerColor(kBlue), MarkerStyle(8), MarkerSize(1.0), Binning(nbins,lowVal,highVal));
      RooFitResult *f = gauss.fitTo(dataset, Save(), Range(mean_jj-15, mean_jj+20));
      gauss.plotOn(jj_frame, LineColor(kRed+2), LineWidth(3));
      RooFitResult * fCorr = gaussCorr.fitTo(dataset, Save(), Range(mean_jjCorr-15, mean_jjCorr+20));
      gaussCorr.plotOn(jjCorr_frame, LineColor(kBlue+2), LineWidth(3));
      RooArgList p = f->floatParsFinal();
      RooArgList pCorr = fCorr->floatParsFinal();
      jjCorr_frame->GetYaxis()->SetTitleOffset(1.25);
      jjCorr_frame->SetTitle("");
      jjCorr_frame->Draw();
      jj_frame->Draw("same");
      double res = sigma_gauss.getVal() / mu_gauss.getVal() * 100.;
      double resErr= res * sqrt(pow(mu_gauss.getError()/mu_gauss.getVal(),2) + pow(sigma_gauss.getError()/sigma_gauss.getVal(),2));
      double resCorr = sigma_gaussCorr.getVal() / mu_gaussCorr.getVal() * 100.;
      double resCorrErr= res * sqrt(pow(mu_gaussCorr.getError()/mu_gaussCorr.getVal(),2) + pow(sigma_gaussCorr.getError()/sigma_gaussCorr.getVal(),2));

      resChangeG = (res-resCorr)/res * 100.;
      resChangeErrG = fabs(resChangeG) * sqrt(pow(resCorrErr/resCorr,2) + pow(resErr/res,2));
      plotParameters( c1, Rmass, &p, &pCorr, 0);
      if(Rmass>0) c1->Print(Form("plots/resImp_gauss_M-%i.pdf",Rmass));
      else c1->Print(Form("plots/resImp_gauss_SM.pdf"));
      c1->Clear();
    }

  double resS = doHistApprox ? sigma_eff(HmassUncorrHist)/mean_jj * 100. : sigma_eff(HmassUncorrArray,sizeOfHmassUncorrArray,mean_jj)/mean_jj * 100.;
  double resSerr = fabs(resS)*sqrt(pow(mean_jjErr/mean_jj,2)+pow(rms_jjErr/(resS*mean_jj/100.),2));
  double resSCorr = doHistApprox ? sigma_eff(HmassCorrHist)/mean_jjCorr * 100. : sigma_eff(HmassCorrArray,sizeOfHmassCorrArray,mean_jjCorr)/mean_jjCorr * 100.;
  double resSCorrerr = fabs(resSCorr)*sqrt(pow(mean_jjCorrErr/mean_jjCorr,2)+pow(rms_jjCorrErr/(resSCorr*mean_jjCorr/100.),2));
  resChangeS = (resS-resSCorr)/resS * 100.;
  resChangeErrS = fabs(resChangeS) * sqrt(pow(resSerr/resS,2) + pow(resSCorrerr/resSCorr,2));

  cout << "\n\n";
  cout << "Resolution change from SigmaEff: " << resChangeS << " +/- " << resChangeErrS << endl; 
  cout << "Resolution change from Voigtian: " << resChangeV << " +/- " << resChangeErrV << endl; 
  cout << "Resolution change from Gaussian: " << resChangeG << " +/- " << resChangeErrG << endl; 
  cout << "Resolution change from Gauss+CB: " << resChangeGCB << " +/- " << resChangeErrGCB << endl; 
  cout << "Distro mean shifted from " << mean_jj << " +/- " << mean_jjErr << " to " << mean_jjCorr << " +/- " << mean_jjCorrErr << endl;
  if(GAUSSCB) cout << "Fitted mean shifted from " << oldMean << " +/- " << oldMeanErr << " to " << newMean << " +/- " << newMeanErr << endl;

  FILE* output = fopen("results.txt","a");
  if(GAUSSCB){
    if(Rmass>0) fprintf(output,"%s %4i %5.2f +/- %4.2f :: %5.1f +/- %3.1f to %5.1f +/- %3.1f \n",treeName,Rmass,resChangeGCB,resChangeErrGCB,oldMean,oldMeanErr,newMean,newMeanErr);
    else fprintf(output,"       %s   SM %5.2f +/- %4.2f :: %5.1f +/- %3.1f to %5.1f +/- %3.1f \n",treeName,resChangeGCB,resChangeErrGCB,oldMean,oldMeanErr,newMean,newMeanErr);
  }
  else{
    if(Rmass>0) fprintf(output,"%s %4i %5.2f +/- %4.2f :: %5.1f +/- %3.1f to %5.1f +/- %3.1f \n",treeName,Rmass,resChangeGCB,resChangeErrGCB,mean_jj,mean_jjErr,mean_jjCorr,mean_jjCorrErr);
    else fprintf(output,"       %s   SM %5.2f +/- %4.2f :: %5.1f +/- %3.1f to %5.1f +/- %3.1f \n",treeName,resChangeGCB,resChangeErrGCB,mean_jj,mean_jjErr,mean_jjCorr,mean_jjCorrErr);
  }
  fclose(output);

  file->Close();

  return;
}

void plotParameters(TCanvas *c, int Rmass, RooArgList *param, RooArgList *paramCorr, int fittype )
{
  c->cd(1);

  double mean1, mean2, sigma1, sigma2, width1, width2, res1, res2;
  double mean1err, mean2err, sigma1err, sigma2err, width1err, width2err, res1err, res2err;

  RooRealVar* obj = new RooRealVar();
  TIterator *it = (TIterator*) param->createIterator();
  if(fittype==2) {it->Next();it->Next();} 
  obj = (RooRealVar*)it->Next();
  mean1 = obj->getVal();
  mean1err = obj->getError();
  if(fittype==2) {it->Next();} 
  obj = (RooRealVar*)it->Next();
  sigma1 = obj->getVal();
  sigma1err = obj->getError();
  if(fittype==1){
    obj = (RooRealVar*)it->Next();
    width1 = obj->getVal();
    width1err = obj->getError();
  }
  it = (TIterator*) paramCorr->createIterator();
  if(fittype==2) {it->Next();it->Next();} 
  obj = (RooRealVar*)it->Next();
  mean2 = obj->getVal();
  mean2err = obj->getError();
  if(fittype==2) {it->Next();} 
  obj = (RooRealVar*)it->Next();
  sigma2 = obj->getVal();
  sigma2err = obj->getError();
  if(fittype==1){
    obj = (RooRealVar*)it->Next();
    width2 = obj->getVal();
    width2err = obj->getError();
  }

  if(fittype==1){
    double fg = 2. * sigma1 * sqrt(2. * log(2));
    double fl = 2. * width1;
    double fv = 0.5346 * fl + sqrt(0.2166 * pow(fl, 2.) + pow(fg, 2.));
    res1 = fv / mean1 * 100. / (2. * sqrt(2. * log(2.)));
    res1err = res1 * sqrt(pow(mean1err/mean1,2) + pow(sigma1err/sigma1,2) + pow(width1err/width1,2));

    fg = 2. * sigma2 * sqrt(2. * log(2));
    fl = 2. * width2;
    fv = 0.5346 * fl + sqrt(0.2166 * pow(fl, 2.) + pow(fg, 2.));
    res2 = fv / mean2 * 100. / (2. * sqrt(2. * log(2.)));
    res2err = res2 * sqrt(pow(mean2err/mean2,2) + pow(sigma2err/sigma2,2) + pow(width2err/width2,2));
  }
  else{
    res1 = sigma1/mean1 * 100.;
    res1err = res1 * sqrt(pow(mean1err/mean1,2) + pow(sigma1err/sigma1,2));
    res2 = sigma2/mean2 * 100.;
    res2err = res2 * sqrt(pow(mean2err/mean2,2) + pow(sigma2err/sigma2,2));
  }

  double resChange = fabs(res1-res2)/res1 * 100.;
  double resChangeErr = resChange * sqrt(pow(res1err/res1,2) + pow(res2err/res2,2));

  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextSize(0.045);
  latexLabel.SetTextAlign(12);
  latexLabel.SetTextColor(kBlack);
  latexLabel.DrawLatex(0.17,0.89,"CMS Simulation");
  latexLabel.SetTextSize(0.035);
  latexLabel.DrawLatex(0.17,0.84,"#sqrt{s} = 8 TeV");
  switch(doAltSample){
  case 0:
    if(Rmass>0){
      latexLabel.DrawLatex(0.17,0.80,"R#rightarrow H(#gamma#gamma)H(b#bar{b})");
      latexLabel.DrawLatex(0.17,0.75,TString::Format("M_{R} = %i GeV",Rmass));
    }
    else{
      latexLabel.DrawLatex(0.17,0.80,"gg#rightarrow H(#gamma#gamma)H(b#bar{b})");
    }
    break;
  case 1: latexLabel.DrawLatex(0.17,0.80,"ZH#rightarrow l#bar{l} b#bar{b}"); break;
  case 2: latexLabel.DrawLatex(0.17,0.80,"HZ#rightarrow #gamma#gamma b#bar{b}"); break;
  case 3: latexLabel.DrawLatex(0.17,0.80,"HZ#rightarrow #nu#bar{#nu} b#bar{b}"); break;
  case 4: latexLabel.DrawLatex(0.17,0.80,"ZZ#rightarrow l#bar{l} b#bar{b}"); break;
  }

  int pos=0;
  float topPos=0.89, spacing=0.03;
  latexLabel.SetTextSize(0.025);
  latexLabel.SetTextColor(kRed+2);
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"Before correction:");
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Mean = (%.1f #pm %.1f) GeV",mean1,mean1err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Sigma = (%.1f #pm %.1f) GeV",sigma1,sigma1err));
  if(fittype==1) latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Width = (%.1f #pm %.1f) GeV",width1,width1err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Res. = %.1f%% #pm %.1f%%",res1,res1err));
  
  pos++;
  latexLabel.SetTextSize(0.025);
  latexLabel.SetTextColor(kBlue+2);
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"After correction:");
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Mean = (%.1f #pm %.1f) GeV",mean2,mean2err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Sigma = (%.1f #pm %.1f) GeV",sigma2,sigma2err));
  if(fittype==1) latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Width = (%.1f #pm %.1f) GeV",width2,width2err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Res. = %.1f%% #pm %.1f%%",res2,res2err));
  
  pos++;
  latexLabel.SetTextSize(0.025);
  latexLabel.SetTextColor(kBlack);
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"Resolution Change:");
  if(res1>res2)
    latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  + %.1f%% #pm %.1f%%",resChange,resChangeErr));
  else
    latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  - %.1f%% #pm %.1f%%",resChange,resChangeErr));

  latexLabel.Draw();

  return;
}

void plotParameters(TCanvas *c, int Rmass, RooArgList *param, RooArgList *paramCorr)
{
  c->cd(1);

  double mean1, mean2, sigma1, sigma2, width1, width2, res1, res2;
  double mean1err, mean2err, sigma1err, sigma2err, width1err, width2err, res1err, res2err;

  RooRealVar* obj = new RooRealVar();
  TIterator *it = (TIterator*) param->createIterator();
  obj = (RooRealVar*)it->Next();
  mean1 = obj->getVal();
  mean1err = obj->getError();
  obj = (RooRealVar*)it->Next();
  sigma1 = obj->getVal();
  sigma1err = obj->getError();
  obj = (RooRealVar*)it->Next();
  width1 = obj->getVal();
  width1err = obj->getError();
  
  it = (TIterator*) paramCorr->createIterator();
  obj = (RooRealVar*)it->Next();
  mean2 = obj->getVal();
  mean2err = obj->getError();
  obj = (RooRealVar*)it->Next();
  sigma2 = obj->getVal();
  sigma2err = obj->getError();
  obj = (RooRealVar*)it->Next();
  width2 = obj->getVal();
  width2err = obj->getError();

  double fg = 2. * sigma1 * sqrt(2. * log(2));
  double fl = 2. * width1;
  double fv = 0.5346 * fl + sqrt(0.2166 * pow(fl, 2.) + pow(fg, 2.));
  res1 = fv / mean1 * 100. / (2. * sqrt(2. * log(2.)));
  res1err = res1 * sqrt(pow(mean1err/mean1,2) + pow(sigma1err/sigma1,2) + pow(width1err/width1,2));

  fg = 2. * sigma2 * sqrt(2. * log(2));
  fl = 2. * width2;
  fv = 0.5346 * fl + sqrt(0.2166 * pow(fl, 2.) + pow(fg, 2.));
  res2 = fv / mean2 * 100. / (2. * sqrt(2. * log(2.)));
  res2err = res2 * sqrt(pow(mean2err/mean2,2) + pow(sigma2err/sigma2,2) + pow(width2err/width2,2));

  double resChange = fabs(res1-res2)/res1 * 100.;
  double resChangeErr = resChange * sqrt(pow(res1err/res1,2) + pow(res2err/res2,2));

  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextSize(0.045);
  latexLabel.SetTextAlign(12);
  latexLabel.SetTextColor(kBlack);
  latexLabel.DrawLatex(0.17,0.89,"CMS Simulation");
  latexLabel.SetTextSize(0.035);
  latexLabel.DrawLatex(0.17,0.84,"#sqrt{s} = 8 TeV");
  switch(doAltSample){
  case 0:
    if(Rmass>0){
      latexLabel.DrawLatex(0.17,0.80,"R#rightarrow H(#gamma#gamma)H(b#bar{b})");
      latexLabel.DrawLatex(0.17,0.75,TString::Format("M_{R} = %i GeV",Rmass));
    }
    else{
      latexLabel.DrawLatex(0.17,0.80,"gg#rightarrow H(#gamma#gamma)H(b#bar{b})");
    }
    break;
  case 1: latexLabel.DrawLatex(0.17,0.80,"ZH#rightarrow l#bar{l} b#bar{b}"); break;
  case 2: latexLabel.DrawLatex(0.17,0.80,"HZ#rightarrow #gamma#gamma b#bar{b}"); break;
  case 3: latexLabel.DrawLatex(0.17,0.80,"HZ#rightarrow #nu#bar{#nu} b#bar{b}"); break;
  case 4: latexLabel.DrawLatex(0.17,0.80,"ZZ#rightarrow l#bar{l} b#bar{b}"); break;
  }

  int pos=0;
  float topPos=0.89, spacing=0.03;
  latexLabel.SetTextSize(0.025);
  latexLabel.SetTextColor(kRed+2);
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"Before correction:");
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Mean = (%.1f #pm %.1f) GeV",mean1,mean1err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Sigma = (%.1f #pm %.1f) GeV",sigma1,sigma1err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Width = (%.1f #pm %.1f) GeV",width1,width1err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Res. = %.1f%% #pm %.1f%%",res1,res1err));
  
  pos++;
  latexLabel.SetTextSize(0.025);
  latexLabel.SetTextColor(kBlue+2);
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"After correction:");
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Mean = (%.1f #pm %.1f) GeV",mean2,mean2err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Sigma = (%.1f #pm %.1f) GeV",sigma2,sigma2err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Width = (%.1f #pm %.1f) GeV",width2,width2err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Res. = %.1f%% #pm %.1f%%",res2,res2err));
  
  pos++;
  latexLabel.SetTextSize(0.025);
  latexLabel.SetTextColor(kBlack);
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"Resolution Change:");
  if(res1>res2)
    latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  + %.1f%% #pm %.1f%%",resChange,resChangeErr));
  else
    latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  - %.1f%% #pm %.1f%%",resChange,resChangeErr));

  latexLabel.Draw();

  return;
}

void plotParameters(TCanvas *c, int Rmass, double *arrayUncorr, int nEventsUncorr, double *arrayCorr, int nEventsCorr, RooRealVar *sig1, RooRealVar *sig2 )
{
  c->cd(1);

  double mean1=-1, mean2=-1, sigma1, sigma2, res1, res2;
  double mean1err, mean2err, sigma1err, sigma2err, res1err, res2err;

  if(doHistApprox){
    sigma1 = sig1->getVal();
    sigma2 = sig2->getVal();
    for(int i=0; i<nEventsUncorr; ++i)
      mean1+=arrayUncorr[i];
    for(int i=0; i<nEventsCorr; ++i)
      mean2+=arrayCorr[i];
    mean1/=nEventsUncorr;
    mean2/=nEventsCorr;
  }
  else{
    sigma1 = sigma_eff(arrayUncorr, nEventsUncorr, mean1);
    sigma2 = sigma_eff(arrayCorr, nEventsCorr, mean2);
  }

  res1 = sigma1/mean1 * 100.;
  res2 = sigma2/mean2 * 100.;

  double tmp1=0.0, tmp2=0.0;
  for(int i=0; i<nEventsUncorr; ++i)
    tmp1+=pow(arrayUncorr[i]-mean1,2);
  for(int i=0; i<nEventsCorr; ++i)
    tmp2+=pow(arrayCorr[i]-mean2,2);
  mean1err = sqrt(tmp1)/(nEventsUncorr-1);
  mean2err = sqrt(tmp2)/(nEventsCorr-1);
  sigma1err = sig1->getError();
  sigma2err = sig2->getError();

  res1err = fabs(res1) * sqrt(pow(sigma1err/sigma1,2) + pow(mean1err/mean1,2));
  res2err = fabs(res2) * sqrt(pow(sigma2err/sigma2,2) + pow(mean2err/mean2,2));

  double resChange = fabs(res1-res2)/res1 * 100.;
  double resChangeErr = resChange * sqrt(pow(res1err/res1,2) + pow(res2err/res2,2));

  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextSize(0.045);
  latexLabel.SetTextAlign(12);
  latexLabel.SetTextColor(kBlack);
  latexLabel.DrawLatex(0.17,0.89,"CMS Simulation");
  latexLabel.SetTextSize(0.035);
  latexLabel.DrawLatex(0.17,0.84,"#sqrt{s} = 8 TeV");
  switch(doAltSample){
  case 0:
    if(Rmass>0){
      latexLabel.DrawLatex(0.17,0.80,"R#rightarrow H(#gamma#gamma)H(b#bar{b})");
      latexLabel.DrawLatex(0.17,0.75,TString::Format("M_{R} = %i GeV",Rmass));
    }
    else{
      latexLabel.DrawLatex(0.17,0.80,"gg#rightarrow H(#gamma#gamma)H(b#bar{b})");
    }
    break;
  case 1: latexLabel.DrawLatex(0.17,0.80,"ZH#rightarrow l#bar{l} b#bar{b}"); break;
  case 2: latexLabel.DrawLatex(0.17,0.80,"HZ#rightarrow #gamma#gamma b#bar{b}"); break;
  case 3: latexLabel.DrawLatex(0.17,0.80,"HZ#rightarrow #nu#bar{#nu} b#bar{b}"); break;
  case 4: latexLabel.DrawLatex(0.17,0.80,"ZZ#rightarrow l#bar{l} b#bar{b}"); break;
  }

  int pos=0;
  float topPos=0.89, spacing=0.03;
  latexLabel.SetTextSize(0.025);
  latexLabel.SetTextColor(kRed+2);
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"Before correction:");
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Mean = (%.1f #pm %.1f) GeV",mean1,mean1err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Sig. Eff. = (%.1f #pm %.1f) GeV",sigma1,sigma1err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Res. = %.1f%% #pm %.1f%%",res1,res1err));
  
  pos++;
  latexLabel.SetTextSize(0.025);
  latexLabel.SetTextColor(kBlue+2);
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"After correction:");
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Mean = (%.1f #pm %.1f) GeV",mean2,mean2err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Sig. Eff. = (%.1f #pm %.1f) GeV",sigma2,sigma2err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Res. = %.1f%% #pm %.1f%%",res2,res2err));
  
  pos++;
  latexLabel.SetTextSize(0.025);
  latexLabel.SetTextColor(kBlack);
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"Resolution Change:");
  if(res1>res2)
    latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  + %.1f%% #pm %.1f%%",resChange,resChangeErr));
  else
    latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  - %.1f%% #pm %.1f%%",resChange,resChangeErr));

  latexLabel.Draw();

  return;
}

void sortArray(double* a, const int len){
  double tmpVal;
  int tmpIndex;

  for(int i=0; i<len; ++i){
    tmpVal=a[i];
    tmpIndex=i;
    for(int j=i+1; j<len; ++j){
      if(a[j]<tmpVal){
	tmpVal=a[j];
	tmpIndex=j;
      }
    }
    if(tmpIndex!=i){
      a[tmpIndex]=a[i];
      a[i]=tmpVal;
    }
  }

  cout<<"!!! "<<len<<' '<<a[0]<<' '<<a[len-1]<<' '<<endl;
}

double sigma_eff(double *data, int nEvents, double &mean){
  cout<<"@@@ "<<nEvents<<' '<<data[0]<<' '<<data[nEvents-1]<<' '<<endl;
  sortArray(data,nEvents);
  cout<<"@@@ "<<nEvents<<' '<<data[0]<<' '<<data[nEvents-1]<<' '<<endl;

  if(mean<0){
    mean=0.0;
    for(int i=0; i<nEvents; ++i) mean+=data[i];
    mean/=nEvents;
  }

  double width=1.0, step=0.01, frac=0.0;
  do{
    width+=step;
    int nEventsPass=0;
    for(int i=0; i<nEvents && data[i]<mean+width/2.0; ++i){
      if(data[i]>mean-width/2.0 && data[i]<mean+width/2.0) ++nEventsPass;
    }
    frac=(double)nEventsPass/nEvents;
  }
  while(frac<0.6827);

  return width/2.0;
}

double sigma_eff(TH1 *hist) {

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    std::cout << "effsigma: Not a valid histo. nbins = " << nb << std::endl;
    return 0.;
  }

  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    std::cout << "effsigma: Not a valid histo. bwid = " << bwid << std::endl;
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
  if(total < 100.) {
    std::cout << "effsigma: Too few entries " << total << std::endl;
    return 0.;
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
  if(ierr != 0) std::cout << "effsigma: Error of type " << ierr << std::endl;

  return widmin;
}

