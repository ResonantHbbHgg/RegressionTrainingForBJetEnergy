// @(#)root/tmva $Id: TMVARegression.C 38475 2011-03-17 10:46:00Z evt $
/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVARegression                                                     *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l TMVARegression.C\(\"LD,MLP\"\)                                      *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set is used.                                     *
 *                                                                                *
 * The output file "TMVAReg.root" can be analysed with the use of dedicated       *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 **********************************************************************************/

#include <cstdlib>
#include <iostream> 
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVARegGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#endif

using namespace TMVA;
   
void TMVARegression( TString Cat="resonant" , TString myMethodList = "" ) 
{

   // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
   // if you use your private .rootrc, or run from a different directory, please copy the 
   // corresponding lines from .rootrc

   // methods to be processed can be given as an argument; use format:
   //
   // mylinux~> root -l TMVARegression.C\(\"myMethod1,myMethod2,myMethod3\"\)
   //

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDEFoam"]         = 0; 
   Use["KNN"]             = 0;
   // 
   // --- Linear Discriminant Analysis
   Use["LD"]		        = 0;
   // 
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   // 
   // --- Neural Network
   Use["MLP"]             = 0; 
   // 
   // --- Support Vector Machine 
   Use["SVM"]             = 0;
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 0;
   Use["BDTG"]            = 1;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVARegression" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Here the preparation phase begins

   // Create a new root output file
   TString outfileName( Form("TMVAReg.root") );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory will
   // then run the performance analysis for you.
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/ 
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in 
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory = new TMVA::Factory( Form("TMVARegression_%s",Cat.Data()), outputFile, 
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Regression" );

   // If you wish to modify default settings 
   // (please check "src/Config.h" to see all available global options)
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
   // You can add so-called "Spectator variables", which are not used in the MVA training, 
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the 
   // input variables, the response values of all trained MVAs, and the spectator variables

   factory->AddVariable( "bjet_pt", "bjet_pt", "units", 'F' );
   factory->AddVariable( "bjet_eta", "bjet_eta", "units", 'F' );
   factory->AddVariable( "bjet_mt", "bjet_mt", "units", 'F' );
   factory->AddVariable( "bjet_phofrac", "bjet_phofrac", "units", 'F' );
   factory->AddVariable( "bjet_nhadfrac", "bjet_nhadfrac", "units", 'F' );
   factory->AddVariable( "bjet_softLeptPt:=(bjet_softLeptIdLooseMu==1 || bjet_softLeptIdEle95==1) ? (bjet_softLeptPt) : (-99)", "bjet_softLeptPt", "units", 'F' );
   factory->AddVariable( "bjet_softLeptPtRel:=(bjet_softLeptIdLooseMu==1 || bjet_softLeptIdEle95==1) ? (bjet_softLeptPtRel) : (-99)", "bjet_softLeptPtRel", "units", 'F' );
   factory->AddVariable( "bjet_secVtxM", "bjet_secVtxM", "units", 'F' );
   factory->AddVariable( "bjet_secVtx3deL", "bjet_secVtx3deL", "units", 'F' );
   factory->AddVariable( "MET", "MET", "units", 'F' );
   factory->AddVariable( "dPhiMETJet:=(abs(bjet_phi-METphi)>3.14159265 ) ? (2*3.14159265-abs(bjet_phi-METphi)) : (abs(bjet_phi-METphi))", "dPhiMETJet", "units", 'F' );
   factory->AddVariable( "bjet_leadTrackPt", "bjet_leadTrackPt", "units", 'F' );
   factory->AddVariable( "bjet_ntot:=bjet_nCharged+bjet_nNeutrals", "bjet_ntot", "units", 'F' );
   factory->AddVariable( "rho", "rho", "units", 'F' );

   factory->AddSpectator( "bjet_mufrac", "bjet_mufrac", "units", 'F' );
   factory->AddSpectator( "bjet_elefrac", "bjet_elefrac", "units", 'F' );
   factory->AddSpectator( "bjet_chadfrac", "bjet_chadfrac", "units", 'F' );
   factory->AddSpectator( "bjet_softLeptDR:=(bjet_softLeptIdLooseMu==1 || bjet_softLeptIdEle95==1) ? (bjet_softLeptDR) : (-99)", "bjet_softLeptDR", "units", 'F' );
   factory->AddSpectator( "bjet_secVtxPt", "bjet_secVtxPt", "units", 'F' );
   factory->AddSpectator( "bjet_secVtx3dL", "bjet_secVtx3dL", "units", 'F' );
   factory->AddSpectator( "bjet_JECUnc", "bjet_JECUnc", "units", 'F' );
   
   // Add the variable carrying the regression target
   factory->AddTarget( "bjet_genPt" ); 

   TChain chainTraining("EventTree");
   TChain chainTesting("EventTree");
   if(Cat.Contains("resonant")){
     chainTraining.Add("files/MSSM_H_m260_8TeV_regInput_evenEvents.root");//128546 events
     chainTraining.Add("files/MSSM_H_m300_8TeV_regInput_evenEvents.root");//134924 events
     chainTraining.Add("files/MSSM_H_m350_8TeV_regInput_evenEvents.root");//141022 events
     chainTesting.Add("files/MSSM_H_m260_8TeV_regInput_oddEvents.root");
     chainTesting.Add("files/MSSM_H_m300_8TeV_regInput_oddEvents.root");
     chainTesting.Add("files/MSSM_H_m350_8TeV_regInput_oddEvents.root");
   }
   else{
     chainTraining.Add("files/ggHH_8TeV_regInput_evenEvents.root");//47405 events, almost an order of magnitude fewer events to the resonant training.
     chainTesting.Add("files/ggHH_8TeV_regInput_oddEvents.root");
   }
   TTree *regTreeTraining = (TTree*) chainTraining;
   TTree *regTreeTesting = (TTree*) chainTesting;

   // global event weights per tree (see below for setting event-wise weights)
   Double_t regWeight  = 1.0;   

   // You can add an arbitrary number of regression trees
   factory->AddTree( regTreeTraining, "Regression", regWeight, "", "train" );
   factory->AddTree( regTreeTesting, "Regression", regWeight, "", "test" );
   if(!Cat.Contains("resonant")) factory->SetWeightExpression( "bjet_ptWeight", "Regression" );

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycut = "bjet_pt>20. && abs(bjet_eta)<2.5 && bjet_genPt>0.";

   // tell the factory to use all remaining events in the trees after training for testing:
   // factory->PrepareTrainingAndTestTree(mycut, "nTrain_Regression=600000:nTest_Regression=600000:SplitMode=Random:NormMode=NumEvents:!V");
   factory->PrepareTrainingAndTestTree(mycut, "V");

   // If no numbers of events are given, half of the events in the tree are used 
   // for training, and the other half for testing:
   //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );  

   // ---- Book MVA methods
   //
   // please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // PDE - RS method
   if (Use["PDERS"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERS", 
                           "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=40:NEventsMax=60:VarTransform=None" );
   // And the options strings for the MinMax and RMS methods, respectively:
   //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );   
   //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );   

   if (Use["PDEFoam"])
       factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoam", 
			    "!H:!V:MultiTargetRegression=F:TargetSelection=Mpv:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Compress=T:Kernel=None:Nmin=10:VarTransform=None" );

   // K-Nearest Neighbour classifier (KNN)
   if (Use["KNN"])
      factory->BookMethod( TMVA::Types::kKNN, "KNN", 
                           "nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

   // Linear discriminant
   if (Use["LD"])
      factory->BookMethod( TMVA::Types::kLD, "LD", 
                           "!H:!V:VarTransform=None" );

	// Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
   if (Use["FDA_MC"]) 
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MC",
                          "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100):FitMethod=MC:SampleSize=100000:Sigma=0.1:VarTransform=D" );
   
   if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options) .. the formula of this example is good for parabolas
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GA",
                           "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100):FitMethod=GA:PopSize=100:Cycles=3:Steps=30:Trim=True:SaveBestGen=1:VarTransform=Norm" );

   if (Use["FDA_MT"]) 
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MT",
                           "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

   if (Use["FDA_GAMT"]) 
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GAMT",
                           "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

   // Neural network (MLP)
   if (Use["MLP"])
      factory->BookMethod( TMVA::Types::kMLP, "MLP", "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator" );

   // Support Vector Machine
   if (Use["SVM"])
      factory->BookMethod( TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );
//100, 5
   // Boosted Decision Tree   //100,5 nCuts=-1
   if (Use["BDT"])
     factory->BookMethod( TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=100:nEventsMin=4:BoostType=AdaBoostR2:SeparationType=RegressionVariance:nCuts=-1:PruneMethod=CostComplexity:PruneStrength=30" );

   int ntrees = 1200;//800 is default
   float shrinkage = 0.1;//between 0.1 and 0.3
   float gradbaggingfrac = 0.7;//between 0.5 and 0.8
   int maxdepth = 3;//between 2 and 4
   if(!Cat.Contains("resonant")){
     ntrees=600; shrinkage=0.1; gradbaggingfrac=0.7; maxdepth=3;
   }

   if (Use["BDTG"]){
     factory->BookMethod( TMVA::Types::kBDT, "BDTG", 
			    Form("!H:!V:NTrees=%i::BoostType=Grad:Shrinkage=%.1f:UseBaggedGrad:GradBaggingFraction=%.1f:nCuts=200:MaxDepth=%i:NegWeightTreatment=IgnoreNegWeights",
				 ntrees,shrinkage,gradbaggingfrac,maxdepth) );
   }

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();    

   // --------------------------------------------------------------
   
   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVARegression is done!" << std::endl;      

   delete factory;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVARegGui( outfileName );
}
