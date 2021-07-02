// (c) MPI 2020
#include "AnalysisEventShapes.h"
 
// // C++ includes
// #include <map>
// #include <set>
#include <vector>

// Root includes
// #include <TH1.h>
// #include <TH2.h>
 
// // H1 includes
#include "H1PhysUtils/H1BoostedJets.h"
#include "H1HadronicCalibration/H1HadronicCalibration.h"

#include "H1Calculator/H1CalcGenericInterface.h"
using namespace H1CalcGenericInterface;

#include "H1Calculator/H1Calculator.h"
#include "H1Calculator/H1CalcTrig.h"
#include "H1Calculator/H1CalcWeight.h"
#include "H1Calculator/H1CalcVertex.h"
#include "H1Calculator/H1CalcEvent.h"
#include "H1Calculator/H1CalcKine.h"
#include "H1Calculator/H1CalcElec.h"
#include "H1Calculator/H1CalcFs.h"
#include "H1Calculator/H1CalcHad.h"
#include "H1Calculator/H1CalcTrack.h"
#include "H1Calculator/H1CalcSystematic.h"
#include "H1Mods/H1PartCandArrayPtr.h"

// #include "H1PhysUtils/H1MakeKine.h"
#include "EventshapeTools.h"
#include "JetTools.h"         // JetsAtHighQ2
#include "H1Mods/H1PartSelTrack.h"
#include "H1JetFinder/H1EventShape.h"
//#include "UsefulTools.h"


// fjcontrib
#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/Centauro.hh"



using namespace std;
using namespace fastjet;

// _______________________________________________________ //
//! Constructor
AnalysisEventShapes::AnalysisEventShapes(TString chain) : AnalysisBase(chain) {
}

// _______________________________________________________ //
AnalysisEventShapes::~AnalysisEventShapes() {
}



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoInitialSettings()
//!
//!  This function is called by the main program
//!  for the very first event once
//!
void AnalysisEventShapes::DoInitialSettings() {
   // nothing todo
   // ... initialize reweightings etc...
}



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoReset()
//!
//!  This function is called by the main program
//!  at the beginning of the event loop
//!
//!  Reset all members.
//!  Note: Also AnalysisBase::DoReset() is called
//!
void AnalysisEventShapes::DoReset() {
   // -- reset event quantites
   fGen     = CrossSectionQuantities();
   fRec     = CrossSectionQuantities();
   fTreeVar = TreeVariables();

}



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoAnalysisCutsGen()
//!
//!  This function is called by the main program
//!  at the beginning of the event loop
//!
//!  Define analysis specific generator level cuts
//!
bool AnalysisEventShapes::DoAnalysisCutsGen() {
   // -- reset event quantites
   fAnalysisCutsGen = true;

   //Exclude MC Events to avoid double counting
   // Q2<4 is covered by the Pythia photo production
   // Q2>60 is covered by Django and Rapgap


   if ( IsBkgMC && fChainName=="DjBkg" ){
      if( gH1Calc->Kine()->GetQ2Gen() > 60 || gH1Calc->Kine()->GetQ2Gen() < 4) {
         fAnalysisCutsGen = false;
      }
   }

   // set fGen.IsGood
   fGen.IsGood  =  fAnalysisCutsGen && fBasicCutsGen;
   return fAnalysisCutsGen;
}

// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoAnalysisCutsRec()
//!
//!  This function is called by the main program
//!  at the beginning of the event loop
//!
//!  Define analysis specific detector level cuts
//!
bool AnalysisEventShapes::DoAnalysisCutsRec() {
   // -- reset event quantites
   fAnalysisCutsRec = true;
   
   // set fRec.IsGood
   fRec.IsGood  =  fAnalysisCutsRec && fBasicCutsRec;

   return fAnalysisCutsRec;
}



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoCrossSectionObservablesGen()
//!
//!  This function is called by the main program
//!  during the event loop
//!
//!  Set observables needed to make the cross section
//!  histograms. Store them in (GenLevQuantities) fGen
//!
void AnalysisEventShapes::DoCrossSectionObservablesGen() {


   static EventshapeTools ESTools;
   // event weight
   // fGen.wgt      = gH1Calc->Weight()->GetWeightGen(); // gH1Calc->GetFloatVariable(Weight_Weight);
   //fGen.BoostToBreit = ESTools.BoostToBreitFrame(fGen.Q2, fGen.Y, fGen.X, ScatElecGen.Phi());
   //fGen.BoostToLab = ESTools.BoostToLabFrame(fGen.BoostToBreit);
   //mini-tree                                                                                                                                                                                            
  

   //  TLorentzVector virtualphoton = ebeam - ScatElec;

   // --------------- general kinematics -------------------                                                                                                                                                
   TLorentzVector ebeam, pbeam;
   H1BoostedJets::Instance()->BeamVectors(&pbeam, &ebeam);

   fTreeVar.event_weight = gH1Calc->Weight()->GetWeightGen();
   TLorentzVector ScatElecGen = gH1Calc->Elec()->GetFirstElectronGen();
   TLorentzVector virtual_photon = ebeam - ScatElecGen;

   fTreeVar.gen_event_Q2 =  gH1Calc->Kine()->GetQ2Gen();
   fTreeVar.gen_event_y  =  gH1Calc->Kine()->GetYGen();
   fTreeVar.gene_px = ScatElecGen.Px();
   fTreeVar.gene_py = ScatElecGen.Py();
   fTreeVar.gene_pz = ScatElecGen.Pz();
   fGen.Q2 = gH1Calc->Kine()->GetQ2Gen();

}



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoCrossSectionObservablesRec()
//!
//!  This function is called by the main program
//!  during the event loop
//!
//!  Set observables needed to make the cross section
//!  histograms. Store them in (RecLevQuantities) fRec
//!
void AnalysisEventShapes::DoCrossSectionObservablesRec() {

   static EventshapeTools ESTools;
   // event weight
   fRec.wgt      = gH1Calc->Weight()->GetWeight();
   // initial and final electron and quark
   TLorentzVector ebeam, pbeam;
   H1BoostedJets::Instance()->BeamVectors(&pbeam, &ebeam);
   TLorentzVector ScatElec = gH1Calc->Elec()->GetFirstElectron();
   TLorentzVector virtual_photon = ebeam - ScatElec;



   vector<H1PartCand*> particlearray = to_vector<H1PartCand*>(H1BoostedJets::Instance()->GetHFSArray());
    
   //minitree

   Float_t hadptda = gH1Calc->Fs()->GetHadPtDa();
   TLorentzVector HFS = gH1Calc->Fs()->GetAllFsLessElectron();


   fTreeVar.event_Q2_e = gH1Calc->Kine()->GetQ2e();
   fTreeVar.event_y_e  = gH1Calc->Kine()->GetYe();
   fTreeVar.event_Q2_sigma = gH1Calc->Kine()->GetQ2s();
   fTreeVar.event_y_sigma  = gH1Calc->Kine()->GetYs();
   fTreeVar.event_Q2_da = gH1Calc->Kine()->GetQ2da();
   fTreeVar.event_y_da  = gH1Calc->Kine()->GetYda();
   fTreeVar.event_Q2_esigma = gH1Calc->Kine()->GetQ2es();
   fTreeVar.event_y_esigma  = gH1Calc->Kine()->GetYes();

   fTreeVar.Empz     = gH1Calc->Fs()->GetEmpz();
   fTreeVar.e_px = ScatElec.Px();
   fTreeVar.e_py = ScatElec.Py();
   fTreeVar.e_pz = ScatElec.Pz();
   fTreeVar.HFS_px = HFS.Px();
   fTreeVar.HFS_py = HFS.Py();
   fTreeVar.HFS_pz = HFS.Pz();

   fTreeVar.ptmiss =  gH1Calc->Fs()->GetPtMiss();
   fTreeVar.pth    =  gH1Calc->Fs()->GetPtCalo();
   fTreeVar.vertex_z = gH1Calc->Vertex()->GetZ();
   fTreeVar.ptratio_da = HFS.Pt()/hadptda;
   fTreeVar.ptratio_ele = HFS.Pt()/ScatElec.Pt();

 
      
}//end main function



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoControlPlotsGen()
//!
//!  This function is called by the main program
//!  during the event loop
//!
void AnalysisEventShapes::DoControlPlotsGen() {
   
  return;
}



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoControlPlotsRec()
//!
//!  This function is called by the main program
//!  during the event loop
//!
void AnalysisEventShapes::DoControlPlotsRec() {
   
   static EventshapeTools ESTools;

   // --- event weight
   double wgt = fRec.wgt;

   return;

}



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoControlPlotsGenRec()
//!
//!  This function is called by the main program
//!  during the event loop
//!
void AnalysisEventShapes::DoControlPlotsGenRec() {
  return;
}



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoCrossSectionsGenRec()
//!
//!  This function is called by the main program
//!  during the event loop
//!
//!  Note: Fill histograms, but make USE ONLY of
//!  observables stored previously in CrossSectionQuantities
//! 
void AnalysisEventShapes::DoCrossSectionsGenRec() {
   

   //Determine acceptance and purity
   //Create 3D Histos with bins in Q2, X and tau_zQ
   
}


