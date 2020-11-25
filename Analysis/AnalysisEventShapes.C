// (c) MPI 2020
#include "AnalysisEventShapes.h"
 
// // C++ includes
// #include <map>
// #include <set>
// #include <vector>

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
// #include "H1PhysUtils/H1MakeKine.h"

#include "JetTools.h"         // JetsAtHighQ2
#include "H1Mods/H1PartSelTrack.h"
#include "H1JetFinder/H1EventShape.h"
//#include "UsefulTools.h"
using namespace std;



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

   // event weight
   fGen.wgt      = gH1Calc->Weight()->GetWeightGen(); // gH1Calc->GetFloatVariable(Weight_Weight);

   // basic DIS quantites
   fGen.Q2       = gH1Calc->Kine()->GetQ2Gen();
   fGen.Y        = gH1Calc->Kine()->GetYGen();
   fGen.X        = gH1Calc->Kine()->GetXGen();

   // event shapes
   fGen.tau_zQ   = 3.1;
   fGen.tau1a    = 3.1;
   fGen.tau_zP   = 3.1;

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

   // event weight
   fRec.wgt      = gH1Calc->Weight()->GetWeight();
   
   // basic DIS observables
   fRec.Q2       = gH1Calc->Kine()->GetQ2();
   fRec.Y        = gH1Calc->Kine()->GetY();
   fRec.X        = gH1Calc->Kine()->GetX();
   
   // event shapes
   fRec.tau_zQ   = 3.1;
   fRec.tau1a    = 3.1;
   fRec.tau_zP   = 3.1;

}



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoControlPlotsGen()
//!
//!  This function is called by the main program
//!  during the event loop
//!
void AnalysisEventShapes::DoControlPlotsGen() {

   const double wgt = fGen.wgt;

   TLorentzVector ebeam, pbeam;
   H1BoostedJets::Instance()->BeamVectors(&pbeam, &ebeam); 
   Float_t EmpzGen = gH1Calc->Fs()->GetEmpzGen(); // = Delta = E - pZ
   double SigmaGen = 0; //SigmaGen ( = sum_{hadron final state} (E-pZ) )

   TLorentzVector ScatElecGen = gH1Calc->Elec()->GetFirstElectronGen();
   TLorentzVector ElecInitGen( 0., 0., -EmpzGen/2, EmpzGen/2);
   vector<TLorentzVector> particlelist;
   vector<H1PartMC*> hadronarray = to_vector<H1PartMC*>(H1BoostedJets::Instance()->GetHadronArray());
   for (long unsigned int ipart=0; ipart<hadronarray.size(); ++ipart){
     H1PartMC* part = static_cast<H1PartMC*>(hadronarray[ipart]);
     if ( !(part->GetE()>0.) ) continue; // reject zeroed-out part cands (badly measured ones, iron muons...)
     particlelist.push_back(TLorentzVector(part->GetPx(),part->GetPy(), part->GetPz(), part->GetE()) );
     SigmaGen += (part->GetE() - std::fabs(part->GetPz()) );
   }
   TLorentzVector virtualphotonGen = ebeam - ScatElecGen;
   TLorentzVector photonISigmaGen = ElecInitGen - ScatElecGen;
   
  //elelctron method
   Float_t YeGen  = GetKine_YeGen();
   Float_t Q2eGen = GetKine_Q2eGen();
   Float_t XeGen = GetKine_XeGen();
   /*
  //hadron method
   Float_t YhGen  = GetKine_YhGen();
   Float_t Q2hGen = GetKine_Q2hGen();
   Float_t XhGen = gH1Calc->Kine()->GetXhGen();
  //eSigma method
   Float_t YesGen = GetKine_YesGen();
   Float_t Q2esGen = GetKine_Q2esGen();
   Float_t XesGen = GetKine_XesGen();
   */
  //Sigma method
   Float_t YsGen = GetKine_YsGen();
   Float_t Q2sGen = GetKine_Q2sGen();
   Float_t XsGen = GetKine_XsGen();
  //ISigma method
   Float_t YisGen = SigmaGen / ( SigmaGen + ScatElecGen.E()*( 1-TMath::Cos(ScatElecGen.Theta()) ) );
   Float_t Q2isGen = TMath::Power( ScatElecGen.E()*TMath::Sin( ScatElecGen.Theta() ) , 2 ) / ( 1-YisGen );
   Float_t XisGen = ( ScatElecGen.E() / pbeam.E() ) * ( TMath::Power( TMath::Cos(ScatElecGen.Theta()/2) , 2 ) / YsGen );

   Float_t Ymin = 0.2;
   Float_t Ymax = 0.7;
   Float_t Q2min = 150;
   Float_t Q2max = 20000;

   H1Boost myboostelectronGen(2*XeGen*pbeam,virtualphotonGen,ebeam,-1.*pbeam);
   H1Boost myboostISGen(2*XisGen*pbeam, photonISigmaGen, ebeam, -1.*pbeam);
   H1Boost test_isigma_boost = CalcBoost(Q2sGen, YsGen, XisGen, ScatElecGen.Phi(), pbeam.E());
   H1Boost test_el_boost = CalcBoost(Q2eGen, YeGen, XeGen, ScatElecGen.Phi(), pbeam.E());
   if ( YeGen > Ymin && YeGen < Ymax ) {
     if ( Q2eGen > Q2min && Q2eGen < Q2max ){
       vector<TLorentzVector> boosted_eGen = CalculateEventShape_tauzQ("Gen_tau_zQ_el", particlelist, myboostelectronGen, Q2eGen);
       vector<TLorentzVector> myboost_eGen = CalculateEventShape_tauzQ("Gen_tau_zQ_el_myboost", particlelist, test_el_boost, Q2eGen);
       PlotKinematicVaribles( "Gen_tau_zQ_el_ycut", Q2eGen, XeGen, YeGen);
       ClassicalEventShapes("eventshapes_eGen", boosted_eGen);
     }
   }
   if ( YsGen > Ymin && YsGen < Ymax ) {
     if ( Q2sGen > Q2min && Q2sGen < Q2max ){
       vector<TLorentzVector> boosted_isGen = CalculateEventShape_tauzQ("Gen_tau_zQ_iSigma", particlelist, myboostISGen, Q2sGen);
       vector<TLorentzVector> myboost_isGen = CalculateEventShape_tauzQ("Gen_tau_zQ_iSigma_myboost", particlelist, test_isigma_boost, Q2sGen);
       PlotKinematicVaribles( "Gen_tau_zQ_iSigma_ycut", Q2sGen, XsGen, YsGen);
       ClassicalEventShapes("eventshapes_isGen", boosted_isGen);
     }
   }
   if ( YsGen > Ymin && YsGen < 1 ) {
     if ( Q2sGen > Q2min && Q2sGen < Q2max ){
       vector<TLorentzVector> boosted_isGen = CalculateEventShape_tauzQ("Gen_tau_zQ_iSigma_highy", particlelist, myboostISGen, Q2sGen);
       PlotKinematicVaribles( "Gen_tau_zQ_iSigma_highy", Q2sGen, XsGen, YsGen);
     }
   }
   if ( YsGen > 0. && YsGen < Ymax ) {
     if ( Q2sGen > Q2min && Q2sGen < Q2max ){
       vector<TLorentzVector> boosted_isGen = CalculateEventShape_tauzQ("Gen_tau_zQ_iSigma_lowy", particlelist, myboostISGen, Q2sGen);
       PlotKinematicVaribles( "Gen_tau_zQ_iSigma_lowy", Q2sGen, XsGen, YsGen);
     }
   }
   //no kinematic cuts
   //vector<TLorentzVector> boosted_eGen_nocut = CalculateEventShape_tauzQ("Gen_tau_zQ_el", particlelist, myboostelectronGen, Q2eGen);
   //vector<TLorentzVector> boosted_isGen_nocut = CalculateEventShape_tauzQ("Gen_tau_zQ_iSigma", particlelist, myboostISGen, Q2sGen);


   TLorentzVector breitpart = myboostelectronGen.Boost(ScatElecGen);
   /*
   cout<<"scat elec"<<endl;
   breitpart.Print();
   cout<<"init elec"<<endl;
   myboostelectronGen.Boost(ebeam).Print();
   cout<<"init proton"<<endl;
   myboostelectronGen.Boost(pbeam).Print();
   cout<<"scat proton"<<endl;
   myboostelectronGen.Boost(H1BoostedJets::Instance()->FillHadronArray()).Print();
   */
}



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoControlPlotsRec()
//!
//!  This function is called by the main program
//!  during the event loop
//!
void AnalysisEventShapes::DoControlPlotsRec() {
   
   // --- event weight
   double wgt = fRec.wgt;

   // fill set of NC DIS control plots
   {
   double Q2 = gH1Calc->Kine()->GetQ2es();
   double Y  = gH1Calc->Kine()->GetQ2es();
   double X  = gH1Calc->Kine()->GetQ2es();
   
   FillBasicNCDISHists("DISEvent", Q2, Y, X, wgt);
   
   if ( Q2 > 800 ) 
      FillBasicNCDISHists("DISEvent highQ2", Q2, Y, X, wgt);
   if ( Y <0.3 ) 
      FillBasicNCDISHists("DISEvent low-Y", Q2, Y, X, wgt);

   if ( X >0.2 ) 
      FillBasicNCDISHists("DISEvent high-X", Q2, Y, X, wgt);
   
   }
   // this is one histmanager, and will become one TDirectory in the file later
   H2020HistManager& hm       = HistMaster::Instance()->GetHistManager("DISEvent");
   // this is is another histmanager
   H2020HistManager& hmtracks = HistMaster::Instance()->GetHistManager("NCDISTracks");

   H2020HistManager& hmGenRatio = HistMaster::Instance()->GetHistManager("GenRatio");

   H2020HistManager& hmtau_zQ_el = HistMaster::Instance()->GetHistManager("tau_zQ_el");

   H2020HistManager& hmtau_zQ_had = HistMaster::Instance()->GetHistManager("tau_zQ_had");

   H2020HistManager& hmRecBoost = HistMaster::Instance()->GetHistManager("RecBoost");

   // put some histogram inside
   hmtracks.Get<TH1D>("N tracks BoostedJets",   ";N tracks Boosedjets;entries", 50, 0, 100)  -> Fill(  H1BoostedJets::Instance()->GetHFSArray()->GetEntries()   , wgt ); 
   hmtracks.Get<TH1D>("N tracks BoostedJets_ly",";N tracks Boosedjets;entries", 50, 0, 100)  -> Fill(  H1BoostedJets::Instance()->GetHFSArray()->GetEntries()   , wgt ); 
   // put tracks into TDirectory: "DISEvents"
   // H2020HistManager& hmtracks2 = HistMaster::Instance()->GetHistManager("NCDISTracks uniquename1","DISEvent");
   // // put some histogram inside
   // hmtracks2.Get<TH1D>("N tracks BoostedJets",   ";N tracks Boosedjets;entries", 50, 0, 100)  -> Fill(  H1BoostedJets::Instance()->GetHFSArray()->GetEntries()   , wgt ); 
   // hmtracks2.Get<TH1D>("N tracks BoostedJets_ly",";N tracks Boosedjets;entries", 50, 0, 100)  -> Fill(  H1BoostedJets::Instance()->GetHFSArray()->GetEntries()   , wgt ); 

   // this is how you make a histogram and fill it...
   
   /*
   double tau_zQ = 0.1;
   hm.Get<TH1D>("tau_zQ",    "tau_zQ;#tau_{zQ};events",  70, -0.3, 1.1)  -> Fill(tau_zQ, wgt); 
   */

   // --------------- general kinematics -------------------
   TLorentzVector ebeam, pbeam;
   H1BoostedJets::Instance()->BeamVectors(&pbeam, &ebeam);
   TLorentzVector ScatElec = gH1Calc->Elec()->GetFirstElectron();
   TLorentzVector HFS = gH1Calc->Fs()->GetFullHadFS();
   TLorentzVector ScatElecGen = gH1Calc->Elec()->GetFirstElectronGen();

   Float_t Empz = gH1Calc->Fs()->GetEmpz(); // = Delta, E-pz for the full final state
   Float_t EmpzGen = gH1Calc->Fs()->GetEmpzGen();
   Float_t Sigma = gH1Calc->Fs()->GetFullHadFsEmpz(); // E-pz for the full hadron final state
   
   Float_t PtMiss = gH1Calc->Fs()->GetPtMiss();
   Float_t prim_zvtx = gH1Calc->Vertex()->GetHatPrimaryVertexZ();
   Float_t optimal_zvtx = gH1Calc->Vertex()->GetZ(H1CalcVertex::vtOptimalNC);
   Bool_t NuclInt = gH1Calc->Vertex()->GetOptimalNCVertexType() == H1CalcVertex::optimalDTNV;
   Float_t gammah_degree = (gH1Calc->Fs()->GetGamma()) * (180.0/TMath::Pi());


   TLorentzVector fullHFS = H1BoostedJets::Instance()->FillHFSArray();
   TObjArray* detectorlevel = H1BoostedJets::Instance()->GetHFSArray();

   vector<TLorentzVector> particlelist;
   vector<H1PartCand*> hfsarray = to_vector<H1PartCand*>(H1BoostedJets::Instance()->GetHFSArray());
   vector<H1PartMC*> hadronarray = to_vector<H1PartMC*>(H1BoostedJets::Instance()->GetHadronArray());
   for (long unsigned int ipart=0; ipart<hfsarray.size(); ++ipart){
     H1PartCand* part = static_cast<H1PartCand*>(hfsarray[ipart]);
     if (part->IsScatElec()) continue;
     if ( !(part->GetE()>0.) ) continue; // reject zeroed-out part cands (badly measured ones, iron muons...)
     particlelist.push_back(TLorentzVector(part->GetPx(),part->GetPy(), part->GetPz(), part->GetE()) );
   }

   //get SigmaGen ( = sum_{hadron final state} (E-pZ) )
   double SigmaGen = 0;
   for (long unsigned int ipart=0; ipart<hadronarray.size(); ++ipart){
     H1PartMC* part = static_cast<H1PartMC*>(hadronarray[ipart]);
     //if (part->IsScatElec()) continue;  // exclude the scattered electron
     if ( !(part->GetE()>0.) ) continue; // reject zeroed-out part cands (badly measured ones, iron muons...)
     SigmaGen += (part->GetE() - part->GetPz() );
   }
   
   TLorentzVector virtualphoton = ebeam - ScatElec;
   TLorentzVector virtualphotonGen = ebeam - ScatElecGen;
   hmGenRatio.Get<TH1D>("Empz_ly", 50, 30, 75., "Empz [GeV]", "Entries")  -> Fill(  Empz   , wgt );
   /*
   Float_t myelectronenergy = ScatElec.Dot(fullHFS)/(ScatElec.E()+ScatElec.Pz()+fullHFS.E()+fullHFS.Pz());
   hmGenRatio.Get<TH1D>("Empz_ly", 50, 30, 75., "Empz [GeV]", "Entries")  -> Fill(  Empz   , wgt );
   hmGenRatio.Get<TH1D>("MyElecEnergy_ly", 50, 30, 75., "2*myelectronenergy [GeV]", "Entries")  -> Fill(  2*myelectronenergy   , wgt );
   hmGenRatio.Get<TH1D>("MyElecEnergyRatio_ly", 50, 0., 5., "2*myelectronenergy/Empz [GeV]", "Entries")  -> Fill(  2*myelectronenergy/Empz   , wgt );
   hmGenRatio.Get<TH1D>("MyElecEnergyRatio", 50, 0., 5., "2*myelectronenergy/Empz [GeV]", "Entries")  -> Fill(  2*myelectronenergy/Empz   , wgt );
   */

   TLorentzVector Pplusq = pbeam + virtualphoton;
   Float_t W   = Pplusq.Mag();
   Float_t Q2  = GetKine_Q2();
   Float_t X   = GetKine_X();
   Float_t Y   = GetKine_Y();

//Detector level
 //elelctron method
   Float_t Ye  = GetKine_Ye();
   Float_t Q2e = GetKine_Q2e();
   Float_t Xe = GetKine_Xe();
 //hadron method
   Float_t Yh  = GetKine_Yh();
   Float_t Q2h = GetKine_Q2h();
   Float_t Xh = gH1Calc->Kine()->GetXh();
   TLorentzVector virtualphoton_hadron = fullHFS - Xh*pbeam;
 //eSigma method
   Float_t Yes = GetKine_Yes();
   Float_t Q2es = GetKine_Q2es();
   Float_t Xes = GetKine_Xes();
 //Sigma method 
   Float_t Ys = GetKine_Ys();
   Float_t Q2s = GetKine_Q2s();
   Float_t Xs = GetKine_Xs();
 //ISigma method
   Float_t Yis = Sigma / ( Sigma + ScatElec.E()*( 1-TMath::Cos(ScatElec.Theta()) ) );
   Float_t Q2is = TMath::Power( ScatElec.E()*TMath::Sin( ScatElec.Theta() ) , 2 ) / ( 1-Yis );
   Float_t Xis = ( ScatElec.E() / pbeam.E() ) * ( TMath::Power( TMath::Cos(ScatElec.Theta()/2) , 2 ) / Yis );

//Particle level
 //elelctron method
   Float_t YeGen  = GetKine_YeGen();
   Float_t Q2eGen = GetKine_Q2eGen();
   Float_t XeGen = GetKine_XeGen();
 //hadron method
   Float_t YhGen  = GetKine_YhGen();
   Float_t Q2hGen = GetKine_Q2hGen();
   Float_t XhGen = gH1Calc->Kine()->GetXhGen();
 //eSigma method
   Float_t YesGen = GetKine_YesGen();
   Float_t Q2esGen = GetKine_Q2esGen();
   Float_t XesGen = GetKine_XesGen();
 //Sigma method
   Float_t YsGen = GetKine_YsGen();
   Float_t Q2sGen = GetKine_Q2sGen();
   Float_t XsGen = GetKine_XsGen();
 //ISigma method
   Float_t YisGen = SigmaGen / ( SigmaGen + ScatElecGen.E()*( 1-TMath::Cos(ScatElecGen.Theta()) ) );
   Float_t Q2isGen = TMath::Power( ScatElecGen.E()*TMath::Sin( ScatElecGen.Theta() ) , 2 ) / ( 1-YisGen );  
   Float_t XisGen = ( ScatElecGen.E() / pbeam.E() ) * ( TMath::Power( TMath::Cos(ScatElecGen.Theta()/2) , 2 ) / YsGen );

   hmGenRatio.Get<TH1D>("RatioYeGen_ly", 50, -0.1, 2., "Ye / YeGen", "Entries")  -> Fill(  Ye/YeGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioYeGen", 50, -0.1, 2., "Ye / YeGen", "Entries")  -> Fill(  Ye/YeGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioQ2eGen_ly", 50, -0.1, 2., "Q2e / Q2eGen", "Entries")  -> Fill(  Q2e/Q2eGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioQ2eGen", 50, -0.1, 2., "Q2e / Q2eGen", "Entries")  -> Fill(  Q2e/Q2eGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioXeGen_ly", 50, -0.1, 2., "Xe / XeGen", "Entries")  -> Fill(  Xe/XeGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioXeGen", 50, -0.1, 2., "Xe / XeGen", "Entries")  -> Fill(  Xe/XeGen   , wgt );

   //Ratio e to isGen

   hmGenRatio.Get<TH1D>("RatioYeToSigma_ly", 50, -0.1, 2., "Ye / Ys", "Entries")  -> Fill(  Ye/Ys   , wgt );
   hmGenRatio.Get<TH1D>("RatioYeToSigma", 50, -0.1, 2., "Ye / Ys", "Entries")  -> Fill(  Ye/Ys   , wgt );
   hmGenRatio.Get<TH1D>("RatioQ2eToSigma_ly", 50, -0.1, 2., "Q2e / Q2s", "Entries")  -> Fill(  Q2e/Q2s   , wgt );
   hmGenRatio.Get<TH1D>("RatioQ2eToSigma", 50, -0.1, 2., "Q2e / Q2s", "Entries")  -> Fill(  Q2e/Q2s   , wgt );
   hmGenRatio.Get<TH1D>("RatioXeToSigma_ly", 50, -0.1, 2., "Xe / Xs", "Entries")  -> Fill(  Xe/Xs   , wgt );
   hmGenRatio.Get<TH1D>("RatioXeToSigma", 50, -0.1, 2., "Xe / Xs", "Entries")  -> Fill(  Xe/Xs   , wgt );

   hmGenRatio.Get<TH1D>("RatioYhGen_ly", 50, -0.1, 2., "Yh / YhGen", "Entries")  -> Fill(  Yh/YhGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioYhGen", 50, -0.1, 2., "Yh / YhGen", "Entries")  -> Fill(  Yh/YhGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioQ2hGen_ly", 50, -0.1, 2., "Q2h / Q2hGen", "Entries")  -> Fill(  Q2h/Q2hGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioQ2hGen", 50, -0.1, 2., "Q2h / Q2hGen", "Entries")  -> Fill(  Q2h/Q2hGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioXhGen_ly", 50, -0.1, 2., "Xh / XhGen", "Entries")  -> Fill(  Xh/XhGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioXhGen", 50, -0.1, 2., "Xh / XhGen", "Entries")  -> Fill(  Xh/XhGen   , wgt );

   hmGenRatio.Get<TH1D>("RatioYesGen_ly", 50, -0.1, 2., "Yes / YesGen", "Entries")  -> Fill(  Yes/YesGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioYesGen", 50, -0.1, 2., "Yes / YesGen", "Entries")  -> Fill(  Yes/YesGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioQ2esGen_ly", 50, -0.1, 2., "Q2es / Q2esGen", "Entries")  -> Fill(  Q2es/Q2esGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioQ2esGen", 50, -0.1, 2., "Q2es / Q2esGen", "Entries")  -> Fill(  Q2es/Q2esGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioXesGen_ly", 50, -0.1, 2., "Xes / XesGen", "Entries")  -> Fill(  Xes/XesGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioXesGen", 50, -0.1, 2., "Xes / XesGen", "Entries")  -> Fill(  Xes/XesGen   , wgt );

   hmGenRatio.Get<TH1D>("RatioYsGen_ly", 50, -0.1, 2., "Ys / YsGen", "Entries")  -> Fill(  Ys/YsGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioYsGen", 50, -0.1, 2., "Ys / YsGen", "Entries")  -> Fill(  Ys/YsGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioQ2sGen_ly", 50, -0.1, 2., "Q2s / Q2sGen", "Entries")  -> Fill(  Q2s/Q2sGen   , wgt );   
   hmGenRatio.Get<TH1D>("RatioQ2sGen", 50, -0.1, 2., "Q2s / Q2sGen", "Entries")  -> Fill(  Q2s/Q2sGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioXsGen_ly", 50, -0.1, 2., "Xs / XsGen", "Entries")  -> Fill(  Xs/XsGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioXsGen", 50, -0.1, 2., "Xs / XsGen", "Entries")  -> Fill(  Xs/XsGen   , wgt );

   hmGenRatio.Get<TH1D>("RatioYisGen_ly", 50, -0.1, 2., "Yis / YisGen", "Entries")  -> Fill(  Yis/YisGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioYisGen", 50, -0.1, 2., "Yis / YisGen", "Entries")  -> Fill(  Yis/YisGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioQ2isGen_ly", 50, -0.1, 2., "Q2is / Q2isGen", "Entries")  -> Fill(  Q2is/Q2isGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioQ2isGen", 50, -0.1, 2., "Q2is / Q2isGen", "Entries")  -> Fill(  Q2is/Q2isGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioXisGen_ly", 50, -0.1, 2., "Xis / XisGen", "Entries")  -> Fill(  Xis/XisGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioXisGen", 50, -0.1, 2., "Xis / XisGen", "Entries")  -> Fill(  Xis/XisGen   , wgt );

   hmGenRatio.Get<TH1D>("RatioXisToSigma_ly", 50, -0.1, 2., "Xis / Xs", "Entries")  -> Fill(  Xis/Xs   , wgt );
   hmGenRatio.Get<TH1D>("RatioXisToSigma", 50, -0.1, 2., "Xis / Xs", "Entries")  -> Fill(  Xis/Xs   , wgt );

   hmGenRatio.Get<TH1D>("RatioEmpzGen_ly", 50, -0.1, 2., "Empz / EmpzGen", "Entries")  -> Fill(  Empz/EmpzGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioEmpzGen", 50, -0.1, 2., "Empz / EmpzGen", "Entries")  -> Fill(  Empz/EmpzGen   , wgt );
   hmGenRatio.Get<TH1D>("Sigma", 50, 0, 50, "Sigma [GeV]", "Entries") -> Fill ( Sigma, wgt);
   hmGenRatio.Get<TH1D>("Sigma_ly", 50, 0, 50, "Sigma [GeV]", "Entries") -> Fill ( Sigma, wgt);
   hmGenRatio.Get<TH1D>("SigmaGen", 50, 0, 50, "SigmaGen [GeV]", "Entries") -> Fill ( SigmaGen, wgt);
   hmGenRatio.Get<TH1D>("SigmaGen_ly", 50, 0, 50, "SigmaGen [GeV]", "Entries") -> Fill ( SigmaGen, wgt);
   hmGenRatio.Get<TH1D>("RatioSigmaGen_ly", 50, -0.1, 2., "Sigma / SigmaGen", "Entries")  -> Fill(  Sigma/SigmaGen   , wgt );
   hmGenRatio.Get<TH1D>("RatioSigmaGen", 50, -0.1, 2., "Sigma / SigmaGen", "Entries")  -> Fill(  Sigma/SigmaGen   , wgt );  

   if ( Sigma/SigmaGen > 1.3 ) {
     hmGenRatio.Get<TH1D>("SigmaRightTailXe", 50, 0., 1., "Xe, #Sigma/#Sigma_{Gen} > 1.3", "Entries")  -> Fill(  Xe   , wgt );
     hmGenRatio.Get<TH1D>("SigmaRightTailXe_lxy", 80, 0.001, 1., "Xe, #Sigma/#Sigma_{Gen} > 1.3", "Entries")  -> Fill(  Xe   , wgt );
     hmGenRatio.Get<TH1D>("SigmaRightTailYe", 50, 0., 1., "Ye, #Sigma/#Sigma_{Gen} > 1.3", "Entries")  -> Fill(  Ye   , wgt );
     hmGenRatio.Get<TH1D>("SigmaRightTailYe_lxy", 80, 0.001, 1., "Ye, #Sigma/#Sigma_{Gen} > 1.3", "Entries")  -> Fill(  Ye   , wgt );
     hmGenRatio.Get<TH1D>("SigmaRightTailXis", 50, 0., 1., "Xis, #Sigma/#Sigma_{Gen} > 1.3", "Entries")  -> Fill(  Xis   , wgt );
     hmGenRatio.Get<TH1D>("SigmaRightTailXis_lxy", 80, 0.001, 1., "Xis, #Sigma/#Sigma_{Gen} > 1.3", "Entries")  -> Fill(  Xis   , wgt );
     hmGenRatio.Get<TH1D>("SigmaRightTailYis", 50, 0., 1., "Yis, #Sigma/#Sigma_{Gen} > 1.3", "Entries")  -> Fill(  Yis   , wgt );
     hmGenRatio.Get<TH1D>("SigmaRightTailYis_lxy", 80, 0.001, 1., "Yis, #Sigma/#Sigma_{Gen} > 1.3", "Entries")  -> Fill(  Yis   , wgt );
   }

   H2020HistManager& hmSelectSigma = HistMaster::Instance()->GetHistManager("SelectSigma");

   if ( Ys > 0.1 && Ys < 0.7 ) {
     if ( Xis > 0. && Xis < 1. ) {
       if ( Q2s > 150 && Q2s < 40000 ) {
	 hmSelectSigma.Get<TH1D>("SigmaGenSelected_ly", 50, -0.1, 2., "Sigma / SigmaGen, el method cuts", "Entries")  -> Fill(  Sigma/SigmaGen   , wgt );
	 hmSelectSigma.Get<TH1D>("SigmaGenSelected", 50, -0.1, 2., "Sigma / SigmaGen, el method cuts", "Entries")  -> Fill(  Sigma/SigmaGen   , wgt );
       }
     }
   }
   
   TLorentzVector ElecInit( 0., 0., -Empz/2, Empz/2); //initial electron, mass neglecte
   TLorentzVector photonISigma = ElecInit - ScatElec; //virtual photon with final state electron from detector
   
   //  TLorentzVector boostISigma = BoostToBreitFrame( Q2is, Yis, Empz/2, ScatElec.Phi() ); //final state elctron calculated

   
   /*
   //Boost vector defined as b = 2xP + q (=0 is Breit frame)
   //Different possibilities to calculate b:
   // 1) b1 = 3 e0 - 3 e' +2 h
   // 2) b2 = 3 xP -h
   // 3) b3 = xP -2 ( e0 -e' ) + h
   // 4) b4 = 2 xP + e0 - e'
   //Write b = vect?1 + vect?2: 

   TLorentzVector vect11 = 3*virtualphoton;
   TLorentzVector vect12 = 2*fullHFS;
   TLorentzVector vect21 = 3*Xes*pbeam;
   TLorentzVector vect22 = -1*fullHFS;
   TLorentzVector vect31 = Xes*pbeam + fullHFS;
   TLorentzVector vect32 = 2*virtualphoton;
   TLorentzVector vect41 = 2*Xes*pbeam;
   TLorentzVector vect42 = virtualphoton;
   */

   //Different options for the boost to Breit frame
   //electron method
   H1Boost myboostelectron(2*Xe*pbeam,virtualphoton,ebeam,-1.*pbeam);
   //hadron method
   H1Boost myboosthadron(2*Xh*pbeam,virtualphoton_hadron,ebeam,-1.*pbeam);
   
   /*
   //eSigma method with b1
   H1Boost myboosteS1(vect11, vect12, ebeam,-1.*pbeam);
   //eSigma method with b2
   H1Boost myboosteS2(vect21, vect22, ebeam,-1.*pbeam);
   //eSigma method with b3
   H1Boost myboosteS3(vect31, vect32, ebeam,-1.*pbeam);
   //eSigma method with b4
   H1Boost myboosteS4(vect41, vect42, ebeam,-1.*pbeam);
   */
   //ISigma method
   H1Boost myboostIS(2*Xis*pbeam, photonISigma, ebeam,-1.*pbeam);
   
   double_t Ymin = 0.2;
   double_t Ymax = 0.7;
   double_t Q2min = 150;
   double_t Q2max = 20000;
   
     if ( Ye > Ymin && Ye < Ymax ) {
       if ( Xe > 0. && Xe < 1. ) {
	 if ( Q2e > Q2min && Q2e < Q2max ) {
	   vector<TLorentzVector> test = CalculateEventShape_tauzQ ("tau_zQ_el_new", particlelist, myboostelectron, Q2e);
	   vector<TLorentzVector> boosted_e = BoostParticleArray("tau_zQ_el", hfsarray, myboostelectron, Q2e, Xe, Ye);
	   ClassicalEventShapes("testshapes_e", boosted_e);
	   //	   hmGenRatio.Get<TH1D>("SigmaGenSelected_ly", 50, -0.1, 2., "Sigma / SigmaGen, el method cuts", "Entries")  -> Fill(  Sigma/SigmaGen   , wgt );
	   //hmGenRatio.Get<TH1D>("SigmaGenSelected", 50, -0.1, 2., "Sigma / SigmaGen, el method cuts", "Entries")  -> Fill(  Sigma/SigmaGen   , wgt );
	 }
       }
     }
     if ( Ys > Ymin && Ys < Ymax ) {
       if ( Xis > 0. && Xis < 1. ) {
	 if ( Q2s > Q2min && Q2s < Q2max ) {
	   vector<TLorentzVector> test_is = CalculateEventShape_tauzQ ("tau_zQ_iSigma_new", particlelist, myboostIS, Q2s);
	   vector<TLorentzVector> boosted_IS = BoostParticleArray("tau_zQ_iSigma", hfsarray, myboostIS, Q2is, Xis, Ys);
	   ClassicalEventShapes("testshapes_is", boosted_IS);
	 }
       }
     }

/*
   vector<TLorentzVector> boosted_h = BoostParticleArray("tau_zQ_had", hfsarray1, myboosthadron, Q2h, Xh, Yh);
   ClassicalEventShapes("testshapes_h", boosted_h);
     
     vector<TLorentzVector> boosted_es1 = BoostParticleArray("tau_zQ_es1", hfsarray1, myboosteS1, Q2es, Xes, Yes);
     vector<TLorentzVector> boosted_es2 = BoostParticleArray("tau_zQ_es2", hfsarray1, myboosteS2, Q2es, Xes, Yes);
     vector<TLorentzVector> boosted_es3 = BoostParticleArray("tau_zQ_es3", hfsarray1, myboosteS3, Q2es, Xes, Yes);
     vector<TLorentzVector> boosted_es4 = BoostParticleArray("tau_zQ_es4", hfsarray1, myboosteS4, Q2es, Xes, Yes);
     */
     //vector<TLorentzVector> boosted_recboost = BoostParticleArray("tau_zQ_recboost", hadronarray, recboost, Q2isGen, XisGen, YsGen);
   
   vector<TLorentzVector> full_event_breit;
   vector<TLorentzVector> full_event_breit_hadron;
   double_t energy_breit = 0;
   double_t current_px_had  = 0;
   double_t current_py_had  = 0;
   double_t current_pz_had  = 0;
   double_t current_px_el  = 0;
   double_t current_py_el  = 0;
   double_t current_pz_el  = 0;
   for (Int_t ipart=0; ipart<detectorlevel->GetEntries(); ++ipart){
     //H1PartMC* part = allparticles[ipart];
     //H1PartCand* parttemp = (H1PartCand*) detectorlevel->At(ipart);
     H1PartCand* part = static_cast<H1PartCand*>(detectorlevel->At(ipart));
     if (part->IsScatElec()) continue;  // exclude the scattered electron
     if ( !(part->GetE()>0.) ) continue; // reject zeroed-out part cands (badly measured ones, iron muons...)
     if ( myboostelectron.IsPhysicalBoost() ) {
       TLorentzVector breitpart = myboostelectron.Boost(part->GetFourVector());
       full_event_breit.push_back(TLorentzVector(breitpart.Px(),breitpart.Py(), breitpart.Pz(), breitpart.E()) );
       //check sign of pz, which eta corresponds to current hemisphere
       //tau for both hemispheres
       if(breitpart.Eta()>0){
	 energy_breit += breitpart.E();
	 current_px_el += breitpart.Px();
	 current_py_el += breitpart.Py();
	 current_pz_el += breitpart.Pz();
       }
     }
     if ( myboosthadron.IsPhysicalBoost() ) {
       TLorentzVector breitpart = myboosthadron.Boost(part->GetFourVector());
       full_event_breit_hadron.push_back(TLorentzVector(breitpart.Px(),breitpart.Py(), breitpart.Pz(), breitpart.E()));
       if(breitpart.Eta()>0){
	 current_px_had += breitpart.Px();
	 current_py_had += breitpart.Py();
	 current_pz_had += breitpart.Pz();
       }
     }
   }

   //UsefulTools Usefultools;
   static JetTools tools;
   auto q2bins  = H2020HistManager::MakeLogBinning(50, 90, 30000.);
   auto xbjbins = H2020HistManager::MakeLogBinning(50, 0.001, 1.);
   hm.Get<TH1D>("Q2_lxy", 50, q2bins, "Q^{2} / GeV^{2}", "Entries")  -> Fill(  Q2   , wgt );
   hm.Get<TH1D>("Q2_lx",  50, q2bins, "Q^{2} / GeV^{2}", "Entries")  -> Fill(  Q2   , wgt );

   hm.Get<TH1D>("X_lxy", 50, xbjbins, "x_{Bj}", "Entries")  -> Fill(  X   , wgt );
   hm.Get<TH1D>("X_lx" , 50, xbjbins, "x_{Bj}", "Entries")  -> Fill(  X   , wgt );

   hm.Get<TH1D>("Y", 50, 0, 1, "y", "Entries")  -> Fill(  Y   , wgt );
   hm.Get<TH1D>("Ye", 50, 0, 1, "y_{e}", "Entries")  -> Fill(  Ye   , wgt );
   hm.Get<TH1D>("Yh", 50, 0, 1, "y_{h}", "Entries")  -> Fill(  Yh   , wgt );

   hm.Get<TH1D>("W", 50, 100, 300, "W / GeV", "Entries")  -> Fill(  W   , wgt );
   hm.Get<TH1D>("W_ly", 50, 100, 300, "W / GeV", "Entries")  -> Fill(  W  , wgt );

   hm.Get<TH1D>("Empz", 50, 30, 80, "E-p_{z} / GeV", "Entries")  -> Fill(  Empz   , wgt );

   hm.Get<TH1D>("Ptmiss",    25, 0, 25, "P_{t} miss / GeV", "Entries")  -> Fill(  PtMiss   , wgt );
   hm.Get<TH1D>("Ptmiss_ly", 25, 0, 25, "P_{t} miss / GeV", "Entries")  -> Fill(  PtMiss   , wgt );

   hm.Get<TH1D>("Zvertex",         50, -50, 50, "Z-Vertex / cm", "Entries")  -> Fill(  prim_zvtx   , wgt );
   hm.Get<TH1D>("Zvertex_ly",     100, -50, 50, "Z-Vertex / cm", "Entries")  -> Fill(  prim_zvtx   , wgt );
   hm.Get<TH1D>("OptimalZvertex",  50, -50, 50, "optimal NC Z-Vertex / cm", "Entries")  -> Fill(  optimal_zvtx   , wgt );
   hm.Get<TH1D>("ZvertexDiff_ly",  50, -100, 100, "optimal - primary z_{vtx} / cm", "Entries")  -> Fill(  optimal_zvtx - prim_zvtx   , wgt );
   hm.Get<TH1D>("NuclIntFound_ly",     2, -0.5, 1.5, "Nuclear Interaction", "Entries")  -> Fill(  NuclInt   , wgt );

   hm.Get<TH1D>("Gammah", 50, 0, 180, "#gamma_{h} [^{o}]", "Entries")  -> Fill(  gammah_degree   , wgt );

   // HFS
   hm.Get<TH1D>("HFS_E",     50, 0, 200,   "Energy of HFS / GeV", "Entries")  -> Fill(  HFS.E()   , wgt );
   hm.Get<TH1D>("HFS_Pz",    50, 0, 200,   "P_{z} of HFS / GeV", "Entries")  -> Fill(  HFS.Pz()   , wgt );
   hm.Get<TH1D>("HFS_Empz",  50, 0, 70,    "E - P_{z} of HFS / GeV", "Entries")  -> Fill(  HFS.E()-HFS.Pz()  , wgt );
   hm.Get<TH1D>("HFS_Pt",    50, 0, 50,    "P_{t} of HFS / GeV", "Entries")  -> Fill(  HFS.Pt()  , wgt );
   hm.Get<TH1D>("HFS_Theta", 50, 0, 180,   "#theta of HFS", "Entries")  -> Fill(   HFS.Theta()*TMath::RadToDeg()  , wgt );
   hm.Get<TH1D>("HFS_Phi",   50, -180, 180,"#phi of HFS", "Entries")  -> Fill(  HFS.Phi()*TMath::RadToDeg()    , wgt );

   // --------------------- electromagnetic/hadronic shower classifiers -------------------------
   // Double_t jetthetabins[] = {7, 10, 15, 30, 55, 80, 110, 135, 155};
   // AddTH2D("NNOutput_AllClusters_ThetaJetBins",            8, jetthetabins, 50, -1.1, 1.1, "#theta_{cluster}", "NN Output")  -> Fill(  14.1   , wgt );

   // general calibration histos
   hm.Get<TH1D>("Pth_Pte", 50, 0, 2,  "P_{T} hadrons / P_{T} electron", "Entries")  -> Fill(  GetFs_HadPtElecPtRatio()  , wgt );
   hm.Get<TH1D>("Pth_Ptda", 50, 0, 2, "P_{T} hadrons / P_{T} double angle", "Entries")  -> Fill(  GetFs_HadPtHadPtDaRatio()  , wgt );
   hm.Get<TH1D>("Pte_Ptda", 50, 0, 2, "P_{T} electron / P_{T} double angle", "Entries")  -> Fill( GetFs_ElecPtHadPtDaRatio()  , wgt );

   // total event weight
   hm.Get<TH1D>("event_weight",    100, 0, 3, "weight of the event", "Entries")  -> Fill(  wgt   , wgt );
   hm.Get<TH1D>("event_weight_ly", 100, 0, 3, "weight of the event", "Entries")  -> Fill(  wgt   , wgt );

   // background
   // todo
   // static int bg1Index = GetHistIndex("BackgroundBits1");
   // static int bg2Index = GetHistIndex("BackgroundBits2");
   // static int bg3Index = GetHistIndex("BackgroundBits3");
   // static H1IntPtr packedbg1("Ibg");
   // static H1IntPtr packedbg2("Ibgam");
   // static H1IntPtr packedbg3("Ibgfm");
   // hm.Get<TH1D>("BackgroundBits1", 10, -0.5, 9.5,  "Non-ep bkg finder 1 bit", " Entries")  -> Fill(  14.1   , wgt );
   // hm.Get<TH1D>("BackgroundBits2", 19, -0.5, 18.5, "Non-ep bkg finder 2 bit", " Entries")  -> Fill(  14.1   , wgt );
   // hm.Get<TH1D>("BackgroundBits3", 17, -0.5, 16.5, "Non-ep bkg finder 3 bit", " Entries")  -> Fill(  14.1   , wgt );


   // ------------------------- electron variables ------------------------
   Float_t e_energy = ScatElec.E();
   Float_t e_theta_degree = (ScatElec.Theta())*(180.0/TMath::Pi());
   Float_t e_phi_degree = (ScatElec.Phi())*(180.0/TMath::Pi());
   Float_t e_pt = ScatElec.Pt();
   Float_t e_zimpact = gH1Calc->Elec()->GetZimpact();

   hm.Get<TH1D>("El_Energy", 50, 0, 50, "Electron Energy / GeV", "Entries")  -> Fill( e_energy   , wgt );
   hm.Get<TH1D>("El_Theta", 50, 0, 180, "#theta_{e} [^{o}]", "Entries")  -> Fill(  e_theta_degree   , wgt );
   hm.Get<TH1D>("El_Phi", 50, -180, 180, "#phi_{e} [^{o}]","Entries")  -> Fill(  e_phi_degree   , wgt );
   hm.Get<TH1D>("El_Pt", 50, 0, 50, "Electron P_{t} / GeV", "Entries")  -> Fill(  e_pt   , wgt );
   hm.Get<TH1D>("El_Zimpact_bwd", 50, -200,   0, "Z_{impact} position of electron", "Entries")  -> Fill(  e_zimpact   , wgt );
   hm.Get<TH1D>("El_Zimpact_fwd", 50,    0, 200, "Z_{impact} position of electron", "Entries")  -> Fill(  e_zimpact   , wgt );
   hm.Get<TH1D>("El_Zimpact",     50, -200, 200, "Z_{impact} position of electron", "Entries")  -> Fill(  e_zimpact   , wgt );
   // technical variables                                                                                                                        
   Float_t e_clusterDTRAdistance = gH1Calc->Elec()->GetTrackClusDistRa();
   Float_t e_clusterDTNVdistance = gH1Calc->Elec()->GetTrackClusDistNv();
   Float_t e_clusterCIPdistance = gH1Calc->Elec()->GetCipdZ();
   Float_t e_NcipHits = gH1Calc->Elec()->GetNCipValidationHits();
   Float_t e_ptratio = gH1Calc->Elec()->GetPtRatio();

   hm.Get<TH1D>("El_PtRatio",  20, 0, 2, "Electron P_{T}(calo)/P_{T}(tracker)", "Entries")  -> Fill(  e_ptratio   , wgt );
   Float_t q2da = GetKine_Q2da();
   Float_t Eebeam = TheConstants().GetElectronBeamEnergy();
   Float_t theta_e = (gH1Calc->Elec()->GetFirstElectron()).Theta();
   Float_t ElEnDa = q2da/(4*Eebeam*TMath::Power(TMath::Cos(theta_e/2.),2));
   hm.Get<TH1D>("El_ERec_EDA", 50, 0, 2, "Electron E_{rec}/E_{da}", "Entries")  -> Fill( e_energy/ElEnDa  , wgt );

   // technical variables
   hm.Get<TH1D>("El_NCIP_hits", 8, -0.5, 7.5, "Number of electron CIP hits", "Entries")  -> Fill(  e_NcipHits   , wgt );
   hm.Get<TH1D>("El_ClusterDTRADistance_ly", 30, 0, 15, "Distance electron cluster and DTRA track", "Entries")   -> Fill(  e_clusterDTRAdistance , wgt );
   hm.Get<TH1D>("El_ClusterDTNVDistance_ly", 30, 0, 15, "Distance electron cluster and DTNV track", "Entries")   -> Fill(  e_clusterDTNVdistance , wgt );
   hm.Get<TH1D>("El_ClusterCIPDistance",  15, 0, 15, "Distance electron cluster and nearest CIP hit", "Entries") -> Fill(  e_clusterCIPdistance  , wgt );
   hm.Get<TH1D>("El_ChargeVerification_ly",    2, -0.5, 1.5, "Elec. Charge Verification", "Entries")  -> Fill(  gH1Calc->Elec()->GetBeamChargeMatch()  , wgt );




   // ---------------------- boost to Breit frame --------------------
   /*
   auto BoostBetaArray = H2020HistManager::MakeLogBinning(50, 0.08, 3.0);
   hmtau_zQ_el.Get<TH1D>("Boost_Beta_el_lxy", 50, BoostBetaArray, "Boost: #beta, el", " Entries")  -> Fill(  myboostelectron.GetBeta()   , wgt );
   hmtau_zQ_el.Get<TH1D>("Boost_Beta_el_lx",  50, BoostBetaArray, "Boost: #beta, el", " Entries")  -> Fill(  myboostelectron.GetBeta()   , wgt );
   hmtau_zQ_el.Get<TH1D>("Boost_Phi_el",  50, -3.1415, 3.1415, "Boost: #phi [^{o}], el", " Entries")  -> Fill(  myboostelectron.GetPhi()   , wgt );
   hmtau_zQ_el.Get<TH1D>("Boost_Theta_el",50, 0, 3.1415, "Boost: #theta [^{o}], el", " Entries")  -> Fill(  myboostelectron.GetTheta()   , wgt );
   
   hm.Get<TH1D>("BreitFramePx_el", 50, -15, 15, "Breit frame: P_{x}^{el}", "Entries")  -> Fill(  current_px_el   , wgt );
   hm.Get<TH1D>("BreitFramePy_el", 50, -15, 15, "Breit frame: P_{y}^{el}", "Entries")  -> Fill(  current_py_el   , wgt );
   hm.Get<TH1D>("BreitFramePt_el", 50,   0, 10, "Breit frame: P_{T}^{el}", "Entries")  -> Fill(  TMath::Sqrt(TMath::Power(current_px_el,2)+TMath::Power(current_py_el,2)), wgt );
   Float_t tau_zQ_el = 1-2*current_pz_el/Q2e;
   hm.Get<TH1D>("tau_zQ_el",    "tau_zQ;#tau_{zQ}, el;events",  70, -0.3, 1.1)  -> Fill(tau_zQ_el, wgt);
   
   hmtau_zQ_had.Get<TH1D>("Boost_Beta_had_lxy", 50, BoostBetaArray, "Boost: #beta, had", " Entries")  -> Fill(  myboosthadron.GetBeta()   , wgt );
   hmtau_zQ_had.Get<TH1D>("Boost_Beta_had_lx",  50, BoostBetaArray, "Boost: #beta, had", " Entries")  -> Fill(  myboosthadron.GetBeta()   , wgt );
   hmtau_zQ_had.Get<TH1D>("Boost_Phi_had",  50, -3.1415, 3.1415, "Boost: #phi [^{o}], had", " Entries")  -> Fill(  myboosthadron.GetPhi()   , wgt );
   hmtau_zQ_had.Get<TH1D>("Boost_Theta_had",50, 0, 3.1415, "Boost: #theta [^{o}], had", " Entries")  -> Fill(  myboosthadron.GetTheta()   , wgt );
   */
   hm.Get<TH1D>("BreitFramePx_had", 50, -15, 15, "Breit frame: P_{x}^{h}", "Entries")  -> Fill(  current_px_had   , wgt );
   hm.Get<TH1D>("BreitFramePy_had", 50, -15, 15, "Breit frame: P_{y}^{h}", "Entries")  -> Fill(  current_py_had   , wgt );
   hm.Get<TH1D>("BreitFramePt_had", 50,   0, 10, "Breit frame: P_{T}^{h}", "Entries")  -> Fill(  TMath::Sqrt(TMath::Power(current_px_had,2)+TMath::Power(current_py_had,2))   , wgt );
   Float_t tau_zQ_had = 1-2*current_pz_had/Q2h;
   hm.Get<TH1D>("tau_zQ_had",    "tau_zQ;#tau_{zQ}, had;events.",  70, -0.3, 1.1)  -> Fill(tau_zQ_had, wgt);



   // ---------------------- jet multiplicity --------------------

   Float_t Ntrackpt3 = 0;
   Float_t Ntrackpt5 = 0;
   Float_t Ntrackpt7 = 0;


   for (Int_t ipart=0; ipart<detectorlevel->GetEntries(); ++ipart){
     H1PartCand* parttemp = (H1PartCand*) detectorlevel->At(ipart);
     H1PartCand* part = static_cast<H1PartCand*>(parttemp);
     if (part->IsScatElec()) continue;
     Float_t en = part->GetE();
     Float_t Delta = -999;
     if (en>0){
       Float_t phi = part->GetPhi();
       Float_t eta = part->GetEta();
       Float_t Dphi = TMath::Pi()/180 * TMath::Abs(ScatElec.Phi() - phi);       //tools.PhiDiff(ScatElec.Phi() - phi);
       Float_t Deta = ScatElec.Eta() - eta;
       Delta = TMath::Sqrt(Dphi*Dphi + Deta*Deta);
     }
     //here add delta plot??
     //FillTH1(Delta, en*Weight, el_hfs_delta_idx);
     Float_t MultEleTracks = 0;
     if (Delta<0.05){
       MultEleTracks = 1;
     }

     hm.Get<TH1D>("Mult_Ele_Tracks",       2, -0.5, 1.5,    "Multiple electron tracks","Entries")  -> Fill(  MultEleTracks   , wgt );

     if (part->IsHFSChargedCls() || part->IsHFSUnidentifiedCls() || part->IsEm()){

       //  FillTH1(Delta, en*Weight, el_clusters_delta_idx);

       // NNOutput
       float NNOutput = part->GetEmHadSepClassifier();
       if (NNOutput < -10) continue;

       // fill histograms
       Float_t theta = part->GetTheta() * TMath::RadToDeg();
       Float_t phi = part->GetPhi() * TMath::RadToDeg();
       Float_t p = part->GetP();
       /*
       FillTH2(theta, NNOutput, Weight, nnoutput_idx);

       FillTH1(theta, p*Weight, theta_clusters_idx);
       FillTH1(phi,   p*Weight, phi_clusters_idx);
       FillTH1(p,     Weight, p_clusters_idx);
       */
     } 
     else {

       // got a track
       const H1PartSelTrack* track = part->GetIDTrack();
       //tracks from boosted jets
       if (!track){
	 // got an iron muon
	 continue;
       }

       //FillTH1(Delta, en*Weight, el_tracks_delta_idx);

       Float_t theta = track->GetTheta() * TMath::RadToDeg();
       Float_t phi = track->GetPhi() * TMath::RadToDeg();
       Float_t p = track->GetP();

       if (track->IsCentralTrk()){
	 
	 hm.Get<TH1D>("theta_ctr_tracks", 50, 0, 180, "#theta ctr tracks [^{o}]", "Entries")  -> Fill(  theta   , wgt );
	 hm.Get<TH1D>("phi_ctr_tracks",   50, 0, 180, "#phi ctr tracks [^{o}]", "Entries")  -> Fill(  phi   , wgt );
	 hm.Get<TH1D>("p_ctr_tracks_ly",  50, 0,  50, "momentum ctr tracks [GeV]", "Entries")  -> Fill(  p   , wgt );
       } else if (track->IsCombinedTrk()){
	 hm.Get<TH1D>("theta_cmb_tracks", 40, 0,  40, "#theta cmb tracks [^{o}]", "Entries")  -> Fill(  theta   , wgt );
	 hm.Get<TH1D>("phi_cmb_tracks",   50, 0, 180, "#phi cmb tracks [^{o}]", "Entries")  -> Fill(  phi   , wgt );
	 hm.Get<TH1D>("p_cmb_tracks_ly",  40, 0,  40, "momentum cmb tracks [GeV]", "Entries")  -> Fill(  p   , wgt );
       } else if (track->IsForwardTrk()){
	 hm.Get<TH1D>("theta_fwd_tracks", 40, 0,  40, "#theta fwd tracks [^{o}]", "Entries")  -> Fill(  theta   , wgt );
	 hm.Get<TH1D>("phi_fwd_tracks",   50, 0, 180, "#phi fwd tracks [^{o}]", "Entries")  -> Fill(  phi   , wgt );
	 hm.Get<TH1D>("p_fwd_tracks_ly",  40, 0,  40, "momentum fwd tracks [GeV]", "Entries")  -> Fill(  p   , wgt );
       }

       if (track->GetPt()>7){
	 Ntrackpt7 +=1;
       } else if (track->GetPt()>5){
	 Ntrackpt5 +=1;
       } else if (track->GetPt()>3){
	 Ntrackpt3 +=1;
       }

     }

   }


   hm.Get<TH1D>("N_Jets", 7, -0.5, 6.5, "Jet Multiplicity, P_{T} > 7 GeV", "Entries")  -> Fill(  Ntrackpt7   , wgt );
   hm.Get<TH1D>("N_Jets_ly", 7, -0.5, 6.5, "Jet Multiplicity, P_{T} > 7 GeV", "Entries")  -> Fill(  Ntrackpt7   , wgt );

   hm.Get<TH1D>("N_Jets_ext5", 7, -0.5, 6.5, "Jet Multiplicity, P_{T} > 5 GeV", "Entries")  -> Fill(  Ntrackpt5   , wgt );
   hm.Get<TH1D>("N_Jets_ext5_ly", 7, -0.5, 6.5, "Jet Multiplicity, P_{T} > 5 GeV", "Entries")  -> Fill(  Ntrackpt5   , wgt );

   hm.Get<TH1D>("N_Jets_ext3", 7, -0.5, 6.5, "Jet Multiplicity, P_{T} > 3 GeV", "Entries")  -> Fill(  Ntrackpt3   , wgt );
   hm.Get<TH1D>("N_Jets_ext3_ly", 7, -0.5, 6.5, "Jet Multiplicity, P_{T} > 3 GeV", "Entries")  -> Fill(  Ntrackpt3   , wgt );

   // ---------------------- multiplicities --------------------

   hm.Get<TH1D>("N_tracks"  , 51,  -0.5,  50.5, "Number of tracks"  , "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("N_em_parts", 10,  -0.5,   9.5, "Number of em-parts", "Entries")  -> Fill(  14.1   , wgt );

   static H1ShortPtr Nfwdtracks("NumForwardTracks");
   static H1ShortPtr Ncmbtracks("NumCombinedTracks");
   static H1ShortPtr Nctrtracks("NumCentralTracks");
   static H1ShortPtr NemParts("NumEmParts");
   // angle between the electron candidates
   if (*NemParts>1){
     TLorentzVector SecElec = gH1Calc->Elec()->GetElectron(H1CalcElec::SecondElectron);
     Float_t phi1 = ScatElec.Phi();
     Float_t phi2 = SecElec.Phi();
     Float_t Dphi = ScatElec.DeltaPhi(SecElec);    //TMath::Pi()/180 * Usefultools.PhiDiff(phi1 - phi2);
     Float_t Deta = ScatElec.Eta() - SecElec.Eta();
     Float_t Delta = TMath::Sqrt(Dphi*Dphi + Deta*Deta);
     Float_t Dtheta = fabs(ScatElec.Theta() - SecElec.Theta());
     Float_t aco = TMath::Cos( TMath::Pi() - TMath::Abs( phi1-phi2 ) );

     hm.Get<TH1D>("Els_Delta_Theta", 50, -90, 90,    "#theta^{1}_{e} - #theta^{2}_{e} [^{o}]", "Entries")  -> Fill(  Dtheta*180/TMath::Pi()   , wgt );
     hm.Get<TH1D>("Els_Delta_Phi",   50, 0, 180,    "#phi^{1}_{e} - #phi^{2}_{e} [^{o}]","Entries")  -> Fill(  (phi1-phi2)*180/TMath::Pi()   , wgt ); //PhiDiff??
     hm.Get<TH1D>("Els_Delta",       50, 0, 5,      "#Delta(El_{1},El_{2})","Entries")  -> Fill(  14.1   , wgt );
     hm.Get<TH1D>("El2_Energy",      50, 0, 50,     "2nd Electron Energy / GeV", "Entries")  -> Fill(  SecElec.E()   , wgt );
     hm.Get<TH1D>("Acoplanarity_ly", 50, 0.5, 1,  "Cos(|#pi - #Delta #phi|)", "Entries")  -> Fill(  aco , wgt );
   }

   hm.Get<TH1D>("N_fwd_tracks_ly", 16, -0.5, 15.5, "N forward tracks",  "Entries")  -> Fill(  *Nfwdtracks   , wgt );
   hm.Get<TH1D>("N_cmb_tracks_ly", 11, -0.5, 10.5, "N combined tracks", "Entries")  -> Fill(  *Ncmbtracks   , wgt );
   hm.Get<TH1D>("N_ctr_tracks",    51, -0.5, 50.5, "N central tracks",  "Entries")  -> Fill(  *Nctrtracks   , wgt );

   // electron-track and HFS object relations
   hm.Get<TH1D>("El_HFS_Delta",         20, 0,  0.5,      "#Delta(El_{1},PartCand)","Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("El_Tracks_Delta",      20, 0,  0.5,      "#Delta(El_{1},Track)","Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("El_Clusters_Delta",    20, 0,  0.5,      "#Delta(El_{1},Cluster)","Entries")  -> Fill(  14.1   , wgt );
   //hm.Get<TH1D>("Mult_Ele_Tracks",       2, -0.5, 1.5,    "Multiple electron tracks","Entries")  -> Fill(  MultEleTracks   , wgt );


   // ---------------------- track quantities --------------------
   /*
   hm.Get<TH1D>("theta_fwd_tracks", 40, 0,  40, "#theta fwd tracks [^{o}]", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("phi_fwd_tracks",   50, 0, 180, "#phi fwd tracks [^{o}]", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("p_fwd_tracks_ly",  40, 0,  40, "momentum fwd tracks [GeV]", "Entries")  -> Fill(  14.1   , wgt );

   hm.Get<TH1D>("theta_cmb_tracks", 40, 0,  40, "#theta cmb tracks [^{o}]", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("phi_cmb_tracks",   50, 0, 180, "#phi cmb tracks [^{o}]", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("p_cmb_tracks_ly",  40, 0,  40, "momentum cmb tracks [GeV]", "Entries")  -> Fill(  14.1   , wgt );
   
   hm.Get<TH1D>("theta_ctr_tracks", 50, 0, 180, "#theta ctr tracks [^{o}]", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("phi_ctr_tracks",   50, 0, 180, "#phi ctr tracks [^{o}]", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("p_ctr_tracks_ly",  50, 0,  50, "momentum ctr tracks [GeV]", "Entries")  -> Fill(  14.1   , wgt );
   */
   hm.Get<TH1D>("theta_clusters_ly",50, 0, 180, "energy weighted #theta clusters [^{o}]", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("phi_clusters",     50, 0, 180, "energy weighted #phi clusters [^{o}]", "Entries")  -> Fill(  14.1   , wgt );
   auto PClusterArray = H2020HistManager::MakeLogBinning(50, 0.1, 500);
   hm.Get<TH1D>("p_clusters_lxy",   50, PClusterArray, "momentum clusters [GeV]", "Entries")  -> Fill(  14.1   , wgt );


   // test the jet calibration
   auto CalFacArray = H2020HistManager::MakeLogBinning(50, 1, 30);
   hm.Get<TH1D>("cls_cal_factor_lxy", 50, CalFacArray, "cluster cal. factor", "Entries")  -> Fill(  14.1   , wgt );

   //auto Q2BinMeasArray = tools.GetQ2BinArray()  -> Fill(  14.1   , wgt );
   // Q2 in measurement bins
   //hm.Get<TH1D>("Q2InMeasBins_lxy", 6, Q2BinMeasArray, "Q^{2} [GeV^{2}]", "Entries");



   // fill some plots with tracks
   FillTrackPlots ("Tracks-BoostedJetsHFSArray", H1BoostedJets::Instance()->GetHFSArray() , wgt);


}



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoControlPlotsGenRec()
//!
//!  This function is called by the main program
//!  during the event loop
//!
void AnalysisEventShapes::DoControlPlotsGenRec() {
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
   
}





// _______________________________________________________ //
//!
//!  AnalysisEventShapes::EventLoop
//!
//!  Run the event shape analysis for one single event
//!  define observables and quantities in fEvent
// //!  
// bool AnalysisEventShapes::EventLoop() {

      
//    // --- set quantities
//    SetGenLevelQuantities(); // generator level
//    if ( fRec.bIsDISEPS && !GenOnly )  
//       SetRecLevelQuantities(); // detector level

//    // --- fill histograms
//    if ( !IsBkgMC ) { // gen-level
//       FillHist_NCDISGen();
//       FillHist_GenLevel();
//    }

//    //esanalysis.PrintEventInfo();                                                                                                                           
//    if ( fRec.bIsDISEPS && !GenOnly ) {
//       FillHist_NCDISall();
//       if ( fRec.X > 0.1 ) FillHist_NCDIS("high-x");
//       if ( fRec.Y > 0.5 ) FillHist_NCDIS("high-y");
//       FillHist_DetOnly();
//       FillHist_CrossSections();
//       FillHist_APS();
//    }
//    FillHist_Unfolding();

//    return true;
// }





// _______________________________________________________ //
//! 
//!  Fill a set of standard DIS control plots
//!  into a histmanager
//! 
void AnalysisEventShapes::FillBasicNCDISHists(const string& hmname, double Q2, double Y, double X, double weight) {

   H2020HistManager& hmNCDIS = HistMaster::Instance()->GetHistManager(hmname);
   hmNCDIS.Get<TH1D>("NCDIS 101 Q2",       "Q2;Q2 [GeV^{2}];events",  50, 100., 30000.) ->Fill(Q2, weight); 
   hmNCDIS.Get<TH1D>("NCDIS 101 Q2_lx",    "Q2;Q2 [GeV^{2}];events",  50, 10., 30000.)  ->Fill(Q2, weight); 
   hmNCDIS.Get<TH1D>("NCDIS 101 Q2_lxy", "Q2;Q2 [GeV^{2}];events",  50, 100., 30000.) ->Fill(Q2, weight); 
   hmNCDIS.Get<TH1D>("NCDIS 102 Y",       "Y;Y;events",  60, -0.1, 1.1)  ->Fill(Y, weight); 
   hmNCDIS.Get<TH1D>("NCDIS 102 Y_lx",    "Y;Y;events",  50, 0.0001, 1.) ->Fill(Y, weight); 
   hmNCDIS.Get<TH1D>("NCDIS 102 Y_lxy", "Y;Y;events",  50, 0.001, 1.)  ->Fill(Y, weight); 
   hmNCDIS.Get<TH1D>("NCDIS 103 X",       "X;X;events",  60, -0.1, 1.1)  ->Fill(X, weight); 
   hmNCDIS.Get<TH1D>("NCDIS 103 X_lx",    "X;X;events",  50, 0.0001, 1.) ->Fill(X, weight); 
   hmNCDIS.Get<TH1D>("NCDIS 103 X_lxy", "X;X;events",  50, 0.001, 1.)  ->Fill(X, weight); 

}

// _______________________________________________________ //
//! 
//!  Fill a set of standard plots for a collection of tracks
//! 
void AnalysisEventShapes::FillTrackPlots(const string& hmname, const std::vector<H1Part*>& parts , double weight) {

   H2020HistManager& hm = HistMaster::Instance()->GetHistManager(hmname);
   hm.Get<TH1D>("Tracks 401 N",      "N;Track multiplicity;entries",  51,-0.5,101.5) -> Fill(parts.size(), weight);
   hm.Get<TH1D>("Tracks 401 N_ly",   "N;Track multiplicity;entries",  51,-0.5,101.5) -> Fill(parts.size(), weight);
   TLorentzVector sum;
   for ( H1Part* part : parts ) {
      TLorentzVector lpart = part->GetFourVector();
      sum += lpart;
      hm.Get<TH1D>("Tracks 401 pt",   ";Track p_{T} [GeV];entries",  40,0,6)       -> Fill(part->GetPt(), weight);
      hm.Get<TH1D>("Tracks 401 pt_lx",   ";Track p_{T} [GeV];entries",  40,0.1,10)       -> Fill(part->GetPt(), weight);
      hm.Get<TH1D>("Tracks 401 pt_lxy",";Track p_{T} [GeV];entries",  40,0.1,20)       -> Fill(part->GetPt(), weight);
      hm.Get<TH1D>("Tracks 401 px",   ";Track p_{X} [GeV];entries",  40,0,50)       -> Fill(lpart.Px(), weight);
      hm.Get<TH1D>("Tracks 401 py",   ";Track p_{Y} [GeV];entries",  40,0,50)       -> Fill(lpart.Py(), weight);
      hm.Get<TH1D>("Tracks 401 pz",   ";Track p_{Z} [GeV];entries",  40,0,4)       -> Fill(part->GetPz(), weight);
      hm.Get<TH1D>("Tracks 401 pz2",   ";Track p_{Z} [GeV];entries",  40,0,10)       -> Fill(part->GetPz(), weight);
      hm.Get<TH1D>("Tracks 401 pz_lx",   ";Track p_{Z} [GeV];entries",  40,0.1,10)       -> Fill(part->GetPz(), weight);
      hm.Get<TH1D>("Tracks 401 pz_lxy",";Track p_{Z} [GeV];entries",  40,0.1,20)       -> Fill(part->GetPz(), weight);
      hm.Get<TH1D>("Tracks 401 eta",  ";Track #eta;entries",         40,-5,5)       -> Fill(part->GetEta(), weight);
      hm.Get<TH1D>("Tracks 401 theta",";Track #theta;entries",       40,-M_PI,M_PI) -> Fill(part->GetTheta(), weight);
      hm.Get<TH1D>("Tracks 401 phi",  ";Track #phi;entries",         40,-M_PI,M_PI) -> Fill(part->GetPhi(), weight);
   }
   hm.Get<TH1D>("Tracks 402 pt",   ";Sum-of-Tracks p_{T} [GeV];entries",  40,0,20)       -> Fill(sum.Pt(), weight);
   hm.Get<TH1D>("Tracks 402 pt",   ";Sum-of-Tracks p_{T} [GeV];entries",  40,0,80)       -> Fill(sum.Pt(), weight);
   hm.Get<TH1D>("Tracks 402 pt_lx",   ";Sum-of-Tracks p_{T} [GeV];entries",  40,0.1,80)       -> Fill(sum.Pt(), weight);
   hm.Get<TH1D>("Tracks 402 pt_lxy",   ";Sum-of-Tracks p_{T} [GeV];entries",  40,0.1,80)       -> Fill(sum.Pt(), weight);
   hm.Get<TH1D>("Tracks 402 px",   ";Sum-of-Tracks p_{X} [GeV];entries",  40,0,20)       -> Fill(sum.Px(), weight);
   hm.Get<TH1D>("Tracks 402 py",   ";Sum-of-Tracks p_{Y} [GeV];entries",  40,0,20)       -> Fill(sum.Py(), weight);
   hm.Get<TH1D>("Tracks 402 pz",   ";Sum-of-Tracks p_{Z} [GeV];entries",  40,-40,120)       -> Fill(sum.Pz(), weight);
   hm.Get<TH1D>("Tracks 402 pz_lx",   ";Sum-of-Tracks p_{Z} [GeV];entries",  40,2,400)       -> Fill(sum.Pz(), weight);
   hm.Get<TH1D>("Tracks 402 pz_lxy",   ";Sum-of-Tracks p_{Z} [GeV];entries",  40,0.1,600)       -> Fill(sum.Pz(), weight);
   hm.Get<TH1D>("Tracks 402 eta",  ";Sum-of-Tracks #eta;entries",         40,-5,5)       -> Fill(sum.Eta(), weight);
   hm.Get<TH1D>("Tracks 402 theta",";Sum-of-Tracks #theta;entries",       40,0,M_PI) -> Fill(sum.Theta(), weight);
   hm.Get<TH1D>("Tracks 402 phi",  ";Sum-of-Tracks #phi;entries",         40,-M_PI,M_PI) -> Fill(sum.Phi(), weight);

}

//! compatibility version
void AnalysisEventShapes::FillTrackPlots(const string& hm, TObjArray* parts , double weight) {
   std::vector<H1Part*> vparts;
   for ( int i = 0  ; i<parts->GetEntries() ; i++ )
      vparts.push_back(static_cast<H1Part*>(parts->At(i)));
   FillTrackPlots(hm,vparts,weight);
}

// _______________________________________________________ //
//! 
//!  Fill a set of standard detector-level H1Calculator control plots
//!  into a histmanager
//! 
void AnalysisEventShapes::FillBasicH1CalculatorRec(const string& hm, double weight) {

   H2020HistManager& hmNCDIS = HistMaster::Instance()->GetHistManager(hm);
   hmNCDIS.Get<TH1D>("h1Rec 201 Empz",       ";empz [GeV];events",  50, 100., 30000.) ->Fill(gH1Calc->Fs()->GetEmpz(), weight); 
   hmNCDIS.Get<TH1D>("h1Rec 211 Elec.E",     ";Elec E [GeV];events",  50, 5., 55.)  ->Fill(gH1Calc->Elec()->GetFirstElectron().E(), weight); 
   hmNCDIS.Get<TH1D>("h1Rec 211 Elec.Theta", ";Elec #theta (deg);events",  50, -180,180)  ->Fill(gH1Calc->Elec()->GetFirstElectron().Theta() * TMath::DegToRad(), weight); 
   hmNCDIS.Get<TH1D>("h1Rec 211 Elec.Eta",   ";Elec #phi;events",  50, -3,3.)  ->Fill(gH1Calc->Elec()->GetFirstElectron().Eta(), weight); 
   hmNCDIS.Get<TH1D>("h1Rec 211 Elec.Phi",   ";Elec #phi;events",  50, -M_PI,M_PI)  ->Fill(gH1Calc->Elec()->GetFirstElectron().Phi(), weight); 

   // and more
   // particle multiplicities, ...

}


// _______________________________________________________ //
//! 
//!  Fill a set of standard detector-level H1Calculator control plots
//!  into a histmanager
//! 
void AnalysisEventShapes::FillBasicH1CalculatorGen(const string& hm, double weight) {

   H2020HistManager& hmNCDIS = HistMaster::Instance()->GetHistManager(hm);
   hmNCDIS.Get<TH1D>("H1Gen 201 Empz",       ";empz [GeV];events",  50, 100., 30000.) ->Fill(gH1Calc->Fs()->GetEmpzGen(), weight); 
   // GetEmpzbal()
   // GetElecPtHadPtDaRatio()
   hmNCDIS.Get<TH1D>("H1Gen 211 Elec.E",     ";Elec E [GeV];events",  50, 5., 55.)  ->Fill(gH1Calc->Elec()->GetFirstElectronGen().E(), weight); 
   //hmNCDIS.Get<TH1D>("211 Elec.Theta", ";Elec #theta (deg);events",  50, -180,180)  ->Fill(gH1Calc->Elec()->Theta() * TMath::DegToRad(), weight); 
   // hmNCDIS.Get<TH1D>("211 Elec.Eta",   ";Elec #phi;events",  50, -3,3.)  ->Fill(gH1Calc->Elec()->Eta(), weight); 
   // hmNCDIS.Get<TH1D>("211 Elec.Phi",   ";Elec #phi;events",  50, -M_PI,M_PI)  ->Fill(gH1Calc->Elec()->Phi(), weight); 

   // and more
   // particle multiplicities, ...

}



//template <class T>
//vector<TLorentzVector> AnalysisEventShapes::BoostParticleArray (const string& hm, const vector<T*> partarray, H1Boost& myboost, double Q2, double X, double Y){

vector<TLorentzVector> AnalysisEventShapes::CalculateEventShape_tauzQ (const string& hm, const vector<TLorentzVector> partarray, H1Boost& myboost, double Q2) {

  H2020HistManager& hmtau_zQ = HistMaster::Instance()->GetHistManager(hm);
  const double wgt = fGen.wgt;
  vector<TLorentzVector> boostedpart_current; //do I need to reset this vector?
  TLorentzVector sum_vector(0., 0., 0., 0.);

  if ( myboost.IsPhysicalBoost() ) {
    for (long unsigned int ipart=0; ipart<partarray.size(); ++ipart){
      TLorentzVector part = static_cast<TLorentzVector>(partarray[ipart]);
      TLorentzVector breitpart = myboost.Boost(part);
      boostedpart_current.push_back(TLorentzVector(breitpart.Px(),breitpart.Py(), breitpart.Pz(), breitpart.E()) );
      if (breitpart.Eta()<0) {
        sum_vector -= breitpart;
      }
    }
    double_t tau = 1-2*sum_vector.Pz()/TMath::Sqrt(Q2);
    auto BoostBetaArray = H2020HistManager::MakeLogBinning(50, 0.08, 3.0);
    hmtau_zQ.Get<TH1D>("Boost_Beta_lxy", 50, BoostBetaArray, "Boost: #beta", "events" )  -> Fill(  myboost.GetBeta()   , wgt );
    hmtau_zQ.Get<TH1D>("Boost_Beta_lx" , 50, BoostBetaArray, "Boost: #beta", "events" )  -> Fill(  myboost.GetBeta()   , wgt );
    hmtau_zQ.Get<TH1D>("Boost_Phi", ";Boost: #phi [^{o}]; events", 60, -180, 180)  -> Fill(  myboost.GetPhi()*180/TMath::Pi()   , wgt );
    hmtau_zQ.Get<TH1D>("Boost_Theta", ";Boost: #theta [^{o}]; events", 60,    0, 180 )  -> Fill(  myboost.GetTheta()*180/TMath::Pi()   , wgt );
    hmtau_zQ.Get<TH1D>("SqrtQ2", ";Q [GeV]; events", 50, 10, 100) -> Fill ( TMath::Sqrt(Q2), wgt);
    hmtau_zQ.Get<TH1D>("SqrtQ2_lxy", ";Q [GeV]; events", 50, 10, 200) -> Fill ( TMath::Sqrt(Q2), wgt);
    hmtau_zQ.Get<TH1D>("current_pz", ";current hem. p_{z}[GeV];events", 50, -1., 80.) -> Fill( sum_vector.Pz(), wgt);
    hmtau_zQ.Get<TH1D>("current_pz_lxy", ";current hem. p_{z}[GeV];events", 60, 0.9, 80.) -> Fill( sum_vector.Pz(), wgt);
    hmtau_zQ.Get<TH1D>("tau_zQ", ";#tau_{zQ};events",     50,  -.1, 1.) -> Fill(tau, wgt);
    hmtau_zQ.Get<TH1D>("tau_zQ_lxy", ";#tau_{zQ};events", 70, 0.001, 1.2) -> Fill(tau, wgt);
    hmtau_zQ.Get<TH1D>("current_sum_px", 50, -15, 15, "Breit frame: P_{x}, #eta<0", "Entries")  -> Fill(  sum_vector.Px()   , wgt );
    hmtau_zQ.Get<TH1D>("current_sum_py", 50, -15, 15, "Breit frame: P_{y}, #eta<0", "Entries")  -> Fill(  sum_vector.Py()   , wgt );
    hmtau_zQ.Get<TH1D>("current_sum_pt", 50,  -1, 30, "Breit frame: P_{T}, #eta<0", "Entries")  -> Fill(sum_vector.Pt(), wgt);
    hmtau_zQ.Get<TH1D>("current_sum_pt_lxy", 50,  0.1, 40, "Breit frame: P_{T}, #eta<0", "Entries")  -> Fill(sum_vector.Pt(), wgt);
  }
  hmtau_zQ.Get<TH1D>("physicalboost", ";Is physical boost;events", 2, -0.5, 1.5) -> Fill(double( myboost.IsPhysicalBoost() ), wgt);
  hmtau_zQ.Get<TH1D>("physicalboost_ly", ";Is physical boost;events", 2, -0.5, 1.5) -> Fill(double( myboost.IsPhysicalBoost() ), wgt);
  return boostedpart_current; 
}
vector<TLorentzVector> AnalysisEventShapes::BoostParticleArray (const string& hm, const vector<H1PartCand*> partarray, H1Boost& myboost, double Q2, double X, double Y){
  //Takes particle list and boost
  //Returns TLorentzVectors of boosted particles
  H2020HistManager& hmtau_zQ = HistMaster::Instance()->GetHistManager(hm);
  const double wgt = fGen.wgt;
  vector<TLorentzVector> boostedpart_current; //do I need to reset this vector?
  TLorentzVector sum_vector(0., 0., 0., 0.);
  if ( myboost.IsPhysicalBoost() ) {
    for (long unsigned int ipart=0; ipart<partarray.size(); ++ipart){
      //H1Part* part = static_cast<H1Part*>(partarray[ipart]);
      H1PartCand* part = static_cast<H1PartCand*>(partarray[ipart]);
      if (part->IsScatElec()) continue;  // exclude the scattered electron
      if ( !(part->GetE()>0.) ) continue; // reject zeroed-out part cands (badly measured ones, iron muons...)
      TLorentzVector breitpart = myboost.Boost(part->GetFourVector());
      boostedpart_current.push_back(TLorentzVector(breitpart.Px(),breitpart.Py(), breitpart.Pz(), breitpart.E()) );
      if (breitpart.Eta()<0) {
	sum_vector -= breitpart;
      }
    }
    double_t tau = 1-2*sum_vector.Pz()/TMath::Sqrt(Q2);
    if ( tau < 0 ) {
      hmtau_zQ.Get<TH1D>("taucurrent_pz", ";current hem. p_{z}[GeV], #tau<0;events", 51, -1., 80.) -> Fill(sum_vector.Pz(), wgt);
      hmtau_zQ.Get<TH1D>("taucurrent_pz_lxy", ";current hem. p_{z}[GeV], #tau<0;events", 60, 0.9, 80.) -> Fill(sum_vector.Pz(), wgt);
      hmtau_zQ.Get<TH1D>("tauSqrtQ2", ";Q [GeV], #tau<0; events", 50, 10, 100) -> Fill ( TMath::Sqrt(Q2), wgt);
      hmtau_zQ.Get<TH1D>("tauSqrtQ2_lxy", ";Q [GeV], #tau<0; events", 50, 10, 200) -> Fill ( TMath::Sqrt(Q2), wgt);
      hmtau_zQ.Get<TH1D>("tauX", ";X, #tau<0; events", 50, 0, 1) -> Fill ( X, wgt);
      hmtau_zQ.Get<TH1D>("tauX_lxy", ";X, #tau<0; events", 50, 0.001, 1) -> Fill ( X , wgt);
      hmtau_zQ.Get<TH1D>("tauY", ";Y, #tau<0; events", 50, 0, 1) -> Fill ( Y, wgt);
      hmtau_zQ.Get<TH1D>("tauY_lxy", ";Y, #tau<0; events", 50, 0.095, 1) -> Fill ( Y , wgt);
    }

    auto BoostBetaArray = H2020HistManager::MakeLogBinning(50, 0.08, 3.0);
    hmtau_zQ.Get<TH1D>("Boost_Beta_lxy", 50, BoostBetaArray, "Boost: #beta", "events" )  -> Fill(  myboost.GetBeta()   , wgt );
    hmtau_zQ.Get<TH1D>("Boost_Beta_lx" , 50, BoostBetaArray, "Boost: #beta", "events" )  -> Fill(  myboost.GetBeta()   , wgt );
    hmtau_zQ.Get<TH1D>("Boost_Phi", ";Boost: #phi [^{o}]; events", 60, -180, 180)  -> Fill(  myboost.GetPhi()*180/TMath::Pi()   , wgt );
    hmtau_zQ.Get<TH1D>("Boost_Theta", ";Boost: #theta [^{o}]; events", 60,    0, 180 )  -> Fill(  myboost.GetTheta()*180/TMath::Pi()   , wgt );
    hmtau_zQ.Get<TH1D>("SqrtQ2", ";Q [GeV]; events", 50, 10, 100) -> Fill ( TMath::Sqrt(Q2), wgt);
    hmtau_zQ.Get<TH1D>("SqrtQ2_lxy", ";Q [GeV]; events", 50, 10, 200) -> Fill ( TMath::Sqrt(Q2), wgt);
    hmtau_zQ.Get<TH1D>("current_pz", ";current hem. p_{z}[GeV];events", 50, -1., 80.) -> Fill(sum_vector.Pz(), wgt);
    hmtau_zQ.Get<TH1D>("current_pz_lxy", ";current hem. p_{z}[GeV];events", 60, 0.9, 80.) -> Fill(sum_vector.Pz(), wgt);
    hmtau_zQ.Get<TH1D>("tau_zQ", ";#tau_{zQ};events",     70,  -3., 1.) -> Fill(tau, wgt);
    hmtau_zQ.Get<TH1D>("tau_zQ_lxy", ";#tau_{zQ};events", 70, 0.001, 1.2) -> Fill(tau, wgt);
    hmtau_zQ.Get<TH1D>("BreitFramePx", 50, -15, 15, "Breit frame: P_{x}, #eta<0", "Entries")  -> Fill(  sum_vector.Px()   , wgt );
    hmtau_zQ.Get<TH1D>("BreitFramePy", 50, -15, 15, "Breit frame: P_{y}, #eta<0", "Entries")  -> Fill(  sum_vector.Py()   , wgt );
    hmtau_zQ.Get<TH1D>("BreitFramePt", 50,  -1, 30, "Breit frame: P_{T}, #eta<0", "Entries")  -> Fill(sum_vector.Pt(), wgt);
    hmtau_zQ.Get<TH1D>("BreitFramePt_lxy", 50,  0.1, 40, "Breit frame: P_{T}, #eta<0", "Entries")  -> Fill(sum_vector.Pt(), wgt);
    hmtau_zQ.Get<TH1D>("BreitX", ";X; events", 50, 0, 1) -> Fill ( X, wgt);
    hmtau_zQ.Get<TH1D>("BreitX_lxy", ";X; events", 50, 0.001, 1) -> Fill ( X , wgt);
    hmtau_zQ.Get<TH1D>("BreitY", ";Y; events", 50, 0, 1) -> Fill ( Y, wgt);
    hmtau_zQ.Get<TH1D>("BreitY_lxy", ";Y; events", 50, 0.095, 1) -> Fill ( Y , wgt);
  }
  hmtau_zQ.Get<TH1D>("physicalboost", ";Is physical boost;events", 2, -0.5, 1.5) -> Fill(double( myboost.IsPhysicalBoost() ), wgt);
  hmtau_zQ.Get<TH1D>("physicalboost_ly", ";Is physical boost;events", 2, -0.5, 1.5) -> Fill(double( myboost.IsPhysicalBoost() ), wgt);
  return boostedpart_current;
}

void AnalysisEventShapes::PlotKinematicVaribles( const string& hm,double Q2, double X, double Y) {
  
  //Plot kinematic variables Q2, x, y ( linear and log scale )

  H2020HistManager& hmKineVaribles = HistMaster::Instance()->GetHistManager(hm);
  const double wgt = fGen.wgt;
  hmKineVaribles.Get<TH1D>("SqrtQ2", ";Q [GeV]; events", 50, 10, 100) -> Fill ( TMath::Sqrt(Q2), wgt);
  hmKineVaribles.Get<TH1D>("SqrtQ2_lxy", ";Q [GeV]; events", 50, 10, 200) -> Fill ( TMath::Sqrt(Q2), wgt);
  hmKineVaribles.Get<TH1D>("X", "; X; events", 50, 0, 1) -> Fill ( X, wgt);
  hmKineVaribles.Get<TH1D>("X_lxy", "; X; events", 50, 0.001, 1) -> Fill ( X , wgt);
  hmKineVaribles.Get<TH1D>("Y", ";Y; events", 50, 0, 1) -> Fill ( Y, wgt);
  hmKineVaribles.Get<TH1D>("Y_lxy", ";Y; events", 50, 0.095, 1) -> Fill ( Y , wgt);

}


void AnalysisEventShapes::ClassicalEventShapes (const string& hm, const vector<TLorentzVector> BoostedHFS) {
  //Calculate and plot classical event shape observables
  //Pass name of HistManager and vector of boosted particles as argument
  //Boosted particles can be obtained from function BoostParticleArray (Scattered electron already excluded)

  //define hist manager and get weights
  H2020HistManager& hmshapes = HistMaster::Instance()->GetHistManager(hm);
  const double wgt = fGen.wgt;

  //get thrust axis
  TObjArray BoostedHFSThreeVector;
  TVector3 mypart3;
  for (long unsigned int ipart=0; ipart<BoostedHFS.size(); ++ipart){
    TLorentzVector part = BoostedHFS[ipart];
    mypart3 = part.Vect();
    BoostedHFSThreeVector.Add(&mypart3);
  }
  H1EventShape h1es;
  h1es.setPartList(&BoostedHFSThreeVector);
  TVector3 thrustaxis = h1es.thrustAxis();

  //Values for event shape observables will be stored here
  Float_t pz = 0;
  Float_t pt = 0;
  Float_t p = 0;
  Float_t thrust = 0;
  Float_t E = 0;
  TVector3 momentum(0., 0., 0.);
  for (long unsigned int ipart=0; ipart<BoostedHFS.size(); ++ipart){
    TLorentzVector part = BoostedHFS[ipart];
    if ( !(part.E()>0.) ) continue; // reject zeroed-out part cands (badly measured ones, iron muons...)
    if (part.Eta()<0){   //only particles in current hemisphere
      pz += std::fabs( part.Pz() );
      p += std::fabs( part.P());
      pt += std::fabs( part.Pt());
      thrust += std::fabs( part.Vect().Dot(thrustaxis));
      E += part.E();
      momentum += part.Vect();
    }
  }

  //Fill histos
  hmshapes.Get<TH1D>("Thrust_C", ";#tau = 1-T_{C};events", 60, -.1, 1.1) -> Fill(1-thrust/p, wgt);
  hmshapes.Get<TH1D>("LogThrust_C", ";log(#tau) = log(1-T_{C});events", 60, -10, 0.1) -> Fill(TMath::Log(1-thrust/p), wgt);
  hmshapes.Get<TH1D>("Thrust_Z", ";#tau = 1-T_{Z};events", 60, -.1, 1.1) -> Fill(1-pz/p, wgt);
  hmshapes.Get<TH1D>("Broadening", ";Broadening B_{C};events", 60, -.1, 1.1) -> Fill(pt/(2*p), wgt);
  hmshapes.Get<TH1D>("Jet_mass", ";Jet mass #rho;events", 60, -.1, 1.1) -> Fill((TMath::Power(E,2)-momentum.Dot(momentum))/(TMath::Power(2*p,2)), wgt);
  //add jet mass from 1999
  //log plots

}



TLorentzVector AnalysisEventShapes::BoostToBreitFrame ( double Q2, double y, double E0, double Phi ) {
  
  //Calculate virtual photon for iSigma method
  //Independant of initial electron energy -> not sensitive to ISR
  //Initial electron energy given by E = Sigma/2
    
  Double_t px = 0;
  Double_t py = 0;
  Double_t pz = 0;
  Double_t Ee = 0;
  
  //Calculate final state electron
  
  if ((Q2>0) && (y>0) && (y<1)){
    
    Ee = Q2 / (4*E0) + E0*(1-y);

    Double_t b = 4*E0*E0*(1-y)/Q2;
    Double_t Theta = TMath::ACos( (1-b)/(1+b) );
  
    px = Ee*TMath::Sin(Theta)*TMath::Cos(Phi);
    py = Ee*TMath::Sin(Theta)*TMath::Sin(Phi);
    pz = Ee*TMath::Cos(Theta);
    
  }

  TLorentzVector ElecScat(px, py, pz, Ee);
  //Double_t Sigma = gH1Calc->Fs()->GetFullHadFsEmpz();
  Double_t Delta = gH1Calc->Fs()->GetEmpz(); // = Sigma + Ee ( 1-cos(theta) ), equal to full state E - pz
  TLorentzVector ElecInit( 0., 0., -Delta/2, Delta/2); //initial electron, mass neglecte
  return ElecInit - ElecScat;

}

H1Boost AnalysisEventShapes::CalcBoost(double q2, double y, double x, double phi, double Ep) {

  double Epxy = Ep*x*y;  // temporary

  double E0 = q2/(4*Epxy); // electron beam energy
  TLorentzVector elec0(0,0,-E0,E0); // beam electron
  TLorentzVector prot0(0, 0, Ep, Ep); // beam proton

  double ElecE     = Epxy + q2*(1-y)/(4*Epxy);
  double ElecPz    = Epxy - q2*(1-y)/(4*Epxy);
  double Theta     = TMath::ACos(ElecPz/ElecE); // temprorary
  double ElecPx    = ElecE*TMath::Sin(Theta)*TMath::Cos(phi);
  double ElecPy    = ElecE*TMath::Sin(Theta)*TMath::Sin(phi);

  TLorentzVector Elec(ElecPx, ElecPy, ElecPz, ElecE); // scattered electron

  return H1Boost(2*x*prot0, elec0-Elec ,elec0, -prot0 );

} 
