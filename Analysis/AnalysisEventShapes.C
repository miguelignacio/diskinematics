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
   fGen = CrossSectionQuantities();
   fRec = CrossSectionQuantities();
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

   H2020HistManager& hmGen = HistMaster::Instance()->GetHistManager("Gen");
   
   const double wgt = fGen.wgt; //todo
   hmGen.Get<TH1D>("GenQ2e", "Q2e;Gen Q2e [GeV^{2}];events",  60, 100., 30000.)->Fill(gH1Calc->Kine()->GetQ2eGen(), fGen.wgt); 
   hmGen.Get<TH1D>("GenQ2s", "Q2s;Gen Q2s [GeV^{2}];events",  60, 100., 30000.)->Fill(gH1Calc->Kine()->GetQ2sGen(), fGen.wgt); 
   hmGen.Get<TH1D>("Q2Gen",  "Q2Gen;Gen Q2Gen [GeV^{2}];events",  60, 100., 30000.)->Fill(gH1Calc->Kine()->GetQ2Gen(), fGen.wgt ); 

   hmGen.Get<TH1D>("GenQ2e_lx", "Q2e;Gen Q2e [GeV^{2}];events",  60, 100., 30000.)->Fill(gH1Calc->Kine()->GetQ2eGen(), fGen.wgt); 
   hmGen.Get<TH1D>("GenQ2s_lx", "Q2s;Gen Q2s [GeV^{2}];events",  60, 100., 30000.)->Fill(gH1Calc->Kine()->GetQ2sGen(), fGen.wgt); 
   hmGen.Get<TH1D>("Q2Gen_lx",  "Q2Gen;Gen Q2Gen [GeV^{2}];events",  60, 100., 30000.)->Fill(gH1Calc->Kine()->GetQ2Gen(), fGen.wgt ); 

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
   // put some histogram inside
   hmtracks.Get<TH1D>("N tracks BoostedJets",   ";N tracks Boosedjets;entries", 50, 0, 100)  -> Fill(  H1BoostedJets::Instance()->GetHFSArray()->GetEntries()   , wgt ); 
   hmtracks.Get<TH1D>("N tracks BoostedJets_ly",";N tracks Boosedjets;entries", 50, 0, 100)  -> Fill(  H1BoostedJets::Instance()->GetHFSArray()->GetEntries()   , wgt ); 

   // this is how you make a histogram and fill it...
   double tau_zQ = 0.1;
   hm.Get<TH1D>("tau_zQ",    "tau_zQ;#tau_{zQ};events",  70, -0.3, 1.1)  -> Fill(tau_zQ, wgt); 


   // --------------- general kinematics -------------------
   TLorentzVector ebeam, pbeam;
   H1BoostedJets::Instance()->BeamVectors(&pbeam, &ebeam);
   TLorentzVector ScatElec = gH1Calc->Elec()->GetFirstElectron();
   TLorentzVector HFS = gH1Calc->Fs()->GetFullHadFS();
   TLorentzVector virtualphoton = ebeam - ScatElec;
   TLorentzVector Pplusq = pbeam + virtualphoton;
   Float_t W   = Pplusq.Mag();
   Float_t Q2  = GetKine_Q2();
   Float_t X   = GetKine_X();
   Float_t Y   = GetKine_Y();
   Float_t Ye  = GetKine_Ye();
   Float_t Yh  = GetKine_Yh();
   Float_t Empz = gH1Calc->Fs()->GetEmpz();
   Float_t PtMiss = gH1Calc->Fs()->GetPtMiss();
   Float_t prim_zvtx = gH1Calc->Vertex()->GetHatPrimaryVertexZ();
   Float_t optimal_zvtx = gH1Calc->Vertex()->GetZ(H1CalcVertex::vtOptimalNC);
   Bool_t NuclInt = gH1Calc->Vertex()->GetOptimalNCVertexType() == H1CalcVertex::optimalDTNV;
   Float_t gammah_degree = (gH1Calc->Fs()->GetGamma()) * (180.0/TMath::Pi());


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

   auto BoostBetaArray = H2020HistManager::MakeLogBinning(50, 0.08, 3.0);
   hm.Get<TH1D>("Boost_Beta_lxy", 50, BoostBetaArray, "Boost: #beta", " Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("Boost_Beta_lx",  50, BoostBetaArray, "Boost: #beta", " Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("Boost_Phi",  50, -3.1415, 3.1415, "Boost: #phi [^{o}]", " Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("Boost_Theta",50, 0, 3.1415, "Boost: #theta [^{o}]", " Entries")  -> Fill(  14.1   , wgt );

   hm.Get<TH1D>("BreitFramePx", 50, -15, 15, "Breit frame: P_{x}^{h}", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("BreitFramePy", 50, -15, 15, "Breit frame: P_{y}^{h}", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("BreitFramePt", 50,   0, 10, "Breit frame: P_{T}^{h}", "Entries")  -> Fill(  14.1   , wgt );


   // ---------------------- jet multiplicity --------------------

   hm.Get<TH1D>("N_Jets", 7, -0.5, 6.5, "Jet Multiplicity, P_{T} > 7 GeV", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("N_Jets_ly", 7, -0.5, 6.5, "Jet Multiplicity, P_{T} > 7 GeV", "Entries")  -> Fill(  14.1   , wgt );

   hm.Get<TH1D>("N_Jets_ext5", 7, -0.5, 6.5, "Jet Multiplicity, P_{T} > 5 GeV", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("N_Jets_ext5_ly", 7, -0.5, 6.5, "Jet Multiplicity, P_{T} > 5 GeV", "Entries")  -> Fill(  14.1   , wgt );

   hm.Get<TH1D>("N_Jets_ext3", 7, -0.5, 6.5, "Jet Multiplicity, P_{T} > 3 GeV", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("N_Jets_ext3_ly", 7, -0.5, 6.5, "Jet Multiplicity, P_{T} > 3 GeV", "Entries")  -> Fill(  14.1   , wgt );

   // ---------------------- multiplicities --------------------

   hm.Get<TH1D>("N_tracks"  , 51,  -0.5,  50.5, "Number of tracks"  , "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("N_em_parts", 10,  -0.5,   9.5, "Number of em-parts", "Entries")  -> Fill(  14.1   , wgt );

   // angle between the electron candidates
   hm.Get<TH1D>("Els_Delta_Theta", 50, -90, 90,    "#theta^{1}_{e} - #theta^{2}_{e} [^{o}]", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("Els_Delta_Phi",   50, 0, 180,    "#phi^{1}_{e} - #phi^{2}_{e} [^{o}]","Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("Els_Delta",       50, 0, 5,      "#Delta(El_{1},El_{2})","Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("El2_Energy",      50, 0, 50,     "2nd Electron Energy / GeV", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("Acoplanarity_ly", 50, 0.5, 1,  "Cos(|#pi - #Delta #phi|)", "Entries")  -> Fill(  14.1   , wgt );

   hm.Get<TH1D>("N_fwd_tracks_ly", 16, -0.5, 15.5, "N forward tracks",  "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("N_cmb_tracks_ly", 11, -0.5, 10.5, "N combined tracks", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("N_ctr_tracks",    51, -0.5, 50.5, "N central tracks",  "Entries")  -> Fill(  14.1   , wgt );

   // electron-track and HFS object relations
   hm.Get<TH1D>("El_HFS_Delta",         20, 0,  0.5,      "#Delta(El_{1},PartCand)","Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("El_Tracks_Delta",      20, 0,  0.5,      "#Delta(El_{1},Track)","Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("El_Clusters_Delta",    20, 0,  0.5,      "#Delta(El_{1},Cluster)","Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("Mult_Ele_Tracks",       2, -0.5, 1.5,    "Multiple electron tracks","Entries")  -> Fill(  14.1   , wgt );


   // ---------------------- track quantities --------------------

   hm.Get<TH1D>("theta_fwd_tracks", 40, 0,  40, "#theta fwd tracks [^{o}]", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("phi_fwd_tracks",   50, 0, 180, "#phi fwd tracks [^{o}]", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("p_fwd_tracks_ly",  40, 0,  40, "momentum fwd tracks [GeV]", "Entries")  -> Fill(  14.1   , wgt );

   hm.Get<TH1D>("theta_cmb_tracks", 40, 0,  40, "#theta cmb tracks [^{o}]", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("phi_cmb_tracks",   50, 0, 180, "#phi cmb tracks [^{o}]", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("p_cmb_tracks_ly",  40, 0,  40, "momentum cmb tracks [GeV]", "Entries")  -> Fill(  14.1   , wgt );

   hm.Get<TH1D>("theta_ctr_tracks", 50, 0, 180, "#theta ctr tracks [^{o}]", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("phi_ctr_tracks",   50, 0, 180, "#phi ctr tracks [^{o}]", "Entries")  -> Fill(  14.1   , wgt );
   hm.Get<TH1D>("p_ctr_tracks_ly",  50, 0,  50, "momentum ctr tracks [GeV]", "Entries")  -> Fill(  14.1   , wgt );

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

