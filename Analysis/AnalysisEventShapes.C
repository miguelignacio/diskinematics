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

   fTreeVar.gen_event_Q2 =  gH1Calc->Kine()->GetQ2esGen();
   fTreeVar.gen_event_x  =  gH1Calc->Kine()->GetXesGen();
   fTreeVar.gen_event_y  =  gH1Calc->Kine()->GetYesGen();
   fTreeVar.gene_px = ScatElecGen.Px();
   fTreeVar.gene_py = ScatElecGen.Py();
   fTreeVar.gene_pz = ScatElecGen.Pz();
   fGen.Q2 = gH1Calc->Kine()->GetQ2esGen();

   //const double etamin = -2.5;
   const double etamax = 2.75;

   vector<PseudoJet> hadron_event_lab;

   TObjArray* hadrons_lab = H1BoostedJets::Instance()->GetHadronArray();
   for (int ipart=0; ipart<hadrons_lab->GetEntries(); ++ipart){
     H1PartMC* part = (H1PartMC*)hadrons_lab->At(ipart);
     if ( part->GetEta() > etamax ) continue;
     fastjet::PseudoJet particle(part->GetPx(),part->GetPy(), part->GetPz(), part->GetE());
     particle.set_user_index(part->GetCharge());
     hadron_event_lab.push_back(particle);
   }

   fastjet::JetDefinition jet_def(fastjet::genkt_algorithm, 1.0,1);
   // --- jets in lab-frame particles (reco level)                                                                                                                                                        
  
    
    ClusterSequence clust_seq_hadlab(hadron_event_lab, jet_def);
   vector<PseudoJet> hadronJetsLab = sorted_by_pt(clust_seq_hadlab.inclusive_jets(3.0));



   fGen.genjets = hadronJetsLab;


   std::vector<float> v_gen_jetpt;                                                                                                                                                                      
   std::vector<float> v_gen_jeteta;                                                                                                                                                                      
   std::vector<float> v_gen_jetphi;                                                                                                                                                                      
   std::vector<float> v_gen_jetz;

   std::vector<float> v_jetpt;
   std::vector<float> v_jeteta;
   std::vector<float> v_jetphi;
   std::vector<float> v_jetz;
   TLorentzVector jet_vec;

 for (unsigned ijet= 0; ijet <  fGen.genjets.size();ijet++) {
     fastjet::PseudoJet genjet =  fGen.genjets[ijet];
     jet_vec.SetPtEtaPhiE(genjet.perp(), genjet.eta(), genjet.phi(), genjet.e());
     float deltaR = 999;
     int matched_index = -999;
     v_gen_jetpt.push_back(genjet.perp());
     v_gen_jeteta.push_back(genjet.eta());
     v_gen_jetphi.push_back(genjet.phi());
     double z = pbeam.Dot(jet_vec)/pbeam.Dot(virtual_photon); 
     v_gen_jetz.push_back(z);
     v_jetpt.push_back(-9999);
     v_jeteta.push_back(-9999);
     v_jetphi.push_back(-9999);
     v_jetz.push_back(-9999);
   }


  fTreeVar.gen_jet_pt = v_gen_jetpt;
  fTreeVar.gen_jet_eta = v_gen_jeteta;
  fTreeVar.gen_jet_phi = v_gen_jetphi;
  fTreeVar.gen_jet_z   = v_gen_jetz;

  
  v_gen_jetpt.clear();
  v_gen_jeteta.clear();
  v_gen_jetphi.clear();
  v_gen_jetz.clear();
  
  fTreeVar.jet_pt = v_jetpt;
  fTreeVar.jet_eta = v_jeteta;
  fTreeVar.jet_phi = v_jetphi;
  fTreeVar.jet_z  = v_jetz;
  v_jetpt.clear();
  v_jeteta.clear();
  v_jetphi.clear();
  v_jetz.clear();

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
   // basic DIS observables, iSigma method
   fRec.Q2       = gH1Calc->Kine()->GetQ2es();
   fRec.Y        = gH1Calc->Kine()->GetYes();
   fRec.X        = ( ScatElec.E() / pbeam.E() ) * ( TMath::Power( TMath::Cos(ScatElec.Theta()/2) , 2 ) / fRec.Y );

   // boost to Breit Frame
   fRec.BoostToBreit = ESTools.BoostToBreitFrame(fRec.Q2, fRec.Y, fRec.X, ScatElec.Phi());
   fRec.BoostToLab = ESTools.BoostToLabFrame(fRec.BoostToBreit);

   // event shapes
   fRec.tau_zQ   =-9;
   fRec.tau_zP   = -9;
   //minitree

   Float_t hadptda = gH1Calc->Fs()->GetHadPtDa();
   TLorentzVector HFS = gH1Calc->Fs()->GetAllFsLessElectron();


   fTreeVar.event_Q2 = gH1Calc->Kine()->GetQ2es();
   fTreeVar.event_y  = gH1Calc->Kine()->GetYes();
   fTreeVar.event_x  = gH1Calc->Kine()->GetXes();
   fTreeVar.Empz     = gH1Calc->Fs()->GetEmpz();
   fTreeVar.e_px = ScatElec.Px();
   fTreeVar.e_py = ScatElec.Py();
   fTreeVar.e_pz = ScatElec.Pz();
   fTreeVar.ptmiss =  gH1Calc->Fs()->GetPtMiss();
   fTreeVar.pth    =  gH1Calc->Fs()->GetPtCalo();
   fTreeVar.vertex_z = gH1Calc->Vertex()->GetZ();
   fTreeVar.ptratio_da = HFS.Pt()/hadptda;
   fTreeVar.ptratio_ele = HFS.Pt()/ScatElec.Pt();

   // fTreeVar.pth    = gH1Calc->Fs()->GetPt();
   fTreeVar.tau1b = -9;//fRec.tau1b;
   //fTreeVar.gen_tau1b = fGen.tau1b;
   fTreeVar.tauzQ = -9;//fRec.tau_zQ;
   //fTreeVar.gen_tauzQ = fGen.tau_zQ;

   //const double etamin = -1.5;
   const double etamax = 2.75;
   // -- define fastjet vectors                                                                                                                                                                           
  
   vector<PseudoJet> full_event_lab;
   //loop over HFS array
   for (long unsigned int ipart=0; ipart<particlearray.size(); ++ipart){
     H1PartCand* part = static_cast<H1PartCand*>(particlearray[ipart]);
     if (part->IsScatElec()) continue;  // exclude the scattered electron                                                                                                                                 
     if ( part->GetEta() > etamax ) continue; // cut on eta, both for lab and breit frame                                                                                      
     fastjet::PseudoJet reco_particle(part->GetPx(),part->GetPy(), part->GetPz(), part->GetE()); 
     reco_particle.set_user_index(part->GetCharge());
     full_event_lab.push_back( reco_particle);
   }
   // --- define fastjet contrib module Centauro                                                                                                                                                          
  
   //fastjet::contrib::CentauroPlugin * centauro_plugin = new fastjet::contrib::CentauroPlugin(1.0);                                                                                                       
   //fastjet::JetDefinition jet_def(centauro_plugin);                                                                                                                                                     
  
   fastjet::JetDefinition jet_def(fastjet::genkt_algorithm, 1.0,1);

   // --- jets in lab-frame particles (reco level)                                                                                                                                                        
  
   ClusterSequence clust_seq_lab(full_event_lab, jet_def);
   vector<PseudoJet> jets_lab = clust_seq_lab.inclusive_jets(5.0);
   vector<PseudoJet> sortedLabJets   = sorted_by_pt(jets_lab);
   

   std::vector<float> v_jetpt; 
   
   std::vector<float> v_jeteta;
   
   std::vector<float> v_jetphi;
   
   std::vector<float> v_jetz;
   // std::vector<float> v_jetcharge;
   //std::vector<float> v_gen_jetcharge;

   TLorentzVector jet_vec;
 
   if( not(fGen.Q2>0)){ //if real data, just fill the reco jets
     //  std::cout<< " REAL DATA " <<std::endl;

       for (unsigned ijet= 0; ijet < sortedLabJets.size();ijet++) {
           fastjet::PseudoJet jet = sortedLabJets[ijet];
           v_jetpt.push_back(jet.perp());
	   v_jeteta.push_back(jet.eta());
	   v_jetphi.push_back(jet.phi());
	   jet_vec.SetPtEtaPhiE(jet.perp(), jet.eta(), jet.phi(), jet.e());
           double z = pbeam.Dot(jet_vec)/pbeam.Dot(virtual_photon);
	   v_jetz.push_back(z);

	   vector<PseudoJet> constituents = jet.constituents();
	   int ncharged = 0;
	   float jet_charge = 0;
	   for (unsigned j = 0; j < constituents.size(); j++) {
	       if(!(constituents[j].user_index()!=0)) continue;
	       ncharged = ncharged +1;
	       TVector3 jetvector(jet.px(), jet.py(), jet.pz());
	       TVector3 hvector(constituents[j].px(), constituents[j].py(), constituents[j].pz());
	       jet_charge += constituents[j].user_index() *pow( hvector.Perp()/jet.perp(), 0.3); 
	   }// end loop over jet constituents
	   //v_jetcharge.push_back(jet_charge);
       }
   }

   
   //MC
   // Loop over generated jets
   for (unsigned ijet= 0; ijet <  fGen.genjets.size();ijet++) {       
       fastjet::PseudoJet genjet =  fGen.genjets[ijet];                                                                                                                                                    
       float deltaR = 999;                                                                                                                                                                                
       int matched_index = -999;                                                                                                                                                                           
      //matching generated jets with reconstructed jets                                                                                                                                                    
       for (unsigned kjet= 0; kjet < sortedLabJets.size();kjet++) {                                                                                                                                      
          fastjet::PseudoJet jet = sortedLabJets[kjet];                                                                                                                                                    
          if(genjet.delta_R(jet) < deltaR)                                                                                                                                                                 
          {                                                                                                                                                                                             
              deltaR = genjet.delta_R(jet);                                                                                                                                                                
              matched_index = kjet;                                                                                                                                                                        
          }                                                                                                                                                                                              
       }//end loop over reconstructed jets                                                                                                                                                               


       if( not(matched_index>-999 and deltaR<0.9)){
           v_jetpt.push_back(-9999);
           v_jeteta.push_back(-9999);
           v_jetphi.push_back(-9999);
           v_jetz.push_back(-9999);
	   //      v_jetcharge.push_back(-9999);
       }
       else{
	 fastjet::PseudoJet matchedjet = sortedLabJets[matched_index];                                                                                                                                    
	 v_jetpt.push_back(matchedjet.perp());
	 v_jeteta.push_back(matchedjet.eta());
	 v_jetphi.push_back(matchedjet.phi());
         TLorentzVector jet_vec;
	 jet_vec.SetPtEtaPhiE(matchedjet.perp(), matchedjet.eta(), matchedjet.phi(), matchedjet.e());
	 double z = pbeam.Dot(jet_vec)/pbeam.Dot(virtual_photon);
	 v_jetz.push_back(z);

       }
   }// loop over gen jets


       /*
       double gen_jetcharge = 0 ;
       //Loop over constituents
       vector<PseudoJet> genconstituents = genjet.constituents();
       for (unsigned j = 0; j < genconstituents.size(); j++) {
           if( not(genconstituents[j].user_index()!=0)) continue;
               TVector3 jetvector(genjet.px(), genjet.py(), genjet.pz());
               TVector3 hvector(genconstituents[j].px(), genconstituents[j].py(), genconstituents[j].pz());
               gen_jetcharge += genconstituents[j].user_index() *pow( hvector.Perp()/genjet.perp(), 0.3); 
       } //end loop over generator-level constituents
       v_gen_jetcharge.push_back(gen_jetcharge);
   }//loop over gen jets
   */
            
   fTreeVar.jet_pt = v_jetpt;
   fTreeVar.jet_eta = v_jeteta;
   fTreeVar.jet_phi = v_jetphi;
   fTreeVar.jet_z  = v_jetz;
   v_jetpt.clear();
   v_jeteta.clear();
   v_jetphi.clear();
   v_jetz.clear();
      
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
