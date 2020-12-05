#include "EventshapeTools.h"
#include "H1PhysUtils/H1BoostedJets.h"

ClassImp(EventshapeTools)


//Float_t JetTools::fJetMatchingRadius = 0.9;

EventshapeTools::EventshapeTools()
{
  // default constructor
}

EventshapeTools::~EventshapeTools()
{
  // default destructor
}

//______________________________________________________

TLorentzVector EventshapeTools::CalcScatElec(double q2, double y, double x, double phi, double Ep)
{
  // Calculate the scattered electron four-vector for the boost 
  // to the Breit frame. Reconstruct the vector from Q2, y, x, Phi_e and E_p,
  // Different reconstruction methods possible (ISigma is default on gen level) 

   double ElecE     = 0;
   double ElecPz    = 0;
   double ElecPy    = 0;
   double ElecPx    = 0;
   if ((q2>0) && (y>0) && (y<1)){
      double Epxy = Ep*x*y;  // temporary
      ElecE     = Epxy + q2*(1-y)/(4*Epxy);
      ElecPz    = Epxy - q2*(1-y)/(4*Epxy);
      double Theta     = TMath::ACos(ElecPz/ElecE); // temporary
      ElecPx    = ElecE*TMath::Sin(Theta)*TMath::Cos(phi);
      ElecPy    = ElecE*TMath::Sin(Theta)*TMath::Sin(phi);
   }
   TLorentzVector Elec(ElecPx, ElecPy, ElecPz, ElecE); // scattered electron

   return Elec;

}

//______________________________________________________


H1Boost EventshapeTools::BoostToBreitFrame(double q2, double y, double x, double phi){
 
   TLorentzVector ebeam, pbeam;
   H1BoostedJets::Instance()->BeamVectors(&pbeam, &ebeam);
   
   // Calculate beam electron
   double Epxy = pbeam.E()*x*y;  // temporary
   double E0 = q2/(4*Epxy); // electron beam energy
   TLorentzVector elec0(0,0,-E0,E0); // beam electron
   
   // Get scattered electron for boost
   TLorentzVector Elec = CalcScatElec(q2, y, x, phi, pbeam.E());
   
   // Boost to reit frame
   H1Boost Boost_To_Breit(2*x*pbeam, elec0-Elec , elec0, -pbeam );
   return Boost_To_Breit;
}

//______________________________________________________


H1Boost EventshapeTools::BoostToLabFrame(H1Boost boost){
   
   // Supply H1Bost as argument
   // Gives the H1Boot back to the initial frame

   // set up Lorentz-vectors of Beam Particles and the scattered electron
   TLorentzVector RestProton(0.,0.,0.,TDatabasePDG::Instance()->GetParticle(2212)->Mass()); // proton mass
   TLorentzVector ZeroLorentzVector(0.,0.,0.,0.);

   TLorentzVector ebeam, pbeam;
   H1BoostedJets::Instance()->BeamVectors(&pbeam, &ebeam);

   // prepare the boost back to the lab-frame:
   TLorentzVector RestProton_bf = boost.Boost(RestProton); // should be at rest again after the boost back
   TLorentzVector xz_lab(1.0,0.0,1.0,TMath::Sqrt(2.0)); // some LorentzVector defining the x-z Plane 
   TLorentzVector xz_bf = boost.Boost(xz_lab); 
   TLorentzVector BeamProt_bf = boost.Boost(pbeam); // defines the positive z-direction
   H1Boost BoostBackToLab(RestProton_bf,ZeroLorentzVector,xz_bf,-BeamProt_bf);

   return BoostBackToLab;

}
