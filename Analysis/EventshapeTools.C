#include "EventshapeTools.h"
#include "H1PhysUtils/H1BoostedJets.h"
#include <H1Calculator/H1CalcEvent.h>
#include <H1Geom/H1Constants.h>
#include <H1Calculator/H1Calculator.h>
#include <TF1.h>
#include <H1Calculator/H1CalcElec.h>
#include <H1Calculator/H1CalcWeight.h>

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

void EventshapeTools::ApplyNCTrackClusterWeight()
{

  // apply the Track-Cluster Weight as determined for a cut on DCA(track,cluster) < 8 cm
  // see HaQ meeting in June/July 2010
  // Taken from JetReweighter.C from the jets at high Q2 analysis by Roman Kogler
  // Calculates new event weight
  // Data not well modeled in MC -> Electron angle theta
  // see PhD Thesis Roman Kogler, section 9.4

   
  H1CalcEvent* Event = H1Calculator::Instance()->Event();
  if (!Event->IsMC()) return;
      
  static H1IntPtr RunNumber("RunNumber");
  static Int_t OldRunNumber = 0;
  static H1Constants* con = H1Constants::Instance();
      
  Float_t TrkClsWeight = 1.0;
  TLorentzVector ScatElec = gH1Calc->Elec()->GetFirstElectron();
  Float_t e_theta_degree = (ScatElec.Theta())*(180.0/TMath::Pi());
  static TF1* MCVtxTrackThetaEff = new TF1("MCVtxTrackThetaEff", "[0]+[1]*x+[2]*x*x+[3]*pow(x,3)", 30, 160);
   
  if(OldRunNumber != *RunNumber) {
    OldRunNumber = *RunNumber;
        
    if     ( con->GetPolPeriod(*RunNumber) >= H1Constants::eEplus06RH1       ){ // 06-07 e+
      MCVtxTrackThetaEff->SetParameters(1.0356, -0.00189382, 2.43179e-05, -9.55152e-08);
    } else if( con->GetPolPeriod(*RunNumber) >= H1Constants::eEminus06LH1    ){ // 06 e-
      MCVtxTrackThetaEff->SetParameters(1.0443, -0.00160364, 1.57696e-05, -5.27797e-08);
    } else if( con->GetPolPeriod(*RunNumber) >= H1Constants::eEminus04RH1    ){ // 04-05 e-
      MCVtxTrackThetaEff->SetParameters(1.02605, -0.00147769, 1.60599e-05, -5.51945e-08);
    } else if( con->GetPolPeriod(*RunNumber) >= H1Constants::eEplus03RH1     ){ // 03-04 e+
      MCVtxTrackThetaEff->SetParameters(1.16988, -0.00611282, 6.26369e-05, -2.03377e-07);
    } else if( con->GetPolPeriod(*RunNumber) >= H1Constants::eEplus9900UNPOL ){ // 99-00 e+
      MCVtxTrackThetaEff->SetParameters(1.0, 0.0, 0.0, 0.0);
    } else if( con->GetPolPeriod(*RunNumber) >= H1Constants::eEminus9899UNPOL){ // 98-99 e-
      MCVtxTrackThetaEff->SetParameters(1.0, 0.0, 0.0, 0.0);
    }

  }
   
  if (e_theta_degree>30){
    TrkClsWeight *= MCVtxTrackThetaEff->Eval(e_theta_degree);
  }

  static Int_t FirstCall = 1;
  if(FirstCall == 1) {
    
    TString PeriodName    = "Unknown Period";
    if     ( con->GetPolPeriod(*RunNumber) >= H1Constants::eEplus06RH1      ) PeriodName = "e+p 06/07";  // 06-07 e+
    else if( con->GetPolPeriod(*RunNumber) >= H1Constants::eEminus06LH1     ) PeriodName = "e-p 06";  // 06 e-
    else if( con->GetPolPeriod(*RunNumber) >= H1Constants::eEminus04RH1     ) PeriodName = "e-p 04/05";  // 04-05 e-
    else if( con->GetPolPeriod(*RunNumber) >= H1Constants::eEplus03RH1      ) PeriodName = "e+p 03/04";  // 03-04 e+
    else if( con->GetPolPeriod(*RunNumber) >= H1Constants::eEplus9900UNPOL  ) PeriodName = "e+p 99/00";  // 99-00 e+
    else if( con->GetPolPeriod(*RunNumber) >= H1Constants::eEminus9899UNPOL ) PeriodName = "e-p 98/99";  // 98-99 e-

    cout << "------------------------------------------------------------" << endl;
    cout << "Apply Vertex-Track-Cluster link weight for " << PeriodName << endl;
    cout << "TrkClsWeight: Parameters of 3rd order polynomial: " << endl;
      for (Int_t i=0;i<4;++i) cout << "              p[" << i << "] = " << MCVtxTrackThetaEff->GetParameter(i) << endl;
      cout << "------------------------------------------------------------" << endl;
      FirstCall = 0;
   }
  
   // apply track-cluster weight
   Double_t OrigWeight = gH1Calc->Weight()->GetWeight();
   //fTrackClusterWeight += TrkClsWeight;
   //fRecWeight *= TrkClsWeight;
   gH1Calc->Weight()->SetWeight(TrkClsWeight*OrigWeight);
   
}
