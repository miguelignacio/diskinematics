#include "JetTools.h"
#include "H1Calculator/H1Calculator.h"
#include "H1Calculator/H1CalcJet.h"
#include "H1Mods/H1PartGenJetArrayPtr.h"
#include "H1PhysUtils/H1BoostedJets.h"
#include "H1Pointers/H1FloatPtr.h"
#include <TRandom.h>


ClassImp(JetTools)


Float_t JetTools::fJetMatchingRadius = 0.9;

JetTools::JetTools()
{
  // default constructor
}

JetTools::~JetTools()
{
  // default destructor
}

//______________________________________________________


H1ArrayI JetTools::MatchModsJets(Float_t R0)
{
  // matches jets between hadron and reconstructed level 
  // with geometrical distance (eta-phi space)
  // the array returned means that array[i] = j 
  // means that hadron jet i matches with reconstructed jet j
  // j = -1 if no matching reconstructed jet was found

  static H1PartGenJetArrayPtr genjets;
  Int_t NumRecJets = gH1Calc->Jet()->GetNumJets();
  
  H1ArrayI MatchingIdx(20);

  for (int i=0; i<20; ++i) MatchingIdx[i] = -1; // initialize the array 

  for (int ihad=0; ihad<genjets.GetEntries(); ++ihad){
    H1PartJet* hadjet = genjets[ihad];

    // test if the generated jet fulfills minimum requirement (Pt>7)
    if (hadjet->GetPt() < 7.) continue;

    TLorentzVector hadjetvec = hadjet->GetFourVector();
    Float_t dist = 0;
    Float_t min_dist = FLT_MAX;
    
    for (int irec=0; irec<NumRecJets; ++irec){
      TLorentzVector recjetvec = gH1Calc->Jet()->GetJet(irec);

      // reconstructed jet should have Pt > 5 GeV
      //if (recjet->GetPt() < 5.) continue;  Don't cut on the rec jet!

      dist = Dist(&hadjetvec, &recjetvec);
      
      if (dist<min_dist && dist<R0){
	min_dist = dist;
	MatchingIdx[ihad] = irec;
      }

    }

  }
  return MatchingIdx;
}


//______________________________________________________


H1ArrayD JetTools::UnmatchedRecJetsPt( H1ArrayI MatchedGenJets )
{

  cout << "ERROR in JetTools: The function JetTools::UnmatchedRecJetsPt is obsolete." << endl;
  return 0;

//   // Method for determining the pt of all 'unmatched' jets in an event.
//   // retruns an H1ArrayD with the pt of the unmatched jets.
//   // use: MatchedGenJets = JetTools::MatchJetsUnfold(...) 

//   H1ArrayD* jetPt = new H1ArrayD();
  
//   // RecJets
//   H1BoostedJets* boostedjetfinder = H1BoostedJets::Instance();
//   TObjArray* boostedjetsRec = boostedjetfinder->GetBoostedJets(H1BoostedJets::eRecLev);
//   //TObjArray* boostedjetsGen = boostedjetfinder->GetBoostedJets(H1BoostedJets::eHadLev);
//   Int_t NumRecJets = boostedjetsRec->GetEntries();

//   for ( Int_t irec = 0 ; irec < NumRecJets ; irec++ ) {

//     Bool_t IsMatched = false;

//     for ( Int_t igen = 0 ; igen < MatchedGenJets.GetEntries() ; igen++ ) {
//       if ( IsMatched && irec == MatchedGenJets[igen] ) cout << "WARNING: One Rec Jet is matched by two Gen Jets. NumRecJets: "  <<  NumRecJets << "\tNGenJets: " << MatchedGenJets.GetEntries()  << endl;
//       if ( irec == MatchedGenJets[igen] ) IsMatched = true;
//     }

//     if ( !IsMatched ) 
//       jetPt->Add( ((H1PartJet*)boostedjetsRec->At(irec))->GetPt() );  
//   }
  
//   return *jetPt; // return array mit den pts jets derjenigen jets die nicht gematched wurden.

}


//______________________________________________________

H1ArrayB JetTools::MatchedRecJetsArray( H1ArrayI MatchedGenJets )
{
  // Method for determining if each rec-Jet was previously matched
  // by an GenJet.
  // use: MatchedGenJets = JetTools::MatchJetsUnfold(...) 

  H1ArrayB IsMatched;// = new H1ArrayB();
  
  // RecJets
  Int_t NumRecJets = H1BoostedJets::Instance()->GetBoostedJets(H1BoostedJets::eRecLev)->GetEntries();

  for ( Int_t irec = 0 ; irec < NumRecJets ; irec++ ) {
    
    IsMatched.Add( false );
    Bool_t matched = false;

    for ( Int_t igen = 0 ; igen < MatchedGenJets.GetEntries() ; igen++ ) {
      if ( matched && irec == MatchedGenJets[igen] ) {
	 cout << "WARNING: One Rec Jet is matched by two Gen Jets. NRecJets: " <<  NumRecJets << "\tNGenJets: " << MatchedGenJets.GetEntries()  << endl;
      }
      else if ( irec == MatchedGenJets[igen] ) { // no double counting (therefore 'else if')
	 matched = true;
	 IsMatched.AddAt(  true , irec );  // overwrite this entry
      }
    }
    
  }
  
  return IsMatched; // return array mit den pts jets derjenigen jets die nicht gematched wurden.

}


//______________________________________________________


H1ArrayI JetTools::MatchJetsUnfold( Bool_t* GenJetsPSArr, Bool_t* RecJetsPSArr ,Int_t DontMatchGenI , Int_t DontMatchRecI) 
{

  // Do not use this function any more
  // use the new jetmatching instead

  cout << " Do not use MatchJetsUnfold, but MatchJetsUnfold2() " << endl;

  // Mathing method to match jets in the unfolding procedure
  //
  // matches jets between generated and reconstructed level 
  // in geometrical distance (eta-phi space) and espc. does not 
  // take pt into account.
  //
  // The array returned means that array[i] = j 
  // means that gen-jet i matches best with reconstructed jet j
  //
  // for no matching:
  // j = -1 if no matching reconstructed jet was found
  // 
  // DontMatchGenI and DontMatchRecI
  // If one knows, that one Rec Jet was matched by two gen jets, then specify those
  // jets-indices (on rec and gen level) and leave them out in the matching procedure.
  //

  // method:
  // 1) loop over gen-jets
  // 2) for each gen-jet look for 'best matching jet'
  //    best matching means:
  //    a) closest jet in r=JetTools::Dist(...)
  //    b) but r<fJetMatchingRadius
  //   
  // ==============================================================


  Bool_t ThisIsAnIteration = false;
  if ( DontMatchGenI != -1 || DontMatchRecI != -1 ) ThisIsAnIteration = true;
  //else ThisIsAnIteration = false;

  //!  initialize the array
  H1ArrayI MatchingIdx;// = new H1ArrayI();

  //!  Jets
  H1BoostedJets* boostedjetfinder = H1BoostedJets::Instance();
  TObjArray* boostedjetsRec = boostedjetfinder->GetBoostedJets(H1BoostedJets::eRecLev);
  TObjArray* boostedjetsGen = boostedjetfinder->GetBoostedJets(H1BoostedJets::eHadLev);
  Int_t NumRecJets = boostedjetsRec->GetEntries();

  //   if ( boostedjetsGen->GetEntries() < NumRecJets ) {
  //     cout << "grepme: LESS gen jets [ " << boostedjetsGen->GetEntries() << " ] than rec jets. rec [ " << NumRecJets << " ]. " << endl;
  //   }
  //   if ( boostedjetsGen->GetEntries() > NumRecJets ) {
  //     cout << "grepme: MORE gen jets [ " << boostedjetsGen->GetEntries() << " ] than rec jets. rec [ " << NumRecJets << " ]. " << endl;
  //   }

  //   if ( ThisIsAnIteration ) cout << "grepme: Iteration with: " << DontMatchGenI <<" " << DontMatchRecI <<" : gen jets [ " << boostedjetsGen->GetEntries() << " ]|  rec jets [ " << NumRecJets << " ]"  << endl;
    
  for ( Int_t igen = 0 ; igen < boostedjetsGen->GetEntries()  ; igen++ ) {
  
    //     if ( ThisIsAnIteration ) cout << " grepme\t" << igen << ":   ";

    MatchingIdx.Add(-1);

    if ( !GenJetsPSArr[igen] ) { 
      //cout << "grepme: testing: this gen jet should be out of the phase space: pt = " <<   hadjet->GetPt() << "\teta = " << hadjet->GetEta()  << endl;
      continue;
    }
    
    //! if iteration
    if ( DontMatchGenI == igen  ) {
      MatchingIdx.AddAt(  DontMatchRecI , igen  );
      //       cout << endl;
      continue;
    }    
    
    H1PartJet* hadjet = (H1PartJet*)boostedjetsGen->At(igen);

    TLorentzVector hadjetvec = hadjet->GetFourVector();
    Float_t dist = -1.;
    Float_t min_dist = 1000; //FLT_MAX
    
    for ( Int_t irec = 0 ; irec < NumRecJets ; irec++ ) {
      
      if ( !RecJetsPSArr[irec] ) {
	//cout << "grepme: testing: this rec jet should be out of the phase space:  pt = " <<  PtDet << "\teta = " << EtaDet  << endl;
	continue;
      }

      if ( DontMatchRecI == irec  ) {
	continue;
      }
      
      H1PartJet* jet_rec = (H1PartJet*)boostedjetsRec->At(irec);
      Double_t PtDet     = jet_rec->GetPt();
      Double_t EtaDet    = jet_rec->GetEta();

      //       if ( ThisIsAnIteration ) cout << irec << " ";
      TLorentzVector recjetvec = ((H1PartJet*)boostedjetsRec->At(irec))->GetFourVector();
      dist = Dist(&hadjetvec, &recjetvec);
      
      //! re-uncomment this, if you are interested in double matched jets
      //       if ( (dist >= min_dist && dist < fJetMatchingRadius) || (min_dist != 1000 && dist < fJetMatchingRadius)  ){
      // 	cout << "grepmegrepme: I could also match a gen jet!  PtDet: " << PtDet << "  EtaDet: " << EtaDet << endl;
      //       }

      if ( dist<min_dist && dist<fJetMatchingRadius ){
	min_dist = dist;
	MatchingIdx.AddAt( irec, igen); // AddAt( number, index )
      }
    }
    //     cout << endl;
  }

  Int_t SortOutGen = -1;
  Int_t SortOutRec = -1;
  Bool_t ReRun     = false;

  //// Int_t leftOver   = -42;

  //!  check if one Rec Jet was matched twice.
  for ( Int_t igen = 0 ; igen < boostedjetsGen->GetEntries()  ; igen++ ) {
    for ( Int_t igen2 = igen+1 ; igen2 < boostedjetsGen->GetEntries()  ; igen2++ ) {
      if ( ( (MatchingIdx)[igen2] == (MatchingIdx)[igen] ) && (MatchingIdx)[igen] != -1 ) {
	ReRun = true;
	//cout << "grepme: Two Gen Jets are matching one Rec Jet. ReRun with this exception the JetMatching." << endl;
	if (ThisIsAnIteration) cout << "grepme: ThisIsAnIteration and there is still one rec jet matched by two gen jets." << endl;
	SortOutRec = (MatchingIdx)[igen];

	// take one Gen-Jet out randomly ( Or the closest ??? )
	static TRandom r;
	if ( r.Rndm() >= 0.5){
	  SortOutGen = igen;
	  //// leftOver   = igen2;
	} else {
	  SortOutGen = igen2;
	  //// leftOver   = igen;
	}
      }
    }
  }
  
  //!  ReRunning
  if ( ReRun && !ThisIsAnIteration ) {
    //cout << "SortOutGen: " << SortOutGen << "\tSortOutRec: " << SortOutRec << endl;
    MatchingIdx = MatchJetsUnfold (   GenJetsPSArr, RecJetsPSArr, SortOutGen , SortOutRec );
    // useless cout: if ( SortOutGen >= MatchingIdx->GetEntries() ) cout << "grepme: ERROR in JetTools!\t SortOutGen: " << SortOutGen << "\tentries: " << MatchingIdx->GetEntries() << endl;
    MatchingIdx.AddAt( SortOutRec ,  SortOutGen );

  }

  //!  check if one Rec Jet is still matched twice (mostly for the case if many_gen_jets and few_rec_jets).
  if ( ReRun && !ThisIsAnIteration ) { // save cpu time
    for ( Int_t k = 0 ; k<2 ; k++ ){ // more tries!
      for ( Int_t igen = 0 ; igen < boostedjetsGen->GetEntries()  ; igen++ ) {
	for ( Int_t igen2 = igen+1 ; igen2 < boostedjetsGen->GetEntries()  ; igen2++ ) {
	  if ( ( (MatchingIdx)[igen2] == (MatchingIdx)[igen] ) && (MatchingIdx)[igen] != -1 ) {
	    // take one Gen-Jet out randomly ( Or the closest ??? )
	    //if ( k == 0 ) cout << "grepme: current debugging: igen: " << igen << "  igen2: " << igen2 << endl;
	    //if (k>0) cout << "grepme: WOW! one rec jet was matched by three (more?) gen jets. I remove one randomly." << endl;
	    static TRandom r;
	    if ( r.Rndm() >= 0.5)
	      (MatchingIdx)[igen2] = -1;
	    else
	      (MatchingIdx)[igen] = -1;
	  }
	}
      }
    }
  }

  return MatchingIdx;

}


//_________________________________________________________________________________________________________________


H1ArrayI JetTools::MatchJetsUnfold2 ( Bool_t* GenJetsPSArr, Bool_t* RecJetsPSArr ) 
{

  //
  // new jet matching!
  //
  // - allow extended phase (pt[,Q2,eta]) space on gen level
  // - look for closest jet within jetmatching radius
  //   
  // The array returned means that array[i] = j 
  // means that gen-jet i matches best with reconstructed jet j
  //
  // for no matching:
  // j = -1 if no matching reconstructed jet was found 
  //


  //!  Jets
  H1BoostedJets* boostedjetfinder = H1BoostedJets::Instance();
  TObjArray* boostedjetsRec = boostedjetfinder->GetBoostedJets(H1BoostedJets::eRecLev);
  TObjArray* boostedjetsGen = boostedjetfinder->GetBoostedJets(H1BoostedJets::eHadLev);
  Int_t NumRecJets = boostedjetsRec->GetEntries();


  //!  initialize the array
  H1ArrayI MatchingIdx ;// = new H1ArrayI(); // n-gen jets. points to the matched rec jet.
  H1ArrayD MatchingDist;// = new H1ArrayD(); 
  

  // --------- loop over gen jets ------------
  for ( Int_t igen = 0 ; igen < boostedjetsGen->GetEntries()  ; igen++ ) {
    
    // ---- extend and initialize arrays ---
    MatchingIdx.Add(-1); 
    MatchingDist.Add(-1);

    // --- skip this jet if it is outside the ps-gen-cuts ---
    if ( !GenJetsPSArr[igen] ) { 
      continue;
    }

    // --- get gen jet properties ---
    TLorentzVector hadjetvec = ((H1PartJet*)boostedjetsGen->At(igen))->GetFourVector();
    Float_t dist = -1.;
    Float_t min_dist = 1000; //FLT_MAX
    Int_t   matchRec = -1;
    

    // ------ for each gen jet loop over rec jet an look for closes rec jet ---------
    for ( Int_t irec = 0 ; irec < NumRecJets ; irec++ ) {
      
      // --- skip this jet if it is outside the ps-rec cuts ---
      if ( !RecJetsPSArr[irec] ) 
	continue;
      
      // --- get rec jet properties ---
      TLorentzVector recjetvec = ((H1PartJet*)boostedjetsRec->At(irec))->GetFourVector();
      
      // --- calc dist between gen and rec jet ---
      dist = Dist(&hadjetvec, &recjetvec);

      if ( dist<min_dist && dist<fJetMatchingRadius ){
	min_dist = dist;
	matchRec = irec;
      }

    } // rec loop
    
    if ( matchRec >= 0 ) {
      MatchingIdx.AddAt( matchRec , igen); // AddAt( number, index )
      MatchingDist.AddAt( min_dist , igen );
    }

  } // gen loop
  

  Bool_t bDoubleMatched = false;
  Int_t  ReUnmatchedGenJet   = -1;
  Int_t  DoubleMatchedRecJet = -1;
  Bool_t bTrippleMatched = false;
  Int_t  TrippleMatchedRecJet = -1;
  Int_t  ReUnmatchedGenJetTripple   = -1;
  
  //! ---------------  check if one Rec Jet was matched twice. ----------------------
  for ( Int_t igen = 0 ; igen < boostedjetsGen->GetEntries()  ; igen++ ) {
     for ( Int_t igen2 = igen+1 ; igen2 < boostedjetsGen->GetEntries()  ; igen2++ ) {

	if ( ( (MatchingIdx)[igen2] == (MatchingIdx)[igen] ) && (MatchingIdx)[igen] != -1 && bDoubleMatched ){ // is this is tripple matched rec-jet?
	   TrippleMatchedRecJet = (MatchingIdx)[igen2];
	   bTrippleMatched = true;
	   cout << "JetTools::JetMatching. Info. One Rec-Jet is matched 'three' times! This is rather uncommmon."<<endl;
	   // --- we take the closest gen jet, that matches this rec jet ---
	   if ( (MatchingDist)[igen2] > (MatchingDist)[igen] ) {
	      (MatchingIdx)[igen2] = -1;
	      ReUnmatchedGenJetTripple = igen2;
	   } else {
	      (MatchingIdx)[igen] = -1;
	      ReUnmatchedGenJetTripple = igen;
	   }
	}
	
	if ( ( (MatchingIdx)[igen2] == (MatchingIdx)[igen] ) && (MatchingIdx)[igen] != -1 ) { // is this a double-matched rec-jet??
	   DoubleMatchedRecJet = (MatchingIdx)[igen2];
	   //cout << "grepme: Two Gen Jets are matching one Rec Jet. ReRun with this exception the JetMatching. ToDo: clearify the double matching." << endl;
	   bDoubleMatched = true;
	   // --- we take the closest gen jet, that matches this rec jet ---
	   if ( (MatchingDist)[igen2] > (MatchingDist)[igen] ) {
	      (MatchingIdx)[igen2] = -1;
	      ReUnmatchedGenJet = igen2;
	   } else {
	      (MatchingIdx)[igen] = -1;
	      ReUnmatchedGenJet = igen;
	   }
	}
    }
  }


  static Int_t countDoubleMatched = 0;
  static Int_t countPossibleRematching = 0;

  if ( bDoubleMatched ){
    //cout << "grepme: testing: we look now, if the reunmatched jet could find another jet, and then, if this other jet would also already be matched." << endl;
    countDoubleMatched++;
    // --------- loop over gen jets ------------
    Int_t igen = ReUnmatchedGenJet;

    // --- get gen jet properties ---
    TLorentzVector hadjetvec = ((H1PartJet*)boostedjetsGen->At(igen))->GetFourVector();
    Float_t dist = -1.;
    Float_t min_dist = 1000; //FLT_MAX
    Int_t   matchRec = -1;

    // ------ for each gen jet loop over rec jet and look for closes rec jet ---------
    for ( Int_t irec = 0 ; irec < NumRecJets ; irec++ ) {
      // --- skip this is the rec jet which is already matched by another (closer) jet ---
      if ( irec == DoubleMatchedRecJet ) {
	continue;  }
      // --- skip this jet if it is outside the ps-rec cuts ---
      if ( !RecJetsPSArr[irec] ) {
	continue;   }
      // --- get rec jet properties ---
      H1PartJet* jet_rec = (H1PartJet*)boostedjetsRec->At(irec);
      TLorentzVector recjetvec = ((H1PartJet*)boostedjetsRec->At(irec))->GetFourVector();
      // --- calc dist between gen and rec jet ---
      dist = Dist(&hadjetvec, &recjetvec);
      if ( dist<min_dist && dist<fJetMatchingRadius ){
	min_dist = dist;
	matchRec = irec;
      }
    } // rec loop
    if ( matchRec >= 0 ) {
      //cout << "greme: well: the reunmatched gen jet could find another rec jet within jm-radius." << endl;
      countPossibleRematching++;

      //! ---------------  check if this second try is also already matched , otherwise assign this match ----------------------
      for ( Int_t igen2 = 0 ; igen2 < boostedjetsGen->GetEntries()  ; igen2++ ) {
	if ( !( (MatchingIdx)[igen2] == matchRec  && (MatchingIdx)[igen2] != -1) ) { 
	  MatchingIdx.AddAt( matchRec , igen ); // AddAt( number, index )
	  MatchingDist.AddAt( min_dist , igen );
	}
      }
      
    }
  }
  
//   // --- just a cout of the status ---
//   if ( bDoubleMatched ){
//     cout << "grepme: Status: Double: "<<countDoubleMatched << "  Dochnicht: " <<countPossibleRematching << "  Complx: " << countMoreComplex << "\tTripleMatching: " <<countTripleMatching <<  endl;
//   } 

  return MatchingIdx;

}


//_________________________________________________________________________________________________________________


Double_t JetTools::Dist(TLorentzVector* jet1, TLorentzVector* jet2)
{
  // return eta-phi distance between jet1 and jet2

  double phi1 = jet1->Phi();
  double phi2 = jet2->Phi();
  double eta1 = jet1->Eta();
  double eta2 = jet2->Eta();
  
  Double_t Dphi = phi1 - phi2;
  Double_t pi   = TMath::Pi();

  if (TMath::Abs(Dphi) > pi){
    if (Dphi>0) Dphi-=2*pi;
    else Dphi+=2*pi;
  }

  //double dist = TMath::Sqrt( TMath::Power( phi1 - phi2, 2) + TMath::Power( eta1 - eta2, 2) );
  Double_t dist = TMath::Sqrt( TMath::Power( Dphi, 2) + TMath::Power( eta1 - eta2, 2) );
 
  return dist;

}


Float_t JetTools::GenElecPhotDist()
{
  // return distance in eta-phi between the generated 
  // electron and the generated QED photon
  // returns 0 if no FSR event

  Float_t dist = 0;
  
  static H1ShortPtr GenRadType("GenRad");
  
  if ((*GenRadType)!=2){
    return FLT_MAX;
  }
  
  // get the photon
  static H1FloatPtr EnPhPtr("GenEnPhoton"); 
  static H1FloatPtr ThPhPtr("GenThPhoton"); 
  static H1FloatPtr PhPhPtr("GenPhPhoton"); 
  Float_t phEn = (*EnPhPtr);
  Float_t phTh = (*ThPhPtr);
  Float_t phPh = (*PhPhPtr);
  
  Float_t phPx = phEn*TMath::Sin(phTh)*TMath::Cos(phPh);
  Float_t phPy = phEn*TMath::Sin(phTh)*TMath::Sin(phPh);
  Float_t phPz = phEn*TMath::Cos(phTh);
  TLorentzVector Photon(phPx, phPy, phPz, phEn);

  // get the generated (uncombined) electron
  static H1FloatPtr EPtr("GenEnElecUncombined");
  static H1FloatPtr ThPtr("GenThElecUncombined");
  static H1FloatPtr PhPtr("GenPhElecUncombined");
  Float_t elecE = (*EPtr);
  Float_t elecTh = (*ThPtr);
  Float_t elecPh = (*PhPtr);
  Float_t elecP  = elecE; // ignore the mass

  Float_t elecPx = elecP*TMath::Sin(elecTh)*TMath::Cos(elecPh); 
  Float_t elecPy = elecP*TMath::Sin(elecTh)*TMath::Sin(elecPh);
  Float_t elecPz = elecP*TMath::Cos(elecTh); 	    
  TLorentzVector Electron(elecPx, elecPy, elecPz, elecE);

  // calculate the distance:
  dist = Dist(&Photon, &Electron);

  return dist;

}


TLorentzVector JetTools::CalcElecForBoost(Float_t q2, Float_t y, Float_t Phi)
{
  // Calculate the scattered electron four-vector for the boost 
  // to the Breit frame. Reconstruct the vector from Q2, y and Phi_e,
  // to be able to take a different reconstruction method

  Double_t E0 = 27.6; // electron beam energy

  Double_t px = 0;
  Double_t py = 0;
  Double_t pz = 0;
  Double_t Ee = 0;
  
  if ((q2>0) && (y>0) && (y<1)){

    Ee = q2 / (4*E0) + E0*(1-y);
    
    Double_t b = 4*E0*E0*(1-y)/q2;
    Double_t Theta = TMath::ACos( (1-b)/(1+b) );
    
    px = Ee*TMath::Sin(Theta)*TMath::Cos(Phi);
    py = Ee*TMath::Sin(Theta)*TMath::Sin(Phi);
    pz = Ee*TMath::Cos(Theta);

  }

  TLorentzVector ElecBoostVec(px, py, pz, Ee);
  
  return ElecBoostVec;

}

Float_t JetTools::GetInvMassFromJetParts(H1PartJet* jet)
{
  // return invariant mass of all jet constituents
  TLorentzVector sum(0,0,0,0);
  for (Int_t i=0; i<jet->GetNumOfParticles();++i){
    const H1Part* part = jet->GetParticle(i);
    sum += part->GetFourVector();
  }
  return sum.Mag();
}
