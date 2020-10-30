#include <iostream>
#include <iomanip>
#include "TDetectQedc.h"
#include "H1Mods/H1PartMCArrayPtr.h"
#include "H1Mods/H1PartMC.h"

using namespace std;

ClassImp(TDetectQedc)

TDetectQedc::TDetectQedc(H1PartMCArrayPtr &mcpart) {
   fIelectron=-2;
   fNisr=0;
   fNfsr=0;
   fIphotonQedc=-2;
   fELIE=kELIE_DEFAULT;
   fAGMA=kAGMA_DEFAULT;
   fAGMI=kAGMI_DEFAULT;
   fELIG=kELIG_DEFAULT;
   fPTMA=kPTMA_DEFAULT;
   fVMIN=kVMIN_DEFAULT;
   fWMIN=kWMIN_DEFAULT;
   fWMAX=kWMAX_DEFAULT;
   fACO =kACO_DEFAULT;
   fNcluster=0;
   fNstring=0;
   fMstring=0.0;
   fMcluster=0.0;
   fMesonPdg=0;
   Int_t vmParent=-1;

   Double_t *eppz=new Double_t[mcpart.GetEntries()];
   Int_t *ieppz=new Int_t[mcpart.GetEntries()];
   Int_t neppz=0;
   for(Int_t i=0;i<mcpart.GetEntries();i++) {
      H1PartMC *p=mcpart[i];
      int pdg=abs(p->GetPDG());
      int status=p->GetStatus();
      if(status<0) continue;
      int parent_status=-2;
      int parent_pdg=-1;
      if(p->GetMother1()>=0) {
         parent_status=mcpart[p->GetMother1()]->GetStatus();
         parent_pdg=abs(mcpart[p->GetMother1()]->GetPDG());
      }
      if(status==201) {
         fEPbeam += p->GetFourVector();
      }
      if(pdg==11) {
         // electron or positron
         if((status==0)&&(parent_status==201) && (fIelectron<0)) {
            fElectron=p->GetFourVector();
            fIelectron=i;
         }
      } else if(pdg==22) {
         // photon
         if((status==202)&&(fIelectron==p->GetMother1())) {
            fPhoton[PHOTON_FSR] += p->GetFourVector();
            fNfsr++;
         } else if((status==202)&&(parent_status==201)) {
            fPhoton[PHOTON_ISR] += p->GetFourVector();
            fNisr++;
         } else if((status==202)||
                   ((status==0) && (parent_pdg<20))) {
            if(fIphotonQedc<0) {
               fPhoton[PHOTON_QEDC]=p->GetFourVector();
               fIphotonQedc=i;
            }
         }
      }
      if(pdg==91) {
         fNcluster++;
         fMcluster += p->GetMass();
      }
      if(pdg==92) {
         fNstring++;
         fMstring += p->GetMass();
      }
      // detect DIFFVM meson
      if((pdg>100)&&(vmParent>=0)&&(vmParent==p->GetMother1())&&
         (pdg==parent_pdg)) {
         fMesonPdg=p->GetPDG();
         fMeson=p->GetFourVector();
      }
      if((pdg>100)&&(parent_pdg==22)) {
         vmParent=i;
      }
      // detect DJANGO meson here...
      // create list of E-Pz of stable particles to detect diffraction
      if(!status) {
         if(fIelectron==i) continue; // exclude scattered electron
         eppz[neppz]=p->GetE()-p->GetPz();
         ieppz[neppz]=i;
         ++neppz;
      }
   }
   // order list in eppz
   for(int inc=(neppz+1)/2;inc;inc /= 2) {
      for(int i=inc;i<neppz;i++) {
         int itmp=ieppz[i];
         int j;
         for(j=i;(j>=inc)&&(eppz[ieppz[j-inc]]>eppz[itmp]);j -=inc) {
            ieppz[j]=ieppz[j-inc];
         }
         ieppz[j]=itmp;
      }
   }
   // add up particles from both sides,
   // such that the sum of both masses is minimized
   Int_t i1=0,i2=neppz-1;
   while(i1<=i2) {
      if((fX+mcpart[ieppz[i1]]->GetFourVector()).M()+fY.M()<
         (fY+mcpart[ieppz[i2]]->GetFourVector()).M()+fX.M()) {
         fX += mcpart[ieppz[i1++]]->GetFourVector();
      } else {
         fY += mcpart[ieppz[i2--]]->GetFourVector();
      }
   }
   //cout<<neppz<<" "<<fX.M()<<" "<<fY.M()<<"\n";
 
   delete [] eppz;
   delete [] ieppz;


   // check cuts of all pairings
   for(Int_t i=0;i<PHOTON_NTYPE;i++) {
      fCutsQedc[i]=GetQedcCuts(fElectron,fPhoton[i]);
   }
   // calculate more 4-vectors
   fW=fEPbeam-fElectron -fPhoton[PHOTON_QEDC] -fPhoton[PHOTON_FSR];

   // mass of all final-state particles
   // which are neither the electron nor the QEDC photon
   fInvis=fW;
   // mass of the outgoing hadronic system
   fW -= fPhoton[PHOTON_ISR];
}

void TDetectQedc::Print(H1PartMCArrayPtr &mcpart) {
   cout<<"TDetectQedc::Print npart="<<mcpart.GetEntries()
       <<" nString="<<fNstring
       <<" nCluster="<<fNcluster
       <<"\n";
   if(fMesonPdg) {
      cout<<"Meson pdg="<<fMesonPdg
          <<" Pt="<<fMeson.Pt()
          <<" Pz="<<fMeson.Pz()
          <<" E="<<fMeson.E()
          <<" M="<<fMeson.M()<<"\n";
   }
   for(int i=0;i<mcpart.GetEntries();i++) {
      H1PartMC *p=mcpart[i];
      int pdg=abs(p->GetPDG());
      int status=p->GetStatus();
      if(status<0) continue;
      int parent_status=-2;
      if(p->GetMother1()>=0) {
         parent_status=mcpart[p->GetMother1()]->GetStatus();
      }
      /* cout<<setw(3)<<i<<" "<<setw(5)<<pdg<<" "<<setw(3)<<status
         <<" "<<setw(3)<<p->GetMother1()<<" "<<setw(3)<<p->GetMother2()
         <<" "<<setw(3)<<p->GetDaughter1()<<" "<<setw(3)<<p->GetDaughter2()
         <<" ("<<p->GetE()<<","<<p->GetPt()<<","<<p->GetPz()<<")"
         <<"\n"; */
      if(i==fIelectron) cout<<" --- Next is eqedc ---\n";
      if(i==fIphotonQedc) cout<<" --- Next is gqedc ---\n";
      if(pdg==22) {
         if((status==202)&&(fIelectron==p->GetMother1())) {
            cout<<" --- Next is FSR ---\n";
         } else if((status==202)&&(parent_status==201)) {
            cout<<" --- Next is ISR ---\n";
         }
      }
      if(!i) {
         p->Print("HEAD");
      }
      cout<<setw(3)<<i<<setw(6)<<p->GetPDG()<<" ";
      p->Print();
   }
}

Bool_t TDetectQedc::IsElectronFound(void) const {
   return fIelectron>=0;
}

Bool_t TDetectQedc::IsPhotonQedcFound(void) const {
   return fIphotonQedc>=0;
}

Float_t const TDetectQedc::kELIE_DEFAULT=1.5;
Float_t const TDetectQedc::kAGMA_DEFAULT=178.;
Float_t const TDetectQedc::kAGMI_DEFAULT=3.6;
Float_t const TDetectQedc::kELIG_DEFAULT=1.5;
Float_t const TDetectQedc::kPTMA_DEFAULT=20.;
Float_t const TDetectQedc::kVMIN_DEFAULT=15.;
Float_t const TDetectQedc::kWMIN_DEFAULT=1.5;
Float_t const TDetectQedc::kWMAX_DEFAULT=310.;
Float_t const TDetectQedc::kACO_DEFAULT=50.;

Int_t TDetectQedc::GetQedcCuts
(TLorentzVector const &electron,TLorentzVector const &photon) const {
   Int_t r=0;
   if(electron.E()<fELIE)  r |= 0x00000001;
   Double_t theta=electron.Theta()*180./M_PI;
   if(theta>fAGMA)         r |= 0x00000002;
   if(theta<fAGMI)         r |= 0x00000004;
   if(photon.E()<fELIG)    r |= 0x00000010;
   theta=photon.Theta()*180./M_PI;
   if(theta>fAGMA)         r |= 0x00000020;
   if(theta<fAGMI)         r |= 0x00000040;
   TLorentzVector s=electron+photon;
   if(s.Pt()>fPTMA)        r |= 0x00000080;
   if(s.E()<fVMIN)         r |= 0x00000100;
   if(s.M()<fWMIN)         r |= 0x00000200;
   if(s.M()>fWMAX)         r |= 0x00000400;
   Double_t dphi=fabs(remainder(electron.Phi()+M_PI-photon.Phi(),2.*M_PI)
                      *180./M_PI);
   if(dphi>fACO)           r |= 0x00000800;
   return r;
}

Bool_t TDetectQedc::IsQedcEvent(void) const {
   for(Int_t i=0;i<PHOTON_NTYPE;i++) {
      // any of the pairs inside all QEDC cuts?
      if(!fCutsQedc[i]) return kTRUE;
   }
   // none of the pairs matches all QEDC cuts
   return kFALSE;
}
