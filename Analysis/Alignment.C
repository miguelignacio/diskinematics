#include "Alignment.h"
#ifdef PROFILE
#include "Profile.h"
#endif

#include <iostream>
#include <TArrayD.h>
#include <cmath>

ClassImp(Alignment)

Alignment::Alignment(void) {
  fPhi=0.0;
  fOctant=0.0;
  fZ0=0.0;
}

void Alignment::SetOrigin(const char *in) {
  TArrayD tmp;
  StringToArray(in,tmp);
  for(int i=0;(i<3)&&(i<tmp.GetSize());i++) {
    fOrigin[i]=tmp[i];
  }
}

void Alignment::SetPhi(Double_t p) {
  fPhi=p;
}

void Alignment::SetOctant(Double_t o) {
  fOctant=o;
}

void Alignment::SetZ0(Double_t z0) {
  fZ0=z0;
}

void Alignment::ImportParameters(Alignment const *a) {
  fOrigin=a->fOrigin;
  fPhi=a->fPhi;
  fOctant=a->fOctant;
  fZ0=a->fZ0;
}

TVector3 Alignment::Transform(TVector3 const &x) const {
#ifdef PROFILE
  Profile::start("Alignment::Transform");
#endif
  TVector3 r=x;
  Double_t phi0=r.Phi();
  r.SetPhi(phi0+fPhi+fOctant*fmod(phi0*180./M_PI+180.,45.));
  r += fOrigin;
#ifdef PROFILE
  Profile::stop("Alignment::Transform");
#endif
  return r;
}

void Alignment::PrintSteer(ostream &out) const {
  out<<"Alignment(\""<<GetName()<<"\") {\n";
  if(fPhi!=0.0) out<<"  fPhi="<<fPhi<<";\n";
  if(fOctant!=0.0) out<<"  fOctant="<<fOctant<<";\n";
  if((fOrigin[0]!=0.0)||(fOrigin[1]!=0.0)||(fOrigin[2]!=0.0)) {
    out<<"  fOrigin=\""<<fOrigin[0]<<","<<fOrigin[1]<<","<<fOrigin[2]<<"\";\n";
  }
  if(fZ0!=0.0) {
    out<<"  fZ0="<<fZ0<<";\n";
  }
  out<<"}\n";
}
