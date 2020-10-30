#include "SpacLinearity.h"
#include "H1Steering/H1SteerManager.h"
#include <TObjString.h>
#include <TArrayD.h>
#include <iostream>
#include <cmath>

using namespace std; 

ClassImp(SpacLinearity)

SpacLinearity::SpacLinearity(void) {
   alreadyPrinted=false;
   for(int i=0;i<5;i++) {
      fAlpha[i]=0.0;
   }
   for(int i=0;i<81;i++) {
      fSuperX[i]=0.;
      fSuperY[i]=0.;
   }
}

void SpacLinearity::SetAlpha(const char *in) {
   TArrayD tmp;
   StringToArray(in,tmp);
   for(int i=0;(i<5)&&(i<tmp.GetSize());i++) {
      fAlpha[i]=tmp[i];
   }
}

void SpacLinearity::SetSuperX(const char *in) {
   TArrayD tmp;
   StringToArray(in,tmp);
   for(int i=0;i<81;i++) {
      fSuperX[i]=tmp[i];
   }
}

void SpacLinearity::SetSuperY(const char *in) {
   TArrayD tmp;
   StringToArray(in,tmp);
   for(int i=0;i<81;i++) {
      fSuperY[i]=tmp[i];
   }
}

void SpacLinearity::SetSuperXY(int i,double x, double y) {
   fSuperX[i] = x;
   fSuperY[i] = y;
}

void SpacLinearity::SetAlphaByIndex(int i, Float_t a) {
   fAlpha[i]=a;
}

TVector3 SpacLinearity::GetPosCorr(TVector3 const &pos,Float_t const *grid,
                                   Float_t const *offset) const {
   TVector3 r(pos);
   // cell-by-cell linearity correction
   double rr=r.Pt();
   for(int i=0;i<2;i++) {
      double arg=(r[i]-offset[i])/grid[i]*2.*M_PI;
      r[i] += fAlpha[0]*(1.+fAlpha[1]*rr+fAlpha[2]*fabs(r[i]))
         *(sin(arg)+fAlpha[3]*sin(2.*arg)+fAlpha[4]*sin(3.*arg));
   }

   // supermodule position correction

   // supermodule number
   int isup[2][2];
   for(int ic=0;ic<2;ic++) {
      for(int ix=0;ix<2;ix++) {
         isup[ic][ix]=(int)(pos[ic]/grid[ic]/4.+4.+ix);
         if(isup[ic][ix]<0) isup[ic][ix]=0;
         if(isup[ic][ix]>8) isup[ic][ix]=8;
      }
   }
   // correct by linear interpolation along both directions
   double f[2];
   for(int ic=0;ic<2;ic++) {
      f[ic]=pos[ic]/grid[ic]/4.+4.-isup[ic][0];
   }
    
   r[0] +=
      (fSuperX[isup[0][0]+9*isup[1][0]]*(1.-f[0])+
       fSuperX[isup[0][1]+9*isup[1][0]]*(f[0]   ))*(1.-f[1])+
      (fSuperX[isup[0][0]+9*isup[1][1]]*(1.-f[0])+
       fSuperX[isup[0][1]+9*isup[1][1]]*(f[0]   ))*f[1];
   r[1] +=
      (fSuperY[isup[0][0]+9*isup[1][0]]*(1.-f[0])+
       fSuperY[isup[0][1]+9*isup[1][0]]*(f[0]   ))*(1.-f[1])+
      (fSuperY[isup[0][0]+9*isup[1][1]]*(1.-f[0])+
       fSuperY[isup[0][1]+9*isup[1][1]]*(f[0]   ))*f[1];

   return r;
}

void SpacLinearity::PrintSteer(ostream &out,bool forcePrint) const {
   if(forcePrint || !alreadyPrinted) {
      out<<"SpacLinearity(\""<<GetName()<<"\") {\n";
      PrintN(out,"Alpha",fAlpha,5);
      PrintN(out,"SuperX",fSuperX,81,9);
      PrintN(out,"SuperY",fSuperY,81,9);
      out<<"}\n\n";
      alreadyPrinted=true;
   }
}

void SpacLinearity::PrintN(ostream &out,char const *name,Float_t const *f,
                           int n,int npl) {
   out<<"  f"<<name<<"=\"";
   for(int i=0;i<n;i++) {
      out<<(f[i]);
      if(i!=(n-1)) {
         out<<",";
         if(npl && ((i+1)%npl ==0)) out<<"\n     ";
      }
   }
   out<<"\";\n";
}
