#ifndef H_ALIGNMENT
#define H_ALIGNMENT

#include "H1Steering/H1Steer.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TLorentzVector.h"

using namespace std;

class Alignment : public H1Steer {
 public:
  Alignment(void);
  void SetOrigin(const char *in);
  void SetPhi(Double_t p);
  void SetOctant(Double_t o);
  void SetZ0(Double_t o);
  inline void SetZOrigin(Double_t const &z) { fOrigin[2]=z; }
  inline void SetOrigin(TVector3 const &pos) { fOrigin=pos; }
  //inline void SetOrigin(Double_t const &x,Double_t const &y,Double_t const &z) {
  //  fOrigin.SetXYZ(x,y,z); }

  void ImportParameters(Alignment const *a);
  //void SetPar(int i,Double_t const &x);
  //Double_t GetPar(int i) const;
  TVector3 Transform(TVector3 const &x) const;
  TVector3 const &GetOrigin(void) const { return fOrigin; }
  Double_t GetZ0(void) const { return fZ0; }
  inline Double_t GetPhi(void) const { return fPhi; }
  // TLorentzVector Rotate(TLorentzVector const &x) const;
  /* enum { SHIFT_X=0,SHIFT_Y=1,SHIFT_Z=2,ROT_PHI=3,ROT_OCTANT=4,
     Z0=5,NPAR=6 }; */

  void PrintSteer(ostream &out) const;

  inline bool IsZero(void) const {
    return (fOrigin[0]==0.0)&&(fOrigin[1]==0.0)&&(fOrigin[2]==0.0)&&
      (fPhi==0.0)&&(fOctant==0.0);
  }
 protected:
  TVector3 fOrigin;
  Double_t fPhi;
  Double_t fOctant;
  Double_t fZ0;
  ClassDef(Alignment,0)
};

#endif
