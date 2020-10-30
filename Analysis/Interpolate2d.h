#include <TVectorD.h>
#include <TH2.h>
#include <iostream>

class Interpolate2d {
 public:
   void ImportHistogram(const TH2 *input,bool useLogZ);
   bool ReadData(const char *coeffFile);
   bool WriteData(const char *coeffFile);
   bool ReadData(std::istream &in);
   bool WriteData(std::ostream &out);
   double Interpolate(double x,double y,
                      bool useLogX,bool useLogY,double extraPol=0.0);
 protected:
   bool fLogZ;
   TVectorD fX,fY;
   TMatrixD fZ;
};
