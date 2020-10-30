// (c) MPI 2020
#include "H2020HistManager.h"
 
// C++ includes
#include <fstream>
#include <set>

// Root includes
 
// H1 includes

using namespace std;


// _______________________________________________________ //
//! Constructor
H2020HistManager::H2020HistManager(string HMname) : fHMname(HMname) {
   TH1::AddDirectory(false);
}

// _______________________________________________________ //E
H2020HistManager::~H2020HistManager() {

}

// _______________________________________________________ //
//! Write histograms to gDirectory
void H2020HistManager::Write() {
   // sort by name
   cout<<"H2020HistManager::Write. Writing histograms into "<<gDirectory->GetPath()<<endl;
   map<string,TH1*> hmap;
   for ( auto [k1, m1] : fHistmap ) {
      for ( auto [k2,th] : m1 ) {
         hmap[th->GetName()] = th;
      }
   }
   for ( auto [k1, th] : hmap ) th->Write();
}
