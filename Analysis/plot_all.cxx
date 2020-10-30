
// ___________________________________________________________ //
enum kCHAIN {
   kData,
   kDjango,
   kRapgap,
   kPythia,
   kDVCS,
   kCompton,
   kCCBkgd,
   kDjBkgd,
   kUndefined = 99,
};
map<kCHAIN,string> gChainMap{
   { kData,    "Data" },
   { kDjango,  "Django" },
   { kRapgap,  "Rapgap" },
   { kPythia,  "Pythia" },
   { kDVCS,    "DVCS" },
   { kCompton, "Compton" },
   { kCCBkgd,  "CCBkgd" },
   { kDjBkgd,  "DjBkgd" }, 
      };

// ___________________________________________________________ //
void SetMainStyle() {
   gStyle->SetTextSize(0.06);
   gStyle->SetLegendTextSize(0.06);

   gStyle->SetPadLeftMargin(0.20);
   gStyle->SetPadRightMargin(0.10);
   gStyle->SetPadTopMargin(0.06);
   gStyle->SetPadBottomMargin(0.19);
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   gROOT->ForceStyle();
}


// ___________________________________________________________ //
//!  Set plotting style for a 1D histogram
void SetStyle_TH1_MC(TH1* hist) {
   hist->SetLineWidth(2);
   hist->SetLineStyle(1);
}

void SetStyle_TH1_Data(TH1* hist) {
   hist->SetMarkerColor(kBlack);
   hist->SetMarkerStyle(8);
   hist->SetMarkerSize(1.0);
}
void SetStyle_TH1(TH1* hist, kCHAIN chain) {
   if ( hist==NULL ) return ;// ausnahmsweise
   if ( chain == kData )   SetStyle_TH1_Data(hist);
   else {
      SetStyle_TH1_MC(hist);
      if ( chain == kDjango )       {
         hist->SetLineColor(kAzure-2);
      }
      else if ( chain == kRapgap )  {
         hist->SetLineColor(kRed);
      }
      else if ( chain == kPythia )  {
         hist->SetLineColor(kGreen+1);
      }
   }
}
void SetStyle_Main(TH1* hist) {
   hist->SetTitleOffset(1.5,"x");
   hist->SetTitleOffset(1.5,"y");
   hist->SetTitleSize(0.06,"x");
   hist->SetTitleSize(0.06,"y");
   hist->SetLabelSize(0.06,"x");
   hist->SetLabelSize(0.06,"y");
   hist->SetLabelOffset(0.006,"x");
   hist->SetLabelOffset(0.008,"y");
   hist->SetStats(0);

}


// ___________________________________________________________ //
//! return the chains in the file:
//! DataEplus0304, Django, Rapgap, etc...
//! 
tuple<string,map<kCHAIN,string> > get_chains(TFile* file, bool DoPrint = false) {
   map<kCHAIN,string> ret;
   const int nEntries = file->GetListOfKeys()->GetEntries();
   for ( int iO = 0 ; iO<nEntries ; iO++ ){
      kCHAIN KKChn = kUndefined;
      string chainname = file->GetListOfKeys()->At(iO)->GetName();
      for ( auto [ kch , sch ] : gChainMap ) {
         if ( chainname.find(sch) == 0 ) KKChn = kch;
      }
      if ( KKChn==kUndefined ) {
         cout<<"ERROR! Cannot identify chain from directory-name: "<<chainname<<endl;
         exit(1);
      }
      ret[KKChn] = chainname;
   }
   if ( ret.empty() ) { 
      cout<<"ERROR! No chains were found. Exiting."<<endl;
      exit(1);
   }
   // get main chain:
   // Django, Rapgap, Data
   string mainchain;
   for ( auto [c,cn] : ret ) 
      if ( cn.find("Django") == 0 ) mainchain = cn; 
   if ( mainchain == "" )
      for ( auto [c,cn] : ret ) 
         if ( cn.find("Rapgap") == 0 ) mainchain = cn; 
   if ( mainchain == "" )
      for ( auto [c,cn] : ret ) 
         if ( cn.find("Data") == 0 ) mainchain = cn; 
   //ret.erase(std::remove(ret.begin(), ret.end(), mainchain), ret.end());
   // some printout
   // if (DoPrint) {
   //    cout<<"\nMain chain:        "<<mainchain<<endl;
   //    cout<<"Further chains: ";
   //    for ( auto c : ret ) cout<<"   "<<c;
   //    cout<<endl<<endl;
   // }
   return {mainchain,ret}; 
}

// ___________________________________________________________ //
//! return list of all directories (aka histmanagers)
vector<string> get_directories(TFile* file) {
   auto [mainchain, chains] = get_chains(file);
   file->cd(mainchain.c_str());
   const int nEntries = gDirectory->GetListOfKeys()->GetEntries();
   cout<<"\nDirectories: "<<endl;
   vector<string> ret;
   for ( int iO = 0 ; iO<nEntries ; iO++ ){
      ret.push_back(gDirectory->GetListOfKeys()->At(iO)->GetName());
      cout<<"  - "<<ret.back()<<endl;
   }
   return ret;
}

// ___________________________________________________________ //
//! return list of all histograms in directory
vector<string> get_histolist(TFile* file, const string& dirname) {
   auto [mainchain, chains] = get_chains(file);
   file->cd((mainchain+"/"+dirname).c_str());
   int nEntries = gDirectory->GetListOfKeys()->GetEntries();
   vector<string> ret;
   for ( unsigned int iO = 0 ; iO<nEntries; iO++ ){
      ret.push_back(gDirectory->GetListOfKeys()->At(iO)->GetName());
   }
   return ret;
}


// ___________________________________________________________ //
//!  helper function to get histogram from fil
//!  implement checks or furhter type-casts if needed
TH1* get_hist_from_file(TFile* file, const string& chain, const string& directory, const string& histname) {
   TH1* hist = file->Get<TH1>((chain+"/"+directory+"/"+"/"+histname).c_str());
   return hist;
}




// ___________________________________________________________ //
//! return list of all histograms in directory
void plot_directory_default(TFile* file, const string& directory, const TString& outps) {

   // set style and all that
   TCanvas* cc = new TCanvas("cc","plots",1200,800);
   SetMainStyle();
   cc->Divide(2,2);
   cc->Print(outps+"["); //// write canvas and keep the ps file open 
   
   

   // names of all histograms, and all chains
   vector<string> allhistos = get_histolist(file,directory);
   auto [mainchain, chains] = get_chains(file);

   int iPad=1;                                                                                                          
   for ( string histname : allhistos ) {
      TH1* hData   = chains.count(kData)   ? get_hist_from_file(file,chains[kData],directory,histname)   :  NULL;
      TH1* hRapgap = chains.count(kRapgap) ? get_hist_from_file(file,chains[kRapgap],directory,histname) :  NULL;
      TH1* hDjango = chains.count(kDjango) ? get_hist_from_file(file,chains[kDjango],directory,histname) :  NULL;

      TH1* hBackground = NULL;// todo, get all background, and also add them to Django and Rapgap

      // get main-histogram
      TH1* hMain = NULL;
      if ( mainchain.find("Data") == 0 ) hMain = hData;
      else if ( mainchain.find("Django") == 0 ) hMain = hDjango;
      else if ( mainchain.find("Rapgap") == 0 ) hMain = hRapgap;
      string title = hMain->GetTitle();
      if ( title == "" ) title = hMain->GetName();
      title = "["+directory+"] "+ title;
      hMain->SetTitle(title.c_str());

      // set style
      SetStyle_Main(hMain);
      SetStyle_TH1(hData,kData);
      SetStyle_TH1(hBackground,kPythia);
      SetStyle_TH1(hDjango,kDjango);
      SetStyle_TH1(hRapgap,kRapgap);


      cc->cd(iPad);
      // plot a 1D-histogram
      if ( hMain->InheritsFrom("TH1D") || hMain->InheritsFrom("TH1F")  || hMain->InheritsFrom("TH1I")  ) {
         double max = hMain->GetMaximum();
         if ( title.find("_lx")!=string::npos || string(hMain->GetName()).find("_lx") !=string::npos ) gPad->SetLogx();
         if ( title.find("_ly")!=string::npos || string(hMain->GetName()).find("_ly") !=string::npos ) gPad->SetLogy();
         if ( gPad->GetLogy() )  hMain->SetMaximum(max*40);
         else                    hMain->SetMaximum(max*1.4);
         if ( gPad->GetLogy() ) hMain->SetMinimum(0.7);
         else                   hMain->SetMinimum(0.);
         

         hMain->DrawClone("hist");
         if ( hBackground ) hBackground->DrawClone("histsame");
         if ( hRapgap ) hRapgap->DrawClone("histsame");
         if ( hDjango ) hDjango->DrawClone("histsame");
         if ( hData   ) hData->DrawClone("PEsame");

         TLegend* l = new TLegend(0.65,0.70,0.98,0.93,"","brNDC");
         l->SetTextSize(0.04);
         if ( hData   )     l->AddEntry(hData,  "Data","PE");
         if ( hDjango )     l->AddEntry(hDjango,"Django","L");
         if ( hRapgap )     l->AddEntry(hRapgap,"Rapgap","L");
         if ( hBackground ) l->AddEntry(hBackground,"Backgrounds","L");
         l->SetFillStyle(0);
         l->SetBorderSize(0);
         l->DrawClone();

         // if ( histname.Contains("_den") ) {
         //    TString numname = histname.ReplaceAll("_den","_num");
         //    ...
         //       }
         gPad->RedrawAxis();
         
      }
      // plot a 1D-histogram
      else if ( hMain->InheritsFrom("TH2D") || hMain->InheritsFrom("TH2F")  || hMain->InheritsFrom("TH2I")  ) {
         gPad->SetLogz();
         cout<<"not yet implemented."<<endl;
         exit(1);
      }
      else {
         cout<<"Warning. Object not recognized! Name: "<<hMain->GetName()<<endl;
         continue;
      }

      // write, clear canvas
      if ( (++iPad)>4 ) {
         iPad=1;
         cc->Print(outps);
         cc->Clear();
         cc->Divide(2,2);
      }

   } // end of loop over all histograms

   cc->Print(outps);
   cc->Print(outps+"]"); // canvas is added to "outps" and ps file is closed
   
}





// ___________________________________________________________ //
//! main
void plot_all(){
   
   string inputfile = "merged.root";
   //string inputfile = "RootOut/DataEplus0304.root";
   TFile* file = TFile::Open(inputfile.c_str(),"READ");

   // print chains in file
   get_chains(file,true);

   // loop over all directories and plot them
   vector<string> directories = get_directories(file);
   for ( string directory : directories ) {
      cout<<" -------------------------------------------------- "<<endl;
      cout<<" [main]  Plotting directory:  "<<directory<<endl;
      cout<<"         using:    plot_directory_default"<<endl;
      string outps = "Plots/"+directory+".ps";
      cout<<"         output:   "<<outps<<endl;
      plot_directory_default(file,directory,outps);
      cout<<" -------------------------------------------------- "<<endl;
   }

   string outps = "plots.ps";

}

int main() { 
   plot_all();
   return 0;
}
