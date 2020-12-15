
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
   kGrape,
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
   { kGrape,   "Grape" },
      };

// ___________________________________________________________ //
void SetMainStyle() {
   gStyle->SetTextSize(0.006);
   gStyle->SetLegendTextSize(0.006);

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
   hist->SetTitleOffset(1.7,"y");
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
      //if ( gDirectory->GetListOfKeys()->At(iO)->InheritsFrom( "TTree" ) ) { 
      //   continue;
      //}
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


void plot_directory_acceptance(TFile* file, const string& directory, const TString& outps) {
   // set style and all that
   TCanvas* cc = new TCanvas("cc","plots",1200,1600);
   SetMainStyle();
   cc->Divide(2,2);
   cc->Print(outps+"["); //// write canvas and keep the ps file open


   // names of all histograms, and all chains
   vector<string> allhistos = get_histolist(file,directory);
   auto [mainchain, chains] = get_chains(file);

   //Get 1D histos
   TH1* HistoGenEvents_1D  = chains.count(kDjango) ? get_hist_from_file(file,chains[kDjango],directory,allhistos.at(0)) : NULL;
   TH1* HistoRecEvents_1D  = chains.count(kDjango) ? get_hist_from_file(file,chains[kDjango],directory,allhistos.at(1)) : NULL;
   TH1* HistoStayEvents_1D = chains.count(kDjango) ? get_hist_from_file(file,chains[kDjango],directory,allhistos.at(2)) : NULL;


   //Get 3D histos
   TH1* HistoGenEvents  = chains.count(kDjango) ? get_hist_from_file(file,chains[kDjango],directory,allhistos.at(3)) : NULL;
   TH1* HistoRecEvents  = chains.count(kDjango) ? get_hist_from_file(file,chains[kDjango],directory,allhistos.at(4)) : NULL;
   TH1* HistoStayEvents = chains.count(kDjango) ? get_hist_from_file(file,chains[kDjango],directory,allhistos.at(5)) : NULL;

   double acceptance, purity;
   int iPad=1;

   TH1D* HistoAcceptance = new TH1D(Form("Acceptance_tau_%0d",i),";Acceptance; bin number",5, 0, 5);
   TH1D* HistoPurity     = new TH1D(Form("Acceptance_tau_%0d",i),";Purity;bin number",    5, 0, 5);
   for (int i = 1; i<6; i++){      //1D taubins
      int binnumber = HistoGenEvents_1D->GetBin(i); //root starts counting bins at 1! 
      // for empty bins, set acceptance and purity to zero
      if(HistoGenEvents_1D->GetBinContent(binnumber) == 0 || HistoRecEvents_1D->GetBinContent(binnumber) == 0){
         acceptance = 0;
         purity = 0;
      }
      else{
         // acceptance = N_rec/N_gen
         acceptance = HistoRecEvents_1D->GetBinContent(binnumber)/HistoGenEvents_1D->GetBinContent(binnumber);
         // purity = N_stay/N_rec
         purity = HistoStayEvents_1D->GetBinContent(binnumber)/HistoRecEvents_1D->GetBinContent(binnumber);
      }
      HistoAcceptance_1D->SetBinContent(binnumber, acceptance);
      HistoPurity_1D->SetBinContent(binnumber, purity);
      //naming convention: first 2 digits -> q2 (=0 -> all bins), third digit -> x (=0 -> all bins), fourth digit -> tau_zQ (=0 -> all bins)
      HistoAcceptance_1D->SetTitle(Form("Acceptance_%04d",binnumber));
      HistoPurity_1D->SetTitle(Form("Purity_%04d",binnumber));
   }
   cc->cd(2*iPad-1);
   HistoAcceptance_1D->SetMarkerColor(kBlue);
   HistoAcceptance_1D->SetMarkerStyle(3);
   HistoAcceptance_1D->SetMarkerSize(1);
   HistoAcceptance_1D->Draw("P");

   cc->cd(2*iPad);
   HistoPurity_1D->SetMarkerColor(kBlue);
   HistoPurity_1D->SetMarkerStyle(3);
   HistoPurity_1D->SetMarkerSize(1);
   HistoPurity_1D->Draw("P");


   for (int i = 1; i<13; i++){      //q2bins
      int BinNumber = 0;
      // create histogams for acceptance and purity plots
      // one plot for each q2 bin -> 25 bins per histo
      TH1D* HistoAcceptance = new TH1D(Form("Acceptance_%0d",i),";Acceptance;global bin number",5*5, 0, 5*5);
      TH1D* HistoPurity     = new TH1D(Form("Acceptance_%0d",i),";Purity;global bin number",    5*5, 0, 5*5);
      for (int j = 1; j<6; j++){    //xbins 
         for (int k = 1; k<6; k++){ //taubins
            int binnumber = HistoGenEvents->GetBin(i,j,k); 
            BinNumber += 1; //bin numbers are the same in each histo
            // for empty bins, set acceptance and purity to zero
            if(HistoGenEvents->GetBinContent(binnumber) == 0 || HistoRecEvents->GetBinContent(binnumber) == 0){
               acceptance = 0;
               purity = 0;
            }
            else{
               // acceptance = N_rec/N_gen
               acceptance = HistoRecEvents->GetBinContent(binnumber)/HistoGenEvents->GetBinContent(binnumber);
               // purity = N_stay/N_rec
               purity = HistoStayEvents->GetBinContent(binnumber)/HistoRecEvents->GetBinContent(binnumber);
            }
            HistoAcceptance->SetBinContent(BinNumber, acceptance);
            HistoPurity->SetBinContent(BinNumber, purity);
            //naming convention: first 2 digits -> q2, third digit -> x (=0 -> all bins), fourth digit -> tau_zQ (=0 -> all bins)
            int histonumber = i*100;
            HistoAcceptance->SetTitle(Form("Acceptance_%04d",histonumber));
            HistoPurity->SetTitle(Form("Purity_%04d",histonumber));
         }
      }

      cc->cd(2*iPad-1);
      HistoAcceptance->SetMarkerColor(kBlue);
      HistoAcceptance->SetMarkerStyle(3);
      HistoAcceptance->SetMarkerSize(1);
      HistoAcceptance->Draw("P");

      cc->cd(2*iPad);
      HistoPurity->SetMarkerColor(kBlue);
      HistoPurity->SetMarkerStyle(3);
      HistoPurity->SetMarkerSize(1);
      HistoPurity->Draw("P");

      //clear canvas
      ++iPad;
      int temp = 2*iPad;
      if ( temp>4 ) {
         iPad=1;
         cc->Print(outps);
         cc->Clear();
         cc->Divide(2,2);
      }

   }


   

   // plot all histos
   for ( string histname : allhistos ) {
      TH1* hDjango = chains.count(kDjango) ? get_hist_from_file(file,chains[kDjango],directory,histname) :  NULL;
      
      TH1* hMain = NULL;
      if ( mainchain.find("Django") == 0 ) hMain = hDjango;
      string title = hMain->GetTitle();
      if ( title == "" ) title = hMain->GetName();
      title = "["+directory+"] "+ title;
      hMain->SetTitle(title.c_str());

      // set style
      SetStyle_Main(hMain);
      SetStyle_TH1(hDjango,kDjango);
      
      cc->cd(iPad);
      

      //1D histos
      if ( hMain->InheritsFrom("TH1D") || hMain->InheritsFrom("TH1F")  || hMain->InheritsFrom("TH1I")  ) {
         continue;
      }

      //3D histos
      if ( hMain->InheritsFrom("TH3D") || hMain->InheritsFrom("TH3F")  || hMain->InheritsFrom("TH3I")  ) {   
         if ( hDjango ) {
            hDjango->GetXaxis()->SetTitle("log_{10}(Q2)");
            hDjango->GetYaxis()->SetTitle("log_{10}(X)");
            hDjango->GetZaxis()->SetTitle("tau_{zQ}");
            hDjango->SetTitleOffset(1.8,"z");
            hDjango->SetTitleSize(0.06,"z");
            hDjango->SetLabelSize(0.06,"z");
            hDjango->SetLabelOffset(0.05,"z");
            
            hDjango->DrawClone("LEGO");
         }
      }
      gPad->RedrawAxis();

   

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
//! return list of all histograms in directory
void plot_directory_default(TFile* file, const string& directory, const TString& outps) {

   // set style and all that
   TCanvas* cc = new TCanvas("cc","plots",1200,1600);
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
      TH1* hCCBkgd = chains.count(kCCBkgd) ? get_hist_from_file(file,chains[kCCBkgd],directory,histname) :  NULL;
      TH1* hDjBkgd = chains.count(kDjBkgd) ? get_hist_from_file(file,chains[kDjBkgd],directory,histname) :  NULL;
      TH1* hBackground = chains.count(kPythia) ? get_hist_from_file(file,chains[kPythia],directory,histname) :NULL;// todo, get all background, and also add them to Django and Rapgap
/*
      if(hBackground && hCCBkgd){
         hBackground->Add(hCCBkgd);
      }
      if(hBackground && hDjBkgd){
         hBackground->Add(hDjBkgd);
      }
*/
      if(hBackground){
         hDjango->Add(hBackground);
         hRapgap->Add(hBackground);
      }

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

         //TCanvas *c1 = new TCanvas("c1","example",600,700);
         if ( hRapgap && hDjango ) {
            TPad *pad1 = new TPad("pad1","pad1",0.,0.35,1.,1.0);
            TPad *pad2 = new TPad("pad2","pad2",0.,0.,1.,0.35);
            if ( title.find("_lx")!=string::npos || string(hMain->GetName()).find("_lx") !=string::npos ) {
               pad1->SetLogx();
               pad2->SetLogx();
            }

            if ( title.find("_ly")!=string::npos || string(hMain->GetName()).find("_ly") !=string::npos ||
                 title.find("_lxy")!=string::npos || string(hMain->GetName()).find("_lxy") !=string::npos ) pad1->SetLogy();

            if ( pad1->GetLogy() )  hMain->SetMaximum(max*40);
            else                    hMain->SetMaximum(max*1.4);
            if ( pad1->GetLogy() ) hMain->SetMinimum(0.7);
            else                   hMain->SetMinimum(0.);
            
            pad2->SetTopMargin(0);
            pad2->SetBottomMargin(0.35);
            pad2->Draw();
            
            pad1->SetBottomMargin(0.);
            pad1->Draw();
            pad1->cd();

            hMain->DrawClone("hist");
//            gStyle->SetPadBottomMargin(0.5);
         //if ( hBackground ) hBackground->DrawClone("histsame");
         //if ( hRapgap ) hRapgap->DrawClone("histsame");
         //if ( hDjango ) hDjango->DrawClone("histsame");
            
            TH1* h3=hRapgap->DrawCopy("histsame");
            //hRapgap->Draw("histsame");
            TH1* h1=hDjango->DrawCopy("histsame");
            //hDjango->Draw("histsame");
            if ( hBackground ) hBackground->DrawClone("histsame");
            if ( hData       ) hData->DrawClone("PEsame");
  
            h1->SetStats(0);
            h3->SetStats(0);

            TLegend* l = new TLegend(0.65,0.70,0.9,0.93,"","brNDC");
            l->SetTextSize(0.04);
            if ( hData   )     l->AddEntry(hData,  "Data","PE");
            if ( hDjango )     l->AddEntry(hDjango,"Django","L");
            if ( hRapgap )     l->AddEntry(hRapgap,"Rapgap","L");
            if ( hBackground ) l->AddEntry(hBackground,"Backgrounds","L");
            l->SetFillStyle(0);
            l->SetBorderSize(1);
            l->DrawClone();

            pad1->RedrawAxis();
            
            //pad2->SetTopMargin(0);
            //pad2->SetBottomMargin(0.2);
            //pad2->Draw();
            pad2->cd();
            hMain->DrawClone("hist");
            //hMain->GetYaxis()->SetRangeUser(0.4, 1.6);
            TH1* h2=hRapgap->DrawCopy();
//            TH1* h4=hDjango->DrawCopy();
            if ( hData ) {
               h2->SetTitle("");
               h2->GetYaxis()->SetTitle("Ratio");
               h2->GetYaxis()->SetTickSize(0.025);
               h2->GetYaxis()->SetNdivisions(205);
               h2->GetYaxis()->SetRangeUser(0.4, 1.6);
               h2->GetXaxis()->SetTickSize(0.04);
               h2->SetTitleOffset(1.5,"x");
               h2->SetTitleOffset(.6,"y");
               h2->SetTitleSize(0.11,"x");
               h2->SetTitleSize(0.11,"y");
               h2->SetLabelSize(0.11,"x");
               h2->SetLabelSize(0.11,"y");
               h2->SetLabelOffset(0.001,"x");
               h2->SetLabelOffset(0.008,"y");

               hDjango->SetStats(0);
               hDjango->Divide(hData);
               hDjango->SetMarkerStyle(21);
               hDjango->SetMarkerSize(.5);
               hDjango->SetMarkerColor(kBlue);
               hDjango->SetLineWidth(1);
               hDjango->SetLineColor(kBlue);
               hDjango->Draw("epsame");

               h2->SetStats(0);
               h2->Divide(hData);
               h2->SetMarkerStyle(21);
               h2->SetMarkerSize(.5);
               h2->SetMarkerColor(kRed);
               h2->SetLineWidth(1);
               h2->SetLineColor(kRed);
               h2->Draw("epsame");
            }
            else {
               hRapgap->SetTitle("");
               hRapgap->GetYaxis()->SetTitle("Ratio");
               hRapgap->GetYaxis()->SetTickSize(0.025);
               hRapgap->GetYaxis()->SetNdivisions(205);
               hRapgap->GetYaxis()->SetRangeUser(0.4, 1.6);
               hRapgap->GetXaxis()->SetTickSize(0.04);
               hRapgap->SetTitleOffset(1.5,"x");
               hRapgap->SetTitleOffset(.6,"y");
               hRapgap->SetTitleSize(0.11,"x");
               hRapgap->SetTitleSize(0.11,"y");
               hRapgap->SetLabelSize(0.11,"x");
               hRapgap->SetLabelSize(0.11,"y");
               hRapgap->SetLabelOffset(0.001,"x");
               hRapgap->SetLabelOffset(0.008,"y");

               hRapgap->SetStats(0);
               hRapgap->Divide(hDjango);
               hRapgap->SetMarkerStyle(21);
               hRapgap->SetMarkerSize(.5);
               hRapgap->SetMarkerColor(kRed);
               hRapgap->SetLineWidth(1);
               hRapgap->SetLineColor(kRed);
               hRapgap->Draw("ep");

            }
            
            pad2->Update();
            
            TLine *line = new TLine(hRapgap->GetXaxis()->GetXmin(),1,hRapgap->GetXaxis()->GetXmax(),1);
            line->SetLineColor(kBlack);
            line->SetLineWidth(1);
            line->SetLineStyle(2);
            line->Draw("same");
            
            pad2->RedrawAxis();
         }
         
         // if ( histname.Contains("_den") ) {
         //    TString numname = histname.ReplaceAll("_den","_num");
         //    ...
         //       }
      }
   

      // plot a 2D-histogram
      else if ( hMain->InheritsFrom("TH2D") || hMain->InheritsFrom("TH2F")  || hMain->InheritsFrom("TH2I")  ) {

         gPad->SetLogz();

         double max = hMain->GetMaximum();
         if ( title.find("_lx")!=string::npos || string(hMain->GetName()).find("_lx") !=string::npos ) gPad->SetLogx();
         if ( title.find("_ly")!=string::npos || string(hMain->GetName()).find("_ly") !=string::npos ||
              title.find("_lxy")!=string::npos || string(hMain->GetName()).find("_lxy") !=string::npos ) gPad->SetLogy();
         if ( gPad->GetLogy() )  hMain->SetMaximum(max*40);
         else                    hMain->SetMaximum(max*1.4);
         if ( gPad->GetLogy() ) hMain->SetMinimum(0.7);
         else                   hMain->SetMinimum(0.);


         hMain->DrawClone();
         //if ( hBackground ) hBackground->DrawClone("histsame");
         //if ( hRapgap ) hRapgap->DrawClone("histsame");
         if ( hDjango ) hDjango->DrawClone("COLZ");
         //if ( hData   ) hData->DrawClone("PEsame");
         //if ( hBackground ) hBackground->DrawClone("histsame");

         gPad->RedrawAxis();
         gPad->SaveAs("test.pdf");


      }


      // plot a 3D-histogram
      else if ( hMain->InheritsFrom("TH3D") || hMain->InheritsFrom("TH3F")  || hMain->InheritsFrom("TH3I")  ) {

         //gPad->SetLogz();

         //double max = hMain->GetMaximum();
         if ( title.find("_lx")!=string::npos || string(hMain->GetName()).find("_lx") !=string::npos ) gPad->SetLogx();
         if ( title.find("_ly")!=string::npos || string(hMain->GetName()).find("_ly") !=string::npos ||
              title.find("_lxy")!=string::npos || string(hMain->GetName()).find("_lxy") !=string::npos ) gPad->SetLogy();
         if ( title.find("_lz")!=string::npos || string(hMain->GetName()).find("_lz") !=string::npos ||
              title.find("_lxz")!=string::npos || string(hMain->GetName()).find("_lxz") !=string::npos ||
              title.find("_lyz")!=string::npos || string(hMain->GetName()).find("_lyz") !=string::npos ||
              title.find("_lxyz")!=string::npos || string(hMain->GetName()).find("_lxyz") !=string::npos ) gPad->SetLogz();
/*
         if ( gPad->GetLogy() )  hMain->SetMaximum(max*40);
         else                    hMain->SetMaximum(max*1.4);
         if ( gPad->GetLogy() ) hMain->SetMinimum(0.7);
         else                   hMain->SetMinimum(0.);
*/

         //hMain->DrawClone();                                                                                                                                                     //if ( hBackground ) hBackground->DrawClone("histsame");
         //if ( hRapgap ) hRapgap->DrawClone("histsame");
         if ( hDjango ) {
            hDjango->GetXaxis()->SetTitle("log_{10}(Q2)");
            hDjango->GetYaxis()->SetTitle("log_{10}(X)");
            hDjango->GetZaxis()->SetTitle("tau_{zQ}");
            hDjango->SetTitleOffset(1.8,"z");
            hDjango->SetTitleSize(0.06,"z");
            hDjango->SetLabelSize(0.06,"z");
            hDjango->SetLabelOffset(0.05,"z");

            hDjango->DrawClone("LEGO");
         }
         //if ( hData   ) hData->DrawClone("PEsame");
         //if ( hBackground ) hBackground->DrawClone("histsame");

         gPad->RedrawAxis();
      
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
void plot_all_ratio(){
   
   string inputfile = "merged.root";
   //string inputfile = "RootOut/DataEplus0304.root";
   TFile* file = TFile::Open(inputfile.c_str(),"READ");

   // print chains in file
   get_chains(file,true);

   // loop over all directories and plot them
   vector<string> directories = get_directories(file);
   for ( string directory : directories ) {
      if ( directory == "minitree" ) continue;
      cout<<" -------------------------------------------------- "<<endl;
      cout<<" [main]  Plotting directory:  "<<directory<<endl;
      cout<<"         using:    plot_directory_default"<<endl;
      string outps = "Plots/"+directory+".ps";
      cout<<"         output:   "<<outps<<endl;
      plot_directory_default(file,directory,outps);
      cout<<" -------------------------------------------------- "<<endl;

      if ( directory == "Acceptance") {
         cout<<"Now analyzing directory Acceptance"<<endl;
         cout<<"Study acceptance and purity"<<endl;
         plot_directory_acceptance(file,directory,outps);
         cout<<" -------------------------------------------------- "<<endl;
      }
   }


   string outps = "plots.ps";

}

int main() { 
   plot_all_ratio();
   return 0;
}
