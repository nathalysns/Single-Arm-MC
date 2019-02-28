{
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Mon Dec  5 13:15:05 2016 by ROOT version5.34/13)
//   from TTree h1/HUT
//   found on file: test.root
//////////////////////////////////////////////////////////


//Reset ROOT and connect tree file
   gROOT->Reset();
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("test.root");
   if (!f) {
      f = new TFile("test.root");
   }
    f->GetObject("h1",tree);

//Declaration of leaves types
   Float_t         hsxfp;
   Float_t         hsyfp;
   Float_t         hsxpfp;
   Float_t         hsypfp;
   Float_t         hsytari;
   Float_t         hsdeltai;
   Float_t         hsyptari;
   Float_t         hsxptari;
   Float_t         hsytar;
   Float_t         hsdelta;
   Float_t         hsyptar;
   Float_t         hsxptar;
   Float_t         fry;
   Float_t         frx;
   Float_t         ok_spec;
   Float_t         stopwhen;
   Float_t         x_stop;
   Float_t         y_stop;

   // Set branch addresses.
   h1->SetBranchAddress("hsxfp",&hsxfp);
   h1->SetBranchAddress("hsyfp",&hsyfp);
   h1->SetBranchAddress("hsxpfp",&hsxpfp);
   h1->SetBranchAddress("hsypfp",&hsypfp);
   h1->SetBranchAddress("hsytari",&hsytari);
   h1->SetBranchAddress("hsdeltai",&hsdeltai);
   h1->SetBranchAddress("hsyptari",&hsyptari);
   h1->SetBranchAddress("hsxptari",&hsxptari);
   h1->SetBranchAddress("hsytar",&hsytar);
   h1->SetBranchAddress("hsdelta",&hsdelta);
   h1->SetBranchAddress("hsyptar",&hsyptar);
   h1->SetBranchAddress("hsxptar",&hsxptar);
   h1->SetBranchAddress("fry",&fry);
   h1->SetBranchAddress("frx",&frx);
   h1->SetBranchAddress("ok_spec",&ok_spec);
   h1->SetBranchAddress("stopwhen",&stopwhen);
   h1->SetBranchAddress("x_stop",&x_stop);
   h1->SetBranchAddress("y_stop",&y_stop);

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// h1->SetBranchStatus("*",0);  // disable all branches
// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

   Long64_t nentries = h1->GetEntries();

   Long64_t nbytes = 0;
//   for (Long64_t i=0; i<nentries;i++) {
//      nbytes += h1->GetEntry(i);
//   }
}
