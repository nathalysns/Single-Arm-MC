{

 TRandom3 r;
   TH1F h1("h1","test1",100,-6,6);
   for (int i=0;i<100000;i++) h1.Fill(r.Gaus(0,1)/pow(2,0.5));
   h1.Draw();
   h1->Fit("gaus");
   new TCanvas;
   TH1F h2("h2","test2",100,-6,6);
   for (int i=0;i<100000;i++) h2.Fill(r.Gaus(0,1/pow(2,0.5)));
   h2.Draw("");
   h2->Fit("gaus");
}
