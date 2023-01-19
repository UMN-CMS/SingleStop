void makeDRPlots(){

  using namespace std;

  string inputDir = "Run3/22-12-12";

  gStyle->SetPalette(kBird);//kPastel);
  gStyle->SetOptTitle(0);
  gStyle->SetNumberContours(256);
  vector<string> compare = {"200_100","300_100","300_200","500_100","500_200","500_400","700_100","700_400","700_600","1000_100","1000_400","1000_900","1500_100","1500_600","1500_1400","2000_100","2000_900","2000_1900","700_200","1000_200","1500_200","1500_400","2000_200","2000_400"};
  vector<string> plotNames = {"passDijet"};//"passDijet"};//"dRmaxchi"}; //"HT","m4","m3NoLead","nJets","nb","pT1","pT2","jetMatch","leadJetMatch","nJetsChiMerged"};
  for(string n:plotNames){
      TCanvas* c1 = new TCanvas();
      c1->SetLeftMargin(0.15);
      c1->SetRightMargin(0.15);
      //c1->SetGrid();
      TH2F* chistop = new TH2F("chistop",Form(";m_{#tilde{t}} [GeV]; m_{#chi^{#pm}} [GeV]; %s","Fraction of events failing dijet search jet reqs."), 20, 50, 2050 , 20, 50, 2050);
      string xTitle, yTitle;
      for(string s:compare){
       
        TFile* file = new TFile(Form("output/%s/signal_%s.root",inputDir.c_str(),s.c_str()));
        TDirectory* plots = file->GetDirectory("plots");
        TH1F* hist = (TH1F*)plots->Get(n.c_str());
 
        ////Find the lower bin we want
        //int b1 = hist->GetXaxis()->FindBin(1);
        ////Find the max bin in the plots (plus the overflow bin)
        //int b2 = hist->GetXaxis()->GetNbins() + 1;
        //// Normalize the histogram
        //hist->Scale(1./hist->GetEntries());
        ////Calculate the integral (fraction of events above b1)
        //double frac = hist->Integral(b1,b2);
      
        int mass_stop = stoi(s.substr(0, s.find("_")));      
        int mass_chargino = stoi(s.substr(s.find("_") + 1));
 
        stringstream plotTitle;
        gStyle->SetPaintTextFormat(".2f");
        chistop -> Fill(mass_stop, mass_chargino, 1 - hist->GetMean());
        plotTitle << "(" << s.substr(0,s.find("_")) << ", " << s.substr(s.find("_") + 1) << ")";
        gStyle->SetOptStat(0);
        chistop->SetTitle(plotTitle.str().c_str());
        xTitle = chistop->GetXaxis()->GetTitle();
        yTitle = chistop->GetYaxis()->GetTitle();
        }

      //chistop->GetXaxis()->SetRangeUser(50,2050);
      //chistop->GetYaxis()->SetRangeUser(50,2050); 
      //chistop->GetYaxis()->SetLimits(50,2000);
      //chistop->GetXaxis()->SetLimits(150,2050);
      chistop->GetXaxis()->SetTickLength(-0.02);
      chistop->GetXaxis()->SetLabelOffset(0.02);
      chistop->GetXaxis()->SetTitleOffset(1.3);
      chistop->GetYaxis()->SetTickLength(-0.02);      
      chistop->GetYaxis()->SetLabelOffset(0.025);
      chistop->GetYaxis()->SetTitleOffset(1.55);
      //chistop->GetYaxis()->SetNdivisions(20);
      //chistop->GetXaxis()->SetNdivisions(20);
      chistop->SetMarkerSize(1.3);
      
      chistop->Draw("COLZ TEXT");

      c1->SaveAs(Form("plots/dRPlots/dRWideJets22-12-12.png"));
  }   
}
