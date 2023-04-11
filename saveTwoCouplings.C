void saveTwoCouplings(){

  int doPlots = 0;
  int doCompare = 1;
  int doTrigPlots = 0;
  int doTrigMassPlots = 0;

  string inputDirs[]= ["313", "312"];
  string saveTo = "313_312";

  if (doPlots){
    vector<string> masses = {"2000_100"};//{"200_100","300_100","300_200","500_100","500_200","500_400","700_100","700_400","700_600","1000_100","1000_400","1000_900","1500_100","1500_600","1500_1400","2000_100","2000_900","2000_1900"};
    for (string inputDir:inputDirs) {
	    for(string s:masses){
		    TFile* file = new TFile(Form("output/%s/hist_%s.root",inputDir.c_str(),s.c_str()));
		    TDirectory* plots = file->GetDirectory("plots");
		    TIter nextkey(plots->GetListOfKeys());
		    TKey* key;
		    while ((key = (TKey*)nextkey())){
			    TH1F* h = (TH1F*)key->ReadObj();
			    TCanvas* c1 = new TCanvas();
			    h->Draw();
			    c1->SaveAs(Form("plots/%s/%s/%s.png",saveTo.c_str(),s.c_str(),h->GetName()));
		    }
	    }
    }
  }


  gStyle->SetPalette(kBlueRedYellow);//kPastel);
  gStyle->SetOptTitle(0);
  if (doCompare){
    vector<string> compare = {"signal_1000_400","signal_1000_900","signal_1500_900","signal_2000_400","signal_2000_1900"};
    vector<string> plotNames = {"nbLoose","nbMedium","nbTight","HT","nJets","mAll","m4","m3","pT1","pT2","pT3","pT4","eta1","eta2","eta3","eta4","m3NoLead","m3NoLeadOrSub","dEta12","dPhi12","dR12","ntLoose","ntMedium","ntTight","MET","dPhiMET1","dPhiMET2","dPhiMET3","dPhiMET4","HTLHE","nQLHE","nGLHE","nJLHE","pT1Gen","pT2Gen","pT3GEn","pT4Gen","eta1Gen","eta2Gen","eta3Gen","eta4Gen","pTStop","pTStopPlus","pTStopMinus","pTChi","pTChiPlus","pTChiMinus","pTBStop","pTBStopPlus","pTBStopMinus","pTBChi1","pTBChi2","pTBChiPlus","pTBChiMinus","etaStop","etaStopPlus","etaStopMinus","etaChi","etaChiPlus","etaChiMinus","etaBStop","etaBStopPlus","etaBStopMinus","etaBChi","etaBChiPlus","etaBChiMinus","dEtaBChi","dPhiBChi","nJetsChiMerged","dRBB","dEtaBB","dPhiBB","passDijet","dEtaWJs","mWJs","dRChiMax","dRBChi"};
    //float yMax = 0;
    for(string n:plotNames){
	    TCanvas* c1 = new TCanvas();
	    THStack* stack = new THStack();
	    TLegend* legend = new TLegend(0.68,0.89,0.89,0.65);
	    string xTitle, yTitle;
	    for(string s:compare){
		for (string inputDir:inputDirs) {
		    TFile* file = new TFile(Form("output/%s/%s.root",inputDir.c_str(),s.c_str()));
		    //TDirectory* plots = file->GetDirectory("plots");
		    TH1D* plot = (TH1D*)file->Get(n.c_str());
		    //TH1D* plot = (TH1D*)plots->Get(n.c_str());
		    stringstream plotTitle;
		    //plotTitle << "(" << s.substr(0,s.find("_")) << ", " << s.substr(s.find("_") + 1) << ")";
		    plotTitle << s;
		    plot->SetTitle(plotTitle.str().c_str());
		    plot->SetLineWidth(3);
		    plot->SetMarkerSize(0);
		    plot->Scale(1. / plot->Integral());
		    legend->AddEntry(plot);
		    stack->Add(plot);
		    xTitle = plot->GetXaxis()->GetTitle();
		    yTitle = plot->GetYaxis()->GetTitle();
		    //plot->Draw("HIST PLC PMC SAME");
		    //if (plot->GetMaximum() > yMax) yMax = plot->GetMaximum();
	    	}
	    } 
	    stack->Draw("HIST PLC NOSTACK");
	    stack->GetXaxis()->SetTitle(xTitle.c_str());
	    //stack->GetYaxis()->SetTitle("Fraction of Events");
	    stack->GetYaxis()->SetTitle("Events");
	    //legend->SetHeader("( m_{#tilde{t}} , m_{#tilde{#chi_{1}^{#pm}}} )","C");
	    //legend->SetTextSize(0.03);
	    legend->Draw();
	    //TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
	    //header->SetTextSize(0.05);
	    gPad->SetLogy();
	    c1->SaveAs(Form("plots/%s/overlaid/%s.pdf",saveTo.c_str(),n.c_str()));
    }
  }

  if (doTrigPlots){
    vector<string> compare = {"500_100","700_100","1000_100","1000_400","1500_100","1500_1400","2000_100","2000_1900"}; //"200_100","300_100","300_200","500_100","500_200","500_400","700_100","700_400","700_600","1000_100","1000_400","1000_900","1500_100","1500_600","1500_1400","2000_100","2000_900","2000_1900"}; //"200_100","2000_100","2000_1900","1000_400"};
    //vector<string> plotNames = {"HLT_AK8PFHT750_TrimMass50","HLT_AK8PFJet360_TrimMass30","HLT_AK8PFJet500","HLT_CaloJet500_NoJetID","HLT_PFHT1050","HLT_PFJet500","HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71","HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94","HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1"};
    vector<string> plotNames = {"HLT_PFHT1050","HLT_AK8PFJet360_TrimMass30","trigOR","HLT_PFHT430"};
    TCanvas* c1 = new TCanvas();
    //TH1F* trigEff = new TH1F("trigEff",";Trigger;Efficiency",9,0,9);
    THStack* stack = new THStack("trigEff",";;Efficiency");
    TLegend* legend = new TLegend(0.7,0.85,0.85,0.45);
    for (string s:compare){
	    TH1F *plotAll = new TH1F("","",9,0,9);
	    plotAll->SetTitle(s.c_str());
	    for (int i = 0; i < plotNames.size(); i++){
		    for (string inputDir:inputDirs) {
			    string n = plotNames[i];
			    TFile* file = new TFile(Form("output/%s/hist_%s.root",inputDir.c_str(),s.c_str()));
			    TDirectory* plots = file->GetDirectory("plots");
			    TH1F* plot = (TH1F*)plots->Get(n.c_str());
			    plotAll->SetBinContent(i+1,plot->GetMean());
			    plotAll->SetLineWidth(2);
		    }
	    }
	    legend->AddEntry(plotAll);
	    stack->Add(plotAll);
    }
  stack->Draw("HIST PLC NOSTACK");
  for(int i = 0; i < plotNames.size(); i++){
	stack->GetXaxis()->SetBinLabel(i+1,plotNames[i].erase(0,4).c_str());
  }  
  stack->GetXaxis()->SetLabelSize(0.028);
  stack->GetXaxis()->LabelsOption("u");
  //gPad->SetLogy();
  legend->Draw();
  gPad->SetBottomMargin(0.3);
  c1->SaveAs(Form("plots/%s/overlaid/%s.png",saveTo.c_str(),"trigEff"));
}

  gStyle->SetPalette(kBlueRedYellow);
  if (doTrigMassPlots){
    int nTot = 10000;
    vector<string> compare = {"200_100","300_100","500_100","700_100","1000_100","1500_100","2000_100"};
    vector<double> mstop = {200,300,500,700,1000,1500,2000};
    vector<string> plotNames = {"HLT_PFHT430","trigOR","HLT_AK8PFJet360_TrimMass30","HLT_PFHT1050"};
    TCanvas* c1 = new TCanvas();
    TMultiGraph* mg = new TMultiGraph();
    TLegend* legend = new TLegend(0.55,0.60,0.85,0.35);
    for (int i = 0; i < plotNames.size(); i++){
	    vector<double> eff, effError;
	    for (string s:compare){
		for (string inputDir:inputDirs) {
		    string n = plotNames[i];
		    TFile* file = new TFile(Form("output/%s/hist_%s.root",inputDir.c_str(),s.c_str()));
		    TDirectory* plots = file->GetDirectory("plots");
		    TH1F* plot = (TH1F*)plots->Get(n.c_str());
		    eff.push_back(plot->GetMean());
		    effError.push_back(sqrt(plot->GetMean()) / nTot);
		}
	    }
	    TGraphErrors* gre = new TGraphErrors(eff.size(),&mstop[0],&eff[0],0,0);//&effError[0]);
	    gre->SetTitle(plotNames[i].c_str());
	    //gre->SetMarkerColor();
	    gre->SetMarkerStyle(20);
	    gre->SetMarkerSize(0.8);
	    mg->Add(gre);
	    legend->AddEntry(gre);
    }
    mg->Draw("ALP PMC PLC");
    mg->GetXaxis()->SetTitle("m_{#tilde{t}} [GeV]");
    mg->GetYaxis()->SetTitle("Trigger Efficiency");
    mg->SetMinimum();
    //for(int i = 0; i < plotNames.size(); i++){
    //  stack->GetXaxis()->SetBinLabel(i+1,plotNames[i].erase(0,4).c_str());
    //}
    //stack->GetXaxis()->SetLabelSize(0.028);
    //stack->GetXaxis()->LabelsOption("u");
    legend->Draw();
    //gPad->SetBottomMargin(0.3);
    c1->SaveAs(Form("plots/%s/overlaid/%s.png",saveTo.c_str(),"trigEffMass"));
    gPad->SetLogy();
    c1->SaveAs(Form("plots/%s/overlaid/%s.png",saveTo.c_str(),"trigEffMassLog"));
    }
}
