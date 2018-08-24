double findCuts (TH2D* [36][2], int);

void makeHisto() {
  int runNumber = 321396;
  TFile* input = TFile::Open(Form("histo_%i.root", runNumber), "READ");
  
  int iHit = 0;
  int iCSC = 0;
  int run_;

  TFile* output = TFile::Open(Form("output_%i.root", runNumber), "RECREATE");
  TTree* runTree = new TTree("runTree", "runTree");
  runTree->Branch("run", &run_, "run/I");

  // Load Geometry Setup
  float stripLength[2][8];
  float stripPitch[2][8];
  float cscArea;
  TTree* t_setup = (TTree*) input->Get("SliceTestAnalysis/Setup");
  t_setup->SetBranchAddress("stripLength", stripLength);
  t_setup->SetBranchAddress("stripPitch", stripPitch);
  t_setup->SetBranchAddress("cscArea", &cscArea);
  t_setup->GetEntry(0);
  float trArea[36][2][8];
  float chArea[36][2];
  float vfatArea[36][2][8][3];
  float totArea;
  
  // Prepare get info about active areas
  TH2D* hAllStrip[36][2];
  TH2D* hMaskedStrip[36][2];
  int nOfActiveStrip[36][2][8];
  int nOfActiveStripVfat[36][2][8][3];
  for (int iChamber = 26; iChamber < 30; iChamber++) {
    for (int iLayer = 0; iLayer < 2; iLayer++) {
      hAllStrip[iChamber][iLayer] = (TH2D*) input->Get(Form("SliceTestAnalysis/allStrips ch %i lay %i", iChamber+1, iLayer+1));
      
      hMaskedStrip[iChamber][iLayer] = new TH2D(Form("maskedStrip ch %i layer %i", iChamber+1, iLayer+1), Form("maskedStrip ch %i layer %i", iChamber+1, iLayer+1), 384, 1, 385, 8, 0.5, 8.5);
      hMaskedStrip[iChamber][iLayer]->GetXaxis()->SetTitle("Strip");
      hMaskedStrip[iChamber][iLayer]->GetYaxis()->SetTitle("iEta");

      for (int iEta = 1; iEta < 9; iEta++) {
        nOfActiveStrip[iChamber][iLayer][iEta-1] = 0;
        for (int iStrip = 1; iStrip < 385; iStrip++) {
          hMaskedStrip[iChamber][iLayer]->Fill(iStrip, iEta);
        }
        for (int iStrip = 0; iStrip < 3; iStrip++) {
          nOfActiveStripVfat[iChamber][iLayer][iEta-1][iStrip] =0;
        }
      }
    }
  }
  
  // Load Event
  int run, lumi, nGEMHits, nCSCHits;
  float instLumi;
  TTree* t_Event = (TTree*) input->Get("SliceTestAnalysis/Event");
  t_Event->SetBranchAddress("run", &run);
  t_Event->SetBranchAddress("lumi", &lumi);
  t_Event->SetBranchAddress("nGEMHits", &nGEMHits);
  t_Event->SetBranchAddress("nCSCHits", &nCSCHits);
  t_Event->SetBranchAddress("instLumi", &instLumi);
  
  int nEvent = t_Event->GetEntries();
  map<int, int> nRunEvent;
  map<int, map<int,int>> nLumiEvent;
  map<int, TH1D*> hInstLumi;
  map<int, TH1D*[36][2]> hHitRateLumi;
  map<int, TH2D*[36][2]> hAllStrips;
  map<int, TH2D*[36][2]> hHitRateVfat;
  map<int, TH1D*[36][2]> hHitRateCSC;
  for (int i = 0; i < nEvent; i++) {
    t_Event->GetEntry(i);
    nRunEvent[run]++;
    nLumiEvent[run][lumi]++;
  }
  for (auto iter = nRunEvent.begin(); iter != nRunEvent.end(); iter++) {
    run_ = iter->first;
    hInstLumi[run_] = new TH1D(Form("instLumi %i", iter->first), "instLumi", 50, 0, 2);
    hInstLumi[run_]->GetXaxis()->SetTitle("instantaneous luminosity[cm^{-2}s^{-1}]");
    hInstLumi[run_]->GetYaxis()->SetTitle("Event");
    
    hHitRateLumi[run_][0][0] = new TH1D(Form("%i hitRate", run_), Form("Run number %i", run_), 50, 0, 2);
    hHitRateLumi[run_][0][0]->GetXaxis()->SetTitle("instantaneous luminosity[cm^{-2}s^{-1}]");
    hHitRateLumi[run_][0][0]->GetYaxis()->SetTitle("Hit rate[Hz/cm^{2}]");

    for (int iChamber = 26; iChamber < 30; iChamber++) {
      for (int iLayer = 0; iLayer < 2; iLayer++) {
        hHitRateLumi[run_][iChamber][iLayer] = new TH1D(Form("%i hitRate ch %i layer %i", run_, iChamber+1, iLayer+1), Form("Run number %i", run_), 50, 0, 2);
        hHitRateLumi[run_][iChamber][iLayer]->GetXaxis()->SetTitle("instantaneous luminosity[cm^{-2}s^{-1}]");
        hHitRateLumi[run_][iChamber][iLayer]->GetYaxis()->SetTitle("Hit rate[Hz/cm^{2}]");
        hHitRateLumi[run_][iChamber][iLayer]->SetMarkerSize(2);

        //hHitRateCSC[run_][iChamber][iLayer] = new TH1D(Form("%i cscRate ch %i layer %i", run_, iChamber+1, iLayer+1), Form("Run number %i", run_), 50, 0, 2);
        //hHitRateCSC[run_][iChamber][iLayer]->GetXaxis()->SetTitle("instantaneous luminosity[cm^{-2}s^{-1}]");
        //hHitRateCSC[run_][iChamber][iLayer]->GetYaxis()->SetTitle("Hit rate[Hz/cm^{2}]");
        //hHitRateCSC[run_][iChamber][iLayer]->SetMarkerSize(2);
        
        hHitRateVfat[run_][iChamber][iLayer] = new TH2D(Form("%i vfat ch %i layer %i", run_, iChamber+1, iLayer+1), Form("ch %i layer %i", iChamber+1, iLayer+1), 3, 0.5, 384.5, 8, 0.5, 8.5);
        hHitRateVfat[run_][iChamber][iLayer]->GetXaxis()->SetTitle("strip");
        hHitRateVfat[run_][iChamber][iLayer]->GetYaxis()->SetTitle("iEta");

        hAllStrips[run_][iChamber][iLayer] = new TH2D(Form("%i allStrips ch %i layer %i", run_, iChamber+1, iLayer+1), Form("ch %i layer %i", iChamber+1, iLayer+1), 384, 1, 385, 8, 0.5, 8.5);
        hAllStrips[run_][iChamber][iLayer]->GetXaxis()->SetTitle("strip");
        hAllStrips[run_][iChamber][iLayer]->GetYaxis()->SetTitle("iEta");
      }
    }
    runTree->Fill();
  }
  for (int i = 0; i < nEvent; i++) {
    t_Event->GetEntry(i);
    hInstLumi[run]->Fill(instLumi);
  }

  // Load RecHits
  int chamber, layer, etaPartition, firstStrip, nStrips;
  TTree* t_Hit = (TTree*) input->Get("SliceTestAnalysis/Hit");
  t_Hit->SetBranchAddress("chamber", &chamber);
  t_Hit->SetBranchAddress("layer", &layer);
  t_Hit->SetBranchAddress("etaPartition", &etaPartition);
  t_Hit->SetBranchAddress("firstStrip", &firstStrip);
  t_Hit->SetBranchAddress("nStrips", &nStrips);

  //int chamber_, layer_, nStrips_;
  //TTree* t_csc = (TTree*) input->Get("SliceTestAnalysis/csc");
  //t_csc->SetBranchAddress("chamber", &chamber_);
  //t_csc->SetBranchAddress("layer", &layer_);
  //t_csc->SetBranchAddress("nStrips", &nStrips_);

  double cut = findCuts(hAllStrip, nEvent);

  for (int iChamber = 26; iChamber < 30; iChamber++) {
    for (int iLayer = 0; iLayer < 2; iLayer++) {
      if (hAllStrip[iChamber][iLayer]->GetEntries() == 0) continue;
      auto xAxis = hAllStrip[iChamber][iLayer]->GetXaxis();
      auto yAxis = hAllStrip[iChamber][iLayer]->GetYaxis();
      for (int iEta = 1; iEta < 9; iEta++) {
        auto yBin = yAxis->FindBin(iEta);
        for (int iStrip = 1; iStrip < 385; iStrip++) {
          auto xBin = xAxis->FindBin(iStrip);
          double val = hAllStrip[iChamber][iLayer]->GetBinContent(xBin, yBin);
          if (val < 1) continue;
          val = double(val/nEvent);
          if (val > cut) continue;
          hMaskedStrip[iChamber][iLayer]->SetBinContent(xBin, yBin, 0);
          nOfActiveStrip[iChamber][iLayer][iEta-1]++;
          nOfActiveStripVfat[iChamber][iLayer][iEta-1][iStrip/128]++;
        }
      }
    }
  }

  totArea = 0;
  for ( int ch = 0; ch < 36; ch++ ) {
    for ( int ly = 0; ly < 2; ly++) {
      chArea[ch][ly] = 0;
      for ( int et = 0; et < 8; et++) {
        trArea[ch][ly][et] = stripLength[(ch+1)%2][et] * stripPitch[(ch+1)%2][et] * (nOfActiveStrip[ch][ly][et]);
        chArea[ch][ly] += trArea[ch][ly][et];
        for ( int ix = 0; ix < 3; ix++) vfatArea[ch][ly][et][ix] = stripLength[(ch+1)%2][et] * stripPitch[(ch+1)%2][et] * (nOfActiveStripVfat[ch][ly][et][ix]);
      }
      totArea += chArea[ch][ly];
    }
  }

  bool flag;
  for (int i = 0; i < nEvent; i++) {
    t_Event->GetEntry(i);
    flag = false;
    int nEvent_ = hInstLumi[run]->GetBinContent(hInstLumi[run]->GetXaxis()->FindBin(instLumi));
    for (int j = 0; j < nGEMHits; j++) {
      t_Hit->GetEntry(iHit);
      iHit++;

      auto etaBin = hAllStrip[chamber-1][layer-1]->GetYaxis()->FindBin(etaPartition);
      
      double normFactor = (chArea[chamber-1][layer-1]*nEvent_*25E-9);

      for (int k = 0; k < nStrips; k++) {
        hAllStrips[run][chamber-1][layer-1]->Fill(firstStrip+k, etaPartition);
        auto xBin = hAllStrip[chamber-1][layer-1]->GetXaxis()->FindBin(firstStrip+k);
        if(hAllStrip[chamber-1][layer-1]->GetBinContent(xBin, etaBin)/nEvent > cut) continue;
        hHitRateLumi[run][chamber-1][layer-1]->Fill(instLumi, 1./(normFactor));
        double normFactorVfat = (vfatArea[chamber-1][layer-1][etaPartition-1][(firstStrip+k)/128]*nEvent_*25E-9);
        hHitRateVfat[run][chamber-1][layer-1]->Fill(firstStrip+k, etaPartition, 1./(normFactorVfat));
        flag = true;
      }
    }
    if (flag) hHitRateLumi[run][0][0]->Fill(instLumi, 1./(totArea*nEvent_*25E-9));
    //for (int j = 0; j < nGEMHits; j++) {
    //  t_csc->GetEntry(iCSC);
    //  iCSC++;

    //  double normFactor = (cscArea*nEvent_*25E-9);

    //  for (int k = 0; k < nStrips_; k++) {
    //    hHitRateCSC[run][chamber_-1][layer_-1]->Fill(instLumi, 1./(normFactor));
    //    flag = true;
    //  }
    //}
  }

  output->Write();
  output->Close();
}


double findCuts(TH2D* hAllStrip[36][2], int nEvent) {
  TH1D* hitHist = new TH1D( "hitHist", "hits per each strip", 500, 0, 4e-4);
  hitHist->GetXaxis()->SetTitle("number of hits per each strip");
  hitHist->GetYaxis()->SetTitle("strips");
  for (int iChamber = 26; iChamber < 30; iChamber++) {
    for (int iLayer = 0; iLayer < 2; iLayer++) {
      if (hAllStrip[iChamber][iLayer]->GetEntries() == 0) continue;
      auto xAxis = hAllStrip[iChamber][iLayer]->GetXaxis();
      auto yAxis = hAllStrip[iChamber][iLayer]->GetYaxis();
      for (int iEta = 1; iEta < 9; iEta++) {
        auto yBin = yAxis->FindBin(iEta);
        for (int iStrip = 1; iStrip < 385; iStrip++) {
          auto xBin = xAxis->FindBin(iStrip);
          double val = hAllStrip[iChamber][iLayer]->GetBinContent(xBin, yBin);
          val = double(val/nEvent);
          if (hAllStrip[iChamber][iLayer]->GetBinContent(xBin, yBin) == 0) continue;
          else if(val >= 4E-4) val = 3.99E-4;
          hitHist->Fill((double)val);
        }
      }
    }
  }
  hitHist->Write();
  double res = (hitHist->GetMean() + 2*hitHist->GetRMS());
  delete hitHist;
  return res; 
}
