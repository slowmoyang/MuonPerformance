array<double,8> findCuts(TH2D*[36][2], TH1D*[36][2][24]);
array<double,8> findCuts(TH2D*[36][2], TH1D*[36][2]);

void rebinHisto(TH1D* h);

void fillHisto() {
  int lumiFill = 7058;
  TFile* input = TFile::Open(Form("data/histo_%i.root", lumiFill), "READ");
  
  int iHit = 0;
  int iCSC = 0;
  int run_;

  TFile* output = TFile::Open(Form("output_%i.root", lumiFill), "RECREATE");
  TTree* runTree = new TTree("runTree", "runTree");

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
  float vfatArea[36][2][24];
  float totArea;
  
  // Load RecHits
  int chamber, layer, etaPartition, firstStrip, nStrips;
  TTree* t_Hit = (TTree*) input->Get("SliceTestAnalysis/Hit");
  t_Hit->SetBranchAddress("chamber", &chamber);
  t_Hit->SetBranchAddress("layer", &layer);
  t_Hit->SetBranchAddress("etaPartition", &etaPartition);
  t_Hit->SetBranchAddress("firstStrip", &firstStrip);
  t_Hit->SetBranchAddress("nStrips", &nStrips);

  // Prepare get info about active areas
  TH2D* hAllStrip[36][2];
  TH2D* hMaskedStrip[36][2];
  int nOfActiveStrip[36][2][8];
  int nOfActiveStripVfat[36][2][24];
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
      }
      for (int vfat = 0; vfat < 24; vfat++) {
        nOfActiveStripVfat[iChamber][iLayer][vfat] =0;
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
  map<int, map<int, int>> nLumiEvent;
  map<int, map<int, TH1D*>> hLumiStatus;
  map<int, TH2D*> hActiveLumi;
  TH1D* hInstLumi   [36][2][24];
  TH1D* hHitRateEta [36][2][8];
  TH1D* hHitRateLumi[36][2];
  TH2D* hAllStrips  [36][2];
  TH2D* hHitRateVfat[36][2];
  TH1D* hHitRateCSC [36][2];

  for (int iChamber = 26; iChamber < 30; iChamber++) {
    for (int iLayer = 0; iLayer < 2; iLayer++) {
      for (int vfat = 0; vfat < 24; vfat++) {
        hInstLumi[iChamber][iLayer][vfat] = new TH1D(Form("instLumi ch %i layer %i vfat %i", iChamber+1, iLayer+1, vfat), "instLumi", 50, 0, 2);
        hInstLumi[iChamber][iLayer][vfat]->GetXaxis()->SetTitle("instantaneous luminosity[cm^{-2}s^{-1}]");
        hInstLumi[iChamber][iLayer][vfat]->GetYaxis()->SetTitle("Event");
      }
    
      hHitRateLumi[iChamber][iLayer] = new TH1D(Form("hitRate ch %i layer %i", iChamber+1, iLayer+1), Form("Run number"), 50, 0, 2);
      hHitRateLumi[iChamber][iLayer]->GetXaxis()->SetTitle("instantaneous luminosity[cm^{-2}s^{-1}]");
      hHitRateLumi[iChamber][iLayer]->GetYaxis()->SetTitle("Hit rate[Hz/cm^{2}]");
      hHitRateLumi[iChamber][iLayer]->SetMarkerSize(2);

      for (int iEta = 0; iEta < 8; iEta++) {
        hHitRateEta[iChamber][iLayer][iEta] = new TH1D(Form("hitRate ch %i layer %i eta %i", iChamber+1, iLayer+1, iEta+1), Form("Run number"), 50, 0, 2);
        hHitRateEta[iChamber][iLayer][iEta]->GetXaxis()->SetTitle("instantaneous luminosity[cm^{-2}s^{-1}]");
        hHitRateEta[iChamber][iLayer][iEta]->GetYaxis()->SetTitle("Hit rate[Hz/cm^{2}]");
        hHitRateEta[iChamber][iLayer][iEta]->SetMarkerSize(2);
      }

      hHitRateVfat[iChamber][iLayer] = new TH2D(Form("vfat ch %i layer %i", iChamber+1, iLayer+1), Form("ch %i layer %i", iChamber+1, iLayer+1), 3, 0.5, 384.5, 8, 0.5, 8.5);
      hHitRateVfat[iChamber][iLayer]->GetXaxis()->SetTitle("strip");
      hHitRateVfat[iChamber][iLayer]->GetYaxis()->SetTitle("iEta");

      hAllStrips[iChamber][iLayer] = new TH2D(Form("allStrips ch %i layer %i", iChamber+1, iLayer+1), Form("ch %i layer %i", iChamber+1, iLayer+1), 384, 1, 385, 8, 0.5, 8.5);
      hAllStrips[iChamber][iLayer]->GetXaxis()->SetTitle("strip");
      hAllStrips[iChamber][iLayer]->GetYaxis()->SetTitle("iEta");
    }
  }
  
  for (int i = 0; i < nEvent; i++) {
    t_Event->GetEntry(i);
    nRunEvent[run]++;
    nLumiEvent[run][lumi]++;
  }
  for (auto iter = nRunEvent.begin(); iter != nRunEvent.end(); iter++) {
    run_ = iter->first;
    hActiveLumi[run_] = new TH2D(Form("%i activeLumi", run_), Form("Run number %i", run_), 1000, 0, 4000, 400, 27, 31);
    hActiveLumi[run_]->GetXaxis()->SetTitle("luminosityBlock");
    hActiveLumi[run_]->GetYaxis()->SetTitle("chamber + layer/2 + vfat/48");
  }
  for (auto iter = nLumiEvent.begin(); iter != nLumiEvent.end(); iter++) {
    run_ = iter->first;
    auto lumiMap = iter->second;
    for (auto iter_ = lumiMap.begin(); iter_ != lumiMap.end(); iter_++) {
      auto lumi_ = iter_->first;
      hLumiStatus[run_][lumi_] = new TH1D(Form("%i %i lumiStatus", run_, lumi_), Form("Run number %i lumi %i", run_, lumi_), 1000, 0, 2);
      hLumiStatus[run_][lumi_]->GetXaxis()->SetTitle("Instantaneous Luminosity [10^{34}cm^{-2}s^{-1}]");
    }
  }

  iHit = 0;
  for (int i = 0; i < nEvent; i++) {
    t_Event->GetEntry(i);
    hLumiStatus[run][lumi]->Fill(instLumi);
    for (int j = 0; j < nGEMHits; j++) {
      t_Hit->GetEntry(iHit);
      auto vfat = (8-etaPartition)+8*((firstStrip-1)/128); 
      iHit++;
      hActiveLumi[run]->Fill(lumi,(chamber+(layer-1)*0.5+vfat/48.));
    }
  }

  for (int i = 0; i < nEvent; i++) {
    t_Event->GetEntry(i);
    if (hLumiStatus[run][lumi]->GetRMS() > 0.01) continue;
    auto lumiBin = hActiveLumi[run]->GetXaxis()->FindBin(lumi);
    for (int iChamber = 27; iChamber < 31; iChamber++) {
      for (int iLayer = 1; iLayer < 3; iLayer++) {
        for (int vfat = 0; vfat < 24; vfat++) {
          //auto fromBin = hActiveLumi[run]->GetYaxis()->FindBin(iChamber+(iLayer-1)*0.5);
          //auto toBin = hActiveLumi[run]->GetYaxis()->FindBin(iChamber+(iLayer-1)*0.5+23/48.);
          auto chanBin = hActiveLumi[run]->GetYaxis()->FindBin(iChamber+(iLayer-1)*0.5+vfat/48.);
          if (hActiveLumi[run]->Integral(lumiBin, lumiBin, chanBin, chanBin) != 0) hInstLumi[iChamber-1][iLayer-1][vfat]->Fill(instLumi);
        }
      }
    }
  }
  
  cout << "test\n";
  auto cut = findCuts(hAllStrip, hInstLumi);

  cout << "test\n";
  for (int iChamber = 26; iChamber < 30; iChamber++) {
    for (int iLayer = 0; iLayer < 2; iLayer++) {
      if (hAllStrip[iChamber][iLayer]->GetEntries() == 0) continue;
      auto xAxis = hAllStrip[iChamber][iLayer]->GetXaxis();
      auto yAxis = hAllStrip[iChamber][iLayer]->GetYaxis();
      for (int iEta = 1; iEta < 9; iEta++) {
        auto yBin = yAxis->FindBin(iEta);
        for (int iStrip = 1; iStrip < 385; iStrip++) {
          auto xBin = xAxis->FindBin(iStrip);
          auto vfat = (8-iEta)+8*((iStrip-1)/128);
          double val = hAllStrip[iChamber][iLayer]->GetBinContent(xBin, yBin);
          if (val < 1) continue;
          int nEvent_ = hInstLumi[iChamber][iLayer][vfat]->GetEntries();
          val = double(val/nEvent_);
          if (val > cut[iEta-1]) continue;
          hMaskedStrip[iChamber][iLayer]->SetBinContent(xBin, yBin, 0);
          nOfActiveStrip[iChamber][iLayer][iEta-1]++;
          nOfActiveStripVfat[iChamber][iLayer][vfat]++;
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
        for ( int ix = 0; ix < 3; ix++) {
          auto vfat = (8-et+1)+8*(ix);
          vfatArea[ch][ly][vfat] = stripLength[(ch+1)%2][et] * stripPitch[(ch+1)%2][et] * (nOfActiveStripVfat[ch][ly][vfat]);
        }
      }
      totArea += chArea[ch][ly];
    }
  }

  cout << "test\n";
  bool flag;
  iHit = 0;
  for (int i = 0; i < nEvent; i++) {
    t_Event->GetEntry(i);
    flag = false;
    auto lumiBin = hActiveLumi[run]->GetXaxis()->FindBin(lumi);
    for (int j = 0; j < nGEMHits; j++) {
      t_Hit->GetEntry(iHit);
      iHit++;
      if(hLumiStatus[run][lumi]->GetRMS() > 0.01) continue;
      auto vfat = (8-etaPartition)+8*((firstStrip-1)/128);
      auto chanBin = hActiveLumi[run]->GetYaxis()->FindBin(chamber+(layer-1)*0.5+vfat/48.);
      if (hActiveLumi[run]->Integral(lumiBin, lumiBin, chanBin, chanBin) < 1) continue;

      double chArea_ = 0;
      for (int vfat_ = 0; vfat_ < 24; vfat_++) {
        auto chanBin = hActiveLumi[run]->GetYaxis()->FindBin(chamber+(layer-1)*0.5+vfat_/48.);
        if (hActiveLumi[run]->Integral(lumiBin, lumiBin, chanBin, chanBin) > 0)
          chArea_ += vfatArea[chamber-1][layer-1][vfat_];
      }
      double etaArea_ = 0;
      for (int ix = 0; ix < 3; ix++) {
        auto vfat_ = (8-etaPartition)+8*ix;
        auto chanBin = hActiveLumi[run]->GetYaxis()->FindBin(chamber+(layer-1)*0.5+vfat_/48.);
        if (hActiveLumi[run]->Integral(lumiBin, lumiBin, chanBin, chanBin) > 0)
          etaArea_ += vfatArea[chamber-1][layer-1][vfat_];
      }

      int nEvent_ = hInstLumi[chamber-1][layer-1][vfat]->GetBinContent(hInstLumi[chamber-1][layer-1][vfat]->GetXaxis()->FindBin(instLumi));
      auto etaBin = hAllStrip[chamber-1][layer-1]->GetYaxis()->FindBin(etaPartition);
      
      double normFactor = (chArea_*nEvent_*25E-9);

      auto xBin = hMaskedStrip[chamber-1][layer-1]->GetXaxis()->FindBin(firstStrip);
      if(hMaskedStrip[chamber-1][layer-1]->GetBinContent(xBin, etaBin) != 0) continue;
      hAllStrips[chamber-1][layer-1]->Fill(firstStrip, etaPartition);
      hHitRateLumi[chamber-1][layer-1]->Fill(instLumi, 1./(normFactor));
      hHitRateEta[chamber-1][layer-1][etaPartition-1]->Fill(instLumi, 1./(etaArea_*nEvent_*25E-9));
      double normFactorVfat = (vfatArea[chamber-1][layer-1][vfat]*nEvent_*25E-9);
      hHitRateVfat[chamber-1][layer-1]->Fill(firstStrip, etaPartition, 1./(normFactorVfat));
      flag = true;
    }
  }
  cout << "test\n";

  output->Write();
  output->Close();
}


array<double,8> findCuts(TH2D* hAllStrip[36][2], TH1D* hInstLumi[36][2]) {
  TH1D* hitHist[8];
  for (int iEta = 0; iEta < 8; iEta++) {
    hitHist[iEta] = new TH1D( Form("hitHist eta %i", iEta+1), 
                             "hits per each strip", 500, 0, 1e-4);
    hitHist[iEta]->GetXaxis()->SetTitle("Average number of hits on Strip");
    hitHist[iEta]->GetYaxis()->SetTitle("strips");
  }
  array<double, 8> res;
  for (int iChamber = 26; iChamber < 30; iChamber++) {
    for (int iLayer = 0; iLayer < 2; iLayer++) {
      if (hAllStrip[iChamber][iLayer]->GetEntries() == 0) continue;
      auto xAxis = hAllStrip[iChamber][iLayer]->GetXaxis();
      auto yAxis = hAllStrip[iChamber][iLayer]->GetYaxis();
      for (int iEta = 1; iEta < 9; iEta++) {
        auto yBin = yAxis->FindBin(iEta);
        for (int iStrip = 1; iStrip < 386; iStrip++) {
          auto nEvent_ = hInstLumi[iChamber][iLayer]->GetEntries();
          auto xBin = xAxis->FindBin(iStrip);
          double val = hAllStrip[iChamber][iLayer]->GetBinContent(xBin, yBin);
          val = double(val/nEvent_);
          if (hAllStrip[iChamber][iLayer]->GetBinContent(xBin, yBin) == 0) continue;
          hitHist[iEta-1]->Fill((double)val);
        }
      }
    }
  }
  for (int iEta = 0; iEta < 8; iEta++) {
    res[iEta] = (hitHist[iEta]->GetMean() + 2*hitHist[iEta]->GetRMS());
  }
  return res; 
}

void rebinHisto(TH1D* h) {
  double bins[19];
  for (int i = 0; i < 6; i++) {
    bins[i] = i * 0.2;
  }
  for (int i = 0; i < 13; i++) {
    bins[i+6] = i * 0.05 + 1.2;
  }
  h->SetBins(18, bins);
}

array<double,8> findCuts(TH2D* hAllStrip[36][2], TH1D* hInstLumi[36][2][24]) {
  TH1D* hitHist[8];
  for (int iEta = 0; iEta < 8; iEta++) {
    hitHist[iEta] = new TH1D( Form("hitHist eta %i", iEta+1), 
                             "hits per each strip", 500, 0, 1e-4);
    hitHist[iEta]->GetXaxis()->SetTitle("Average number of hits on Strip");
    hitHist[iEta]->GetYaxis()->SetTitle("strips");
  }
  array<double,8> res;
  for (int iChamber = 26; iChamber < 30; iChamber++) {
    for (int iLayer = 0; iLayer < 2; iLayer++) {
      if (hAllStrip[iChamber][iLayer]->GetEntries() == 0) continue;
      auto xAxis = hAllStrip[iChamber][iLayer]->GetXaxis();
      auto yAxis = hAllStrip[iChamber][iLayer]->GetYaxis();
      for (int iEta = 1; iEta < 9; iEta++) {
        auto yBin = yAxis->FindBin(iEta);
        for (int iStrip = 1; iStrip < 385; iStrip++) {
          auto vfat = (8-iEta)+8*((iStrip-1)/128);
          auto nEvent_ = hInstLumi[iChamber][iLayer][vfat]->GetEntries();
          auto xBin = xAxis->FindBin(iStrip);
          double val = hAllStrip[iChamber][iLayer]->GetBinContent(xBin, yBin);
          val = double(val/nEvent_);
          if (hAllStrip[iChamber][iLayer]->GetBinContent(xBin, yBin) == 0) continue;
          hitHist[iEta-1]->Fill((double)val);
        }
      }
    }
  }
  for (int iEta = 0; iEta < 8; iEta++) {
    res[iEta] = (hitHist[iEta]->GetMean() + 2*hitHist[iEta]->GetRMS());
  }
  return res; 
}

