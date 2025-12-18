// plot_reco.C
// To run: root -l plot_reco.C

void plot_reco() {
    // Open the simulation output
    TFile* file = new TFile("../buildPos5cm/B4.root", "READ");
    if (!file || file->IsZombie()) {
        cout << "Error: Could not open B4.root" << endl;
        return;
    }
    
    // Get the reconstructed positions tree
    TTree* tree = (TTree*)file->Get("ReconstructedPositions");
    if (!tree) {
        cout << "Error: Could not find ReconstructedPositions tree" << endl;
        file->Close();
        return;
    }
    
    cout << "Found tree with " << tree->GetEntries() << " entries" << endl;
    
    // Set branch addresses
    Float_t xReco, yReco, nTotal;
    Int_t eventID;
    
    tree->SetBranchAddress("EventID", &eventID);
    tree->SetBranchAddress("X_reconstructed", &xReco);
    tree->SetBranchAddress("Y_reconstructed", &yReco);
    tree->SetBranchAddress("N_total", &nTotal);
    
    // Create histograms
    TH1D* hX = new TH1D("hX", "Reconstructed X Position;X (mm);Counts", 100, -50, 50);
    TH1D* hY = new TH1D("hY", "Reconstructed Y Position;Y (mm);Counts", 100, -50, 50);
    TH2D* hXY = new TH2D("hXY", "2D Reconstructed Positions;X (mm);Y (mm)", 
                         100, -50, 50, 100, -50, 50);
    TH1D* hN = new TH1D("hN", "Photons per Event;N_{photons};Counts", 100, 0, 500);
    
    // Fill histograms
    Long64_t nEntries = tree->GetEntries();
    Long64_t nEventsWithPhotons = 0;
    
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        
        if (nTotal > 0) {
            hX->Fill(xReco);
            hY->Fill(yReco);
            hXY->Fill(xReco, yReco);
            hN->Fill(nTotal);
            nEventsWithPhotons++;
        }
        
        if (i % 10000 == 0) {
            cout << "Processing event " << i << " of " << nEntries << endl;
        }
    }
    
    cout << "\n=== Reconstruction Statistics ===" << endl;
    cout << "Total events in file: " << nEntries << endl;
    cout << "Events with photons: " << nEventsWithPhotons << endl;
    cout << "X position - Mean: " << hX->GetMean() << " mm, RMS: " << hX->GetRMS() << " mm" << endl;
    cout << "Y position - Mean: " << hY->GetMean() << " mm, RMS: " << hY->GetRMS() << " mm" << endl;
    cout << "Average photons per event: " << hN->GetMean() << endl;
    
    // Set style
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    gStyle->SetPalette(kRainBow);
    
    // Create main canvas
    TCanvas* c1 = new TCanvas("c1", "Reconstruction Results", 1400, 1000);
    c1->Divide(3, 2);
    
    // Plot 1: X distribution with Gaussian fit
    c1->cd(1);
    hX->SetFillColor(kBlue);
    hX->SetFillStyle(3001);
    hX->Draw();
    TF1* fgausX = new TF1("fgausX", "gaus", -50, 50);
    hX->Fit(fgausX, "RQ");
    fgausX->SetLineColor(kRed);
    fgausX->Draw("same");
    
    // Plot 2: Y distribution with Gaussian fit
    c1->cd(2);
    hY->SetFillColor(kGreen);
    hY->SetFillStyle(3001);
    hY->Draw();
    TF1* fgausY = new TF1("fgausY", "gaus", -50, 50);
    hY->Fit(fgausY, "RQ");
    fgausY->SetLineColor(kRed);
    fgausY->Draw("same");
    
    // Plot 3: 2D histogram
    c1->cd(3);
    hXY->Draw("COLZ");
    
    // Plot 4: Photon count distribution
    c1->cd(4);
    hN->SetFillColor(kRed);
    hN->SetFillStyle(3001);
    hN->Draw();
    hN->Fit("expo", "RQ", "", hN->GetMean()/2, hN->GetMean()*2);
    
    // Plot 5: Profile histogram
    c1->cd(5);
    TProfile2D* profXY = new TProfile2D("profXY", "Profile XY;X (mm);Y (mm);Mean N_{photons}", 
                                        25, -50, 50, 25, -50, 50);
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        if (nTotal > 0) {
            profXY->Fill(xReco, yReco, nTotal);
        }
    }
    profXY->Draw("COLZ");
    
    // Plot 6: Scatter plot (first 1000 events)
    c1->cd(6);
    TGraph* gr = new TGraph();
    gr->SetTitle("Scatter Plot (first 1000 events);X (mm);Y (mm)");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.5);
    gr->SetMarkerColor(kBlue);
    
    Int_t nPoints = TMath::Min(nEventsWithPhotons, (Long64_t)1000);
    for (Int_t i = 0; i < nPoints; i++) {
        tree->GetEntry(i);
        if (nTotal > 0) {
            gr->SetPoint(gr->GetN(), xReco, yReco);
        }
    }
    gr->Draw("AP");
    gr->GetHistogram()->SetMinimum(-50);
    gr->GetHistogram()->SetMaximum(50);
    
    c1->SaveAs("reconstruction_summary.png");
    
    // Create detailed 2D plot
    TCanvas* c2 = new TCanvas("c2", "Detailed 2D Reconstruction", 800, 800);
    hXY->Draw("COLZ");
    hXY->SetContour(20);
    hXY->Draw("CONTZ SAME");
    
    // Add statistics box
    TPaveText* stats = new TPaveText(0.15, 0.85, 0.45, 0.95, "NDC");
    stats->SetFillColor(0);
    stats->SetBorderSize(1);
    stats->AddText(Form("Events: %lld", nEventsWithPhotons));
    stats->AddText(Form("X: %.2f #pm %.2f mm", hX->GetMean(), hX->GetRMS()));
    stats->AddText(Form("Y: %.2f #pm %.2f mm", hY->GetMean(), hY->GetRMS()));
    stats->AddText(Form("<N>: %.1f photons", hN->GetMean()));
    stats->Draw();
    
    c2->SaveAs("detailed_2d_reconstruction.png");
    
    // Create resolution plot
    TCanvas* c3 = new TCanvas("c3", "Resolution Analysis", 1200, 400);
    c3->Divide(3, 1);
    
    c3->cd(1);
    TH1D* hResX = new TH1D("hResX", "X Position Residuals;X (mm);Counts", 100, -10, 10);
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        if (nTotal > 0) {
            hResX->Fill(xReco - hX->GetMean());
        }
    }
    hResX->SetFillColor(kBlue);
    hResX->SetFillStyle(3001);
    hResX->Draw();
    hResX->Fit("gaus", "RQ");
    
    c3->cd(2);
    TH1D* hResY = new TH1D("hResY", "Y Position Residuals;Y (mm);Counts", 100, -10, 10);
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        if (nTotal > 0) {
            hResY->Fill(yReco - hY->GetMean());
        }
    }
    hResY->SetFillColor(kGreen);
    hResY->SetFillStyle(3001);
    hResY->Draw();
    hResY->Fit("gaus", "RQ");
    
    c3->cd(3);
    TH2D* hResXY = new TH2D("hResXY", "2D Residuals;X (mm);Y (mm)", 
                           50, -10, 10, 50, -10, 10);
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        if (nTotal > 0) {
            hResXY->Fill(xReco - hX->GetMean(), yReco - hY->GetMean());
        }
    }
    hResXY->Draw("COLZ");
    
    c3->SaveAs("resolution_analysis.png");
    
    // Save histograms to file
    TFile* outFile = new TFile("reconstruction_results.root", "RECREATE");
    hX->Write();
    hY->Write();
    hXY->Write();
    hN->Write();
    hResX->Write();
    hResY->Write();
    hResXY->Write();
    profXY->Write();
    c1->Write();
    c2->Write();
    c3->Write();
    outFile->Close();
    
    cout << "\n=== Analysis Complete ===" << endl;
    cout << "Plots saved as:" << endl;
    cout << "1. reconstruction_summary.png - Multi-panel summary" << endl;
    cout << "2. detailed_2d_reconstruction.png - Detailed 2D plot" << endl;
    cout << "3. resolution_analysis.png - Resolution analysis" << endl;
    cout << "4. reconstruction_results.root - ROOT file with histograms" << endl;
    
    // Don't close file if you want to keep working with it
    // file->Close();
}
