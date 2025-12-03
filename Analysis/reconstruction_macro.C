// reconstruct_from_txt.C
// Run with: root -l reconstruct_from_txt.C
{
    // Open the position data file
    ifstream infile("../buildPosRec/position_data.txt");
    if (!infile.is_open()) {
        cout << "Error: Could not open position_data.txt" << endl;
        return;
    }
    
    // Output histograms
    TH1D* histReconX = new TH1D("histReconX", "Reconstructed X Position", 100, -50, 50);
    TH1D* histReconY = new TH1D("histReconY", "Reconstructed Y Position", 100, -50, 50);
    TH2D* histReconXY = new TH2D("histReconXY", "Reconstructed X vs Y", 100, -50, 50, 100, -50, 50);
    
    // For final resolution calculation
    vector<double> allReconX, allReconY;
    int totalEvents = 0;
    int reconstructedEvents = 0;
    
    cout << "Processing position data..." << endl;
    
    // Variables to read
    int eventID, nBottom, nTop, nLeft, nRight;
    double value;
    
    // Read each event from the file
    while (infile >> eventID >> nBottom >> nTop >> nLeft >> nRight) {
        totalEvents++;
        
        // Create vectors for this event
        vector<double> xBottom, xTop, yLeft, yRight;
        
        // Read bottom positions
        for (int i = 0; i < nBottom; i++) {
            if (!(infile >> value)) break;
            xBottom.push_back(value);
        }
        
        // Read top positions
        for (int i = 0; i < nTop; i++) {
            if (!(infile >> value)) break;
            xTop.push_back(value);
        }
        
        // Read left positions
        for (int i = 0; i < nLeft; i++) {
            if (!(infile >> value)) break;
            yLeft.push_back(value);
        }
        
        // Read right positions
        for (int i = 0; i < nRight; i++) {
            if (!(infile >> value)) break;
            yRight.push_back(value);
        }
        
        // Skip events without enough data
        if (xBottom.empty() || xTop.empty() || yLeft.empty() || yRight.empty()) {
            continue;
        }
        
        // ========== PERFORM RECONSTRUCTION ==========
        
        // Create temporary histograms for fitting
        TH1D* histXBottom = new TH1D("histXBottom_temp", "X Bottom Array", 50, -50, 50);
        TH1D* histXTop = new TH1D("histXTop_temp", "X Top Array", 50, -50, 50);
        TH1D* histYLeft = new TH1D("histYLeft_temp", "Y Left Array", 50, -50, 50);
        TH1D* histYRight = new TH1D("histYRight_temp", "Y Right Array", 50, -50, 50);
        
        // Fill histograms
        for (double x : xBottom) histXBottom->Fill(x);
        for (double x : xTop) histXTop->Fill(x);
        for (double y : yLeft) histYLeft->Fill(y);
        for (double y : yRight) histYRight->Fill(y);
        
        // Fit Gaussian distributions
        TF1* gaussFitXBottom = new TF1("gaussFitXBottom_temp", "gaus", -50, 50);
        TF1* gaussFitXTop = new TF1("gaussFitXTop_temp", "gaus", -50, 50);
        TF1* gaussFitYLeft = new TF1("gaussFitYLeft_temp", "gaus", -50, 50);
        TF1* gaussFitYRight = new TF1("gaussFitYRight_temp", "gaus", -50, 50);
        
        histXBottom->Fit(gaussFitXBottom, "QN");
        histXTop->Fit(gaussFitXTop, "QN");
        histYLeft->Fit(gaussFitYLeft, "QN");
        histYRight->Fit(gaussFitYRight, "QN");
        
        // Get fit parameters
        double meanXBottom = gaussFitXBottom->GetParameter(1);
        double sigmaXBottom = gaussFitXBottom->GetParameter(2);
        double ampXBottom = gaussFitXBottom->GetParameter(0);
        
        double meanXTop = gaussFitXTop->GetParameter(1);
        double sigmaXTop = gaussFitXTop->GetParameter(2);
        double ampXTop = gaussFitXTop->GetParameter(0);
        
        double meanYLeft = gaussFitYLeft->GetParameter(1);
        double sigmaYLeft = gaussFitYLeft->GetParameter(2);
        double ampYLeft = gaussFitYLeft->GetParameter(0);
        
        double meanYRight = gaussFitYRight->GetParameter(1);
        double sigmaYRight = gaussFitYRight->GetParameter(2);
        double ampYRight = gaussFitYRight->GetParameter(0);
        
        // Use weighted mean calculation for reconstruction
        double reconX = calculateWeightedMeanX(meanXBottom, ampXBottom, sigmaXBottom, 
                                              meanXTop, ampXTop, sigmaXTop);
        double reconY = calculateWeightedMeanY(meanYLeft, ampYLeft, sigmaYLeft,
                                              meanYRight, ampYRight, sigmaYRight);
        
        // Check if reconstruction is valid
        if (reconX > -50 && reconX < 50 && reconY > -50 && reconY < 50) {
            histReconX->Fill(reconX);
            histReconY->Fill(reconY);
            histReconXY->Fill(reconX, reconY);
            allReconX.push_back(reconX);
            allReconY.push_back(reconY);
            reconstructedEvents++;
        }
        
        // Cleanup temporary objects
        delete histXBottom;
        delete histXTop;
        delete histYLeft;
        delete histYRight;
        delete gaussFitXBottom;
        delete gaussFitXTop;
        delete gaussFitYLeft;
        delete gaussFitYRight;
        
        // Progress indicator
        if (totalEvents % 100 == 0) {
            cout << "Processed " << totalEvents << " events, reconstructed " 
                 << reconstructedEvents << " events" << endl;
        }
    }
    
    infile.close();
    
    cout << "\n==========================================" << endl;
    cout << "PROCESSING COMPLETE" << endl;
    cout << "==========================================" << endl;
    cout << "Total events in file: " << totalEvents << endl;
    cout << "Successfully reconstructed: " << reconstructedEvents << endl;
    cout << "Reconstruction efficiency: " 
         << (100.0 * reconstructedEvents / totalEvents) << "%" << endl;
    
    // Calculate final resolution if we have data
    if (!allReconX.empty()) {
        calculateFinalResolution(allReconX, allReconY);
    }
    
    // Save results to ROOT file
    TFile* outFile = new TFile("reconstruction_from_txt.root", "RECREATE");
    histReconX->Write();
    histReconY->Write();
    histReconXY->Write();
    outFile->Close();
    
    // Create visualization
    TCanvas* c1 = new TCanvas("c1", "Position Reconstruction Results", 1200, 800);
    c1->Divide(2, 2);
    
    c1->cd(1);
    histReconX->SetLineColor(kBlue);
    histReconX->SetLineWidth(2);
    histReconX->GetXaxis()->SetTitle("X Position (mm)");
    histReconX->GetYaxis()->SetTitle("Counts");
    histReconX->Draw();
    
    c1->cd(2);
    histReconY->SetLineColor(kRed);
    histReconY->SetLineWidth(2);
    histReconY->GetXaxis()->SetTitle("Y Position (mm)");
    histReconY->GetYaxis()->SetTitle("Counts");
    histReconY->Draw();
    
    c1->cd(3);
    histReconXY->GetXaxis()->SetTitle("X Position (mm)");
    histReconXY->GetYaxis()->SetTitle("Y Position (mm)");
    histReconXY->Draw("colz");
    
    c1->cd(4);
    // Add statistics text
    TPaveText* stats = new TPaveText(0.1, 0.1, 0.9, 0.9);
    stats->AddText(Form("Total Events: %d", totalEvents));
    stats->AddText(Form("Reconstructed: %d", reconstructedEvents));
    stats->AddText(Form("Efficiency: %.1f%%", 100.0 * reconstructedEvents / totalEvents));
    if (!allReconX.empty()) {
        // Calculate mean positions
        double meanX = 0, meanY = 0;
        for (size_t i = 0; i < allReconX.size(); i++) {
            meanX += allReconX[i];
            meanY += allReconY[i];
        }
        meanX /= allReconX.size();
        meanY /= allReconY.size();
        stats->AddText(Form("Mean X: %.2f mm", meanX));
        stats->AddText(Form("Mean Y: %.2f mm", meanY));
    }
    stats->Draw();
    
    c1->SaveAs("reconstruction_results.png");
    
    cout << "\nResults saved to:" << endl;
    cout << "  - reconstruction_from_txt.root (ROOT file with histograms)" << endl;
    cout << "  - reconstruction_results.png (visualization)" << endl;
    
    // Cleanup
    delete c1;
}

// ========== YOUR ORIGINAL FUNCTIONS ==========

double calculateWeightedMeanX(double P_x1, double N_x1, double sigma_x1, 
                             double P_x2, double N_x2, double sigma_x2) {
    if (sigma_x1 == 0 || sigma_x2 == 0) {
        return -1000; // Invalid value
    }
    double numerator = (P_x1 * N_x1 / sigma_x1) + (P_x2 * N_x2 / sigma_x2);
    double denominator = (N_x1 / sigma_x1) + (N_x2 / sigma_x2);
    return numerator / denominator;
}

double calculateWeightedMeanY(double P_y1, double N_y1, double sigma_y1, 
                             double P_y2, double N_y2, double sigma_y2) {
    if (sigma_y1 == 0 || sigma_y2 == 0) {
        return -1000; // Invalid value
    }
    double numerator = (P_y1 * N_y1 / sigma_y1) + (P_y2 * N_y2 / sigma_y2);
    double denominator = (N_y1 / sigma_y1) + (N_y2 / sigma_y2);
    return numerator / denominator;
}

void calculateFinalResolution(const vector<double>& allReconX, 
                            const vector<double>& allReconY) {
    if (allReconX.empty() || allReconY.empty()) {
        cout << "No data for final resolution calculation!" << endl;
        return;
    }
    
    // Create histograms from all reconstructed positions
    TH1D* histFinalX = new TH1D("histFinalX", "Final Reconstructed X", 100, -50, 50);
    TH1D* histFinalY = new TH1D("histFinalY", "Final Reconstructed Y", 100, -50, 50);
    
    for (size_t i = 0; i < allReconX.size(); i++) {
        histFinalX->Fill(allReconX[i]);
        histFinalY->Fill(allReconY[i]);
    }
    
    // Perform Gaussian fits
    TF1* gaussFitX = new TF1("gaussFitX", "gaus", -50, 50);
    TF1* gaussFitY = new TF1("gaussFitY", "gaus", -50, 50);
    
    histFinalX->Fit(gaussFitX, "Q");
    histFinalY->Fit(gaussFitY, "Q");
    
    // Get resolution parameters
    double meanX = gaussFitX->GetParameter(1);
    double sigmaX = gaussFitX->GetParameter(2);
    double meanY = gaussFitY->GetParameter(1);
    double sigmaY = gaussFitY->GetParameter(2);
    
    // Calculate FWHM and spatial resolution
    double fwhmX = 2.355 * sigmaX;
    double fwhmY = 2.355 * sigmaY;
    double spatialResolution = TMath::Sqrt(sigmaX * sigmaX + sigmaY * sigmaY);
    double spatialResolutionFWHM = 2.355 * spatialResolution;
    
    cout << "\n=============================================" << endl;
    cout << "      FINAL SPATIAL RESOLUTION RESULTS" << endl;
    cout << "=============================================" << endl;
    cout << "Based on " << allReconX.size() << " reconstructed events" << endl;
    cout << "X Position: Mean = " << meanX << " mm, Sigma = " << sigmaX << " mm" << endl;
    cout << "Y Position: Mean = " << meanY << " mm, Sigma = " << sigmaY << " mm" << endl;
    cout << "X FWHM: " << fwhmX << " mm" << endl;
    cout << "Y FWHM: " << fwhmY << " mm" << endl;
    cout << "Spatial Resolution (sigma): " << spatialResolution << " mm" << endl;
    cout << "Spatial Resolution (FWHM): " << spatialResolutionFWHM << " mm" << endl;
    cout << "=============================================" << endl;
    
    // Save final resolution histograms
    TFile* resFile = new TFile("final_resolution.root", "RECREATE");
    histFinalX->Write();
    histFinalY->Write();
    gaussFitX->Write();
    gaussFitY->Write();
    resFile->Close();
    
    // Cleanup
    delete histFinalX;
    delete histFinalY;
    delete gaussFitX;
    delete gaussFitY;
}
