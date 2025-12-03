// proper_reconstruction.C
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <utility>

using namespace std;

// Forward declarations
double calculateMean(const vector<double>& values);
pair<double, double> calculateStatistics(const vector<double>& values);
double calculateWeightedMeanX(double P_x1, double N_x1, double sigma_x1, 
                             double P_x2, double N_x2, double sigma_x2);

void proper_reconstruction() {
    ifstream infile("../buildPosRec/position_data.txt");
    if (!infile.is_open()) {
        cout << "Error opening file!" << endl;
        return;
    }
    
    // Output histograms
    TH1D* histReconX = new TH1D("histReconX", "Reconstructed X Position", 100, -50, 50);
    TH1D* histReconY = new TH1D("histReconY", "Reconstructed Y Position", 100, -50, 50);
    TH2D* histReconXY = new TH2D("histReconXY", "Reconstructed X vs Y", 100, -50, 50, 100, -50, 50);
    
    // For tracking
    vector<double> allReconX, allReconY;
    int totalEvents = 0;
    int eventsWithX = 0;
    int eventsWithY = 0;
    int eventsWithBoth = 0;
    
    cout << "Starting proper reconstruction..." << endl;
    
    int eventID, nBottom, nTop, nLeft, nRight;
    double value;
    
    while (infile >> eventID >> nBottom >> nTop >> nLeft >> nRight) {
        totalEvents++;
        
        // Read positions
        vector<double> xBottom, xTop, yLeft, yRight;
        
        for (int i = 0; i < nBottom; i++) {
            infile >> value;
            xBottom.push_back(value);
        }
        for (int i = 0; i < nTop; i++) {
            infile >> value;
            xTop.push_back(value);
        }
        for (int i = 0; i < nLeft; i++) {
            infile >> value;
            yLeft.push_back(value);
        }
        for (int i = 0; i < nRight; i++) {
            infile >> value;
            yRight.push_back(value);
        }
        
        // ========== PROPER RECONSTRUCTION ==========
        
        double reconX = -999;
        double reconY = -999;
        
        // Reconstruct X coordinate
        if (!xBottom.empty() && !xTop.empty()) {
            // Both X arrays available - use simple mean for now
            reconX = (calculateMean(xBottom) + calculateMean(xTop)) / 2.0;
            eventsWithX++;
        }
        else if (!xBottom.empty()) {
            // Only bottom array
            reconX = calculateMean(xBottom);
            eventsWithX++;
        }
        else if (!xTop.empty()) {
            // Only top array
            reconX = calculateMean(xTop);
            eventsWithX++;
        }
        
        // Reconstruct Y coordinate
        if (!yLeft.empty() && !yRight.empty()) {
            // Both Y arrays available - use simple mean for now
            reconY = (calculateMean(yLeft) + calculateMean(yRight)) / 2.0;
            eventsWithY++;
        }
        else if (!yLeft.empty()) {
            // Only left array
            reconY = calculateMean(yLeft);
            eventsWithY++;
        }
        else if (!yRight.empty()) {
            // Only right array
            reconY = calculateMean(yRight);
            eventsWithY++;
        }
        
        // Fill histograms if we have at least one coordinate
        if (reconX > -900) {
            histReconX->Fill(reconX);
            allReconX.push_back(reconX);
        }
        if (reconY > -900) {
            histReconY->Fill(reconY);
            allReconY.push_back(reconY);
        }
        if (reconX > -900 && reconY > -900) {
            histReconXY->Fill(reconX, reconY);
            eventsWithBoth++;
        }
        
        // Progress
        if (totalEvents % 1000 == 0) {
            cout << "Processed " << totalEvents << " events..." << endl;
        }
    }
    
    infile.close();
    
    cout << "\n==========================================" << endl;
    cout << "PROPER RECONSTRUCTION RESULTS" << endl;
    cout << "==========================================" << endl;
    cout << "Total events: " << totalEvents << endl;
    cout << "Events with X reconstruction: " << eventsWithX 
         << " (" << (100.0 * eventsWithX / totalEvents) << "%)" << endl;
    cout << "Events with Y reconstruction: " << eventsWithY
         << " (" << (100.0 * eventsWithY / totalEvents) << "%)" << endl;
    cout << "Events with BOTH coordinates: " << eventsWithBoth
         << " (" << (100.0 * eventsWithBoth / totalEvents) << "%)" << endl;
    
    if (!allReconX.empty()) {
        pair<double, double> xStats = calculateStatistics(allReconX);
        cout << "\nX Resolution:" << endl;
        cout << "  Mean: " << xStats.first << " mm" << endl;
        cout << "  Sigma: " << xStats.second << " mm" << endl;
        cout << "  FWHM: " << 2.355 * xStats.second << " mm" << endl;
    }
    
    if (!allReconY.empty()) {
        pair<double, double> yStats = calculateStatistics(allReconY);
        cout << "\nY Resolution:" << endl;
        cout << "  Mean: " << yStats.first << " mm" << endl;
        cout << "  Sigma: " << yStats.second << " mm" << endl;
        cout << "  FWHM: " << 2.355 * yStats.second << " mm" << endl;
    }
    
    if (!allReconX.empty() && !allReconY.empty()) {
        pair<double, double> xStats = calculateStatistics(allReconX);
        pair<double, double> yStats = calculateStatistics(allReconY);
        double spatialRes = sqrt(xStats.second * xStats.second + yStats.second * yStats.second);
        cout << "\nOverall Spatial Resolution: " << spatialRes << " mm (sigma)" << endl;
        cout << "Overall Spatial Resolution: " << 2.355 * spatialRes << " mm (FWHM)" << endl;
    }
    
    // Gaussian fit for publication-quality results
    if (!allReconX.empty()) {
        TF1* gaussX = new TF1("gaussX", "gaus", -50, 50);
        histReconX->Fit(gaussX, "Q");
        cout << "\nX Gaussian Fit:" << endl;
        cout << "  Mean: " << gaussX->GetParameter(1) << " ± " << gaussX->GetParError(1) << " mm" << endl;
        cout << "  Sigma: " << gaussX->GetParameter(2) << " ± " << gaussX->GetParError(2) << " mm" << endl;
    }
    
    if (!allReconY.empty()) {
        TF1* gaussY = new TF1("gaussY", "gaus", -50, 50);
        histReconY->Fit(gaussY, "Q");
        cout << "\nY Gaussian Fit:" << endl;
        cout << "  Mean: " << gaussY->GetParameter(1) << " ± " << gaussY->GetParError(1) << " mm" << endl;
        cout << "  Sigma: " << gaussY->GetParameter(2) << " ± " << gaussY->GetParError(2) << " mm" << endl;
    }
    
    // Save results
    TFile* outFile = new TFile("proper_reconstruction.root", "RECREATE");
    histReconX->Write();
    histReconY->Write();
    histReconXY->Write();
    outFile->Close();
    
    // Draw
    TCanvas* c1 = new TCanvas("c1", "Proper Reconstruction", 1200, 800);
    c1->Divide(2, 2);
    
    c1->cd(1);
    histReconX->SetLineColor(kBlue);
    histReconX->SetLineWidth(2);
    histReconX->Draw();
    
    c1->cd(2);
    histReconY->SetLineColor(kRed);
    histReconY->SetLineWidth(2);
    histReconY->Draw();
    
    c1->cd(3);
    histReconXY->Draw("colz");
    
    c1->SaveAs("proper_reconstruction.png");
    
    cout << "\nResults saved to proper_reconstruction.root and proper_reconstruction.png" << endl;
}

// Helper functions
double calculateMean(const vector<double>& values) {
    if (values.empty()) return 0;
    double sum = 0;
    for (double v : values) sum += v;
    return sum / values.size();
}

pair<double, double> calculateStatistics(const vector<double>& values) {
    if (values.empty()) return make_pair(0, 0);
    
    double sum = 0, sum2 = 0;
    for (double v : values) {
        sum += v;
        sum2 += v * v;
    }
    
    double mean = sum / values.size();
    double sigma = sqrt((sum2 / values.size()) - (mean * mean));
    
    return make_pair(mean, sigma);
}

// Your original weighted mean functions (for when you have both arrays)
double calculateWeightedMeanX(double P_x1, double N_x1, double sigma_x1, 
                             double P_x2, double N_x2, double sigma_x2) {
    if (sigma_x1 == 0 || sigma_x2 == 0) {
        return (P_x1 + P_x2) / 2.0; // Fallback to simple mean
    }
    double numerator = (P_x1 * N_x1 / sigma_x1) + (P_x2 * N_x2 / sigma_x2);
    double denominator = (N_x1 / sigma_x1) + (N_x2 / sigma_x2);
    return numerator / denominator;
}
