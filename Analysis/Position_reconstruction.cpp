// position_reconstruction_weighted_only.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>

// ROOT includes
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TProfile2D.h"
using namespace std;

// Structure to hold event data
struct EventData {
    int eventID;
    vector<double> xBottom;  // positions from bottom array
    vector<double> xTop;     // positions from top array  
    vector<double> yLeft;    // positions from left array
    vector<double> yRight;   // positions from right array
    
    EventData(int id) : eventID(id) {}
    
    void addXBottom(double x) { xBottom.push_back(x); }
    void addXTop(double x) { xTop.push_back(x); }
    void addYLeft(double y) { yLeft.push_back(y); }
    void addYRight(double y) { yRight.push_back(y); }
    
    void getArrayStats(const vector<double>& positions, 
                       double& mean, double& sigma, double& N) {
        N = positions.size();
        if (N == 0) {
            mean = 0.0;
            sigma = 1.0;
            return;
        }
        
        double sum = 0.0;
        for (double pos : positions) sum += pos;
        mean = sum / N;
        
        double sumSq = 0.0;
        for (double pos : positions) sumSq += (pos - mean) * (pos - mean);
        sigma = sqrt(sumSq / N);
        
        if (sigma < 0.1) sigma = 0.1;
    }
    
    void getGaussianStats(const vector<double>& positions,
                         double& mean, double& sigma, double& N) {
        N = positions.size();
        if (N < 5) {
            getArrayStats(positions, mean, sigma, N);
            return;
        }
        
        double minVal = *min_element(positions.begin(), positions.end());
        double maxVal = *max_element(positions.begin(), positions.end());
        double range = maxVal - minVal;
        
        minVal -= 0.1 * range;
        maxVal += 0.1 * range;
        
        string histName = "h_temp_" + to_string(eventID);
        TH1D* hist = new TH1D(histName.c_str(), "Temp histogram", 100, minVal, maxVal);
        
        for (double pos : positions) hist->Fill(pos);
        
        TF1* gauss = new TF1("gauss", "gaus", minVal, maxVal);
        gauss->SetParameters(N/10.0, (minVal + maxVal)/2.0, range/4.0);
        hist->Fit(gauss, "QN");
        
        mean = gauss->GetParameter(1);
        sigma = gauss->GetParameter(2);
        
        delete hist;
        delete gauss;
        
        if (sigma <= 0 || !isfinite(sigma)) {
            getArrayStats(positions, mean, sigma, N);
        }
    }
    
    // Weighted reconstruction as per paper formula
    void reconstructWeighted(double& x, double& y) {
        double Px1, sigmaX1, Nx1;  // Bottom array
        double Px2, sigmaX2, Nx2;  // Top array
        double Py1, sigmaY1, Ny1;  // Left array
        double Py2, sigmaY2, Ny2;  // Right array
        
        getGaussianStats(xBottom, Px1, sigmaX1, Nx1);
        getGaussianStats(xTop, Px2, sigmaX2, Nx2);
        getGaussianStats(yLeft, Py1, sigmaY1, Ny1);
        getGaussianStats(yRight, Py2, sigmaY2, Ny2);
        
        // Paper's formula for X
        double numeratorX = 0.0;
        double denominatorX = 0.0;
        
        if (Nx1 > 0 && sigmaX1 > 0) {
            numeratorX += Px1 * Nx1 / sigmaX1;
            denominatorX += Nx1 / sigmaX1;
        }
        if (Nx2 > 0 && sigmaX2 > 0) {
            numeratorX += Px2 * Nx2 / sigmaX2;
            denominatorX += Nx2 / sigmaX2;
        }
        x = (denominatorX > 0) ? numeratorX / denominatorX : 0.0;
        
        // Paper's formula for Y
        double numeratorY = 0.0;
        double denominatorY = 0.0;
        
        if (Ny1 > 0 && sigmaY1 > 0) {
            numeratorY += Py1 * Ny1 / sigmaY1;
            denominatorY += Ny1 / sigmaY1;
        }
        if (Ny2 > 0 && sigmaY2 > 0) {
            numeratorY += Py2 * Ny2 / sigmaY2;
            denominatorY += Ny2 / sigmaY2;
        }
        y = (denominatorY > 0) ? numeratorY / denominatorY : 0.0;
    }
};

vector<EventData> readPositionData(const string& filename) {
    vector<EventData> events;
    ifstream infile(filename);
    
    if (!infile.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return events;
    }
    
    string line;
    while (getline(infile, line)) {
        istringstream iss(line);
        int eventID, nBottom, nTop, nLeft, nRight;
        
        if (!(iss >> eventID >> nBottom >> nTop >> nLeft >> nRight)) {
            cerr << "Error reading event header for line: " << line << endl;
            continue;
        }
        
        EventData event(eventID);
        
        for (int i = 0; i < nBottom; i++) {
            double pos; iss >> pos;
            event.addXBottom(pos);
        }
        for (int i = 0; i < nTop; i++) {
            double pos; iss >> pos;
            event.addXTop(pos);
        }
        for (int i = 0; i < nLeft; i++) {
            double pos; iss >> pos;
            event.addYLeft(pos);
        }
        for (int i = 0; i < nRight; i++) {
            double pos; iss >> pos;
            event.addYRight(pos);
        }
        
        events.push_back(event);
    }
    
    infile.close();
    cout << "Read " << events.size() << " events from " << filename << endl;
    return events;
}

void createPlots(const vector<double>& x_vals, const vector<double>& y_vals) {
    // Set style
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    gStyle->SetPalette(kRainBow);
    
    // Create canvas
    TCanvas* c1 = new TCanvas("c1", "Paper Formula Reconstruction", 1200, 400);
    c1->Divide(3, 1);
    
    // Plot 1: X distribution
    c1->cd(1);
    TH1D* h_x = new TH1D("h_x", "X Position (Paper Formula);X (mm);Counts", 
                         100, -50, 50);
    for (double x : x_vals) h_x->Fill(x);
    h_x->SetFillColor(kBlue);
    h_x->SetFillStyle(3001);
    h_x->Draw();
    
    // Fit Gaussian to X distribution
    TF1* gauss_x = new TF1("gauss_x", "gaus", -50, 50);
    gauss_x->SetParameters(h_x->GetMaximum(), h_x->GetMean(), h_x->GetRMS());
    h_x->Fit(gauss_x, "RQ");
    gauss_x->SetLineColor(kRed);
    gauss_x->Draw("same");
    
    // Plot 2: Y distribution
    c1->cd(2);
    TH1D* h_y = new TH1D("h_y", "Y Position (Paper Formula);Y (mm);Counts", 
                         100, -50, 50);
    for (double y : y_vals) h_y->Fill(y);
    h_y->SetFillColor(kGreen);
    h_y->SetFillStyle(3001);
    h_y->Draw();
    
    // Fit Gaussian to Y distribution
    TF1* gauss_y = new TF1("gauss_y", "gaus", -50, 50);
    gauss_y->SetParameters(h_y->GetMaximum(), h_y->GetMean(), h_y->GetRMS());
    h_y->Fit(gauss_y, "RQ");
    gauss_y->SetLineColor(kRed);
    gauss_y->Draw("same");
    
    // Plot 3: 2D XY distribution
    c1->cd(3);
    TH2D* h_xy = new TH2D("h_xy", "2D Position (Paper Formula);X (mm);Y (mm)", 
                          100, -50, 50, 100, -50, 50);
    for (size_t i = 0; i < x_vals.size(); i++) {
        h_xy->Fill(x_vals[i], y_vals[i]);
    }
    h_xy->Draw("COLZ");
    
    // Save canvas
    c1->SaveAs("paper_formula_reconstruction.png");
    
    // Create another canvas for detailed 2D plot
    TCanvas* c2 = new TCanvas("c2", "Detailed 2D Reconstruction", 800, 800);
    h_xy->Draw("COLZ");
    
    // Add contour lines
    h_xy->SetContour(20);
    h_xy->Draw("CONTZ SAME");
    
    // Create profile
    TProfile2D* prof_xy = new TProfile2D("prof_xy", "Profile XY;X (mm);Y (mm)", 
                                         50, -50, 50, 50, -50, 50);
    for (size_t i = 0; i < x_vals.size(); i++) {
        prof_xy->Fill(x_vals[i], y_vals[i], 1.0);
    }
    prof_xy->SetMarkerStyle(20);
    prof_xy->SetMarkerSize(0.5);
    prof_xy->SetMarkerColor(kBlack);
    prof_xy->Draw("P SAME");
    
    // Add statistics box
    TPaveText* stats = new TPaveText(0.15, 0.85, 0.45, 0.95, "NDC");
    stats->SetFillColor(0);
    stats->SetBorderSize(1);
    stats->AddText(Form("Events: %lu", x_vals.size()));
    stats->AddText(Form("X mean: %.2f #pm %.2f mm", h_x->GetMean(), h_x->GetRMS()));
    stats->AddText(Form("Y mean: %.2f #pm %.2f mm", h_y->GetMean(), h_y->GetRMS()));
    stats->Draw();
    
    c2->SaveAs("paper_formula_2d_detailed.png");
    
    cout << "\n=== Reconstruction Results ===" << endl;
    cout << "Total reconstructed events: " << x_vals.size() << endl;
    cout << "X position: mean = " << h_x->GetMean() << " mm, sigma = " << h_x->GetRMS() << " mm" << endl;
    cout << "Y position: mean = " << h_y->GetMean() << " mm, sigma = " << h_y->GetRMS() << " mm" << endl;
    cout << "\nPlots saved as:" << endl;
    cout << "1. paper_formula_reconstruction.png - 1D distributions" << endl;
    cout << "2. paper_formula_2d_detailed.png - 2D detailed plot" << endl;
}

int main() {
    // Read data
    string filename = "../buildPosRec/position_data.txt";
    vector<EventData> events = readPositionData(filename);
    
    if (events.empty()) {
        cerr << "No events to process!" << endl;
        return 1;
    }
    
    // Store reconstructed positions
    vector<double> x_positions, y_positions;
    vector<int> eventIDs;
    vector<double> N_totals;
    
    // Create ROOT file for storage
    TFile* outputFile = new TFile("weighted_reconstruction.root", "RECREATE");
    TTree* tree = new TTree("reco", "Weighted Reconstruction Results");
    
    int eventID;
    double x_weighted, y_weighted;
    double N_total;
    tree->Branch("eventID", &eventID);
    tree->Branch("x_weighted", &x_weighted);
    tree->Branch("y_weighted", &y_weighted);
    tree->Branch("N_total", &N_total);
    
    // Process events
    int processed = 0;
    for (auto& event : events) {
        eventID = event.eventID;
        
        // Calculate total photons
        N_total = event.xBottom.size() + event.xTop.size() + 
                  event.yLeft.size() + event.yRight.size();
        
        // Skip events with too few photons
        if (N_total < 1) continue;
        
        // Reconstruct using paper's formula
        event.reconstructWeighted(x_weighted, y_weighted);
        
        // Store results
        x_positions.push_back(x_weighted);
        y_positions.push_back(y_weighted);
        eventIDs.push_back(eventID);
        N_totals.push_back(N_total);
        
        // Fill tree
        tree->Fill();
        processed++;
        
        // Print progress
        if (processed % 100 == 0) {
            cout << "Processed " << processed << " events" << endl;
        }
    }
    
    // Write tree to file
    tree->Write();
    
    // Create and save plots
    if (!x_positions.empty()) {
        createPlots(x_positions, y_positions);
        
        // Save histograms to file
        TH1D* h_x_file = new TH1D("h_x_weighted", "X Position;X (mm);Counts", 
                                  100, -50, 50);
        TH1D* h_y_file = new TH1D("h_y_weighted", "Y Position;Y (mm);Counts", 
                                  100, -50, 50);
        TH2D* h_xy_file = new TH2D("h_xy_weighted", "2D Position;X (mm);Y (mm)", 
                                   100, -50, 50, 100, -50, 50);
        
        for (size_t i = 0; i < x_positions.size(); i++) {
            h_x_file->Fill(x_positions[i]);
            h_y_file->Fill(y_positions[i]);
            h_xy_file->Fill(x_positions[i], y_positions[i]);
        }
        
        h_x_file->Write();
        h_y_file->Write();
        h_xy_file->Write();
    }
    
    outputFile->Close();
    delete outputFile;
    
    cout << "\n=== Analysis Complete ===" << endl;
    cout << "Processed " << processed << " out of " << events.size() << " events" << endl;
    cout << "Results saved to weighted_reconstruction.root" << endl;
    cout << "Plots saved as PNG files" << endl;
    
    // Open ROOT file to view results
    cout << "\nTo view results in ROOT:" << endl;
    cout << "  root -l weighted_reconstruction.root" << endl;
    cout << "  TBrowser tb" << endl;
    
    return 0;
}
