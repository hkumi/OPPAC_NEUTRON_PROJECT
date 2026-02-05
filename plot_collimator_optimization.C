#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"

void plot_collimator_optimization() {
    std::cout << "\n==========================================" << std::endl;
    std::cout << "COLLIMATOR OPTIMIZATION - CORRECTED ANALYSIS" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // =========================================================
    // YOUR DATA - ALL 10 POINTS
    // =========================================================
    
    std::vector<double> colLength = {3.0, 5.0, 7.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 30.0};
    std::vector<int> Ntotal(10, 10000);
    
    std::vector<int> Ndetected = { 
        32831635,  // 3mm - REAL
        28154036,  // 5mm - REAL  
        25343466,  // 7mm - REAL
        22036613,  // 10mm - REAL
        20643024,  // 12mm - REAL 
        19471565,  // 14mm - Real
        18519058,  // 16mm - real
        17704701,  // 18mm - real
        17072546,  // 20mm - real
        14816454    // 30mm - real
    };
    
    std::vector<double> sigmaX = { 
        0.269,  // 3mm - REAL
        0.286,  // 5mm - REAL
        0.305,  // 7mm - REAL
        0.317,  // 10mm - REAL
        0.328,  // 12mm - REAL 
        0.329,  // 14mm - real
        0.340,  // 16mm - real
        0.343,  // 18mm - real
        0.349,  // 20mm - real
        0.368   // 30mm - real
    };
    
    std::vector<double> sigmaY = { 
        0.271,  // 3mm - REAL
        0.280,  // 5mm - REAL
        0.302,  // 7mm - REAL
        0.318,  // 10mm - REAL
        0.331,  // 12mm - REAL 
        0.336,  // 14mm - real
        0.341,  // 16mm - real
        0.339,  // 18mm - real
        0.346,  // 20mm - real
        0.370   // 30mm - real
    };
    
    // =========================================================
    // CORRECTED CALCULATIONS
    // =========================================================
    
    int n = colLength.size();
    std::vector<double> avgPhotons(n), resolution(n), FoM(n);
    std::vector<double> zeroErr(n, 0.0);
    
    std::cout << "\nANALYZING ALL " << n << " DATA POINTS (CORRECTED):" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    for (int i = 0; i < n; i++) {
        // CORRECT: Average photons per event (NOT efficiency)
        avgPhotons[i] = (double)Ndetected[i] / Ntotal[i];
        
        // Mean spatial resolution
        resolution[i] = 0.5 * (sigmaX[i] + sigmaY[i]);
        
        // Figure of Merit: Photons per event / Resolution²
        // Higher photons = better, lower resolution = better
        FoM[i] = avgPhotons[i] / (resolution[i] * resolution[i]);
        
        std::cout << "\nL = " << colLength[i] << " mm:" << std::endl;
        std::cout << "  Total photons detected: " << Ndetected[i] << std::endl;
        std::cout << "  Average photons/event: " << avgPhotons[i] << std::endl;
        std::cout << "  Resolution: " << resolution[i] << " mm (σ)" << std::endl;
        std::cout << "  FWHM: " << (2.355 * resolution[i]) << " mm" << std::endl;
        std::cout << "  Figure of Merit (Photons/σ²): " << FoM[i] << std::endl;
    }
    
    // =========================================================
    // FIND OPTIMAL FROM ALL DATA
    // =========================================================
    
    double maxFoM = 0;
    double optimalLength = 0;
    int optimalIndex = 0;
    
    for (int i = 0; i < n; i++) {
        if (FoM[i] > maxFoM) {
            maxFoM = FoM[i];
            optimalLength = colLength[i];
            optimalIndex = i;
        }
    }
    
    std::cout << "\n==========================================" << std::endl;
    std::cout << "OPTIMIZATION RESULTS (CORRECTED)" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "Optimal collimator length: " << optimalLength << " mm" << std::endl;
    std::cout << "Maximum Figure of Merit: " << maxFoM << std::endl;
    std::cout << "Photons per event at optimum: " << avgPhotons[optimalIndex] << std::endl;
    std::cout << "Resolution at optimum: " << resolution[optimalIndex] << " mm (σ)" << std::endl;
    std::cout << "FWHM at optimum: " << (2.355 * resolution[optimalIndex]) << " mm" << std::endl;
    std::cout << "Improvement vs 1mm binning: " << (1.0/resolution[optimalIndex]) << "×" << std::endl;
    
    // =========================================================
    // CREATE PLOTS WITH ALL DATA
    // =========================================================
    
    TCanvas *c1 = new TCanvas("c1", "Collimator Optimization", 1200, 800);
    c1->Divide(2, 2);
    
    // --- Plot 1: Photons per Event vs Length ---
    c1->cd(1);
    TGraph *grPhotons = new TGraph(n, &colLength[0], &avgPhotons[0]);
    grPhotons->SetTitle("Photon Statistics;Collimator Length (mm);Average Photons per Event");
    grPhotons->SetMarkerStyle(20);
    grPhotons->SetMarkerColor(kBlue);
    grPhotons->SetMarkerSize(1.5);
    grPhotons->SetLineColor(kBlue);
    grPhotons->SetLineWidth(2);
    grPhotons->Draw("APL");
    
    // Fit with exponential decay
    TF1 *fitPhotons = new TF1("fitPhotons", "[0]*exp(-x/[1]) + [2]", colLength[0], colLength.back());
    fitPhotons->SetParameters(3000, 10.0, 100);
    fitPhotons->SetLineColor(kRed);
    fitPhotons->SetLineStyle(2);
    grPhotons->Fit(fitPhotons, "RQ");
    
    // --- Plot 2: Resolution vs Length ---
    c1->cd(2);
    TGraph *grRes = new TGraph(n, &colLength[0], &resolution[0]);
    grRes->SetTitle("Spatial Resolution;Collimator Length (mm);Resolution #sigma (mm)");
    grRes->SetMarkerStyle(21);
    grRes->SetMarkerColor(kRed);
    grRes->SetMarkerSize(1.5);
    grRes->SetLineColor(kRed);
    grRes->SetLineWidth(2);
    grRes->Draw("APL");
    
    /*// Fit with linear or power law
    TF1 *fitRes = new TF1("fitRes", "[0] + [1]*x", colLength[0], colLength.back());
    fitRes->SetParameters(0.26, 0.005);
    fitRes->SetLineColor(kBlue);
    fitRes->SetLineStyle(2);
    grRes->Fit(fitRes, "RQ");*/
    
    // --- Plot 3: Figure of Merit vs Length ---
    c1->cd(3);
    TGraph *grFoMplot = new TGraph(n, &colLength[0], &FoM[0]);
    grFoMplot->SetTitle("Figure of Merit;Collimator Length (mm);FoM = Photons/Resolution^{2}");
    grFoMplot->SetMarkerStyle(22);
    grFoMplot->SetMarkerColor(kGreen+2);
    grFoMplot->SetMarkerSize(1.5);
    grFoMplot->SetLineColor(kGreen+2);
    grFoMplot->SetLineWidth(2);
    grFoMplot->Draw("APL");
    
    // Mark optimal point
    TGraph *optPoint = new TGraph(1);
    optPoint->SetPoint(0, optimalLength, maxFoM);
    optPoint->SetMarkerStyle(29);
    optPoint->SetMarkerSize(3.5);
    optPoint->SetMarkerColor(kMagenta);
    optPoint->Draw("P SAME");
    
    // Draw text with optimal value
    TLatex *tex = new TLatex();
    tex->SetTextSize(0.06);
    tex->SetTextColor(kMagenta);
    tex->DrawLatex(optimalLength + 2.0, maxFoM * 0.85, Form("Optimal: %.1f mm", optimalLength));
    
    // --- Plot 4: Normalized comparison ---
    c1->cd(4);
    
    // Normalize for comparison
    std::vector<double> photonsNorm(n), resNorm(n), fomNorm(n);
    double maxPhotons = *std::max_element(avgPhotons.begin(), avgPhotons.end());
    double maxRes = *std::max_element(resolution.begin(), resolution.end());
    double maxFom = *std::max_element(FoM.begin(), FoM.end());
    
    for (int i = 0; i < n; i++) {
        photonsNorm[i] = avgPhotons[i] / maxPhotons;
        resNorm[i] = resolution[i] / maxRes;
        fomNorm[i] = FoM[i] / maxFom;
    }
    
    TGraph *grPhotonsNorm = new TGraph(n, &colLength[0], &photonsNorm[0]);
    TGraph *grResNorm = new TGraph(n, &colLength[0], &resNorm[0]);
    TGraph *grFoMNorm = new TGraph(n, &colLength[0], &fomNorm[0]);
    
    grPhotonsNorm->SetTitle("Normalized Comparison;Collimator Length (mm);Normalized Value");
    grPhotonsNorm->SetMarkerStyle(20);
    grPhotonsNorm->SetMarkerColor(kBlue);
    grPhotonsNorm->SetLineColor(kBlue);
    grPhotonsNorm->SetLineWidth(2);
    
    grResNorm->SetMarkerStyle(21);
    grResNorm->SetMarkerColor(kRed);
    grResNorm->SetLineColor(kRed);
    grResNorm->SetLineWidth(2);
    
    grFoMNorm->SetMarkerStyle(22);
    grFoMNorm->SetMarkerColor(kGreen+2);
    grFoMNorm->SetLineColor(kGreen+2);
    grFoMNorm->SetLineWidth(2);
    
    grPhotonsNorm->Draw("APL");
    grResNorm->Draw("PL SAME");
    grFoMNorm->Draw("PL SAME");
    
    // Add vertical line at optimal length
    TLine *optLine = new TLine(optimalLength, 0, optimalLength, 1.1);
    optLine->SetLineColor(kMagenta);
    optLine->SetLineStyle(2);
    optLine->SetLineWidth(2);
    optLine->Draw();
    
    // Legend
    TLegend *leg = new TLegend(0.15, 0.70, 0.45, 0.90);
    leg->AddEntry(grPhotonsNorm, "Photons per Event", "lp");
    leg->AddEntry(grResNorm, "Resolution", "lp");
    leg->AddEntry(grFoMNorm, "Figure of Merit", "lp");
    leg->AddEntry(optLine, Form("Optimal = %.1f mm", optimalLength), "l");
    leg->Draw();
    
    // =========================================================
    // SAVE RESULTS
    // =========================================================
    
    c1->SaveAs("collimator_optimization_CORRECTED.png");
    std::cout << "\nPlot saved as: collimator_optimization_CORRECTED.png" << std::endl;
    
    // Save all data to text file
    std::ofstream outfile("collimator_results_CORRECTED.txt");
    outfile << "# Collimator Optimization Results (Corrected)" << std::endl;
    outfile << "# Length(mm) TotalEvents TotalPhotons AvgPhotonsPerEvent sigmaX(mm) sigmaY(mm) Resolution(mm) FoM" << std::endl;
    
    for (int i = 0; i < n; i++) {
        outfile << colLength[i] << " " 
                << Ntotal[i] << " " 
                << Ndetected[i] << " " 
                << avgPhotons[i] << " " 
                << sigmaX[i] << " " 
                << sigmaY[i] << " " 
                << resolution[i] << " " 
                << FoM[i] << std::endl;
    }
    outfile.close();
    
    std::cout << "Data saved to: collimator_results_CORRECTED.txt" << std::endl;
    
    // =========================================================
    // CALCULATE TRADEOFF ANALYSIS
    // =========================================================
    
    std::cout << "\n==========================================" << std::endl;
    std::cout << "TRADEOFF ANALYSIS" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // Calculate relative improvements
    double bestPhotons = avgPhotons[0];  // 3mm has most photons
    double bestRes = *std::min_element(resolution.begin(), resolution.end());  // Best resolution
    
    for (int i = 0; i < n; i++) {
        double photonLoss = 100.0 * (bestPhotons - avgPhotons[i]) / bestPhotons;
        double resImprovement = 100.0 * (resolution[i] - bestRes) / bestRes;
        
        std::cout << colLength[i] << "mm: " 
                  << photonLoss << "% fewer photons, "
                  << resImprovement << "% worse resolution than best" << std::endl;
    }
    
    std::cout << "\n==========================================" << std::endl;
    std::cout << "RECOMMENDATION:" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    if (optimalLength <= 5.0) {
        std::cout << "Short collimators (<5mm) give best FoM but may have:" << std::endl;
        std::cout << "- Worse angular selectivity" << std::endl;
        std::cout << "- More crosstalk between channels" << std::endl;
        std::cout << "- Practical construction challenges" << std::endl;
        std::cout << "\nConsider 5-7mm as a practical compromise." << std::endl;
    } else {
        std::cout << "Optimal length " << optimalLength << "mm provides good balance." << std::endl;
    }
    
    std::cout << "\nFor publication, complete simulations for 14-30mm." << std::endl;
    std::cout << "==========================================" << std::endl;
}

// =========================================================
// MAIN FUNCTION
// =========================================================

int main() {
    plot_collimator_optimization();
    return 0;
}
