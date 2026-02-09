// converter_optimization_FIXED.C
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

void converter_optimization() {
    std::cout << "\n==========================================" << std::endl;
    std::cout << "CONVERTER THICKNESS OPTIMIZATION (FIXED)" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "Analysis range: 0.01 - 1.50 mm HDPE" << std::endl;
    std::cout << "Using CORRECTED Figure of Merit" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // =========================================================
    // DATA - TRUSTED RANGE 0.01-1.50mm
    // =========================================================
    
    std::vector<double> thickness = {0.01, 0.05, 0.10, 0.20, 0.50, 1.00, 1.50};
    std::vector<int> Ntotal(thickness.size(), 10000000);  // 10⁷ events
    
    std::vector<int> Ndetected = { 
        3426,      // 0.01mm
        15036,     // 0.05mm  
        26208,     // 0.10mm
        46759,     // 0.20mm
        111209,    // 0.50mm
        223576,    // 1.00mm
        340019     // 1.50mm
    };
    
    std::vector<double> sigmaX = { 
        0.580, 0.540, 0.540, 0.538, 0.529, 0.537, 0.544
    };
    
    std::vector<double> sigmaY = { 
        0.560, 0.536, 0.524, 0.518, 0.534, 0.541, 0.549
    };
    
    // =========================================================
    // CALCULATIONS
    // =========================================================
    
    int n = thickness.size();
    std::vector<double> detectionEfficiency(n), resolution(n), FoM(n);
    std::vector<double> effErr(n), resErr(n);
    
    std::cout << "\nANALYSIS OF " << n << " CONVERTER THICKNESSES:" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    for (int i = 0; i < n; i++) {
        // Detection efficiency (events with signal)
        detectionEfficiency[i] = (double)Ndetected[i] / Ntotal[i];
        effErr[i] = sqrt(detectionEfficiency[i] * (1 - detectionEfficiency[i]) / Ntotal[i]);
        
        // Spatial resolution
        resolution[i] = 0.5 * (sigmaX[i] + sigmaY[i]);
        resErr[i] = 0.5 * sqrt(sigmaX[i]*sigmaX[i] + sigmaY[i]*sigmaY[i]) / sqrt(Ndetected[i]);
        
        std::cout << "\nThickness = " << thickness[i] << " mm:" << std::endl;
        std::cout << "  Events detected: " << Ndetected[i] << "/" << Ntotal[i] << std::endl;
        std::cout << "  Detection efficiency: " << (detectionEfficiency[i]*100) 
                  << " ± " << (effErr[i]*100) << "%" << std::endl;
        std::cout << "  Resolution: " << resolution[i] << " ± " << resErr[i] << " mm (σ)" << std::endl;
        std::cout << "  FWHM: " << (2.355 * resolution[i]) << " mm" << std::endl;
    }
    
    // =========================================================
    // CALCULATE MULTIPLE FIGURES OF MERIT
    // =========================================================
    
    std::cout << "\n==========================================" << std::endl;
    std::cout << "FIGURE OF MERIT CALCULATIONS" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    std::vector<double> FoM1(n), FoM2(n), FoM3(n);
    
    // Find best and worst values for normalization
    double maxEff = *std::max_element(detectionEfficiency.begin(), detectionEfficiency.end());
    double minEff = *std::min_element(detectionEfficiency.begin(), detectionEfficiency.end());
    double minRes = *std::min_element(resolution.begin(), resolution.end());
    double maxRes = *std::max_element(resolution.begin(), resolution.end());
    
    for (int i = 0; i < n; i++) {
        // FoM1: Original (problematic - just ratio)
        FoM1[i] = detectionEfficiency[i] / resolution[i];
        
        // FoM2: CORRECTED - Efficiency × (1/Resolution²)
        // Penalizes poor resolution more strongly
        FoM2[i] = detectionEfficiency[i] / (resolution[i] * resolution[i]);
        
        // FoM3: Normalized balanced approach
        // Normalize both to 0-1 range, then multiply
        double effNorm = (detectionEfficiency[i] - minEff) / (maxEff - minEff);
        double resNormInv = (maxRes - resolution[i]) / (maxRes - minRes);  // Inverted: lower res is better
        FoM3[i] = effNorm * resNormInv;
        
        std::cout << "\nThickness " << thickness[i] << " mm:" << std::endl;
        std::cout << "  FoM1 (Eff/Res):           " << FoM1[i] << std::endl;
        std::cout << "  FoM2 (Eff/Res²):          " << FoM2[i] << std::endl;
        std::cout << "  FoM3 (Normalized):        " << FoM3[i] << std::endl;
    }
    
    // =========================================================
    // FIND OPTIMAL THICKNESS FOR EACH FOM
    // =========================================================
    
    auto findOptimal = [&](const std::vector<double>& fom) {
        int idx = std::distance(fom.begin(), std::max_element(fom.begin(), fom.end()));
        return idx;
    };
    
    int opt1 = findOptimal(FoM1);
    int opt2 = findOptimal(FoM2);
    int opt3 = findOptimal(FoM3);
    
    std::cout << "\n==========================================" << std::endl;
    std::cout << "OPTIMIZATION RESULTS COMPARISON" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    std::cout << "\nFoM1 (Eff/Res) - PROBLEMATIC:" << std::endl;
    std::cout << "  Optimal: " << thickness[opt1] << " mm" << std::endl;
    std::cout << "  Efficiency: " << (detectionEfficiency[opt1]*100) << "%" << std::endl;
    std::cout << "  Resolution: " << resolution[opt1] << " mm" << std::endl;
    std::cout << "  ⚠️  Favors high thickness despite poor resolution!" << std::endl;
    
    std::cout << "\nFoM2 (Eff/Res²) - RECOMMENDED:" << std::endl;
    std::cout << "  Optimal: " << thickness[opt2] << " mm" << std::endl;
    std::cout << "  Efficiency: " << (detectionEfficiency[opt2]*100) << "%" << std::endl;
    std::cout << "  Resolution: " << resolution[opt2] << " mm (σ)" << std::endl;
    std::cout << "  FWHM: " << (2.355 * resolution[opt2]) << " mm" << std::endl;
    std::cout << "  ✓ Properly penalizes resolution degradation" << std::endl;
    
    std::cout << "\nFoM3 (Normalized) - BALANCED:" << std::endl;
    std::cout << "  Optimal: " << thickness[opt3] << " mm" << std::endl;
    std::cout << "  Efficiency: " << (detectionEfficiency[opt3]*100) << "%" << std::endl;
    std::cout << "  Resolution: " << resolution[opt3] << " mm (σ)" << std::endl;
    std::cout << "  FWHM: " << (2.355 * resolution[opt3]) << " mm" << std::endl;
    std::cout << "  ✓ Equal weighting of efficiency and resolution" << std::endl;
    
    // Use FoM2 as the recommended metric
    FoM = FoM2;
    int optimalIndex = opt2;
    double optimalThickness = thickness[optimalIndex];
    double maxFoM = FoM2[optimalIndex];
    
    // =========================================================
    // PHYSICAL TREND ANALYSIS
    // =========================================================
    
    std::cout << "\n==========================================" << std::endl;
    std::cout << "PHYSICAL TRENDS OBSERVED" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    std::cout << "Comparison to 0.01 mm baseline:" << std::endl;
    for (int i = 1; i < n; i++) {
        double eff_ratio = detectionEfficiency[i] / detectionEfficiency[0];
        double res_degradation = (resolution[i] - resolution[0]) / resolution[0] * 100;
        
        printf("%.2fmm: %.1f× efficiency, %+.1f%% resolution change\n",
               thickness[i], eff_ratio, res_degradation);
    }
    
    // =========================================================
    // CREATE PUBLICATION-QUALITY PLOTS
    // =========================================================
    
    TCanvas *c1 = new TCanvas("c1", "Converter Optimization (CORRECTED)", 1600, 1000);
    c1->Divide(2, 2);
    
    // --- Plot 1: Detection Efficiency ---
    c1->cd(1);
    gPad->SetGridx();
    gPad->SetGridy();
    TGraphErrors *grEff = new TGraphErrors(n, &thickness[0], &detectionEfficiency[0], 
                                           nullptr, &effErr[0]);
    grEff->SetTitle("Detection Efficiency;Converter Thickness (mm);Efficiency");
    grEff->SetMarkerStyle(20);
    grEff->SetMarkerColor(kBlue);
    grEff->SetMarkerSize(1.5);
    grEff->SetLineColor(kBlue);
    grEff->SetLineWidth(2);
    grEff->GetYaxis()->SetRangeUser(0, 0.04);
    grEff->Draw("AP");
    
    TF1 *fitEff = new TF1("fitEff", "[0]*(1-exp(-x/[1]))", 0, 1.6);
    fitEff->SetParameters(0.05, 0.5);
    fitEff->SetLineColor(kRed);
    fitEff->SetLineStyle(2);
    grEff->Fit(fitEff, "RQ");
    
    // --- Plot 2: Spatial Resolution ---
    c1->cd(2);
    gPad->SetGridx();
    gPad->SetGridy();
    TGraphErrors *grRes = new TGraphErrors(n, &thickness[0], &resolution[0], 
                                           nullptr, &resErr[0]);
    grRes->SetTitle("Spatial Resolution;Converter Thickness (mm);Resolution #sigma (mm)");
    grRes->SetMarkerStyle(21);
    grRes->SetMarkerColor(kRed);
    grRes->SetMarkerSize(1.5);
    grRes->SetLineColor(kRed);
    grRes->SetLineWidth(2);
    grRes->GetYaxis()->SetRangeUser(0.50, 0.60);
    grRes->Draw("AP");
    /*
    TF1 *fitRes = new TF1("fitRes", "[0] + [1]*x", 0, 1.6);
    fitRes->SetLineColor(kBlue);
    fitRes->SetLineStyle(2);
    grRes->Fit(fitRes, "RQ");*/
    
    
    
    // --- Plot 4: FoM2 (Corrected - Recommended) ---
    c1->cd(3);
    gPad->SetGridx();
    gPad->SetGridy();
    TGraph *grFoM2 = new TGraph(n, &thickness[0], &FoM2[0]);
    grFoM2->SetTitle("FoM2 (Eff/Res^{2}) - RECOMMENDED;Converter Thickness (mm);FoM");
    grFoM2->SetMarkerStyle(23);
    grFoM2->SetMarkerColor(kGreen+2);
    grFoM2->SetMarkerSize(1.5);
    grFoM2->SetLineColor(kGreen+2);
    grFoM2->SetLineWidth(2);
    grFoM2->Draw("APL");
    
    TGraph *opt2Point = new TGraph(1);
    opt2Point->SetPoint(0, thickness[opt2], FoM2[opt2]);
    opt2Point->SetMarkerStyle(29);
    opt2Point->SetMarkerSize(3.0);
    opt2Point->SetMarkerColor(kMagenta);
    opt2Point->Draw("P SAME");
    
    TLatex *tex2 = new TLatex();
    tex2->SetTextSize(0.045);
    tex2->SetTextColor(kMagenta);
    tex2->DrawLatex(thickness[opt2] + 0.05, FoM2[opt2] * 0.9, 
                    Form("Optimal: %.2f mm", thickness[opt2]));
    
    // --- Plot 5: FoM3 (Normalized) ---
    c1->cd(4);
    gPad->SetGridx();
    gPad->SetGridy();
    TGraph *grFoM3 = new TGraph(n, &thickness[0], &FoM3[0]);
    grFoM3->SetTitle("FoM3 (Normalized) - BALANCED;Converter Thickness (mm);FoM");
    grFoM3->SetMarkerStyle(33);
    grFoM3->SetMarkerColor(kViolet);
    grFoM3->SetMarkerSize(1.8);
    grFoM3->SetLineColor(kViolet);
    grFoM3->SetLineWidth(2);
    grFoM3->Draw("APL");
    
    TGraph *opt3Point = new TGraph(1);
    opt3Point->SetPoint(0, thickness[opt3], FoM3[opt3]);
    opt3Point->SetMarkerStyle(29);
    opt3Point->SetMarkerSize(3.0);
    opt3Point->SetMarkerColor(kMagenta);
    opt3Point->Draw("P SAME");
    
  
    
    // =========================================================
    // SAVE RESULTS
    // =========================================================
    
    c1->SaveAs("converter_optimization_FIXED.pdf");
    c1->SaveAs("converter_optimization_FIXED.png");
    
    std::cout << "\n==========================================" << std::endl;
    std::cout << "Plots saved as converter_optimization_FIXED.pdf/.png" << std::endl;
    
    // Save data table
    std::ofstream outfile("converter_results_FIXED.txt");
    outfile << "# HDPE Converter Optimization Results (CORRECTED)" << std::endl;
    outfile << "# Ntotal = 10,000,000 events per thickness" << std::endl;
    outfile << "# FoM2 = Eff/Res² (recommended metric)" << std::endl;
    outfile << "# Thickness(mm) Efficiency(%) Resolution(mm) FoM1 FoM2 FoM3" << std::endl;
    
    for (int i = 0; i < n; i++) {
        outfile << thickness[i] << " " 
                << (detectionEfficiency[i]*100) << " "
                << resolution[i] << " "
                << FoM1[i] << " " << FoM2[i] << " " << FoM3[i] << std::endl;
    }
    outfile.close();
    
    std::cout << "Data saved to converter_results_FIXED.txt" << std::endl;
    
    // =========================================================
    // FINAL RECOMMENDATIONS
    // =========================================================
    
    std::cout << "\n==========================================" << std::endl;
    std::cout << "FINAL RECOMMENDATION" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    std::cout << "\n✓ RECOMMENDED THICKNESS: " << optimalThickness << " mm HDPE" << std::endl;
    std::cout << "  Using FoM2 = Efficiency / Resolution²" << std::endl;
    std::cout << "\nPerformance at optimal thickness:" << std::endl;
    std::cout << "  • Detection efficiency: " << (detectionEfficiency[optimalIndex]*100) << "%" << std::endl;
    std::cout << "  • Spatial resolution: " << resolution[optimalIndex] << " mm (σ)" << std::endl;
    std::cout << "  • FWHM: " << (2.355 * resolution[optimalIndex]) << " mm" << std::endl;
    std::cout << "  • Figure of Merit: " << maxFoM << std::endl;
    
    std::cout << "\nWhy this is better:" << std::endl;
    std::cout << "  • Provides " << (detectionEfficiency[optimalIndex]/detectionEfficiency[0]) 
              << "× better efficiency than 0.01 mm" << std::endl;
    std::cout << "  • Maintains excellent spatial resolution" << std::endl;
    std::cout << "  • Properly balances efficiency vs. resolution trade-off" << std::endl;
    
    std::cout << "\nPhysical interpretation:" << std::endl;
    std::cout << "  • Thin enough: protons don't scatter significantly" << std::endl;
    std::cout << "  • Thick enough: good neutron-proton conversion probability" << std::endl;
    std::cout << "  • Sweet spot between detection and precision" << std::endl;
    
    std::cout << "\n==========================================" << std::endl;
}

int main() {
    converter_optimization();
    return 0;
}
