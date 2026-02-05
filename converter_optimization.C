// converter_optimization_WITH_CONSTRAINTS.C
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
    std::cout << "CONVERTER OPTIMIZATION WITH CONSTRAINTS" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // =========================================================
    // DATA
    // =========================================================
    
    std::vector<double> thickness = {0.01, 0.05, 0.10, 0.20, 0.50, 1.00, 1.50};
    std::vector<int> Ntotal(thickness.size(), 10000000);
    
    std::vector<int> Ndetected = { 
        3426, 15036, 26208, 46759, 111209, 223576, 340019
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
    std::vector<double> efficiency(n), resolution(n);
    std::vector<double> effErr(n), resErr(n);
    
    for (int i = 0; i < n; i++) {
        efficiency[i] = (double)Ndetected[i] / Ntotal[i];
        effErr[i] = sqrt(efficiency[i] * (1 - efficiency[i]) / Ntotal[i]);
        resolution[i] = 0.5 * (sigmaX[i] + sigmaY[i]);
        resErr[i] = 0.5 * sqrt(sigmaX[i]*sigmaX[i] + sigmaY[i]*sigmaY[i]) / sqrt(Ndetected[i]);
    }
    
    std::cout << "\nDATA SUMMARY:" << std::endl;
    std::cout << "==========================================" << std::endl;
    for (int i = 0; i < n; i++) {
        printf("%.2f mm: Eff=%.3f%%, Res=%.3f mm (FWHM=%.2f mm)\n",
               thickness[i], efficiency[i]*100, resolution[i], resolution[i]*2.355);
    }
    
    // =========================================================
    // MULTIPLE OPTIMIZATION STRATEGIES
    // =========================================================
    
    std::cout << "\n==========================================" << std::endl;
    std::cout << "OPTIMIZATION STRATEGIES" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // Strategy 1: Resolution Constraint
    std::cout << "\n--- STRATEGY 1: Resolution Constraint ---" << std::endl;
    double maxAllowedResolution = 0.530;  // Max acceptable resolution (mm)
    std::cout << "Constraint: Resolution must be ≤ " << maxAllowedResolution << " mm" << std::endl;
    
    int opt_strategy1 = -1;
    double maxEff_constrained = 0;
    for (int i = 0; i < n; i++) {
        if (resolution[i] <= maxAllowedResolution && efficiency[i] > maxEff_constrained) {
            maxEff_constrained = efficiency[i];
            opt_strategy1 = i;
        }
    }
    
    if (opt_strategy1 >= 0) {
        std::cout << "Optimal: " << thickness[opt_strategy1] << " mm" << std::endl;
        std::cout << "  Efficiency: " << (efficiency[opt_strategy1]*100) << "%" << std::endl;
        std::cout << "  Resolution: " << resolution[opt_strategy1] << " mm" << std::endl;
    }
    
    // Strategy 2: Weighted FoM (user-defined weights)
    std::cout << "\n--- STRATEGY 2: Weighted FoM ---" << std::endl;
    double alpha = 5.0;  // Resolution weight (higher = more emphasis on resolution)
    std::cout << "Using FoM = Efficiency / Resolution^α with α = " << alpha << std::endl;
    
    std::vector<double> FoM_weighted(n);
    for (int i = 0; i < n; i++) {
        FoM_weighted[i] = efficiency[i] / pow(resolution[i], alpha);
    }
    
    int opt_strategy2 = std::distance(FoM_weighted.begin(), 
                                      std::max_element(FoM_weighted.begin(), FoM_weighted.end()));
    
    std::cout << "Optimal: " << thickness[opt_strategy2] << " mm" << std::endl;
    std::cout << "  Efficiency: " << (efficiency[opt_strategy2]*100) << "%" << std::endl;
    std::cout << "  Resolution: " << resolution[opt_strategy2] << " mm" << std::endl;
    std::cout << "  FoM: " << FoM_weighted[opt_strategy2] << std::endl;
    
    // Strategy 3: Normalized product with equal weights
    std::cout << "\n--- STRATEGY 3: Normalized Product ---" << std::endl;
    
    double maxEff = *std::max_element(efficiency.begin(), efficiency.end());
    double minRes = *std::min_element(resolution.begin(), resolution.end());
    double maxRes = *std::max_element(resolution.begin(), resolution.end());
    
    std::vector<double> FoM_normalized(n);
    for (int i = 0; i < n; i++) {
        double effNorm = efficiency[i] / maxEff;
        double resScore = 1.0 - (resolution[i] - minRes) / (maxRes - minRes);
        FoM_normalized[i] = effNorm * resScore;
    }
    
    int opt_strategy3 = std::distance(FoM_normalized.begin(), 
                                      std::max_element(FoM_normalized.begin(), FoM_normalized.end()));
    
    std::cout << "Optimal: " << thickness[opt_strategy3] << " mm" << std::endl;
    std::cout << "  Efficiency: " << (efficiency[opt_strategy3]*100) << "%" << std::endl;
    std::cout << "  Resolution: " << resolution[opt_strategy3] << " mm" << std::endl;
    std::cout << "  FoM: " << FoM_normalized[opt_strategy3] << std::endl;
    
    // Strategy 4: Knee point detection (diminishing returns)
    std::cout << "\n--- STRATEGY 4: Diminishing Returns Analysis ---" << std::endl;
    std::cout << "Finding where efficiency gain per resolution loss is maximized" << std::endl;
    
    std::vector<double> efficiency_gain_rate(n-1);
    std::vector<double> resolution_loss_rate(n-1);
    std::vector<double> trade_off_ratio(n-1);
    
    for (int i = 1; i < n; i++) {
        efficiency_gain_rate[i-1] = (efficiency[i] - efficiency[i-1]) / (thickness[i] - thickness[i-1]);
        resolution_loss_rate[i-1] = (resolution[i] - resolution[i-1]) / (thickness[i] - thickness[i-1]);
        
        // Avoid division by zero
        if (fabs(resolution_loss_rate[i-1]) > 1e-6) {
            trade_off_ratio[i-1] = efficiency_gain_rate[i-1] / fabs(resolution_loss_rate[i-1]);
        } else {
            trade_off_ratio[i-1] = efficiency_gain_rate[i-1] / 1e-6;
        }
        
        printf("  %.2f→%.2f mm: Eff gain rate=%.4f%%/mm, Res loss rate=%.4f mm/mm, Ratio=%.2f\n",
               thickness[i-1], thickness[i], 
               efficiency_gain_rate[i-1]*100,
               resolution_loss_rate[i-1],
               trade_off_ratio[i-1]);
    }
    
    int max_ratio_idx = std::distance(trade_off_ratio.begin(), 
                                      std::max_element(trade_off_ratio.begin(), trade_off_ratio.end()));
    int opt_strategy4 = max_ratio_idx;  // Index before the step with max ratio
    
    std::cout << "\nBest trade-off occurs at: " << thickness[opt_strategy4] << " mm" << std::endl;
    std::cout << "  Efficiency: " << (efficiency[opt_strategy4]*100) << "%" << std::endl;
    std::cout << "  Resolution: " << resolution[opt_strategy4] << " mm" << std::endl;
    
    // Strategy 5: Target efficiency threshold
    std::cout << "\n--- STRATEGY 5: Efficiency Threshold ---" << std::endl;
    double targetEfficiency = 0.005;  // 0.5% target
    std::cout << "Target: Achieve ≥" << (targetEfficiency*100) << "% efficiency with best resolution" << std::endl;
    
    int opt_strategy5 = -1;
    double bestRes_aboveThreshold = 1000;
    for (int i = 0; i < n; i++) {
        if (efficiency[i] >= targetEfficiency && resolution[i] < bestRes_aboveThreshold) {
            bestRes_aboveThreshold = resolution[i];
            opt_strategy5 = i;
        }
    }
    
    if (opt_strategy5 >= 0) {
        std::cout << "Optimal: " << thickness[opt_strategy5] << " mm" << std::endl;
        std::cout << "  Efficiency: " << (efficiency[opt_strategy5]*100) << "%" << std::endl;
        std::cout << "  Resolution: " << resolution[opt_strategy5] << " mm" << std::endl;
    }
    
    // =========================================================
    // SUMMARY AND RECOMMENDATION
    // =========================================================
    
    std::cout << "\n==========================================" << std::endl;
    std::cout << "SUMMARY OF ALL STRATEGIES" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    printf("Strategy 1 (Resolution ≤%.3f mm): %.2f mm\n", maxAllowedResolution, thickness[opt_strategy1]);
    printf("Strategy 2 (Weighted FoM, α=%.1f): %.2f mm\n", alpha, thickness[opt_strategy2]);
    printf("Strategy 3 (Normalized Product): %.2f mm\n", thickness[opt_strategy3]);
    printf("Strategy 4 (Diminishing Returns): %.2f mm\n", thickness[opt_strategy4]);
    if (opt_strategy5 >= 0) {
        printf("Strategy 5 (Eff ≥%.1f%%): %.2f mm\n", targetEfficiency*100, thickness[opt_strategy5]);
    }
    
    // Pick recommended strategy (let's use Strategy 1 or 2)
    int recommended = opt_strategy1;
    
    std::cout << "\n==========================================" << std::endl;
    std::cout << "RECOMMENDED CONFIGURATION" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "\n✓ RECOMMENDED: " << thickness[recommended] << " mm HDPE" << std::endl;
    std::cout << "  (Based on resolution constraint strategy)" << std::endl;
    std::cout << "\nPerformance:" << std::endl;
    std::cout << "  • Efficiency: " << (efficiency[recommended]*100) << "%" << std::endl;
    std::cout << "  • Resolution: " << resolution[recommended] << " mm (σ)" << std::endl;
    std::cout << "  • FWHM: " << (resolution[recommended]*2.355) << " mm" << std::endl;
    
    std::cout << "\nRationale:" << std::endl;
    std::cout << "  • Maintains excellent spatial resolution (<0.54 mm)" << std::endl;
    std::cout << "  • Provides " << (efficiency[recommended]/efficiency[0]) 
              << "× improvement over thinnest converter" << std::endl;
    std::cout << "  • Best efficiency without compromising resolution" << std::endl;
    
    // =========================================================
    // CREATE 2x2 CANVAS PLOTS
    // =========================================================
    
    TCanvas *c1 = new TCanvas("c1", "Converter Optimization", 1200, 900);
    c1->Divide(2, 2);
    
    // ---------------------------------------------------------
    // Plot 1: Efficiency vs Thickness (top-left) - NO MARKER
    // ---------------------------------------------------------
    c1->cd(1);
    gPad->SetGridx();
    gPad->SetGridy();
    TGraphErrors *grEff = new TGraphErrors(n, &thickness[0], &efficiency[0], nullptr, &effErr[0]);
    grEff->SetTitle("Detection Efficiency;Thickness (mm);Efficiency");
    grEff->SetMarkerStyle(20);
    grEff->SetMarkerColor(kBlue);
    grEff->SetMarkerSize(1.5);
    grEff->SetLineColor(kBlue);
    grEff->SetLineWidth(2);
    grEff->Draw("APL");
    
    // NO MARKER ADDED HERE (Removed marker1)
    
    // ---------------------------------------------------------
    // Plot 2: Resolution vs Thickness (top-right) - NO MARKER
    // ---------------------------------------------------------
    c1->cd(2);
    gPad->SetGridx();
    gPad->SetGridy();
    TGraphErrors *grRes = new TGraphErrors(n, &thickness[0], &resolution[0], nullptr, &resErr[0]);
    grRes->SetTitle("Spatial Resolution;Thickness (mm);Resolution (mm)");
    grRes->SetMarkerStyle(21);
    grRes->SetMarkerColor(kRed);
    grRes->SetMarkerSize(1.5);
    grRes->SetLineColor(kRed);
    grRes->SetLineWidth(2);
    grRes->GetYaxis()->SetRangeUser(0.51, 0.58);
    grRes->Draw("APL");
    
   /* // Add resolution constraint line
    TLine *constraintLine = new TLine(0.0, maxAllowedResolution, 
                                      thickness.back() + 0.1, maxAllowedResolution);
    constraintLine->SetLineColor(kGreen+2);
    constraintLine->SetLineStyle(2);
    constraintLine->SetLineWidth(2);
    constraintLine->Draw("SAME");
    
    TLatex *texConst = new TLatex();
    texConst->SetTextSize(0.03);
    texConst->SetTextColor(kGreen+2);
    texConst->DrawLatex(thickness.back() * 0.7, maxAllowedResolution + 0.002, 
                        Form("Constraint: ≤%.3f mm", maxAllowedResolution));*/
    
    // NO MARKER ADDED HERE (Removed marker2)
    
    // ---------------------------------------------------------
    // Plot 3: Efficiency-Resolution Trade-off (bottom-left)
    // ---------------------------------------------------------
    c1->cd(3);
    gPad->SetGridx();
    gPad->SetGridy();
    
    std::vector<double> resInv(n);
    for (int i = 0; i < n; i++) {
        resInv[i] = 1.0 / resolution[i];
    }
    TGraph *grCompRes = new TGraph(n, &resInv[0], &efficiency[0]);
    
    grCompRes->SetTitle("Efficiency-Resolution Trade-off;1/Resolution (mm^{-1});Efficiency");
    grCompRes->SetMarkerStyle(20);
    grCompRes->SetMarkerSize(1.5);
    grCompRes->SetLineColor(kBlue);
    grCompRes->SetLineWidth(2);
    grCompRes->Draw("APL");
    
    // Mark recommended point
    TGraph *recPoint = new TGraph(1);
    recPoint->SetPoint(0, 1.0/resolution[recommended], efficiency[recommended]);
    recPoint->SetMarkerStyle(29);
    recPoint->SetMarkerSize(3.5);
    recPoint->SetMarkerColor(kMagenta);
    recPoint->Draw("P SAME");
    
    TLatex *texRec = new TLatex();
    texRec->SetTextSize(0.04);
    texRec->SetTextColor(kMagenta);
    texRec->DrawLatex(1.0/resolution[recommended] - 0.02, efficiency[recommended] * 1.1, 
                      Form("Optimal: %.2f mm", thickness[recommended]));
    
    // Add labels for each thickness on trade-off plot
    TLatex *texLabels = new TLatex();
    texLabels->SetTextSize(0.03);
    texLabels->SetTextColor(kBlue);
    for (int i = 0; i < n; i++) {
        if (i != recommended) { // Don't label the optimal point twice
            texLabels->DrawLatex(1.0/resolution[i] + 0.01, efficiency[i], 
                                Form("%.2f mm", thickness[i]));
        }
    }
    
    // ---------------------------------------------------------
    // Plot 4: Figure of Merit Comparison (bottom-right)
    // ---------------------------------------------------------
    c1->cd(4);
    gPad->SetGridx();
    gPad->SetGridy();
    
    // Create a plot for Weighted FoM
    TGraph *grFoMWeighted = new TGraph(n, &thickness[0], &FoM_weighted[0]);
    grFoMWeighted->SetTitle("Figure of Merit (Weighted);Thickness (mm);FoM = Eff/Res^{#alpha}");
    grFoMWeighted->SetMarkerStyle(22);
    grFoMWeighted->SetMarkerColor(kOrange+2);
    grFoMWeighted->SetMarkerSize(1.5);
    grFoMWeighted->SetLineColor(kOrange+2);
    grFoMWeighted->SetLineWidth(2);
    grFoMWeighted->Draw("APL");
    
    /*// Mark optimal point from Strategy 2
    TGraph *markerFoM = new TGraph(1);
    markerFoM->SetPoint(0, thickness[opt_strategy2], FoM_weighted[opt_strategy2]);
    markerFoM->SetMarkerStyle(29);
    markerFoM->SetMarkerSize(3.0);
    markerFoM->SetMarkerColor(kRed);
    markerFoM->Draw("P SAME");
    
    TLatex *texFoM = new TLatex();
    texFoM->SetTextSize(0.035);
    texFoM->SetTextColor(kRed);
    texFoM->DrawLatex(thickness[opt_strategy2] + 0.02, FoM_weighted[opt_strategy2] * 0.95, 
                      Form("Max FoM: %.2f mm", thickness[opt_strategy2]));*/
    
    // Add annotation for alpha value
    TLatex *texAlpha = new TLatex();
    texAlpha->SetTextSize(0.04);
    texAlpha->SetTextColor(kOrange+2);
    texAlpha->DrawLatex(0.15, FoM_weighted[opt_strategy2] * 0.7, 
                        Form("#alpha = %.1f", alpha));
    
    c1->Update();
    c1->SaveAs("converter_optimization_CONSTRAINED.pdf");
    c1->SaveAs("converter_optimization_CONSTRAINED.png");
    
    std::cout << "\nPlots saved as converter_optimization_CONSTRAINED.pdf/.png" << std::endl;
    
    // Save results
    std::ofstream outfile("converter_results_CONSTRAINED.txt");
    outfile << "# Converter Optimization Results with Multiple Strategies" << std::endl;
    outfile << "# Thickness(mm) Efficiency(%) Resolution(mm) FoM_weighted FoM_normalized" << std::endl;
    for (int i = 0; i < n; i++) {
        outfile << thickness[i] << " " 
                << (efficiency[i]*100) << " "
                << resolution[i] << " "
                << FoM_weighted[i] << " " 
                << FoM_normalized[i] << std::endl;
    }
    outfile << "\n# Recommended: " << thickness[recommended] << " mm" << std::endl;
    outfile.close();
    
    std::cout << "Data saved to converter_results_CONSTRAINED.txt" << std::endl;
    std::cout << "\n==========================================" << std::endl;
}

int main() {
    converter_optimization();
    return 0;
}
