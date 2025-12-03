#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TVirtualFitter.h>

struct SensorHit {
    int eventID;
    int sensorID;
    double posX;
    double posY;
    double energy;
    double time;
    int photonCount;
};

class PositionReconstructor {
private:
    std::vector<SensorHit> allHits;
    std::map<int, std::vector<SensorHit>> events;
    std::map<int, std::pair<double, double>> sensorMap;
    
    // Results storage
    std::vector<double> reconstructedX;
    std::vector<double> reconstructedY;
    std::vector<double> trueX;
    std::vector<double> trueY;
    std::vector<double> residualsX;
    std::vector<double> residualsY;
    std::vector<double> residualsR; // Radial residuals

public:
    PositionReconstructor() {
        InitializeSensorMap_1mm_54cells();
    }
    
    void InitializeSensorMap_1mm_54cells() {
        // 1mm sensor configuration with 54 cells per array
        int sensorsPerArray = 54;
        double cellSize = 1.84; // mm (1mm sensor + 0.84mm gap)
        double arrayLength = sensorsPerArray * cellSize; // ~100 mm
        double startPos = -(sensorsPerArray - 1) * cellSize / 2.0;
        
        std::cout << "Initializing 1mm sensors with " << sensorsPerArray 
                  << " cells per array" << std::endl;
        std::cout << "Cell size: " << cellSize << " mm, Array length: " 
                  << arrayLength << " mm" << std::endl;
        
        // Bottom array (0-53) - sensors below active area
        for (int i = 0; i < sensorsPerArray; i++) {
            double x = startPos + i * cellSize;
            double y = -50.0; // Fixed Y position for bottom array
            sensorMap[i] = std::make_pair(x, y);
        }
        
        // Left array (54-107) - sensors left of active area
        for (int i = 0; i < sensorsPerArray; i++) {
            double x = -50.0; // Fixed X position for left array
            double y = startPos + i * cellSize;
            sensorMap[i + sensorsPerArray] = std::make_pair(x, y);
        }
        
        // Top array (108-161) - sensors above active area
        for (int i = 0; i < sensorsPerArray; i++) {
            double x = startPos + i * cellSize;
            double y = 50.0; // Fixed Y position for top array
            sensorMap[i + 2*sensorsPerArray] = std::make_pair(x, y);
        }
        
        // Right array (162-215) - sensors right of active area
        for (int i = 0; i < sensorsPerArray; i++) {
            double x =  50.0; // Fixed X position for right array
            double y = startPos + i * cellSize;
            sensorMap[i + 3*sensorsPerArray] = std::make_pair(x, y);
        }
        
        std::cout << "Sensor map initialized. Total sensors: " << sensorMap.size() << std::endl;
    }
    
    void LoadData(const std::string& filename) {
        TFile* file = new TFile(filename.c_str(), "READ");
        TTree* tree = (TTree*)file->Get("SensorData");
        
        if (!tree) {
            std::cerr << "Error: Cannot find SensorData tree in file " << filename << std::endl;
            return;
        }
        
        SensorHit hit;
        tree->SetBranchAddress("EventID", &hit.eventID);
        tree->SetBranchAddress("SensorID", &hit.sensorID);
        tree->SetBranchAddress("PosX", &hit.posX);
        tree->SetBranchAddress("PosY", &hit.posY);
        tree->SetBranchAddress("Energy", &hit.energy);
        tree->SetBranchAddress("Time", &hit.time);
        tree->SetBranchAddress("PhotonCount", &hit.photonCount);
        
        int nEntries = tree->GetEntries();
        std::cout << "Loading " << nEntries << " sensor hits from " << filename << std::endl;
        
        for (int i = 0; i < nEntries; i++) {
            tree->GetEntry(i);
            allHits.push_back(hit);
            events[hit.eventID].push_back(hit);
        }
        
        file->Close();
        delete file;
        
        std::cout << "Loaded data for " << events.size() << " events" << std::endl;
        
        // Print some statistics
        int totalPhotons = 0;
        for (const auto& hit : allHits) {
            totalPhotons += hit.photonCount;
        }
        std::cout << "Total photons detected: " << totalPhotons << std::endl;
    }
    
    std::pair<double, double> ReconstructPosition_Weighted(int eventID) {
        if (events.find(eventID) == events.end()) {
            return std::make_pair(0.0, 0.0);
        }
        
        const auto& hits = events[eventID];
        int sensorsPerArray = 54;
        
        // Group hits by array
        std::map<std::string, std::vector<SensorHit>> arrayHits;
        for (const auto& hit : hits) {
            if (hit.sensorID < sensorsPerArray) 
                arrayHits["bottom"].push_back(hit);
            else if (hit.sensorID < 2*sensorsPerArray) 
                arrayHits["left"].push_back(hit);
            else if (hit.sensorID < 3*sensorsPerArray) 
                arrayHits["top"].push_back(hit);
            else 
                arrayHits["right"].push_back(hit);
        }
        
        // Paper's algorithm: weighted by N/σ
        double xPos = 0.0, yPos = 0.0;
        double totalWeightX = 0.0, totalWeightY = 0.0;
        
        // Process horizontal arrays for X coordinate
        for (const auto& array : {"bottom", "top"}) {
            if (arrayHits[array].empty()) continue;
            
            double arrayPos = 0.0, arrayWeight = 0.0, arraySigma = 0.0;
            
            // Calculate weighted mean position
            for (const auto& hit : arrayHits[array]) {
                double weight = hit.photonCount;
                double pos = sensorMap[hit.sensorID].first;
                arrayPos += pos * weight;
                arrayWeight += weight;
            }
            
            if (arrayWeight > 0) {
                arrayPos /= arrayWeight;
                
                // Calculate standard deviation of the distribution
                for (const auto& hit : arrayHits[array]) {
                    double weight = hit.photonCount;
                    double pos = sensorMap[hit.sensorID].first;
                    arraySigma += weight * std::pow(pos - arrayPos, 2);
                }
                arraySigma = std::sqrt(arraySigma / arrayWeight);
                
                // Paper's weighting: weight = N/σ
                double weight = arrayWeight / (arraySigma + 1e-6); // Avoid division by zero
                xPos += arrayPos * weight;
                totalWeightX += weight;
            }
        }
        
        // Process vertical arrays for Y coordinate
        for (const auto& array : {"left", "right"}) {
            if (arrayHits[array].empty()) continue;
            
            double arrayPos = 0.0, arrayWeight = 0.0, arraySigma = 0.0;
            
            for (const auto& hit : arrayHits[array]) {
                double weight = hit.photonCount;
                double pos = sensorMap[hit.sensorID].second;
                arrayPos += pos * weight;
                arrayWeight += weight;
            }
            
            if (arrayWeight > 0) {
                arrayPos /= arrayWeight;
                
                // Calculate standard deviation
                for (const auto& hit : arrayHits[array]) {
                    double weight = hit.photonCount;
                    double pos = sensorMap[hit.sensorID].second;
                    arraySigma += weight * std::pow(pos - arrayPos, 2);
                }
                arraySigma = std::sqrt(arraySigma / arrayWeight);
                
                double weight = arrayWeight / (arraySigma + 1e-6);
                yPos += arrayPos * weight;
                totalWeightY += weight;
            }
        }
        
        // Final weighted averages
        if (totalWeightX > 0) xPos /= totalWeightX;
        else xPos = 0.0;
        
        if (totalWeightY > 0) yPos /= totalWeightY;
        else yPos = 0.0;
        
        return std::make_pair(xPos, yPos);
    }
    
    void ReconstructAllPositions() {
        reconstructedX.clear();
        reconstructedY.clear();
        trueX.clear();
        trueY.clear();
        residualsX.clear();
        residualsY.clear();
        residualsR.clear();
        
        std::cout << "Reconstructing positions for " << events.size() << " events..." << std::endl;
        
        int reconstructedEvents = 0;
        for (const auto& [eventID, hits] : events) {
            // Skip events with too few photons for reliable reconstruction
            int totalPhotons = 0;
            for (const auto& hit : hits) {
                totalPhotons += hit.photonCount;
            }
            if (totalPhotons < 1) continue; // Minimum threshold
            
            auto pos = ReconstructPosition_Weighted(eventID);
            reconstructedX.push_back(pos.first);
            reconstructedY.push_back(pos.second);
            
            // For centered beam, true position is (0,0)
            trueX.push_back(0.0);
            trueY.push_back(0.0);
            
            residualsX.push_back(pos.first);
            residualsY.push_back(pos.second);
            residualsR.push_back(std::sqrt(pos.first*pos.first + pos.second*pos.second));
            
            reconstructedEvents++;
        }
        
        std::cout << "Successfully reconstructed " << reconstructedEvents << " events" << std::endl;
    }
    
    double CalculateFWHM(TH1D* histogram) {
        // Fit Gaussian to get proper FWHM
        histogram->Fit("gaus", "Q");
        TF1* gaussFit = histogram->GetFunction("gaus");
        
        if (gaussFit) {
            double sigma = gaussFit->GetParameter(2);
            double fwhm = 2.355 * sigma;
            
            // Calculate fit quality
            double chi2 = gaussFit->GetChisquare();
            int ndf = gaussFit->GetNDF();
            double chi2_ndf = (ndf > 0) ? chi2/ndf : 0;
            
            std::cout << "Gaussian fit: μ = " << gaussFit->GetParameter(1) 
                      << " mm, σ = " << sigma << " mm" << std::endl;
            std::cout << "Fit quality: χ²/NDF = " << chi2_ndf << std::endl;
            
            return fwhm;
        }
        
        // Fallback: estimate FWHM from histogram
        double max = histogram->GetMaximum();
        double halfMax = max / 2.0;
        
        // Find bins where histogram crosses half-maximum
        int nBins = histogram->GetNbinsX();
        double leftEdge = 0, rightEdge = 0;
        bool foundLeft = false;
        
        for (int i = 1; i <= nBins; i++) {
            if (histogram->GetBinContent(i) >= halfMax && !foundLeft) {
                leftEdge = histogram->GetBinCenter(i);
                foundLeft = true;
            }
            if (histogram->GetBinContent(i) < halfMax && foundLeft) {
                rightEdge = histogram->GetBinCenter(i-1);
                break;
            }
        }
        
        return (rightEdge - leftEdge);
    }
    
    void CalculateResolution() {
        if (residualsX.empty()) {
            std::cout << "No data to calculate resolution" << std::endl;
            return;
        }
        
        // Create histograms for fitting
        TH1D* hResidualX = new TH1D("hResidualX_temp", "X Residuals", 100, -10, 10);
        TH1D* hResidualY = new TH1D("hResidualY_temp", "Y Residuals", 100, -10, 10);
        TH1D* hResidualR = new TH1D("hResidualR_temp", "Radial Residuals", 100, 0, 15);
        
        for (size_t i = 0; i < residualsX.size(); i++) {
            hResidualX->Fill(residualsX[i]);
            hResidualY->Fill(residualsY[i]);
            hResidualR->Fill(residualsR[i]);
        }
        
        std::cout << "\n=== HIGH-RESOLUTION RESULTS (1mm sensors, 54 cells) ===" << std::endl;
        
        // Calculate FWHM using Gaussian fits
        double fwhmX = CalculateFWHM(hResidualX);
        double fwhmY = CalculateFWHM(hResidualY);
        
        std::cout << "\nFINAL RESOLUTION RESULTS:" << std::endl;
        std::cout << "X-direction FWHM: " << fwhmX << " mm" << std::endl;
        std::cout << "Y-direction FWHM: " << fwhmY << " mm" << std::endl;
        std::cout << "Average FWHM: " << (fwhmX + fwhmY)/2.0 << " mm" << std::endl;
        
        // Calculate RMS for comparison
        double rmsX = hResidualX->GetRMS();
        double rmsY = hResidualY->GetRMS();
        double rmsR = hResidualR->GetRMS();
        
        std::cout << "\nRMS values:" << std::endl;
        std::cout << "X RMS: " << rmsX << " mm" << std::endl;
        std::cout << "Y RMS: " << rmsY << " mm" << std::endl;
        std::cout << "Radial RMS: " << rmsR << " mm" << std::endl;
        
        delete hResidualX;
        delete hResidualY;
        delete hResidualR;
    }
    
    void GeneratePlots(const std::string& outputName = "results_1mm_highres.root") {
        TFile* outputFile = new TFile(outputName.c_str(), "RECREATE");
        
        // Create comprehensive histograms
        TH1D* hResidualX = new TH1D("hResidualX", "X Position Residuals; #Delta X (mm); Counts", 
                                   100, -10, 10);
        TH1D* hResidualY = new TH1D("hResidualY", "Y Position Residuals; #Delta Y (mm); Counts", 
                                   100, -10, 10);
        TH1D* hResidualR = new TH1D("hResidualR", "Radial Position Residuals; R (mm); Counts", 
                                   100, 0, 15);
        TH2D* hReconstructed = new TH2D("hReconstructed", "Reconstructed Positions; X (mm); Y (mm)", 
                                       100, -50, 50, 100, -50, 50);
        TH2D* hResidual2D = new TH2D("hResidual2D", "2D Position Residuals; #Delta X (mm); #Delta Y (mm)", 
                                    100, -10, 10, 100, -10, 10);

        // NEW: Separate histograms for reconstructed X and Y
        TH1D* hReconstructedX = new TH1D("hReconstructedX", "Reconstructed X Positions; X (mm); Counts", 
                                    100, -50, 50);
        TH1D* hReconstructedY = new TH1D("hReconstructedY", "Reconstructed Y Positions; Y (mm); Counts", 
                                    100, -50, 50);
        
        // Fill histograms
        for (size_t i = 0; i < residualsX.size(); i++) {
            hResidualX->Fill(residualsX[i]);
            hResidualY->Fill(residualsY[i]);
            hResidualR->Fill(residualsR[i]);
            hReconstructed->Fill(reconstructedX[i], reconstructedY[i]);
            hResidual2D->Fill(residualsX[i], residualsY[i]);
            hReconstructedX->Fill(reconstructedX[i]);  // NEW: Fill X positions
            hReconstructedY->Fill(reconstructedY[i]);  // NEW: Fill Y positions
        }
        
        // Create detailed plots
        TCanvas* c1 = new TCanvas("c1", "Position Reconstruction - 1mm Sensors", 1200, 800);
        c1->Divide(2, 2);
        
        c1->cd(1);
        hReconstructed->Draw("colz");
        hReconstructed->SetTitle("Reconstructed Positions (1mm sensors); X (mm); Y (mm)");
        
        c1->cd(2);
        hResidual2D->Draw("colz");
        hResidual2D->SetTitle("2D Position Residuals; #Delta X (mm); #Delta Y (mm)");
        
        c1->cd(3);
        hResidualX->Draw();
        hResidualX->SetLineColor(kBlue);
        hResidualX->SetLineWidth(2);
        hResidualX->Fit("gaus");
        hResidualX->SetTitle("X Residuals with Gaussian Fit; #Delta X (mm); Counts");
        
        c1->cd(4);
        hResidualY->Draw();
        hResidualY->SetLineColor(kRed);
        hResidualY->SetLineWidth(2);
        hResidualY->Fit("gaus");
        hResidualY->SetTitle("Y Residuals with Gaussian Fit; #Delta Y (mm); Counts");
        
        c1->SaveAs("highres_reconstruction_1mm.png");

        // NEW: Canvas for separate X and Y reconstruction
    TCanvas* c3 = new TCanvas("c3", "Separate X and Y Reconstruction", 1200, 600);
    c3->Divide(2, 2);
    
    c3->cd(1);
    hReconstructedX->Draw();
    hReconstructedX->SetLineColor(kBlue);
    hReconstructedX->SetLineWidth(2);
    hReconstructedX->Fit("gaus");
    hReconstructedX->SetTitle("Reconstructed X Positions; X (mm); Counts");
    
    c3->cd(2);
    hReconstructedY->Draw();
    hReconstructedY->SetLineColor(kRed);
    hReconstructedY->SetLineWidth(2);
    hReconstructedY->Fit("gaus");
    hReconstructedY->SetTitle("Reconstructed Y Positions; Y (mm); Counts");

    c3->cd(3);
        hReconstructed->Draw("colz");
        hReconstructed->SetTitle("Reconstructed Positions (1mm sensors); X (mm); Y (mm)");

    
    c3->SaveAs("separate_XY_reconstruction_1mm.png");
    
        
        // Additional canvas for radial distribution
        TCanvas* c2 = new TCanvas("c2", "Radial Distribution", 600, 400);
        hResidualR->Draw();
        hResidualR->SetLineColor(kGreen);
        hResidualR->SetLineWidth(2);
        hResidualR->SetTitle("Radial Position Residuals; R (mm); Counts");
        c2->SaveAs("radial_distribution_1mm.png");
        
        // Write everything to file
        hResidualX->Write();
        hResidualY->Write();
        hResidualR->Write();
        hReconstructed->Write();
        hResidual2D->Write();
        c1->Write();
        c2->Write();
        
        outputFile->Close();
        delete outputFile;
        
        std::cout << "\nPlots saved to:" << std::endl;
        std::cout << " - highres_reconstruction_1mm.png" << std::endl;
        std::cout << " - radial_distribution_1mm.png" << std::endl;
        std::cout << " - " << outputName << " (ROOT file with all histograms)" << std::endl;
    }
};

int main(int argc, char* argv[]) {
    std::string inputFile = "B4.root";
    if (argc > 1) {
        inputFile = argv[1];
    }
    
    std::cout << "==========================================" << std::endl;
    std::cout << "HIGH-RESOLUTION POSITION RECONSTRUCTION" << std::endl;
    std::cout << "Configuration: 1mm sensors, 54 cells" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    PositionReconstructor reconstructor;
    
    std::cout << "Loading data from: " << inputFile << std::endl;
    reconstructor.LoadData(inputFile);
    
    std::cout << "\nReconstructing positions using weighted algorithm..." << std::endl;
    reconstructor.ReconstructAllPositions();
    
    std::cout << "\nCalculating spatial resolution..." << std::endl;
    reconstructor.CalculateResolution();
    
    std::cout << "\nGenerating analysis plots..." << std::endl;
    reconstructor.GeneratePlots();
    
    std::cout << "\n==========================================" << std::endl;
    std::cout << "ANALYSIS COMPLETE" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    return 0;
}
