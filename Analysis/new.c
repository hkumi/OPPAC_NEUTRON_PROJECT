{
    // --- Setup input
    std::ifstream infile("../buildPosRec/position_data.txt");
    if (!infile.is_open()) {
        std::cout << "Error: Could not open position_data.txt" << std::endl;
        return;
    }

    // --- Output histograms
    TH1D* histReconX = new TH1D("histReconX", "Reconstructed X Position", 100, -50, 50);
    TH1D* histReconY = new TH1D("histReconY", "Reconstructed Y Position", 100, -50, 50);
    TH2D* histReconXY = new TH2D("histReconXY", "Reconstructed X vs Y", 100, -50, 50, 100, -50, 50);

    std::vector<double> allReconX, allReconY;
    int totalEvents = 0;
    int reconstructedEvents = 0;

    std::cout << "Processing position data..." << std::endl;

    // Temporary histograms created once and reused (avoids name collisions)
    TH1D* histXBottom = new TH1D("histXBottom_temp", "X Bottom Array", 50, -50, 50);
    TH1D* histXTop    = new TH1D("histXTop_temp",    "X Top Array",    50, -50, 50);
    TH1D* histYLeft   = new TH1D("histYLeft_temp",   "Y Left Array",   50, -50, 50);
    TH1D* histYRight  = new TH1D("histYRight_temp",  "Y Right Array",  50, -50, 50);

    // Variables to read
    int eventID, nBottom, nTop, nLeft, nRight;
    double value;

    while (infile >> eventID >> nBottom >> nTop >> nLeft >> nRight) {
        totalEvents++;

        // read vectors for this event
        std::vector<double> xBottom, xTop, yLeft, yRight;
        bool read_ok = true;

        for (int i = 0; i < nBottom; ++i) {
            if (!(infile >> value)) { read_ok = false; break; }
            xBottom.push_back(value);
        }
        for (int i = 0; i < nTop; ++i) {
            if (!(infile >> value)) { read_ok = false; break; }
            xTop.push_back(value);
        }
        for (int i = 0; i < nLeft; ++i) {
            if (!(infile >> value)) { read_ok = false; break; }
            yLeft.push_back(value);
        }
        for (int i = 0; i < nRight; ++i) {
            if (!(infile >> value)) { read_ok = false; break; }
            yRight.push_back(value);
        }

        if (!read_ok) {
            std::cerr << "Warning: incomplete event " << eventID << " — skipping." << std::endl;
            // If read failed mid-event, best to stop (file might be truncated) or skip safely:
            // break; // <-- uncomment to stop processing completely
            continue;
        }

        // require at least some points in each side (tweak thresholds as needed)
        if (xBottom.empty() || xTop.empty() || yLeft.empty() || yRight.empty()) {
            continue;
        }

        // Reset temp histograms
        histXBottom->Reset();
        histXTop->Reset();
        histYLeft->Reset();
        histYRight->Reset();

        for (double x : xBottom) histXBottom->Fill(x);
        for (double x : xTop)    histXTop->Fill(x);
        for (double y : yLeft)   histYLeft->Fill(y);
        for (double y : yRight)  histYRight->Fill(y);

        // Require enough statistics to attempt a fit
        const int minEntriesForFit = 3;
        bool goodFit = true;

        TF1* fitXBottom = nullptr;
        TF1* fitXTop    = nullptr;
        TF1* fitYLeft   = nullptr;
        TF1* fitYRight  = nullptr;

        double meanXBottom=0, sigmaXBottom=0, ampXBottom=0;
        double meanXTop=0,    sigmaXTop=0,    ampXTop=0;
        double meanYLeft=0,   sigmaYLeft=0,   ampYLeft=0;
        double meanYRight=0,  sigmaYRight=0,  ampYRight=0;

        if (histXBottom->GetEntries() >= minEntriesForFit &&
            histXTop->GetEntries()    >= minEntriesForFit &&
            histYLeft->GetEntries()   >= minEntriesForFit &&
            histYRight->GetEntries()  >= minEntriesForFit) {

            // create unique fit names to avoid collisions
            fitXBottom = new TF1(Form("gaussXBottom_evt%d", eventID), "gaus", -50, 50);
            fitXTop    = new TF1(Form("gaussXTop_evt%d",    eventID), "gaus", -50, 50);
            fitYLeft   = new TF1(Form("gaussYLeft_evt%d",   eventID), "gaus", -50, 50);
            fitYRight  = new TF1(Form("gaussYRight_evt%d",  eventID), "gaus", -50, 50);

            // Do fits quietly and request TFitResultPtr to check status
            TFitResultPtr r1 = histXBottom->Fit(fitXBottom, "QSR"); // Q: quiet, S: return result, R: use range
            TFitResultPtr r2 = histXTop->Fit(fitXTop, "QSR");
            TFitResultPtr r3 = histYLeft->Fit(fitYLeft, "QSR");
            TFitResultPtr r4 = histYRight->Fit(fitYRight, "QSR");

            if ( (int)r1 != 0 || (int)r2 != 0 || (int)r3 != 0 || (int)r4 != 0 ) {
                // non-zero status often indicates fit problems
                goodFit = false;
            } else {
                // read parameters & extra checks
                ampXBottom   = fitXBottom->GetParameter(0);
                meanXBottom  = fitXBottom->GetParameter(1);
                sigmaXBottom = fitXBottom->GetParameter(2);

                ampXTop   = fitXTop->GetParameter(0);
                meanXTop  = fitXTop->GetParameter(1);
                sigmaXTop = fitXTop->GetParameter(2);

                ampYLeft   = fitYLeft->GetParameter(0);
                meanYLeft  = fitYLeft->GetParameter(1);
                sigmaYLeft = fitYLeft->GetParameter(2);

                ampYRight   = fitYRight->GetParameter(0);
                meanYRight  = fitYRight->GetParameter(1);
                sigmaYRight = fitYRight->GetParameter(2);

                // basic sanity checks on sigma/amp
                if (!(sigmaXBottom > 0 && sigmaXTop > 0 && sigmaYLeft > 0 && sigmaYRight > 0 &&
                      ampXBottom > 0 && ampXTop > 0 && ampYLeft > 0 && ampYRight > 0)) {
                    goodFit = false;
                }
            }
        } else {
            goodFit = false;
        }

        if (!goodFit) {
            // Optionally print debug info for a few events
            // std::cerr << "Event " << eventID << " fit failed or low statistics. Skipping." << std::endl;
            if (fitXBottom) delete fitXBottom;
            if (fitXTop)    delete fitXTop;
            if (fitYLeft)   delete fitYLeft;
            if (fitYRight)  delete fitYRight;
            continue;
        }

        // Weighted mean (assumes you have these functions)
        double reconX = calculateWeightedMeanX(meanXBottom, ampXBottom, sigmaXBottom,
                                               meanXTop,   ampXTop,   sigmaXTop);
        double reconY = calculateWeightedMeanY(meanYLeft, ampYLeft, sigmaYLeft,
                                               meanYRight, ampYRight, sigmaYRight);

        if (std::isfinite(reconX) && std::isfinite(reconY) &&
            reconX > -50 && reconX < 50 && reconY > -50 && reconY < 50) {
            histReconX->Fill(reconX);
            histReconY->Fill(reconY);
            histReconXY->Fill(reconX, reconY);
            allReconX.push_back(reconX);
            allReconY.push_back(reconY);
            reconstructedEvents++;
        }

        // cleanup fit objects (temp histograms are reused)
        if (fitXBottom) delete fitXBottom;
        if (fitXTop)    delete fitXTop;
        if (fitYLeft)   delete fitYLeft;
        if (fitYRight)  delete fitYRight;

        if (totalEvents % 100 == 0) {
            std::cout << "Processed " << totalEvents << " events, reconstructed " 
                      << reconstructedEvents << " events" << std::endl;
        }
    } // end while

    infile.close();

    std::cout << "\n==========================================" << std::endl;
    std::cout << "PROCESSING COMPLETE" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "Total events in file: " << totalEvents << std::endl;
    std::cout << "Successfully reconstructed: " << reconstructedEvents << std::endl;
    std::cout << "Reconstruction efficiency: " 
              << (totalEvents ? 100.0 * reconstructedEvents / totalEvents : 0.0) << "%" << std::endl;

    if (!allReconX.empty()) {
        calculateFinalResolution(allReconX, allReconY); // assumes you implemented it
    }

    TFile* outFile = new TFile("reconstruction_from_txt.root", "RECREATE");
    histReconX->Write();
    histReconY->Write();
    histReconXY->Write();
    outFile->Close();

    // Visualization (unchanged; color calls are fine in ROOT macros)
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
    TPaveText* stats = new TPaveText(0.1, 0.1, 0.9, 0.9);
    stats->AddText(Form("Total Events: %d", totalEvents));
    stats->AddText(Form("Reconstructed: %d", reconstructedEvents));
    stats->AddText(Form("Efficiency: %.1f%%", (totalEvents?100.0 * reconstructedEvents / totalEvents:0.0)));
    stats->Draw();

    c1->SaveAs("reconstruction_summary.png");

    // delete temp histograms if you want to free memory
    delete histXBottom;
    delete histXTop;
    delete histYLeft;
    delete histYRight;

    std::cout << "Saved results to reconstruction_from_txt.root and reconstruction_summary.png" << std::endl;
}
