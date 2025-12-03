// check_positions.C
{
    ifstream infile("../buildPosRec/position_data.txt");
    if (!infile.is_open()) {
        cout << "Error opening file!" << endl;
        return;
    }
    
    TH1D* histAllX = new TH1D("histAllX", "All X Positions", 200, -100, 100);
    TH1D* histAllY = new TH1D("histAllY", "All Y Positions", 200, -100, 100);
    TH2D* histXY = new TH2D("histXY", "All X vs Y", 200, -100, 100, 200, -100, 100);
    
    int eventID, nBottom, nTop, nLeft, nRight;
    double value;
    int totalEvents = 0;
    int eventsWithData = 0;
    
    while (infile >> eventID >> nBottom >> nTop >> nLeft >> nRight) {
        totalEvents++;
        
        // Read and store all positions
        for (int i = 0; i < nBottom; i++) {
            infile >> value;
            histAllX->Fill(value);  // Bottom array measures X
        }
        for (int i = 0; i < nTop; i++) {
            infile >> value;
            histAllX->Fill(value);  // Top array measures X
        }
        for (int i = 0; i < nLeft; i++) {
            infile >> value;
            histAllY->Fill(value);  // Left array measures Y
        }
        for (int i = 0; i < nRight; i++) {
            infile >> value;
            histAllY->Fill(value);  // Right array measures Y
        }
        
        if (nBottom + nTop + nLeft + nRight > 0) {
            eventsWithData++;
        }
    }
    
    infile.close();
    
    cout << "\n==========================================" << endl;
    cout << "DATA ANALYSIS" << endl;
    cout << "==========================================" << endl;
    cout << "Total events: " << totalEvents << endl;
    cout << "Events with any data: " << eventsWithData << endl;
    cout << "Percentage with data: " << (100.0 * eventsWithData / totalEvents) << "%" << endl;
    
    cout << "\nX Position Statistics:" << endl;
    cout << "  Mean: " << histAllX->GetMean() << " mm" << endl;
    cout << "  RMS: " << histAllX->GetRMS() << " mm" << endl;
    cout << "  Entries: " << histAllX->GetEntries() << endl;
    
    cout << "\nY Position Statistics:" << endl;
    cout << "  Mean: " << histAllY->GetMean() << " mm" << endl;
    cout << "  RMS: " << histAllY->GetRMS() << " mm" << endl;
    cout << "  Entries: " << histAllY->GetEntries() << endl;
    
    // Draw
    TCanvas* c1 = new TCanvas("c1", "Position Analysis", 1200, 800);
    c1->Divide(2, 2);
    
    c1->cd(1);
    histAllX->Draw();
    
    c1->cd(2);
    histAllY->Draw();
    
    c1->cd(3);
    histXY->Draw("colz");
    
    c1->cd(4);
    TPaveText* stats = new TPaveText(0.1, 0.1, 0.9, 0.9);
    stats->AddText(Form("Total events: %d", totalEvents));
    stats->AddText(Form("Events with data: %d", eventsWithData));
    stats->AddText(Form("X mean: %.2f mm", histAllX->GetMean()));
    stats->AddText(Form("X RMS: %.2f mm", histAllX->GetRMS()));
    stats->AddText(Form("Y mean: %.2f mm", histAllY->GetMean()));
    stats->AddText(Form("Y RMS: %.2f mm", histAllY->GetRMS()));
    stats->Draw();
    
    c1->SaveAs("position_analysis.png");
}
