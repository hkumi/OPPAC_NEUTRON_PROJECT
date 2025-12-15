// debug_reconstruction.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>

using namespace std;

struct EventData {
    int eventID;
    vector<double> xBottom, xTop, yLeft, yRight;
    
    EventData(int id) : eventID(id) {}
    
    void addXBottom(double x) { xBottom.push_back(x); }
    void addXTop(double x) { xTop.push_back(x); }
    void addYLeft(double y) { yLeft.push_back(y); }
    void addYRight(double y) { yRight.push_back(y); }
    
    void getStats(const vector<double>& pos, double& mean, double& sigma, double& N) {
        N = pos.size();
        if (N == 0) { 
            mean = 0.0; 
            sigma = 1.0; 
            return; 
        }
        
        double sum = 0.0;
        for (double p : pos) sum += p;
        mean = sum / N;
        
        double sumSq = 0.0;
        for (double p : pos) sumSq += (p - mean) * (p - mean);
        sigma = sqrt(sumSq / N);
        if (sigma < 0.1) sigma = 0.1;
    }
};

vector<EventData> readData(const string& filename) {
    vector<EventData> events;
    ifstream infile(filename);
    
    if (!infile.is_open()) {
        cerr << "Cannot open " << filename << endl;
        return events;
    }
    
    string line;
    while (getline(infile, line)) {
        istringstream iss(line);
        int id, nb, nt, nl, nr;
        if (!(iss >> id >> nb >> nt >> nl >> nr)) continue;
        
        EventData event(id);
        double pos;
        for (int i = 0; i < nb; i++) { iss >> pos; event.addXBottom(pos); }
        for (int i = 0; i < nt; i++) { iss >> pos; event.addXTop(pos); }
        for (int i = 0; i < nl; i++) { iss >> pos; event.addYLeft(pos); }
        for (int i = 0; i < nr; i++) { iss >> pos; event.addYRight(pos); }
        
        events.push_back(event);
    }
    
    infile.close();
    return events;
}

int main() {
    string filename = "../buildPosRec/position_data.txt";
    vector<EventData> events = readData(filename);
    
    if (events.empty()) {
        cerr << "No events read!" << endl;
        return 1;
    }
    
    cout << "Total events: " << events.size() << endl;
    
    // Analyze first 10 events in detail
    for (int i = 0; i < min(10, (int)events.size()); i++) {
        const auto& ev = events[i];
        int total = ev.xBottom.size() + ev.xTop.size() + ev.yLeft.size() + ev.yRight.size();
        
        cout << "\n=== Event " << ev.eventID << " ===" << endl;
        cout << "Total photons: " << total << endl;
        cout << "Bottom array: " << ev.xBottom.size() << " photons" << endl;
        cout << "Top array: " << ev.xTop.size() << " photons" << endl;
        cout << "Left array: " << ev.yLeft.size() << " photons" << endl;
        cout << "Right array: " << ev.yRight.size() << " photons" << endl;
        
        // Check if any photons exist
        if (ev.xBottom.size() > 0) {
            cout << "  Bottom positions range: " 
                 << *min_element(ev.xBottom.begin(), ev.xBottom.end()) << " to "
                 << *max_element(ev.xBottom.begin(), ev.xBottom.end()) << " mm" << endl;
        }
        if (ev.xTop.size() > 0) {
            cout << "  Top positions range: " 
                 << *min_element(ev.xTop.begin(), ev.xTop.end()) << " to "
                 << *max_element(ev.xTop.begin(), ev.xTop.end()) << " mm" << endl;
        }
    }
    
    // Count statistics
    int zeroPhotons = 0;
    int lessThan5 = 0;
    int lessThan10 = 0;
    int validEvents = 0;
    
    for (const auto& ev : events) {
        int total = ev.xBottom.size() + ev.xTop.size() + ev.yLeft.size() + ev.yRight.size();
        if (total == 0) zeroPhotons++;
        if (total < 5) lessThan5++;
        if (total < 10) lessThan10++;
        if (total >= 10) validEvents++;
    }
    
    cout << "\n=== Statistics ===" << endl;
    cout << "Events with 0 photons: " << zeroPhotons << endl;
    cout << "Events with <5 photons: " << lessThan5 << endl;
    cout << "Events with <10 photons: " << lessThan10 << endl;
    cout << "Events with >=10 photons: " << validEvents << endl;
    
    return 0;
}
