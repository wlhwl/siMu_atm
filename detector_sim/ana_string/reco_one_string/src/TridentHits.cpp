#define TridentHits_cxx
#include "trident/TridentHits.hpp"
#include "trident/DataFormats.hpp"
#include <cmath>
#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "Utils.hpp"

// Logging macro to include file, line number, and function name
#define LOG_INFO(message, ...) Logger::getInstance().info("[{}:{}] " message, Logger::getRelativePath(__FILE__), __LINE__, ##__VA_ARGS__)
#define LOG_WARN(message, ...) Logger::getInstance().warn("[{}:{}] " message, Logger::getRelativePath(__FILE__), __LINE__, ##__VA_ARGS__)

const float speed_of_light = 0.3; // m/ns

// Define a helper function to calculate the sign of a value
template<typename T>
int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

void TridentHits::Loop()
{
//   In a ROOT session, you can do:
//      root> .L TridentHits.C
//      root> TridentHits t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TH1F* h_score_sig =new TH1F("score_sig", "", 50,-20,0);
   TH1F* h_score_bkg =new TH1F("score_bkg", "", 50,-20,0);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // Find prime DOM
      // Merge hits into dom hits
      std::map<int, std::vector<Hit>> dom_hits_map;
      std::map<int, std::vector<Hit>> dom_hits_map_sig;
      for (size_t j = 0; j < t0->size(); ++j) {
          Hit hit = {(*t0)[j], (*DomId)[j], (*x0)[j], (*y0)[j], (*z0)[j], (*Type)[j]};
          dom_hits_map[(*DomId)[j]].push_back(hit);
          if ((*Type)[j] == 1) dom_hits_map_sig[(*DomId)[j]].push_back(hit);
      }

      // Get coordiante of prime dom
      // Variables to track the maximum size and corresponding DomId
      int max_dom_id = -1;
      int max_dom_id_sig = -1;
      size_t max_size = 0;
      size_t max_size_sig = 0;
      float prime_x, prime_y, prime_z, prime_t = 0;

      // Iterate through the map to find the DomId with the largest vector size
      for (const auto& pair : dom_hits_map) {
          if (pair.second.size() > max_size) {
              max_size = pair.second.size();
              max_dom_id = pair.first;
              prime_x = pair.second[0].x0;
              prime_y = pair.second[0].y0;
              prime_z = pair.second[0].z0;
              prime_t = pair.second[0].t0;
          }
      }
      TLorentzVector prime_p4 = TLorentzVector();
      prime_p4.SetXYZT(prime_x, prime_y, prime_z, prime_t);

      for (const auto& pair : dom_hits_map_sig) {
          if (pair.second.size() > max_size_sig) {
              max_size_sig = pair.second.size();
              max_dom_id_sig = pair.first;
          }
      }

      // Some basic checks
      if (ientry % 1000 == 0) {
          auto result_hDOM = Utils::CountUniqueElements(*DomId);
          int nhDOMs_have_hits = result_hDOM.first;
          std::map<int, int> hDOM_frequencyMap = result_hDOM.second;
          int nhDOM_passL1 = Utils::PrintElementsWithMinFrequency(hDOM_frequencyMap, 2);

          auto result_string = Utils::CountUniqueElements( Utils::DivideElementsBy(*DomId, 20) );
          int nStrings_have_hits = result_string.first;
          LOG_INFO("Event {} has {} Hits and {} hDOMs have hits, {} hDOMs passing L1 trigger, {} strings has hits", 
                  jentry, PmtId->size(), nhDOMs_have_hits, nhDOM_passL1, nStrings_have_hits); 

          // Output the result
          if (max_dom_id != -1) {
              LOG_INFO("Event {}: DomId with the largest hit number: {}, number of hits: {}, number of signal hits {}",
                      jentry, max_dom_id, max_size, max_size_sig);
          } else {
              std::cout << "The map is empty." << std::endl;
          }
      }

      // Calculate score for each hit, loop over all hits
      // Plot score for both signal and background
      for (const auto& pair : dom_hits_map) {
          if (pair.first != max_dom_id) {
              vector<Hit> hits = pair.second;
              prime_x = hits[0].x0;
              prime_y = hits[0].y0;
              prime_z = hits[0].z0;
              prime_t = hits[0].t0;
              TLorentzVector hit_p4 = TLorentzVector();
              hit_p4.SetXYZT(hits[0].x0, hits[0].y0, hits[0].z0, hits[0].t0);
              TLorentzVector dR = hit_p4 - prime_p4;
              float mag2 = dR.Mag2();
              float score = log(1+abs(mag2))*sign(mag2);

              /*
              log_LorentzMag2 = (np.log(1 + (relative_rt.LorentzMag2).abs()) * np.sign(relative_rt.LorentzMag2)).to_numpy()
              type_hits = relative_rt.Type.to_numpy()
              relative_rt['score'] = -np.abs(log_LorentzMag2)
              */
              for (auto hit : hits) {
                if (hit.Type == 1) h_score_sig->Fill(score);
                if (hit.Type == 0) h_score_bkg->Fill(score);
              }
          }
      }


   }
   // Create a ROOT file
    //TFile *file = new TFile(outputfilename + ".root", "RECREATE");
    TFile *file = new TFile("outputfilename.root", "RECREATE");

    // Write the histogram to the ROOT file
    h_score_sig->Write();
    h_score_bkg->Write();

    // Close the file
    file->Close();
}

void TridentHits::SetOutputname(std::string outfilename){ 
    outputfilename =  outfilename;
}
