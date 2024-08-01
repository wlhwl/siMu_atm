#define TridentHits_cxx
#include "trident/TridentHits.hpp"
#include "trident/DataFormats.hpp"
#include <cmath>
#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "Utils.hpp"
#include <TMatrixD.h>
#include <TDecompSVD.h>

// Logging macro to include file, line number, and function name
#define LOG_INFO(message, ...) Logger::getInstance().info("[{}:{}] " message, Logger::getRelativePath(__FILE__), __LINE__, ##__VA_ARGS__)
#define LOG_WARN(message, ...) Logger::getInstance().warn("[{}:{}] " message, Logger::getRelativePath(__FILE__), __LINE__, ##__VA_ARGS__)
#define LOG_DEBUG(message, ...) Logger::getInstance().debug("[{}:{}] " message, Logger::getRelativePath(__FILE__), __LINE__, ##__VA_ARGS__)

const float speed_of_light = 0.3; // m/ns

// Define a helper function to calculate the sign of a value
template<typename T>
int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

// Function to perform linear fit
#include <vector>
#include <iostream>
#include <TMatrixD.h>
#include <TDecompSVD.h>
#include <TLorentzVector.h>

void PrintTMatrix(std::string name, TMatrixD matrix) {
// Print the matrix
    LOG_INFO("Name of matrix is {}: ", name);
    for (int i = 0; i < matrix.GetNrows(); ++i) {
        for (int j = 0; j < matrix.GetNcols(); ++j) {
            std::cout << matrix(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

#include <TMatrixD.h>
#include <TDecompSVD.h>
#include <TVectorD.h>
#include <iostream>

std::vector<double> LinearFit (const std::vector<TLorentzVector>& points) {
    // Check if the input vector is not empty
    if (points.empty()) {
        throw std::invalid_argument("The vector of TLorentzVector points is empty.");
    }

    // Calculate the means of x, y, z
    double x_mean = 0;
    double y_mean = 0;
    double z_mean = 0;

    for (const auto& point : points) {
        x_mean += point.X();
        y_mean += point.Y();
        z_mean += point.Z();
    }

    x_mean /= points.size();
    y_mean /= points.size();
    z_mean /= points.size();

    // Calculate the covariance matrix
    TMatrixD covarianceMatrix(3, 3);
    for (const auto& point : points) {
        double x = point.X() - x_mean;
        double y = point.Y() - y_mean;
        double z = point.Z() - z_mean;

        covarianceMatrix(0, 0) += x * x;
        covarianceMatrix(0, 1) += x * y;
        covarianceMatrix(0, 2) += x * z;
        covarianceMatrix(1, 0) += y * x;
        covarianceMatrix(1, 1) += y * y;
        covarianceMatrix(1, 2) += y * z;
        covarianceMatrix(2, 0) += z * x;
        covarianceMatrix(2, 1) += z * y;
        covarianceMatrix(2, 2) += z * z;
    }

    // Regularization: Add a small value to the diagonal elements to improve conditioning
    double regularizationFactor = 1e-5;
    for (int i = 0; i < 3; ++i) {
        covarianceMatrix(i, i) += regularizationFactor;
    }

    // Perform eigen decomposition
    TVectorD eigenValues(3);
    TMatrixD eigenVectors = covarianceMatrix.EigenVectors(eigenValues);

    // The eigenvector corresponding to the largest eigenvalue is the direction vector
    TVector3 direction(eigenVectors(0, 2), eigenVectors(1, 2), eigenVectors(2, 2));
    TVector3 pointOnLine(x_mean, y_mean, z_mean);

    // Return the point on the line and the direction vector
    return {eigenVectors(0, 2), eigenVectors(1, 2), eigenVectors(2, 2)};
}

// Function to rotate points
TVector3 rotatePoint(const TVector3& point, const TVector3& direction) {
    TVector3 zAxis(0, 0, 1);
    TVector3 dir = direction.Unit(); // Normalize the direction vector
    TVector3 v = zAxis.Cross(dir);
    double s = v.Mag();
    double c = zAxis.Dot(dir);

    if (s == 0) {
        return point;  // Already aligned
    }

    // Construct skew-symmetric matrix kmat
    TMatrixD kmat(3, 3);
    kmat(0, 0) = 0;     kmat(0, 1) = -v.Z(); kmat(0, 2) = v.Y();
    kmat(1, 0) = v.Z(); kmat(1, 1) = 0;     kmat(1, 2) = -v.X();
    kmat(2, 0) = -v.Y(); kmat(2, 1) = v.X(); kmat(2, 2) = 0;

    // Calculate the rotation matrix
    TMatrixD I(3, 3);
    I.UnitMatrix();
    TMatrixD rotationMatrix = I + kmat + kmat * kmat * ((1 - c) / (s * s));
    TVector3 rotatedPoint;
    for (int i = 0; i < 3; ++i) {
        rotatedPoint[i] = rotationMatrix(i, 0) * point.X() + rotationMatrix(i, 1) * point.Y() + rotationMatrix(i, 2) * point.Z();
    }
    return rotatedPoint;
}

vector<Hit> RotateAll(vector<Hit> hits, TVector3 direction) {

     vector<TLorentzVector> points;
     for (auto hit : hits) {
     
          float x = hit.x0;
          float y = hit.y0;
          float z = hit.z0;
          float t = hit.t0;
          TLorentzVector p4 = TLorentzVector();
          p4.SetXYZT(x, y, z, t*speed_of_light);
          points.push_back(p4);
     }

     // Rotate all points
     std::vector<TLorentzVector> rotatedPoints;
     for (const auto& point : points) {
         TVector3 rotatedPoint = rotatePoint(TVector3(point.X(), point.Y(), point.Z()), direction);
         rotatedPoints.emplace_back(rotatedPoint.X(), rotatedPoint.Y(), rotatedPoint.Z(), point.T());
     }

     // Determine the positive z-axis direction
     auto minMaxZ = std::minmax_element(rotatedPoints.begin(), rotatedPoints.end(),
                                        [](const TLorentzVector& a, const TLorentzVector& b) {
                                            return a.Z() < b.Z();
                                        });
     TLorentzVector minZPoint = *minMaxZ.first;
     TLorentzVector maxZPoint = *minMaxZ.second;

     if (minZPoint.T() > maxZPoint.T()) {
         for (auto& point : rotatedPoints) {
             point.SetZ(-point.Z());
         }
     }

     vector<Hit> rotatedHits;
     for (int i = 0; i < rotatedPoints.size(); i++) {
         float t = hits[i].t0;
         int DomId = hits[i].DomId;
         float x = rotatedPoints[i].X();
         float y = rotatedPoints[i].Y();
         float z = rotatedPoints[i].Z();
         int Type = hits[i].Type;
         Hit hit = {t,DomId,x,y,z,Type};
         rotatedHits.push_back(hit);
     
     }
     return rotatedHits;
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
   TH1F* h_causality_sig =new TH1F("causality_sig", "", 50,-20,20);
   TH1F* h_causality_bkg =new TH1F("causality_bkg", "", 50,-20,20);
   vector<float> v_mag2_3D_sig, v_mag2_3D_bkg;
   vector<float> v_ct_sig, v_ct_bkg;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=40000; jentry<nentries;jentry++) {

      if (jentry > 50000) continue; // TODO make this configurable

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //if (t0->size() > 100000) LOG_DEBUG("This entry {} has {} hits", jentry, t0->size());
      //if (t0->size() > 100000) continue;

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
      std::vector<TLorentzVector> hits_for_fit; 

      for (const auto& pair : dom_hits_map) {
          if (pair.second.size() > max_size) {
              max_size = pair.second.size();
              max_dom_id = pair.first;
              prime_x = pair.second[0].x0;
              prime_y = pair.second[0].y0;
              prime_z = pair.second[0].z0;
              prime_t = pair.second[0].t0;
          }

          if (pair.second.size() > 2) {
             TLorentzVector fithit_p4 = TLorentzVector();
             Hit firsthit = pair.second[0];
             fithit_p4.SetXYZT(firsthit.x0, firsthit.y0, firsthit.z0, firsthit.t0*speed_of_light);
             hits_for_fit.push_back(fithit_p4);
          }
      }

      // Perform linear fit
      //LOG_DEBUG("Number of hits used for fit: {}", hits_for_fit.size());
      /*
      TVector3 direction(0,0,0);
      if (hits_for_fit.size() > 2) {
          std::vector<double> fitParams = LinearFit(hits_for_fit);
          direction.SetXYZ(fitParams[0], fitParams[1], fitParams[2]);
      } 
      */

      //LOG_DEBUG("direction: {}, {}, {}", direction.X(), direction.Y(), direction.Z());

      TLorentzVector prime_p4 = TLorentzVector();
      prime_p4.SetXYZT(prime_x, prime_y, prime_z, prime_t*speed_of_light);

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
      vector<Hit> allhits;
      for (const auto& pair : dom_hits_map) {
          int domID = pair.first;
          vector<Hit> hits = pair.second;

          if (domID != max_dom_id) {
              for (auto hit : hits) {
                  allhits.push_back(hit);
                  float x = hit.x0;
                  float y = hit.y0;
                  float z = hit.z0;
                  float t = hit.t0;
                  TLorentzVector hit_p4 = TLorentzVector();

                  hit_p4.SetXYZT(x,y,z,t*speed_of_light);
                  TLorentzVector dR = hit_p4 - prime_p4;
                  float mag2 = dR.Mag2();
                  //float score = log(1+abs(mag2))*sign(mag2);
                  float score = -log(1+abs(mag2));
                  float mag2_3D = (hit_p4.Vect() - prime_p4.Vect()).Mag2();

                  if (hit.Type == 1) {
                      h_score_sig->Fill(score);
                      if (v_mag2_3D_sig.size() < 100000) { 
                          v_mag2_3D_sig.push_back(log(mag2_3D));
                          v_ct_sig.push_back( log(pow(t*speed_of_light, 2)) );
                      }
                  }
                  if (hit.Type == 0) {
                      h_score_bkg->Fill(score);
                      if (v_mag2_3D_bkg.size() < 100000) { 
                          v_mag2_3D_bkg.push_back(log(mag2_3D));
                          v_ct_bkg.push_back( log(pow(t*speed_of_light, 2)) );
                      }
                  }
              }// loop over all hits from one dom

          }


          } // loop over hit map

          //LOG_DEBUG("direction: {}, {}, {}", direction.X(), direction.Y(), direction.Z());
          /*
          std::vector<Hit> rotatedHits = RotateAll(allhits, direction);
          for (auto hit : rotatedHits) {

              int domID = hit.DomId;
              if (domID == max_dom_id) continue;
              float x = hit.x0;
              float y = hit.y0;
              float z = hit.z0;
              float t = hit.t0;
              TLorentzVector p4 = TLorentzVector();
              float cherenkovFactor = pow(TMath::Tan(42 * TMath::DegToRad()), 2);
              p4.SetXYZT(sqrt(cherenkovFactor)*x,sqrt(cherenkovFactor)*y,z,t*speed_of_light);
              TLorentzVector prime_p4_mod = TLorentzVector();
              prime_p4_mod.SetXYZT(sqrt(cherenkovFactor)*prime_p4.X(), sqrt(cherenkovFactor)*prime_p4.Y(), prime_p4.Z(), prime_p4.T());
              float causality = (p4-prime_p4_mod).Mag2();

              //float causality = pow( ((t-prime_p4.T())*speed_of_light - (z-prime_p4.Z())) , 2) 
              //    - pow(TMath::Tan(42 * TMath::DegToRad()), 2) * ( pow(x-prime_p4.X(),2) + pow(y-prime_p4.Y(),2));

              int type = hit.Type;
              causality = TMath::Sign(1.0, causality) * TMath::Log(TMath::Abs(causality));
              //causality = 1-causality/20;
              if (type == 1) h_causality_sig->Fill(causality);
              if (type == 0) h_causality_bkg->Fill(causality);
          }*/

   }

   TGraph *graph_3D_t_sig = new TGraph(v_mag2_3D_sig.size(), v_mag2_3D_sig.data(), v_ct_sig.data());
   TGraph *graph_3D_t_bkg = new TGraph(v_mag2_3D_bkg.size(), v_mag2_3D_bkg.data(), v_ct_bkg.data());
   graph_3D_t_sig->SetMarkerColor(2);
   graph_3D_t_sig->SetName("graph_sig");
   graph_3D_t_bkg->SetMarkerColor(4);
   graph_3D_t_bkg->SetName("graph_bkg");

   // Create a ROOT file
    //TFile *file = new TFile(outputfilename + ".root", "RECREATE");
    TFile *file = new TFile("outputfilename.root", "RECREATE");

    // Write the histogram to the ROOT file
    h_score_sig->Write();
    h_score_bkg->Write();
    h_causality_sig->Write();
    h_causality_bkg->Write();
    graph_3D_t_sig->Write();
    graph_3D_t_bkg->Write();


    // Close the file
    file->Close();
}

void TridentHits::SetOutputname(std::string outfilename){ 
    outputfilename =  outfilename;
}
