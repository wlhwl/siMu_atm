//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May 31 15:00:43 2024 by ROOT version 6.30/06
// from TTree Hits/Hits
// found on file: muon.root
//////////////////////////////////////////////////////////

#ifndef TridentHits_h
#define TridentHits_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "Logging.hpp"

// Header file for the classes stored in the TTree if any.
#include "vector"

using namespace std;

class TridentHits {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *t0;
   vector<int>     *PmtId;
   vector<int>     *DomId;
   vector<float>   *x0;
   vector<float>   *y0;
   vector<float>   *z0;
   vector<int>     *Type;

   // List of branches
   TBranch        *b_t0;   //!
   TBranch        *b_PmtId;   //!
   TBranch        *b_DomId;   //!
   TBranch        *b_x0;   //!
   TBranch        *b_y0;   //!
   TBranch        *b_z0;   //!
   TBranch        *b_Type;   //!

   std::string outputfilename;

   TridentHits(TTree *tree=0);
   virtual ~TridentHits();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void SetOutputname(std::string outfilename);
};

#endif

#ifdef TridentHits_cxx
TridentHits::TridentHits(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("muon.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("muon.root");
      }
      f->GetObject("Hits",tree);

   }
   Init(tree);
}

TridentHits::~TridentHits()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TridentHits::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TridentHits::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TridentHits::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   t0 = 0;
   PmtId = 0;
   DomId = 0;
   x0 = 0;
   y0 = 0;
   z0 = 0;
   Type = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("t0", &t0, &b_t0);
   fChain->SetBranchAddress("PmtId", &PmtId, &b_PmtId);
   fChain->SetBranchAddress("DomId", &DomId, &b_DomId);
   fChain->SetBranchAddress("x0", &x0, &b_x0);
   fChain->SetBranchAddress("y0", &y0, &b_y0);
   fChain->SetBranchAddress("z0", &z0, &b_z0);
   fChain->SetBranchAddress("Type", &Type, &b_Type);
   Notify();
}

Bool_t TridentHits::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TridentHits::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TridentHits::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TridentHits_cxx
