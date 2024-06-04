import ROOT
from tdrstyle import *

setTDRStyle()

# Load the ROOT file containing the histograms
file = ROOT.TFile("../build/outputfilename.root", "READ")

# Get the histograms from the ROOT file
h_score_sig = file.Get("score_sig")
h_score_bkg = file.Get("score_bkg")

# Create a canvas
c = ROOT.TCanvas("c", "c", 800, 600)

# Set colors for the histograms
h_score_sig.SetLineColor(ROOT.kBlue)
h_score_bkg.SetLineColor(ROOT.kRed)

h_score_sig.Scale(1/h_score_sig.Integral())
h_score_bkg.Scale(1/h_score_bkg.Integral())

# Draw the histograms on the same canvas
h_score_sig.Draw("hist")
h_score_bkg.Draw("same hist")

# Create a legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(h_score_sig, "Signal", "l")
legend.AddEntry(h_score_bkg, "Background", "l")
legend.Draw()

# Save the plot as an image file
c.SaveAs("overlay_histograms.png")

file.Close()

