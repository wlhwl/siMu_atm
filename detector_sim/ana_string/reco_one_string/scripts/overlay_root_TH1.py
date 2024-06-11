import ROOT
from tdrstyle import *

setTDRStyle()

def PrintHist(hist):
    # Get the number of bins
    nbins = hist.GetNbinsX()

    # Print bin contents
    for bin in range(1, nbins + 1):
        bin_center = hist.GetBinCenter(bin)
        content = hist.GetBinContent(bin)
        print(f"Bin {bin_center}: {content}")

tag = "u"
#tag = "e"

cut = 12 
if tag == "u":
    cut = 11

# Load the ROOT file containing the histograms
file = ""
if tag == "e":
    file = ROOT.TFile("../build/outputfilename_ne_big.root", "READ")
if tag == "u":
    file = ROOT.TFile("../build/outputfilename_nu_big.root", "READ")

# Get the histograms from the ROOT file
h_score_sig = file.Get("score_sig")
h_score_bkg = file.Get("score_bkg")
#h_score_sig = file.Get("causality_sig")
#h_score_bkg = file.Get("causality_bkg")

# Create a canvas
c = ROOT.TCanvas("c", "c", 800, 600)

# Set colors for the histograms
h_score_sig.SetLineColor(ROOT.kBlue)
h_score_bkg.SetLineColor(ROOT.kRed)

h_score_sig.Scale(1/h_score_sig.Integral())
h_score_bkg.Scale(1/h_score_bkg.Integral())

h_score_sig.Rebin(2)
h_score_bkg.Rebin(2)
# Draw the histograms on the same canvas
h_score_sig.Draw("hist text")
h_score_sig.GetXaxis().SetTitle("Score")
h_score_sig.GetYaxis().SetTitle("Fraction of hits")
h_score_bkg.Draw("hist text same")

print ("sig")
#PrintHist(h_score_sig)
print (h_score_sig.Integral(cut, 50))

print ("bkg")
#PrintHist(h_score_bkg)
print (1-h_score_bkg.Integral(cut, 50))
# Create a legend
legend = ROOT.TLegend(0.75, 0.7, 0.95, 0.9)
if tag == "e":
    legend.AddEntry(h_score_sig, "Signal (#nu_{e})", "l")
if tag == "u":
    legend.AddEntry(h_score_sig, "Signal (#nu_{#mu})", "l")
legend.AddEntry(h_score_bkg, "Background", "l")
legend.Draw()

# Save the plot as an image file
if tag == "e":
    c.SaveAs("overlay_histograms_e.png")
    c.SaveAs("overlay_histograms_e.pdf")
if tag == "u":
    c.SaveAs("overlay_histograms_u.png")
    c.SaveAs("overlay_histograms_u.pdf")
#c.SaveAs("overlay_causality_histograms.png")

file.Close()

