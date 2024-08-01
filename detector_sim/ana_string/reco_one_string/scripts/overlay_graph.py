import ROOT

def plot_combined_graphs(root_file, graph1_name, graph2_name, output_image):
    # Open the ROOT file
    file = ROOT.TFile(root_file, "READ")
    
    # Get the graphs from the file
    graph1 = file.Get(graph1_name)
    graph2 = file.Get(graph2_name)
    
    if not graph1 or not graph2:
        print("Error: One or both graphs not found in the file.")
        return
    
    # Create a canvas
    canvas = ROOT.TCanvas("canvas", "Combined Graphs", 800, 600)
    
    # Draw the first graph
    graph2.SetMarkerStyle(1)
    graph2.SetMarkerColor(ROOT.kBlue)
    graph2.SetLineColor(ROOT.kBlue)
    graph2.Draw("AP")
    
    # Draw the second graph on the same canvas
    graph1.SetMarkerStyle(1)
    graph1.SetMarkerColor(ROOT.kRed)
    graph1.SetLineColor(ROOT.kRed)
    graph1.Draw("P SAME")
    
    graph1.SetMarkerColorAlpha(4, 0.5)
    graph2.SetMarkerColorAlpha(2, 0.5)

    # Add a legend
    legend = ROOT.TLegend(0.1, 0.7, 0.3, 0.9)
    legend.AddEntry(graph1, "Graph 1", "P")
    legend.AddEntry(graph2, "Graph 2", "P")
    legend.Draw()
    
    # Save the canvas as an image
    canvas.SaveAs(output_image)
    
    # Clean up
    file.Close()

# Usage example
plot_combined_graphs("../build/outputfilename.root", "graph_sig", "graph_bkg", "combined_graph.png")

