from ROOT import TFile, TH1F
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

# Switch to the 'Agg' backend for matplotlib
plt.switch_backend('Agg')

# Enable LaTeX in matplotlib
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

def calculate_percentile(hist, percentile):
    cumulative_sum = np.cumsum([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)])
    total_entries = cumulative_sum[-1]
    threshold = total_entries * (percentile / 100.0)
    for i, value in enumerate(cumulative_sum):
        if value >= threshold:
            return hist.GetBinLowEdge(i + 1)
            #return i + 1

def calculate_fraction_below(hist, x_value):
    cumulative_sum = np.cumsum([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)])
    bin_index = hist.FindBin(x_value)
    #bin_index = x_value
    print (bin_index, cumulative_sum, len(cumulative_sum))
    if bin_index < 1:
        return 0.0
    return cumulative_sum[bin_index - 1] / cumulative_sum[-1]

def plot_histograms(root_file_path, histogram1_name, histogram2_name, output_base_name, legend1_name, legend2_name, additional_text, percentage):
    # Open the ROOT file
    root_file = TFile.Open(root_file_path)

    # Extract the histograms
    hist1 = root_file.Get(histogram1_name).Clone()
    hist2 = root_file.Get(histogram2_name).Clone()
    hist1.Rebin(2)
    hist2.Rebin(2)
    hist1.Sumw2()
    hist2.Sumw2()

    # normalize the histograms
    hist1.Scale(1.0 / hist1.Integral())
    hist2.Scale(1.0 / hist2.Integral())

    # extract the data from hist1
    hist1_values = [hist1.GetBinContent(i) for i in range(1, hist1.GetNbinsX() + 1)]
    hist1_errors = [hist1.GetBinError(i) for i in range(1, hist1.GetNbinsX() + 1)]
    hist1_edges = [hist1.GetBinLowEdge(i) for i in range(1, hist1.GetNbinsX() + 2)]

    # extract the data from hist2
    hist2_values = [hist2.GetBinContent(i) for i in range(1, hist2.GetNbinsX() + 1)]
    hist2_errors = [hist2.GetBinError(i) for i in range(1, hist2.GetNbinsX() + 1)]
    hist2_edges = [hist2.GetBinLowEdge(i) for i in range(1, hist2.GetNbinsX() + 2)]

    # calculate the 99th percentile x-value for hist1
    x_percentile = calculate_percentile(hist1, percentage)
    
    # calculate the fraction of entries for hist2 below this x-value
    fraction_below = calculate_fraction_below(hist2, x_percentile)

    # create the plot
    fig, ax = plt.subplots()

    # plot the first histogram
    ax.step(hist1_edges[:-1], hist1_values, where='mid', label=legend1_name, color='cornflowerblue')
    ax.fill_between(hist1_edges[:-1], hist1_values, step='mid', color='cornflowerblue', alpha=0.3)
    #ax.errorbar((np.array(hist1_edges[:-1]) + np.array(hist1_edges[1:])) / 2, hist1_values, yerr=hist1_errors, fmt='.', color='cornflowerblue')

    # plot the second histogram
    ax.step(hist2_edges[:-1], hist2_values, where='mid', label=legend2_name, color='orange')
    ax.fill_between(hist2_edges[:-1], hist2_values, step='mid', color='orange', alpha=0.3)
    #ax.errorbar((np.array(hist2_edges[:-1]) + np.array(hist2_edges[1:])) / 2, hist2_values, yerr=hist2_errors, fmt='.', color='orange')

    # customize the plot
    ax.set_xlabel('score')  # update with the appropriate label
    ax.set_ylabel('fraction of hits')
    #ax.set_title('overlay of two normalized histograms')
    ax.legend()
    ax.set_xlim(-20, 0)  # set x-axis range
    ax.set_ylim(0)  # set x-axis range
    ax.set_xticks(np.linspace(-20, 0, 11))
    #ax.yaxis.set_maijor_locator(ticker.loglocator(base=10, numticks=10))
    #ax.yaxis.set_minor_locator(ticker.loglocator(base=10, subs=np.arange(2, 10) * .1, numticks=10))

 # add additional text at the top left of the canvas
    additional_text_with_fraction = f"{additional_text}\nfraction below {x_percentile:.2f} in {legend2_name}: {fraction_below:.2%}"
    # add additional text at the top left of the canvas
    ax.text(0.02, 1.05, additional_text, transform=ax.transAxes, fontsize=12, verticalalignment='top', horizontalalignment='left',color='black')
    #ax.text(0.02, 1.10, additional_text_with_fraction, transform=ax.transaxes, fontsize=12, verticalalignment='top', horizontalalignment='left',color='black')
    print (additional_text_with_fraction)

    # save the plot as png and pdf
    plt.savefig(f'{output_base_name}.png')
    plt.savefig(f'{output_base_name}.pdf')

    # show the plot (commented out since we are using the 'agg' backend)
    # plt.show()

    # clean up root objects
    root_file.Close()

# example usage
plot_histograms(
    '../build/outputfilename_nu_big_correct.root',
    'score_sig',
    'score_bkg',
    'overlay_nu_mu_big',
    r'$\nu_\mu$',
    'noise',
    'TRIDENT preliminary',
    1 
)

plot_histograms(
    '../build/outputfilename_ne_big_correct.root',
    'score_sig',
    'score_bkg',
    'overlay_nu_e_big',
    r'$\nu_e$',
    'noise',
    'TRIDENT Preliminary',
    1 
)

# example usage
#plot_histograms('../build/outputfilename.root', 'score_sig', 'score_bkg', 'overlay_test')
#plot_histograms('../build/outputfilename_nu_big_correct.root', 'score_sig', 'score_bkg', 'overlay_nu_mu_big')
#plot_histograms('../build/outputfilename_ne_big.root', 'score_sig', 'score_bkg', 'overlay_nu_e_big')


'''
def plot_histograms(root_file_path, histogram1_name, histogram2_name, output_base_name):
    # Open the ROOT file
    root_file = TFile.Open(root_file_path)
    
    # Extract the histograms
    hist1 = root_file.Get(histogram1_name).Clone()
    hist2 = root_file.Get(histogram2_name).Clone()
    
    # Normalize the histograms
    hist1.Scale(1.0 / hist1.Integral())
    hist2.Scale(1.0 / hist2.Integral())
    
    # Extract the data from hist1
    hist1_values = [hist1.GetBinContent(i) for i in range(1, hist1.GetNbinsX() + 1)]
    hist1_errors = [hist1.GetBinError(i) for i in range(1, hist1.GetNbinsX() + 1)]
    hist1_edges = [hist1.GetBinLowEdge(i) for i in range(1, hist1.GetNbinsX() + 2)]
    
    # Extract the data from hist2
    hist2_values = [hist2.GetBinContent(i) for i in range(1, hist2.GetNbinsX() + 1)]
    hist2_errors = [hist2.GetBinError(i) for i in range(1, hist2.GetNbinsX() + 1)]
    hist2_edges = [hist2.GetBinLowEdge(i) for i in range(1, hist2.GetNbinsX() + 2)]

    # Create the plot
    fig, ax = plt.subplots()

    # Plot the first histogram
    ax.step(hist1_edges[:-1], hist1_values, where='mid', label=histogram1_name, color='cornflowerblue')
    ax.fill_between(hist1_edges[:-1], hist1_values, step='mid', color='cornflowerblue', alpha=0.3)
    ax.errorbar((np.array(hist1_edges[:-1]) + np.array(hist1_edges[1:])) / 2, hist1_values, yerr=hist1_errors, fmt='.', color='cornflowerblue')

    # Plot the second histogram
    ax.step(hist2_edges[:-1], hist2_values, where='mid', label=histogram2_name, color='orange')
    ax.fill_between(hist2_edges[:-1], hist2_values, step='mid', color='orange', alpha=0.3)
    ax.errorbar((np.array(hist2_edges[:-1]) + np.array(hist2_edges[1:])) / 2, hist2_values, yerr=hist2_errors, fmt='.', color='orange')

    # Customize the plot
    ax.set_xlabel('Score')  # Update with the appropriate label
    ax.set_ylabel('Fraction of hits')
    #ax.set_title('Overlay of Two Normalized Histograms')
    ax.legend()
    ax.set_xlim(-20, 0)  # Set x-axis range
    ax.set_xticks(np.linspace(-20, 0, 5))
    #ax.yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=10))
    #ax.yaxis.set_minor_locator(ticker.LogLocator(base=10, subs=np.arange(2, 10) * .1, numticks=10))

    # Save the plot as PNG and PDF
    plt.savefig(f'{output_base_name}.png')
    #plt.savefig(f'{output_base_name}.pdf')

    # Show the plot
    #plt.show()
    # Close the ROOT file
    root_file.Close()

# Example usage
#plot_histograms('../build/outputfilename.root', 'score_sig', 'score_bkg', 'overlay_test')
plot_histograms('../build/outputfilename_nu_big_correct.root', 'score_sig', 'score_bkg', 'overlay_nu_mu_big')
#plot_histograms('../build/outputfilename_ne_big.root', 'score_sig', 'score_bkg', 'overlay_nu_e_big')
'''

