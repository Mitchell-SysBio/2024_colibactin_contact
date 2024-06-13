# 2024_colibactin_contact
Figure 1

example_metadata.xlsx is an example metadata file to use as input for both flowGate.m and flowPlot.m to sort and keep track of .fcs files in the analysis.

flowGate.m reads in .fcs files to gate on the cell population, single cell population, YFP positive population, CFP positive population, and mCherry positive population. It requires a metadata spreadsheet to sort files and extract representative samples to perform the gating on and saves the coordinates of each gate to use in downstream analyses (xyCoordinate.mat).

flowPlot.m uses xyCoordinate.mat to extract fluorescent signal for cells within each gate defined in flowGate.m to analyze positive reporter cells.

xyCoordinate.mat is an example output of flowGate.m to be used as an input for flowPlot.m

Figure 2

quant_backgroundSignal.m is a function to draw an roi line on a control, non-fluorescent lawn of cells to get background signal for each fluorescent channel used in an experiment. The median signal of the roi line is averaged across all technical replicates. This is used in the plot_lawnBackground.m wrapper script

plot_lawnBackground.m is a wrapper script that uses the quant_backgroundSignal.m function to calculate background signal on a bacterial lawn. The initial output of the function is saved as
bkgrnd.mat. The function output is then plotted over time and the background signal by channel and time is saved in a new output called backgroundSignal.mat

bkgrnd.mat is the saved output after running the quant_backgroundSignal.m function

backgroundSignal.mat is the output of the plot_lawnBackground.m wrapper script and is used as an input to plot_lawn_recASignal.m

segment_lawn_recASignal.m is a function used by the plot_lawn_recASignal.m wrapper script. The function allows the user to draw an roi line from a colony to a reporter lawn and determine signal intensity along the line.

plot_lawn_recASignal.m is a wrapper script that uses the segment_lawn_recASignal.m function to quantify fluorescent signal intensity from a producer colony to the reporter lawn. The output from the function is saved as data.mat if data needs to be plotted in the future, to save time re-drawing roi lines. This script produces several plots looking at decay curves for each measurement and the average across all measurements. Decay curve plots can be saved. Parameters such as max signal and dist50 are saved in an output file at the end called rep1.mat (or other replicate number).

data.mat is an example file that is generated by the segment_lawn_recASignal function. It can be used as an input in the plot_lawn_recASignal.m script to re-plot data.

rep1.mat, rep2.mat & rep3.mat are example replicate output files from plot_lawn_recASignal.m. These files are inputs to the compareReplicates.m script that allows for comparing measurements across biological replicates

compareReplicates.m uses output files from plot_lawn_recASignal.m to plot and compare data across biological replicates. This script plots the max YFP signal, the dist50, and the decay curves for each biological replicate over time.

Figures 3 and 4

segment_recASignal.m is a function used to calculate fluorescent profiles along a cross-section line between two colonies. A czi file is read into the function and each scene within the file is opened. The user marks a line between two colonies and additional lines mark the edge of each colony. The position of each colony border and the fluorescent signal along the line cross-sectioning the two colonies is saved for downstream analysis. recA_reps.m is a wrapper script for applying this function and visualizing the outputs.

recA_reps.m is a wrapper script that applies the segment_recASignal function to one or several czi microscopy files. The script iterates through all files, grouping by bacteria strain and condition, to extract distance between colonies and fluorescent profiles of colonies. The signal decay is plotted for each condition and parameters such as max signal and distance to 50% of the max signal response are reported and saved in output files for further analysis comparing experiments from separate dates using compare_reps.m. Plots of signal decay curves for each replicate can be saved from this script.

BAC.mat & EcN.mat are example output files from recA_reps.m of data produced from the segment_recASignal.m function. These include measurements of fluorescence and distance and colony edge markers for each colony imaged in the experiment. They are broken down into two cells, one containing pks+ data and one with pks- data.

rep1.mat, rep2.mat & rep4.mat are example output files from the end of recA_reps.mat. They have saved variables like the average signal decay, dist50, and max signal to plot and run statistical tests in compare_reps.m

compare_reps.m takes outputs from recA_reps.m (i.e. rep1.mat) to plot parameters across biological replicates and run statistical tests. Plots of signal parameters can be saved from this script.

plotColonyTraces.m reads in a czi file and allows you to select an index within the file to plot fluorescent profiles for. The user draws a line crossing the entirety of both colonies and plots the fluorescent profiles of all fluorophores (min-max normalized) across the line. 



![image](https://github.com/Mitchell-SysBio/2024_colibactin_contact/assets/154358451/6f7e0295-89e3-4ee3-8d02-4e37100a2e09)
