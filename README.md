This is a model of auditory cortex that was used in the manuscript "Mathematical framework for place coding in the auditory system" submitted to PLoS Computational biology and uploaded into BioRxiv. It is a 2 dimensional network with pyramidal cells (henceforth termed E cells) and fast spiking interneurons (I cells).  The firing and synaptic properties and the network architecture are based on experimental data from Levy and Reyes, J. Neurosci 2012 and the model is a simplified version of the model in Levy and Reyes, PLoS Comp. 2011 (only feedforward configuration where E and I cells receive external inputs and the I cells send inhibitory inputs to E cells). The spatial profiles of connections between E and I cells are Gaussian (standard deviations given in line 52 of audCTX.m)

There are 2 Matlab files.  audCTX.m is the actual simulation and AudCTX_analysis.m is for analyzing the data.  

Running audCTX.m
1. Create a file that contains the connection matrices between E and I cells. To do this, specify the number of cells (line 5; total network size is ncells x ncells). Start off with ncells = 100 until you get used to the program as it may take a long time the network is too large.  Set RW (line 21) to 1 and specify a name for the file (line 26).  Run it once to create the file.  If you want to use the same connection schemes for subsequent simulations (it's faster if the network isn't made eac time), set RW to 0 and it will read the file specified on line 26.

2. Setting parameters: 

a. Specify a file name that will contain the data (line 25).

b. "repetition" is number of sweeps.  For each run, the repetition number will be added to the file name.

c. External input: The external input is a set of impulse trains converted to synaptic barrages and delivered to a set of E and I neurons.  The input is Gaussian distributed in space (X-Y plane; cells in the center of the Gaussian receives the largest input). 

  Two sets of excitatory inputs may be delivered to the network with standard deviations inSigmaA and inSigmaB (lies 10,11) with separation sepX (line 12).  The units are in cell number. To deliver a single input, set inSigma B to 0; doing so delivers a single input to the center of the network with standard deviaion inSigmaA. 
  The inhibitory input spread is specified by inSigmaI (line 17). There are two modes:
  
a. In the 'normal' mode (inhMode=0, line 18), the inhibition is aligned with the excitatory inputs.

b. In the 'independent' mode (inhMode=1), the location of the inhibition can be uncoupled from the excitation (as in Fig. 7 of manuscript)
    1. To be able to vary the effective inhibitory sigma to the E cells, the I cells are bypassed and instead the thalamic inputs to I cells are converted to 
    inhibitory inputs to E cells.
    2. In this mode, inSigmaB should be set to zero. sepX(line 12) then sets the separaiton between inSigmaA and inSigmaI.
    3. To reproduce Fig. 8, you will need to manually set sepX to zero in line 293 and set inSigmaB to some non-zero value.
    
Running the analysis program AudCTX_analysis(baseName) where 'baseName' is the file name specified in line 25 of audCTX.m, without the repeition number ("myFile")
1. The analysis program runs immediately after the simulation (line 30 of audCTX.m) but can be commented out and run separately.
2. The program generates 3 graphs:

a. Figure 1 shows the locations of neurons that fired (blue dots), the outermost boundary of activated cells (red), and the fitted circle to the boundary (black).

b. Figures 2 and 3 shows the excitatory and inhibitory conductances generated in the E cells.

3. The program calculates the mean (+/- SD) number of active spikes (meanArea, stdArea on line 12, 13), diameter of fitted circle of active cells (lines 14,15), diameter of the summed excitatory/inhibitory inputs (lines 16,17) that exceeded rheobase current (cThold, line 3, see Manuscript for details), the projection (lines 18,19) of the active cells on tonotopic axis (X axis), and the projection of the synaptic current that exceeded rheobase (lines 20,21).

5. These data are stored in a file of a given name (line 23).
