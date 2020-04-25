# Modeling-and-Analyzing-Neural-Signals-with-Phase-Variability-using-the-Fisher-Rao-Registration
This project is to compare our Fisher-Rao registration (FRR) framework to the Dynamic Time Warping (DTW)  to analyze neural signals such as EEG and fMRI where phase variability plays an important role in the data. The code provided here illustrate comparisons with simulation examples (simulated data are generated in each corresponding matlab code). 


1. DynamicProgrammingQ.c -- c code for implementation of FRR method to align q2 to q1 (SRVFs) via dynamic programming. 
2. DynamicProgrammingQ.mexmacci64 and DynamicProgrammingQ.mexw64 -- call c code running on Matlab on specific system MacOS and Windows. 
3. fig1_illstration.m -- matlab to output results in Figure 1 in the paper to show the algorithm of DTW. 
4. fig2_alignment_properness.m -- matlab to output results in Figure 2 in the paper to show the comparison of the two alignment methods (FRR and DTW) by two simple sequences. 
5. fig3_feature_preservation.m -- matlab to output results in Figure 3 in the paper to show the comparison of the two alignment methods (FRR and DTW) in terms of feature preservation. 
6. fig4_fpca.m -- matlab to output results in Figure 4 in the paper to show the performance on functional PCA adopted on the amplitude and phase results obtained by FRR. 
