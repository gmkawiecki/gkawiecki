# gkawiecki
SUMMARY

The document "Development_of_Efficient_Vertical_Axis_Wind_Turbines" presents the result of a comparative experimental study of two basic Vertical Axis Wind Turbine (VAWT) configurations: a classical one with fixed blade pitch and a novel configuration with passive blade pitch control using a linkage mechanism with a center of rotation eccentric with respect to turbine’s vertical axis (Benedict and Chopra, 2016; Erickson et al, 2011). To the best of my knowledge, this is is the first study to compare experimentally both configurations, using the same basic test specimen. I carried out this work in my garage during the pandemic lockdown.

Files for VAWT performance estimation:

1. forMST.f95: slightly modified FORTRAN code published by J. H. Strickland (DART), SAND75-0431 report.

2. STRIDATA_ORIG.txt: DART input data published in Strickland's report. Used in forMST for validation.

3. forMST_output.txt: forMST.f95 output.

4. DART_digitized_output.csv: DART output.

5. pyMSTv1.py: Python code based on Strickland's and Sanyer's work.
  
6. modMST.py: module with methods needed to run pyMSTv1.py

7. STRIDATA_ORIG.csv: pyMSTv1 input.

8. pyMSTv1_output.txt: pyMSTv1 output, for validation against Strickland results.

9. cP_py_Vs_Strickland.png: plot presenting pyMSTv1-produced results against Strickland's results.

10. fit_Strickland_data.png: plot for validating the interpolation method for Strickland's input data. 

Conclusions

Obtained results indicate that the variable blade pitch configuration of a VAWT offers far superior performance for low and moderate wind velocities and for low load resistance values.  Those are precisely the typical conditions for urban and residential environments and applications of growing interest, such as car or renewable energy production system batteries.

The agreement of modMST.py results with original Strickland's results is satisfactory. Differences showing in cP_py_Vs_Strickland.png may be attributed to different methods for cN, cT vs angle of attack interpolation methods or to a, yet undiscovered, bug. This is work in progress. 

References

Benedict, M. and Chopra I., 2016, “Aerodynamics of a Small-Scale Vertical-Axis Wind Turbine with Dynamic Blade Pitching,” AIAA Journal, Vol. 54, No. 3.

Erickson D. W., Wallace J. J. and Peraire J., 2011, " Performance Characterization of Cyclic Blade Pitch Variation on a Vertical Axis Wind Turbine," AIAA Paper 2011-638, Proceedings of 49th AIAA Aerospace Sciences Meeting including the New Horizons Forum and Aerospace Exposition.

Strickland J. H., 1975, "The Darrieus Turbine: A Performance Prediction Model Using Multiple Streamtubes," SAND75-0431, Sandia Laboratories energy report.
Sanyer W. E., 2011, "The Development of a Wind Turbine for Residential Use," M.Sc. thesis, Dept. of Mech. Eng., North Carolina State University.
