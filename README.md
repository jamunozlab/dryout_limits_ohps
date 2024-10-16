# Performance Dryout Limits of Oscillating Heat Pipes

The following repository contains the Python codes developed to calculate the performance dryout limits of oscillating heat pipes, based on the analytical model developed by Drolen and Smoot. The main methods for limits calculations are contained in the module called OHP_Performance_Limits_Calculator_Methods, which dependent on the application can be assumed under a constant adiabatic temperature approach or a constant cold-plate approach. This respective module and the functions contained need to be imported into a secondary Python script, where they can be called to generate the limits plots and calculations. Examples of this can be found in the Plotting Scripts directory, which are also the codes that generate the performance limits plots in Diaz-Caraveo et. al. (Fig. 8).

Any questions about the code can be addressed at cesardc1@stanford.edu.

References:
- Diaz-Caraveo, C., Wolk, K., Miesner, S., Montemayor, M., Rodriguez, A., Kumar, V., Muñoz, J.A., Daimaru, T., Furst, B.I., & Roberts, S. N. Performance-Dryout Limits of Oscillating Heat Pipes: A Comprehensive Theoretical and Experimental Determination. Journal of Thermophysics and Heat Transfer, 1-11 (2024).
- Drolen, B., Smoot, C. Performance limits of oscillating heat pipes: Theory and validation. Journal of Thermophysics and Heat Transfer, 31(4) (2017), 920-936
