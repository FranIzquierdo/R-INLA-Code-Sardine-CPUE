# R-INLA-Code-Sardine-CPUE

R Code corresponding to the manuscript: Bayesian spatio-temporal CPUE standardization: study case of European sardine *(Sardina pilchardus)* in the western coast of Portugal

The exploratory analysis code is not included, Protocol from Zuur et al. (2016) was applied.

Script 1_Linear_model_comparison: We can find a loop function to perform automatically all possible combination of variables through an INLA linear model. This was done in order to select the best effort variable and also to have an idea about the effect (linear effects) of the environmental variables over the response.

Script 2_Selection_of_models: We start constructing the model from the effort variable. In each of the steps we select and add a new component to the model. Model selection in terms of WAIC, DIC and LCPO. See the begining of the script for more info.

Script 3_Best_model\_&\_prediction: Best model formula, prediction and plot effects.

References:

Advanced Spatial Modeling with Stochastic Partial Differential Equations Using R and INLA. Elias T. Krainski, Virgilio Gómez-Rubio, Haakon Bakka, Amanda Lenzi, Daniela Castro-Camilo, Daniel Simpson, Finn Lindgren and Håvard Rue. CRC Press/Taylor and Francis Group, 2019 (<https://becarioprecario.bitbucket.io/spde-gitbook/ch-stapp.html#sec:hgst>)

Zuur, A. F., Ieno, E. N., & Saveliev, A. A. (2017). Spatial, Temporal and Spatial-Temporal Ecological Data Analysis with R-INLA. Highland Statistics Ltd, 1
