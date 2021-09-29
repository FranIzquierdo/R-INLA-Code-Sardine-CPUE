# R-INLA-Code-Sardine-CPUE
R Code corresponding to the manuscript: Bayesian spatio-temporal CPUE standardization: study case of European sardine (Sardina pilchardus) in the western coast of Portugal

The exploratory analysis code is not included, Protocol from Zuur et al. (2016) was applied.

Script 1_Linear_model_comparison: We can find a loop function to perform automatically all possible combination of variables through an INLA linear model. We did this approach to select the best effort variable and also to have an idea about the effect (linear effects) of the environmental variables over the response.

Script 2_Selection_of_models: We start constructing the model from the effort variable. In each of the steps we select and add a new component to the model. Model selection in terms of WAIC, DIC and LCPO. 0) Select the best effort variable (linear or Rw2 effects). 1) Select the best  ID variable (iid effects). 234) try different spatio-temporal effects (progressive structure spde + AR1 month, year iid and month rw2 cyclic). 3) Add the environmental variables. As the envars were not relevant, 4) we included the environmental variables without spatio-temporal component in order to check if there was a confounding. This step was repeated also without the rest of model terms.

Script 3_:_Best_model_&_prediction: Best model, prediction and plot effects. 


References:

Advanced Spatial Modeling with Stochastic Partial Differential Equations Using R and INLA. Elias T. Krainski, Virgilio Gómez-Rubio, Haakon Bakka, Amanda Lenzi, Daniela Castro-Camilo, Daniel Simpson, Finn Lindgren and Håvard Rue. CRC Press/Taylor and Francis Group, 2019. (https://becarioprecario.bitbucket.io/spde-gitbook/ch-stapp.html#sec:hgst)

Zuur, A. F., Ieno, E. N., & Saveliev, A. A. (2017). Spatial, Temporal and Spatial-Temporal Ecological Data Analysis with R-INLA. Highland Statistics Ltd, 1.
