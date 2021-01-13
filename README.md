# R-INLA-Code-Sardine-CPUE
R Code corresponding to the manuscript: Bayesian spatio-temporal CPUE standardization: studycase of European sardine (Sardina pilchardus)commercial data in the western coast of Portugal

The exploratory analysis code is not included, Protocol from Zuur et al. (2016) was applied.

Script 1: we can find a loop function to perform automatically all possible combination of variables through an INLA linear model. We did this approach to select the best effort variable.

Script 2: Random effects (iid) model comparison for Vessel, Harbour, Region and the possible nested combinations. We selected Vessel ID as the best variable.

Script 3: Base model comparison including vessel characteristics iid term, spatio-temporal (AR1) structure term, year (iid) and month (RW2) temporal terms. Once the base model was reached, 5 environmental variables as smoothed (RW2) terms were tried. 


References:

Advanced Spatial Modeling with Stochastic Partial Differential Equations Using R and INLA. Elias T. Krainski, Virgilio Gómez-Rubio, Haakon Bakka, Amanda Lenzi, Daniela Castro-Camilo, Daniel Simpson, Finn Lindgren and Håvard Rue. CRC Press/Taylor and Francis Group, 2019. (https://becarioprecario.bitbucket.io/spde-gitbook/ch-stapp.html#sec:hgst)

Zuur, A. F., Ieno, E. N., & Saveliev, A. A. (2017). Spatial, Temporal and Spatial-Temporal Ecological Data Analysis with R-INLA. Highland Statistics Ltd, 1.
