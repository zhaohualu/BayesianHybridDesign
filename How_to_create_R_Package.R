#######################################

#USE THE FOLLOWING FOR A NEW FOLDER

rm(list = ls())
remove.packages("BayesianHybridDesign")

devtools::create("C:\\External\\Bayesian Hybrid Design\\BayesianHybridDesign")
devtools::document(pkg="C:\\External\\Bayesian Hybrid Design\\BayesianHybridDesign")

devtools::install_github("phe9480/BayesianHybridDesign")


#------------------



library(BayesianHybridDesign)

help(package="BayesianHybridDesign")


install.packages("remotes")

remotes::install_github("keaven/nphsim")

