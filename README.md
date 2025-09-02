# COINT
A collection of procedures for unit root test with structural breaks and fully-modified estimator of cointegrating equations. For example, KPSS/ADF test with 1- and 2- break; CCR, FM-OLS, FM-GMM, and FM-VAR etc.. 
In R environment, you can implement 

          devtools::install_github("tsungwu/COINT") 

to install COINT.

COINT includes procedures include Phillips' (1995) Fully-Modofied VAR and OLS <https://doi.org/10.2307/2171721>,  Kitamura and Phillips' (1997) Fully-Modofied GMM <https://doi.org/10.1016/S0304-4076(97)00004-3>, and Park's (1992) CCR(Canonical Cointegrating Regression) <https://doi.org/10.2307/2951679>, and so on. Tests with 1 or 2 structural breaks include Gregory and Hansen (1996) <https://doi.org/10.1016/0304-4076(69)41685-7>, ADF with 1 break of Zivot and Andrews (1992) <https://doi.org/10.2307/1391541>, and KPSS with 1 break by Kurozumi (2002) <https://doi.org/10.1016/S0304-4076(01)00106-3>.
