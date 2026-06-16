**DES** is supported by **StreamObs**.

# Available releases
## DES Y6 Gold

We have added the DES Y6 Gold as a supported survey to use with streamobs. 
The survey dataset is described in [Bechtol et al. 2025](https://arxiv.org/abs/2501.05739), and the catalogs can are documented/publically available from [DESDM](https://des.ncsa.illinois.edu/releases). 
The maglim, completeness, and photoerror files should be downloaded and placed in the `data/surveys/des_y6/` folder and loaded in the following manner:
```
des_y6= surveys.Survey.load(survey = 'des', release='y6')
```
Any questions about the creation of these survey specific files can be addressed to Peter Ferguson. 