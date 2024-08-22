# DNFSC
The DNFSC package implements stability-based criteria for determining the number of factors in linear factor models, as proposed by Lee et al. (2024+).

To install this package in R, run the following commands:  

```R
library(devtools) 
install_github("Arthurlee51/DNFSC")
```

## Usage 
dnfsc_sc(X, K.cand)


## Brief description
The DNFSC package implements the stability-based criteria proposed by Lee et al. (2024+) for determining the number of factors. The main function, dnfsc_sc(), computes three criteria—SC1, SC2, and SC3—along with the estimated number of factors and the loading instability measure.  For detailed parameter descriptions, usage instructions, and examples, execute the following after installation:

```R
?dnfsc_sc
```

## Reference 
Sze Ming Lee and Yunxiao Chen. **Determining number of factors under stability considerations**. 2024+. (manuscript)
