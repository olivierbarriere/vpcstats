ComputePI
========

[![Travis-CI Build Status](https://travis-ci.org/smouksassi/ComputePI.svg?branch=master)](https://travis-ci.org/smouksassi/ComputePI)

### Installation and Running information
```
# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("smouksassi/ComputePI")

```

### Usage
```
require(computePI)

VPCDATA<- computePI(
  obsdata = exampleobs, simdata = examplesim, stratify = ~ISM,
  NBINS = NULL, LLOQ = LLOQ)
  
```