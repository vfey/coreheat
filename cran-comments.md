---
title: CRAN package coreheat
---

## Resubmission 2026-02-04
This is resubmission of the package. The version was increased to 0.3.0 after addressing the comments by Prof Brian Ripley regarding failure when accessing an internet resource.  
The fix required a change in the dependency 'convertid', which is resubmitted at the same time. Changes in this package are hence minor and include mainly a number of error handling iterations to allow the function to continue if the online resource is not responding.

### Test environments (2026-02-04 - )
* local OS X install: aarch64-apple-darwin25.2.0, R 4.5.2
* win-builder (devel, release and oldrelease)
* Red Hat Enterprise Linux release 9.7 (Plow), R 4.5.2

## Resubmission 2021-09-17
This is resubmission of the package. The version was increased to 0.2.2 after addressing the comments by CRAN member Gregor Seyer:  

* replaced `cat` with `message` to make verbose output be suppressable
* added argument `verbose` defaulting to `FALSE` for convenience
* The word "Ensembl" is not misspelled but refers to the Ensembl project. It is used, e.g., in the DESCRIPTION of the `biomaRt` package.  

In addition, the following corrections/additions were made:  

* extended correlation map filtering capabilities
* fixed and now exporting function to automatically split clusters based on noise level and hierarchy for filtering
* added argument controlling plot sub-title
* corrected bug and fixed behavior when input is a correlation matrix
* adapted some defaults
* minor corrections to documentation (misspelled words etc)

## Resubmission 2021-09-03
This is resubmission of the package. The version was increased to 0.2.1 after checking and slightly modifying the vignette to address the re-building error given in the __win-builder Debian__ check:

```
Error(s) in re-building vignettes:
  ...
--- re-building ‘coreheat-vignette.Rmd’ using rmarkdown
Quitting from lines 110-111 (coreheat-vignette.Rmd) 
Error: processing vignette 'coreheat-vignette.Rmd' failed with diagnostics:
Values argument contains no data.
--- failed re-building ‘coreheat-vignette.Rmd’
```

Please note that I could not reproduce the error on any of the test systems. I do not have a Debian system, I tested on a Docker image running Ubuntu 18.04.5 LTS.

## Test environments
* local OS X install: x86_64-apple-darwin17.0, R 4.0.2
* local Ubuntu Docker image: 18.04.5 LTS, R 4.0.3
* win-builder (devel and release)
* CentOS Linux release 7.9.2009 (Core) [:core-4.1-amd64:core-4.1-noarch], R 4.0.4

## R CMD check results
There were no ERRORs or WARNINGs.  

There was 1 NOTE when testing in the Mac and CentOS environments:  

```
* checking installed package size ... NOTE
  installed size is  5.9Mb
  sub-directories of 1Mb or more:
    doc       2.4Mb
    extdata   3.4Mb
```
