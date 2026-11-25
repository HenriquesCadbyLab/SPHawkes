## Submission

This is a new submission.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

## Test environments

* local: Ubuntu 24.04, R 4.3.3, `R CMD check --as-cran`, plus a check with
  only the packages in Depends and Imports installed
* GitHub Actions, `R CMD check --as-cran`: macOS (R release), Windows
  (R release), Ubuntu (R release and R oldrel-1)
* win-builder: R-devel (results to be added before submitting)

## Notes for the reviewer

* spHawkes started from the R code accompanying D'Angelo et al. (2024), but
  has since been rewritten, with bug fixes and new methods, so its authors
  are not listed as contributors. The README and package documentation
  acknowledge the original code.
* Words flagged as possibly misspelled in DESCRIPTION (e.g. Hawkes, Mateu,
  Zhuang, declustering) are names or established technical terms.
