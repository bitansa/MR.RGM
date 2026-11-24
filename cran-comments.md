## Submission

This is an update of MR.RGM from 0.1.3, the version on CRAN, to 0.1.4. The changes are listed in NEWS.md.

Version 0.1.4 sets an excluded instrument effect exactly to zero when the instruments are selected (a point-mass
spike in place of a narrow normal spike). This is the model of the accompanying methods paper. Fits with a known instrument map are unchanged. We are sorry for the short interval since
0.1.3. The update is needed so that the version on CRAN is the version described in the paper.

## Test environments

* local macOS (aarch64-apple-darwin23), R 4.6.1
* win-builder, R-devel (2026-11-17 r90926 ucrt)
* win-builder, R-release (4.6.2)

## R CMD check results

0 errors | 0 warnings | 0 notes
