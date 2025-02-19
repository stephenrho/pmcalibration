# pmcalibration 0.1.0

* Initial CRAN submission.

# pmcalibration 0.1.1

* Fixed bug when smooth = "none"
* predictions (p) that are <= 0 or >= 1 are now excluded and produce a warning
* pmcalibration output now includes 'p_c' slot containing the value of the calibration curve for each observation
* changed default smooth to "gam"
* improved plot function to show risk distribution (for binary outcomes)
* pmcalibration produces plot by default
* moved vignette material to readme
