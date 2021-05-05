## Purpose
Fixes missing rmarkdown dependency issue which was brought to my attention by
Kurt Hornik in a 4/16/21 email.


## Test environments
I tested on the default environments specified in
rhub::check_for_cran  (Windows 2008 server, Ubuntu, Fedora)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking Rd cross-references ... NOTE
  Package unavailable to check Rd xrefs: ‘competitiontoolbox’

The competitiontoolbox package is available on CRAN. The link also works
on my local install and running devtools::check(cran=TRUE) does not produce
this note.

