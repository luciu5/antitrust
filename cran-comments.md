## Purpose
Update maintainer email address


## Test environments
I tested on the default environments specified in
rhub::check_for_cran  (Windows 2008 server, Ubuntu, Fedora)
https://win-builder.r-project.org/

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... [13s] NOTE
Maintainer: 'Charles Taragin <ctaragin+antitrustr@gmail.com>'

New maintainer:
  Charles Taragin <ctaragin+antitrustr@gmail.com>
Old maintainer(s):
  Charles Taragin <ctaragin@ftc.gov>

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1111/1756-2171.12385
    From: man/SupplyChain-Functions.Rd
    Status: 503
    Message: Service Unavailable
  URL: https://www.jstor.org/stable/1806699
    From: inst/doc/Reference.html
    Status: 403
    Message: Forbidden
  URL: https://www.jstor.org/stable/2950522
    From: inst/doc/Reference.html
    Status: 403
    Message: Forbidden



I am switching jobs and will shortly no longer have access to ctaragin@ftc.gov.

I have checked the doi URL and it works, but with a redirect. 
I also checked the two jStor URLs by copying and pasting into two different browsers 
and they appear to work without a redirect.

