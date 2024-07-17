## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results
❯ On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release)
  checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Martin Law <martin.law@mrc-bsu.cam.ac.uk>'
  
  Found the following (possibly) invalid DOIs:
    DOI: 10.1002/pst.2067
      From: DESCRIPTION
      Status: Service Unavailable
      Message: 503

❯ On windows-x86_64-devel (r-devel)
  checking for detritus in the temp directory ... NOTE

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

* DOI has been checked at doi.org and is valid.
* Detritus in temp directory is file/directory 'lastMiKTeXException'. As noted in R-hub issue #503, this could be due to a bug/crash in MiKTeX and can likely be ignored.
