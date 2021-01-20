## Resubmission
This is a resubmission.

There was one comment to address:

> Found the following (possibly) invalid URLs:  
>      URL: https://www.rstudio.com/ide/docs/packages/prerequisites (moved to https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites)
>        From: README.md  
>        Status: 200  
>        Message: OK  
> 
> Please change http --> https, add trailing slashes, or follow moved
content as appropriate.

The offending URLs have been removed from the README, as they are not 
essential for CRAN. They merely provide advice on using the development 
version.

## Comments

First release on CRAN.

Note: this package implements a model from a paper that is currently under 
review. There is no DOI or arXiv link at this time. The package will be 
updated with citation information and a DOI in the Description field upon 
publication.

## Test environments
* Windows 10, R 4.0.3
* Fedora 33, R 4.0.3
* winbuilder (devel)

## R CMD check results

0 errors v | 0 warnings v | 2 notes x

There are 2 notes:

> * checking CRAN incoming feasibility ... NOTE
> Maintainer: ‘Neil Marchant <ngmarchant@gmail.com>’
> 
> New submission
> 
> Possibly mis-spelled words in DESCRIPTION:  
>   Marchant (22:8)  
>   al (22:20)  
>   et (22:17)  

The words are not mis-spelled

> The Date field is over a month old.

The package was ready for submission prior to the Christmas/NY break.

> * checking installed package size ... NOTE  
>   installed size is  9.8Mb  
>   sub-directories of 1Mb or more:  
>     libs   9.3Mb  

The large size of the compiled code seems unavoidable

