<div class="container-fluid main-container">

<div id="header">

<div class="btn-group pull-right float-right">

Code <span class="caret"></span>

-   <a href="#" id="rmd-show-all-code">Show All Code</a>
-   <a href="#" id="rmd-hide-all-code">Hide All Code</a>
-   
-   <a href="#" id="rmd-download-source">Download Rmd</a>

</div>

# Figure 1.Luminal Study Groups

</div>

<div id="r-markdown" class="section level2">

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. For more details on
using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that
includes both content as well as the output of any embedded R code
chunks within the document. You can embed an R code chunk like this:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

</div>

<div id="including-plots" class="section level2">

## Including Plots

You can also embed plots, for example:

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.

</div>

<div id="rmd-source-code">

LS0tCnRpdGxlOiAiRmlndXJlIDEuTHVtaW5hbCBTdHVkeSBHcm91cHMiCmdlb21ldHJ5OiAibGVmdD0yY20scmlnaHQ9MmNtLHRvcD0yY20sYm90dG9tPTJjbSIKb3V0cHV0OiAKICBwZGZfZG9jdW1lbnQ6CiAgICBmaWdfY2FwdGlvbjogeWVzCiAgICBlY2hvOiBubwogICAgZGV2OgogICAgICAtIHBkZgogICAgICAtIHBuZwpoZWFkZXItaW5jbHVkZXM6IAotIFx1c2VwYWNrYWdle2Zsb2F0fQplZGl0b3Jfb3B0aW9uczogCiAgCiAgY2h1bmtfb3V0cHV0X3R5cGU6IGNvbnNvbGUKa25pdDogKGZ1bmN0aW9uKGlucHV0RmlsZSwgZW5jb2RpbmcsIC4uLikgewogICAgc291cmNlKCIuLi9iaW4vY3VzdG9tX2tuaXRfZnVuY3Rpb25zLlIiKTsKICAgIGtuaXRfZ2l0aHViKGlucHV0RmlsZSkKICAgIH0pCi0tLQoKYGBge3Igc2V0dXAsIGluY2x1ZGU9RkFMU0V9CmlmKGtuaXRyOjppc19sYXRleF9vdXRwdXQoKSl7CiAgcHJpbnRfc2Vzc2lvbl9pbmZvID0gRkFMU0UKICBrbml0cjo6b3B0c19jaHVuayRzZXQoCiAgICBlY2hvPUZBTFNFLAogICAgZGV2PWMoJ3BkZicsICdwbmcnKSwKICAgIGRwaSA9IDMwMCkKfQoKaWYoa25pdHI6OmlzX2h0bWxfb3V0cHV0KCkpewogIHByaW50X3Nlc3Npb25faW5mbyA9IFRSVUUKICBrbml0cjo6b3B0c19jaHVuayRzZXQoCiAgICBlY2hvPVRSVUUsCiAgICBkZXY9YygncGRmJywgJ3BuZycpLAogICAgdG9jID0gRkFMU0UsCiAgICB0b2NfZmxvYXQgPSBGQUxTRSwKICAgIGNvZGVfZm9sZGluZyA9ICJzaG93IiwKICAgIGNvZGVfZG93bmxvYWQgPSBUUlVFKQp9Cgprbml0cjo6b3B0c19jaHVuayRzZXQoCiAgZmlnLndpZHRoID0gNi42OTI5MTMzODU4LAogIGZpZy5wYXRoPSIuL0ZpZ3VyZXMvIiwKICBmaWcuYWxpZ24gPSAiY2VudGVyIiwKICAjZGV2PWMoJ3BkZicsICdwbmcnKSwgZHBpPTMwMCwKICBmaWcucHJvY2VzcyA9IGZ1bmN0aW9uKGZpbGVuYW1lKXsKICAgIG5ld19maWxlbmFtZSA8LSBzdHJpbmdyOjpzdHJfcmVtb3ZlKHN0cmluZyA9IGZpbGVuYW1lLCAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHBhdHRlcm4gPSAiLTEiKQogICAgZnM6OmZpbGVfbW92ZShwYXRoID0gZmlsZW5hbWUsIG5ld19wYXRoID0gbmV3X2ZpbGVuYW1lKQogICAgaWZlbHNlKGZzOjpmaWxlX2V4aXN0cyhuZXdfZmlsZW5hbWUpLCBuZXdfZmlsZW5hbWUsIGZpbGVuYW1lKQogICAgCiAgICAgICMgQ3JlYXRlIFBORwogICAgICAjIHBhdGgyIDwtIHhmdW46OndpdGhfZXh0KG5ld19maWxlbmFtZSwgInBuZyIpCiAgICAgICMgaW1nIDwtIG1hZ2ljazo6aW1hZ2VfcmVhZChuZXdfZmlsZW5hbWUpCiAgICAgICMgbWFnaWNrOjppbWFnZV93cml0ZShpbWcsIHBhdGgyLCBmb3JtYXQgPSAicG5nIikKfSkKIyBzZXR3ZCgiL1VzZXJzL3ZpbGthbC93b3JrL0Jyb2xpZGVuc193b3JrL1Byb2plY3RzL0hJVl9wcm9qZWN0L3JlcG9ydHMvbWFudXNjcmlwdCIpCgojIyMjIEZpZ3VyZSByZXF1aXJlbWVudHMgTWljcm9iaW9tZTogIyMjIwogICMgLSB3aWR0aCBvZiA4NSBtbSBmb3IgaGFsZiBwYWdlIHdpZHRoIGZpZ3VyZQogICMgLSB3aWR0aCBvZiAxNzAgbW0gZm9yIGZ1bGwgcGFnZSB3aWR0aCBmaWd1cmUKICAjIC0gbWF4aW11bSBoZWlnaHQgb2YgMjI1IG1tIGZvciBmaWd1cmUgYW5kIGxlZ2VuZAogICMgLSBpbWFnZSByZXNvbHV0aW9uIG9mIGFwcHJveGltYXRlbHkgMzAwIGRwaSAoZG90cyBwZXIgaW5jaCkgYXQgdGhlIGZpbmFsIHNpemUKYGBgCgojIyBSIE1hcmtkb3duCgpUaGlzIGlzIGFuIFIgTWFya2Rvd24gZG9jdW1lbnQuIE1hcmtkb3duIGlzIGEgc2ltcGxlIGZvcm1hdHRpbmcgc3ludGF4IGZvciBhdXRob3JpbmcgSFRNTCwgUERGLCBhbmQgTVMgV29yZCBkb2N1bWVudHMuIEZvciBtb3JlIGRldGFpbHMgb24gdXNpbmcgUiBNYXJrZG93biBzZWUgPGh0dHA6Ly9ybWFya2Rvd24ucnN0dWRpby5jb20+LgoKV2hlbiB5b3UgY2xpY2sgdGhlICoqS25pdCoqIGJ1dHRvbiBhIGRvY3VtZW50IHdpbGwgYmUgZ2VuZXJhdGVkIHRoYXQgaW5jbHVkZXMgYm90aCBjb250ZW50IGFzIHdlbGwgYXMgdGhlIG91dHB1dCBvZiBhbnkgZW1iZWRkZWQgUiBjb2RlIGNodW5rcyB3aXRoaW4gdGhlIGRvY3VtZW50LiBZb3UgY2FuIGVtYmVkIGFuIFIgY29kZSBjaHVuayBsaWtlIHRoaXM6CgpgYGB7ciBjYXJzfQpzdW1tYXJ5KGNhcnMpCmBgYAoKIyMgSW5jbHVkaW5nIFBsb3RzCgpZb3UgY2FuIGFsc28gZW1iZWQgcGxvdHMsIGZvciBleGFtcGxlOgoKYGBge3IgcHJlc3N1cmUsIGVjaG89RkFMU0V9CnBsb3QocHJlc3N1cmUpCmBgYAoKTm90ZSB0aGF0IHRoZSBgZWNobyA9IEZBTFNFYCBwYXJhbWV0ZXIgd2FzIGFkZGVkIHRvIHRoZSBjb2RlIGNodW5rIHRvIHByZXZlbnQgcHJpbnRpbmcgb2YgdGhlIFIgY29kZSB0aGF0IGdlbmVyYXRlZCB0aGUgcGxvdC4K

</div>

</div>
