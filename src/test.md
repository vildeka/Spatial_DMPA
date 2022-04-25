## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. For more details on
using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that
includes both content as well as the output of any embedded R code
chunks within the document. You can embed an R code chunk like this:

``` {.r}
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

``` {.r}
##################
# LOAD LIBRARIES #
##################
library(ggplot2)

#############
# LODA DATA #
#############
#use color brewer as default discrete colors
scale_colour_discrete <- function(...) scale_color_brewer(palette="Set1", ...)
scale_fill_discrete <- function(...) scale_fill_brewer(palette="Set1", ...)

data('mtcars')
# create factors with value labels
mtcars$gear <- factor(mtcars$gear,levels=c(3,4,5),
     labels=c("3gears","4gears","5gears"))
mtcars$am <- factor(mtcars$am,levels=c(0,1),
     labels=c("Automatic","Manual"))
mtcars$cyl <- factor(mtcars$cyl,levels=c(4,6,8),
   labels=c("4cyl","6cyl","8cyl"))
head(mtcars)
```

    ##                    mpg  cyl disp  hp drat    wt  qsec vs        am   gear carb
    ## Mazda RX4         21.0 6cyl  160 110 3.90 2.620 16.46  0    Manual 4gears    4
    ## Mazda RX4 Wag     21.0 6cyl  160 110 3.90 2.875 17.02  0    Manual 4gears    4
    ## Datsun 710        22.8 4cyl  108  93 3.85 2.320 18.61  1    Manual 4gears    1
    ## Hornet 4 Drive    21.4 6cyl  258 110 3.08 3.215 19.44  1 Automatic 3gears    1
    ## Hornet Sportabout 18.7 8cyl  360 175 3.15 3.440 17.02  0 Automatic 3gears    2
    ## Valiant           18.1 6cyl  225 105 2.76 3.460 20.22  1 Automatic 3gears    1

``` {.r}
#################
# COLOR PALETTS #
#################
```

## Including Plots

You can also embed plots, for example:

`<img src="./Figures/pressure.png" style="display: block; margin: auto;" />`{=html}

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
