Class 19: Investigating Pertussis Resurgence
================
Kira

## 1. Investigating Pertussis Cases by Year:

The CDC tracks cases of Pertussis in the US. We can get their data via
web-scraping.

> Q1. With the help of the R “addin” package datapasta assign the CDC
> pertussis case number data to a data frame called cdc and use ggplot
> to make a plot of cases numbers over time:

Now we can make a plot of cases vs. year:

``` r
library(ggplot2)
baseplot <- ggplot(cdc) + aes(Year,Cases) + geom_point() + geom_line() + labs(title="Pertussis Cases by Year 1922-2019", subtitle="Data from the CDC") 
baseplot
```

![](class19_pertussis_files/figure-commonmark/unnamed-chunk-2-1.png)

## 2. A Tale of Two Vaccines (wP & aP):

> Q2. Using the ggplot geom_vline() function add lines to your previous
> plot for the 1946 introduction of the wP vaccine and the 1996 switch
> to aP vaccine (see example in the hint below). What do you notice?:

``` r
baseplot + geom_vline(xintercept=1946, col="purple", linetype=2) + geom_vline(xintercept=1996, col="orange", linetype=2)
```

![](class19_pertussis_files/figure-commonmark/unnamed-chunk-3-1.png)

Cases dropped after the initial vaccine was introduced but have begun to
increase again recently.

> Q3. Describe what happened after the introduction of the aP vaccine?
> Do you have a possible explanation for the observed trend?

After the introduction of the aP vaccine, cases of pertussis remained
low for a short time but then began to rise again. A possible
explanation for this is that the bacteria is evolving to evade vaccine
immunity.

## 3. Exploring CMI-PB Data:

The CMI-PB API (like most APIs) sends responses in JSON format. Briefly,
JSON data is formatted as a series of key-value pairs, where a
particular word (“key”) is associated with a particular value.

``` r
# Allows us to read, write and process JSON data
library(jsonlite)
```

``` r
subject <- read_json("https://www.cmi-pb.org/api/subject", simplifyVector = TRUE) 

head(subject)
```

      subject_id infancy_vac biological_sex              ethnicity  race
    1          1          wP         Female Not Hispanic or Latino White
    2          2          wP         Female Not Hispanic or Latino White
    3          3          wP         Female                Unknown White
    4          4          wP           Male Not Hispanic or Latino Asian
    5          5          wP           Male Not Hispanic or Latino Asian
    6          6          wP         Female Not Hispanic or Latino White
      year_of_birth date_of_boost      dataset
    1    1986-01-01    2016-09-12 2020_dataset
    2    1968-01-01    2019-01-28 2020_dataset
    3    1983-01-01    2016-10-10 2020_dataset
    4    1988-01-01    2016-08-29 2020_dataset
    5    1991-01-01    2016-08-29 2020_dataset
    6    1988-01-01    2016-10-10 2020_dataset

> Q4. How may aP and wP infancy vaccinated subjects are in the dataset?

``` r
table(subject$infancy_vac)
```


    aP wP 
    47 49 

There are 47 aP infancy vaccinated individuals and 49 wP infancy
vaccinated individuals in this dataset.

> Q5. How many Male and Female subjects/patients are in the dataset?

``` r
table(subject$biological_sex)
```


    Female   Male 
        66     30 

There are 66 female and 30 male subjects in this dataset.

> Q6. What is the breakdown of race and biological sex (e.g. number of
> Asian females, White males etc…)?

``` r
table(subject$race, subject$biological_sex)
```

                                               
                                                Female Male
      American Indian/Alaska Native                  0    1
      Asian                                         18    9
      Black or African American                      2    0
      More Than One Race                             8    2
      Native Hawaiian or Other Pacific Islander      1    1
      Unknown or Not Reported                       10    4
      White                                         27   13

> Q7. Using this approach determine (i) the average age of wP
> individuals, (ii) the average age of aP individuals; and (iii) are
> they significantly different?

``` r
library(lubridate)
```


    Attaching package: 'lubridate'

    The following objects are masked from 'package:base':

        date, intersect, setdiff, union

Calculate the age in years of all subjects:

``` r
age_days <- today()-ymd(subject$year_of_birth)
age_years <- time_length(age_days, "years")
subject$age <- age_years
```

Now I want to calculate the average age of all individuals:

``` r
mean(subject$age)
```

    [1] 31.05627

Now I can calculate the average age of all aP vs. wP vaccinated
subjects:

``` r
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
ap.age <- filter(subject, infancy_vac == "aP")$age
wp.age <- filter(subject, infancy_vac == "wP")$age
mean(ap.age)
```

    [1] 25.52108

``` r
mean(wp.age)
```

    [1] 36.36553

The average age of aP vaccinated subjects is 25.5 and the average age of
wP vaccinated subjects is 36.3.

``` r
t.test(ap.age, wp.age)
```


        Welch Two Sample t-test

    data:  ap.age and wp.age
    t = -12.092, df = 51.082, p-value < 2.2e-16
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -12.644857  -9.044045
    sample estimates:
    mean of x mean of y 
     25.52108  36.36553 

Since the p-value is less than 0.05, the difference is statistically
significant.

> Q8. Determine the age of all individuals at time of boost?

``` r
int <- ymd(subject$date_of_boost) - ymd(subject$year_of_birth)
age_at_boost <- time_length(int, "year")
head(age_at_boost)
```

    [1] 30.69678 51.07461 33.77413 28.65982 25.65914 28.77481

> Q9. With the help of a faceted boxplot (see below), do you think these
> two groups are significantly different?

``` r
ggplot(subject) +
  aes(age,
      fill=as.factor(infancy_vac)) +
  geom_histogram(show.legend=FALSE) +
  facet_wrap(vars(infancy_vac), nrow=2) 
```

    `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](class19_pertussis_files/figure-commonmark/unnamed-chunk-15-1.png)

## Joining multiple tables:

> Q9. Complete the code to join specimen and subject tables to make a
> new merged data frame containing all specimen records along with their
> associated subject details:

``` r
# Complete the API URLs...
specimen <- read_json("http://cmi-pb.org/api/specimen", simplifyVector = TRUE) 
titer <- read_json("http://cmi-pb.org/api/ab_titer", simplifyVector = TRUE) 
```

``` r
meta <- inner_join(specimen,subject, by=join_by("subject_id"))
dim(meta)
```

    [1] 729  14

> Q10. Now using the same procedure join meta with titer data so we can
> further analyze this data in terms of time of visit aP/wP, male/female
> etc.

``` r
abdata <- inner_join(titer,meta, by=join_by("specimen_id"))
dim(abdata)
```

    [1] 32675    21

> Q11. How many specimens (i.e. entries in abdata) do we have for each
> isotype?

``` r
table(abdata$isotype)
```


     IgE  IgG IgG1 IgG2 IgG3 IgG4 
    6698 1413 6141 6141 6141 6141 

> Q12. What do you notice about the number of visit 8 specimens compared
> to other visits?

``` r
table(abdata$visit)
```


       1    2    3    4    5    6    7    8 
    5795 4640 4640 4640 4640 4320 3920   80 

There are a lot fewer samples for visit 8 compared with the other
visits. This is because the study is still ongoing.

## 4. Examine IgG1 Ab Titer Levels:

Now using our joined/merged/linked abdata dataset filter() for IgG1
isotype and exclude the small number of visit 8 entries.

``` r
ig1 <- abdata %>% filter(isotype == "IgG1", visit!=8)
head(ig1)
```

      specimen_id isotype is_antigen_specific antigen        MFI MFI_normalised
    1           1    IgG1                TRUE     ACT 274.355068      0.6928058
    2           1    IgG1                TRUE     LOS  10.974026      2.1645083
    3           1    IgG1                TRUE   FELD1   1.448796      0.8080941
    4           1    IgG1                TRUE   BETV1   0.100000      1.0000000
    5           1    IgG1                TRUE   LOLP1   0.100000      1.0000000
    6           1    IgG1                TRUE Measles  36.277417      1.6638332
       unit lower_limit_of_detection subject_id actual_day_relative_to_boost
    1 IU/ML                 3.848750          1                           -3
    2 IU/ML                 4.357917          1                           -3
    3 IU/ML                 2.699944          1                           -3
    4 IU/ML                 1.734784          1                           -3
    5 IU/ML                 2.550606          1                           -3
    6 IU/ML                 4.438966          1                           -3
      planned_day_relative_to_boost specimen_type visit infancy_vac biological_sex
    1                             0         Blood     1          wP         Female
    2                             0         Blood     1          wP         Female
    3                             0         Blood     1          wP         Female
    4                             0         Blood     1          wP         Female
    5                             0         Blood     1          wP         Female
    6                             0         Blood     1          wP         Female
                   ethnicity  race year_of_birth date_of_boost      dataset
    1 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    2 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    3 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    4 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    5 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    6 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
           age
    1 37.20192
    2 37.20192
    3 37.20192
    4 37.20192
    5 37.20192
    6 37.20192

> Q13. Complete the following code to make a summary boxplot of Ab titer
> levels for all antigens:

``` r
ggplot(ig1) +
  aes(MFI, antigen) +
  geom_boxplot() + 
  facet_wrap(vars(visit), nrow=2)
```

![](class19_pertussis_files/figure-commonmark/unnamed-chunk-22-1.png)

> Q14. What antigens show differences in the level of IgG1 antibody
> titers recognizing them over time? Why these and not others?

Antigens such as FIM 2/3, PT, and FHA show change over time. These
antibody levels change because they represent components of the
Pertussis vaccine.

We can attempt to examine differences between wP and aP here by setting
color and/or facet values of the plot to include infancy_vac status:

``` r
ggplot(ig1) +
  aes(MFI, antigen, col=infancy_vac ) +
  geom_boxplot(show.legend = FALSE) + 
  facet_wrap(vars(visit), nrow=2) +
  theme_bw()
```

![](class19_pertussis_files/figure-commonmark/unnamed-chunk-23-1.png)

Another version of this plot adding infancy_vac to the faceting:

``` r
ggplot(ig1) +
  aes(MFI, antigen, col=infancy_vac ) +
  geom_boxplot(show.legend = FALSE) + 
  facet_wrap(vars(infancy_vac, visit), nrow=2)
```

![](class19_pertussis_files/figure-commonmark/unnamed-chunk-24-1.png)

> Q15. Filter to pull out only two specific antigens for analysis and
> create a boxplot for each. You can chose any you like. Below I picked
> a “control” antigen (“Measles”, that is not in our vaccines) and a
> clear antigen of interest (“FIM2/3”, extra-cellular fimbriae proteins
> from B. pertussis that participate in substrate attachment):

``` r
filter(ig1, antigen=="Measles") %>%
  ggplot() +
  aes(MFI, col=infancy_vac) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(vars(visit)) +
  theme_bw()
```

![](class19_pertussis_files/figure-commonmark/unnamed-chunk-25-1.png)

``` r
filter(ig1, antigen=="FIM2/3") %>%
  ggplot() +
  aes(MFI, col=infancy_vac) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(vars(visit)) +
  theme_bw()
```

![](class19_pertussis_files/figure-commonmark/unnamed-chunk-26-1.png)

``` r
filter(ig1, antigen=="FHA") %>%
  ggplot() +
  aes(MFI, col=infancy_vac) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(vars(visit)) +
  theme_bw()
```

![](class19_pertussis_files/figure-commonmark/unnamed-chunk-27-1.png)

> Q16. What do you notice about these two antigens time course and the
> FIM2/3 data in particular?

The levels of the antigen control I chose (Measles) did not
significantly change throughout the course of the 7 visits. In
comparison, the levels of FIM2/3 increased throughout the first 5 visits
and then began to decrease again.

> Q17. Do you see any clear difference in aP vs. wP responses?

There is some differences but it’s unclear if this is statistically
significant from the data we have so far.

## 5. Obtaining CMI-PB RNASeq Data:

For RNA-Seq data the API query mechanism quickly hits the web browser
interface limit for file size. We will present alternative download
mechanisms for larger CMI-PB datasets in the next section.

``` r
url <- "https://www.cmi-pb.org/api/v2/rnaseq?versioned_ensembl_gene_id=eq.ENSG00000211896.7"

rna <- read_json(url, simplifyVector = TRUE) 

head(rna)
```

      versioned_ensembl_gene_id specimen_id raw_count      tpm
    1         ENSG00000211896.7         344     18613  929.640
    2         ENSG00000211896.7         243      2011  112.584
    3         ENSG00000211896.7         261      2161  124.759
    4         ENSG00000211896.7         282      2428  138.292
    5         ENSG00000211896.7         345     51963 2946.136
    6         ENSG00000211896.7         244     49652 2356.749

``` r
ssrna <- inner_join(rna, meta)
```

    Joining with `by = join_by(specimen_id)`

> Q18. Make a plot of the time course of gene expression for IGHG1 gene
> (i.e. a plot of visit vs. tpm).

``` r
ggplot(ssrna) +
  aes(visit, tpm, group=subject_id) +
  geom_point() +
  geom_line(alpha=0.2)
```

![](class19_pertussis_files/figure-commonmark/unnamed-chunk-30-1.png)

> Q19. What do you notice about the expression of this gene (i.e. when
> is it at it’s maximum level)?

The expression of this gene is at its maximum level at visit 4.

> Q20. Does this pattern in time match the trend of antibody titer data?
> If not, why not?

This doesn’t exactly match the trends of antibody data which peak around
visit 5 and drop off more slowly. The expression of the IGHG1 gene peaks
at visit 4 and drops off rather quickly.

We can dig deeper and color and/or facet by infancy_vac status:

``` r
ggplot(ssrna) +
  aes(tpm, col=infancy_vac) +
  geom_boxplot() +
  facet_wrap(vars(visit))
```

![](class19_pertussis_files/figure-commonmark/unnamed-chunk-31-1.png)

There is no obvious wP vs. aP differences here even if we focus in on a
particular visit:

``` r
ssrna %>%  
  filter(visit==4) %>% 
  ggplot() +
    aes(tpm, col=infancy_vac) + geom_density() + 
    geom_rug() 
```

![](class19_pertussis_files/figure-commonmark/unnamed-chunk-32-1.png)
