# Statistical Analysis Results for 'Groups'

Analysis performed on: 2025-07-09 15:10:42.288291

```R
Performing Kruskal-Wallis tests for Groups 
[1] "Observed Richness (Kruskal-Wallis): p-value = 0.00073971697080019"
Observed Richness (Kruskal-Wallis): p-value = 0.00073971697080019 
  Performing Dunn's post-hoc test for Observed Richness:
  Kruskal-Wallis rank sum test

data: x and groups
Kruskal-Wallis chi-squared = 14.4185, df = 2, p-value = 0


                           Comparison of x by groups                           
                                 (Bonferroni)                                  
Col Mean-|
Row Mean |         CD    Healthy
---------+----------------------
 Healthy |  -3.675457
         |    0.0004*
         |
      UC |  -0.867958   2.527285
         |     0.5781    0.0172*

alpha = 0.05
Reject Ho if p <= alpha/2 
[1] "Chao1 Richness (Kruskal-Wallis): p-value = 0.0194795861896877"
Chao1 Richness (Kruskal-Wallis): p-value = 0.0194795861896877 
  Performing Dunn's post-hoc test for Chao1 Richness:
  Kruskal-Wallis rank sum test

data: x and groups
Kruskal-Wallis chi-squared = 7.8768, df = 2, p-value = 0.02


                           Comparison of x by groups                           
                                 (Bonferroni)                                  
Col Mean-|
Row Mean |         CD    Healthy
---------+----------------------
 Healthy |  -2.792998
         |    0.0078*
         |
      UC |  -1.544718   1.030681
         |     0.1836     0.4540

alpha = 0.05
Reject Ho if p <= alpha/2 
[1] "ACE Richness (Kruskal-Wallis): p-value = 0.013118205807352"
ACE Richness (Kruskal-Wallis): p-value = 0.013118205807352 
  Performing Dunn's post-hoc test for ACE Richness:
  Kruskal-Wallis rank sum test

data: x and groups
Kruskal-Wallis chi-squared = 8.6675, df = 2, p-value = 0.01


                           Comparison of x by groups                           
                                 (Bonferroni)                                  
Col Mean-|
Row Mean |         CD    Healthy
---------+----------------------
 Healthy |  -2.940436
         |    0.0049*
         |
      UC |  -1.240263   1.473121
         |     0.3223     0.2111

alpha = 0.05
Reject Ho if p <= alpha/2 
```

--- End of Analysis ---


Here's a short explanation of the results for each richness variable:

* **Observed Richness:**
    * The Kruskal-Wallis test was highly significant (p < 0.001), indicating overall differences among the "Groups."
    * Dunn's post-hoc test shows that **Healthy** is significantly different from **CD** (p = 0.0004), and **UC** is significantly different from **Healthy** (p = 0.0172). 
      There was no significant difference between UC and CD.

* **Chao1 Richness:**
    * The Kruskal-Wallis test was significant (p = 0.019), suggesting overall differences.
    * Dunn's post-hoc test indicates that **Healthy** is significantly different from **CD** (p = 0.0078).
      No other pairwise comparisons (UC vs CD, UC vs Healthy) showed significant differences.

* **ACE Richness:**
    * The Kruskal-Wallis test was significant (p = 0.013), suggesting overall differences.
    * Dunn's post-hoc test reveals that **Healthy** is significantly different from **CD** (p = 0.0049). 
      No other pairwise comparisons (UC vs CD, UC vs Healthy) showed significant differences.
