Right now I have two versions of permutation test.
Version 1 is fast but the result may not be accurate, Version2 is very slow, for permutation=1e4 it takes 2 hours to test a single gene and for permutation=1e5 it takes more than one day.

I am interested to see if version1 and version2 gives the similar result and if so we could probably use version1.

In the github, I provide result for score test for all genes and permutation test version1 for all genes and permutation test version2 for the first 50 genes(becasues it takes long time to run).

I perform a two sample t-test on SNP to compare if the p-values are similar 
The result is that 
score test vs permutation test version1    significantly different  with  p=3.54e-140
score test vs permutation test version2    significantly different  with  p=4.0182e-08
permutation test1 and permutation test2    significantly different  with  p=9.4527e-11

My conclusion up to now : permutation test version 1 and 2 are not similar and thus permutation test version1 result is not reliable.

Although score test and permutation test version2 show different result too, it could be the reason that #permutations is very low so the smallest p-value we can reach is  9.9900e-04. We need to solve this by speeding up the permutation test version2.