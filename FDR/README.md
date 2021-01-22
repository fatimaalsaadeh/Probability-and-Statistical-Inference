1.Sort P values

2.Count tests (m)

3.Set Q

4.Plot sorted Pvalues(smallest to largest) vs line Q*c(1:m)/m

5.(If not independent, Q*c(1:m)/(m*(sum(1/i)i=1,…m)))

6.Find P* =largest P value <line

7.Every P<=P* is “interesting”

8.Assignment 8: code this to find all interesting tests from an input vector of P values. Output includes list of hypothesis which are interesting by number in the original unsorted list of p values. input is q values if data is independent and the vector of original p values 

Test on a test vector v1<-c(1e-5*runif(100),runif(900)) use Q=0.05,
