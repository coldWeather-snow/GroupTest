# Optimal Group Testing Design

2025-04-19

## Description

This repo holds the code to reproduce the manuscript entitled **Single
and multi-objective optimal designs for group testing experiments** for
computing optimal group testing design

## File Explain

Under the [Code/](Code/) directory, you will find the the files that are
used to reproduce the Tables and Figures in the manuscript.

- **Table1.m**: used to generate Table 1 in the manuscript. It is about
  single objective optimal approximate design with $M=61$ and
  $\mathbf{\theta}=(0.07,0.093,0.96)$. Table 2 can be easily obtained by
  changing the value of $M$ to $150$.

- **Table3.m**: used to generate Table 3 in the manuscript. It is about
  multi-objective optimal approximate design with $M=61,150$ and
  $\mathbf{\theta}=(0.07,0.093,0.96)$.

- **Table4.m**: used to generate Table 4 in the manuscript. It is about
  single objective exact design with $M=61,150$,
  $\mathbf{\theta}=(0.07,0.093,0.96)$ with different budget
  $C=100,500,10000$ for the cost situation.

- **Table5.m**: used to generate Table 5 in the manuscript. It is about
  single objective exact design with $M=61,150$,
  $\mathbf{\theta}=(0.07,0.093,0.96)$ with different budget
  $C=100,500,10000$ for the cost-free situation.

- **Table6_budget.m** and **Table6_budgetfree.m**: used to generate
  Table 6 in the manuscript. It is about multi-objective exact design
  with $M=61,150$, $\mathbf{\theta}=(0.07,0.093,0.96)$ with different
  small budget $C=100,500$ or run $n=10,25,50$ cost and cost-free
  situations, respectively.

- **Fig1.m**: used to generate Figure 1 in the manuscript. It is the
  equivalence theorem plot for $D$-$A$ criterion with $M=61$ and
  $q=0.8$.

- **Fig2.m**: used to generate Figure 2 in the manuscript. It is the
  equivalence theorem plot for $D$-$A$-$D_s$ criterion with $M=150$ and
  $q=0.2$.

## References

- Huang, S.-H., Lo Huang, M.-N., Shedeen, K. & Wong, W.K. (2017).  
  [Optimal group testing designs for estimating prevalence with
  uncertain testing
  errors](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssb.12223).
  *JRSSB*, **79**, 1547–1563.

- Huang, S.-H., Huang, M.-N.L. & Shedeen, K. (2021).  
  [Optimal group testing designs for prevalence estimation combining
  imperfect and gold standard
  assays](https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-15/issue-1/Optimal-group-testing-designs-for-prevalence-estimation-combining-imperfect-and/10.1214/20-EJS1786.full).
  *Electronic Journal of Statistics*, **15**, 630–649.

- Huang, S.-H., Huang, M.-N.L. & Shedeen, K. (2020).  
  [Cost considerations for efficient group testing
  studies](https://www.jstor.org/stable/pdf/26892784.pdf). *Statistica
  Sinica*, **30**, 285–302.
