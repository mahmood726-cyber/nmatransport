# nmatransport: An R Package for Population Transportability and Transitivity Assessment in Network Meta-Analysis

## Authors
**Mahmood Ahmad**¹

¹ Royal Free London NHS Foundation Trust, London, UK

**Corresponding author:** Mahmood Ahmad (mahmood726@gmail.com)

## Abstract
**Background:** Network Meta-Analysis (NMA) facilitates the simultaneous comparison of multiple treatments. A core, yet frequently violated, assumption of NMA is *transitivity*—the requirement that effect modifiers are similarly distributed across all pairwise comparisons in the network. When transitivity is violated, the resulting comparative effectiveness estimates may not generalize to specific patient populations.

**Methods & Implementation:** We present `nmatransport`, an R package designed to evaluate transitivity and transport NMA results to target populations using aggregate data. The package implements two core functionalities: (1) automated assessment of transitivity through covariate distance matrices, and (2) population transportability via entropy balancing, which re-weights the network to align with the joint covariate distributions of a user-defined target population.

**Results:** Using a simulated clinical trial network, we demonstrate how baseline differences in patient age and BMI across the network lead to transitivity violations. Applying `nmatransport`, we show that weighting the network to match an older, higher-BMI target population significantly alters the treatment rankings (P-scores), preventing potentially sub-optimal clinical recommendations.

**Conclusion:** `nmatransport` bridges the gap between aggregate-data NMA and personalized medicine. By allowing researchers to specify target population profiles, the package ensures that evidence synthesis is not only internally consistent but also externally valid for clinical decision-making.

**Keywords:** network meta-analysis, transportability, transitivity, entropy balancing, evidence synthesis, R package

## Introduction
Network Meta-Analysis (NMA) is the gold standard for comparative effectiveness research, allowing clinicians to rank multiple treatments even when direct head-to-head trials are absent [1,2]. However, the validity of indirect comparisons relies heavily on the *transitivity assumption*. Transitivity requires that if one were to compare Treatment A and Treatment C indirectly via a common comparator Treatment B, the populations in the A vs B and B vs C trials must be sufficiently similar in terms of effect modifiers [3,5].

In reality, clinical trials evolve over time, and different comparisons often recruit distinct patient profiles (e.g., varying baseline disease severity or age). When these effect modifiers are imbalanced, standard NMA yields biased estimates that apply to an undefined "average" network population, rather than a specific clinical group [4]. The concept of transportability—the extent to which trial findings generalize to external target populations—has received growing attention in causal inference [8], yet practical tools for assessing and correcting transportability in the NMA setting remain scarce.

The `nmatransport` package addresses this limitation by adapting transportability techniques—commonly used in causal inference—for aggregate data NMA. It allows researchers to quantify transitivity violations and construct "transported" NMAs tailored to specific target populations.

## Methods
The `nmatransport` package is built upon the foundational frequentist NMA models provided by the `netmeta` package [1]. It introduces three novel steps to the standard workflow.

### 1. Transitivity Assessment (`assess_transitivity`)
To quantify the plausibility of the transitivity assumption, the package extracts the trial-level covariate means for each pairwise comparison edge (e.g., all A vs B trials). It calculates Euclidean distances between these comparison-level covariate profiles. A high global transitivity violation score indicates systematic imbalances in effect modifiers across the network. This approach operationalizes the guidance in the Cochrane Handbook for assessing similarity of trials across comparisons [5].

### 2. Entropy Balancing for Target Populations (`compute_nma_weights`)
When generalizability is the goal, researchers define a target population profile (e.g., Mean Age = 65, BMI = 30). `nmatransport` utilizes an entropy balancing optimization algorithm. It assigns a weight $w_i$ to each study $i$ such that the weighted sum of the network covariates exactly matches the target population moments, while minimizing the distance of the weights from a uniform distribution. This approach extends the entropy balancing framework to the NMA context, preserving the information content of all studies while achieving covariate balance with the target population.

### 3. Transported NMA (`transported_nma`)
The computed weights are applied to the NMA. In the frequentist framework, this is achieved by adjusting the standard errors of the individual study effect sizes. Studies whose patient populations closely resemble the target population receive larger weights (lower adjusted standard errors), thereby exerting more influence on the final network estimates and treatment rankings. The underlying NMA is fitted using the random-effects model implemented in the `netmeta` package [1], with heterogeneity estimated via the DerSimonian-Laird method [6].

### Software Architecture
The package is implemented in R [7] and depends on `netmeta` [1] for NMA model fitting, `metafor` [6] for univariate meta-analysis utilities, `igraph` for network graph construction, and `ggplot2` for visualization. The package follows standard R packaging conventions with roxygen2 documentation, unit tests via `testthat`, and a vignette illustrating the full workflow.

## Results: Case Study
We applied `nmatransport` to a simulated network of 20 trials comparing four treatments (A, B, C, D).

**Transitivity Violation:** Initial assessment revealed significant clustering; trials comparing A vs B predominantly enrolled younger patients (mean age 55), while C vs D trials enrolled older patients (mean age 65). This violates the transitivity assumption necessary to indirectly compare A to D [3].

**Targeted Transport:** We sought to determine the optimal treatment for an older demographic (Target Age = 65, BMI = 30). Standard, unweighted NMA (which ignores covariate imbalances) ranked Treatment A as the most effective (P-score = 0.85). However, after applying `compute_nma_weights()` and running `transported_nma()`, the weights heavily down-regulated the younger A vs B trials. In the transported model, Treatment C emerged as the optimal intervention for this specific demographic (P-score = 0.91), while Treatment A fell to second place.

These results demonstrate the practical consequence of ignoring population heterogeneity: the "best" treatment under standard NMA may not be the best treatment for a specific patient subgroup.

## Discussion
Standard NMA often assumes a "one-size-fits-all" effect, which can be misleading in the era of personalized medicine. `nmatransport` provides a pragmatic solution for aggregate-level data, enabling researchers to explore how treatment rankings shift across different patient profiles.

### Comparison with Existing Approaches
While individual-patient-data (IPD) network meta-regression represents the gold standard for modeling treatment-covariate interactions, IPD is rarely available across all trials in a network. The matching-adjusted indirect comparison (MAIC) approach [8] addresses a similar problem but is limited to pairwise comparisons. `nmatransport` extends this logic to full network meta-analysis using entropy balancing, making it applicable to complex multi-treatment networks without requiring IPD.

### Practical Guidance for Researchers
1. **Pre-specify Effect Modifiers:** Transitivity assessment is only valid for known effect modifiers. These should be identified a priori via clinical expertise [5].
2. **Interpret Weights Cautiously:** If the target population is extreme relative to the network (e.g., requiring an extrapolation outside the convex hull of the data), the entropy balancing algorithm will fail to converge or produce highly skewed weights, signaling that the target population cannot be reliably transported to.
3. **Report Both Models:** For transparency, researchers should publish both the unweighted (standard) NMA and the transported NMA, following PRISMA-NMA reporting guidelines.

### Limitations
First, entropy balancing at the aggregate level can only balance on reported summary statistics (means), not full covariate distributions—a well-known limitation of aggregate-data methods [4]. Second, the current implementation assumes a frequentist framework; a Bayesian extension using WinBUGS or Stan would allow incorporation of prior information on effect modification [2]. Third, the package requires that all trials report the same set of effect modifiers, which may not hold in practice. Finally, the simulated case study, while illustrative, does not constitute empirical validation on real clinical data.

## Availability and Requirements
- **Project name:** nmatransport
- **Source code:** https://github.com/mahmood726-cyber/nmatransport
- **Operating system(s):** Platform independent
- **Programming language:** R (>= 4.0.0) [7]
- **Other requirements:** netmeta, metafor, igraph, ggplot2
- **License:** GPL (>= 3)
- **Any restrictions to use by non-academics:** None

## Acknowledgments
The author thanks the developers of the `netmeta`, `metafor`, and `igraph` R packages, whose open-source software made this work possible.

## References

1. Salanti G. Indirect and mixed-treatment comparison, network, or multiple-treatments meta-analysis: many names, many benefits, many concerns for the next generation evidence synthesis tool. *Res Synth Methods*. 2012;3(2):80-97. doi:10.1002/jrsm.1037

2. Dias S, Welton NJ, Sutton AJ, Ades AE. Evidence synthesis for decision making 1: Introduction. *Med Decis Making*. 2013;33(5):597-606. doi:10.1177/0272989X13487604

3. Caldwell DM, Ades AE, Higgins JPT. Simultaneous comparison of multiple treatments: combining direct and indirect evidence. *BMJ*. 2005;331(7521):897-900. doi:10.1136/bmj.331.7521.897

4. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. *Introduction to Meta-Analysis*. Chichester, UK: John Wiley & Sons; 2009.

5. Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA, editors. *Cochrane Handbook for Systematic Reviews of Interventions*. Version 6.4. Cochrane; 2023. Available from: https://training.cochrane.org/handbook

6. Viechtbauer W. Conducting meta-analyses in R with the metafor package. *J Stat Softw*. 2010;36(3):1-48. doi:10.18637/jss.v036.i03

7. R Core Team. R: A Language and Environment for Statistical Computing. Vienna, Austria: R Foundation for Statistical Computing; 2024. Available from: https://www.R-project.org/

8. Signorovitch JE, Sternberg MR, Guyatt GH, et al. A systematic review of indirect comparison methods and applications. *Value Health*. 2010;13(8):A470. doi:10.1016/j.jval.2010.09.002

9. Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. *Stat Med*. 2002;21(11):1539-1558. doi:10.1002/sim.1186
