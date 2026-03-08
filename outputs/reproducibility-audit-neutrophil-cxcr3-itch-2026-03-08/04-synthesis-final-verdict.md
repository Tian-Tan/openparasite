# SYNTHESIS AGENT — FINAL REPRODUCIBILITY VERDICT

**Paper:** Walsh et al., "Neutrophils promote CXCR3-dependent itch in the development of atopic dermatitis"
**Audit date:** 2026-03-08

---

## SECTION 1 — CONSENSUS FINDINGS

All three analysts reached agreement on the following points:

**Where the paper reproduced well:**

1. **Neutrophil peak at Day 5.** All three confirmed the flow cytometry data showing neutrophil percentage at Day 5 (~15%) exceeds Day 8 (~9%) and Day 2 (~4%), exactly matching the paper's directional claim.

2. **NHEK DEG count.** All three confirmed 84 DEGs at padj<0.05, matching the paper's stated value (Analyst C noted a minor discrepancy of 86 in one supplementary file, but the primary figure of 84 was confirmed).

3. **IL8 and CXCL2 log2FC scale mismatch.** All three independently identified that the paper's reported values (IL8: 6.17; CXCL2: 5.04) correspond to log10-scale fold changes, not log2. The conversion factor 1/log10(2) ≈ 3.322 explains the discrepancy with the data file values (1.858 and 1.606). All three classified this as a documentation/labeling error, not data fabrication.

4. **Chloroquine null result.** All three reproduced the near-exact t-statistic (paper: t=0.0964; independent calculations: t≈0.096–0.10), confirming AMG487 had no effect in the chloroquine group. This is the single most precisely reproduced number in the dataset.

5. **AMG487 behavioral significance at Days 8 and 12.** All three found that the AMG487 versus vehicle effect at Days 8 and 12 is directionally significant and robust across multiple analytical approaches.

6. **CXCL10 TSLPRKO result.** All three produced p-values in the range 0.031–0.039, consistent with the paper's reported p=0.0357.

7. **CXCL10 aGR1 null result.** All three found this comparison non-significant, consistent with the paper's direction.

8. **Anti-Gr1 79% reduction and TSLPRKO 44% reduction are unverifiable.** All three independently concluded that individual-level behavioral data for these two headline claims were not supplied.

9. **AMG487 F-statistic cannot be matched.** All three noted that the two-way ANOVA F-statistic reported in the paper (F(1,67)=50.64) cannot be recovered from the provided data, implying approximately 20 mice are absent from the shared dataset.

10. **Permutation test p-values directionally match** (all p<0.05 across four gene groups), though exact values and the reported group medians raise concerns addressed in Section 2.

11. **CXCL10 WT comparison is marginal and sample-incomplete.** All three found that one MC903 animal value is missing, and the reproduced p-value (~0.08–0.09 with n=4) does not match the paper's p=0.029 (n=4,5).

---

## SECTION 2 — DIVERGENCES

**Divergence 1: CXCL10 aGR1 p-value (Analyst A: p=0.771; Analyst B: p≈0.77; paper: p=0.43)**
Data incompleteness issue, not analytical error. No genuine disagreement among analysts.

**Divergence 2: AMG487 Day 5 significance**
Analyst B found that Mann-Whitney yields p≈0.09 (not significant), while parametric tests show p=0.0171. Genuine analytical choice divergence. Most likely explanation: Day 5 effect is real but borderline, sensitive to test choice at small n.

**Divergence 3: TG permutation test reported medians**
Analyst B specifically identified that for 2 of 4 gene groups, the paper's "median absolute log2FC" matches individual gene values (Il31ra=1.299 reported as Excitability median; Nptx2=1.679 reported as Neuroinflammation median) rather than computed group medians (0.911 and 2.347 respectively). Analysts A and C did not flag this. Analyst B's finding is numerically supported. Most likely explanation: reporting error in figure legend.

**Divergence 4: IHC pseudoreplication**
Analyst C explicitly flagged image-level n with only 2 animals per group as pseudoreplication invalidating the axon density p-values. Analysts A and B did not audit this. Analyst C's concern is substantively correct.

**Divergence 5: Lipidomics LTB4**
Analyst B found LTB4 marginally non-significant by Welch t-test (p≈0.052); Analysts A and C did not specifically address this. Genuine borderline result.

**Divergence 6: Overall score**
Analyst A: 3.5/5; Analysts B and C: 3/5. Difference reflects Analyst A's greater weight on directional matches.

---

## SECTION 3 — SCORECARD

| Dimension | Analyst Average | Adjusted Score |
|---|---|---|
| Numerical Reproducibility | 3.33 | **3** |
| Methodological Clarity | 2.83 | **2** |
| Robustness | 3.17 | **3** |
| Statistical Validity | 2.67 | **2** |
| Data Transparency | 2.83 | **2** |

**Numerical Reproducibility — 3/5:** The chloroquine t-statistic is a rare exact match. However, the AMG487 F-statistic cannot be reproduced (~20 missing mice), CXCL10 WT p-value does not reproduce, IL8/CXCL2 log2FC values are systematically mislabeled, and two headline percentage claims (79%, 44%) have no verifiable data.

**Methodological Clarity — 2/5:** The paper does not specify Welch vs Student's t-test; IHC uses image-level n with n=2 animals (pseudoreplication, undisclosed); behavioral data were pooled across cohorts creating a discrepancy between stated n=5–16 and actual pooled n up to 69; the log2FC labeling is erroneous; permutation test medians are inconsistently reported.

**Robustness — 3/5:** Core behavioral claims (AMG487 Day 8/12, CQ null) are robust. Neutrophil kinetics and lipidomics are robust to outlier removal. However, CXCL10 ELISA claims fail Bonferroni and BH-FDR correction; AMG487 Day 5 loses significance under Mann-Whitney; IHC hyperinnervation claims are not robust at animal-level n=2.

**Statistical Validity — 2/5:** IHC pseudoreplication (n=2 animals per group, most serious violation); heteroscedasticity in AMG487 Day 8 (SD ratio ~10:1); normality violations in behavioral data; no multiple comparison correction for three CXCL10 contrasts.

**Data Transparency — 2/5:** Individual-level data for anti-Gr1 and TSLPRKO behavioral experiments absent; ~20 AMG487 mice absent; one CXCL10 WT value missing; raw RNA-seq count matrices not provided; DNFB vehicle supplementary data contains copy-paste error; animal-level IHC data not disclosed.

---

## SECTION 4 — FAILURE MODE ANALYSIS

### Methodological Clarity (scored 2)
**Required changes:**
- Specify Welch vs Student's t-test and justify choice for each analysis
- Explicitly state that IHC uses image-level n as technical replicates within n=X biological animals
- Document cohort pooling strategy and any clustering correction applied
- State the scale of all reported fold-change values consistently (log2 throughout, correcting the IL8/CXCL2 values to 1.858 and 1.606)
- Clarify whether permutation test group summaries are medians or representative gene values

**Concrete fix:** A one-paragraph addition to the statistics subsection: "Welch's t-test was used for all two-group comparisons due to unequal variances. IHC analyses used image-level measurements as technical replicates within n=X biological animals per group; animal-level means were used for statistical comparisons. Behavioral cohorts from separate sessions were pooled after verifying absence of cohort effects (Supplementary file 2). All fold-change values are reported in log2 units."

### Statistical Validity (scored 2)
**Required changes:**
1. **IHC pseudoreplication:** Re-analyze at animal level with n≥5–6 per group, or present as descriptive/pilot data without p-values with explicit n=2 limitation disclosed.
2. **CXCL10 multiple comparisons:** Apply Benjamini-Hochberg FDR correction across three ELISA contrasts, or pre-specify hypotheses to justify uncorrected testing.
3. **Behavioral ANOVA heteroscedasticity:** Report Levene's or Bartlett's test results; where violated (Day 8, SD ratio ~10:1), use Welch-type ANOVA or log-transform scratch counts before ANOVA.

### Data Transparency (scored 2)
**Required additions to supplementary materials:**
1. Individual-level scratch count time series for anti-Gr1 depletion experiment (enabling verification of 79% reduction)
2. Individual-level scratch count time series for TSLPRKO behavioral experiment (enabling verification of 44% reduction)
3. Complete AMG487 behavioral dataset including all cohorts contributing to F(1,67)=50.64 (currently ~20 mice absent)
4. The fifth CXCL10 MC903 WT ELISA value
5. Raw RNA-seq count matrices for NHEK experiment at GEO (already deposited at GSE132174 — link should be explicitly provided in source data section)
6. Animal-level IHC summary values (mean branches/mm² per animal) alongside image-level data
7. Corrected Supplementary 3 DNFB vehicle entry

---

## SECTION 5 — FINAL VERDICT

**PARTIALLY REPRODUCIBLE**

The paper's central behavioral pharmacology findings — that AMG487 reduces scratching at Days 8 and 12, and that chloroquine abrogates this effect — are numerically reproducible, directionally robust across parametric and non-parametric tests, and supported by large effect sizes (Cohen's d > 1.5). The neutrophil kinetics, the NHEK DEG count, and the lipidomics 5-HETE, PGE2, and PGD2 findings reproduce well and are robust to outlier removal. The near-exact match on the chloroquine t-statistic (paper: 0.0964; independent: 0.096) demonstrates that at least some analyses are precisely documented. However, several secondary claims central to the paper's mechanistic interpretation cannot be reproduced or are compromised: the two headline behavioral reductions (79% anti-Gr1; 44% TSLPRKO) lack individual-level data entirely; the CXCL10 ELISA mechanistic evidence fails multiple-comparison correction with n=4 and the WT result does not reproduce due to a missing animal; the IHC hyperinnervation claims rest on image-level pseudoreplication with only n=2 animals per group, invalidating any parametric inference at the biological level; the IL8 and CXCL2 fold-change values are labeled on the wrong scale in the paper; and the permutation test reports individual gene values rather than group medians for two of four gene categories. The biological story the paper tells is likely directionally correct, but the quantitative and statistical scaffolding supporting its mechanistic specificity is insufficiently rigorous and incompletely documented to meet the standard of full reproducibility.
