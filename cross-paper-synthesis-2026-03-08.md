# CROSS-PAPER REPRODUCIBILITY SYNTHESIS
**Papers audited:**
- **Paper 1:** Walsh et al., "Neutrophils promote CXCR3-dependent itch in the development of atopic dermatitis," *eLife* 2019;8:e48448. Audit verdict: **PARTIALLY REPRODUCIBLE** (3/5 overall).
- **Paper 2:** Tavares-Ferreira et al., "Single-nucleus transcriptomic analysis of human dorsal root ganglion neurons," *eLife*, GSE168243. Audit verdict: **PARTIALLY REPRODUCIBLE** (2.4/5 overall).

Both audits were conducted 2026-03-08 using a four-stage isolated-analyst protocol (Data Collector, Methodology Extractor, three independent Analysts, and a Synthesis Agent).

---

## PART 1 — CLAIM STABILITY INVENTORY

### PAPER 1: Walsh et al. (atopic dermatitis / CXCR3-itch)

**TIER 1 — STABLE AND VERIFIED**

**P1-C1. Neutrophil infiltration peaks at Day 5 (~15%) and declines by Day 8 (~9%) and Day 2 (~4%).**
All three analysts independently confirmed this from flow cytometry data; the directional ordering reproduces exactly with no analytical ambiguity.

**P1-C2. NHEK cells treated with neutrophil supernatant show 84 differentially expressed genes at padj < 0.05.**
All three analysts confirmed 84 DEGs (Analyst C noted 86 in one supplementary file, but the primary figure was verified). The result is documented in a deposited GEO dataset (GSE132174) and reproduced exactly.

**P1-C3. AMG487 (CXCR3 antagonist) reduces scratching behavior at Days 8 and 12; chloroquine (histamine-dependent itch) is unaffected (t ≈ 0.0964).**
The chloroquine t-statistic is the single most precisely reproduced number in the dataset (paper: 0.0964; independent: 0.096–0.10). AMG487 Days 8 and 12 effects are directionally significant and robust across parametric and non-parametric tests with Cohen's d > 1.5. These are the paper's core pharmacological claims and they hold.

---

**TIER 2 — DIRECTIONALLY SUPPORTED**

**P1-C4. Lipidomics: 5-HETE, PGE2, and PGD2 are elevated in MC903-treated skin.**
Directionally reproduced and robust to outlier removal; no analyst challenged the direction. Exact magnitudes and test details are underspecified (Welch vs. Student's t-test not stated), introducing uncertainty around precise p-values.

**P1-C5. CXCL10 is elevated in TSLPRKO versus WT comparisons (p = 0.0357).**
Reproduced in the range p = 0.031–0.039 across all three analysts, consistent with the stated value. However, this comparison is one of three ELISA contrasts conducted without multiple-comparison correction; under Bonferroni or BH-FDR, it would not survive. The direction is real; the unadjusted p-value should not be taken at face value.

**P1-C6. CXCL10 aGR1 null result (neutrophil-depleted animals do not differ in CXCL10).**
Directionally confirmed by all three analysts. The paper's p-value (0.43) does not reproduce (analysts obtain ~0.77) due to incomplete data, but the direction — non-significance — is consistent and uncontested.

**P1-C7. Permutation test across four neural gene groups: all groups show p < 0.05.**
All four groups are directionally significant across all three analyses. However, Analyst B established that for two of four groups (Excitability and Neuroinflammation), the paper reports individual gene values (Il31ra = 1.299; Nptx2 = 1.679) as group medians rather than computed medians (0.911 and 2.347 respectively). The significance is real; the descriptive statistics are mislabeled.

---

**TIER 3 — CONTESTED**

**P1-C8. IL8 log2FC = 6.17 and CXCL2 log2FC = 5.04 in NHEK RNA-seq.**
All three analysts independently established that these values correspond to log10-scale, not log2-scale fold changes. The conversion factor 1/log10(2) ≈ 3.322 explains the discrepancy; the actual log2 fold-changes are 1.858 (IL8) and 1.606 (CXCL2), representing approximately 3.6-fold and 3.1-fold upregulation respectively — biologically meaningful but far less dramatic than the paper implies. This is a documentation/labeling error, not fabrication, but it directly distorts the reader's interpretation of NHEK activation magnitude.

**P1-C9. CXCL10 WT comparison (p = 0.029).**
Reproduced as p ≈ 0.08–0.09 (n = 4, with one MC903 WT ELISA value missing from the deposited data). The result does not reproduce, and the missing value makes the stated n = 4, 5 impossible to verify. Contested.

**P1-C10. AMG487 F-statistic F(1,67) = 50.64 from two-way ANOVA.**
Cannot be recovered from the deposited data; approximately 20 mice are absent from the shared dataset. The behavioral effect of AMG487 is well-supported (Tier 1), but the specific test statistic and degrees of freedom cannot be independently confirmed.

**P1-C11. AMG487 reduces scratching at Day 5.**
Analyst B found p ≈ 0.09 (Mann-Whitney, not significant at conventional threshold); Analysts A and C using parametric tests found p = 0.0171. Genuine analytical choice divergence. The Day 5 effect, if real, is borderline and sensitive to test choice at small n. The paper's claim of significance here is contested.

**P1-C12. LTB4 elevated in MC903 lipidomics.**
Borderline non-significant by Welch t-test in Analyst B's analysis (p ≈ 0.052); Analysts A and C did not specifically address LTB4. This claim should not be treated as established.

**P1-C13. IHC: increased axon density / hyperinnervation in AD skin.**
Analyst C explicitly identified pseudoreplication: statistical analysis was conducted at the image level with only n = 2 biological animals per group. Any parametric p-values derived from this design are invalid at the animal-comparison level. The directional observation (more nerve fibers visible) may be real, but no valid statistical inference is possible from the reported data.

---

**TIER 4 — UNVERIFIABLE FROM AVAILABLE DATA**

**P1-C14. Anti-Gr1 depletion reduces scratching by 79%.**
Individual-level behavioral time-series data for this experiment were not deposited. The claim cannot be evaluated. This is a documentation failure.

**P1-C15. TSLPRKO mice show 44% reduction in scratching.**
Same failure mode. Individual-level behavioral data absent. Claim cannot be evaluated.

---

### PAPER 2: Tavares-Ferreira et al. (human DRG snRNA-seq)

**TIER 1 — STABLE AND VERIFIED**

**P2-C1. 1,837 neuronal nuclei across six preparations with per-preparation counts of 212/152/770/281/80/342.**
Confirmed exactly and without ambiguity by all three analysts via direct file inspection. This is a data-level fact, independent of any analytical choice. The most robustly reproducible finding in the audit.

**P2-C2. TMEM100 is near-absent in human DRG (<0.5% of cells positive).**
Confirmed across all QC parameterizations (0.33–0.39% across strict, baseline, and lenient filtering). This finding is entirely independent of normalization method, clustering algorithm, or feature selection, as it reflects raw UMI count data.

---

**TIER 2 — DIRECTIONALLY SUPPORTED**

**P2-C3. Human DRG neurons segregate into approximately a dozen transcriptomically distinct classes (H1–H15, ~12–15).**
Silhouette analysis by Analyst B yields optimal k = 12; the qualitative claim of "approximately a dozen" is consistently supported across all analytical parameterizations (PCs 20–50, HVGs 1,000–2,000, normalization methods). However, the mean ARI across parameterizations is 0.547, meaning specific cell-to-cluster assignments change moderately with analytical choices. The count is directionally correct; the precise identities of H1–H15 are unverifiable.

**P2-C4. Multigene ISH with NEFH, TAC1, and OSMR labels essentially every cell.**
95.2% of filtered cells in the snRNA-seq data express at least one of the three probe genes (Analyst B), consistent with the paper's claim. The snRNA-seq data systematically underestimates expression relative to ISH (cytoplasmic vs. nuclear mRNA), so the true ISH coverage is likely higher. Directionally supported; exact coverage unconfirmable from count matrices.

**P2-C5. H5 neurons are CALCA-negative, NTRK1-negative, and NTRK2-positive.**
Global marker expression is quantifiable (NTRK2: 23.44%; CALCA: 49.7%; NTRK1: 30.87%), and the pattern of NTRK2 co-occurring with CALCA absence is visible in the raw data. However, the H5 cluster label cannot be mapped to a specific k-means partition without the original Seurat object. The profile is directionally plausible; the cluster-specific attribution is unconfirmable.

---

**TIER 3 — CONTESTED**

**P2-C6. H10 and H11 together account for approximately 20% of neurons.**
Analyst B's two smallest clusters sum to 4.5% under baseline parameterization (range: 1.2–7.5% across alternatives). The ~20% figure is not reproduced under any tested configuration. The most likely explanation is label-mapping uncertainty (which unlabeled clusters correspond to H10/H11), but the gap is too large to dismiss. Contested.

**P2-C7. S1PR3 is "not detected" in human DRG.**
This claim is directly and consistently contradicted by the raw count data across all three independent analyses: ~40–60 cells express S1PR3 across multiple preparations (1.77% post-filtering, Analyst B; ~6.6% in prep1 pre-filter, Analyst A; ~2–3% overall, Analyst C). The categorical statement "not detected" is incorrect as written. The correct characterization is "detected at very low frequency."

---

**TIER 4 — UNVERIFIABLE FROM AVAILABLE DATA**

**P2-C8. Spatial clustering is statistically significant (p ≤ 6.96 × 10⁻⁴², Mann-Whitney U).**
ISH tissue section XY coordinates were not deposited in GEO (GSE168243). The exact p-value cannot be reproduced. A transcriptomic-space proxy test supports the qualitative direction (Analyst B: 10/40 neighbor-count values significant), but the specific statistical result is unverifiable.

**P2-C9. H9 shows only weak similarity to any mouse neuron class (KL divergence).**
Mouse reference single-cell data were not deposited alongside the human data. KL divergence parameters (gene set, pseudocount, directionality, distribution construction method) are not specified. Entirely untestable from available data.

**P2-C10. No human DRG class shows appreciable similarity to mouse cLTMRs.**
Same failure mode as P2-C9. The comparative biology claim is scientifically important but completely unverifiable without the missing mouse data and KL implementation.

---

## PART 2 — CROSS-PAPER SYNTHESIS

### A. CONVERGING CLAIMS

**Convergence 1: Molecular specificity in itch, not generic inflammation.**
Paper 1 demonstrates that itch in MC903 atopic dermatitis requires a specific chemokine pathway (CXCL10/CXCR3), not a nonspecific inflammatory cascade — the null chloroquine result (P1-C3, Tier 1) is precisely the control demonstrating this. Paper 2 demonstrates that the sensory neurons that would receive such signals are not a monolithic population but approximately a dozen molecularly distinct classes (P2-C3, Tier 2). Together, these papers push toward the same conceptual position: itch is a molecularly specific, circuit-level phenomenon localized to discrete populations, not an emergent property of generalized inflammation or generic nociceptor activation.

The confidence this convergence adds is moderate, not high. Paper 1's pharmacological specificity is well-supported (Tier 1 data), but Paper 2's subtype-resolution claims rest on Tier 2 evidence with unresolved batch effects. The convergence is conceptually coherent but cannot yet be grounded numerically in the direct claim "CXCR3-expressing neurons belong to human DRG subtype X."

**Convergence 2: The OSMR/IL-31 receptor axis as a candidate skin-to-neuron signaling pathway in AD.**
Paper 1's NHEK RNA-seq found Il31ra among the differentially upregulated genes (verified, Tier 1 by NHEK DEG count). IL-31 receptor function requires co-assembly of IL31RA with OSMR. Paper 2 directly quantifies OSMR expression in human DRG neurons: 25.8% of all filtered cells are OSMR-positive. These two findings, from different experimental systems (mouse keratinocytes and human sensory neurons), independently implicate the IL-31/OSMR axis as a plausible conduit through which epithelial signals in atopic dermatitis reach DRG neurons. This is a scientifically meaningful convergence, but an important caveat applies: Paper 2's OSMR data is Tier 2, and neither paper tested IL-31 signaling functionally at the neuron level. The convergence is suggestive, not mechanistically established.

**Convergence 3: Both papers share a common methodological failure mode — pseudoreplication.**
Paper 1 conducts IHC with n = 2 animals per group analyzed at image level. Paper 2 uses data from the same donor twice (preparations 3 and 5, 850/1,837 cells = 46.3%) without acknowledging or correcting for this statistical dependency. This is not a scientific convergence but a convergence in study design weakness that should inform how results from both papers are interpreted: neither paper's statistics fully account for biological replication at the animal or individual level.

---

### B. DIVERGING CLAIMS

**Divergence 1: Mouse versus human DRG neuron type correspondence.**
Paper 1 implicitly treats the mouse DRG as a valid model for human itch circuitry. Paper 2's central scientific contribution, if the cross-species comparison were verifiable, is that the correspondence between mouse and human DRG neuron types is imperfect — H9 shows weak similarity to all mouse classes (Tier 4), and no human class aligns with mouse cLTMRs (Tier 4). The tension is this: Paper 1 builds its mechanistic story entirely on mouse biology; Paper 2 suggests the mouse-to-human translation of DRG subtype biology may be non-trivial.

The most plausible explanation is not contradiction but complementarity at different levels. Pharmacological targets like CXCR3 are conserved proteins; the molecular pathways Paper 1 identifies likely operate in both species even if the subtype landscape differs. However, if human DRG genuinely lacks direct mouse equivalents for some subtypes, the precise cell population responsible for CXCR3-dependent itch in mice may not map cleanly to a human counterpart. This tension cannot be resolved with available data because Paper 2's cross-species claims are Tier 4.

**Divergence 2: The extent to which IL-31 signaling involves specific DRG subtypes.**
Paper 1 implies (via NHEK Il31ra upregulation and TSLPRKO behavioral phenotype) that epithelial IL-31/TSLP signals act on sensory neurons in a pathway-specific way. Paper 2's OSMR data (25.8% of all DRG neurons) suggests the receptor is broadly, not narrowly, distributed — which would predict broad rather than subtype-specific sensitization. This is not a strong contradiction (Il31ra and OSMR have different expression patterns, and co-expression is what matters), but it creates uncertainty about whether the TSLP/OSMR axis is a fine-grained subtype-specific route or a broad sensitization mechanism. Resolving this requires querying both IL31RA and OSMR co-expression within specific DRG subtypes from Paper 2's deposited data — a HIGH feasibility analysis (see Part 3).

---

### C. COMPLEMENTARY GAPS

**Gap 1 (most important): Which human DRG subtype expresses CXCR3?**
Paper 1 establishes CXCL10/CXCR3 as a functional itch axis in mouse (Tier 1). Paper 2 provides an atlas of 1,837 human DRG neuronal nuclei with expression data for 58,188 genes. CXCR3 is one of those genes. Paper 1 cannot answer "which neuron type responds to CXCL10" because it lacks single-cell transcriptomics. Paper 2 could answer "which human DRG subtype expresses CXCR3" but doesn't query it. A combined analysis could close this gap entirely using existing deposited data.

**Gap 2: IL31RA co-expression with OSMR in human DRG subtypes.**
Paper 1 shows IL31RA is upregulated in keratinocytes; Paper 2 shows OSMR is present in 25.8% of human DRG neurons. Neither paper establishes whether both receptors are co-expressed in the same human neuronal subtype. Joint query of IL31RA and OSMR in GSE168243 could identify the specific DRG subtype(s) competent for IL-31 signaling — directly connecting Paper 1's mechanistic inference to a human substrate.

**Gap 3: Mouse DRG CXCR3 expression by subtype.**
Paper 1's pharmacology is all in mouse, but the paper never establishes which mouse DRG neuron type expresses CXCR3. Published mouse DRG single-cell datasets exist (external, MEDIUM feasibility) and could fill this gap, providing the link between Paper 1's behavioral pharmacology and a specific cell-type identity in the same species.

**Gap 4: Longitudinal neuronal changes during AD development.**
Paper 1 has a time-course (Days 2, 5, 8, 12) for behavioral and immune data but no single-cell neuronal profiling. Paper 2 has a detailed snapshot of DRG neuron types but no disease condition or time-course. Together, a time-resolved transcriptomics study of DRG neurons during MC903-induced AD — or during human AD flares — would answer whether itch-related DRG subtypes undergo dynamic transcriptional changes that Paper 1's pharmacology is targeting and that Paper 2's atlas could predict.

---

## PART 3 — FOLLOW-UP ANALYSIS RECOMMENDATIONS

---

**RECOMMENDATION 1 (HIGH feasibility)**

**QUESTION:** Does CXCR3 expression in human DRG neurons localize to a specific transcriptomic subclass, or is it broadly distributed?

**MOTIVATION:** Paper 1's Tier 1 finding establishes CXCL10/CXCR3 as the primary itch-driving axis in mouse AD. The identities of the neurons that receive this signal are unknown. Paper 2's deposited data (GSE168243) contains expression measurements for all 58,188 genes — including CXCR3 — in 1,837 neuronal nuclei. Closing this gap would directly translate Paper 1's pharmacological finding to a specific human cell population.

**WHAT DATA IS NEEDED:** GSE168243 (already available locally in `geo_data/`). Query CXCR3 (gene symbol: CXCR3) across all cells. Cross-tabulate CXCR3 detection with existing cluster assignments — if the authors' Seurat object is unavailable, use any reproducible clustering (e.g., the k=12 configuration validated by Analyst B) as a proxy. Also query CXCL11 (the other CXCR3 ligand) and CCR3, CCR4 for completeness.

**SUGGESTED METHOD:** Single-gene query on the raw count matrices: binarize detection at UMI ≥ 1 (as done for TMEM100 and S1PR3 in the audit). If CXCR3 is detected in < 5% of cells (as S1PR3 was), cluster-level enrichment tests become underpowered; use co-expression with established subtype markers (NTRK1, NTRK2, CALCA, NEFH, TAC1) instead of cluster labels to infer subtype identity. Fisher's exact test for enrichment in marker-positive vs. marker-negative cells.

**FEASIBILITY:** HIGH. All required data are deposited and verified. Can be executed with existing code from Analyst B's pipeline (Step 7 gene expression queries). No new data collection required.

**RISK:** CXCR3 is a low-abundance transcript in snRNA-seq data generally (G-protein-coupled receptor; low nuclear mRNA). It may be below detection in the majority of cells, making subtype-level conclusions impossible. If CXCR3 is detected in fewer than ~20 cells, the cluster enrichment analysis will be underpowered and results will be directional at best.

---

**RECOMMENDATION 2 (HIGH feasibility)**

**QUESTION:** Are IL31RA and OSMR co-expressed in the same human DRG neuronal subtype, and if so, which?

**MOTIVATION:** Both Paper 1 (NHEK Il31ra upregulation, Tier 1 DEG count) and Paper 2 (OSMR present in 25.8% of human DRG neurons, Tier 2) implicate the IL-31 receptor complex in the skin-to-neuron itch pathway. The IL-31 receptor requires IL31RA and OSMR co-assembly. Co-expression in a specific subtype would identify the human DRG cell population most likely targeted by IL-31 signals from inflamed keratinocytes, directly integrating the two papers' findings.

**WHAT DATA IS NEEDED:** GSE168243 (`geo_data/`). Query IL31RA and OSMR expression. OSMR detection rate (25.8%) is confirmed by the audit. IL31RA detection rate was not explicitly reported by any analyst but can be directly extracted.

**SUGGESTED METHOD:** Generate a 2×2 co-expression table (IL31RA positive/negative × OSMR positive/negative) across all 1,837 cells. Chi-squared test for non-independence. Then stratify by neuronal marker (NTRK1, NTRK2, CALCA, TAC1, NEFH, TRPV1) to identify which cell class drives co-expression. This is a straightforward extension of Analyst B's Step 7 queries.

**FEASIBILITY:** HIGH. Data in hand, method requires minimal additional code.

**RISK:** IL31RA may be very sparsely detected in snRNA-seq data (the nuclear mRNA for surface receptors is often near the detection limit). If both IL31RA and OSMR are low-abundance, co-expression will appear rare by chance regardless of biology, and results will not be interpretable. Absence of co-expression in snRNA-seq does not imply absence of protein co-expression.

---

**RECOMMENDATION 3 (HIGH feasibility)**

**QUESTION:** Does batch correction for 10x Genomics v2 versus v3 chemistry alter the cluster identities or proportions reported in Paper 2?

**MOTIVATION:** Paper 2's Tier 3 claim for H10+H11 proportions (~20% vs. the audited 4.5%) and all Tier 2 cluster-identity claims are threatened by the finding that three of fifteen clusters are 100% chemistry-pure (Analyst B, chi² p = 3.2 × 10⁻¹¹⁶). This is the most serious internal threat to the paper's biological interpretation. A batch-corrected re-analysis would establish whether the claimed neuron subtypes persist after technical confound removal, or whether some are artifacts of sequencing chemistry.

**WHAT DATA IS NEEDED:** GSE168243 (`geo_data/`). Chemistry labels are available in the data collector report (prep1–3: v2; prep4–6: v3). No new data required.

**SUGGESTED METHOD:** Apply Harmony (R/Python) or ComBat-seq with chemistry (v2/v3) as the batch covariate. Re-cluster after correction using the same parameter space Analyst B tested (k = 12–15, 30 PCs, 2,000 HVGs). Compute ARI between pre- and post-correction cluster assignments. Specifically identify whether the three 100%-pure chemistry clusters persist, merge, or dissolve post-correction. Report proportion of chemistry-pure clusters before and after.

**FEASIBILITY:** HIGH. All data and a tested analysis pipeline are in hand. Harmony is a standard, publicly available package. The analysis requires approximately 1–2 hours of implementation and computation.

**RISK:** Harmony performs integration in PCA space, not in count space; it may over-correct, removing real biological variation that correlates with chemistry version (e.g., if v2 and v3 samples were from donors with different biological characteristics). The two donors whose samples were split across chemistry versions (preparations 3 and 5 from the same donor, both female age 34) provide an internal control: if batch correction works correctly, these two preparations should merge. If they do not, over-correction is occurring.

---

**RECOMMENDATION 4 (HIGH feasibility)**

**QUESTION:** What is the corrected biological interpretation of IL8 and CXCL2 upregulation in NHEK cells, given that the reported log2FC values are actually in log10 scale?

**MOTIVATION:** Paper 1's P1-C8 claim is Tier 3: the fold-change values for IL8 (reported 6.17, actual 1.858 log2FC) and CXCL2 (reported 5.04, actual 1.606 log2FC) are mislabeled. The correct values represent ~3.6-fold and ~3.1-fold upregulation — meaningfully elevated but far less dramatic than the ~73-fold and ~32-fold implied by the paper. The NHEK RNA-seq data is deposited at GSE132174. A corrected analysis could either validate or re-contextualize the claim that NHEK cells in neutrophil-conditioned medium undergo strong pro-nociceptive gene upregulation.

**WHAT DATA IS NEEDED:** GEO accession GSE132174 (NHEK RNA-seq data deposited by the authors; not currently among the locally available files but was referenced by analysts as accessible). The corrected log2FC values (1.858 for IL8; 1.606 for CXCL2) are already established by the audit.

**SUGGESTED METHOD:** Re-run the differential expression analysis from raw counts using a standard pipeline (DESeq2 or edgeR). Generate a complete results table with correct fold-change units throughout. Specifically assess whether the permutation test gene group medians also require correction — Analyst B identified two of four group median values as individual gene values rather than computed medians, and this correction should be done consistently.

**FEASIBILITY:** HIGH for the re-labeling and re-interpretation; MEDIUM for the full re-analysis from raw counts (GSE132174 access required but it is a public GEO dataset).

**RISK:** The corrected analysis may reveal that the NHEK transcriptional response is modest (3–4-fold upregulation), which would weaken the paper's claim that neutrophil-conditioned medium strongly activates a pruritic gene program. This is the scientifically honest outcome — it changes the interpretation of the experiment but does not invalidate the DEG count (Tier 1).

---

**RECOMMENDATION 5 (MEDIUM feasibility)**

**QUESTION:** Do the three CXCL10 ELISA contrasts in Paper 1 survive multiple-comparison correction, and if not, which claims remain valid?

**MOTIVATION:** Paper 1 reports three ELISA comparisons using CXCL10 levels (WT vs. TSLPRKO, with and without neutrophil depletion) without applying FDR or Bonferroni correction. P1-C5 (CXCL10 TSLPRKO, p = 0.031–0.039) and P1-C9 (CXCL10 WT, p ≈ 0.08–0.09 with missing data) are both Tier 2 or Tier 3 precisely because of this issue. The mechanistic interpretation of the paper depends on CXCL10 being specifically regulated by neutrophils via a TSLP-dependent pathway — a claim that requires all three contrasts to hold.

**WHAT DATA IS NEEDED:** The five ELISA data values for the CXCL10 WT MC903 group (one value is missing from the deposited dataset; it must be requested from the corresponding author). The remaining values are in the deposited data files.

**SUGGESTED METHOD:** Apply Bonferroni correction (multiply each p-value by 3) and BH-FDR correction across the three contrasts. Report corrected p-values for each comparison. If the missing WT MC903 value cannot be obtained, perform a sensitivity analysis bounding the corrected p-values under all plausible values for the missing observation (e.g., from the minimum to maximum of observed values in that group).

**FEASIBILITY:** MEDIUM. Requires the missing data point, which means contacting the corresponding author. The authors may provide it; the analysis itself takes minutes once the data are in hand.

**RISK:** After correction, the TSLPRKO contrast may remain significant (p ≈ 0.09–0.12 after Bonferroni) or may not, depending on the exact pre-correction values. If the WT contrast already fails to reproduce at p ≈ 0.08 without correction, the mechanistic CXCL10 story loses its primary comparator. The risk is that the chemokine-level evidence for the neutrophil-CXCL10-CXCR3 mechanism weakens substantially, leaving the behavioral pharmacology (strong, Tier 1) poorly explained at the molecular level.

---

**RECOMMENDATION 6 (MEDIUM feasibility)**

**QUESTION:** Which mouse DRG neuron subtypes express CXCR3, and do they correspond to the human subtype(s) identified in Recommendation 1?

**MOTIVATION:** Paper 1 establishes CXCL10/CXCR3 as the itch-driving axis in mouse but never characterizes which mouse DRG neurons express CXCR3. Published mouse DRG single-cell atlases exist (from laboratories including Usoskin et al., Zeisel et al., and subsequent comprehensive atlases). Querying CXCR3 in an existing public mouse DRG dataset, then comparing CXCR3-positive mouse neurons to the human DRG subtypes in Paper 2, would begin to answer whether Paper 1's pharmacological target acts on a conserved or divergent cell population across species.

**WHAT DATA IS NEEDED:** A publicly available mouse DRG snRNA-seq dataset with per-cell gene expression matrices. Multiple such datasets exist in GEO; the most comprehensive current atlases have been deposited by groups studying pain and itch. No new mouse experiments required. The human comparison data is already in `geo_data/`.

**SUGGESTED METHOD:** Identify CXCR3-positive cells in the mouse dataset by UMI ≥ 1 threshold. Characterize their marker expression profile (Trpm8, Calca, Ntrk1, Ntrk2, Mrgprd, Th, Nefh, etc.) against established mouse DRG subtype definitions (NP1–3, PEP1–2, NF1–5, etc.). Compare the marker profile of mouse CXCR3+ cells to the profile of human CXCR3+ cells (from Recommendation 1). Pearson or Spearman correlation of top marker genes as a cross-species alignment metric.

**FEASIBILITY:** MEDIUM. Data likely publicly available in GEO; requires identifying and downloading the appropriate mouse reference dataset. No new experiments required.

**RISK:** Mouse DRG atlases vary substantially in QC, cell-type annotation granularity, and completeness. CXCR3 is a surface receptor expressed at low levels; detection in snRNA-seq depends heavily on sequencing depth, and many cells may register as CXCR3-negative even if they express the protein. False-negative rate for low-abundance transcripts in snRNA-seq is high. Results may indicate "not detected in mouse snRNA-seq" without ruling out CXCR3 expression at the protein level.

---

**RECOMMENDATION 7 (MEDIUM feasibility)**

**QUESTION:** Do Paper 1's behavioral claims survive re-analysis with the complete AMG487 dataset including all ~20 missing animals?

**MOTIVATION:** Paper 1's Tier 3 claim P1-C10 documents that F(1,67) = 50.64 cannot be recovered from the deposited data because approximately 20 mice are absent. The behavioral effect at Days 8 and 12 is Tier 1 and robust, but the exact statistical test statistic — and therefore the complete sample size — cannot be verified. For regulatory or clinical translation purposes, the precise n and the pooling strategy (multiple cohorts) must be verified.

**WHAT DATA IS NEEDED:** The complete AMG487 behavioral dataset from all pooled cohorts, including the ~20 animals not deposited. This requires a data request to the corresponding author.

**SUGGESTED METHOD:** Upon receiving the full dataset, reproduce the two-way repeated-measures ANOVA (treatment × day) and verify F(1,67) = 50.64. Additionally conduct a cohort effect test (add cohort as a blocking factor) to determine whether the pooling across sessions was valid. Report mixed-effects model results alongside ANOVA for transparency.

**FEASIBILITY:** MEDIUM. The data almost certainly exists in the authors' lab and the request is specific and justified. Author compliance is not guaranteed but the request is reasonable.

**RISK:** Authors may be unable or unwilling to locate raw behavioral data from a 2019 paper. If the full dataset is obtained and the F-statistic still cannot be matched, this would suggest more serious data integrity concerns. The risk to the overall behavioral conclusion is low (the effect is real, Tier 1), but the documentation failure could have implications for the paper's official record.

---

**RECOMMENDATION 8 (LOW feasibility)**

**QUESTION:** Do anti-Gr1 neutrophil depletion (79% scratching reduction) and TSLPRKO (44% scratching reduction) results replicate in independent experiments with individual-level data fully documented?

**MOTIVATION:** These are the two headline claims of the paper (P1-C14 and P1-C15, both Tier 4) — they are what makes the title claim "neutrophils promote CXCR3-dependent itch" testable at a mechanistic level. Without individual animal data, neither claim can be audited, and neither should be cited as established. The behavioral pharmacology (Tier 1) supports the conclusion that CXCR3 mediates itch, but the neutrophil-to-CXCL10 upstream link depends entirely on these unverifiable results.

**WHAT DATA IS NEEDED:** New or replicated animal experiments using the MC903 model with anti-Gr1 antibody depletion (or CXCL10-knockout mice, or a complementary depletion strategy) and TSLPRKO scratch behavior, with individual animal scratch count time-series recorded and deposited.

**SUGGESTED METHOD:** Pre-registered replication design with n ≥ 10 per group (power analysis based on the effect size reported in the paper, if the exact size could be estimated from published figures). Use video-based automated scratch counting (standardized across cohorts). Report individual animal data in supplementary materials. Separately, if the original individual-level data can be obtained from the authors, retrospective analysis is sufficient.

**FEASIBILITY:** LOW. New animal experiments are required if original data cannot be retrieved. This involves animal ethics approval, significant time, and substantial cost.

**RISK:** If replication fails — if the 79% and 44% reductions are not reproduced — the mechanistic case for neutrophil-driven CXCL10 signaling weakens substantially, though the CXCR3 pharmacology (Tier 1) would remain unaffected. The risk of not pursuing this analysis is that influential claims about specific biology (neutrophil necessity) continue to be cited in the literature without verifiable support.

---

> **Blocked recommendation:** A follow-up analysis of Paper 2's cross-species KL divergence claims (P2-C9 and P2-C10) is blocked because both are Tier 4: the mouse reference dataset, the KL implementation, and the co-clustering method are all absent. No meaningful follow-up can be built on these claims until the methods and data are disclosed.

---

## PART 4 — RECOMMENDED READING MAP

### (a) Papers that could immediately unlock HIGH-feasibility analyses

**1. Published human DRG snRNA-seq papers with analysis code deposited (to unlock Recommendation 3 and contextualize Recommendation 1).**
What to look for: Studies that processed GSE168243 or a similar human DRG atlas using Seurat with fully documented parameters, preferably with Harmony batch correction for multi-chemistry datasets. An ideal paper would include the resolution parameter, PC count, HVG selection method, and integration algorithm used — providing the missing parameter chain for Paper 2.
What it contributes: Directly provides the undisclosed Seurat parameters needed to reproduce Paper 2's H1–H15 classification, unlocking Recommendation 3 (batch correction) and enabling proper cluster label mapping for Recommendations 1 and 2 (CXCR3 and IL31RA/OSMR queries).
Most directly supports: **Recommendation 3.**

**2. Papers establishing CXCR3 protein expression in DRG neurons by immunohistochemistry or single-cell proteomics (to validate the low-abundance snRNA-seq signal in Recommendation 1).**
What to look for: IHC or flow cytometry studies quantifying CXCR3 surface expression in DRG tissue, either mouse or human, particularly from itch or pain research groups. Studies using CXCR3 reporter mice (CXCR3-EGFP) would be especially informative.
What it contributes: Provides protein-level validation of the transcript-level query from Recommendation 1, mitigating the low-detection risk identified as the primary methodological concern.
Most directly supports: **Recommendation 1.**

---

### (b) Papers that would resolve the most important contradiction between the papers

**3. Comparative transcriptomic studies directly aligning human and mouse DRG neuron types (to resolve Part 2B Divergence 1).**
What to look for: Papers explicitly mapping human DRG neuron subtypes (from Paper 2's atlas or similar datasets) to mouse subtype classifications using published cross-species alignment methods (SAMap, LIGER, or direct KL divergence with disclosed parameters). Particularly relevant are papers that conclude whether a human equivalent to mouse NP1 (the primary CXCR3-candidate in mouse) exists.
What it contributes: Resolves whether the translatability concern is real or overstated — whether Paper 1's CXCR3 finding has a direct human correlate.
Most directly supports: **Recommendation 6, and Part 2B Divergence 1.**

**4. Papers reporting the mouse DRG cell type identity of neurons expressing CXCR3, CXCL10, or responding to CXCL10 by calcium imaging or electrophysiology.**
What to look for: Studies combining single-cell transcriptomics with functional assays in mouse DRG, specifically identifying which neuronal subtype(s) respond to chemokine receptor activation. Papers from labs working at the intersection of neuroimmunology and DRG physiology in the context of skin inflammation.
What it contributes: Provides the cell-type identity missing from Paper 1's pharmacological study, grounding its behavioral results in a specific circuit element.
Most directly supports: **Recommendations 1 and 6.**

---

### (c) Foundational methodological papers that would improve the analytical approach

**5. Statistical guidelines for pseudoreplication in IHC and behavioral neuroscience (to address Paper 1 P1-C13).**
What to look for: Methodological papers or statistical guidelines specifically addressing image-level vs. animal-level inference in histology, nested data structures in behavioral experiments, and appropriate use of mixed-effects models when pooling across cohorts. The work of Lazic et al. on experimental design in neuroscience is the canonical starting point.
What it contributes: Provides the methodological framework for re-analyzing Paper 1's IHC data if additional animals are available, and for designing the replication experiment in Recommendation 8.
Most directly supports: **Recommendation 8.**

**6. Benchmarking papers for single-cell batch correction methods (Harmony, ComBat-seq, scVI, RPCA) specifically tested on multi-chemistry snRNA-seq data.**
What to look for: Papers comparing batch correction performance on snRNA-seq datasets with mixed 10x Genomics v2 and v3 chemistry, ideally from neural tissue. Papers that quantify over-correction risk and provide guidance on choosing the right method for datasets where chemistry and biology partially co-vary.
What it contributes: Provides the methodological justification and parameter choices for Recommendation 3's batch correction analysis.
Most directly supports: **Recommendation 3.**

**7. Papers establishing QC thresholds and analysis standards for low-n human DRG snRNA-seq (to contextualize Paper 2's pseudo-replication problem).**
What to look for: Studies grappling with DRG snRNA-seq from small human donor cohorts (n = 3–6 donors), specifically addressing how to handle repeated measures from the same donor and how to apply mixed-effects models or donor-level random effects in scRNA-seq clustering.
What it contributes: Provides the analytical framework needed to address the single most serious statistical concern in Paper 2 — that 46.3% of the cell population comes from one donor — and for designing any follow-up human DRG study.
Most directly supports: **Recommendations 3 and 2.**

---

## PART 5 — HONEST LIMITATIONS STATEMENT

**What this audit process can tell us:**
The reproducibility audit protocol used here — independent data collection, methodology extraction, and three isolated analyst perspectives — is substantially more rigorous than informal replication attempts. When analysts converge on the same numbers (the chloroquine t-statistic, the cell counts, the TMEM100 percentage), those findings deserve genuine confidence. The protocol is well-calibrated for detecting numerical irreproducibility, missing data, and explicit statistical assumption violations.

**What this audit process cannot tell us:**
The audit cannot establish biological truth. A claim can be numerically unverifiable (Tier 4) and still be scientifically correct — both headline behavioral claims in Paper 1 (79% and 44% reductions) may be accurate, but we have no way to know. Conversely, a Tier 1 claim is reproducible within the assumptions of the original design, but those assumptions (mouse model validity, the MC903 protocol as a faithful AD model, human DRG from 5 donors representing the population) are not themselves audited. The audit evaluates the internal consistency and documentability of the analysis; it says nothing about external validity.

**Where the synthesis in Parts 2 and 3 is speculative:**
The convergence on the OSMR/IL-31 axis (Part 2A, Convergence 2) is a hypothesis generated by juxtaposing two findings from different experimental systems — one a Tier 1 DEG count, one a Tier 2 expression percentage. No experiment in either paper directly tested whether OSMR-expressing human DRG neurons respond to the IL-31 pathway activated in Paper 1's NHEK model. This is a literature-level inference, not an audited finding. Similarly, the divergence discussion about mouse-to-human translatability (Part 2B, Divergence 1) rests partly on Paper 2's cross-species claims, which are Tier 4. If those claims are wrong — if the mouse-human subtype correspondence is in fact closer than Paper 2 implies — the divergence largely disappears.

**What readers should be most cautious about when acting on these recommendations:**

The recommendations in Part 3 are only as strong as the underlying claims they are designed to test or extend. Recommendations 1 and 2 (querying CXCR3 and IL31RA/OSMR in Paper 2's data) are HIGH feasibility, but they rest on Paper 2's cluster assignments, which carry a mean ARI of 0.547 and an unaddressed batch confound. A CXCR3 subtype-enrichment result from this data should be treated as a preliminary observation requiring validation by orthogonal methods (protein staining, functional assay), not as an established cell-type identification.

Recommendation 3 (batch correction) is technically rigorous but carries the specific risk noted: over-correction may remove genuine biology. If a batch-corrected analysis dissolves claimed subtypes, this could reflect either the removal of a technical artifact (good) or the removal of real donor-level biological variation (problematic). The result of Recommendation 3 must be interpreted cautiously and cannot by itself establish or refute the existence of any particular DRG subtype.

Most critically: both papers underlying this synthesis are themselves PARTIALLY REPRODUCIBLE. A synthesis built on two partially reproducible papers inherits all of their uncertainties. The integrative picture offered in Part 2 — a CXCL10/CXCR3 signal from inflamed keratinocytes reaching a specific DRG subtype — is scientifically plausible and worth investigating, but it is a narrative constructed from Tier 1 and Tier 2 evidence in two different species, without any single experiment that tests the complete mechanism. The reader should treat it as a well-motivated hypothesis, not a synthesis of established facts.
