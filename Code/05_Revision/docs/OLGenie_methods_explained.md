# The OLGenie site-classification approach to dN/dS in overlapping frames

**What this document is for.** You raised the objection that a dN/dS computed on the
alternate reading frame of a DCR cannot be independent of the canonical frame, because
purifying selection on the canonical protein constrains codon positions 1 and 2, and those
same nucleotides are the alternate frame's synonymous and nonsynonymous sites. You are
right, and this document explains exactly why, what the standard fix is, how I implemented
it, and what the resulting numbers mean.

---

## 1. Why a naive alternate-frame dN/dS is uninterpretable

Consider a DCR where a second transcript reads the same nucleotides shifted by +1
(the case OLGenie calls `ss13`). Line up the two frames:

```
nucleotide position   1   2   3   4   5   6   7   8   9
canonical frame       |---codon A---| |---codon B---| |--
                       p1  p2  p3      p1  p2  p3
alternate frame (+1)      |---codon a---| |---codon b---|
                           p1  p2  p3      p1  p2  p3
```

Reading the mapping off the diagram, for a +1 shift:

| canonical position | alternate position |
|---|---|
| p1 (highly constrained) | p3 (mostly degenerate) |
| p2 (most constrained)   | p1 |
| p3 (wobble, degenerate) | p2 |

Now ask what happens if the canonical protein is under strong purifying selection.
Substitutions are **suppressed** at canonical p1 and p2, and **tolerated** at canonical p3.
But canonical p3 is alternate p2 — the alternate frame's *most* constrained position. So the
only positions that freely accumulate substitutions are precisely the ones that look
nonsynonymous in the alternate frame.

The consequence: the alternate frame's apparent dN is inflated (its substitutions are
concentrated at alternate-p2), and its apparent dS is deflated (alternate-p3 sits on
canonical-p1, where change is suppressed). **Naive omega_alt is biased upward.** A DCR under
genuine alternate-frame purifying selection could easily read omega_alt near or above 1, and
you would wrongly conclude there is no constraint — or, in the other direction, mistake the
canonical frame's own constraint for alternate-frame signal.

For a +2 shift (`ss12`) the position mapping is different again, so the bias has a
different magnitude. This is why the two frame relationships must never be pooled.

This is not a subtlety I invented: it is the motivating problem of the overlapping-gene
dN/dS literature. Standard methods do not account for the fact that a change which is
synonymous in one frame is usually nonsynonymous in the other, which is why they fail to
detect purifying selection, or erroneously report positive selection, when applied to
overlapping reading frames.

---

## 2. The fix: cross-classify every site by its effect in BOTH frames

Instead of asking "is this substitution synonymous?" (which has no single answer when two
frames overlap), ask **"what is it in each frame?"** Every site gets a two-letter label:

- first letter = effect in the **reference** (canonical) frame
- second letter = effect in the **alternate** frame

giving four classes:

| class | canonical frame | alternate frame |
|---|---|---|
| **NN** | nonsynonymous | nonsynonymous |
| **NS** | nonsynonymous | synonymous |
| **SN** | **synonymous** | nonsynonymous |
| **SS** | **synonymous** | synonymous |

From these we form four normalised divergences dNN, dNS, dSN, dSS — each is
(observed differences of that class) / (number of sites of that class).

Selection on the alternate frame is then measured by a ratio that **holds the canonical
letter constant**:

```
omega_alt  =  dSN / dSS        <-- both canonical-synonymous
omega_alt  =  dNN / dNS        <-- both canonical-nonsynonymous
```

and, as a control, selection on the canonical frame by holding the *alternate* letter fixed:

```
omega_ref  =  dNN / dSN        <-- both alternate-nonsynonymous
```

### Why this answers your objection exactly

Look at `dSN/dSS`. Every single change counted in the numerator **and** every single change
counted in the denominator is **synonymous in the canonical frame**. Whatever purifying
selection is acting on the canonical protein applies identically to both, so it **cancels in
the ratio by construction**. The ratio is then free to reflect only what differs between the
numerator and denominator — namely, whether the change alters the alternate protein.

`dNN/dNS` does the same thing at canonical-nonsynonymous sites. It has far more sites
available (canonical-nonsynonymous sites vastly outnumber canonical-synonymous ones), so it
is better powered in short regions; `dSN/dSS` is the conceptually cleanest answer to your
objection. I report both. They agree, which is reassuring — if they disagreed, that would
itself be diagnostic of a residual frame-specific artefact.

---

## 3. How I implemented it

**Sequence pairs.** Human DCR intervals and their mouse orthologs, lifted with UCSC
`liftOver` (hg38 -> mm39). Only **indel-free** pairs were used (human and mouse interval the
same length), so codons stay in register without needing a gapped alignment. This is a real
filter: 12.1% of lifted DCRs and 15.0% of lifted control intervals were dropped for indels —
comparably, so it does not bias the contrast.

**Reading frames.** *Both* DCR frames are taken from the pipeline's own protein-level
evidence, not from GENCODE. The authoritative source is the filtered UniProtKB isoform
models (`FINAL_Filtered_UniprotKB_isoforms_mapping_toDCR_*.extBed.bed`, 36,733 isoform
models). These BED12 records were exploded into CDS blocks in 5'->3' order and a phase
assigned to each block as the cumulative preceding CDS length modulo 3; the phase was then
converted to an offset within the DCR interval accounting for strand. Intersecting the DCR
intervals against these models recovers **>=2 distinct reading frames for all 1,741 intervals**
(1,718 with exactly 2, 23 with 3), which is what a protein-evidence-based DCR call
guarantees by construction. Frame relationship splits 661 `ss13` / 567 `ss12` among scored
intervals.

An earlier version of this analysis derived the DCR frames from GENCODE v48 CDS annotation
instead, and consequently classified 266 intervals as having "only one coding frame" because
their second isoform overlaps the interval without a GENCODE CDS annotation there. That was
an error of authority: the DCR calls come from protein evidence, so the interval is known to
be coding in both frames regardless of whether GENCODE annotates it. Those 266 have been
restored with their correct protein-derived frames (see section 5). GENCODE phase is still
used for the *control* set, where it is the appropriate source: those intervals are ordinary
annotated CDS with a detected latent ORF as the alternate frame.

**Site and difference counting.** For each nucleotide position, both the canonical codon and
the alternate codon containing it are extracted. Each of the three possible alternative bases
is evaluated in both codons simultaneously and the position contributes 1/3 of a site to
whichever of the four classes applies. Observed human-mouse differences are classified the
same way. Counts are pooled across intervals, converted to proportions, and
Jukes-Cantor-corrected before ratios are taken.

**A bug this design caught.** My first implementation converted the annotated codon position
to a frame offset as `offset = cpos`, when the correct conversion is `offset = (3 - cpos) % 3`.
The error was invisible in the DCR set (both of its frames are genuinely coding, so a
scrambled offset still lands on a real reading frame) but catastrophic in the control set. It
surfaced through the built-in positive control: omega_ref for the controls read 0.98
(indistinguishable from neutral) when it must be strongly purifying, because these are
annotated CDS. After the fix, omega_ref = 0.110. **This is why the positive control is not
optional** — it is the only thing that detects a scrambled frame.

---

## 4. Why an absolute omega_alt still would not be enough

Even after canonical constraint cancels, the sites carry other non-neutral pressures:
codon usage bias, exonic splicing enhancers, and the GC enrichment we documented separately
for DCRs. Any of these could push omega_alt below 1 without any alternate-frame protein
existing.

So the claim is **never** "omega_alt < 1, therefore selection." The claim is the **contrast
against a matched empirical null**: intervals inside non-DCR exons that carry a latent
dual-coding ORF but are not annotated as dual-coding. These are matched on length and
stratified by frame relationship, and they experience the same codon-usage, splicing, and
composition pressures. The null answers the question "what does omega_alt look like when
there is dual-coding *potential* but no dual-coding *function*?"

The answer: **omega_alt = 0.984, 95% CI [0.951, 1.017]** — the null sits exactly at
neutrality. That is the single most important number in this analysis, because it validates
that the design has successfully removed the confounds. The contrast is therefore meaningful.

---

## 5. Results

| | real DCRs (n=1,228) | latent null (n=21,340) |
|---|---|---|
| omega_alt = dNN/dNS | **0.244** [0.223, 0.264] | 0.798 [0.787, 0.809] |
| omega_alt = dSN/dSS | **0.446** [0.376, 0.529] | 0.984 [0.951, 1.017] |
| omega_ref = dNN/dSN (control) | 0.260 | 0.110 |

Bootstrap over intervals, 2000 replicates. Difference in omega_alt between sets: 0.554
[0.530, 0.576] under dNN/dNS and 0.536 [0.449, 0.613] under dSN/dSS, p < 0.001 for both.

**Controls all pass:**
- canonical-frame control strongly purifying in both sets (0.260, 0.110) — frames are correct
- holds separately in `ss13` (D 0.443 vs N 0.827) and `ss12` (D 0.451 vs N 1.196)
- survives length matching (D 0.439 vs N 1.034) — the control that erased the earlier
  potential-retention signal
- DCR and null are separated in **every** length bin from 30 nt to 300 nt, with no bin where
  the DCR estimate approaches the null

**Ascertainment.** All 1,741 DCR intervals carry >=2 protein-derived frames, so no interval is
excluded for biological reasons. The 513 not scored were lost to purely technical filters:
245 shorter than 30 nt, 201 with indels in the human-mouse lift, 49 that failed to lift. The
indel filter is symmetric between arms (12.1% of lifted DCRs, 15.0% of lifted controls), so it
cannot manufacture the contrast.

The 266 intervals that the earlier GENCODE-based version misclassified were rescored with
their correct protein-derived frames and show **purifying selection, not neutrality**:
dSN/dSS = 0.506 [0.359, 0.726], with a canonical control of 0.277. Under the earlier
mistaken frame assignment they had read 1.138. They are statistically indistinguishable
from the other 956 (0.429), which is the expected result if they are ordinary DCRs that
simply lack a GENCODE CDS annotation on the second isoform.

## 6. What this licenses you to claim, and what it does not

**Supported:** the alternate reading frames of annotated DCRs are under purifying selection,
significantly more so than sequence-matched intervals with equivalent dual-coding potential
but no annotated dual-coding function. Since dual-coding potential is near-universal in
coding exons (X = 0.994), this selection signal — not the mere existence of potential — is
what distinguishes real DCRs.

**Not supported:** (i) any per-DCR claim. A single human-mouse pairwise comparison over a
median 72 nt yields too few substitutions to rank individual DCRs; this is a pooled,
genome-wide result only. (ii) any claim about which specific alternate proteins are made —
selection on a reading frame is evidence of function, not proof of translation.

**Natural extensions if a reviewer pushes:** multi-species alignments (UCSC 100-way) would
give per-DCR resolution and could identify a strongly-constrained subset; human polymorphism
(gnomAD) would test present-day constraint with far more variants. Both were scoped and are
straightforward; neither is needed for the global claim.

---

## 7. Reference

Nelson CW, Ardern Z, Wei X. *OLGenie: Estimating Natural Selection to Predict Functional
Overlapping Genes.* Molecular Biology and Evolution 37(8):2440-2449, 2020.
The dNN/dNS/dSN/dSS cross-classification and the ratio estimators used here are theirs; the
implementation in this analysis is my own, written to work directly on the lifted interval
pairs, with the reference-frame estimator retained as a positive control.
