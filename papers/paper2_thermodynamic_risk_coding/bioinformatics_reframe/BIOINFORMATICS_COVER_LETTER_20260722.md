# Exact reversible choice coding with retrospective assay-guided selection in PCR-tested oligonucleotide libraries

## Bioinformatics cover letter — author-completion draft — 22 July 2026

Dear Editors,

We submit the Original Paper entitled “Exact reversible choice coding with retrospective assay-guided selection in PCR-tested oligonucleotide libraries” for consideration in *Bioinformatics*.

The manuscript addresses a computational molecular-design problem: how to connect an assay-specific sequence scorer to an exact constrained oligonucleotide codec without making payload recovery depend on that scorer. The interface assigns each payload a disjoint, fixed-size choice fiber, incurs an exact (r)-bit rate cost and recovers the payload from every valid fiber member without the scorer. Optional canonical-output verification is a separate task. A scorer may therefore be lost, upgraded or replaced without compromising recovery of previously encoded payloads.

The contribution is not rank/unrank, affine permutation or representative selection in isolation. It is the combined contract of full-domain exact reversibility, fixed candidate budget, scorer-independent inversion, scorer-replacement invariance, separate canonical verification and identical-fiber benchmarking. We evaluate this interface retrospectively using real PCR measurements, including an external-laboratory pool reported in the same source publication.

The evidence boundaries are explicit. This is not prospective wet-lab validation of codec-emitted sequences: generated candidates were computationally verified and scored but were not synthesized or assayed. We do not claim an assay-independent sequence-quality function, universal improvement, general superiority over released 1D-CNN outputs or coverage of all DNA-storage failure modes. The released-CNN comparisons are endpoint-dependent, 17.2--30.9% of matched-condition fibers have negative measured gain and the secondary exploratory Q5 analysis reverses direction. Figure 2 presents the matched positive and altered-protocol negative intervals on the same zero-effect axis.

We believe the manuscript is suited to *Bioinformatics* because it links an exact algorithmic interface for constrained DNA sequence design to measured biological outcomes, identical-fiber comparators, conditional two-stage uncertainty analysis, explicit transfer failure and a reproducible implementation. The OUP preflight is seven pages.

The manuscript is original, is not under consideration elsewhere and has been approved by all authors. **[AUTHOR CONFIRMATION REQUIRED]**

Code, tests, frozen derived data and reproduction commands are public at https://github.com/nessajzhang/exact-reversible-choice-coding under the BSD-3-Clause licence. **[INSERT VERIFIED RELEASE TAG, COMMIT AND ARCHIVE DOI AFTER CREATOR APPROVAL]**

During author-directed revision, OpenAI Codex (GPT-5 model family) assisted with code review and modification, test generation, literature-metadata checking, internal-consistency auditing, presentation review and language suggestions. All reported numerical results were produced by the versioned analysis scripts; the tool generated no experimental data. The authors retain responsibility for the mathematical arguments, code, references, analyses and final wording. **[ALL-AUTHOR VERIFICATION REQUIRED]**

Sincerely,

**[Corresponding author]**  
**[Affiliation]**  
**[Email]**
