# Paper 2 PDF visual-quality audit

Audit date: 2026-07-22

Status: `PASS_RENDERED_PAGE_QC; OUP_PREFLIGHT_7_PAGES`

## Files inspected

| Artifact | Pages | Page format | Result |
|---|---:|---|---|
| `build/main.pdf` | 9 | A4 review layout | PASS |
| `build/supplementary_codec_evidence.pdf` | 17 | A4 review layout | PASS |
| `oup_preflight/build/main_oup_preflight.pdf` | 7 | OUP Bioinformatics preflight | PASS |

All pages were rendered with Poppler and inspected as full-document contact sheets. The two figure-bearing OUP pages and the reference/availability page were additionally inspected at page resolution.

Inspected final SHA-256 values:

- `build/main.pdf`: `095ba48e4b0d6d86560af4f679dd43743fd7c98a2f18b3d9f04292cd28025dbb`
- `build/supplementary_codec_evidence.pdf`: `0246c7586b568e85304d713d0c8d47e9891027e3edc83de2807a2bd75f3b50be`
- `oup_preflight/build/main_oup_preflight.pdf`: `6002c6bf7d72cc77a2ef8faef0535ccfccfff9ff7f5423c7feccecc02e15fa98`

## Findings

- No page, figure, formula, table, caption or reference block is visibly clipped.
- Figure 1 and Figure 2 remain within the OUP text area and their panel labels, axes and legends are legible at the compiled page size.
- Figure 2d places cross-pool, external-KAPA and secondary exploratory Q5 two-stage intervals on one zero-effect axis. Positive matched-condition points and negative Q5 points remain distinguishable by marker shape and colour; no interval or legend is clipped.
- The large `Overfull \\hbox` messages emitted while the OUP output routine places full-width figures reflect the double-column float mechanism; visual inspection found no content outside the page or figure bounds.
- The main review PDF and Supplementary Information show no overlapping text, missing glyphs, unresolved-reference markers or abnormal blank pages.
- The OUP preflight is exactly seven pages without changing the class font size, margins or line spacing. Six recent graph/channel/recovery references remain fully cited in the Supplementary prior-art matrix rather than expanding the main-paper bibliography beyond the page limit.
- The revised 149-word abstract sentence about inversion and payload recovery is legible in both the review and OUP layouts and introduces no line collision or clipping.
- The revised persistent-archive wording was re-rendered and inspected on main-paper page 7, Supplementary page 17 and OUP page 7; it remains within the text area and does not imply that a DOI has already been created.

## Boundary

This is a rendered-layout audit, not a scientific-result audit or a guarantee that the publisher production system will make identical line and float breaks. The source manuscript, figures and bibliography remain the controlling submission artifacts.
