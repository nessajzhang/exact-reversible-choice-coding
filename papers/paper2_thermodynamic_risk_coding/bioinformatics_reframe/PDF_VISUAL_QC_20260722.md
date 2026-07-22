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

- `build/main.pdf`: `74b8c746497c96ae37cddae52dedf08a0867dd5393b4c0817d3e0de6961e4c5c`
- `build/supplementary_codec_evidence.pdf`: `7f607d790c5d7da16dbbaf5a615988e7a2816a6e3cc909115857c1a484854873`
- `oup_preflight/build/main_oup_preflight.pdf`: `21df701b282672871a6d8402afd2e9023613dec941ec15acd139459d3383e0b4`

## Findings

- No page, figure, formula, table, caption or reference block is visibly clipped.
- Figure 1 and Figure 2 remain within the OUP text area and their panel labels, axes and legends are legible at the compiled page size.
- Figure 2d places cross-pool, external-KAPA and secondary exploratory Q5 two-stage intervals on one zero-effect axis. Positive matched-condition points and negative Q5 points remain distinguishable by marker shape and colour; no interval or legend is clipped.
- The large `Overfull \\hbox` messages emitted while the OUP output routine places full-width figures reflect the double-column float mechanism; visual inspection found no content outside the page or figure bounds.
- The main review PDF and Supplementary Information show no overlapping text, missing glyphs, unresolved-reference markers or abnormal blank pages.
- The OUP preflight is exactly seven pages without changing the class font size, margins or line spacing. Six recent graph/channel/recovery references remain fully cited in the Supplementary prior-art matrix rather than expanding the main-paper bibliography beyond the page limit.

## Boundary

This is a rendered-layout audit, not a scientific-result audit or a guarantee that the publisher production system will make identical line and float breaks. The source manuscript, figures and bibliography remain the controlling submission artifacts.
