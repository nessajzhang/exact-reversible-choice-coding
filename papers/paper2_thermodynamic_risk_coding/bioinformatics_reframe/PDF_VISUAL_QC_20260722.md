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

- `build/main.pdf`: `58c0552ef0456ab5a4983be3a7504ee321a077baadebd6722756d91865b4078a`
- `build/supplementary_codec_evidence.pdf`: `c9da4ff78af0a585f855a509f1744f2b32bdc6be214698b420d5d6ea6909cb65`
- `oup_preflight/build/main_oup_preflight.pdf`: `e1840d67cb0582d02b006138e6767ae6273ebfa3fd6b379a88e6037f90f4d229`

## Findings

- No page, figure, formula, table, caption or reference block is visibly clipped.
- Figure 1 and Figure 2 remain within the OUP text area and their panel labels, axes and legends are legible at the compiled page size.
- The large `Overfull \\hbox` messages emitted while the OUP output routine places full-width figures reflect the double-column float mechanism; visual inspection found no content outside the page or figure bounds.
- The main review PDF and Supplementary Information show no overlapping text, missing glyphs, unresolved-reference markers or abnormal blank pages.
- The OUP preflight is exactly seven pages without changing the class font size, margins or line spacing.

## Boundary

This is a rendered-layout audit, not a scientific-result audit or a guarantee that the publisher production system will make identical line and float breaks. The source manuscript, figures and bibliography remain the controlling submission artifacts.
