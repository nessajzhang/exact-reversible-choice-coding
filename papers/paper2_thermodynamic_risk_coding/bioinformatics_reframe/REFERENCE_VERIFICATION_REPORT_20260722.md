# Exact reversible choice coding — prior-art verification report (2026-07-22)

## Scope and method

The six references added during the final pre-submission audit were checked by exact DOI lookup in the configured academic-search service, Crossref REST metadata, DOI content negotiation, and, where an article number was absent from Crossref, the publisher landing page. Titles, author order, journal, year, volume/issue, article or page number, and DOI were compared field by field. No citation was inferred from a search-result title alone.

## Verified records

| BibTeX key | Verified record | DOI | Cross-check outcome |
|---|---|---|---|
| `dou2024explorer` | Dou C, Yang Y, Zhu F, Li B, Duan Y. *Explorer: efficient DNA coding by De Bruijn graph toward arbitrary local and global biochemical constraints*. **Briefings in Bioinformatics** 25(5), bbae363 (2024). | 10.1093/bib/bbae363 | PASS; DOI, author order, volume/issue and e-locator agree across Crossref and OUP metadata. |
| `press2020hedges` | Press WH, Hawkins JA, Jones SK Jr, Schaub JM, Finkelstein IJ. *HEDGES error-correcting code for DNA storage corrects indels and allows sequence constraints*. **PNAS** 117(31), 18489–18496 (2020). | 10.1073/pnas.2004821117 | PASS; Crossref and DOI citation metadata agree. |
| `zhao2024compositehedges` | Zhao X, Li J, Fan Q, Dai J, Long Y, Liu R, Zhai J, Pan Q, Li Y. *Composite Hedges Nanopores codec system for rapid and portable DNA data readout with high INDEL-Correction*. **Nature Communications** 15, 9395 (2024). | 10.1038/s41467-024-53455-3 | PASS; article number 9395 and publication date verified on the publisher page. |
| `volkel2023framed` | Volkel KD, Lin KN, Hook PW, Timp W, Keung AJ, Tuck JM. *FrameD: framework for DNA-based data storage design, verification, and validation*. **Bioinformatics** 39(10), btad572 (2023). | 10.1093/bioinformatics/btad572 | PASS; OUP page confirms issue, e-locator, date and author order. |
| `gimpel2026codecbenchmark` | Gimpel AL, Remschak A, Stark WJ, Heckel R, Grass RN. *Comparison of state-of-the-art error-correction coding for sequence-based DNA data storage*. **Nature Communications** 17, 3963 (2026). | 10.1038/s41467-026-70548-3 | PASS; article number 3963 and 14 March 2026 publication date verified on the publisher page. |
| `zhang2026gungnir` | Zhang J, Chen L, Sun J, Li S, Zhou Y, Wu Z, Li C, Zheng Z, Luo R. *Gungnir codec enabling high error-tolerance and low-redundancy DNA storage through substantial computing power*. **Nature Communications** 17, 4828 (2026). | 10.1038/s41467-026-71485-x | PASS; Crossref metadata and publisher page agree on title, authors, volume and article number. |

## Claim-use audit

- Explorer is cited as constrained De Bruijn-graph coding, not as an assay-guided fixed-fiber method.
- HEDGES, Composite Hedges Nanopores and Gungnir are cited as channel-error/recovery systems, not as direct PCR-ranking comparators.
- FrameD is cited as system-level simulation and verification infrastructure.
- The 2026 Nature Communications study is cited for standardized codec/error/dropout benchmarking and the separation between a shaping layer and a complete recovery stack.
- None of these scope citations is represented as a numerical comparator run on the present PCR fibers.

The executable numerical comparator remains the released Gimpel 1D-CNN output and published hard rule evaluated on identical fibers; the new records refine prior-art positioning only.
