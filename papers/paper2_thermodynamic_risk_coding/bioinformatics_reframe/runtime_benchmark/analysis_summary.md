# Reversible choice-codec runtime benchmark

Status: `COMPLETE_MACHINE_SPECIFIC_RUNTIME_AUDIT`

Cold exact completion-table construction: 0.100001 s; cache entries after all checks: 41,101.

Median timings:

| Operation | r | Candidates | Items | Median ms/item | Peak process RSS (MiB) |
|---|---:|---:|---:|---:|---:|
| P5_choice_encode | 0 | 1 | 128 | 1.705808 | 194.703 |
| P5_choice_encode | 2 | 4 | 128 | 6.842102 | 194.773 |
| P5_choice_encode | 4 | 16 | 128 | 28.914271 | 195.613 |
| base_rank | 0 | 1 | 4096 | 0.188970 | 194.516 |
| base_unrank | 0 | 1 | 4096 | 0.191709 | 194.516 |
| choice_decode | 0 | 1 | 128 | 0.178113 | 194.703 |
| choice_decode | 2 | 4 | 128 | 0.180324 | 194.773 |
| choice_decode | 4 | 16 | 128 | 0.173910 | 195.613 |

All timed encode/decode and base rank/unrank checks passed. Timings are machine-specific, single-process implementation measurements; they are not biological replicates, a physical-density result or an end-to-end system comparison.
