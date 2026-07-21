# Reversible choice-codec runtime benchmark

Status: `COMPLETE_MACHINE_SPECIFIC_RUNTIME_AUDIT`

Cold exact completion-table construction: 0.432223 s; cache entries after all checks: 41,101.

Median timings:

| Operation | r | Candidates | Items | Median ms/item | Peak process RSS (MiB) |
|---|---:|---:|---:|---:|---:|
| P5_choice_encode | 0 | 1 | 128 | 8.576748 | 141.801 |
| P5_choice_encode | 2 | 4 | 128 | 42.591569 | 142.012 |
| P5_choice_encode | 4 | 16 | 128 | 154.996404 | 142.066 |
| base_rank | 0 | 1 | 4096 | 0.840850 | 138.582 |
| base_unrank | 0 | 1 | 4096 | 0.927683 | 138.582 |
| choice_decode | 0 | 1 | 128 | 0.914650 | 141.801 |
| choice_decode | 2 | 4 | 128 | 0.974828 | 142.012 |
| choice_decode | 4 | 16 | 128 | 1.021818 | 142.066 |

All timed encode/decode and base rank/unrank checks passed. Timings are machine-specific, single-process implementation measurements; they are not biological replicates, a physical-density result or an end-to-end system comparison.
