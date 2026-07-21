# External KAPA mapping sensitivity

Status: `POST_HOC_EXPLORATORY_ALGORITHMIC_SENSITIVITY`

The analysis-plan-locked mapping remains primary. Thirty-two additional namespaces depend only on a fixed index and sequence SHA-256. Ranges below are deterministic min--median--max summaries, not confidence intervals or biological repeats.

- GCall, r=2: FullContext gain 0.00244063--0.00265519--0.00298752 (32/32 positive); FullContext-minus-AssayContext 0.00047662--0.00080106--0.00134060 (32/32 positive); FullContext-minus-CNN 0.00076106--0.00113616--0.00156641 (32/32 positive).
- GCall, r=4: FullContext gain 0.00285422--0.00378461--0.00440120 (32/32 positive); FullContext-minus-AssayContext 0.00019314--0.00132695--0.00212793 (32/32 positive); FullContext-minus-CNN 0.00097495--0.00205135--0.00313682 (32/32 positive).
- GCfix, r=2: FullContext gain 0.00233205--0.00253388--0.00290979 (32/32 positive); FullContext-minus-AssayContext 0.00045301--0.00080821--0.00138723 (32/32 positive); FullContext-minus-CNN 0.00067983--0.00104349--0.00151975 (32/32 positive).
- GCfix, r=4: FullContext gain 0.00295867--0.00353341--0.00412392 (32/32 positive); FullContext-minus-AssayContext 0.00074027--0.00142503--0.00216994 (32/32 positive); FullContext-minus-CNN 0.00097015--0.00192148--0.00370313 (32/32 positive).

Evidence boundary: post hoc deterministic mapping sensitivity on public external-KAPA measurements; mappings are not biological repeats, do not define a confidence interval and do not replace the analysis-plan-locked primary mapping.
