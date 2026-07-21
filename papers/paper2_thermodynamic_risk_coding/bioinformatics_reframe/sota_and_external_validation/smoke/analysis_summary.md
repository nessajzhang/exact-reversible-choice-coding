# Paper 2 SOTA and external-validation analysis

Run class: SMOKE

## Fixed-model identical-fiber results

### Gimpel2025_cross_pool_public_codebook

- GCall_to_GCfix, r=2, P5_combined_context: gain=0.00234903 (0.2503 SD), 95% fixed-model fiber CI [0.00186977, 0.00280474], payload retention=1.0000.
- GCall_to_GCfix, r=2, published_1dcnn: gain=0.00100442 (0.1070 SD), 95% fixed-model fiber CI [0.00049683, 0.00155504], payload retention=1.0000.
- GCall_to_GCfix, r=2, published_sota_hard_filter: gain=0.00001564 (0.0017 SD), 95% fixed-model fiber CI [-0.00045650, 0.00056003], payload retention=1.0000.
- GCall_to_GCfix, r=4, P5_combined_context: gain=0.00318967 (0.3398 SD), 95% fixed-model fiber CI [0.00264833, 0.00373975], payload retention=1.0000.
- GCall_to_GCfix, r=4, published_1dcnn: gain=0.00123227 (0.1313 SD), 95% fixed-model fiber CI [0.00014719, 0.00222812], payload retention=1.0000.
- GCall_to_GCfix, r=4, published_sota_hard_filter: gain=0.00064176 (0.0684 SD), 95% fixed-model fiber CI [-0.00007099, 0.00126365], payload retention=1.0000.
- GCfix_to_GCall, r=2, P5_combined_context: gain=0.00213546 (0.2744 SD), 95% fixed-model fiber CI [0.00171119, 0.00242468], payload retention=1.0000.
- GCfix_to_GCall, r=2, published_1dcnn: gain=0.00041859 (0.0538 SD), 95% fixed-model fiber CI [0.00007883, 0.00082483], payload retention=1.0000.
- GCfix_to_GCall, r=2, published_sota_hard_filter: gain=-0.00062404 (-0.0802 SD), 95% fixed-model fiber CI [-0.00122350, -0.00002172], payload retention=1.0000.
- GCfix_to_GCall, r=4, P5_combined_context: gain=0.00296666 (0.3812 SD), 95% fixed-model fiber CI [0.00252262, 0.00372064], payload retention=1.0000.
- GCfix_to_GCall, r=4, published_1dcnn: gain=0.00100211 (0.1288 SD), 95% fixed-model fiber CI [0.00006368, 0.00188299], payload retention=1.0000.
- GCfix_to_GCall, r=4, published_sota_hard_filter: gain=-0.00046809 (-0.0602 SD), 95% fixed-model fiber CI [-0.00178064, 0.00050731], payload retention=1.0000.

### Gimpel2025_external_laboratory_Taq

- GCall_to_external_Taq, r=2, P5_combined_context: gain=0.00245844 (0.2345 SD), 95% fixed-model fiber CI [0.00212107, 0.00294796], payload retention=1.0000.
- GCall_to_external_Taq, r=2, published_1dcnn: gain=0.00161693 (0.1542 SD), 95% fixed-model fiber CI [0.00110964, 0.00214498], payload retention=1.0000.
- GCall_to_external_Taq, r=4, P5_combined_context: gain=0.00408539 (0.3896 SD), 95% fixed-model fiber CI [0.00330682, 0.00493276], payload retention=1.0000.
- GCall_to_external_Taq, r=4, published_1dcnn: gain=0.00160625 (0.1532 SD), 95% fixed-model fiber CI [0.00090999, 0.00230299], payload retention=1.0000.
- GCfix_to_external_Taq, r=2, P5_combined_context: gain=0.00239766 (0.2287 SD), 95% fixed-model fiber CI [0.00180172, 0.00273383], payload retention=1.0000.
- GCfix_to_external_Taq, r=2, published_1dcnn: gain=0.00161057 (0.1536 SD), 95% fixed-model fiber CI [0.00088198, 0.00213671], payload retention=1.0000.
- GCfix_to_external_Taq, r=4, P5_combined_context: gain=0.00343257 (0.3274 SD), 95% fixed-model fiber CI [0.00248389, 0.00421343], payload retention=1.0000.
- GCfix_to_external_Taq, r=4, published_1dcnn: gain=0.00138308 (0.1319 SD), 95% fixed-model fiber CI [0.00027755, 0.00216627], payload retention=1.0000.

### Gimpel2025_external_laboratory_Q5_sensitivity

- GCall_to_external_Q5, r=2, P5_combined_context: gain=-0.00571101 (-0.2512 SD), 95% fixed-model fiber CI [-0.00769253, -0.00347599], payload retention=1.0000.
- GCall_to_external_Q5, r=2, published_1dcnn: gain=-0.00023701 (-0.0104 SD), 95% fixed-model fiber CI [-0.00243449, 0.00206683], payload retention=1.0000.
- GCall_to_external_Q5, r=4, P5_combined_context: gain=-0.01859621 (-0.8179 SD), 95% fixed-model fiber CI [-0.02468242, -0.01154691], payload retention=1.0000.
- GCall_to_external_Q5, r=4, published_1dcnn: gain=0.00226958 (0.0998 SD), 95% fixed-model fiber CI [-0.00354270, 0.00641297], payload retention=1.0000.
- GCfix_to_external_Q5, r=2, P5_combined_context: gain=-0.00378485 (-0.1665 SD), 95% fixed-model fiber CI [-0.00639778, -0.00151593], payload retention=1.0000.
- GCfix_to_external_Q5, r=2, published_1dcnn: gain=0.00282629 (0.1243 SD), 95% fixed-model fiber CI [0.00049398, 0.00517744], payload retention=1.0000.
- GCfix_to_external_Q5, r=4, P5_combined_context: gain=-0.01186458 (-0.5218 SD), 95% fixed-model fiber CI [-0.01804461, -0.00418834], payload retention=1.0000.
- GCfix_to_external_Q5, r=4, published_1dcnn: gain=0.00208809 (0.0918 SD), 95% fixed-model fiber CI [-0.00370649, 0.00591210], payload retention=1.0000.

## Two-stage source-plus-target intervals

- GCall_to_GCfix, r=2: mean=0.00226031, 95% two-stage interval [0.00177024, 0.00273948], fraction above zero=1.0000.
- GCall_to_GCfix, r=4: mean=0.00309344, 95% two-stage interval [0.00242689, 0.00364912], fraction above zero=1.0000.
- GCall_to_external_Q5, r=2: mean=-0.00532892, 95% two-stage interval [-0.00821252, -0.00268344], fraction above zero=0.0000.
- GCall_to_external_Q5, r=4: mean=-0.01753839, 95% two-stage interval [-0.02286074, -0.01311146], fraction above zero=0.0000.
- GCall_to_external_Taq, r=2: mean=0.00250059, 95% two-stage interval [0.00218264, 0.00287240], fraction above zero=1.0000.
- GCall_to_external_Taq, r=4: mean=0.00398450, 95% two-stage interval [0.00301400, 0.00502811], fraction above zero=1.0000.
- GCfix_to_GCall, r=2: mean=0.00216357, 95% two-stage interval [0.00177414, 0.00247132], fraction above zero=1.0000.
- GCfix_to_GCall, r=4: mean=0.00304018, 95% two-stage interval [0.00248251, 0.00367269], fraction above zero=1.0000.
- GCfix_to_external_Q5, r=2: mean=-0.00335295, 95% two-stage interval [-0.00595912, -0.00121691], fraction above zero=0.0000.
- GCfix_to_external_Q5, r=4: mean=-0.01312118, 95% two-stage interval [-0.01825587, -0.00638703], fraction above zero=0.0000.
- GCfix_to_external_Taq, r=2: mean=0.00246764, 95% two-stage interval [0.00205172, 0.00284271], fraction above zero=1.0000.
- GCfix_to_external_Taq, r=4: mean=0.00354498, 95% two-stage interval [0.00261608, 0.00427038], fraction above zero=1.0000.

## Evidence boundary

retrospective measured public-sequence selection and an exact computational choice-codec interface; not prospective codec-output wet lab, mechanism, NUPACK, material, sequencing, recovery, full 110-nt weighted-language capacity, achieved storage density, end-to-end superiority, or unpublished-sister-manuscript evidence.
