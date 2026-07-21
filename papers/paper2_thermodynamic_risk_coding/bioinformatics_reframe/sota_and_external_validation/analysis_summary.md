# Paper 2 SOTA and external-validation analysis

Run class: FROZEN FULL

## Fixed-model identical-fiber results

### Gimpel2025_cross_pool_public_codebook

- GCall_to_GCfix, r=2, P5_combined_context: gain=0.00234903 (0.2503 SD), 95% fixed-model fiber CI [0.00191334, 0.00281727], payload retention=1.0000.
- GCall_to_GCfix, r=2, published_1dcnn: gain=0.00100442 (0.1070 SD), 95% fixed-model fiber CI [0.00044925, 0.00154628], payload retention=1.0000.
- GCall_to_GCfix, r=2, published_sota_hard_filter: gain=0.00001564 (0.0017 SD), 95% fixed-model fiber CI [-0.00057285, 0.00058053], payload retention=1.0000.
- GCall_to_GCfix, r=4, P5_combined_context: gain=0.00318967 (0.3398 SD), 95% fixed-model fiber CI [0.00252222, 0.00388960], payload retention=1.0000.
- GCall_to_GCfix, r=4, published_1dcnn: gain=0.00123227 (0.1313 SD), 95% fixed-model fiber CI [0.00005988, 0.00219037], payload retention=1.0000.
- GCall_to_GCfix, r=4, published_sota_hard_filter: gain=0.00064176 (0.0684 SD), 95% fixed-model fiber CI [-0.00031639, 0.00148880], payload retention=1.0000.
- GCfix_to_GCall, r=2, P5_combined_context: gain=0.00213546 (0.2744 SD), 95% fixed-model fiber CI [0.00166899, 0.00261517], payload retention=1.0000.
- GCfix_to_GCall, r=2, published_1dcnn: gain=0.00041859 (0.0538 SD), 95% fixed-model fiber CI [-0.00005728, 0.00091215], payload retention=1.0000.
- GCfix_to_GCall, r=2, published_sota_hard_filter: gain=-0.00062404 (-0.0802 SD), 95% fixed-model fiber CI [-0.00135301, 0.00003188], payload retention=1.0000.
- GCfix_to_GCall, r=4, P5_combined_context: gain=0.00296666 (0.3812 SD), 95% fixed-model fiber CI [0.00230641, 0.00358140], payload retention=1.0000.
- GCfix_to_GCall, r=4, published_1dcnn: gain=0.00100211 (0.1288 SD), 95% fixed-model fiber CI [-0.00000743, 0.00190630], payload retention=1.0000.
- GCfix_to_GCall, r=4, published_sota_hard_filter: gain=-0.00046809 (-0.0602 SD), 95% fixed-model fiber CI [-0.00174509, 0.00071138], payload retention=1.0000.

### Gimpel2025_external_laboratory_Taq

- GCall_to_external_Taq, r=2, P5_combined_context: gain=0.00245844 (0.2345 SD), 95% fixed-model fiber CI [0.00194486, 0.00301679], payload retention=1.0000.
- GCall_to_external_Taq, r=2, published_1dcnn: gain=0.00161693 (0.1542 SD), 95% fixed-model fiber CI [0.00110589, 0.00219093], payload retention=1.0000.
- GCall_to_external_Taq, r=4, P5_combined_context: gain=0.00408539 (0.3896 SD), 95% fixed-model fiber CI [0.00317002, 0.00496631], payload retention=1.0000.
- GCall_to_external_Taq, r=4, published_1dcnn: gain=0.00160625 (0.1532 SD), 95% fixed-model fiber CI [0.00081720, 0.00236793], payload retention=1.0000.
- GCfix_to_external_Taq, r=2, P5_combined_context: gain=0.00239766 (0.2287 SD), 95% fixed-model fiber CI [0.00188645, 0.00296526], payload retention=1.0000.
- GCfix_to_external_Taq, r=2, published_1dcnn: gain=0.00161057 (0.1536 SD), 95% fixed-model fiber CI [0.00103984, 0.00218680], payload retention=1.0000.
- GCfix_to_external_Taq, r=4, P5_combined_context: gain=0.00343257 (0.3274 SD), 95% fixed-model fiber CI [0.00236462, 0.00442358], payload retention=1.0000.
- GCfix_to_external_Taq, r=4, published_1dcnn: gain=0.00138308 (0.1319 SD), 95% fixed-model fiber CI [0.00026556, 0.00240606], payload retention=1.0000.

### Gimpel2025_external_laboratory_Q5_sensitivity

- GCall_to_external_Q5, r=2, P5_combined_context: gain=-0.00571101 (-0.2512 SD), 95% fixed-model fiber CI [-0.00817519, -0.00324061], payload retention=1.0000.
- GCall_to_external_Q5, r=2, published_1dcnn: gain=-0.00023701 (-0.0104 SD), 95% fixed-model fiber CI [-0.00285817, 0.00224254], payload retention=1.0000.
- GCall_to_external_Q5, r=4, P5_combined_context: gain=-0.01859621 (-0.8179 SD), 95% fixed-model fiber CI [-0.02485317, -0.01214992], payload retention=1.0000.
- GCall_to_external_Q5, r=4, published_1dcnn: gain=0.00226958 (0.0998 SD), 95% fixed-model fiber CI [-0.00260955, 0.00682206], payload retention=1.0000.
- GCfix_to_external_Q5, r=2, P5_combined_context: gain=-0.00378485 (-0.1665 SD), 95% fixed-model fiber CI [-0.00640725, -0.00118606], payload retention=1.0000.
- GCfix_to_external_Q5, r=2, published_1dcnn: gain=0.00282629 (0.1243 SD), 95% fixed-model fiber CI [0.00035891, 0.00512888], payload retention=1.0000.
- GCfix_to_external_Q5, r=4, P5_combined_context: gain=-0.01186458 (-0.5218 SD), 95% fixed-model fiber CI [-0.01908351, -0.00492572], payload retention=1.0000.
- GCfix_to_external_Q5, r=4, published_1dcnn: gain=0.00208809 (0.0918 SD), 95% fixed-model fiber CI [-0.00315853, 0.00678960], payload retention=1.0000.

## Two-stage source-plus-target intervals

- GCall_to_GCfix, r=2: mean=0.00226033, 95% two-stage interval [0.00175426, 0.00276960], fraction above zero=1.0000.
- GCall_to_GCfix, r=4: mean=0.00314820, 95% two-stage interval [0.00246768, 0.00383601], fraction above zero=1.0000.
- GCall_to_external_Q5, r=2: mean=-0.00580474, 95% two-stage interval [-0.00879021, -0.00299124], fraction above zero=0.0000.
- GCall_to_external_Q5, r=4: mean=-0.01791040, 95% two-stage interval [-0.02522317, -0.01049660], fraction above zero=0.0000.
- GCall_to_external_Taq, r=2: mean=0.00249308, 95% two-stage interval [0.00194047, 0.00303514], fraction above zero=1.0000.
- GCall_to_external_Taq, r=4: mean=0.00398696, 95% two-stage interval [0.00300744, 0.00495754], fraction above zero=1.0000.
- GCfix_to_GCall, r=2: mean=0.00214864, 95% two-stage interval [0.00168748, 0.00259633], fraction above zero=1.0000.
- GCfix_to_GCall, r=4: mean=0.00303532, 95% two-stage interval [0.00233101, 0.00376684], fraction above zero=1.0000.
- GCfix_to_external_Q5, r=2: mean=-0.00395548, 95% two-stage interval [-0.00694793, -0.00113082], fraction above zero=0.0010.
- GCfix_to_external_Q5, r=4: mean=-0.01279110, 95% two-stage interval [-0.02135270, -0.00493600], fraction above zero=0.0005.
- GCfix_to_external_Taq, r=2: mean=0.00245358, 95% two-stage interval [0.00192216, 0.00301382], fraction above zero=1.0000.
- GCfix_to_external_Taq, r=4: mean=0.00349036, 95% two-stage interval [0.00235988, 0.00450761], fraction above zero=1.0000.

## Statistical boundary

Source bootstrap copies retained fixed original-sequence CV folds, so no original sequence occurred in both training and validation within a replicate. Intervals condition on the published sequence-level outcomes and do not propagate technical-replicate, batch, normalization or source-assay measurement uncertainty.

## Evidence boundary

retrospective measured public-sequence selection and an exact computational choice-codec interface; not prospective codec-output wet lab, mechanism, NUPACK, material, sequencing, recovery, full 110-nt weighted-language capacity, achieved storage density, end-to-end superiority, or unpublished-sister-manuscript evidence.
