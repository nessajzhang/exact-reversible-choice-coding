# Paper 2 clean-Linux evidence

This directory contains the downloaded text evidence for the successful Ubuntu 24.04 portability run of commit `7d192be6a2a143750805444de3a0e3d92283976e`.

- Workflow run: `https://github.com/nessajzhang/exact-reversible-choice-coding/actions/runs/29942185585`
- Frozen-release job: `88998562611`, started `2026-07-22T17:24:24Z`, completed `2026-07-22T17:26:28Z`, conclusion `success`.
- Public-input job: `88998562688`, started `2026-07-22T17:24:24Z`, completed `2026-07-22T17:36:24Z`, conclusion `success`.
- Frozen artifact ID: `8538691007`.
- Public-input artifact ID: `8538980112`.

The logs are the stdout artifacts uploaded by the workflow. The GitHub run remains the authoritative record for combined stdout, stderr, step status and annotations. The downloaded files preserve the contract-test, integrity, frozen-reproduction, compilation, consistency, PDF metadata/font and full public-input completion evidence used in the portability audit.

The preceding probe at run `29941571798` is intentionally not accepted as public-input evidence: its missing-script failure was masked by a pipeline that did not propagate the producer's exit status. The corrected workflow uses `pipefail`, and the successful run above executed for 12 minutes and reached the declared full-audit completion markers.
