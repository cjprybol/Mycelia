# Momentum Fork-Resolution Weight Ablation

Synthetic repeat benchmark comparison for uniform, linear, saturating, and LLR/mapping-quality weighting.

Synthetic NGA50 is a repeat-resolution contiguity proxy derived from the benchmark's aligned-span model.

| Weight mode | alpha=beta | Fork accuracy | Correct resolution rate | Misassembly rate | Synthetic NGA50 (bp) | Synthetic NGA50 / reference | Resolution rate | Mean observations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| linear | 0.2 | 0.836 | 0.833 | 0.164 | 480.0 | 0.837 | 0.997 | 1.34 |
| llr_mapping_quality | 0.2 | 0.903 | 0.679 | 0.073 | 480.0 | 0.837 | 0.752 | 4.18 |
| saturating | 0.2 | 0.845 | 0.782 | 0.144 | 480.0 | 0.837 | 0.926 | 2.31 |
| uniform | 0.2 | 1.0 | 0.019 | 0.0 | 240.0 | 0.419 | 0.019 | 6.63 |
