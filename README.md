# Beta-cell-functional-network


This repository contains MATLAB code for constructing temporal functional networks from calcium imaging signals using wavelet coherence in the fast frequency band.

## Input data

The script expects:

- `data_jki_2_ds_5.txt`: matrix of calcium signals with dimensions `[time x cells]`
- `koordinate_jki_2.txt`: matrix of cell coordinates `[cells x 2]`
- `sampling.txt`: scalar containing the original sampling interval

## Main parameters

- analyzed cells: `2:169`
- global interval: `1000-2800 s`
- fast band: `0.02-0.15 Hz`
- window length: `60 s`
- step size: `30 s`
- retained connections per window: top `3%`
- persistence threshold: `0.20`
