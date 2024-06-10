# fp-trace-analysis

A small python program designed to read pClamp .abf files, and automatically perform corrections and analysis as performed in Baird et al., 2023 (https://doi.org/10.1007/s00213-023-06340-8), in an attempt to simplify the workflow of analyzing these data.

Used for dLight 1.1 and GRABDA2m with experimental signal at 470 and isosbestic signal at 405 to control for movement artifacts.
Signals are recorded in 15-second "sweeps," to reduce intra-trial bleaching. Each sweep contained 50,000 points.

In essence, these corrections include:
1. Calculating the average baseline fluorescence of the fiberoptic cables and subtracting that from the experimental signal.
2. Gaussian filtering the 405 nm signal to reduce noise when the ratio is taken.
3. Finding the ratio of the 470/405 nm signals to eliminate movement artifacts.
4. Gaussian filtering the combined signal to reduce noise.

Current Analyses Include:
1. Finding the average fluorescence value of each 15s trace, then calculating the ΔF/F ((trace average - pre injection average)/pre injection average)
2. Frequency of events in each trace (Hz)
3. Duration of event decay (ms)
4. Event height or amplitude (V)
5. Event width at 50% peak (ms)

Events are defined as any peak in a trace with a prominence of 0.05, with a maximum window length of 10,000 points.