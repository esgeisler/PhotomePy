# fp-trace-analysis

A small python program designed to read pClamp .abf files, and automatically perform corrections and analysis as performed in Baird et al., 2023 (https://doi.org/10.1007/s00213-023-06340-8)

Used for dLight 1.1 and GRABDA2m with experimental signal at 470 and isosbestic signal at 405 to control for movement artifacts.

In essence, these corrections include:
1. Calculating the average baseline fluorescence of the fiberoptic cables and subtracting that from the experimental signal.
2. Gaussian filter the 405 nm signal.
3. Find ratio of the 470/405 nm signals.
4. Gaussian filter the combined signal.

Current Analyses Include:
1. Finding the average of each 15s trace, then calculating the ΔF/F ((trace avg - pre inj avg)/pre inj avg)

Planned Analysis Include:
1. Frequency of events
2. Duration of event decay
3. Event height
4. Event width
5. Area under the curve of the ΔF/F