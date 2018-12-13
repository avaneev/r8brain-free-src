## Benchmark ##

This folder includes test sample rate conversion tool for Win64 and test
`DrumsSrc.wav` and `DrumsDst96.wav` files. The `DrumsDst96.wav` was produced
by the author via the `r8bfreesrc.exe` tool from the source `DrumsSrc.wav`
file. This converted file can be used for comparison with your own sample rate
converter implementation using this library.

The `rmscompare` tool can be used to compare differences between two WAV
files. E.g. two 24-bit files with the magnitude of differences of -141 dB or
below can be considered equal.
