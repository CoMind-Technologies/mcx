# COMIND SPECIFIC README 
 - Author: Alex Antrobus
 - Last Update: 21/04/2021
This document lays out which of these directories is specific to COMIND, namely:
* `/scripts` - all scripts for building volumes, parameter sets, configuration files, etc.

* `/volumes` - .mat format volume files for mcx. These should all be 3-dimensional cell arrays with uintX (for X=8,16,32 etc.).

* `/configs` - .mat configuration files for simulations. Ideally, these should be self-contained and fully define a simulation.
	That said, the largest aspect of these will generally be volumetric data, and in some cases (e.g. replay simulations) fluence data detected photon data,
	all of which roughly scale with volume and duration.
	For these reasons, you may wish to design simulations to call volume files separately from the config files.

## Getting mcxlab working
Building mcxlab from source (this source) has proved somewhat of a nightmare, and it is generally easier to simply download a pre-compiled build for your OS from [SourceForge](https://sourceforge.net/projects/mcx/files/). For safety, you may want to pick the newest build preceding date-stamp on this README. _(Currently *2020Furious Fermion*_).

Once you've got the binaries, simply point your path to both them AND the relevant directories in this repo (e.g. ./mcxlab/). Note that this is true for your Matlab installation too.
