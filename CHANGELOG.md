
# Change-Log

### 2025/12/94
* Allow serial number of equipment to be numeric only. 
* Merge dev into master


### 2025/12/01
* Major improvements to power measurement to allow users to specify which beams to calibrate, if there are multiple beams configured. Sumiya Kuroda [Issues #108](https://github.com/SWC-Advanced-Microscopy/multiphoton-qc/issues/108).

### 2025/06/18
* UPDATE: Convert power measurement code to a class and add features to enable automation 
     of power measurements. Streamline how the interface to the power meter is handled. RAAC


### 2025/03/07
* NEW: Add facility for checking calibrated power at the sample plane automatically. Isabell Whiteley [PR #95]([https://github.com/SWC-Advanced-Microscopy/multiphoton-qc/pull/95).



### 2025/02/27
* Run lens paper and standard source without ScanImage offset subtraction. 


### 2025/02/24
* NEW: filter uniform slide more aggressively, which is needed for cases with uniform images
* BUGFIX: correctly display photon counts in lens paper data where multiple channels were acquired at once.
* BUGFIX: if weighted fit fails do a Huber fit.
* BUGFIX: correctly loads ScanImage TIFFs from systems with 2 beams
* BUGFIX: deal nicely with low photon counts.
* Minor improvements to lens paper recorder.


### 2025/02/17
Merge of a bunch of recent changes by Rob Campbell
* BREAKING CHANGE! Renamed repo to multiphoton-qc, mpsf -> mpqc. Your settings file will
  now need to be renamed `MPQC_SystemSettings.yml`.
* General doc text tidying. 
* Add `CONTRIBUTION_GUIDELINES.md`
* Made a standalone `CHANGELOG.md` and tidied it and improve formatting.
* Lens paper function now saves all channels.
* Standard source and lens paper produce more similar data (same pixel size)
* Dark noise is longer recorded along with electrical noise.
* Electrical noise function now saves all available channels automatically. 
* Added simple code in `tools` for converting the standard source to photons. 


### 2024/07/19
* NEW FEATURE: `mpsf.record` functions now all accept parameter/value pairs via standard interface.
Inputs that are required not supplied when the function is called are requested interactively at the CLI. Isabell Whiteley [PR #70](https://github.com/SWC-Advanced-Microscopy/measurePSF/pull/70).


### 2024/07/05
* Updates to standard light source. Plotting of said. Bugfixes.

### 2024/06/14
* NEW FEATURE: Add standard light source function. Add dark noise to electrical noise.

### 2024/05/23
* Implement a more elaborate microscope settings (parameters) system.

### 2023/07/31
* NEW FEATURE: Integrate functionality of making PDF reports, uniform slide analyses, and plots of lens paper. 

### 2022/08/02
* NEW FEATURE: Add function for imaging electrical noise and document protocol.

### 2022/08/01
* NEW FEATURE: Add functions for recording lens paper and uniform slides.

### 2020/02/19
* Add tiff stack name to title of top right plot.

### 2020/02/18
* Tidy `measurePSF` PDF and add dummy values to demo mode.

### 2020/02/18
* Bug fixes, check coarse z acquisition works, add PDF saving.

### 2020/02/12
* Bugfixes

### 2020/01/30
* Add `mpsf_tools.meanFrame` for displaying a rolling frame average.

### 2020/01/14
* Add button that allows the current image to be saved to the desktop.

### 2020/01/14
* Add edit boxes and checkboxes to allow the user to modify on the fly what would otherwise have been input arguments.

### 2020/01/14
* Get voxel size from ScanImage TIFF header.

### 2020/01/14
* If no input args to `measurePSF`, bring up the load GUI.

### 2020/01/13
* Convert Grid2MicsPerPixel to a class and add buttons to interact with SI.

### 2020/01/08
* Grid2MicsPerPixel optionally can extract the grid image directly from ScanImage.

### 2018/11/09
* Add `record.PSF`.

### 2017/11/28
* Simple GUI for interactive cropping of a desired bead.

### 2017/11/28
* Improve output data and don't display FWHM for directions in which the user defined no microns per pixel.

### 2017/11/27
* Convert `measurePSF` to a class so adding new features is easier.



