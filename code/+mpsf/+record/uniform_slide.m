function uniform_slide(varargin)
    % Set up ScanImage to acquire images of a uniform Chroma slide
    %
    % function mpsf.record.uniform_slide('param1','val1')
    %
    % Purpose
    % The uniform Chroma slide is a bright, uniform, fluorescent target that can be used
    % assess image brightness across the field of view. Typically the Fov is inhomogeneous,
    % with the brightest signal at the centre and dimming towards the edges and, in
    % particular, the corners. This effect is of course most visible at the maximum
    % field of view and rapidly disappears with zoom. The dimming at the edges is most
    % likely due to vignetting of the excitation beam, but it could also be due to
    % increased aberrations of the PSF at the field edges or to vignetting within the
    % detection path.
    %
    % One indication of a misaligned microscope is that the uniform slide produces an
    % image where the maximum brightness peak is far from the centre. A microscope with
    % poorly conjugated scan optics can produce very oblate intensity pattern.
    %
    %
    % Inputs
    % The user must supply the laser power in mW and the laser wavelength in nm as
    % parameter/value pairs: 'wavelength' and 'power'. This is for logging purposes: the
    % values are used to build the file name. If the user fails to do this they are
    % prompted at CLI.
    %
    %
    % Example
    % 1. Select one channel to save in ScanImage
    % 2. Choose a low laser power and ensure you get a clear image of the slide.
    % 3. Run this function as shown in the examples below
    %
    % >> mpsf.record.uniform_slide
    % >> mpsf.record.uniform_slide('wavelength',920,'power',10)
    % >> mpsf.record.uniform_slide('power',7,'wavelength',800)
    %
    %
    % Rob Campbell, first committed SWC 2022
    % Updated: Isabell Whiteley, SWC 2024


    % Parse the input arguments
    out =  parseInputVariable(varargin{:});
    laser_wavelength=out.wavelength;
    laser_power_in_mW = out.power;


    % Connect to ScanImage using the linker class
    API = sibridge.silinker;

    if API.linkSucceeded == false
        return
    end


    if ~strcmp(API.scannerType,'resonant')
        fprintf('%s is only safe with a resonant scanner. Quitting. \n',mfilename)
        return
    end

    % Errors if the user has selected more than one channel to save
    saveChanName = API.getSaveChannelName;
    if isempty(saveChanName)
        % message produced by API.getSaveChannelName;
        return
    end


    % Create 'diagnostic' directory in the user's desktop
    saveDir = mpsf.tools.makeTodaysDataDirectory;
    if isempty(saveDir)
        return
    end

    %Record the state of all ScanImage settings we will change so we can change them back
    settings = mpsf.tools.recordScanImageSettings(API);


    %Apply common setting
    API.setZSlices(1)

    API.hSI.hStackManager.framesPerSlice=20; % We will record multiple frames
    API.hSI.hRoiManager.pixelsPerLine=256;

    API.hSI.hScan2D.logAverageFactor = 1;
    API.hSI.hDisplay.volumeDisplayStyle='Current';

    API.hSI.hRoiManager.scanZoomFactor = 1; % Set zoom

    API.hSI.hChannels.loggingEnable=true;


    % Set file name and save dir
    SETTINGS=mpsf.settings.readSettings;
    fileStem = sprintf('%s_uniform_slide_zoom__%s_%dnm_%dmW_%s__%s', ...
        SETTINGS.microscope.name, ...
        obj.returnZoomFactorAsString, ...
        laser_wavelength, ...
        laser_power_in_mW, ...
        saveChanName, ...
        datestr(now,'yyyy-mm-dd_HH-MM-SS'));

    API.hSI.hScan2D.logFileStem=fileStem;
    API.hSI.hScan2D.logFilePath=saveDir;
    API.hSI.hScan2D.logFileCounter=1;

    API.acquireAndWait;


    mpsf.tools.reapplyScanImageSettings(API,settings);

    % Report where the file was saved
    mpsf.tools.reportFileSaveLocation(saveDir,fileStem)

    % Save system settings to this location
    settingsFilePath = mpsf.settings.findSettingsFile;
    copyfile(settingsFilePath, saveDir)
