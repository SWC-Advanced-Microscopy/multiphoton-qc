function lens_paper(varargin)
    % Image lens paper to generate raw data for assessing overall image quality
    %
    % function mpqc.record.lens_paper('power',value, 'wavelength', value, ...)
    %
    % Purpose
    % Lens paper (e.g. purchased from ThorLabs) is useful standard sample. It is
    % consistent from batch to batch, it contains structure, it bleaches slowly. Images
    % of lens paper show a sparse structure of fibres that do not dim rapidly with depth.
    % depth and anyway it is intrinsically pretty thin. Lens paper can easily be sealed
    % under a coverslip but note the paper must remain dry for imaging. Lens paper samples
    % can be assessed in two ways:
    % 1. Qualitatively. Here the user learns by experience that, say, a barely reasonable
    %   image of the lens tissue paper can be achieved at a power of 15 mW at the sample.
    %   Failure to get a decent image at that power would indicate a potential problem.
    % 2. Quantitatively. Here we would use a higher laser power, such as 50 to 100 mW, and
    %   analyse the data post-acquisition to calculate the mean number of photons emitted
    %   from the sample per pixel. If the same laser power is used each time, a deviation
    %   in mean photon count may indicate a problem with the microscope. A higher laser
    %   power is needed to provide a larger SNR for the measurement.
    %
    %
    % Inputs (optional param/val pairs. If not defined, a CLI prompt appears)
    %  'wavelength' - Excitation wavelength of the laser. Defined in nm.
    %  'power' - Power at the sample. Defined in mW.
    %
    % Inputs (optional, exact & case sensitive,  param/val pairs)
    %  'nFames' - The number of frames to acquire. If not defined a reasonable value is
    %           chosen: the larger of 40 frames or 5 seconds of data.
    %  'numGains' - The number of PMT gains at which to record the data. If not defined,
    %               the current PMT gains are used and one set of images obtained only. If
    %               "numGains" is non empty and >1 then a series of recordings at
    %               different gains is made. NOTE: set numGains to -1 to use the default.
    %
    %
    % The channels to image is determined based on the channel selected to be saved
    % within ScanImage.
    %
    % Running the function
    % The wavelength and power inputs are required for the function to run. The user
    % must supply power and wavelengths as integer for logging purposes. If the user
    % does not define these as parameter/value pairs, they are prompted at CLI to do so.
    % The order of the two arguments does not matter.
    %
    %
    % e.g. In the following example the user is imaging at 920 nm and 10 mW
    % at the sample.
    %
    % >> mpqc.record.lens_paper('wavelength',920,'power',10)
    %
    % e.g. If no inputs are given, user will be prompted for values
    % >> mpqc.record.lens_paper
    %
    % Rob Campbell, SWC AMF, initial commit 2022


    %%
    % process inputs
    out =  parseInputVariable(varargin{:});

    % Extract critical input variable
    laser_wavelength=out.wavelength;
    laser_power_in_mW = out.power;

    % Extract optional inputs
    numFramesToAcquire = [];
    if isfield(out,'nFrames')
        numFramesToAcquire = out.nFrames;
    end

    numGains = [];
    if isfield(out,'mumGains')
        mumGains = out.mumGains;
    end

    % Connect to ScanImage using the linker class
    API = sibridge.silinker;

    if API.linkSucceeded == false
        return
    end


    % Create 'diagnostic' directory in the user's desktop
    saveDir = mpqc.tools.makeTodaysDataDirectory;
    if isempty(saveDir)
        return
    end


    % Get measurePSF settings for use later
    SETTINGS=mpqc.settings.readSettings;

    %Record the state of all ScanImage settings we will change so we can change them back
    initialSettings = mpqc.tools.recordScanImageSettings(API);


    %Apply common settings
    API.setZSlices(1)

    API.hSI.hRoiManager.pixelsPerLine=256;

    % If the user did not define how many frames to acquire, we will acquire 5 seconds
    % of data or 40 frames, whichever is larger.
    if isempty(numFramesToAcquire)
        minFrames = 40;
        fiveSecondsOfFrames = ceil(5 / API.hSI.hRoiManager.scanFramePeriod);
        if fiveSecondsOfFrames < 40
            numFramesToAcquire = minFrames;
        else
            numFramesToAcquire = fiveSecondsOfFrames;
        end
    end

    % TODO -- make it acquire on all available PMTS?

    % Do not subtract channel offsets. This probably matters for using lens paper data
    % to calibrate the standard source.
    API.disableChannelOffsetSubtraction;

    API.hSI.hStackManager.framesPerSlice=numFramesToAcquire;

    API.hSI.hScan2D.logAverageFactor = 1;
    API.hSI.hDisplay.volumeDisplayStyle='Current';

    API.hSI.hRoiManager.scanZoomFactor = 2; % To keep vignetting down

    API.hSI.hChannels.loggingEnable=true;


    API.hSI.acqsPerLoop=1;

    if isempty(numGains)
        % If no gain was defined, we simply acquire with existing
        % PMT gains

        % TODO -- warn if gains across detectors are different. Give option to bail.
        % TODO -- read PMT gains and set file name according to the first active PMT's gain.
        % TODO -- maybe prompt to set the gain?

        % Set file name and save dir
        fileStem = sprintf('%s_lens_paper_%dnm_%dmW__%s', ...
            SETTINGS.microscope.name, ...
            laser_wavelength, ...
            laser_power_in_mW, ...
            datestr(now,'yyyy-mm-dd_HH-MM-SS'));

        API.hSI.hScan2D.logFileStem=fileStem;
        API.hSI.hScan2D.logFilePath=saveDir;
        API.hSI.hScan2D.logFileCounter=1;

        % Start the acquisition
        API.acquireAndWait;
    else
        % If specific gains were requested, loop through these
        % one at a time.

        gainsToTest = getPMTGainsToTest(numGains);

        for ii=1:length(gainsToTest)

            % In case a zero gain is specified across all detectors, do not record it.
            if sum(gainsToTest(:,ii)) == 0
                continue
            end

            % Set file name and save dir
            fileStem = sprintf('%s_lens_paper_%dV_%dnm_%dmW__%s', ...
                SETTINGS.microscope.name, ...
                gainsToTest(1,ii), ...
                laser_wavelength, ...
                laser_power_in_mW, ...
                datestr(now,'yyyy-mm-dd_HH-MM-SS'));

            API.hSI.hScan2D.logFileStem=fileStem;
            API.hSI.hScan2D.logFilePath=saveDir;
            API.hSI.hScan2D.logFileCounter=1;

            API.setPMTgains(gainsToTest(:,ii)); % Set gain
            pause(1) % Wait for settling

            API.acquireAndWait;
        end

    end


    % Report saved file location and copy mpqc settings there
    postAcqTasks(saveDir,fileStem)

    mpqc.tools.reapplyScanImageSettings(API,initialSettings);

end
