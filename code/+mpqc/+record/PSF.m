function varargout = PSF(varargin)
    % Record a z-stack of beads to later generate a PSF
    %
    %  pathToTiff = mpqc.record.PSF('param1','val1', ...)
    %
    % Purpose
    % This function simplifies setting up ScanImage to acquire a z-stack of nano beads.
    % The user defines the depth of the stack in microns and the separation between
    % planes within this function and ScanImage is set up appropriately. The image stack
    % is averaged using the number of frames entered in the ScanImage IMAGE CONTROLS
    % window. The path to the saved data is displayed to screen and can optionally be
    % returned as an output argument.
    %
    % Inputs
    % All inputs are optional parameter/value pairs. The user will be interactively
    % prompted at the CLI to fill in any undefined values.
    % 'depthMicrons'      - Total depth to image in microns
    % 'stepSizeInMicrons' - Number of microns between each z step
    % 'wavelength'        - Excitation wavelength in nm
    % 'power'             - Excitation power at the sample in mW
    %
    %
    % Outputs
    % Optionally returns path to the TIFF.
    %
    %
    % Examples
    % 1) Record a 12 micron stack, and prompts user for other inputs.
    % >> record.PSF('depthMicrons',12)
    %
    % 2) Record a 20 micron stack every 0.5 microns and return path to tiff
    % >> F=record.PSF('depthMicrons',20,'stepSizeInMicrons',0.5);
    %
    %
    % Rob Campbell, SWC AMF, initial commit Nov 2018



    out =  parseInputVariable(varargin{:});
    laser_wavelength=out.wavelength;
    laser_power_in_mW = out.power;
    micronsToImage = out.depthMicrons;
    stepSizeInMicrons = out.stepSize;


    % Connect to ScanImage using the linker class
    API = sibridge.silinker;

    if API.linkSucceeded == false
        return
    end


    if length(API.hSI.hChannels.channelSave) > 1
        fprintf('Select just one channel to save\n')
        return
    end

    % Create 'PSF' directory in the user's desktop
    saveDir = mpqc.tools.makeTodaysDataDirectory;
    if isempty(saveDir)
        return
    end


    %Record the state of all ScanImage settings we will change so we can change them back
    initialSettings = mpqc.tools.recordScanImageSettings(API);

    % We will set up ScanImage to acquire the z-stack
    framesToAverage = API.hSI.hDisplay.displayRollingAverageFactor;
    numSlices = round(micronsToImage/stepSizeInMicrons);
    fileStem = sprintf('PSF__%d_nm__%d_mW__%s', ....
            laser_wavelength, ...
            laser_power_in_mW, ...
            datestr(now,'yyyy-mm-dd_HH-MM-SS'));

    try

        if API.versionGreaterThan('2020')
            API.hSI.hStackManager.closeShutterBetweenSlices = false;
            API.hSI.hStackManager.numVolumes = 1;
            API.hSI.hStackManager.stackActuator = 'fastZ';
            API.hSI.hStackManager.centeredStack = 0;
            API.hSI.hStackManager.enable = true;
        else
            API.hSI.hFastZ.enable=false;
            API.hSI.hFastZ.numVolumes=1;
            API.hSI.hStackManager.stackStartCentered=false; %TODO: SI doesn't work correctly when true
            API.hSI.hStackManager.shutterCloseMinZStepSize=stepSizeInMicrons+1;
            API.hSI.hStackManager.slowStackWithFastZ=true;
        end
        %%
        API.hSI.hFastZ.waveformType='step';


        API.hSI.hStackManager.numSlices=numSlices;
        API.hSI.hStackManager.stackZStepSize=stepSizeInMicrons;

        API.hSI.hChannels.loggingEnable=true;

        API.hSI.hScan2D.logFileStem=fileStem;
        API.hSI.hScan2D.logFilePath=saveDir;
        API.hSI.hScan2D.logFileCounter=1;

        API.hSI.hStackManager.framesPerSlice = framesToAverage;
        API.hSI.hScan2D.logAverageFactor = framesToAverage;
        API.hSI.hScan2D.bidirectional = 0; % disable bidi to avoid bidi artifacts in bead stacks

        API.hSI.hDisplay.volumeDisplayStyle='Current';
    catch ME
        %If something went wrong we revert the scan settings
        fprintf('Failed to set scan settings\n')
        fprintf(ME.message)
        return
    end

    % Start the acquisition and wait for it to finish
    API.acquireAndWait;

    % Report saved file location and copy mpqc settings there
    postAcqTasks(saveDir,fileStem)




    if nargout>0
        varargout{1} = pathToTiff;
    end

    mpqc.tools.reapplyScanImageSettings(API,initialSettings);

end % record.PSF
