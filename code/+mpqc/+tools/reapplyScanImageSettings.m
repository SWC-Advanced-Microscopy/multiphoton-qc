function settings = reapplyScanImageSettings(API,settings,verbose)
% Record settings from ScanImage
%
% function settings = reapplyScanImageSettings(API,settings)
%
%
% Inputs
% API - silinker object
% settings - output for recordScanImageSettings
%
% Output
% settings - structure of settings
%
% See also
% recordScanImageSettings
%
%
% Rob Campbell, SWC AMF


    if nargin<3
        verbose=false;
    end

    if verbose
        fprintf('Re-applying ScanImage Settings.\n')
    end


    if API.versionGreaterThan('2020')
        API.hSI.hStackManager.closeShutterBetweenSlices = settings.stackManShutterClose;
        API.hSI.hStackManager.numVolumes = settings.numVolumes;
        API.hSI.hStackManager.stackActuator = settings.stackActuator;
        API.hSI.hStackManager.enable = settings.stackManEnable;
        API.hSI.hStackManager.centeredStack = settings.stackManCentr;
    else
        API.hSI.hFastZ.enable = settings.fastZEnable;
        API.hSI.hFastZ.numVolumes = settings.fastZNumVolumes;
        API.hSI.hStackManager.stackStartCentered = settings.stackManCentr;
        API.hSI.hStackManager.shutterCloseMinZStepSize =  settings.stackManShutterClose;
        API.hSI.hStackManager.slowStackWithFastZ = settings.stackManSlowFastZ;
    end

    API.hSI.hRoiManager.scanZoomFactor = settings.zoomFactor;
    API.hSI.hRoiManager.pixelsPerLine = settings.pixelsPerLine;

    API.hSI.hFastZ.waveformType = settings.fastZwaveform;
    API.hSI.hStackManager.numSlices = settings.stackManNumSlices;
    API.hSI.hStackManager.stackZStepSize = settings.stackManStepSize;

    API.hSI.hChannels.loggingEnable = settings.loggingEnabled;
    API.hSI.hScan2D.logAverageFactor = settings.logAveFact;
    API.hSI.hStackManager.framesPerSlice = settings.framesPerSlice;
    API.hSI.acqsPerLoop = settings.acqsPerLoop;

    API.hSI.hScan2D.logFileStem = settings.logFileStem;
    API.hSI.hScan2D.logFilePath = settings.logFilePath;
    API.hSI.hScan2D.logFileCounter = settings.logFileCounter;

    API.hSI.hBeams.powers = settings.laserPower;

    API.hSI.hPmts.gains = settings.pmtGains;

    API.hSI.extTrigEnable = settings.extTrigEnable;

    API.hSI.hChannels.channelSubtractOffset = settings.channelSubtractOffset;
