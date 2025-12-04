function settings = recordScanImageSettings(API)
% Record settings from ScanImage
%
% function settings = recordScanImageSettings(API)
%
%
% Inputs
% API - silinker object
%
% Output
% settings - structure of settings
%
% See also
% reapplyScanImageSettings
%
% Rob Campbell, SWC AMF



    if API.versionGreaterThan('2020')
        settings.stackManShutterClose = API.hSI.hStackManager.closeShutterBetweenSlices;
        settings.numVolumes = API.hSI.hStackManager.numVolumes;
        settings.stackActuator = API.hSI.hStackManager.stackActuator;
        settings.stackManCentr = API.hSI.hStackManager.centeredStack;
        settings.stackManEnable = API.hSI.hStackManager.enable;
    else
        settings.fastZEnable = API.hSI.hFastZ.enable;
        settings.fastZNumVolumes = API.hSI.hFastZ.numVolumes;
        settings.stackManCentr = API.hSI.hStackManager.stackStartCentered;
        settings.stackManShutterClose = API.hSI.hStackManager.shutterCloseMinZStepSize;
        settings.stackManSlowFastZ = API.hSI.hStackManager.slowStackWithFastZ;  %This will be enabled
    end

    settings.zoomFactor = API.hSI.hRoiManager.scanZoomFactor;
    settings.pixelsPerLine = API.hSI.hRoiManager.pixelsPerLine;

    settings.fastZwaveform = API.hSI.hFastZ.waveformType;
    settings.stackManNumSlices = API.hSI.hStackManager.numSlices;
    settings.stackManStepSize = API.hSI.hStackManager.stackZStepSize;

    settings.loggingEnabled = API.hSI.hChannels.loggingEnable;
    settings.logAveFact = API.hSI.hScan2D.logAverageFactor;
    settings.framesPerSlice = API.hSI.hStackManager.framesPerSlice;
    settings.acqsPerLoop = API.hSI.acqsPerLoop;

    settings.logFileStem = API.hSI.hScan2D.logFileStem;
    settings.logFilePath = API.hSI.hScan2D.logFilePath;
    settings.logFileCounter = API.hSI.hScan2D.logFileCounter;

    settings.laserPower = API.hSI.hBeams.powers;

    settings.pmtGains = API.hSI.hPmts.gains;

    settings.extTrigEnable = API.hSI.extTrigEnable;

    settings.channelSubtractOffset = API.hSI.hChannels.channelSubtractOffset;
