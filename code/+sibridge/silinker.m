classdef silinker < handle
    % Linking to the ScanImage API
    %
    % Performs some useful operations as well as exposing the ScanImage API
    %
    %
    % sitools.silinker

    properties (Hidden)
        scanimageObjectName = 'hSI' % If connecting to ScanImage look for this variable in the base workspace
        hSI % The ScanImage API attaches here
        listeners = {} % Reserved for listeners we might make
        linkSucceeded % true if SI connected
    end % Close hidden properties


    methods

        function obj = silinker(connectToSI)
            % By default connect to ScanImage on startup
            if nargin<1
                connectToSI=true;
            end

            if connectToSI
                obj.linkToScanImageAPI;
            end

        end % Constructor


        function delete(obj)
            obj.hSI=[];
        end % Destructor

        function success = linkToScanImageAPI(obj)
            % Link to ScanImage API by importing from base workspace and
            % copying handling to obj.hSI

            success=false;


            API = sibridge.getSIobject;
            if isempty(API)
                obj.linkSucceeded = false;
                return
            end

            obj.hSI=API; % Make composite object
            obj.linkSucceeded = true;
            success=true;
        end % linkToScanImageAPI


        function reportError(~,ME)
            % Reports error from error structure, ME
            fprintf('ERROR: %s\n',ME.message)
            for ii=1:length(ME.stack)
                 fprintf(' on line %d of %s\n', ME.stack(ii).line,  ME.stack(ii).name)
            end
            fprintf('\n')
        end % reportError


        function isGreater = versionGreaterThan(obj,verToTest)
            % Return true if the current ScanImage version is newer than that defined by string verToTest
            %
            % SIBT.versionGreaterThan(obj,verToTest)
            %
            % Inputs
            % verToTest - should be in the format '5.6' or '5.6.1' or
            % '2020.0'
            %
            % Note: this method does not know what to do with the update
            % mumber from SI Basic. So 2020.1 is OK but 2020.1.4 won't
            % produce correct results

            isGreater = nan;
            if ~ischar(verToTest)
                return
            end

            % Add '.0' if needed
            if length(strfind(verToTest,'.'))==0
                verToTest = [verToTest,'.0'];
            end

            % Turn string into a number
            verToTestAsNum = str2num(strrep(verToTest,'.',''));

            % Current version
            curVersion = [obj.hSI.VERSION_MAJOR,obj.hSI.VERSION_MINOR];
            if ischar(curVersion(1))
                % Likely this a free release
                curVersionAsNum = str2num(strrep(curVersion,'.',''));
            else
                % Likely this is Basic or Premium
                curVersionAsNum = curVersion(1)*10 + curVersion(2);
            end

            isGreater = curVersionAsNum>verToTestAsNum;
        end % versionGreaterThan


        function scannerType = scannerType(obj)
            % Since SI 5.6, scanner type "resonant" is returned as "rg"
            % This method returns either "resonant" or "linear"
            scannerType = lower(obj.hSI.hScan2D.scannerType);
            if contains(scannerType,'rg') || strcmp('resonant',scannerType)
                scannerType = 'resonant';
            elseif strcmp('gg',scannerType)
                scannerType='linear';
            else
                fprintf('Unknown scanner type %s\n', scannerType)
            end
        end % scannerType


        function acquireAndWait(obj,block)
            % Start a Grab acquisition and block until SI completes it.

            if nargin<2
                block=true;
            end

            obj.hSI.startGrab % Acquire

            if ~block
                return
            end
            while 1
                if strcmp(obj.hSI.acqState,'idle') %Break when finished
                    break
                end
                pause(0.5)
            end % while
        end % acquireAndWait


        function setZSlices(obj,nSlices)
            % Set the number of slices to acquire in a z-stack
            % Handles differences across versions of SI.
            if obj.versionGreaterThan('2020')
                obj.hSI.hStackManager.numSlices=nSlices;
                obj.hSI.hStackManager.numVolumes = 1;
            else
                obj.hSI.hStackManager.numSlices=nSlices;
                obj.hSI.hFastZ.numVolumes=nSlices;
            end
        end % setZSlices


        function chanName = getSaveChannelName(obj)
            % Return the name of the channel being saved as a string
            %
            % Purpose
            % We want to log to the file name the channel name being saved.
            % If more than one channel has been selected for saving we will
            % return empty and prompt the user to select only one channel
            % to save.
            %
            % Outputs
            % chanName - string defining the name of the channel to save.
            %       If more than one channel is being saved it returns empty.

            if length(obj.hSI.hChannels.channelSave) > 1
                fprintf('Select just one channel to save\n')
                chanName = [];
                return
            end
            chanName = obj.hSI.hChannels.channelName{obj.hSI.hChannels.channelSave};
            chanName = strrep(chanName,' ', '_');
        end % getSaveChannelName


        function turnOffAllPMTs(obj)
            % Turn off all PMTs
            obj.hSI.hPmts.powersOn = obj.hSI.hPmts.powersOn*0;
        end % turnOffAllPMTs


        function turnOnAllPMTs(obj)
            % Turn on all PMTs
            obj.hSI.hPmts.powersOn = ones(1,length(obj.hSI.hPmts.powersOn));
        end % turnOffAllPMTs


        function saveAllChannels(obj)
            % Enable all available channels for saving
            obj.hSI.hChannels.channelSave = 1:obj.numberOfAvailableChannels;
        end % saveAllChannels


        function numChans = numberOfAvailableChannels(obj)
            % Return the number of available channels as an integer
            % These are all the channels that the microscope system can possibly acquire.
            % They may not all have a connected PMT.
            % TODO -- I think there is no way to know whether one is connected
            numChans = obj.hSI.hChannels.channelsAvailable;
        end % numberOfAvailableChannels


        function numPMTs = numberOfAvailablePMTs(obj)
            % Return the number of PMTs with an connected DAQ line as an integer
            numPMTs =  cellfun(@(x) ~isempty(x.hAOGain), obj.hSI.hPmts.hPMTs)
        end % numberOfAvailableChannels


        function setPMTgains(obj,gain)
            % Set gains of all PMTs
            %
            % Inputs
            % gain - If a scalar, the same gain is applied to all PMTs. If a vector
            %      with the same length as the number of PMTs, we set each PMT gain
            %      to the value it corresponds to in the vector.

            if isempty(gain)
                return
            end

            if length(gain)==1
                obj.hSI.hPmts.gains = repmat(gain,1,4);
            elseif length(obj.hSI.hPmts.gains) == length(gain)
                obj.hSI.hPmts.gains = gain(:)';
            end
        end % turnOffAllPMTs


        function zFactStr = returnZoomFactorAsString(obj)
            % Return the zoom factor as a neatly formatted string for file names.
            %
            % Inputs
            % none
            %
            % Outputs
            % Returns a string specifying the current ScanImage zoom factor. The string
            % is used for building file names so the '.' is replaced with '-'

            zFactStr = strrep(num2str(obj.hSI.hRoiManager.scanZoomFactor),'.','-');
        end % returnZoomFactorAsString


        function pointBeam(obj)
            obj.hSI.scanPointBeam
        end


        function parkBeam(obj)
            % Park the beam (abort scanning)
            %

            obj.hSI.abort
        end % parkBeam


        function setLaserPower(obj,powerFraction,beamIndex)
            % Set the laser power as a fraction
            %
            % Inputs
            % powerFraction - the power fraction to which the laser power should be set
            % beamIndex - 1 by default. Defines the beam to change.
            %
            % Outputs
            % none

            if nargin<3
                beamIndex = 1;
            end

            obj.hSI.hBeams.hBeams{beamIndex}.setPowerFraction(powerFraction)
        end % setLaserPower


        function powerIn_mW = powerPercent2Watt(obj,powerFraction,beamIndex)
            % Convert a laser power fraction value to mW
            %
            % Inputs
            % powerFraction - a laser power fraction (0 to 1)
            % beamIndex - 1 by default. Defines the beam to change.
            %
            % Outputs
            % powerIn_mW - the expected number of mW at this power fraction

            if nargin<3
                beamIndex = 1;
            end

            if powerFraction<0 || powerFraction>1
                powerIn_mW = [];
                return
            end

            powerIn_mW = obj.hSI.hBeams.hBeams{beamIndex}.convertPowerFraction2PowerWatt(powerFraction);
        end % powerPercent2Watt


        function setBeamMinMaxPowerInW(obj,minMaxW,beamIndex)
            % Set the min and max power of the beam 
            %
            % Purpose
            % The Machine Data File of SI via the Settings panel determines the maximum and
            % minimum power. The user can change this by editing the settings box manually.
            % This method does the same thing for a defined beam.
            %
            % Inputs
            % minMaxW - Vector of length 2 that defines [minW,maxW] of the laser
            % beamIndex - 1 by default. Defines the beam to change.
            %
            % Outputs
            % none


            if nargin<3
                beamIndex = 1;
            end

            if isempty(minMaxW) || ~isvector(minMaxW) || length(minMaxW)~=2
                return
            end

            obj.hSI.hBeams.hBeams{beamIndex}.powerFraction2PowerWattLut(:,2) = minMaxW;

        end % setBeamMinMaxPowerInW


        function numBeams = numberOfAvailableBeams(obj)
            % Return the number of available beams as an integer
            numBeams = numel(obj.hSI.hBeams.hBeams);
        end % numberOfAvailableBeams

        function disableChannelOffsetSubtraction(obj)
            % Disable the offset subtraction for the PMT inputs
            %
            % Purpose
            % For some data we would like to calculate the offset ourselves rather than
            % relying on ScanImage's offset subtraction.
            %
            % Inputs
            % none
            %
            % Outputs
            % None

            obj.hSI.hChannels.channelSubtractOffset(:)=0;
        end % disableChannelOffsetSubtraction


    end % Close methods


end % Close classdef
