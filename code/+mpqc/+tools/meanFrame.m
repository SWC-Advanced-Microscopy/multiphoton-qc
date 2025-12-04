classdef meanFrame < handle
    % Monitor open ScanImage channel windows and plot a rolling average trace.
    %
    % To use this function simply run
    % "mpqc.tools.meanFrame" then press "Focus" in ScanImage
    % The current frame mean is present in the figure title along with a rolling average.
    %
    %
    % Advanced usage:
    % M=mpqc.tools.meanFrame
    %
    % % Change channel colors
    % M.channelColors='bkrr';
    % M.channelColors='rgkb';
    %
    % % Change the minimum update interval
    % M.minUpdateInterval=2; %To once every two seconds
    % M.minUpdateInterval=0.1; %To every 100 ms (assuming your FPS is fast enough)
    %
    % % Change the number of seconds to display
    % M.secondsToDisplay=5
    % M.secondsToDisplay=3
    %
    %
    % Rob Campbell, SWC AMF


    properties

        hSI % Refence to ScanImage object

        listener_acqState
        listener_frameDone

        hFig
    end

    properties (SetObservable)
        % minUpdateInterval - the minimum time in seconds between plot updates
        % channelColors - a character array of length 4 that defines the channel colors
        % secondsToDisplay - how many seconds worth of data to display on the screen
        minUpdateInterval = 0.05 %The minimum time to wait in seconds between updates of the plot
        channelColors = 'rbgk'
        secondsToDisplay = 30

    end

    properties (Hidden)
        maxChannels = 4
        figureTag = 'meanTraceFig' %to help find the figure window
        axTag = 'meanTraceAx'
        plotDataH
        numDataPoints = 50 % Just a default. This will be modified automatically at run time
    end

    methods
        %constructor
        function obj = meanFrame
            % function mpfstools.meanFrame

            API = sibridge.getSIobject;
            if isempty(API)
                fprintf('ScanImage not found, unable to link to it.\n')
                return
            end

            obj.hSI=API; % Make composite object

            %obj.listener_acqState = addlistener(obj.hSI,'acqState', 'PostSet',  @obj.startStopAcqWithScanImage);

            obj.listener_frameDone = addlistener(obj.hSI.hUserFunctions,'frameAcquired', @obj.updatePlot);



            obj.makeFigWindow


        end %constructor

        function delete(obj)
            delete(obj.listener_acqState)
            delete(obj.listener_frameDone)
            delete(obj.hFig)
        end


        function startStopAcqWithScanImage(obj,~,~)
            % If ScanImage is connected and it starts imaging then
            % acquisition starts. If a file is being saved in ScanImage
            % then this causes this class to save data to disk
            if isempty(obj.hSI)
                fprintf('No ScanImage connected to meanFrame\n')
                return
            end

            switch obj.hSI.acqState
                case {'grab','loop'}
                    % Do nothing because we are in data acquisition mode
                    fprintf('Stopping listener\n')
                    obj.listener_frameDone.Enable=0;
                case 'focus'
                    % Ensure we have a figure window open
                    fprintf('Starting listener\n')
                    obj.listener_frameDone.Enable=1;
                case 'idle'
                    fprintf('Stopping listener\n')
                    obj.listener_frameDone.Enable=0;
            end
        end % startStopAcqWithScanImage


        function updatePlot(obj,~,~)
            framePeriod = obj.hSI.hRoiManager.scanFramePeriod; %frame period in seconds
            obj.numDataPoints=round(obj.secondsToDisplay/framePeriod); %max number of data points


            %Set up the figure for the first time if it does not currently exist
            %obj.makeFigWindow % This really is overly careful...

            % scanimage stores image data in a data structure called 'stripeData'
            lastStripe = obj.hSI.hDisplay.stripeDataBuffer{obj.hSI.hDisplay.stripeDataBufferPointer}; % get the pointer to the last acquired stripeData
            channels =  lastStripe.channelNumbers;


            %Extract means and store in figure's UserData
            timeSinceLastUpdate = (now-obj.hFig.UserData.lastUpdate)*24*60^2;
            meansToDisplay=nan(1,length(obj.channelColors));
            for ii = 1:length(obj.channelColors)
                if ~any(channels==ii) %then it's not being recorded
                    mu=nan;
                else
                    f=find(channels==ii);
                    imData=lastStripe.roiData{1}.imageData{f}{1};
                    mu=mean(imData(:));
                    meansToDisplay(ii)=mu;
                end

                if obj.hFig.UserData.pointsPlotted<obj.numDataPoints
                    obj.hFig.UserData.plotData{ii}(obj.hFig.UserData.pointsPlotted)=mu;
                else %start scrolling once all data points have been filled
                    obj.hFig.UserData.plotData{ii}(end+1)=mu;
                    obj.hFig.UserData.plotData{ii}(1)=[];
                    xlim([1,obj.numDataPoints])
                end

            end

            %Restrict the rate of plot update.
            if (now-obj.hFig.UserData.lastUpdate)*24*60^2 > obj.minUpdateInterval
                titleStr='';
                for ii=1:4
                    obj.plotDataH(ii).YData = obj.hFig.UserData.plotData{ii};
                    obj.plotDataH(ii).Color = obj.channelColors(ii);
                    if ~isnan(meansToDisplay(ii))
                        titleStr = [titleStr,sprintf('CH%d=%0.2f ',ii,meansToDisplay(ii))];
                    end
                end
                obj.hFig.UserData.lastUpdate=now;
                title(titleStr)
                drawnow
            end


            if obj.hFig.UserData.pointsPlotted<obj.numDataPoints
                obj.hFig.UserData.pointsPlotted = obj.hFig.UserData.pointsPlotted+1;
            end


        end %updatePlot


        function makeFigWindow(obj)
            % Build the figure window if it is not already there

            axH = findobj('Tag',obj.axTag);

            if isempty(axH)
                obj.hFig = figure('Tag',obj.figureTag);
                figure(obj.hFig);
                plotData=nan(obj.numDataPoints,obj.maxChannels);
                obj.plotDataH = plot(plotData);
                axH = gca;
                axH.Tag = obj.axTag;

                data.pointsPlotted=1; %Store the number of points plotted to date in the userdata
                data.lastUpdate=now; %Time of last plot update
                for ii=1:length(obj.channelColors)
                    data.plotData{ii}=plotData(:,ii);
                end
                obj.hFig.UserData=data;
                grid on
                set(obj.hFig,'Name','Mean traces')
                xlabel('Frame #')
                ylabel('mean pixel value')
                obj.hFig.CloseRequestFcn = @obj.windowCloseFcn; %So closing the window triggers the destructor
            else
                obj.hFig = findobj('Tag',obj.figureTag);
                obj.plotDataH = axH.Children;
            end
        end % makeFigWindow


        function settings = returnSettings(obj)
            settings.channelColors = obj.channelColors;
            settings.minUpdateInterval = obj.minUpdateInterval;
        end %Returns the settings that will be fed to the mean plotting user functions







        %Getters and setters for the settings
        function set.channelColors(obj,colors)
            if ~ischar(colors)
                fprintf('colors must be a character array\n')
                return
            end
            if length(colors)~=4
                fprintf('colors must be a character array of length 4\n')
                return
            end
            obj.channelColors=colors;
        end

        function set.secondsToDisplay(obj,secondsToDisplay)
            if ~isnumeric(secondsToDisplay)
                fprintf('secondsToDisplay must be a number\n')
                return
            end
            if ~isscalar(secondsToDisplay)
                fprintf('secondsToDisplay must be a scalar\n')
                return
            end

            obj.listener_frameDone.Enabled=0;
            obj.secondsToDisplay=secondsToDisplay;
            obj.wipegraph
            obj.listener_frameDone.Enabled=1;
        end

        function set.minUpdateInterval(obj,minUpdateInterval)
            if ~isnumeric(minUpdateInterval)
                fprintf('minUpdateInterval must be a number\n')
                return
            end
            if ~isscalar(minUpdateInterval)
                fprintf('minUpdateInterval must be a scalar\n')
                return
            end
            obj.minUpdateInterval=minUpdateInterval;
        end


    end %methods



    methods (Hidden)
        function windowCloseFcn(obj,~,~)
            % This runs when the user closes the figure window.
            obj.delete % simply call the destructor
        end %close windowCloseFcn

        function wipegraph(obj)
            % Used when changing the number of displayed points
            framePeriod = obj.hSI.hRoiManager.scanFramePeriod; %frame period in seconds
            obj.numDataPoints=round(obj.secondsToDisplay/framePeriod); %max number of data points

            plotData=nan(obj.numDataPoints,obj.maxChannels);
            obj.plotDataH = plot(plotData);

            data.pointsPlotted=1; %Store the number of points plotted to date in the userdata
            data.lastUpdate=now; %Time of last plot update
            for ii=1:length(obj.channelColors)
                data.plotData{ii}=plotData(:,ii);
            end
            obj.hFig.UserData=data;

            xlim([1,obj.numDataPoints])

            for ii=1:4
                obj.plotDataH(ii).YData = obj.hFig.UserData.plotData{ii};
                obj.plotDataH(ii).Color = obj.channelColors(ii);
            end
            obj.hFig.UserData.lastUpdate=now;
            title('')
        end

    end %hidden methods


end % classdef meanFrame
