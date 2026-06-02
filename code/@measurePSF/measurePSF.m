classdef measurePSF < handle
    % Display PSF and measure its size in X, Y, and Z
    %
    % measurePSF(PSFstack,varargin)
    %
    % USAGE
    % Fit and display a PSF. Reports FWHM to on-screen figure with simple GUI elements
    % for interacting with it. The brightest bead is automatically found and used to estimate
    % the PSF, however the function does not zoom in to the bead. To do this, click on
    % "Select bead", draw a box around the desired bead in the bottom/left plot, then
    % double-click to select. The other plots will update. Zoom back out again by clicking
    % "Reset view".
    %
    %
    % To run a demo: measurePSF('demo')
    %
    %
    % INPUTS (optional)
    % PSFstack  - Either: A 3-D array (imagestack). First layer should be that nearest the objective
    %                 OR: Path to a file containing a TIFF stack with one channel.
    %                 OR: If empty or no input arguments then a file load GUI is presented.
    %
    % INPUTS (optional param/val pairs)
    % micsPixZ - number of microns per pixel in Z (i.e. distance between adjacent Z planes)
    % micsPixXY - number of microns per pixel in X and Y
    % useZmax - [false by default] if true we use the max intensity projection
    %                   for the Z PSFs. This is likely necessary if the PSF is very tilted.
    % zFitOrder - [1 by default]. Number of Gaussians to use for the fit of the Z PSF
    % medFiltSize - [1 by default -- no filtering]. If more than one performs a median filtering
    %               operation on each slice with a filter of this size.
    %
    %
    % OUTPUTS
    % Press the "PSF data to workspace" button to return a structure containing useful
    % data to the base workspace. The "PSFstats" structure contains:
    %
    %
    %          fitStats: [1x1 struct]
    %          PSFstack: [256x256x50 double] % This is the image stack used for the plots in the GUI
    %    micsPerPixelXY: 0.0100              % Microns per pixel supplied by user
    %     micsPerPixelZ: 0.0500              % Microns per pixel supplied by user
    %    figureSnapshot: [800x800x3 uint8]   % Snapshot of figure window which you may view with "imshow"
    %
    % PSFstats.fitStats contains:
    %
    %                 X: [1x1 struct]
    %                 Y: [1x1 struct]
    %                ZX: [1x1 struct]
    %                ZY: [1x1 struct]
    %            FWHMxy: 0.3100 % This of the average of the two fits along the rows and columns
    %             FWHMz: 3.4500 % This is the average of the ZX and ZY fits
    %
    % Where the first four fields have detailed fit statistics for each of axes where
    % a FWHM was estimated. e.g. PSFstats.fitStats.X contains:
    %
    %     fit: [1x1 cfit]
    %    data: [1x1 struct]
    %
    %
    %
    % EXAMPLES:
    % One: bring up file-load GUI and extract voxel size automatically if this is a ScanImage TIFF
    % >> measurePSF;
    %
    % Two: Feed in a matrix and define the voxel size
    % >> T=mpqc.tools.load3Dtiff('PSF_2019-59-15_16-11-55_00001.tif');
    % >> measurePSF(T, 'micsPixZ',0.5, 'micsPixXY',0.1);
    %
    % Three: load a specific file from disk at the command line and also manually specify Z voxel size
    % >> measurePSF('PSF_2019-59-15_16-11-55_00001.tif', 'micsPixZ',0.5);
    %
    % Four: demo mode
    % >> measurePSF('demo');
    %
    %
    % Rob Campbell, Basel Biozentrum, initial commit 2016
    %
    % Requires:
    % Curve-Fitting Toolbox, Image Processing Toolbox


    properties (SetAccess=protected)
        PSFstats  %All stats relating to the PSF are stored here. These can be exported to the base workspace using a button
    end


    properties (SetObservable)
        PSFstack %The image stack containing the PSF


        useMaxIntensityForZpsf
        zFitOrder
        medFiltSize % size of the median filter used to clean up the image stack


        % The following are default values. If they aren't changed
        % (along with the "reportFWHM" properties being true) then
        % the associated FWHM value is not reported.
        micsPerPixelXY=1 % Number of microns per pixel in X/Y
        micsPerPixelZ=5  % Number of microns per pixel in Z
    end


    properties (SetObservable, SetAccess=protected, Hidden)
        % This is the z-plane that contains the brightest region of the bead--we will use this to fit the x/y PSF
        maxZplane
        maxZplaneForFit %These are the data used for fitting. They will always have at least a moderate median filter applied

        %The following are fit-related parameters
        psfCenterInX
        psfCenterInY
        psfCenterInZ
        badFit

        zoomedArea % Area to zoom into if the user has drawn a rectangle over the bottom/left image.
                   % Should be [bottom_pos, left_pos, width, height]
                   % See "areaSelector" and "resetView" callbacks.
    end

    properties (Hidden, SetAccess=protected)
        PSFstack_Orig %The original (un-zoomed) stack
    end

    properties (Hidden)
        hFig % Handle containing the figure
        hPSF_XYmidpointImageAx % The PSF at the estimated mid-point (bottom left image)
        hPSF_XYmidpointImageIM % Handle to the image object
        hPS_midPointImageXhairs % The dashed cross-hairs
        hPSF_midPointText

        hxSectionRowsAx
        hxSectionColsAx

        hPSF_ZXax
        hPSF_ZX_currentZplane % Plot handle for current z-plane line
        hPSF_ZX_fitAx

        hPSF_ZYax
        hPSF_ZY_currentZplane % Plot handle for current z-plane line
        hPSF_ZY_fitAx

        hUserSelectedPlaneAx % Handle to the axis showing the user-selected plane within the PSF (top right)
        hUserSelectedPlaneIM
        hUserSelectedPlaneTitle


        hSlider % Slider handle

        drawBox_PushButton
        reset_PushButton
        fitToBaseWorkSpace_PushButton
        saveImage_PushButton

        useMaxIntensityForZpsf_checkBox
        zFitOrder_editBox
        medFiltSize_editBox
        textMedian
        textZfit


        showHelpTextIfTooFewArgsProvided=false
        listeners={}

        % If user does not supply a pixel size then the associated FWHM value will not be reported to screen
        reportFWHMxy=false
        reportFWHMz=false

        verbose = false % If true we report to screen debugging information of various sorts

        fname % the name of the tiff stack loaded by the user
        dataFolder % The folder that contains the PSF data
    end




    methods
        function obj=measurePSF(inputPSFstack,varargin)

            % Check if required toolboxes are installed
            toolboxesPresent = true;
            toolboxesToTest = {'Image Processing Toolbox','Curve Fitting Toolbox'};
            for ii = 1:length(toolboxesToTest)
                if ~contains(struct2array(ver),toolboxesToTest{ii})
                    fprintf('** measurePSF requires %s **\n', toolboxesToTest{ii})
                    toolboxesPresent = false;
                end
            end

            % Bail out if required toolbox not present (message to user supplied above)
            if toolboxesPresent == false
                return
            end


            % If no input arguments are provided, we bring up the load GUI
            if nargin==0
                [fname, pname]=uigetfile({'*.tif';'*.tiff'},'Select PSF stack','MultiSelect','off');
                if ~ischar(fname) && fname==0
                    return
                end
                inputPSFstack = fullfile(pname,fname);
            end

            if nargin==1 && ischar(inputPSFstack) && strcmp(inputPSFstack,'demo')
                if obj.showHelpTextIfTooFewArgsProvided
                    help(mfilename)
                end
                inputPSFstack=[];
                demoMode=true;
                % We will display the default PSF and use the default values for XY and Z pixel size
                % as defined in the properties section
                % The default PSF is loaded at the end of the constructor
            else
                demoMode=false;
            end

            % Parse optional param/val pairs

            params = inputParser;
            params.CaseSensitive = false;
            params.addParamValue('useZmax', 1, @(x) islogical(x) || x==0 || x==1);
            params.addParamValue('zFitOrder', 1, @(x) isnumeric(x) && isscalar(x));
            params.addParamValue('medFiltSize', 1, @(x) isnumeric(x) && isscalar(x));
            params.addParamValue('micsPixZ', [], @(x) isnumeric(x) && isscalar(x));
            params.addParamValue('micsPixXY',[], @(x) isnumeric(x) && isscalar(x));

            params.parse(varargin{:});

            obj.useMaxIntensityForZpsf = params.Results.useZmax;
            obj.zFitOrder = params.Results.zFitOrder;
            obj.medFiltSize = params.Results.medFiltSize;
            micsPerPixelZ = params.Results.micsPixZ;
            micsPerPixelXY = params.Results.micsPixXY;


            if demoMode
                micsPerPixelZ=0.5;
                micsPerPixelXY=0.1;
            end

            % Load PSF stack if it was provided as a file
            if ischar(inputPSFstack)

                obj.fname=inputPSFstack;

                [inputPSFstack,metadata] = mpqc.tools.scanImage_stackLoad(obj.fname);


                % Allow for user-specified values to over-ride what the header returns
                if ~exist('micsPerPixelZ','var') || isempty(micsPerPixelZ)
                    micsPerPixelZ = metadata.stackZStepSize;
                end
                if ~exist('micsPerPixelXY','var') || isempty(micsPerPixelXY)
                    micsPerPixelXY = metadata.micsPerPixelXY;
                end

                % Get the folder in which the file is situated
                obj.dataFolder = fileparts(which(obj.fname));
            else
                % The "dataFolder" is otherwise the current directory
                obj.dataFolder = pwd;
            end % ischar(inputPSFstack)


            if exist('micsPerPixelZ','var') && isnumeric(micsPerPixelZ) && isscalar(micsPerPixelZ)
                obj.micsPerPixelZ = micsPerPixelZ;
                obj.reportFWHMz=true;
            else
                fprintf('No valid value found for FWHM in Z. Not reporting it.\n')
            end

            if exist('micsPerPixelXY','var') && isnumeric(micsPerPixelXY) && isscalar(micsPerPixelXY)
                obj.micsPerPixelXY = micsPerPixelXY;
                obj.reportFWHMxy=true;
            else
                fprintf('No valid value found for FWHM in XY. Not reporting it.\n')
            end

            % Make empty data and generate empty plots using these
            obj.PSFstack = zeros(2^8);
            obj.maxZplane = zeros(2^8);
            obj.maxZplaneForFit = zeros(2^8);
            obj.psfCenterInX=2^7;
            obj.psfCenterInZ=2^7;
            obj.psfCenterInZ=1;


            obj.setUpFigureWindow; % Make the axes and so forth

            % Set up listeners that will update the plots when the PSF stack is modified or
            % other relevant properties are changed.
            obj.listeners{end+1} = addlistener(obj, 'PSFstack', 'PostSet', @obj.plotNewImageStack);
            obj.listeners{end+1} = addlistener(obj, 'maxZplaneForFit', 'PostSet', @obj.fitPSFandUpdateSlicePlots);
            obj.listeners{end+1} = addlistener(obj, 'useMaxIntensityForZpsf', 'PostSet', @obj.redrawGUI);
            obj.listeners{end+1} = addlistener(obj, 'medFiltSize', 'PostSet', @obj.redrawGUI);
            obj.listeners{end+1} = addlistener(obj, 'zFitOrder', 'PostSet', @obj.redrawGUI);

            % If no PSF stack was provided, we loads the default
            if demoMode
                P = load('PSF');
                obj.fname = 'DEMO_SIMULATED_PSF';
                obj.addNewStack(P.PSF)
            else
                obj.addNewStack(inputPSFstack);
            end

        end %Close constructor

        function addNewStack(obj,newStack)
            newStack = double(newStack);
            newStack = newStack - min(newStack(:));  %needed in case the amps are offset from zero
            % Replace the current image stack with a new one (or add a stack on startup)
            obj.zoomedArea = [1,1,size(newStack,1)-1,size(newStack,2)-1];
            obj.PSFstack_Orig = newStack;
            obj.PSFstack = newStack;
            obj.updateUserSelected; %Ensure the white lines showing the current user z-plane are correct
        end


        function delete(obj)
            cellfun(@delete,obj.listeners)
            if ~isempty(obj.hFig)
                obj.hFig.delete
            end
        end % Close destructor


        function windowCloseFcn(obj,~,~)
            % This runs when the user closes the figure window
            obj.delete % simply call the destructor
        end % Close windowCloseFcn



        function reportMethodEntry(obj)
            % Print to screen the name of the method that was just entered
            if ~obj.verbose
                return
            end
            ST=dbstack;
            fprintf('Entered method %s\n',ST(2).name)

        end % Close reportMethodEntry


        %-----------------------------------------------------------------------------
        % Callback functions follow
        function plotNewImageStack(obj,~,~)
            obj.reportMethodEntry
            %This is run when the PSFstack property is changed
            s=size(obj.PSFstack);
            obj.hFig.Name = sprintf('Image size: %d x %d',s(1:2));
            %Clean the stack and find the mid-point and produce a filtered image plane at this point
            obj.denoiseImStackAndFindPSFcenterInZ;

            %Plot the mid-point of the stack
            obj.hPSF_XYmidpointImageIM.CData = obj.maxZplane;

            %Update overlay text
            obj.hPSF_midPointText.Position(1) = size(obj.PSFstack,1)*0.025;
            obj.hPSF_midPointText.Position(2) = size(obj.PSFstack,2)*0.04;
            obj.hPSF_midPointText.String = sprintf('PSF center at slice #%d',obj.psfCenterInZ);


            %Optionally, show the axes. Right now, I don't think we want this at all so it's not an input argument
            showAxesInMainPSFplot=0;
            if showAxesInMainPSFplot
                Xtick = linspace(1,size(obj.maxZplane,1),8);
                Ytick = linspace(1,size(obj.maxZplane,2),8);
                set(obj.PSF_XYmidpointImageAx,'XTick',Xtick,'XTickLabel',round(Xtick*obj.micsPerPixelXY,2),...
                        'YTick',Ytick,'YTickLabel',round(Ytick*obj.micsPerPixelXY,2));
            else
                set(obj.hPSF_XYmidpointImageAx,'XTick',[],'YTick',[], ...
                    'XLim',[0.5,size(obj.PSFstack,1)], ...
                    'YLim',[0.5,size(obj.PSFstack,2)])
            end

            % Place image into the top/right plot and update the slider
            obj.hUserSelectedPlaneIM.CData = obj.maxZplane;
            obj.hUserSelectedPlaneAx.XLim = [0.5,size(obj.PSFstack,1)];
            obj.hUserSelectedPlaneAx.YLim = [0.5,size(obj.PSFstack,2)];
            set(obj.hSlider, 'Max',size(obj.PSFstack,3), 'Value',obj.psfCenterInZ)


            % Modify the lines to show where we are slicing it to take the cross-sections
            obj.setCrossSectionLinesInMainPSFImage

            obj.hUserSelectedPlaneTitle.String = sprintf('Slice #%d', obj.psfCenterInZ);
        end % Close plotNewImageStack


        function setCrossSectionLinesInMainPSFImage(obj,~,~)
            %Set new cross-hair location on the bottom/left image based upon the
            %x/y centroid of the bead (PSF)
            obj.reportMethodEntry
            set(obj.hPS_midPointImageXhairs(1), 'XData', obj.hPSF_XYmidpointImageAx.XLim, ...
                'YData', [obj.psfCenterInY,obj.psfCenterInY])
            set(obj.hPS_midPointImageXhairs(2), 'XData', [obj.psfCenterInX,obj.psfCenterInX], ...
                'YData', obj.hPSF_XYmidpointImageAx.YLim)
        end % Close setCrossSectionLinesInMainPSFImage


        function fitPSFandUpdateSlicePlots(obj,~,~)
            % This callback is run whenever the raw data are updated or
            % whenever properties that might affect the fit are updated.
            obj.reportMethodEntry

            userZdepth = round(get(obj.hSlider,'Value')); %So we can plot the currently shown z-plane

             % Find the X/Y max location and update the properties psfCenterInX and psfCenterInY
            obj.findPSF_centreInXY(obj.maxZplaneForFit);

            OUT=obj.updateXYfits;

            % Obtain images showing the PSF's extent in Z
            % We do this by taking maximum intensity projections or slices through the maximum

            %PSF in Z/X (upper panel)
            if obj.useMaxIntensityForZpsf
                PSF_ZX=squeeze(max(obj.PSFstack,[],1));
            else
                PSF_ZX=squeeze(obj.PSFstack(obj.psfCenterInY,:,:));
            end

            imagesc(obj.hPSF_ZXax, PSF_ZX)


            obj.hPSF_ZXax.NextPlot='Add';
            obj.hPSF_ZX_currentZplane = plot(obj.hPSF_ZXax,[userZdepth, userZdepth], obj.hPSF_ZXax.YLim, ':w');
            obj.hPSF_ZXax.NextPlot='Replace';

            Ytick = linspace(1,size(PSF_ZX,1),3);
            set(obj.hPSF_ZXax,'XAxisLocation','Top',...
                    'CLim', [min(PSF_ZX(:)),max(PSF_ZX(:))], ...
                    'XTick',[],...
                    'YTick',Ytick,'YTickLabel',round(Ytick*obj.micsPerPixelXY,2));

            %This is the fitted Z/Y PSF with the FWHM
            axes(obj.hPSF_ZX_fitAx)
            cla
            maxPSF_ZX = max(PSF_ZX,[],1);
            baseline = sort(maxPSF_ZX);
            baseline = mean(baseline(1:5));
            maxPSF_ZX = maxPSF_ZX-baseline;

            fitZX = obj.fit_Intensity(maxPSF_ZX, obj.micsPerPixelZ);
            x = (1:length(maxPSF_ZX))*obj.micsPerPixelZ;
            [OUT.ZX.FWHM,OUT.ZX.fitPlot_H] = obj.plotCrossSectionAndFit(x, maxPSF_ZX, fitZX, obj.micsPerPixelZ/4, 0, 'XZ');
            set(obj.hPSF_ZX_fitAx,'XAxisLocation','Top')

            %Suppress title with FWHM estimate if no mics per pixel was provided
            if ~obj.reportFWHMz
                title('XZ: pizel size not provided')
            end
            obj.PSFstats.ZX.im = maxPSF_ZX;
            obj.PSFstats.ZX.fit = fitZX;

            % PSF in Z/Y (panel on the right on the right)
            if obj.useMaxIntensityForZpsf
                PSF_ZY=squeeze(max(obj.PSFstack,[],2));
            else
                PSF_ZY=squeeze(obj.PSFstack(:,obj.psfCenterInX,:));
            end

            PSF_ZY=rot90(PSF_ZY,3);
            imagesc(obj.hPSF_ZYax, PSF_ZY)

            obj.hPSF_ZYax.NextPlot='Add';
            obj.hPSF_ZY_currentZplane = plot(obj.hPSF_ZYax, obj.hPSF_ZYax.XLim, [userZdepth, userZdepth], ':w');
            obj.hPSF_ZYax.NextPlot='Replace';

            Xtick = linspace(1,size(PSF_ZY,2),3);
            set(obj.hPSF_ZYax,'YAxisLocation','Right', ...
                    'CLim', [min(PSF_ZY(:)),max(PSF_ZY(:))], ...
                    'XTick',Xtick,'XTickLabel',round(Xtick*obj.micsPerPixelXY,2), ...
                    'YTick',[])

            %This is the fitted Z/X PSF with the FWHM
            axes(obj.hPSF_ZY_fitAx);
            cla
            maxPSF_ZY = max(PSF_ZY,[],2);
            baseline = sort(maxPSF_ZY);
            baseline = mean(baseline(1:5));
            maxPSF_ZY = maxPSF_ZY-baseline;

            fitZY = obj.fit_Intensity(maxPSF_ZY, obj.micsPerPixelZ);
            x = (1:length(maxPSF_ZY))*obj.micsPerPixelZ;
            [OUT.ZY.FWHM, OUT.ZY.fitPlot_H] = obj.plotCrossSectionAndFit(x,maxPSF_ZY,fitZY,obj.micsPerPixelZ/4,1,'YZ');
            set(obj.hPSF_ZY_fitAx,'XAxisLocation','Top')

            %Suppress title with FWHM estimate if no mics per pixel was provided
            if ~obj.reportFWHMz
                title('YZ: pizel size not provided')
            end
            obj.PSFstats.ZY.im = maxPSF_ZY;
            obj.PSFstats.ZY.fit = fitZY;

            %Add PSF FWHM to stats if appropriate
            if obj.reportFWHMxy
                obj.PSFstats.FWHMxy = mean([OUT.XYcols.FWHM,OUT.XYrows.FWHM]);
            else
                obj.PSFstats.FWHMxy = nan;
            end

            if obj.reportFWHMz
                obj.PSFstats.FWHMz = mean([OUT.ZY.FWHM,OUT.ZX.FWHM]);
            else
                obj.PSFstats.FWHMz = nan;
            end

        end % Close fitPSFandUpdateSlicePlots


        function updateUserSelected(obj,~,~)
            % Runs when the user moves the slider
            obj.reportMethodEntry
            thisSlice = round(get(obj.hSlider,'Value'));
            obj.hUserSelectedPlaneIM.CData = obj.PSFstack(:,:,thisSlice);

            caxis([min(obj.PSFstack(:)), max(obj.PSFstack(:))])

            obj.hUserSelectedPlaneTitle.String = sprintf('Slice #%d', thisSlice);

            %Move the dashed lines on the cross-section plots
            obj.hPSF_ZX_currentZplane.XData = [thisSlice,thisSlice];
            obj.hPSF_ZY_currentZplane.YData = [thisSlice,thisSlice];

            obj.hUserSelectedPlaneAx.CLim = obj.hPSF_XYmidpointImageAx.CLim; %Don't change the lookup table

        end %Close updateUserSelected


        function copyFitToBaseWorkSpace(obj,~,~)
            % Copy the PSFstats property to the base workspace.
            obj.reportMethodEntry
            if isempty(obj.PSFstats)
                return
            end

            varName = 'PSFstats';
            fprintf('Copying PSF fit to base work space as variable "%s"\n', varName)

            OUT.fitStats = obj.PSFstats;
            OUT.PSFstack = obj.PSFstack;

            if obj.reportFWHMxy
                OUT.micsPerPixelXY = obj.micsPerPixelXY;
            else
                OUT.micsPerPixelXY = nan;
            end

            if obj.reportFWHMz
                OUT.micsPerPixelZ = obj.micsPerPixelZ;
            else
                OUT.micsPerPixelZ = nan;
            end

            %Take snapshot of figure window
            snap = getframe(obj.hFig);
            OUT.figureSnapshot = snap.cdata;

            assignin('base',varName, OUT)
        end % Close copyFitToBaseWorkSpace


        function areaSelector(obj,~,~)
            obj.reportMethodEntry
            %select a sub-region of the bottom left plots
            h = imrect(obj.hPSF_XYmidpointImageAx);
            rect_pos = wait(h);
            obj.zoomedArea = round([rect_pos(1:2), mean(rect_pos(3:4)), mean(rect_pos(3:4))]);
            delete(h)
            za = obj.zoomedArea;

            obj.PSFstack = obj.PSFstack(za(2):za(2)+za(3), za(1):za(1)+za(4), :);

            obj.updateUserSelected

        end % Close areaSelector


        function resetView(obj,~,~)
            % Callback that resets (zooms back out) the bottom left view
            obj.reportMethodEntry
            % Un-zoom other panels
            resetSize = [1,1,size(obj.PSFstack_Orig,1)-1,size(obj.PSFstack_Orig,2)-1];

            %Only apply if different to avoid hitting listeners
            if ~isequal(obj.zoomedArea,resetSize)
                obj.zoomedArea=resetSize;
                obj.PSFstack = obj.PSFstack_Orig;
                obj.updateUserSelected
            end
        end % Close areaSelector


        function redrawGUI(obj,~,~)
            % Used to apply changes to setting such as median filter size
            % This is called by a listener on the properties themselves.
            obj.PSFstack = obj.PSFstack_Orig;
        end % Close redrawGUI


        function maxIntCallback(obj,~,~)
            obj.useMaxIntensityForZpsf = obj.useMaxIntensityForZpsf_checkBox.Value;
        end % Close maxIntCallback


        function medFiltSizeCallback(obj,~,~)
            % Callback that runs a median filter on the stack
            newVal = str2double(obj.medFiltSize_editBox.String);
            if isnan(newVal) || newVal<=0
                obj.medFiltSize_editBox.String = num2str(obj.medFiltSize);
            else
                obj.medFiltSize=newVal;
            end
        end % Close medFiltSizeCallback


        function zFitOrderCallback(obj,~,~)
            newVal = str2double(obj.zFitOrder_editBox.String);
            if isnan(newVal) || newVal<=0
                obj.zFitOrder_editBox.String = num2str(obj.zFitOrder);
            else
                obj.zFitOrder=newVal;
            end
        end % Close medFiltSizeCallback


        function saveImage(obj,~,~)
            % Callback that saves the figure window to a PDF
            obj.toggleUIelements('off')
            origString = '';
            origFontSize= [];
            try
                % Get the first part of the name
                [a,b]=regexp(obj.fname,'PSF_.*20\d{2}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');
                im_name = obj.fname(a:b);

                % Get the centre of the ROI in order to create a unique file name
                ROI_centre = round(obj.zoomedArea(1:2) + obj.zoomedArea(3:4)/2);
                im_name = sprintf('Bead_%s__ROI_%d_%d__Zfwhm_%d-%d', ...
                                im_name, ...
                                ROI_centre, ...
                                floor(obj.PSFstats.FWHMz), ...
                                round(rem(obj.PSFstats.FWHMz, floor(obj.PSFstats.FWHMz)),1)*10 );

                fname = fullfile(obj.dataFolder,im_name);

                % Show the file name of the tiff stack (if available) on screen.
                origString = obj.hUserSelectedPlaneTitle.String;
                origFontSize = obj.hUserSelectedPlaneTitle.FontSize;
                if ~isempty(obj.fname)
                    obj.hUserSelectedPlaneTitle.String = strrep(obj.fname,'_','\_');
                    obj.hUserSelectedPlaneTitle.FontSize = 9;
                end
            catch ME
                obj.toggleUIelements('on')
                if ~isempty(origString)
                    obj.hUserSelectedPlaneTitle.String = origString;
                end
                if ~isempty(origFontSize)
                    obj.hUserSelectedPlaneTitle.FontSize = origFontSize;
                end
                rethrow(ME)
            end

            %warning('off','MATLAB:handle_graphics:exceptions:SceneNode')
            print('-dpdf','-bestfit',[fname,'.pdf']) % PDF

            obj.hFig.UserData = obj.PSFstats;
            saveas(obj.hFig, [fname,'.fig']) % Figure file to incorporate into a report
            obj.hFig.UserData = [];

            fprintf('Saved image to: %s\n',fname)
            obj.toggleUIelements('on')
            obj.hUserSelectedPlaneTitle.String = origString;
            obj.hUserSelectedPlaneTitle.FontSize = origFontSize;
        end % Close saveImage

    end % close methods


    methods (Hidden)

        function toggleUIelements(obj,toggleState)
            % Toggle UI elemnents to allow for prettier saved images
            % toggleState should be the string 'on' or 'off'
            obj.drawBox_PushButton.Visible = toggleState;
            obj.reset_PushButton.Visible = toggleState;
            obj.fitToBaseWorkSpace_PushButton.Visible = toggleState;
            obj.saveImage_PushButton.Visible = toggleState;
            obj.useMaxIntensityForZpsf_checkBox.Visible=toggleState;
            obj.zFitOrder_editBox.Visible=toggleState;
            obj.medFiltSize_editBox.Visible=toggleState;
            obj.textZfit.Visible=toggleState;
            obj.textMedian.Visible=toggleState;
        end %toggleUIelements

    end %Hidden methods

end % close class
