function varargout = lens_paper(fname,aveBy)
    % Display lens paper images
    %
    % plot.lens_paper(fname)
    %
    % Purpose
    % Display lens paper images. These are used for qualitative comparison only.
    % If the file contains more than one saved channel, one figure is made per
    % channel.
    %
    % Inputs
    % fname - relative or absolute path to tiff containing the data
    % aveBy - 1 by default. If >1 average by this much to simulate a slower scanner.
    %
    %
    % Outputs
    % params - optionally return key imaging parameters as a structure per channel
    % txt - optionally return legend text for report as a cell array, one per
    %       channel
    % figHandles - optionally return the figure handles, one per channel
    %
    % Rob Campbell, SWC AMF, initial commit 2022



    % bring up a file picker UI if the user did not input a file
    if nargin<1 || isempty(fname)
        [t_file, t_path] = uigetfile('*.tif', 'Select a file');
        if isequal(t_file, 0)
            return
        else
            fname = fullfile(t_path, t_file);
        end
    end

    if nargin<2
        aveBy = 1;
    end

    [fullstack,metadata] = mpqc.tools.scanImage_stackLoad(fname);
    if isempty(fullstack)
        return
    end

    micsPerPixelXY = metadata.micsPerPixelXY;
    nChansSaved = length(metadata.channelSave);

    % Trick to slice fullstack with multiple PMTs
    h = sibridge.readTifHeader(fname);
    nPMTs = length(h.gains);
    savedPMTs = 1:min(nChansSaved,nPMTs);

    figHandles = gobjects(1,length(savedPMTs));

    for kk = 1:length(savedPMTs)
        pmt = savedPMTs(kk);
        selectedChan = metadata.channelSave(pmt);

        % select just this channel out of the interleaved stack
        imstack = fullstack(:,:,pmt:nChansSaved:end);

        % try averaging to simulate a slower scanner
        if aveBy>1
            n=floor(size(imstack,3)/aveBy);
            t=ones(size(imstack,1),size(imstack,2),n);
            ind = 1;
            for ii = 1:aveBy:size(imstack,3)-aveBy+1
                t(:,:,ind) = mean(imstack(:,:,ii:ii+aveBy-1),3);
                ind = ind+1;
            end

            imstack=t;
        end



        % Make a new figure or return a plot handle as appropriate
        fig = mpqc.tools.returnFigureHandleForFile([fname,mfilename,sprintf('_ch%d',selectedChan)]);
        figHandles(kk) = fig;

        im_mu = mean(imstack,3);


        subplot(2,2,1)

        imagesc(im_mu)
        axis equal tight
        colormap gray
        cMax = getColorScaleLim(im_mu,0.005);
        caxis([0,cMax])
        colorbar

        mpqc.tools.add_scale_axis_tick_labels(gca,micsPerPixelXY)
        title(sprintf('Mean lens paper image (Channel %d)', selectedChan))


        subplot(2,2,2)
        imagesc(imstack(:,:,1))
        axis equal tight
        colormap gray
        cMax = getColorScaleLim(im_mu,0.001);
        caxis([0,cMax])
        colorbar

        mpqc.tools.add_scale_axis_tick_labels(gca,micsPerPixelXY)
        title(sprintf('Single frame (Channel %d)', selectedChan))



        subplot(2,2,3:4)
        hist(im_mu(:),1000)
        ax = gca;
        set(ax.XAxis,'Scale','Log')
        xlabel('Log mean pixel intensity')
        ylabel('#')


        % Optionally return key parameters as a structure
        if nargout>0
            out(kk).laser_power_in_mw = mpqc.report.laser_power_from_fname(fname);
            out(kk).laser_wavelength_in_nm = mpqc.report.laser_wavelength_from_fname(fname);

            out(kk).PMT_gain_in_V = h.gains(pmt);
            out(kk).input_range = h.channelsInputRanges{pmt};
            out(kk).PMT_name = h.names{pmt};
            out(kk).channel = selectedChan;
            out(kk).nChans = length(savedPMTs);
            out(kk).image_size = size(im_mu);
            out(kk).averagFrames = aveBy;
        end

        if nargout>1
            [~,main_fname,ext] = fileparts(fname);
            txt{kk} = sprintf(['%s (Channel %d)\nLens paper imaged at %d mW at %d nm. ', ...
                'Using %s at %dV. Input range %d/%d V. Acquired at %d x %d at %d FPS. '...
                ], ...
                [main_fname,ext], ...
                out(kk).channel, ...
                out(kk).laser_power_in_mw, ...
                out(kk).laser_wavelength_in_nm, ...
                out(kk).PMT_name, ...
                out(kk).PMT_gain_in_V, ...
                out(kk).input_range, ...
                metadata.pixelsPerLine, ...
                metadata.linesPerFrame, ...
                round(metadata.scanFrameRate) );
        end

    end % kk

    if nargout>0
        varargout{1} = out;
    end

    if nargout>1
        varargout{2} = txt;
    end

    if nargout>2
        varargout{3} = figHandles;
    end



function colorScaleLim = getColorScaleLim(im,clip_prop)
    % return the maximum color value to plot such that we are not clipping
    % "prop" proportion of the values. e.g. prop of about 0.9 should work.

    if nargin<2
        clip_prop = 0.01;
    end

    sortedVals = sort(im(:),'descend');
    f = round(length(sortedVals)*clip_prop);
    colorScaleLim = sortedVals(f);
