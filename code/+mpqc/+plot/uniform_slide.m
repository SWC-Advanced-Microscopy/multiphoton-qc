function uniform_slide(fname,varargin)
    % Make nice plots of the uniform slide data to explore illumination
    %
    % mpqc.plot.uniform_slide(fname,'param1','val1','param2','val2',...)
    %
    % Purpose
    % Make a plot of field homogeneity based on a uniform fluorescent slide or solution.
    % Assumes data obtained from ScanImage (ideally at Zoom 1).
    %
    % Inputs [required]
    % fname - relative or absolute path to image stack
    %
    % Inputs [optional]
    % overlayZoom - Vector indicating which zoom values to overlay as boxes.
    %               A reasonable selection chosen by default.
    % crossSections - Which directions the image cross-sections should run. 'diagonal' or 'scanner'.
    %              If 'scanner', the lines run through the centre parallel with the scan axes.
    %              This is the default. If 'diagonal' they run from the image corners, which are the
    %              darkest parts of the field of view.
    %
    % Examples
    % Plot with no zoom boxes overlaid
    % mpqc.plot.uniform_slide('uniform_slide_zoom_1_920nm_5mW__2022-08-02_10-09-33_00001.tif','overlayZoom',[])
    %
    %
    % Rob Campbell, SWC AMF


    %Parse optional arguments
    params = inputParser;
    params.CaseSensitive = false;
    params.addParamValue('overlayZoom', [1.2,2,4], @(x) isnumeric(x) && isscalar(x) || isvector(x) || isempty(x));
    params.addParamValue('crossSections', 'scanner', @(x) isstr(x) && (strcmpi(x,'scanner') || strcmpi(x,'diagonal')));
    params.parse(varargin{:});

    overlayZoom = params.Results.overlayZoom;
    crossSections = params.Results.crossSections;

    [imstack,metadata] = mpqc.tools.scanImage_stackLoad(fname);
    if isempty(imstack)
        return
    end

    micsPerPixelXY = metadata.micsPerPixelXY;


    % Make a new figure or return a plot handle as appropriate
    fig = mpqc.tools.returnFigureHandleForFile([fname,mfilename]);

    subplot(1,2,1)
    plotData = mean(imstack,3);

    % Smooth it a bit before plotting. Contours and cross sections are further
    % smoothed on top of this (see below).
    % The imresize along rows removes artifacts caused by amplifier ringing
    plotData = imresize(plotData,[round(size(plotData,1)*0.75), size(plotData,2)]);
    plotData = imresize(plotData,size(imstack,[1,2]));
    plotData = medfilt2(plotData,[7,7]); %filter heavily

    imagesc(plotData)
    axis equal tight
    colormap gray

    mpqc.tools.add_scale_axis_tick_labels(gca,micsPerPixelXY)

    hold on
    nContours = 10;
    colThreshold = 0.6;
    contourf(plotData, nContours, 'LineColor', 'none');
    if median(plotData(:)/max(plotData(:))) > colThreshold
        caxis([0,colThreshold*max(plotData(:))])
    else
        caxis([0, max(plotData(:))])
    end
    color_tmp = gray(nContours+2);
    colormap(color_tmp(1:end-1,:)) % use grayscale wuthout white


    % Add diagonal lines which we will use later to associate with the next plot
    switch crossSections
        case 'diagonal'
            plot([1,size(plotData,1)], [1,size(plotData,2)], '-r', 'linewidth',2)
            plot([1,size(plotData,1)], [size(plotData,2),1], '-c', 'linewidth',2)
        case 'scanner'
            plot([1,size(plotData,1)], [size(plotData,2)/2,size(plotData,2)/2], '-r', 'linewidth',2) %Y
            plot([size(plotData,2)/2,size(plotData,2)/2], [1,size(plotData,1)], '-c', 'linewidth',2) %X
    end


    % Diagnostic/extra information feature: show boxes indicating the FOV
    % at different zooms.
    if ~isempty(overlayZoom)

        zoom_cols = parula(length(overlayZoom));

        for ii=1:length(overlayZoom)
            L=length(imstack);
            newSize = L/overlayZoom(ii);
            offset = (L-newSize)/2;
            r(ii) = rectangle('Position', [offset,offset,newSize,newSize], ...
                                'EdgeColor', zoom_cols(ii,:));
            text_zoom(ii) = text(offset+1,offset+4, ...
                            sprintf('Zoom %0.1f', overlayZoom(ii)), ...
                            'Color', zoom_cols(ii,:), ...
                            'FontSize',12, ...
                            'FontWeight','Bold');
        end

    end

    hold off
    set(gca,'FontSize',12)


    % Plot intensity cross-sections along the red/cyan lines
    subplot(1,2,2)

    normPlotData = plotData/max(plotData(:));
    switch crossSections
        case 'diagonal'
            micsPerDataPoint = sqrt(2*micsPerPixelXY^2);

            f_diag = eye(length(plotData));
            yData = normPlotData(find(f_diag));
            xData = (1:length(yData)) * micsPerDataPoint;
            xData = xData - mean(xData);

            hXsection1 = plot(xData,yData,'-r','linewidth',2);

            hold on

            yData = normPlotData(find(rot90(f_diag)));

            hXsection2 = plot(xData,yData,'-c','linewidth',2);

        case 'scanner'
            micsPerDataPoint = micsPerPixelXY;

            xSectionX = normPlotData(:, round(size(normPlotData,2)/2));
            xSectionY = normPlotData(round(size(normPlotData,2)/2),:);

            xData = (1:length(xSectionY)) * micsPerDataPoint;
            xData = xData - mean(xData);
            hXsection1 = plot(xData,xSectionY,'-r','linewidth',2);
            hold on
            hXsection2 = plot(xData,xSectionX,'-c','linewidth',2);
    end


    xlim([xData(1),xData(end)])
    ylim([0,1])

    % Add tick labels
    xticks = round(linspace(xData(1)+0.5,xData(end)-0.5,5));
    set(gca, 'Xtick', xticks)
    grid on

    xlabel('microns')
    ylabel('normalized intensity')

    % Overlay lines corresponding with the zooms
    if ~isempty(overlayZoom)
        imSizeInMicrons = length(xData)*micsPerDataPoint;
        L = length(hXsection1.XData);
        for ii=1:length(overlayZoom)
            thisZoom = imSizeInMicrons / overlayZoom(ii);
            newSize = L/overlayZoom(ii);
            offset = (L-newSize)/2;
            plot([thisZoom/2,thisZoom/2], ylim, 'Color', zoom_cols(ii,:))
            plot(-[thisZoom/2,thisZoom/2], ylim, 'Color', zoom_cols(ii,:))
        end
    end

    % Tweak plot properties
    set(gca,'Color',[1,1,1]*0.7, ...
        'FontSize',12)

    set(gcf,'InvertHardcopy','off', 'Color','w')
    % Nicely scale the plot window so the two figure are sized well with respect to each other
    fig.Position(3) = fig.Position(4)*2.3;
