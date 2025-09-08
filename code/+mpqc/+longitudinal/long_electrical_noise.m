function varargout = long_electrical_noise(data_dir,varargin)
    % Longitudinal electrical noise plots
    %
    % mpqc.longitudinal.long_electrical_noise(maintenace_folder_path, varargin)
    % Optional inputs: Starting date- month-year
    % TODO - change to varargin so can give a date range to use
    %
    % Purpose
    % Plots of the maximum value and FWHM of electrical noise for each channel 
    % with PMTs off. If there is significant change in electrical noise, the FWHM 
    % distributions will likely increase in value.
    %
    %
    % Outputs
    % out (optional) - structure containing key information and data.
    %
    % Isabell Whiteley, SWC AMF
    
    if nargin<1
        data_dir = pwd;
    end


    debugPlots = true;
    
    d = dir(fullfile(data_dir,'\**\*.tif')); % cTODO -- change "d" to a more descriptive name
    n=1;
    
        for ii=1:length(d)
            tmp = d(ii);

            if contains(tmp.name,'electrical_noise')
                plotting_template(n) = generic_generator_template(tmp);
                plotting_template(n).type = 'electrical_noise';
                plotting_template(n).plotting_func = @mpqc.plot.electrical_noise;
                plotting_template(n).date = string(datetime(regexp(tmp.name, '(\d{4}-\d{2}-\d{2})_(\d{2}-\d{2}-\d{2})','match'),'InputFormat','yyyy-MM-dd_HH-mm-ss'));
                [pathstr,plotting_template(n).name,ext] = fileparts(tmp.name);
                n=n+1;
                % else
                %     TODO add else for it no files found
            end
        end
if ~exist('plotting_template','var')% plotting_template == 0;
disp('No electrical noise files found')
varargout{1} = [];
return
end

    % sort plotting_template data by the date/time
    date_list = [plotting_template.date];
    [~,order] = sort(datenum(date_list,'dd-mm-yyyy hh:MM:ss'),1,'ascend');
    plotting_template = plotting_template(order);

    if nargin > 1
        % take out dates from the order based on the given input date
        startDate = datetime(varargin{1});
        
        % This will represent the index at which our (sorted) data is >=
        % the given threshold value
        startIndex = 1;

        % While the date at the current index is BEFORE our minimum date
        % move the index along as we haven't found the threshold data point
        % The loop will stop at the point at which our order[startIndex]
        % is GREATER THAN OR EQUAL TO our given minimum date threshold
        while [plotting_template(startIndex).date] < startDate
            % if datetime(order) => datetime(startDate) % is equal to or later than startDate keep
                % order = order(withoutdates);
            startIndex = startIndex + 1;
        end

        % Remove the data points before our minimum date by slicing
        % from the found startIndex to the end of the list
        plotting_template = plotting_template(startIndex:end);
    end
    
    for ii = 1:length(plotting_template)
        if contains(plotting_template(ii).full_path_to_data, '.tif')
            noiseData(:,:,:,ii) = mpqc.tools.scanImage_stackLoad(plotting_template(ii).full_path_to_data);
        end
    end
    
    noiseData = single(noiseData);

    for q = 1:size(noiseData,4) % each date

        if debugPlots
            fig = mpqc.tools.returnFigureHandleForFile(sprintf('%s_%02d',mfilename,q));
        end
        for t = 1:size(noiseData,3) % each PMT

            % Extract data
            t_im = noiseData(:,:,t,q);
            [n,x] = hist(t_im(:),100); % plots all data as histogram
            m = smoothdata(n,'gaussian',5);
            detail = interp1(x,m,[1:1000]);

            maxVal(t,q) = max(detail(:));
            halfMaxVal = maxVal(t,q)/2;
            leftIndex = find(detail(:) >= halfMaxVal, 1, 'first');
            rightIndex = find(detail(:) >= halfMaxVal, 1, 'last');
            fwhm(t,q) = rightIndex -leftIndex;

            % Optionally plot
            if debugPlots
                subplot(2,2,t)        
                a=area(n);
                a.EdgeColor=[0,0,0.75];
                a.FaceColor=[0.5,0.5,1];
                a.LineWidth=2;
                hold on
                b = plot(m);
                b.LineWidth = 2;
                sgtitle(plotting_template(q).date)
                title(['PMT # ',num2str(t)])
                hold off
            end
        end
    end
    
    for ii = 1:size(noiseData,3) % plotting FWHM and max value over time for each PMT
        xlabels = {plotting_template.date};
        fig = mpqc.tools.returnFigureHandleForFile(sprintf('%s_%02d',mfilename,ii));
        subplot(2,1,1)
        plot(maxVal(ii,:))
        xticks(1:length(xlabels))
        xticklabels(xlabels)
        title('Max value')
    
        subplot(2,1,2)
        plot(fwhm(ii,:))
        xticks(1:length(xlabels))
        xticklabels(xlabels)
        title('FWHM')
        sgtitle(['PMT # ',num2str(ii)])
    end

    
    % Output of the main function 
    if nargout>0
        out.fileName = {plotting_template(:).name};
        out.noiseData = noiseData; 
        out.fwhm = fwhm;
        out.maxValues = maxVal;
        out.date ={plotting_template(:).date};
        varargout{1} = out;
    end

end % close main funtion



function out = generic_generator_template(t_dir)
    out.full_path_to_data = fullfile(t_dir.folder,t_dir.name);
    out.type = [];
    out.plotting_func = [];
    out.name = [];
    out.date = [];
end