function varargout = uniform_slide(data_dir,varargin)
% Longitudinal uniform slide plots showing cross sections of the FOV 
%
% mpqc.longitudinal.uniform_slide(maintenace_folder_path, varargin)
%
% Optional inputs: 'startDate', 'year-month-day'
% mpqc.longitudinal.uniform_slide(maintenance_folder_path, 'startDate', '2024-06-20')
% Plots all data from given day forward
%
% Purpose
% Displays the central cross section in x and y of a uniform FOV. Can be
% used to identify any beam movement/misalignment over time
%
%
% Outputs
% out (optional) - structure containing key information and data.
%
%
% Isabell Whiteley, SWC AMF 2026

if nargin<1
    data_dir = pwd;
end

inputOptions = parseLongitudinalInputVariable(varargin{:});

debugPlots = false;

maintenanceFiles = dir(fullfile(data_dir,'**','*.tif'));
n=1;

for ii=1:length(maintenanceFiles)
    tmp = maintenanceFiles(ii);

    if contains(tmp.name,'uniform_slide')
        plotting_template(n).full_path_to_data = fullfile(tmp.folder,tmp.name);
        plotting_template(n).type = 'uniform_slide';
        plotting_template(n).plotting_func = @mpqc.plot.uniform_slide;
        plotting_template(n).date = datetime(regexp(tmp.name, '(\d{4}-\d{2}-\d{2})_(\d{2}-\d{2}-\d{2})','match'),'InputFormat','yyyy-MM-dd_HH-mm-ss');
        [pathstr,plotting_template(n).name,ext] = fileparts(tmp.name);
        n=n+1;
    end
end
if ~exist('plotting_template','var')
    disp('No uniform slide files found')
    varargout{1} = [];
    return
end

% sort plotting_template data by the date/time
date_list = [plotting_template.date];
[~,order] = sort(date_list,'ascend');
plotting_template = plotting_template(order);

if ~isempty(inputOptions.startDate) % Optional variable for selecting starting date
    startDate = inputOptions.startDate;
    startIndex = 1;

    while startIndex <= numel(plotting_template) && [plotting_template(startIndex).date] < startDate
        startIndex = startIndex + 1;
    end

    plotting_template = plotting_template(startIndex:end);
end

% Longitudinal analysis

profile_x = nan(length(plotting_template),256); % hard coded because record.uniform_slide sets the number of pixels
profile_y = nan(length(plotting_template),256);

for ii = 1:length(plotting_template)
    if ~contains(plotting_template(ii).full_path_to_data, '.tif')
        continue
    end

    [slideData,metaData] = mpqc.tools.scanImage_stackLoad(plotting_template(ii).full_path_to_data);
    meanSlide = mean(slideData,[3]);
    [nx, ny] = size(meanSlide);
    micsPerPixelXY = metaData.micsPerPixelXY;

    % Utilising renderer.uniform_slide to take the cross sections and
    % plotting info
    [profile_x(ii,:),profile_y(ii,:),xData] = mpqc.plot.renderer.uniform_slide(meanSlide,micsPerPixelXY,[],'plotting',false);

    fig = mpqc.tools.returnFigureHandleForFile(['long_',mfilename]);
    subplot(1,2,1)
    hold on
    for ii = 1:size(profile_x,1)
        plot(xData,profile_x(ii,:))
    end
    hold off
    legend(string([plotting_template.date]),'location','best')
    title('X-axis')
    xlim([xData(1),xData(end)])
    ylim([0,1])
    ylabel('normalized intensity')
    xlabel('microns')

    subplot(1,2,2)
    hold on
    for ii = 1:size(profile_y,1)
        plot(xData,profile_y(ii,:))
    end
    hold off
    legend(string([plotting_template.date]),'location','best')
    title('Y-axis')
    xlim([xData(1),xData(end)])
    ylim([0,1])
    ylabel('normalized intensity')
    xlabel('microns')

    % Output of the main function
    if nargout>0
        out.fileName = {plotting_template(:).name};
        out.profile_y = profile_y;
        out.profile_x = profile_x;
        out.date ={plotting_template(:).date};
        out.FOVsize = xData;
        varargout{1} = out;
    end

end