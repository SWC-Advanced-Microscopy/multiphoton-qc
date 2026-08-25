function varargout = uniform_slide(data_dir,varargin)
% Longitudinal uniform slide plots showing FWHM in x and y
%
% mpqc.longitudinal.uniform_slide(maintenace_folder_path, varargin)
%
% Optional inputs: 'startDate', 'year-month-day'
% mpqc.longitudinal.uniform_slide(maintenance_folder_path, 'startDate', '2024-06-20')
% Plots all data from given day forward
%
% Purpose
%
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

    % Data smoothed and normalised using same method as plot.renderer.uniform_slide

     % Smooth it a bit before plotting. Contours and cross sections are further
    % smoothed on top of this (see below).
    % The imresize along rows removes artifacts caused by amplifier ringing
    meanSlide = imresize(meanSlide,[round(size(meanSlide,1)*0.75), size(meanSlide,2)]);
    meanSlide = imresize(meanSlide,size(slideData,[1 2]));
    fSize=round(size(meanSlide,1)/10);
    meanSlide = medfilt2(meanSlide,[fSize,fSize],'symmetric'); %filter heavily

    % data normalised using same method as plot.renderer.uniform_slide
    meanSlide = meanSlide-min(meanSlide(:));
    meanSlide = meanSlide/max(meanSlide(:));


    cx = round(nx/2);
    cy = round(ny/2);

    % Central x and y profiles
    profile_x(ii,:) = meanSlide(cx,:);
    profile_y(ii,:) = meanSlide(:,cy)';

end


% fig = mpqc.tools.returnFigureHandleForFile(sprintf('%s_%02d',mfilename,q));
subplot(1,2,1)
hold on 
for ii = 1:size(profile_x,1)
    plot(profile_x(ii,:))
end
hold off
legend(string([plotting_template.date]),'location','best') 
% xlim([profile_x(1,1),profile_x(1,end)])
ylim([0,1])

subplot(1,2,2)
hold on 
for ii = 1:size(profile_y,1)
    plot(profile_y(ii,:))
end
hold off
legend(string([plotting_template.date]),'location','best')
% xlim([profile_y(1),profile_y(end)])
ylim([0,1])

end