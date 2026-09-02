function varargout = electrical_noise(data_dir,varargin)
% Longitudinal electrical noise plots showing FWHM and max values over time
%
% mpqc.longitudinal.electrical_noise(maintenace_folder_path, varargin)
%
% Optional inputs: 'startDate', 'year-month-day'
% mpqc.longitudinal.electrical_noise(maintenace_folder_path, 'startDate', '2024-06-20')
% Plots all data from given day forward
%
% Purpose
% Plots of the pixel value at two standard deviations of electrical noise for each channel
% with PMTs off. If there is significant change in electrical noise, the
% pixel value will increase.
%
%
% Outputs
% out (optional) - structure containing key information and data.
%
%
% Isabell Whiteley, SWC AMF 2025

if nargin<1
    data_dir = pwd;
end

inputOptions = parseLongitudinalInputVariable(varargin{:});

debugPlots = false;

maintenanceFiles = dir(fullfile(data_dir,'**','*.tif'));
n=1;

for ii=1:length(maintenanceFiles)
    tmp = maintenanceFiles(ii);

    if contains(tmp.name,'electrical_noise')
        plotting_template(n).full_path_to_data = fullfile(tmp.folder,tmp.name);
        plotting_template(n).type = 'electrical_noise';
        plotting_template(n).plotting_func = @mpqc.plot.electrical_noise;
        plotting_template(n).date = datetime(regexp(tmp.name, '(\d{4}-\d{2}-\d{2})_(\d{2}-\d{2}-\d{2})','match'),'InputFormat','yyyy-MM-dd_HH-mm-ss');
        [pathstr,plotting_template(n).name,ext] = fileparts(tmp.name);
        n=n+1;
    end
end
if ~exist('plotting_template','var')
    disp('No electrical noise files found')
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

% Each file is reduced to two scalars per channel: the maximum pixel value and twice
% the standard deviation. The image stacks themselves are not needed beyond that, so a
% file is processed then discarded rather than all of them being held in memory.
chanSave = cell(1,length(plotting_template)); % Saved channels of each file
allChan = []; % Union of the saved channels over all files
maxVal = [];
im_2SD = [];

for q = 1:length(plotting_template) % each date

    if ~contains(plotting_template(q).full_path_to_data, '.tif')
        continue
    end

    [noiseData,metaData] = mpqc.tools.scanImage_stackLoad(plotting_template(q).full_path_to_data);
    noiseData = single(noiseData);
    chanSave{q} = metaData.channelSave;
    allChan = union(allChan,chanSave{q});

    if isempty(maxVal)
        % Results are indexed by hardware channel number rather than by position within
        % the saved channels. scanImage_stackLoad returns the third dimension in
        % channelSave order, so a system saving channels 2 and 4 would otherwise report
        % them as 1 and 2. Channels that were never saved stay NaN. The arrays cannot be
        % sized until the first file has been read, since the channel count comes from
        % its metadata.
        channelName = metaData.channelName;
        maxVal = nan(numel(channelName),length(plotting_template));
        im_2SD = nan(numel(channelName),length(plotting_template));
    end

    if debugPlots
        fig = mpqc.tools.returnFigureHandleForFile(sprintf('%s_%02d',mfilename,q));
    end

    for t = 1:length(chanSave{q}) % each PMT

        % Extract data
        tChan = chanSave{q}(t);
        t_im = noiseData(:,:,t);

        maxVal(tChan,q) = max(t_im(:));
        im_2SD(tChan,q) = std(t_im(:))*2;

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
            sgtitle(string(plotting_template(q).date))
            title(['PMT # ',num2str(t)])
            hold off
        end
    end
end

fig = mpqc.tools.returnFigureHandleForFile(sprintf('%s_%02d',mfilename,q));
xlabels = string([plotting_template.date]);

hold on
for ii = 1:length(allChan)
    tChan = allChan(ii);
    plot(im_2SD(tChan,:),  'DisplayName', channelName{tChan})
end

hold off
xticks(1:length(xlabels))
xticklabels(xlabels)
title('Two SD')
ylabel('Pixel value')
legend('Location','best')


% Output of the main function
if nargout>0
    out.fileName = {plotting_template(:).name};
    out.twoSD = im_2SD;
    out.maxValues = maxVal;
    out.date ={plotting_template(:).date};
    out.channelSave = allChan; % Rows of twoSD/maxValues are indexed by hardware channel
    out.channelName = channelName;
    varargout{1} = out;
end

end % close main funtion
