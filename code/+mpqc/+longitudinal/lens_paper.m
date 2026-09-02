function varargout = lens_paper(data_dir,varargin)
% Longitudinal lens paper plots showing mean photons per pixel over time
%
% mpqc.longitudinal.lens_paper(maintenace_folder_path, varargin)
%
% Optional inputs: 'startDate', 'year-month-day'
% mpqc.longitudinal.lens_paper(maintenace_folder_path, 'startDate', '2024-06-20')
% Plots all data from given day forward
% 'skipStandardSource', true
% Excludes standard light source data when calculating photons per pixel
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
skipStandardSource = inputOptions.skipStandardSource;

maintenanceFiles = dir(fullfile(data_dir,'**','*.tif'));
n=1;

for ii=1:length(maintenanceFiles)
    tmp = maintenanceFiles(ii);

    if contains(tmp.name,'lens_paper')
        plotting_template(n).full_path_to_data = fullfile(tmp.folder,tmp.name);
        plotting_template(n).type = 'lens_paper';
        plotting_template(n).plotting_func = @mpqc.plot.lens_paper;
        plotting_template(n).wavelength = str2num(cell2mat(regexp(tmp.name,'\d*(?=nm)','match')));
        plotting_template(n).power = str2num(cell2mat(regexp(tmp.name,'\d*(?=mW)','match')));
        plotting_template(n).date = datetime(regexp(tmp.name, '(\d{4}-\d{2}-\d{2})_(\d{2}-\d{2}-\d{2})','match'),'InputFormat','yyyy-MM-dd_HH-mm-ss');
        [pathstr,plotting_template(n).name,ext] = fileparts(tmp.name);
        n=n+1;
    end
end

if ~exist('plotting_template','var')
    disp('No lens paper files found')
    varargout{1} = [];
    return
end

% sort plotting_template data by the date/time
date_list = [plotting_template.date];
[~,order] = sort(date_list,'ascend');
plotting_template = plotting_template(order);


if ~isempty(inputOptions.startDate)
    startDate = inputOptions.startDate;
    startIndex = 1;
    while startIndex <= numel(plotting_template) && [plotting_template(startIndex).date] < startDate
        startIndex = startIndex + 1;
    end
    plotting_template = plotting_template(startIndex:end);
end

% If only one wavelength is measured and the power is within 20mW
if isequal(plotting_template(:).wavelength) &&  max([plotting_template.power]) - min([plotting_template.power]) <= 20

    powerRange = cell(1,1);
    legendLabels = cell(1,1);
    photonsPerPixel = zeros(1,length(plotting_template));

    for ii = 1:length(plotting_template)
        if contains(plotting_template(ii).full_path_to_data, '.tif')

            % calculate photons per pixel
            data = mpqc.analyse.get_quantalsize_from_file(plotting_template(ii).full_path_to_data,[],skipStandardSource);
            photonsPerPixel(ii) = data.mean_photons_per_pixel;
            powerRange = [min([plotting_template.power]), max([plotting_template.power])];
            legendLabels = sprintf('Power between  %g-%g mW', ...
                min([plotting_template.power]), max([plotting_template.power]));
        end
    end

    fig = mpqc.tools.returnFigureHandleForFile(['long_',mfilename]);
    plot(photonsPerPixel)
    xlabels = string([plotting_template.date]);
    xticks(1:length(xlabels))
    xticklabels(xlabels)
    title('Photons per Pixel')
    ylabel('Photons')
    legend(legendLabels,'Location','best')

    if nargout>0
        out.fileName = {plotting_template(:).name};
        out.photonsPerPixel = photonsPerPixel;
        out.date ={plotting_template(:).date};
        out.powerRanges = powerRange;
        varargout{1} = out;
    end

else
    % Group measurements by power within 20mW
    [sortedPower, sortIdx] = sort([plotting_template.power]);

    groups = {};
    currentGroup = sortIdx(1);
    groupMin = sortedPower(1);

    for k = 2:numel(sortedPower)
        if sortedPower(k) - groupMin <= 20
            currentGroup(end+1) = sortIdx(k);
        else
            groups{end+1} = currentGroup;
            currentGroup = sortIdx(k);
            groupMin = sortedPower(k);
        end
    end
    groups{end+1} = currentGroup;


    fig = mpqc.tools.returnFigureHandleForFile(['long_',mfilename]);

    legendLabels = cell(1,numel(groups));
    powerRange = cell(1,numel(groups));
    photonsPerPixel = cell(1,numel(groups));

    for g = 1:numel(groups)
        idx = groups{g};
        groupPowers = [plotting_template(idx).power];

        legendLabels{g} = sprintf('Power between  %g-%g mW', ...
            min(groupPowers), max(groupPowers));
        powerRange{g}= [min(groupPowers), max(groupPowers)];
        photonsPerPixel{g} = nan(1,length(plotting_template));
        for ii = 1:length(idx)
            data = mpqc.analyse.get_quantalsize_from_file( ...
                plotting_template(idx(ii)).full_path_to_data,[],skipStandardSource);
            photonsPerPixel{g}(idx(ii)) = data.mean_photons_per_pixel;
        end

        plot(photonsPerPixel{g})
        hold on
    end

    hold off
    xlabels = string([plotting_template.date]);
    xticks(1:length(xlabels))
    xticklabels(xlabels)
    legend(legendLabels,'Location','best')
    title('Photons per Pixel')
    ylabel('Photons')


    % TODO If more than one wavelength has been used

    if nargout>0
        out.fileName = {plotting_template(:).name};
        out.photonsPerPixel = photonsPerPixel;
        out.date ={plotting_template(:).date};
        out.powerRanges = powerRange;
        varargout{1} = out;
    end
end

end
