function generateDashboardData(data_dir)
% Creates machine readable data for html dashboard
%
% mpqc.report.generateDashboardData(maintenance_folder_path)
%
% Purpose
% Outputs files for longitudinal tracking of power, electrical noise, standard light
% source data, and photon counting from lens paper data. These files will be used to
% generate files for an HTML dashboard
%
%
% Isabell Whiteley, SWC AMF 2026


% Note what figures are open by user
figsBefore = findall(groot,'Type','figure');

dashboardData = struct();

% System information, read from the most recent settings file in the data directory
dashboardData.system = systemInfoFromSettings(data_dir);

% Data section
dashboardData.metrics = struct();

% Electrical noise
disp('Searching for electrical noise data')
en = mpqc.longitudinal.electrical_noise(data_dir);
if isempty(en)
    disp('No electrical noise data')
else
    % Rows of twoSD are indexed by hardware channel, so unsaved channels are NaN.
    % Report only the channels this system actually acquires, named by their hardware
    % channel number rather than by their position in the array.
    chan = en.channelSave;
    varNames = ['date', cellstr("PMT" + chan(:)')];
    twoSD = en.twoSD(chan,:)';
    % The longitudinal functions return dates as datetime, so no parsing is needed.
    % Going via the default display format would tie this to the machine's locale.
    dates = [en.date{:}]';
    datesISO = cellstr(string(dates, "yyyy-MM-dd'T'HH:mm:ss"));
    n = numel(datesISO);
    data = cell(n, 1);
    for ii = 1:n
        data{ii} = [{datesISO{ii}}, num2cell(twoSD(ii, :))];
    end
    dashboardData.metrics.electricalNoise.label = 'Electrical Noise';
    dashboardData.metrics.electricalNoise.units = 'pixel value';
    dashboardData.metrics.electricalNoise.variable_names = varNames;
    dashboardData.metrics.electricalNoise.data = data;
end

% Power
disp('Searching for power data')
pow = mpqc.longitudinal.power(data_dir);
if isempty(pow)
    disp('No power data')
else
    % TO DO add in wavelength
    varNames = {'date','maxPower_mW','percentAt100mW'};
    dates = [pow.date{:}]';
    datesISO = cellstr(string(dates, "yyyy-MM-dd'T'HH:mm:ss"));
    n = numel(datesISO);
    power = cell(n, 1);
    for ii = 1:n
        power{ii} = {datesISO{ii}, pow.maxPower(ii), pow.percentAt100mW(ii)};
    end
    dashboardData.metrics.laserPower.label = 'Laser Power';
    dashboardData.metrics.laserPower.units = 'mW';
    dashboardData.metrics.laserPower.variable_names = varNames;
    dashboardData.metrics.laserPower.data = power;
end

% Standard light source
disp('Searching for standard light source data')
stdLight = mpqc.longitudinal.standard_light_source(data_dir);
if isempty(stdLight)
    disp('No standard light source data')
else
    % As for electrical noise, meanValue rows are indexed by hardware channel
    chan = stdLight.channelSave;
    varNames = ['date', cellstr("ch" + chan(:)')];
    dates = [stdLight.maxDate{:}]';
    datesISO = cellstr(string(dates, "yyyy-MM-dd'T'HH:mm:ss"));
    n = numel(datesISO);
    stdLightData = cell(n, 1);
    for ii = 1:n
        stdLightData{ii} = [{datesISO{ii}}, num2cell(stdLight.meanValue(chan, ii)')];
    end
    dashboardData.metrics.standardLight.label = 'Standard Light Source';
    dashboardData.metrics.standardLight.units = 'mean pixel value';
    dashboardData.metrics.standardLight.variable_names = varNames;
    dashboardData.metrics.standardLight.data = stdLightData;
end

% Lens paper photons per pixel
disp('Searching for photon counting data')
photons = mpqc.longitudinal.lens_paper(data_dir,'skipStandardSource',true);
if isempty(photons)
    disp('No photon counting data')
else
    % lens_paper returns cell arrays when the data are split into power groups but bare
    % numeric arrays when there is only one group. Normalise to cells so that both cases
    % produce one named series per power range.
    powerRanges = photons.powerRanges;
    perPixel = photons.photonsPerPixel;
    if ~iscell(powerRanges)
        powerRanges = {powerRanges};
        perPixel = {perPixel};
    end

    % The dashboard needs series names as strings, not as [min,max] numeric pairs
    varNames = ['date', cellfun(@powerRangeName, powerRanges, 'UniformOutput', false)];
    dates = [photons.date{:}]';
    datesISO = cellstr(string(dates, "yyyy-MM-dd'T'HH:mm:ss"));
    n = numel(datesISO);
    photonsPerPixel = cell(n, 1);
    for ii = 1:n
        values = cellfun(@(x) x(ii), perPixel);
        photonsPerPixel{ii} = [{datesISO{ii}}, num2cell(values)];
    end
    dashboardData.metrics.photonsPerPixel.label = 'Photons Per Pixel';
    dashboardData.metrics.photonsPerPixel.units = 'photons/pixel';
    dashboardData.metrics.photonsPerPixel.variable_names = varNames;
    dashboardData.metrics.photonsPerPixel.data = photonsPerPixel;
end

% Uniform slide
disp('Searching for uniform slide data')
xSections = mpqc.longitudinal.uniform_slide(data_dir);
if isempty(xSections)
    disp('No cross section data')
else
    varNames = {'date','X-axis','Y-axis', 'FOV size'};
    dates = [xSections.date{:}]';
    datesISO = cellstr(string(dates, "yyyy-MM-dd'T'HH:mm:ss"));
    n = numel(datesISO);
    xSections.FOVsize = abs(xSections.FOVsize(1))+abs(xSections.FOVsize(end));
    FOV = cell(n, 1);
    for ii = 1:n
        FOV{ii} = {datesISO{ii}, xSections.profile_x(ii,:), xSections.profile_y(ii,:),xSections.FOVsize(1)};
    end
    dashboardData.metrics.uniformSlide.label = 'Field Uniformity';
    dashboardData.metrics.uniformSlide.units = 'um';
    dashboardData.metrics.uniformSlide.variable_names = varNames;
    dashboardData.metrics.uniformSlide.data = FOV;
end

outfile = fullfile(data_dir, 'longitudinalDashboardData.json');
fid = fopen(outfile, 'w');
if fid == -1
    error('mpqc:report:generateDashboardData:FileOpenFailed', ...
        'Could not open %s for writing.', outfile);
end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, jsonencode(dashboardData,PrettyPrint=true));
clear cleanup

% closes any figures opened by code
figsAfter = findall(groot,'Type','figure');
newFigs = setdiff(figsAfter, figsBefore);
delete(newFigs)

end % generateDashboardData



function name = powerRangeName(range)
% Turn a [min,max] power range in mW into a series name for the dashboard
range = double(range(:))';
if isscalar(range) || range(1) == range(end)
    name = sprintf('power_%gmW', range(1));
else
    name = sprintf('power_%gto%gmW', min(range), max(range));
end
end % powerRangeName



function sys = systemInfoFromSettings(data_dir)
% Build the "system" section of the dashboard data from the settings YAML
%
% The settings file lives in each acquisition sub-directory. We use the most recent
% one, since that describes the current state of the microscope. If no settings file
% can be found or read we return an empty structure rather than failing: the metrics
% are still worth writing out.

sys = struct();

d = dir(fullfile(data_dir, '**', '*_SystemSettings.yml'));
if isempty(d)
    fprintf('No settings file found in %s. Dashboard "system" section will be empty.\n', ...
        data_dir)
    return
end
[~, ind] = max([d.datenum]);
settingsFile = fullfile(d(ind).folder, d(ind).name);

try
    y = mpqc.yaml.ReadYaml(settingsFile);
catch ME
    fprintf('Could not read %s: %s\n', settingsFile, ME.message)
    return
end

if isfield(y, 'microscope')
    sys.microscope = nullIfEmpty(getFieldOrEmpty(y.microscope, 'name'));
    sys.roomNumber = nullIfEmpty(getFieldOrEmpty(y.microscope, 'roomNumber'));
end

if isfield(y, 'objective')
    sys.objective = struct( ...
        'name', nullIfEmpty(getFieldOrEmpty(y.objective, 'name')), ...
        'serialNumber', nullIfEmpty(getFieldOrEmpty(y.objective, 'serialNumber')));
end

% All four PMTs are reported, including those that are not fitted, so that the
% channel numbering in the dashboard lines up with the hardware.
PMTs = {};
for ii = 1:4
    tPMT = sprintf('PMT_%d', ii);
    if ~isfield(y, tPMT)
        continue
    end
    PMTs{end+1} = struct( ...
        'id', tPMT, ...
        'model', nullIfEmpty(getFieldOrEmpty(y.(tPMT), 'model')), ...
        'channelName', nullIfEmpty(getFieldOrEmpty(y.(tPMT), 'microscopeChannelName')), ...
        'bandPassFilter', nullIfEmpty(getFieldOrEmpty(y.(tPMT), 'bandPassFilter')));
end
sys.PMTs = PMTs;

% Only the first configured laser is reported
for ii = 1:3
    tLaser = sprintf('imagingLaser_%d', ii);
    if isfield(y, tLaser) && ~isempty(getFieldOrEmpty(y.(tLaser), 'model'))
        sys.imagingLaser = struct( ...
            'model', y.(tLaser).model, ...
            'serialNumber', nullIfEmpty(getFieldOrEmpty(y.(tLaser), 'serialNumber')));
        break
    end
end

if isfield(y, 'QC')
    sourceIDs = getFieldOrEmpty(y.QC, 'sourceIDs');
    if isempty(sourceIDs)
        sourceIDs = {};
    elseif ~iscell(sourceIDs)
        sourceIDs = {sourceIDs};
    end
    sys.QC = struct('sourceIDs', {sourceIDs});
end

end % systemInfoFromSettings



function val = getFieldOrEmpty(s, fieldName)
if isstruct(s) && isfield(s, fieldName)
    val = s.(fieldName);
else
    val = [];
end
end % getFieldOrEmpty



function val = nullIfEmpty(val)
% jsonencode writes NaN as null but writes [] as an empty array. An unset setting is
% a missing scalar value, so null is the correct representation for it.
if isempty(val)
    val = NaN;
end
end % nullIfEmpty
