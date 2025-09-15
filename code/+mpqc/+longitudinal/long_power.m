function varargout = long_power(data_dir,varargin)
% Function to track the changes in laser power over time
%
% mpqc.longitudinal.long_power(maintenace_folder_path, varargin)
%
% Optional inputs: 
% 'startDate', month-year
% 'wavelength, value
%
% Purpose
% Plots the power at the objective from 0-100% and compares maximum output
% power over time. Used to monitor the health of a laser
%
%                
% Outputs
% out (optional) - structure containing key information and data.
%
%
%
% Isabell Whiteley, SWC AMF 2025

if nargin<1
    data_dir = pwd;
end

% For elseif statement of wavelength given in varargin
 % % process inputs 
 %    % addpath("+tools");
 %    optInputs =  parseInputVariable(varargin{:});
 % 
 %    % Extract critical input variable
 %    selectWavelength = optInputs.wavelength;
 %   % TODO Add to parseInputVariable startDate

debugPlots = true;

maintenanceFiles = dir(fullfile(data_dir,'\**\*.mat')); 
n=1;

for ii=1:length(maintenanceFiles)
    tmp = maintenanceFiles(ii);

    if contains(tmp.name,'power')
        plotting_template(n) = generic_generator_template(tmp);
        plotting_template(n).type = 'power';
        plotting_template(n).plotting_func = @mpqc.plot.power;
        plotting_template(n).date = string(datetime(regexp(tmp.name, '(\d{4}-\d{2}-\d{2})_(\d{2}-\d{2}-\d{2})','match'),'InputFormat','yyyy-MM-dd_HH-mm-ss'));
        plotting_template(n).wavelength = str2num(cell2mat(regexp(tmp.name,'\d*(?=nm)','match')));
        [pathstr,plotting_template(n).name,ext] = fileparts(tmp.name);
        n=n+1;
    end
end
if ~exist('plotting_template','var')
    disp('No power measurement files found')
    varargout{1} = [];
    return
end

% sort plotting_template data by the date/time
date_list = [plotting_template.date];
[~,order] = sort(datenum(date_list,'dd-mm-yyyy hh:MM:ss'),1,'ascend');
plotting_template = plotting_template(order);

if nargin > 1 % Optional variable for selecting starting date
    startDate = datetime(varargin{1});
    startIndex = 1;

    while [plotting_template(startIndex).date] < startDate
        startIndex = startIndex + 1;
    end

    plotting_template = plotting_template(startIndex:end);
end

if isequal(plotting_template(:).wavelength) % if wavelength is the same
    for ii = 1:length(plotting_template)
        if contains(plotting_template(ii).full_path_to_data, '.mat')
            % If power data is found, load it and find max value
            powerData(ii) = load(plotting_template(ii).full_path_to_data);
            maxPower(ii) = powerData(ii).powerMeasurements.observedPower(end);

            % Plot the power curves for each date
            hold on
            subplot(2,1,1)
            plot([0:5:100],powerData(ii).powerMeasurements.observedPower,'.')
            legend(plotting_template.date,'location', 'Northwest')
            title(cell2mat(['Power at ', string(plotting_template(1).wavelength), 'nm']))
            xlabel('Percent power')
            ylabel('Power (mW)')
            hold off
        end
    end
    % Plot the maximum laser output over time
    subplot(2,1,2)
    plot(maxPower, '-*')
    title('Maximum laser power')
    xlabels = {plotting_template.date};
    xticks(1:length(xlabels))
    xticklabels(xlabels)
    ylabel('Maximum power (mW)')

% elseif isequal(plotting_template(:).wavelength,selectWavelength)
% elseif ~isempty(selectWavelength)
%     tempStruc = find(plotting_template(:).wavelength == selectWavelength);
% % if wavelength given in varargin, plot only that wavelength and dates 

else
    % disp('Different wavelengths found')
    % powerData = [];
    % maxPower = [];
    for i = 1:length(plotting_template)
        allWave(i) = plotting_template(i).wavelength;
    end
    wavelengthVals = unique(allWave);
    numWavelength = length(wavelengthVals);
% maxPower = cell(length(plotting_template),length(wavelengthVals)); % Need to be cells because they will have empty values
% powerData = cell(length(plotting_template),length(wavelengthVals));
 for jj = 1:numWavelength
     figure;
     for ii = 1:length(plotting_template)
     
        if contains(plotting_template(ii).full_path_to_data, '.mat') && isequal(plotting_template(ii).wavelength,wavelengthVals(jj))
            % If power data is found, load it and find max value
            powerData(ii,jj) = load(plotting_template(ii).full_path_to_data);
            maxPower(ii,jj) = powerData(ii).powerMeasurements.observedPower(end);


              hold on
            % subplot(2,1,1)
            plot([0:5:100],powerData(ii,jj).powerMeasurements.observedPower,'.')
            legend(plotting_template(ii).date,'location', 'Northwest') % incorrect
            title(cell2mat(['Power at ', string(wavelengthVals(jj)), 'nm']))
            xlabel('Percent power')
            ylabel('Power (mW)')
            hold off
        end
     end
 end
    % for ii = 1:length(numWavelength)       
    % end
    % separate plots for each wavelength
    % tried groupcounts
    % in R2025a numunique function may work

    % plot wavelengths separately, but keep all dates - so assisgn NaN to
    % unused dates. Make temporary structures so values are not permanently
    %  replaced
end



% identify wavelength - done
% plot curves on same plot - done
% separate plots based on wavelength
% max value to see decay of laser - done






% Output of the main function
if nargout>0
    out.fileName = {plotting_template(:).name};
    out.date ={plotting_template(:).date};
    out.powerData = powerData;
    out.maxPower = maxPower;
    out.wavelength = {plotting_template(:).wavelength};
    varargout{1} = out;
end

end

function out = generic_generator_template(t_dir)
out.full_path_to_data = fullfile(t_dir.folder,t_dir.name);
out.type = [];
out.plotting_func = [];
out.name = [];
out.date = [];
out.wavelength = [];
end