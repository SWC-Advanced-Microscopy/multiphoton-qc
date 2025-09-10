function varargout = long_power(data_dir,varargin)
% Function to track the changes in laser power over time
%
%
%
if nargin<1
    data_dir = pwd;
end


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



% Output of the main function
if nargout>0
    out.fileName = {plotting_template(:).name};
    out.date ={plotting_template(:).date};
    varargout{1} = out;
end

end

function out = generic_generator_template(t_dir)
out.full_path_to_data = fullfile(t_dir.folder,t_dir.name);
out.type = [];
out.plotting_func = [];
out.name = [];
out.date = [];
end