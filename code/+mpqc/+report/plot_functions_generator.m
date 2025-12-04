function out = plot_functions_generator(data_dir)
% Look in a directory containing maintenance data and generate a structure
% that can be used to produce the plots.
%
% function out = mpqc.report.plot_functions_generator(data_dir)
%
% Purpose
% This function generates a structure that contains information about which
% data exist in the defined folder and which functions will plot these. The
% structure can then be passed to a downstream plotting functions to actually
% make the plots. This could be simply making and saving the plots as separate
% PDFs or using the MATLAB report generator to make a nicely formatted report.
%
%
% Inputs
% data_dir - optional, if nothing is provided the current directory is used
%
%
% Also see:
% mpqc.report.plotAllBasic
%
%
% Rob Campbell, SWC AMF, initial commit 2022


if nargin<1
    data_dir = pwd;
end


d = dir(data_dir);


n=1;
for ii=1:length(d)
    tmp = d(ii);
    if contains(tmp.name,'.tif')
        if contains(tmp.name,'electrical_noise')
            out(n) = generic_generator_template(tmp);
            out(n).type = 'electrical_noise';
            out(n).plotting_func = @mpqc.plot.electrical_noise;
            n=n+1;
        elseif contains(tmp.name,'uniform_slide_')
            out(n) = generic_generator_template(tmp);
            out(n).type = 'uniform_slide';
            out(n).plotting_func = @mpqc.plot.uniform_slide;
            n=n+1;

            out(n) = generic_generator_template(tmp);
            out(n).type = 'laser_stability';
            out(n).plotting_func = @mpqc.plot.uniform_slide_laser_stability;
            n=n+1;
        elseif contains(tmp.name,'lens_paper_')
            out(n) = generic_generator_template(tmp);
            out(n).type = 'lens_paper';
            out(n).plotting_func = @mpqc.plot.lens_paper;
            n=n+1;
        elseif contains(tmp.name,'standard_source')
            out(n) = generic_generator_template(tmp);
            out(n).type = 'standard_source';
            out(n).plotting_func = @mpqc.plot.standard_light_source;
            n=n+1;
        end
    end
end

% If the user has acquired PSFs but not selected beads and saved the results, we warn them.
PSF_files = dir('PSF_*');
bead_files = dir('Bead_PSF_*.fig');

if ~isempty(PSF_files)>0 && isempty(bead_files)
    fprintf('\nYou have acquired PSF data but not selected any beads and saved the images.\n')
    fprintf('You should do this if you want to incorporate bead PSFs in the summary graphics.\n\n')
end

% Otherwise if we have found bead PSF figure files, we can add them to the list
if ~isempty(bead_files)
    for ii = 1:length(bead_files)
        tmp = bead_files(ii);
        out(n) = generic_generator_template(tmp);
        out(n).type = 'bead_psf';
        out(n).plotting_func = @mpqc.report.loadBeadPSF;
        n=n+1;
    end
end


if ~exist('out','var')
    fprintf('No data found in directory %s\n', data_dir)
    out = [];
end


% Internal functions follow
function out = generic_generator_template(t_dir)
    out.full_path_to_data = fullfile(t_dir.folder,t_dir.name);
    out.type = [];
    out.plotting_func = [];
    out.laser_wavelength = mpqc.report.laser_wavelength_from_fname(t_dir.name); %get laser wavelength
    out.laser_power = mpqc.report.laser_power_from_fname(t_dir.name); %get laser power

