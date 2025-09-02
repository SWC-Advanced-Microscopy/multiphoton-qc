function [d,out,noiseData,fwhm] = long_electrical_noise(data_dir)
% Longitudinal electrical noise plots
%
% mpqc.longitudinal.long_electrical_noise
%
% Purpose
% Plots of the maximum value and FWHM of electrical noise for each channel 
% with PMTs off. If there is significant change in electrical noise, the FWHM 
% distributions will likely increase in value.
%
%
% Isabell Whiteley, SWC AMF

if nargin<1
    data_dir = pwd;
end

d = dir(fullfile(data_dir,'\**\*.tif'));
n=1;
for ii=1:length(d)
    tmp = d(ii);

    if contains(tmp.name,'electrical_noise')
        out(n) = generic_generator_template(tmp);
        out(n).type = 'electrical_noise';
        out(n).plotting_func = @mpqc.plot.electrical_noise;
        out(n).date = string(tmp.date);
        [pathstr,out(n).name,ext] = fileparts(tmp.name);
        n=n+1;
        % else
        %     disp('No files')
    end
end

    function out = generic_generator_template(t_dir)
        out.full_path_to_data = fullfile(t_dir.folder,t_dir.name);
        out.type = [];
        out.plotting_func = [];
        out.name = [];
        out.date = [];
    end

for i = 1:length(out)
    if contains(out(i).full_path_to_data, '.tif')
        % noiseData(:,:,i) = imread(out(i).full_path_to_data);
        noiseData(:,:,:,i) = mpqc.tools.scanImage_stackLoad(out(i).full_path_to_data);
    end
end

noiseData = single(noiseData);
fwhm = 1; % TEMPORARY to allow script to run
for q = 1:size(noiseData,4) % each date
    figure;
    for t = 1:size(noiseData,3) % each PMT
        %
        
        subplot(2,2,t)
        t_im = noiseData(:,:,t,q);
        [n,x] = hist(t_im(:),100); % plots all data as histogram
        
        a=area(n);
        a.EdgeColor=[0,0,0.75];
        a.FaceColor=[0.5,0.5,1];
        a.LineWidth=2;
        hold on
        m = smoothdata(n,'gaussian',5);
        detail = interp1(x,m,[1:1000]);
        b = plot(m);
        b.LineWidth = 2;
        sgtitle(out(q).date)
        title(['PMT # ',num2str(t)])

        hold off
        % maxVal(t,q) = max(detail(:));
        % halfMaxVal = maxVal(t,q)/2;
        % leftIndex = find(detail(:) >= halfMaxVal, 1, 'first');
        % rightIndex = find(detail(:) >= halfMaxVal, 1, 'last');
        % fwhm(t,q) = rightIndex -leftIndex;
    end
end

% xlabels = {out.date};
% figure;
% subplot(2,1,1)
% plot(maxVal)
% xticks(1:length(xlabels))
% xticklabels(xlabels)
% title('Max value')
% 
% subplot(2,1,2)
% plot(fwhm)
% xticks(1:length(xlabels))
% xticklabels(xlabels)
% title('FWHM')
end