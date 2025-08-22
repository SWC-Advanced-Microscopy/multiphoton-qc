function [d,out,noiseData,fwhm] = long_electrical_noise(data_dir)



if nargin<1
    data_dir = pwd;
end

d = dir(fullfile(data_dir,'\**'));
n=1;
for ii=1:length(d)
    tmp = d(ii);

    if contains(tmp.name,'electrical_noise')
        out(n) = generic_generator_template(tmp);
        out(n).type = 'electrical_noise';
        out(n).plotting_func = @mpqc.plot.electrical_noise;
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
end
for i = 1:length(out)
    if contains(out(i).full_path_to_data, '.tif')
        noiseData(:,:,i) = imread(out(i).full_path_to_data);
    end
end
noiseData = single(noiseData);
% fwhm = 1;
for t = 1:size(noiseData,3)
% 
    t_im = noiseData(:,:,t);
    % t_im = noiseData(:,:,2);
%     % figure;
    [n,x] = hist(t_im(:),100); % plots all data as histogram
    figure;
    a=area(x,n);
    a.EdgeColor=[0,0,0.75];
    a.FaceColor=[0.5,0.5,1];
    a.LineWidth=2;
% 
    % halfMaxVal(t) = max(n(:))/2;
    halfMaxVal = max(n(:))/2;
    leftIndex = find(n(:) >= halfMaxVal, 1, 'first'); 
    rightIndex = find(n(:) >= halfMaxVal, 1, 'last');
    fwhm(t) = n(rightIndex) - n(leftIndex);
    % fwhm = n(rightIndex) - n(leftIndex);
end
figure; plot(fwhm)
end