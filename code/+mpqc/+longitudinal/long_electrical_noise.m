function [d,out,noiseData,fwhm] = long_electrical_noise(data_dir)



if nargin<1
    data_dir = pwd;
end

d = dir(fullfile(data_dir,'\**\*.tif'));
n=1;
for ii=1:length(d)
    tmp = d(ii);

    if contains(tmp.name,'electrical_noise')% && contains(tmp.name)
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
    [n,x] = hist(t_im(:),100); % plots all data as histogram
    figure;
    a=area(n);
    a.EdgeColor=[0,0,0.75];
    a.FaceColor=[0.5,0.5,1];
    a.LineWidth=2;
    hold on
    m = smoothdata(n,'gaussian',5);
    detail = interp1(x,m,[1:1000]);
    b = plot(m);
    b.LineWidth = 2;

    hold off
    maxVal(t) = max(detail(:));
    halfMaxVal = maxVal(t)/2;
    leftIndex = find(detail(:) >= halfMaxVal, 1, 'first'); 
    rightIndex = find(detail(:) >= halfMaxVal, 1, 'last');
    fwhm(t) = rightIndex -leftIndex;
end
figure; 
subplot(2,1,1)
plot(maxVal)
subplot(2,1,2)
plot(fwhm)
end