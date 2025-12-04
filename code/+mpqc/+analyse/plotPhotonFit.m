function plotPhotonFit(STATS)
% Diagnostic plots for the compute_sensitivity function.
%
% function mpqc.analyse.plotPhotonFit(STATS)
%
% Purpose
% The function "compute_sensitivity" takes as input a series of image frames from a
% scanning microscope and returns an array of statistics as output. The key statistic
% is "quantal_size", which defines the number or raw data units per photon. The quantal
% size is the slope of variance as a function of mean over pixels. If the fit is poor, or
% the data have an unusual distribution, then the estimate of the slope can be compromised.
% This function allows the user to check visually if the fit is of reasonable quality.
%
% Inputs (Required)
% STATS - The output of compute_sensitivity
%
%
% Example
% OUT = mpqc.tools_get_quantal_size_from_file(fname);
% mpqc.analyse.plotPhotonFit(OUT)
%
%
%
% Rob Campbell, SWC AMF, initial commit February 2025
%
%
% Acknowledgements
% This function is based on this work:
% https://github.com/datajoint/anscombe-numcodecs by Dimitri Yatsenko


tFname = STATS.filename;

if isempty(tFname)
    % If the file name is empty, then probably we are feeding in data manually at the
    % CLI and this function was not called by get_quantalsize_from_file. So make up any
    % filename so figure generation proceeds.
    tFname = 'photon_fit';
end

[~,fname,ext] = fileparts(tFname);
fname = [fname,ext];

fig = mpqc.tools.returnFigureHandleForFile(fname);
fig.Name = fname;



%%
% The photon fit

ax1=subplot(2,3,1);
% Plot the variance as a function of the mean and overlay the fit line
intensity = STATS.min_intensity:STATS.max_intensity - 1;
plot(intensity,STATS.variance,'.', 'color',[0.35,0.35,1],'markersize',4)

hold on
y = STATS.model(1) + intensity*STATS.model(2);
plot(intensity,y,'-r','LineWidth',2)
hold off

xlabel('Intensity')
ylabel('Variance')
grid on
axis tight
% Add some key statistics to the plot title
title(sprintf('Quantal size: %0.1f. Mean photons per pixel: %0.2f', ...
    STATS.quantal_size, STATS.photons_per_pixel))


ax2=subplot(2,3,4);
copyobj(ax1.Children,ax2);
set(gca,'YScale','log')
set(gca,'XScale','log')
xlabel('Intensity')
ylabel('Variance')
grid on
box on
axis tight
% Add some key statistics to the plot title
title(sprintf('Quantal size: %0.1f. Mean photons per pixel: %0.2f', ...
    STATS.quantal_size, STATS.photons_per_pixel))


%%
% The converted image
subplot(2,3,[2,5])
[imStack,metadata]=mpqc.tools.scanImage_stackLoad(fname,false); %Do not subtract the offset

% Correctly index this channel
if ~isempty(STATS.channel)
    chanInd = find(metadata.channelSave == STATS.channel);
    nChans = length(metadata.channelSave);
else
    % Assume we have just one channel
    fprintf('No channel information found in metadata. Assuming there is just one.\n')
    chanInd = 1;
    nChans = 1;
end

% Get this channel
imStack = imStack(:,:,chanInd:nChans:end);
muIm = mean(imStack,3);
muIm_p = mpqc.analyse.convertImageToPhotons(muIm, STATS);
imagesc(floor(muIm_p))
set(gca,'ColorScale','log')
title('log(photon) mean image')
colorbar

%%
% The histogram of photon counts


maxVal = max(muIm_p(:));
if maxVal > 20
    nBins =  round(( maxVal/2.5)/10)*10;
else
    nBins = 10;
end

subplot(2,3,3)
[n,x] = hist(muIm_p(:),nBins);

if isempty(x)
    fprintf('Image conversion failed\n')
    return
end

plot(x,n,'LineWidth',3)
xlabel('Intensity [photons]')
ylabel('log(n)')
set(gca,'YScale','log')
set(gca,'XScale','log')

grid on
xlim([0,max(x)])

hold on
mu = mean(muIm_p(:));
plot([mu,mu], ylim,'r-')
hold off

title(sprintf('log(photons), mean = %0.2f', mu))



nBins =  round(( max(muIm_p(:)) /1.5)/10)*10;
[n,x] = hist(muIm_p(:),nBins);

if isempty(x)
    fprintf('Image conversion failed\n')
    return
end

subplot(2,3,6)
plot(x,n,'LineWidth',3)
xlabel('Intensity [photons]')
ylabel('n')
grid on

hold on
plot([mu,mu], ylim,'r-')
hold off

% set xlim to encompass 99% of values
CLIP=0.99;
p = cumsum(n)/sum(n);
f=find(p>CLIP);
xlim([0,x(f(1))])

title(sprintf('photons clipped at %d%% of max', CLIP*100))
