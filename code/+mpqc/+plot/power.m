function power(fname,varargin)
% Plots of power curve compared to scanimage percent power
%
% mpqc.plot.power(fname,varargin)
%
% Purpose
% Plots of the comparison of power at the sample plane versus scan image
% percent power and estimated poer. Curves should be close to identical, if
% not the power scanimage reports is not the true power at the sample
% plane and should be updated to reflect the true power.
%
% power measurement structure loads:
% observedPower_mW
% SIpower_mW
% powerSeriesPercent
% currentTime
% laserWavelength
% fittedMinAndMax
% sensorName
%
% Isabell Whiteley, SWC AMF, initial commit 2025

measurements = []; % not needed if just using the loaded powerMeasurements

if ~exist(fname,'file')
    fprintf('%s does not exist. mpqc.plot.power will not load it.\n',fname)
    return
else 
    load(fname, "-mat")
end

% Make a new figure or return a plot handle as appropriate
fig = mpqc.tools.returnFigureHandleForFile([fname,mfilename]);

%%
figure
subplot(1,2,1)
title(sprintf('Wavelength = %d nm', powerMeasurements.laserWavelength))
P.observed = plot(repmat(powerMeasurements.powerSeriesPercent,1,...
    size(powerMeasurements.observedPower_mW,2))',powerMeasurements.observedPower_mW(:),'.k');
hold on
P.meanObs =  plot(powerMeasurements.powerSeriesPercent, mean(powerMeasurements.observedPower_mW,2),'-r');
P.SIpower = plot(powerMeasurements.powerSeriesPercent, powerMeasurements.SIpower_mW, '-b'); % removed the x1000 from SIpower
P.fitObs = plot([0,100],[powerMeasurements.fittedMinAndMax(1),powerMeasurements.fittedMinAndMax(2)], '-r','LineWidth',2); % hack method
legend({'','Mean Observed Power', 'SI Power', 'Observed Power fit'}, 'Location', 'NorthWest');
xlabel('Percent Power')
ylabel('Power (mW)')
hold off

subplot(1,2,2)
 plot(mean(powerMeasurements.observedPower_mW,2), (powerMeasurements.SIpower_mW - mean(powerMeasurements.observedPower_mW,2)),...
 'ok')%, 'MarkerFaceColor'), [1,1,1]*0.5,)
 xlabel('Observed Power (mW)')
ylabel('SI\_Power - Observed\_Power (mW)')