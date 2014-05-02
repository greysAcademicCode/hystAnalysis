%analysis for eva's perovskite hysteresis

clear all
%==============EDIT HERE===================
%if you're using the new data file format, this will be overwritten using the value in the file
myArea = 0.12; %cm^2 

%set this to true if you don't want any of the output figures to be
%displayed while the script is running (they'll still be saved as .png)
dontDrawFigures = true;

%skip analysis of the first X segments
segmentsToSkip = 1;

%shows the analysis plot for each voltage step (useful to check if the fits are good)
showAnalysisPlots = true;

%set to true if you would like interpolation to be used between datapoints
%to find the maximum power of the device for each simulated IV sweep
%condition (slow)
generateInterpolatedPowerMap = false;

%choose data file sweep direction here
%if you're using the new data file format, this will be overwritten using the value in the file
sweepUp = false;

%========STOP EDITING NOW PROBABLY=========
%read in data and get it ready
[file, path] = uigetfile('*.*');
raw = importdata([path file]);
if isfield(raw,'textdata')
    nHeaderRows = size(raw.textdata,1);
    for i = 1:nHeaderRows
        thisLine = raw.textdata{i,1};
        
        %pick out area
        if any(strfind(thisLine,'Area'))
            splitString = strsplit(thisLine,' ');
            myArea = str2double(splitString{4});
        end
        
        %pick out sweep direction
        if any(strfind(thisLine,'sweepUp'))
            splitString = strsplit(thisLine,' ');
            sweepUp = strcmp(splitString{4},'1');
        end
        
        
    end
end
if isfield(raw,'data') %new format file
    V = raw.data(:,1); % in volts
    I = raw.data(:,2) * 1000/myArea; %let's do current in mA/cm^2
    t = raw.data(:,3) - raw.data(1,3); %in seconds
    status = raw.data(:,4);%TODO: prune data with bad status bits
else %old format file
    V = raw(:,1); % in volts
    I = raw(:,2) * 1000/myArea; %let's do current in mA/cm^2
    t = raw(:,3) - raw(1,3); %in seconds
    status = raw(:,4);%TODO: prune data with bad status bits
end

dir = [path [file '.outputs']];
[~,~,~] = mkdir(dir);

%change current sign if needed (because I hate when I-V curves are
%upside-down)
if ((V(1)>V(end)) && (I(1)>I(end))) || ((V(1) < V(end)) && (I(1)<I(end)))
    I = I*-1;
end

%sgment the data at the voltage steps
dV = diff(V);
pd = fitdist(dV,'Normal');
nSigmas = 6;
%assume steps more than 6 sigma away from the mean are keithley voltage
%changes
boolStep = (dV < pd.mu-nSigmas*pd.sigma) | (dV > pd.mu+nSigmas*pd.sigma);
boolStep(1) = true;%put marker at start
boolStep(end) = true;%put marker at end
iStep = find(boolStep);

voltageStepsTaken = length(iStep)-1;
averageDwellTime = mean(diff(t(iStep))); %in seconds

%plot up the raw data
if dontDrawFigures
    f = figure('Visible','off');
else
    f = figure;
end

hold on
[AX,H1,H2] = plotyy(t,I,t,V);
[hAx,hLine1,hLine2] = plotyy(t(iStep),I(iStep),t(iStep),V(iStep));
set(hLine1,'Marker','o','LineStyle','none')
set(hLine2,'Marker','o','LineStyle','none')
hold off
set(get(AX(1),'Ylabel'),'String','Current Density [mA/cm^2]')
set(get(AX(2),'Ylabel'),'String','Voltage  [V]')
h = title(file);
set(h,'interpreter','none')
grid on
xlabel('Time [s]')
print(f,'-dpng',[dir filesep 'data.png'])

iStart = 1 + segmentsToSkip;

%this is the equation of the line we'll fit the tail to
line = @(m,b,x) m*x+b;
%this is the equation of the line plus the exponential decay
%p(1) is the line slope
%p(2) is the y offset
%p(3) is tau
%p(4) is the time delay variable
f = @(p,x) p(1)*(x-p(4))+p(2)+exp(-1/p(3)*(x+p(4)));

%probably only need to get the signs right here for the fit to converge...
if sweepUp
    initialGuess = [-1 0 1 1];
else %then sweep down
    initialGuess = [-1 0 -1 1];
end

%set some fit options here
options =optimset('TolFun',1e-19,'TolX',0,'Algorithm','levenberg-marquardt','Display','off');

%filterWindow = 21;
%golayOrder = 3;

powerMapResolution = 100; %total power data points will be this number^2

%here we set up the bounds for the parameters of the simulated IV
%measurement system
minDelay = 0;
maxDelay = averageDwellTime*0.95/2;
delays = linspace(minDelay,maxDelay,powerMapResolution);

minWindowLength = 0.001;
maxWindowLength = averageDwellTime*0.95/2;
windows = linspace(minWindowLength,maxWindowLength,powerMapResolution);

%preallocate this guy because he could get big:
apparentCurrent = zeros(powerMapResolution,powerMapResolution,voltageStepsTaken);

%we'll assume that the simulated IV measurement system samples this fast
assummedSamplingFrequency = 1000; %in Hz

%get lets get the voltages even for skipped steps
for i = 1:voltageStepsTaken
    si = iStep(i)+1;%segment start index
    ei = iStep(i+1);%segment end index
    thisV = V(si:ei);
    thisVoltage(i) = mean(thisV);
end

%analyze each segment of the curve
for i = iStart:voltageStepsTaken
    fprintf('Fitting segment %i of %i\n',i,voltageStepsTaken)
    si = iStep(i)+1;%segment start index
    ei = iStep(i+1);%segment end index
    thist = t(si:ei);
    thisStartT = thist(1);
    thisEndT = thist(end);
    thisI = I(si:ei);
    
    newGuess = initialGuess;
    %need to do a bit better on the delay variable guess for the fit to
    %converge
    newGuess(4) = initialGuess(4) - thisStartT;
    
    %fit the data for this segment
    f1=lsqcurvefit(f,newGuess,thist,thisI,[],[],options);
    thisF = @(x)f(f1,x);
    
    %now we can extract the fit parameters from the (hopefully) successful
    %fit
    tau(i) = f1(3);
    b = f1(2) - f1(1)*f1(4);%y-intercept
    m(i) = f1(1);%slope
    thisLine = @(x)line(m(i),b,x);
    
    %intigrate under the analytical expressions as found by the fits
    lineArea = integral(thisLine,thisStartT,thisEndT);
    curveArea = integral(thisF,thisStartT,thisEndT);
    
    %probably no real reason to smooth out the noise in the data...
    %thisISmooth = sgolayfilt(thisI,golayOrder,filterWindow);
    
    %ensuring that these match is a nice sanity check
    %if the numerical intigration roughly matches the analytical
    %intigration that's great, otherwise something is botched
    qAnalytical(i) = curveArea - lineArea; %in mili-columbs
    qNumerical = trapz(thist,thisI) - lineArea; %in mili-columbs
    
    if showAnalysisPlots
        figure
        hold on
        %plot(thist,thisI,'.',thist,f(f1,thist),'r',thist,thisISmooth,'g')
        plot(thist,thisI,'.',thist,thisF(thist),'r',thist,thisLine(thist),'g')
        myxLim = xlim;
        myyLim = ylim;
        if thisI(end) < 0
            h2 = area(thist,thisLine(thist));
            set(h2,'FaceColor','red','LineStyle','none')
            h1 = area(thist,thisF(thist));
            set(h1,'FaceColor','white','LineStyle','none')
        else
            h2 = area(thist,thisF(thist));
            set(h2,'FaceColor','red','LineStyle','none')
            h1 = area(thist,thisLine(thist));
            set(h1,'FaceColor','white','LineStyle','none')
        end
        
        xlim(myxLim)
        ylim(myyLim)
        plot(thist,thisI,'.',thist,thisF(thist),'r',thist,thisLine(thist),'g')
        words = sprintf('Voltage constant at %0.2f V',thisVoltage(i));
        title(words)
        xlabel('Time [s]')
        ylabel('Current [mA/cm^2]')
        hold off
    end
    
    %build up the data needed to calculate apparent cell power
    for iw = 1:powerMapResolution
        for id = 1:powerMapResolution
            sampleStartTime = thisStartT + delays(id);
            sampleEndTime = thisStartT + delays(id) + windows(iw);
            apparentCurrent(iw,id,i) = mean(thisF(sampleStartTime:1/assummedSamplingFrequency:sampleEndTime));
        end
    end
end

%voltage step
voltageStepSize = abs(thisVoltage(1) - thisVoltage(end))/voltageStepsTaken;

%put NaNs in the proper places if we did not do the analysis on the first
%segments
if segmentsToSkip > 0
    for i =1:segmentsToSkip
        tau(i) = nan;
        m(i) = nan;
        qAnalytical(i) = nan;
        %thisVoltage(i) = nan;
        apparentCurrent(:,:,i) = nan;
    end
end

apparentCurrent = reshape(apparentCurrent,[],voltageStepsTaken);
apparentPower = bsxfun(@times,apparentCurrent,thisVoltage);
pce = nanmax(apparentPower,[],2);
[bestPow,bestPowI] = max(pce);
[worstPow,worstPowI] = min(pce);

%plot all the possible IV curves on top of eachother
if dontDrawFigures
    f = figure('Visible','off');
else
    f = figure;
end
plot(thisVoltage,apparentCurrent)
h = title(file);
set(h,'interpreter','none')
xlabel('Voltage [V]')
ylabel('Current [mA/cm^2]')
grid on
print(f,'-dpng',[dir filesep 'all_iv_curves.png'])

%generate the power map here
if generateInterpolatedPowerMap
    fprintf('Generating interpolated power map\r',i,voltageStepsTaken)
    ft = fittype( 'smoothingspline' );
    
    pce = zeros(1,powerMapResolution^2);
    warning off
    pSize = powerMapResolution^2;
    for i = 1:pSize
        fprintf('Finding max power %i of %i\n',i,pSize)
        % Fit model to data.
        [p,in] = max(thisVoltage.*apparentCurrent(i,:));
        [xData, yData] = prepareCurveData( thisVoltage, apparentCurrent(i,:) );
        [fitresult, gof] = fit( xData,yData , ft );
        invPower = @(x) fitresult(x)*x*-1;
        x0 = thisVoltage(in);
        vMax = fminsearch(invPower,x0);
        pce(i) = vMax*fitresult(vMax);
    end
    warning on
end

pceMin = min(pce);
pceMax = max(pce);

pce = reshape(pce,powerMapResolution,powerMapResolution);

if dontDrawFigures
    f = figure('Visible','off');
else
    f = figure;
end
imagesc(delays,windows,pce)
axis square
xlabel('Delay [s]')
ylabel('Averaging window length [s]')
myTitle = sprintf('Maximum PCE = %0.2f%% Minimum PCE = %0.2f%%',pceMax,pceMin);
title(myTitle)
colorbar
print(f,'-dpng',[dir filesep 'powerMap.png'])

if dontDrawFigures
    f = figure('Visible','off');
else
    f = figure;
end
plot(thisVoltage,tau)
h = title(file);
set(h,'interpreter','none')
xlabel('Voltage [V]')
ylabel('\tau [s]')
set(gca,'xdir','reverse')
grid on
print(f,'-dpng',[dir filesep 'tau.png'])

if dontDrawFigures
    f = figure('Visible','off');
else
    f = figure;
end
plot(thisVoltage,qAnalytical)
xlabel('Voltage [V]')
ylabel('Charge Stored [mC/cm^2]')
set(gca,'xdir','reverse')
h = title(file);
set(h,'interpreter','none')
grid on
print(f,'-dpng',[dir filesep 'q.png'])

if dontDrawFigures
    f = figure('Visible','off');
else
    f = figure;
end
plot(thisVoltage,m)
xlabel('Voltage [V]')
ylabel('Linear Decay [mA/cm^2/s]')
h = title(file);
set(h,'interpreter','none')
set(gca,'xdir','reverse')
grid on
print(f,'-dpng',[dir filesep 'decayRate.png'])

%write out data file
outFileName = [dir filesep 'numericalOutputs.csv'];
header = {'Voltage [V]', 'Best Current Density [mA/cm^2]', 'Worst Current Density [mA/cm^2]', 'Tau [s]', 'Q [mC/cm^2]','Linear Decay Rate [mA/cm^2/s]'};
headerLine = sprintf('%s,',header{:});
dlmwrite(outFileName,headerLine,'');
outData = [thisVoltage' apparentCurrent(bestPowI,:)' apparentCurrent(worstPowI,:)' tau' qAnalytical' m'];
dlmwrite(outFileName,outData,'-append','delimiter',',');

fprintf('Done\n')
