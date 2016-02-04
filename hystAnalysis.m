%analysis for perovskite hysteresis
% written by Grey Christoforo <first name [at] last name [dot] net>

clear all
%==============EDIT HERE===================
%if you're using the new data file format, this will be overwritten using the value in the file
myArea = 0.12; %cm^2

%set this to true if you don't want any of the output figures to be
%displayed while the script is running (they'll still be saved as .png)
dontDrawFigures = false;

%skip analysis of the first X segments
segmentsToSkip = 0;

%shows the analysis plot for each voltage step (useful to check if the fits are good)
showAnalysisPlots = true;

%set to true if you would like interpolation to be used between datapoints
%to find the maximum power of the device for each simulated IV sweep
%condition (slow)
generateInterpolatedPowerMap = false;

%function describing assumed decay shape
%variables named 'm' and 'b' an 'x0' can not be used.
%'x' must be used as the independant variable
%note that the entry here will be added to the equation of a line and
%approritely time shifted
%assumedDecayFunction = 'c/(x+d)';
%assumedDecayFunction = 'c*exp(-1/tau*x)';
%assumedDecayFunction = 'c*exp(-1/tau*x)+d*exp(-1/tau2*x)';
assumedDecayFunction = 'c*exp(-1/tau*x)+d*exp(-1/tau2*x)';


%you must provide guesses for your variables for the fit algorithm
%normally this doesn't matter too much
%list the guesses ordered alphabetically by variable name used in
%assumedDecayFunction above
%variableGuesses = [0, 1e-3];
%variableGuesses = [0, 1];
variableGuesses = [0, 0, 0.1, 2];

%========STOP EDITING NOW PROBABLY=========
release = version('-release'); %check for matlab version
if str2num(release(1:end-1)) < 2013
    fprintf('Sorry, your MATLAB is too old'\n)
end
[files, path] = uigetfile('*.*','Select data file(s) to load','MultiSelect', 'on');
if isequal(files,0)
    return
end
if ~iscell(files)
    files = {files};
else%if we're in batch mode(multiple files), let's not draw any pots
    fprintf('Batch mode activated!\n')
    showAnalysisPlots = false;
    dontDrawFigures = true;
    segmentsToSkip = 0;
end
for nFile=1:length(files)
    file = files{nFile};
    fprintf('Processing "%s"\n',file);
    
    %read in data and get it ready
    raw = importdata([path file]);
    
    %do we have a new enough version of matlab to support strfind?
    newMatlab = exist('strfind');
    
    if isfield(raw,'textdata') && newMatlab
        nHeaderRows = size(raw.textdata,1);
        for i = 1:nHeaderRows
            thisLine = raw.textdata{i,1};
            
            %pick out area
            if any(strfind(thisLine,'Area'))
                splitString = strsplit(thisLine,' ');
                myArea = str2double(splitString{4});
            end
        end
    end
    
    if isfield(raw,'data') %new format file
        V = raw.data(:,1); % in volts
        I = raw.data(:,2) * 1000/myArea; %let's do current in mA/cm^2
        t = raw.data(:,3);
        status = raw.data(:,4);%TODO: prune data with bad status bits
    else %old format file
        V = raw(:,1); % in volts
        I = raw(:,2) * 1000/myArea; %let's do current in mA/cm^2
        t = raw(:,3);
        status = raw(:,4);%TODO: prune data with bad status bits
    end
    
    %let's make sure the data is in chronological order
    [t,i] = sort(t);
    V = V(i);
    I = I(i);
    t = t - t(1); %time should start from zero
    
    dir = [path [file '.outputs']];
    [~,~,~] = mkdir(dir);
    log=fopen([dir filesep 'log.txt'],'w');
    fprintf(log,'Processing "%s"\n',file);
    segDir = [dir filesep 'segmentPlots'];
    [~,~,~] = mkdir(segDir);
    
    %change current sign if needed (because I hate when I-V curves are
    %upside-down)
    if ((V(1)>V(end)) && (I(1)>I(end))) || ((V(1) < V(end)) && (I(1)<I(end)))
        I = I*-1;
    end
    
    %plot up the raw data
    if dontDrawFigures
        f = figure('Visible','off');
    else
        f = figure;
    end
    
    [AX,H1,H2] = plotyy(t,I,t,V);
    set(get(AX(1),'Ylabel'),'String','Current Density [mA/cm^2]')
    set(get(AX(2),'Ylabel'),'String','Voltage  [V]')
    h = title(file);
    set(h,'interpreter','none')
    grid on
    xlabel('Time [s]')
    print(f,'-dpng',[dir filesep 'data.png'])
    
    if ~newMatlab %end here if we have an old version of matlab
        fprintf('your version of matlab is not new enough to continue')
        return
    end
    
    %sgment the data at the voltage steps
    dV = diff(V);
    pd = fitdist(dV,'Normal');
    nSigmas = 6;
    %assume steps more than 6 sigma away from the mean are keithley voltage
    %changes
    boolStep = (dV < pd.mu-nSigmas*pd.sigma) | (dV > pd.mu+nSigmas*pd.sigma);
    iStep = find(boolStep);
    stepValues = dV(boolStep);
    averageStepSize = mean(stepValues);
    problemSteps = find(averageStepSize*stepValues < 0);
    V(iStep(problemSteps+1)) = [];
    I(iStep(problemSteps+1)) = [];
    dV = diff(V);
    pd = fitdist(dV,'Normal');
    boolStep = (dV < pd.mu-nSigmas*pd.sigma) | (dV > pd.mu+nSigmas*pd.sigma);
    iStep = find(boolStep);
    
    voltageStepsTaken = length(iStep)-1;
    averageDwellTime = mean(diff(t(iStep))); %in seconds
    
    iStart = 1 + segmentsToSkip;
    
    %this is the equation of the line we'll fit the tail to
    lineString = 'm*(x-x0)+b';
    lineSym = sym(lineString);
    line = @(m,b,x0,x) m*(x-x0)+b;
    variableString = '@(m,b,';
    
    decaySym = sym(assumedDecayFunction);
    decaySym = subs(decaySym,'x','x-x0');
    symVars = symvar(decaySym);
    for i = 1:length(symVars)
        charVar = char(symVars(i));
        if ~(strcmp(charVar,'x') || strcmp(charVar,'x0'))
            variableString = [variableString charVar ','];
        end
    end
    variableString = [variableString 'x0,x) '];
    funcString = ['m*(x-x0)+b+' char(decaySym)];
    
    %elementwise e'r'thing:
    funcStringE = strrep(funcString,'*','.*');
    funcStringE = strrep(funcStringE,'/','./');
    funcStringE = strrep(funcStringE,'^','.^');
    
    
    f = str2func([variableString funcStringE]);
    
    %initial guesses for fit variables
    intialGuess = [0 0 variableGuesses];
    
    lineStr = func2str(line);
    lineScanned = textscan(func2str(line),'@(%[^)]) %s');
    lineVarNames = strsplit(lineScanned{1}{1},',');
    lineStr=lineScanned{2}{1};
    funcStr = func2str(f);
    funcScanned = textscan(func2str(f),'@(%[^)]) %s');
    funcVarNames = strsplit(funcScanned{1}{1},',');
    funcStr=funcScanned{2}{1};
    
    %setup the fits
    ft = fitoptions('Method','NonlinearLeastSquares','TolFun',1e-19,'TolX',1e-19,'Display','notify','MaxFunEvals',1e5,'MaxIter',1e4,'Algorithm','Levenberg-Marquardt');
    %ft = fitoptions('Method','NonlinearLeastSquares','TolFun',1e-19,'TolX',1e-9,'Display','notify','MaxFunEvals',1e5,'MaxIter',1e4,'Algorithm','Levenberg-Marquardt');
    %ft = fitoptions('TolFun',1e-12,'TolX',0,'Method','NonlinearLeastSquares','Display','notify','MaxFunEvals',1e5,'MaxIter',1e4,'Algorithm','Levenberg-Marquardt','Normalize','on');
    %ft = fitoptions('Method','NonlinearLeastSquares','TolFun',1e-19,'TolX',1e-19,'Display','notify','MaxFunEvals',1e5,'MaxIter',1e4,'Algorithm','Trust-Region');
    %ft = fitoptions('Method','NonlinearLeastSquares','Algorithm','Trust-Region');
    %ft = fitoptions('Method','NonlinearLeastSquares','Algorithm','Levenberg-Marquardt');
    
    fFitter = fittype(f,'Problem','x0');
    %fFitter = fittype(f);
    lineFitter = fittype(line,'Problem','x0');
    
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
    
    %lets get the voltages even for skipped steps
    for i = 1:voltageStepsTaken
        si = iStep(i)+1;%segment start index
        ei = iStep(i+1);%segment end index
        thisV = V(si:ei);
        thisI = I(si:ei);
        thisVoltage(i) = mean(thisV);
    end
    
    %analyze each segment of the curve
    for i = iStart:voltageStepsTaken
        fprintf(log,'Starting segment %i of %i',i,voltageStepsTaken);
        si = iStep(i)+1;%segment start index
        ei = iStep(i+1);%segment end index
        thist = t(si:ei);
        thisStartT = thist(1);
        x0 = thisStartT;
        thisEndT = thist(end);
        thisI = I(si:ei);
        assumedLinearI = thisI(round(4*end/5):end);
        assumedLineart = thist(round(4*end/5):end);
        
        %print the time (for debugging)
        %x0
        
        %fit a line to the assumed linear region:
        lineGuess = [intialGuess(1),intialGuess(2)];
        
        ft.StartPoint = lineGuess;
        [earlyLine,gof,output] = fit(assumedLineart,assumedLinearI,lineFitter,'problem',x0,ft);
        
        %fit the data for this segment
        fGuess = intialGuess;
        fGuess(1) = earlyLine.m;
        %fGuess(1) = 0;
        fGuess(2) = earlyLine.b;
        %fGuess = [fGuess x0];
        ft.StartPoint = fGuess;
        goodFit = false;
        nTry = 0;
        clear rSquares someFits outputs
        maxTrys = 5;
        while ~goodFit && nTry <=maxTrys
            nTry = nTry+1;
            try
                fittedFunction = [];
                gof = [];
                output = [];
                [fittedFunction,gof,output] = fit(thist,thisI,fFitter,'problem',x0,ft);
                rSquares(nTry) = gof.rsquare;
                someFits{nTry} = fittedFunction;
                outputs{nTry} = output;
                if gof.rsquare > 0.99
                    goodFit = true;
                end
            catch ME %if the fit was not good, try again with a slightly different guess
                rSquares(nTry) = 0;
                someFits{nTry} = [];
                outputs{nTry} = 'fail';
            end
            if ~goodFit
                fprintf(log,'Trying again (R^2 was %0.6f)\n',rSquares(nTry));
                fprintf('Trying again (R^2 was %0.6f)\n',rSquares(nTry));
                fGuess = fGuess+(rand(1,length(fGuess))-0.5)*0.01;
                ft.StartPoint = fGuess;
            end
        end
        [a,b] = max(rSquares);
        gof.rsquare = a;
        fittedFunction = someFits{b};
        output = outputs{b};
        
        %[fittedFunction,gof,output] = fit(thist,thisI,fFitter,'problem',x0,'StartPoint',fGuess,'Algorithm','Levenberg-Marquardt'); %ft

        %now we can extract the fit parameters from the (hopefully) successful
        %fit
        %tmp = coeffvalues(fittedFunction);
        %fitParams(i,:) = tmp(1:end-1);
        fitParams(i,:) = coeffvalues(fittedFunction);
        
        %this assumes the first two terms of the fit function are line slope and intercept
        fittedLine = @(x) line(fitParams(i,1),fitParams(i,2),x0,x);
        
        symbolEqn = decaySym;
        for j = 3:(length(funcVarNames)-2)
            symbolEqn = subs(symbolEqn,funcVarNames{j},fitParams(i,j));
        end
        symbolEqn = subs(symbolEqn,'x0',0);
        qAnalytical(i) = double(int(symbolEqn,0,inf));
        
        R2(i) = gof.rsquare;
        badFitThreshold = 0;
        if R2(i) < badFitThreshold
            reportString = 'bad fit';
            badFit = true;
            fitParams(i,:) = nan(1,size(fitParams,2));
            qAnalytical(i) = nan;
        else %the fit was "good"
            reportString = 'good fit';
            badFit = false;
        end
        reportString = sprintf('Segment %i of %i: Voltage constant at %0.2f V, R^2=%0.5f, %s\n',i,voltageStepsTaken,thisVoltage(i),R2(i),reportString);
        fprintf(reportString)
        fprintf(log,reportString);
        fprintf(log,evalc('disp(output)'));
        %output
        
        if showAnalysisPlots
            f = figure;
        else
            f = figure('Visible','off');
        end
        hold on
        plot(thist,thisI,'.',thist,fittedFunction(thist),'r',thist,fittedLine(thist),'g')
        myxLim = xlim;
        myyLim = ylim;
        if abs(fittedLine(x0)) > abs(fittedFunction(x0))
            h2 = area(thist,fittedLine(thist));
            set(h2,'FaceColor','red','LineStyle','none')
            h1 = area(thist,fittedFunction(thist));
            set(h1,'FaceColor','white','LineStyle','none')
        else
            h1 = area(thist,fittedFunction(thist));
            set(h1,'FaceColor','red','LineStyle','none')
            h2 = area(thist,fittedLine(thist));
            set(h2,'FaceColor','white','LineStyle','none')
        end
        plot(thist,thisI,'.',thist,fittedFunction(thist),'r',thist,fittedLine(thist),'g')
        xlim(myxLim)
        ylim(myyLim)
        
        fLegend = sym(funcString);
        lineLegend = lineSym;
        fLegend = subs(fLegend,'x','t');
        lineLegend = subs(lineLegend,'x','t');
        fLegend = subs(fLegend,'x0',x0);
        lineLegend = subs(lineLegend,'x0',x0);
        
        for j =1:length(fitParams(i,:))
             fLegend = subs(fLegend,funcVarNames{j},fitParams(i,j));
             lineLegend = subs(lineLegend,funcVarNames{j},fitParams(i,j));
        end
        digits(3)
        legend('Segment Data',['y=' char(vpa(fLegend))],['y=' char(vpa(lineLegend))],'Location','SouthOutside')
        
        title(reportString)
        xlabel('Time [s]')
        ylabel('Current [mA/cm^2]')
        hold off
        segName = sprintf('%iof%i.png',i,voltageStepsTaken);
        print(f,'-dpng',[segDir filesep segName])
        
        %build up the data needed to calculate apparent cell power
        for iw = 1:powerMapResolution
            for id = 1:powerMapResolution
                sampleStartTime = thisStartT + delays(id);
                sampleEndTime = thisStartT + delays(id) + windows(iw);
                apparentCurrent(iw,id,i) = mean(fittedFunction(sampleStartTime:1/assummedSamplingFrequency:sampleEndTime));
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
    
    
    header = {'Voltage [V]', 'Best Current Density [mA/cm^2]', 'Worst Current Density [mA/cm^2]'};
    outData = [thisVoltage' apparentCurrent(bestPowI,:)' apparentCurrent(worstPowI,:)'];
    for j =1:length(fitParams(i,:))
        if dontDrawFigures
            f = figure('Visible','off');
        else
            f = figure;
        end
        plot(thisVoltage,fitParams(:,j))
        h = title(file);
        %set(h,'interpreter','none')
        xlabel('Voltage [V]')
        ylabel(funcVarNames{j})
        set(gca,'xdir','reverse')
        grid on
        print(f,'-dpng',[dir filesep funcVarNames{j} '.png'])
        header{j+3} = funcVarNames{j};
        outData(:,j+3) = fitParams(:,j)';
    end
    header{j+4} = 'Q [mC/cm^2]';
    outData(:,j+4) = qAnalytical';
    header{j+5} = 'R^2';
    outData(:,j+5) = R2';
    
    
    
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
    
    
    
    %write out data file
    outFileName = [dir filesep 'numericalOutputs.csv'];
    headerLine = sprintf('%s,',header{:});
    dlmwrite(outFileName,headerLine,'');
    dlmwrite(outFileName,outData,'-append','delimiter',',');
    
    fprintf('Mean R^2 = %0.5f\n',mean(R2));
    fprintf(log,'Mean R^2 = %0.5f\n',mean(R2));
    fclose(log);
    fprintf('Done\n')
end