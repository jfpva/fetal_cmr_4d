function S = cardsync_intraloc( S, varargin )
%CARDSYNC_INTRALOC  estimate heart rate within each slice-location
%

%   jfpva (joshua.vanamerom@kcl.ac.uk)


%% Optional Input Argument Default Values

default.hrRange         = [105,180]; % heart rate range (bpm)
default.resultsDir      = pwd;      % path of directory to save results
default.isVerbose       = true;


%% Parse Input

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end


addRequired(   p, 'S', ... 
    @(x) validateattributes( x, {'struct'}, {'vector'}, mfilename) );

add_param_fn(   p, 'hrrange', default.hrRange, ...
    @(x) validateattributes( x, {'numeric'}, {'positive','numel',2}, mfilename) );

add_param_fn(   p, 'resultsdir', default.resultsDir, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(   p, 'verbose', default.isVerbose, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

parse( p, S, varargin{:} );

hrRange         = p.Results.hrrange;
resultsDir      = p.Results.resultsdir;
isVerbose       = p.Results.verbose;


%% Setup

nStack = numel( S );


%% Estimate Heart Rates

for iStk = 1:nStack
 
    S(iStk).hrEstRelPeakProminence = nan(1,S(iStk).nLoc);
    S(iStk).hrEstPeakWidth = nan(1,S(iStk).nLoc);
    S(iStk).hrEstPeakLocDiff = nan(1,S(iStk).nLoc);
     
    for iLoc = 1:S(iStk).nLoc
 
        % Estimate heart rate
        if any( S(iStk).H{iLoc}(:) )
            [ S(iStk).tRR{iLoc}, tRTrigger, S(iStk).hrEstRelPeakProminence(iLoc), S(iStk).hrEstPeakWidth(iLoc), S(iStk).hrEstPeakLocDiff(iLoc) ] = estimate_heartrate_xf( permute( S(iStk).Y{iLoc}, [1,2,4,3] ), S(iStk).frameDuration, 'roi', S(iStk).H{iLoc}, 'hrRange', hrRange, 'verbose', true, 'visible', 'off' );
            S(iStk).tRTrigger{iLoc} = tRTrigger + S(iStk).tFrame{iLoc}(1);  % offset to time of first frame in slice
            hFig = gcf;
            hFig.Name = sprintf( '%ssl%02i_%s', S(iStk).desc, iLoc, hFig.Name );
            title( sprintf( '%ssl%02i', S(iStk).desc, iLoc ) )
            save_figs( pngfigDir, hFig, matfigDir )
            S(iStk).figs(iLoc).heartrateEstimate = strcat( hFig.Name, '.png' );
            close( hFig ),
        else
            S(iStk).tRR{iLoc} = NaN;
            S(iStk).tRTrigger{iLoc} = NaN;
            S(iStk).figs(iLoc).heartrateEstimate = '';
        end
         
        % Calculate cardiac phase 
        if ~any( isnan( S(iStk).tRTrigger{iLoc} ) )
            [ ~, cardPhaseFraction ] = calc_cardiac_timing( S(iStk).tFrame{iLoc}, S(iStk).tRTrigger{iLoc} );
            S(iStk).thetaFrame{iLoc} = 2 * pi * cardPhaseFraction;
        else
            S(iStk).thetaFrame{iLoc} = nan( size( S(iStk).tFrame{iLoc} ) );
        end
         
    end
     
    S(iStk).tRRest = S(iStk).tRR;
     
end


%% Find and Replace Outliers

% Find Outliers Using Peak Height
hrEstRelPeakProminence = [S.hrEstRelPeakProminence];
indVal1 = find(~isnan( hrEstRelPeakProminence ));
[isOutlierPeakHeightIndVal,peakHeightThreshUpper,peakHeightThreshLower,peakHeightCentre] = isoutlier( 1./hrEstRelPeakProminence( indVal1 ) );
indOutlierPeakHeight = indVal1( isOutlierPeakHeightIndVal );  % outliers have relative peak prominence^-1 > 3 * median absolute deviations (MAD) above the median
                                                       % MAD = median(abs(hrEstPeakWidth(indVal1)-median(hrEstPeakWidth(indVal1))))
indOutlierPeakHeight = indOutlierPeakHeight( 1./hrEstRelPeakProminence( indOutlierPeakHeight ) > peakHeightCentre );  % only use estimates with large peak heights
isOutlierPeakHeight  = false( size( hrEstRelPeakProminence ) );
isOutlierPeakHeight( indOutlierPeakHeight ) = true;


% Find Outliers Using Peak Width
hrEstPeakWidth = [S.hrEstPeakWidth];
indVal2 = setdiff( find(~isnan(hrEstRelPeakProminence )), indOutlierPeakHeight );
[isOutlierPeakWidthIndVal,peakWidthThreshLower,peakWidthThreshUpper,peakWidthCentre] = isoutlier( hrEstPeakWidth( indVal2 ) );
indOutlierPeakWidth = indVal2( isOutlierPeakWidthIndVal );  % outliers have peak width > 3 * median absolute deviations (MAD) above the median                       
indOutlierPeakWidth = indOutlierPeakWidth( hrEstPeakWidth( indOutlierPeakWidth ) > peakWidthCentre );  % use estimates with small peak widths
isOutlierPeakWidth  = false( size( hrEstPeakWidth ) );
isOutlierPeakWidth( indOutlierPeakWidth ) = true;

% Combine Outliers
indOutlier = union( indOutlierPeakHeight, indOutlierPeakWidth );

% Identify Inliers
isInlier  = true( size( [S.hrEstPeakWidth] ) );
isInlier(indOutlier) = false;

% Replace Outliers with Interpolation of Inliers and Update Plot
for iStk = 1:nStack
    tRR             = cell2mat(S(iStk).tRRest);
    isInlierLoc     = ~isnan( tRR ) & isInlier( sum([S(1:(iStk-1)).nLoc]) + (1:S(iStk).nLoc) );
    indInlierLoc    = find( isInlierLoc );
    tRR(~isInlierLoc) = NaN;
    tRR = interp1( indInlierLoc, tRR(indInlierLoc), 1:S(iStk).nLoc, 'linear' );  % linear interpolation for slice-locations between inlier heart rate estimates
    indInlierLoc2   = find( ~isnan( tRR ) );
    tRR = interp1( indInlierLoc2, tRR(indInlierLoc2), 1:S(iStk).nLoc, 'nearest', 'extrap' ); % nearest-neighbour interpolation to extrapolate edge slice-locations
    S(iStk).hrOutlierPeakHeight = isOutlierPeakHeight( sum([S(1:(iStk-1)).nLoc]) + (1:S(iStk).nLoc) );
    S(iStk).hrOutlierPeakWidth  = isOutlierPeakWidth( sum([S(1:(iStk-1)).nLoc]) + (1:S(iStk).nLoc) );
    S(iStk).hrOutlier = S(iStk).hrOutlierPeakHeight | S(iStk).hrOutlierPeakWidth; 
    for iLoc = find( S(iStk).hrOutlier )
        % Update Timing
        rrInterval  = tRR(iLoc);    
        nTrigger    = ceil( numel(S(iStk).tFrame{iLoc}) * S(iStk).frameDuration / rrInterval );
        tRTrigger   = rrInterval * (0:nTrigger);
        S(iStk).tRR{iLoc} = rrInterval;
        S(iStk).tRTrigger{iLoc} = tRTrigger + S(iStk).tFrame{iLoc}(1);
        [ ~, cardPhaseFraction ] = calc_cardiac_timing( S(iStk).tFrame{iLoc}, S(iStk).tRTrigger{iLoc} );
        S(iStk).thetaFrame{iLoc} = 2 * pi * cardPhaseFraction;
    end
end


%% Plot R-R v Time

if ( isVerbose )
    
    % Get Frame Timing
    i=1;
    tStart = nan( size( [S.tRR] ) );
    tEnd = nan( size( [S.tRR] ) );
    for iStk = 1:nStack
        for iLoc = 1:S(iStk).nLoc  
            tStart(i) = S(iStk).tFrame{iLoc}(1) / 60;  % min
            tEnd(i)   = ( S(iStk).tFrame{iLoc}(end)+S(iStk).frameDuration ) / 60;  % min 
            i = i + 1;
        end
    end
    t0 = 0.5;
    tOffset = min(tStart) - t0;
    tStart = tStart - tOffset;
    tEnd   = tEnd - tOffset;
    xLim = [ 0, max(tEnd)+t0 ];
    
    % Colours and Markers
    nLocMax = max([S.nLoc]);
    rgbPeakHeight = [171,217,233]/255; 
    rgbPeakWidth  = [44,123,182]/255; 
    rgbOut        = [215,25,28]/255; 
    rgbIn         = 0.5*ones(nLocMax+1,3); 
    markerIconOut = '+';
    markerSizeOut = 5;

    % Figure
    hFig = figure( 'Name', 'rr_peak_height_width_v_time' );
    hFig.Position(1) = 0;
    hFig.Position(3) = .8*hFig.Position(3);
    hFig.Position(4) = .8*hFig.Position(4);
    hFig.Position(3) = range(xLim)/10 * hFig.Position(3);
    
    % Plot Estimated HR
    subplot(7,1,1:3)
    hold on 
    i=1;
    for iStk = 1:nStack
        for iLoc = 1:S(iStk).nLoc  
            hr     = 0.06./S(iStk).tRRest{iLoc};
            if isOutlierPeakHeight(i)
                h = plot( mean([tStart(i),tEnd(i)]), hr*1000, markerIconOut, 'Color', rgbOut, 'LineWidth', 2, 'MarkerSize', markerSizeOut  );
                h.DisplayName = 'outlier due to peak height';
            elseif isOutlierPeakWidth(i)
                h = plot( mean([tStart(i),tEnd(i)]), hr*1000, markerIconOut, 'Color', rgbOut, 'LineWidth', 2, 'MarkerSize', markerSizeOut  );
                h.DisplayName = 'outlier due to peak width';
            end
            i=i+1;
        end
    end
    hold off
    hAx1 = gca;
    hAx1.XLim = xLim;
    hAx1.YLim = hrRange;  %set(hAx1,'YLim',flip(60000./hrRange))
    hAx1.XTick = 0:1:tEnd(end);
    hAx1.XTickLabel = {};
    ylabel( sprintf('Heart Rate\n[bpm]'), 'FontSize', 14 )  %ylabel( 'R-R Interval [ms]', 'FontSize', 14 )
    hAx1.XAxis.MinorTickValues = setdiff( 0:(1/6):xLim(2), 0:1:xLim(2) );
    hAx1.YAxis.MinorTickValues = setdiff( 0:5:250, hAx1.YAxis.TickValues );
    grid on
    grid minor
    box on

    % Plot Corrected HR
    i=1;
    hold on
    for iStk = 1:nStack
        for iLoc = 1:S(iStk).nLoc
            hr     = 0.06./S(iStk).tRR{iLoc};
            h = plot( [tStart(i),tEnd(i)], hr*1000*[1,1], 'LineWidth', 4, 'Color', rgbIn(iLoc,:) );
            h.DisplayName = sprintf( 'slice %02i used', iLoc );
            i=i+1;
        end
    end
    hold off
    
    % Plot Peak Prominence
    subplot(7,1,4:5)
    hold on
    i=1;
    peakH = nan( size( [S.tRR] ) );
    plot( xLim, 1/peakHeightThreshLower*[1,1], 'Color', 0.5*[1,1,1], 'LineWidth', 1 )
    for iStk = 1:nStack
        for iLoc = 1:S(iStk).nLoc  
            peakH(i)  = S(iStk).hrEstRelPeakProminence(iLoc);
            if isOutlierPeakHeight(i)
                plot( mean([tStart(i),tEnd(i)]), peakH(i), markerIconOut, 'Color', rgbOut, 'LineWidth', 2, 'MarkerSize', markerSizeOut );
                
            else
                plot( [tStart(i),tEnd(i)], peakH(i)*[1,1], 'LineWidth', 4, 'Color', rgbPeakHeight );
            end
            i=i+1;
        end
    end
    hold off
    hAx2 = gca;
    set(hAx2,'XLim',xLim)
    hAx2.XTick = hAx1.XTick;
    hAx2.YLim(1) = 0;
    hAx2.XAxis.MinorTickValues = hAx1.XAxis.MinorTickValues;
    box on
    grid on
    grid minor
    ylabel(sprintf('Peak Height\n[a.u.]'),'FontSize',14)
    set(gca,'XTickLabel',{})

    % Plot Peak Width
    subplot(7,1,6:7)
    hold on
    plot( xLim, 60*peakWidthThreshUpper*[1,1], 'Color', 0.5*[1,1,1], 'LineWidth', 1 )
    i=1;
    peakW = nan( size( [S.tRR] ) );
    for iStk = 1:nStack
        for iLoc = 1:S(iStk).nLoc  
            peakW(i)  = S(iStk).hrEstPeakWidth(iLoc);
            if isOutlierPeakWidth(i)
                plot( mean([tStart(i),tEnd(i)]), 60*peakW(i), markerIconOut, 'Color', rgbOut, 'LineWidth', 2, 'MarkerSize', markerSizeOut );
            else
                plot( [tStart(i),tEnd(i)], 60*peakW(i)*[1,1], 'LineWidth', 4, 'Color', rgbPeakWidth );
            end
            i = i+1;
        end
    end
    hold off
    hAx3 = gca;
    set(hAx3,'XLim',xLim)
    hAx3.XTick = hAx1.XTick;
    XTickLabel = hAx3.XTickLabel;
    for iT = 2:2:numel(XTickLabel)
        XTickLabel{iT} = '';
    end
    hAx3.XTickLabel = XTickLabel;
    hAx3.YLim(1) = 0;
    hAx2.XAxis.MinorTickValues = hAx1.XAxis.MinorTickValues;
    box on
    grid on
    grid minor
    ylabel(sprintf('Peak Width\n[bpm]'), 'FontSize', 14)
    xlabel( 'Acquisition Time [min]', 'FontSize', 14 )
    
    % Align YLabels
    yLabelXOffset = min( [ hAx1.YLabel.Position(1) hAx2.YLabel.Position(1) hAx3.YLabel.Position(1) ] );
    hAx1.YLabel.Position(1) = yLabelXOffset;
    hAx2.YLabel.Position(1) = yLabelXOffset;
    hAx3.YLabel.Position(1) = yLabelXOffset;

end


%% Save Figure

if exist( 'hFig', 'var' )
    pngfigDir   = fullfile( resultsDir, 'figs', 'png' );
    matfigDir   = fullfile( resultsDir, 'figs', 'fig' );
    epsfigDir   = fullfile( resultsDir, 'figs', 'eps' );
    if ~exist( fullfile( resultsDir, 'figs' ), 'dir' )
        mkdir( fullfile( resultsDir, 'figs' ) )
    end
    save_figs( pngfigDir, hFig, matfigDir, epsfigDir )
end


end  % cardsync_intraloc(...)

