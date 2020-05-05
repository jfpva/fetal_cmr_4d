function H = plot_displacement_and_outrej_v_time( I, D, reconDir, resultsDir, showImages, showOutrej )
%PLOT_DISPLACEMENT_AND_OUTREJ_V_TIME  plots volumetric reconstruction transformation and outlier rejection results
%
%   Example
%
%       infoFile   = '/path/to/info.tsv';
%       reconDir   = '/path/to/recon';
%       resultsDir = '/path/to/save/results';
%
%       showImages = true;
%       showOutrej = true;
%
%       I = read_info_tsv( infoFile );
%       D = summarise_transformations( reconDir, I );
%       H = plot_displacement_and_outrej_v_time( I, D, reconDir, resultsDir, showImages, showOutrej );
%
%   See also read_info_tsv, summarise_transformations.

%   jfpva (joshua.vanamerom@kcl.ac.uk) 


%% Optional Input Arguments

if ~exist( 'showImages', 'var' )
    showImages = false;
end

if ~exist( 'showOutrej', 'var' )
    showOutrej = true;
end


%% Pre-Process

stackID = unique( I.StackIndex );

% stack description
stackDesc = cell( size( stackID ) );
for iSt = stackID(:)'+1
    stackFilePath = I.File(find(I.StackIndex==stackID(iSt),1));  
    [ ~, desc, ~ ] = fileparts( stackFilePath.char );
    while any( desc == '.' )
        [ ~, desc, ~ ] = fileparts( desc );
    end
    stackDesc{iSt} = desc;
end

% slice-location IDs
sliceID   = cell( size( stackID ) );
for iSt = stackID(:)'+1
    sliceID{iSt} = unique( I.StackLocIndex( I.StackIndex == stackID(iSt) ) );
end


%% Load Images and Masks 

if ( showImages )   
    rltFile = unique( I.File );
    M = struct( 'im', [], 'mask', [], 'maskCenter', [], 'displayWindow', [], 'displayXLim', [], 'displayYLim', [], 'displaySIRange', [] );
    for iStack = 1:numel(stackID)
        dcFile = char( fullfile( fileparts( reconDir ), strrep( strrep( rltFile(iStack), '_rlt_', '_dc_' ), '../', '' ) ) );
        N = load_untouch_nii( dcFile );
        M(iStack).im = N.img;
        maskFile = char( fullfile( fileparts( reconDir ), strrep( strrep( rltFile(iStack), '_rlt_ab', '_mask_heart' ), '../data', 'mask' ) ) );
        N = load_untouch_nii( maskFile );
        M(iStack).mask = N.img;
        ind = find( M(iStack).mask );
        [rowMask,colMask,~] = ind2sub( size(M(iStack).mask), ind );
        M(iStack).maskCenter = round( [ mean(colMask), mean(rowMask) ] );
        M(iStack).displayWindow = [ 76, 60 ];
        M(iStack).displayXLim = M(iStack).maskCenter(1) + M(iStack).displayWindow(1)/2 * [-1,+1];
        M(iStack).displayYLim = M(iStack).maskCenter(2) + M(iStack).displayWindow(2)/2 * [-1,+1];
        vec = @(x) x(:);
        imVec = vec( M(iStack).im(M(iStack).displayYLim(1):M(iStack).displayYLim(2),M(iStack).displayXLim(1):M(iStack).displayXLim(2),:) );
        M(iStack).displaySIRange = [ prctile( imVec, 1 ), prctile( imVec, 99.9 ) ];
    end
end


%% Setup Plot

showAlParam     = true;

SubFigSizeX = 150;
SubFigSizeY = 125;

setlocation = @( h ) setappdata( h, 'SubplotDefaultAxesLocation', [0.05, 0.025, 0.925, 0.95] );
setinset    = @( h ) setappdata( h, 'SubplotDefaultInset', [0, 0, 0, 0]);

nRow = 3;

if showImages
    nRow = nRow + 1;
end

if showOutrej
    nRow = nRow + 1;
end

limDisplacement = 10;  % mm
limTranslation  = limDisplacement;
limRotation     = limDisplacement;

cX              = [223,194,125]/255;
cY              = [166,97,26]/255;
cZ              = [1,133,113]/255;
cWeight         = [253,174,97]/255; 
cDisp           = [128,205,193]/255; 
cDev            = [31,120,180]/255; 

lineWidth0      = 1.5;

calc_dist = @(A1,A2,worldCoordinates) sqrt( sum( ( A1*worldCoordinates-A2*worldCoordinates ).^2, 1 ) );
calc_disp = @(dist) mean( dist );

plot_Ak = @(t,v,c,l) plot(t,v,'-', 'Color',c,'LineWidth',l*lineWidth0); 
plot_AkX = @(t,v) plot_Ak(t,v,cX,1.25); 
plot_AkY = @(t,v) plot_Ak(t,v,cY,1); 
plot_AkZ = @(t,v) plot_Ak(t,v,cZ,1); 

plot_Al = @(t,v,c,l) plot(t,v*ones(size(t)),':','Color',c,'LineWidth',l*lineWidth0);
plot_AlX = @(t,v) plot_Al(t,v,hsv2rgb(rgb2hsv(cX).^[1,1,0.5]),1.375);
plot_AlY = @(t,v) plot_Al(t,v,hsv2rgb(rgb2hsv(cY).^[1,1,0.5]),1.375);
plot_AlZ = @(t,v) plot_Al(t,v,hsv2rgb(rgb2hsv(cZ).^[1,1,0.5]),1.375);

xMajorTickIncrement = 2;
xMinorTickIncrement = 1;


%% Get SliceTiming

rltFiles        = I.File(I.StackLocIndex==0&I.StackDynIndex==0);
rlt2paramFile   = @(rltFile) fullfile( reconDir, strrep( strrep(rltFile,'ab.nii.gz','parameters.mat'), '/data/', '/ktrecon/' ) );
sliceStartTime  = cell( size( stackID ) );
sliceStartTime0 = inf;
for iSt = stackID(:)'+1
    load( rlt2paramFile(rltFiles(iSt)), 'PARAM' );
    sliceStartTime{iSt} = PARAM.Timing.sliceTime;
    sliceStartTime0 = min( [ sliceStartTime0 sliceStartTime{iSt} ] );
end
for iSt = stackID(:)'+1
    sliceStartTime{iSt} = sliceStartTime{iSt} - sliceStartTime0;
end
clear PARAM


%% Create Figures

for iSt = stackID(:)'+1

    sliceNo = sliceID{iSt}' + 1;
    nSlice  = numel(sliceNo) + 1;

    hFig = figure('Name',[stackDesc{iSt},'_motion_correction'],'Position',[0,0,nSlice*SubFigSizeX,nRow*SubFigSizeY]);

    H(iSt) = hFig;

    setlocation( hFig );
    setinset( hFig );

    t0 = sliceStartTime{iSt};

    hAxRef = subplot(nRow,nSlice,(nRow-1)*nSlice+1);
    hAxXOffset = hAxRef.Position(1) * 1.5;
    if nRow == 3
        hAxYOffset = hAxRef.Position(2) * 4.0;
    elseif nRow == 4
        hAxYOffset = hAxRef.Position(2) * 3.0;
    elseif nRow == 5
        hAxYOffset = hAxRef.Position(2) * 2.4;
    else
        hAxYOffset = hAxRef.Position(2) * 1.5;
    end
    hAxWidth   = hAxRef.Position(3);
    hAxHeight  = hAxRef.Position(4);
    hAxXSpace  = (1-(hAxWidth*nSlice+2*hAxXOffset))/(nSlice-1)/2;
    hAxYSpace  = (1-(hAxHeight*nRow+1.5*hAxYOffset))/(nRow-1);
    delete(hAxRef)
    set_ax_position = @(h,iRow,iCol) set(h,'Position',[hAxXOffset+(iCol-1)*(hAxWidth+hAxXSpace),hAxYOffset+(nRow-iRow)*(hAxHeight+hAxYSpace),hAxWidth,hAxHeight]);

    hAx = cell(nRow,nSlice);
    hLgnd = cell(nRow,1);

    for iSlice = [ sliceNo 1 ]
        ind = I.StackIndex == stackID(iSt) & I.StackLocIndex == sliceID{iSt}(iSlice);
        w   = I.Weight( ind );
        tx  = I.TranslationX( ind );
        ty  = I.TranslationY( ind );
        tz  = I.TranslationZ( ind );
        rx  = I.RotationX( ind );
        ry  = I.RotationY( ind );
        rz  = I.RotationZ( ind );
        dx  = I.MeanDisplacementX( ind );
        dy  = I.MeanDisplacementY( ind );
        dz  = I.MeanDisplacementZ( ind );
        disp = I.MeanDisplacement( ind );
        dev = D.devIm(ind);
        worldCoord = cell2mat({D.worldCoord{ind}});
        Aslice2vol = D.Aslice2vol{unique(I.LocIndex(ind))+1};
        Pslice2vol = convert_transformation_matrix_to_parameters( Aslice2vol );
        dist = calc_dist( eye(4), Aslice2vol, worldCoord );
        dispAl = calc_disp( dist );
        dw  = I.WeightedMeanDisplacement( ind );
        tre = I.TRE( ind );
        iCol = iSlice-min(sliceNo)+1;
        t = I.Time( ind ) + t0(iCol);
        iRow = 1;
        xTick = floor(min(t)):xMajorTickIncrement:max(t); 
        % Images
        hAx{iRow,iCol} = subplot(nRow,nSlice,(iRow-1)*nSlice+iCol);
        set_ax_position(hAx{iRow,iCol},iRow,iCol);
        B = bwboundaries( M(iSt).mask(:,:,iSlice), 'noholes' );  % NOTE: 'noholes' option added to avoid crashing in MATLAB2017a
        imshow( M(iSt).im(:,:,iSlice), M(iSt).displaySIRange )
        hold on
        for iB = 1:numel(B)
            line(B{iB}(:,2),B{iB}(:,1),'LineWidth',1,'Color','c')
        end
        hold off
        hAx{iRow,iCol}.XLim = M(iSt).displayXLim + 0.5;
        hAx{iRow,iCol}.YLim = M(iSt).displayYLim + 0.5;
        iRow = iRow + 1;
        % Translation
        hAx{iRow,iCol} = subplot(nRow,nSlice,(iRow-1)*nSlice+iCol);
        set_ax_position(hAx{iRow,iCol},iRow,iCol);
        hold on
        if (showAlParam)
            hAlTx = plot_AlX(t,Pslice2vol.tx);
            hAlTy = plot_AlY(t,Pslice2vol.ty);
            hAlTz = plot_AlZ(t,Pslice2vol.tz);
        end
        hAkTx = plot_AkX(t,tx);
        hAkTy = plot_AkY(t,ty);
        hAkTz = plot_AkZ(t,tz);
        hold off
        hAx{iRow,iCol}.XLim = [min(t),max(t)];
        hAx{iRow,iCol}.XTick = xTick;
        hAx{iRow,iCol}.XTickLabel = '';
        hAx{iRow,iCol}.YLim = limTranslation*[-1,+1];
        hAx{iRow,iCol}.YTick = hAx{iRow,iCol}.YLim(2)*[-1,-0.5,0,0.5,1];
        hAx{iRow,iCol}.XAxis.MinorTickValues = setdiff( floor(min(t)):xMinorTickIncrement:hAx{iRow,iCol}.XLim(2), hAx{iRow,iCol}.XAxis.TickValues );
        hAx{iRow,iCol}.YAxis.MinorTickValues = setdiff( hAx{iRow,iCol}.YLim(2)*(-1:0.25:1), hAx{iRow,iCol}.YAxis.TickValues );
        grid on
        grid minor
        box on
        if iCol == 1
            ylabel( sprintf('Translation\n[mm]'), 'FontSize', 14 )
        else
            hAx{iRow,iCol}.YTickLabel = '';
        end
        if iCol == nSlice-1
            hLgnd{iRow} = legend( [hAkTx,hAkTy,hAkTz], {'\it{tx}','\it{ty}','\it{tz}'}, 'Orientation', 'vertical', 'Location', 'West', 'FontSize', 12 );
            hAx{iRow,iCol+1} = subplot(nRow,nSlice,(iRow-1)*nSlice+iCol+1);
            set_ax_position(hAx{iRow,iCol+1},iRow,iCol+1);
            hLgnd{iRow}.Position(1) = hAx{iRow,iCol+1}.Position(1);
            hAx{iRow,iCol+1}.delete;
        end
        iRow = iRow + 1;
        % Rotation
        hAx{iRow,iCol} = subplot(nRow,nSlice,(iRow-1)*nSlice+iCol);
        set_ax_position(hAx{iRow,iCol},iRow,iCol);
        hold on
        if (showAlParam)
            hAlRx = plot_AlX(t,Pslice2vol.rx);
            hAlRy = plot_AlY(t,Pslice2vol.ry);
            hAlRz = plot_AlZ(t,Pslice2vol.rz);
        end
        hAkRx = plot_AkX(t,rx);
        hAkRy = plot_AkY(t,ry);
        hAkRz = plot_AkZ(t,rz);
        hold off
        hAx{iRow,iCol}.XLim = [min(t),max(t)];
        hAx{iRow,iCol}.XTick = xTick;
        hAx{iRow,iCol}.XTickLabel = '';
        hAx{iRow,iCol}.YLim = limRotation*[-1,+1];
        hAx{iRow,iCol}.YTick = hAx{iRow,iCol}.YLim(2)*[-1,-0.5,0,0.5,1];
        hAx{iRow,iCol}.XAxis.MinorTickValues = setdiff( floor(min(t)):xMinorTickIncrement:hAx{iRow,iCol}.XLim(2), hAx{iRow,iCol}.XAxis.TickValues );
        hAx{iRow,iCol}.YAxis.MinorTickValues = setdiff( hAx{iRow,iCol}.YLim(2)*(-1:0.25:1), hAx{iRow,iCol}.YAxis.TickValues );
        grid on
        grid minor
        box on
        if iCol == 1
            ylabel( sprintf('Rotation\n[degrees]'), 'FontSize', 14 )
        else
            hAx{iRow,iCol}.YTickLabel = '';
        end
        if iCol == nSlice-1
            hLgnd{iRow} = legend( [hAkRx,hAkRy,hAkRz], {'\it{rx}','\it{ry}','\it{rz}'}, 'Orientation', 'vertical', 'Location', 'West', 'FontSize', 12 );
            hAx{iRow,iCol+1} = subplot(nRow,nSlice,(iRow-1)*nSlice+iCol+1);
            set_ax_position(hAx{iRow,iCol+1},iRow,iCol+1);
            hLgnd{iRow}.Position(1) = hAx{iRow,iCol+1}.Position(1);
            hAx{iRow,iCol+1}.delete;
        end
        iRow = iRow + 1;
        % Displacement and Deviation
        hAx{iRow,iCol} = subplot(nRow,nSlice,(iRow-1)*nSlice+iCol);
        set_ax_position(hAx{iRow,iCol},iRow,iCol);
        hold on
        hAlDisp = plot_Al(t,dispAl,hsv2rgb(rgb2hsv(cDisp).^[1,1,0.75]),1.375);
        hAkDisp = plot_Ak(t,disp,cDisp,1.25);
        hAkDev  = plot_Ak(t,dispAl+dev,cDev,1);
        hold off
        hAx{iRow,iCol}.XLim = [min(t),max(t)];
        hAx{iRow,iCol}.XTick = xTick;
        if  ( iRow~=nRow )
            hAx{iRow,iCol}.XTickLabel = '';
        elseif ( iCol == round( nSlice / 2 ) )
            xlabel( sprintf('Acquisition Time [s]'), 'FontSize', 14 )
        end
        hAx{iRow,iCol}.YLim = 2*limDisplacement*[0,+1];
        hAx{iRow,iCol}.YTick = hAx{iRow,iCol}.YLim(2)*[0,0.25,0.5,0.75,1];
        hAx{iRow,iCol}.XAxis.MinorTickValues = setdiff( floor(min(t)):xMinorTickIncrement:hAx{iRow,iCol}.XLim(2), hAx{iRow,iCol}.XAxis.TickValues );
        hAx{iRow,iCol}.YAxis.MinorTickValues = setdiff( hAx{iRow,iCol}.YLim(2)*(-1:0.125:1), hAx{iRow,iCol}.YAxis.TickValues );
        grid on
        grid minor
        box on
        if iCol == 1
            ylabel( sprintf('Displacement\n[mm]'), 'FontSize', 14 )
        else
            hAx{iRow,iCol}.YTickLabel = '';
        end
        if iCol == nSlice-1
            hLgnd{iRow} = legend( [hAkDisp,hAlDisp,hAkDev], {'\it{disp}\rm(A_{\it{k}}\rm)','\it{disp}\rm(A_{\it{l}}\rm)','\it{dev}\rm(A_{\it{k}}\rm)'}, 'Orientation', 'vertical', 'Location', 'West', 'FontSize', 12 );
            hAx{iRow,iCol+1} = subplot(nRow,nSlice,(iRow-1)*nSlice+iCol+1);
            set_ax_position(hAx{iRow,iCol+1},iRow,iCol+1);
            hLgnd{iRow}.Position(1) = hAx{iRow,iCol+1}.Position(1);
            hAx{iRow,iCol+1}.delete;
        end
        iRow = iRow + 1;
        % Frame-Wise Probability
        if showOutrej
            hAx{iRow,iCol} = subplot(nRow,nSlice,(iRow-1)*nSlice+iCol);
            set_ax_position(hAx{iRow,iCol},iRow,iCol);
            hAkPframe = plot_Ak(t,w,cWeight,1.5);
            hAx{iRow,iCol}.XLim = [min(t),max(t)];
            hAx{iRow,iCol}.XTick = xTick;
            if  ( iRow~=nRow )
                hAx{iRow,iCol}.XTickLabel = '';
            elseif ( iCol == round( nSlice / 2 ) )
                xlabel( sprintf('Acquisition Time [s]'), 'FontSize', 14 )
            end
            hAx{iRow,iCol}.YLim = [0,1.02];
            hAx{iRow,iCol}.YTick = hAx{iRow,iCol}.YLim(2)*(0:0.25:1);
            hAx{iRow,iCol}.YTickLabel = {'0','','0.5','','1'};
            hAx{iRow,iCol}.XAxis.MinorTickValues = setdiff( floor(min(t)):xMinorTickIncrement:hAx{iRow,iCol}.XLim(2), hAx{iRow,iCol}.XAxis.TickValues );
            hAx{iRow,iCol}.YAxis.MinorTickValues = setdiff( hAx{iRow,iCol}.YLim(2)*(0:0.125:1), hAx{iRow,iCol}.YAxis.TickValues );
            grid on
            grid minor
            box on
            if iCol == 1
                ylabel( sprintf('Posterior\nProbability'), 'FontSize', 14 )
            else
                hAx{iRow,iCol}.YTickLabel = '';
            end
            if iCol == nSlice-1
                hLgnd{iRow} = legend( {'\it{p^{frame}}'}, 'Orientation', 'vertical', 'Location', 'West', 'FontSize', 12 );
                hAx{iRow,iCol+1} = subplot(nRow,nSlice,(iRow-1)*nSlice+iCol+1);
                set_ax_position(hAx{iRow,iCol+1},iRow,iCol+1);
                hLgnd{iRow}.Position(1) = hAx{iRow,iCol+1}.Position(1);
                hAx{iRow,iCol+1}.delete;
            end
        end
    end

    % Adjust Legend Box Widths
    hLgndWidth = 0;
    for iRow = 1:nRow
        if ~isempty( hLgnd{iRow} )
            hLgndWidth = max( hLgndWidth, hLgnd{iRow}.Position(3) );
        end
    end
    for iRow = 1:nRow
        if ~isempty( hLgnd{iRow} )
            hLgnd{iRow}.Position(3) = hLgndWidth;
        end
    end

end


%% Save Figure

pngDir     = fullfile( resultsDir );
figDir     = fullfile( resultsDir, 'fig' );
epsDir     = fullfile( resultsDir, 'eps' );

save_figs( pngDir, H, figDir, epsDir );


end  % plot_displacement_and_probability_v_time(...)


function P = convert_transformation_matrix_to_parameters( A )
%CONVERT_TRANSFORMATION_MATRIX_TO_PARAMETERS  
%
%   P = CONVERT_TRANSFORMATION_MATRIX_TO_PARAMETERS( A ) converts 
%   transformation matrix A to parameters structure P with field tx, ty,
%   tz, rx, ry and rz.

% jfpva (joshua.vanamerom@kcl.ac.uk)

P.tx = A(1,4);
P.ty = A(2,4);
P.tz = A(3,4);

tol = 0.000001;
ry = asin(-1 * A(1, 3) );

%   // asin returns values for tmp in range -pi/2 to +pi/2, i.e. cos(tmp) >=
%   // 0 so the division by cos(tmp) in the first part of the if clause was
%   // not needed.
if (abs(cos(ry)) > tol)
    P.rx = atan2(A(2,3), A(3,3));
    P.ry = ry;
    P.rz = atan2(A(1,2), A(1,1));
else
    % A(1,3) is close to +1 or -1
    P.rx = atan2(-A(1,3)*A(2,1), -A(1,3)*A(3,1));
    P.ry = ry;
    P.rz = 0;
end

%   // Convert to degrees.
P.rx = P.rx * 180 / pi;
P.ry = P.ry * 180 / pi;
P.rz = P.rz * 180 / pi;

%{

% IRTK method, for reference: 

void irtkRigidTransformation::Matrix2Parameters(irtkMatrix m, double* params) const
{
  double tmp;
  double TOL = 0.000001;

  params[TX] = A(0, 3);
  params[TY] = A(1, 3);
  params[TZ] = A(2, 3);

  tmp = asin(-1 * A(0, 2));

  // asin returns values for tmp in range -pi/2 to +pi/2, i.e. cos(tmp) >=
  // 0 so the division by cos(tmp) in the first part of the if clause was
  // not needed.
  if (fabs(cos(tmp)) > TOL) {
    params[RX] = atan2(A(1,2), A(2,2));
    params[RY] = tmp;
    params[RZ] = atan2(A(0,1), A(0,0));
  } else {
    //A(0,2) is close to +1 or -1
    params[RX] = atan2(-1.0*A(0,2)*A(1,0), -1.0*A(0,2)*A(2,0));
    params[RY] = tmp;
    params[RZ] = 0;
  }

  // Convert to degrees.
  params[RX] *= 180.0/M_PI;
  params[RY] *= 180.0/M_PI;
  params[RZ] *= 180.0/M_PI;

}
%}

end  % convert_transformation_matrix_to_parameters(...)