%% Stress2Grid v1.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These scripts/data are freely available under a Creative Commons Attribution 4.0     % 
% International (CC-BY 4.0) Licence. When using the scripts/data please cite them as:  %
%                                                                                      %
% Ziegler, Moritz; Heidbach, Oliver (2017): Matlab script Stress2Grid. GFZ German      % 
% Research Center for Geosciences. http://doi.org/10.5880/wsm.2017.002                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all
addpath('routines') % The path to the auxiliary files
addpath('data') % Path to the data files
% Input variables
% In this section all input variables which are required can be altered.

% Input files and settings
%input_xls = 'wsm2016.xlsx'; % The WSM database (*.xlsx - File)
input_csv= 'wsm2016.csv'; % The WSM database (*.csv -File)

plate_boundaries = 'plates_Bird_2002_stress2grid.dat'; % Provide the file of the plate boundaries
euler_poles = 'HS3_NUVEL1A.txt'; % Provide the file for the Euler poles for the different plates
num_plates = 15; % The number of plates with Euler poles

polygon_exclude = 'yes'; % If you want to exclude the evaluation for a specific polygon.
exclude_poly = 'PB2002_orogens.txt'; % Provide a file with a polygon from which you want to exclude all data.

% Grid parameters
gridsize = 2;             % gridsize in degree
west = 0;                 % minimum and Maximum Lat and Lon
east = 20;
south = 40; 
north = 54; 

% Data processing
% WSM-Quality weighting
%         A      B     C    D  E
qual = [ 1/15, 1/20, 1/25, 1/40, 0 ];
apply_qw = 'yes'; % If a quality weighting should be applied: 'yes'

% WSM method weighting (from 0 to 5)
%       FMS FMF BO DIF HF GF GV OC, rest
methd = [ 4, 5, 5, 5,  4, 5, 4, 2, 1 ];
apply_mw = 'yes'; % If a method weighting should be applied: 'yes'
m_exclude = {}; % Certain stress indicators can be completely removed from the analysis.

pbe_exclude = 'yes'; % Exclude PBE-flagged events from the algorithm.
pb_dist_exclude = 0; % Exclude data records in a given distance from the next plate boundary.

dist_weight = 'linear';    % method of distance weighting (linear, inverse, or none)
dist_threshold = 0.1;       % distance weight to prevent overweight of data nearby (0 to 1)

% Initialise script functions
R_range = [50:50:1000];
min_data = 3;                % minimum number of data per bin
threshold = 25;             % threshold for deviation of orientation
plate_affil = 'yes';   % Only consider data from the same plate as the gridpoint
arte_thres = 200; % Maximum distance (in km) of the gridpoint to the next datapoint.
compare_pm = 'yes'; % Set to 'yes' if plate motion should be compared.

% Output
output = 'both'; % Define output type: excel, gmt or both.
arrowlength = 0.1 ; % Define arrowlength for GMT.

% Plots:
plot_output = 'yes'; % Specify whether the results should be plotted immediately

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% This is the end of the user input %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-Processing
% In the following the data is prepared, the grid is generated
% and several variables are assigned.
% Several fixed variables are assigned and data is loaded:
deg2rad = pi/180;       % degree to rad
R_earth = 6371;         % radius of the earth

num_r = length(R_range);

h = waitbar(0,'Loading and rating input data...');

exls = exist('input_xls','var');
ecsv = exist('input_csv','var');

% Load stress data according to specified input
if exls > 0
    [~, ~, alldata] = xlsread(input_xls);
    alldata(1,:) = [];

    lat_data0 = round(cell2mat(alldata(:,3)),3);
    lon_data0 = round(cell2mat(alldata(:,4)),3);
    orient_data0 = round(cell2mat(alldata(:,5)));
    method = alldata(:,6);
    depth_data = round(cell2mat(alldata(:,7)),3);
    quality = alldata(:,8);
    regime = alldata(:,9);
    plate = alldata(:,55);
    pbe = alldata(:,58);
    pb_dist = round(cell2mat(alldata(:,57)),2);
    
elseif ecsv > 0
    fileID = fopen(input_csv);
    wsmformat = ('%*s %*s %.6f %.6f %f %s %.2f %s %s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %s %*s %f %s %*[^\n]');
    alldata = textscan(fileID,wsmformat,'HeaderLines',1,'Delimiter',',');
    fclose(fileID);

    lat_data0 = round(alldata{:,1},3);
    lon_data0 = round(alldata{:,2},3);
    orient_data0 = round(alldata{:,3});
    method = alldata{:,4};
    depth_data = round(alldata{:,5},3);
    quality = alldata{:,6};
    regime = alldata{:,7};
    plate = alldata{:,8};
    pb_dist = round(alldata{:,9},2);
    pbe = alldata{:,10};
    
else
    disp(' ')
    disp('No input file specified.')
    return
end

% Check input data
[ec, se] = check_input(lat_data0,lon_data0,orient_data0,method,depth_data,quality,plate,...
    pb_dist,pbe,regime);

disp(' ')
if ec > 0
    disp('Fatal error in the input file. Check messages.')
    disp(' ')
    close(h)
    return
elseif se > 0
    disp('Some non fatal errors were detected in the input file. The script resumes. It is recmommended to check the messages and alter the input accordingly.')
    disp(' ')
else
    disp('Input file checked')
    disp(' ')
end
clear ec se

% Exclude certain stress indicators
if length(m_exclude) > 0
    ind_bad = 0;
    for i = 1:length(m_exclude)
        for j = 1:length(lat_data0)
            if isequal(m_exclude{i},method{j})
                ind_bad = [ind_bad, j];
            end
        end
    end
    ind_bad(1) = [];
    
    lat_data0(ind_bad) = [];
    lon_data0(ind_bad) = [];
    orient_data0(ind_bad) = [];
    method(ind_bad) = [];
    depth_data(ind_bad) = [];
    quality(ind_bad) = [];
    regime(ind_bad) = [];
    plate(ind_bad) = [];
    pbe(ind_bad) = [];
    pb_dist(ind_bad) = [];
    
    disp([num2str(length(ind_bad)),' stress data records removed'])
    disp(' ')
    
    clear s ind_bad
end

% The data is sorted and weighted
if isequal(apply_qw,'yes')
    qw0 = qualityweight(quality,qual(:,1),qual(:,2),qual(:,3),qual(:,4),qual(:,5));
else
    qw0 = ones(length(lat_data0),1);
end

if isequal(apply_mw,'yes')
    ind_method0 = methodweight(method,methd(:,1),methd(:,2),methd(:,3),methd(:,4),methd(:,5),methd(:,6),methd(:,7),methd(:,8),methd(:,9));
    ind_method0 = ind_method0 / 5;
else
    ind_method0 = ones(length(lat_data0),1);
end

datas0_ = [lon_data0,lat_data0,sin(2*orient_data0*deg2rad),cos(2*orient_data0*deg2rad),ind_method0,orient_data0,qw0,pb_dist];

% Exclude deselected data
ind_bad = find(datas0_(:,5) == 0);
if isempty(ind_bad) ~= 1      % deselected qualities
    datas0_(ind_bad,:) = [];
    plate(ind_bad,:) = [];
    regime(ind_bad,:) = [];
    pbe(ind_bad,:) = [];
    orient_data0(ind_bad) = [];
end

ind_bad = find(datas0_(:,7) == 0);
if isempty(ind_bad) ~= 1      %% deselected methods
    datas0_(ind_bad,:) = [];
    plate(ind_bad,:) = [];
    regime(ind_bad,:) = [];
    pbe(ind_bad,:) = [];
    orient_data0(ind_bad) = [];
end

if isequal(pbe_exclude,'yes')
    ind_bad = find(ismember(pbe,'PBE'));
    if isempty(ind_bad) ~= 1      % deselected events
        datas0_(ind_bad,:) = [];
        plate(ind_bad,:) = [];
        regime(ind_bad,:) = [];
        pbe(ind_bad,:) = [];
        orient_data0(ind_bad) = [];
    end
end

ind_bad = find(orient_data0 == 999);
if isempty(ind_bad) ~= 1      % Events with SHmax = 999
        datas0_(ind_bad,:) = [];
        plate(ind_bad,:) = [];
        regime(ind_bad,:) = [];
end

ind_bad = find(datas0_(:,8) < pb_dist_exclude);
if isempty(ind_bad) ~= 1      % exclude events close to plate boundaries
        datas0_(ind_bad,:) = [];
        plate(ind_bad,:) = [];
        regime(ind_bad,:) = [];
end

close(h);
%
init_1 = datas0_(:,1);
init_2 = datas0_(:,2);
init_orient = datas0_(:,6);

datas_ = datas0_;

disp('Input data processed')
disp(' ')

clear input alldata Qual methd plot_input...
    depth_data ind_bad ind_method0 lat_data0 lon_data0 method...
    orient_data0 quality qw0 regime datas0_ h area a c h numdat...
    dens apply_mw apply_qw pbe_exclude pb_dist_exclude pb_dist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regular grid
[XG,YG] = meshgrid(west:gridsize:east,south:gridsize:north);
XG = reshape(XG,[],1);
YG = reshape(YG,[],1);

if isequal(polygon_exclude,'yes')
    % Read the polygon data from specified file
    fid = fopen(exclude_poly);
    X = 0;
    Y = 0;
    ExPoly_num = 0;
    
    tline = fgetl(fid);
    while ischar(tline)
        if isletter(tline)
            ExPoly_num = ExPoly_num + 1;
        elseif isequal(tline,'*** end of line segment ***');
            X(1) = [];
            Y(1) = [];
            assignin ('base',strcat('ExPoly_',num2str(ExPoly_num)), [X Y]);
            X = 0;
            Y = 0;
        else
            L = sscanf(tline,' %f,%f');
            X = [ X; L(1) ];
            Y = [ Y; L(2) ];
        end
        tline = fgetl(fid);
    end
    fclose(fid);

    % Exclude the gridpoints from within polygons 
    h = waitbar(0,'Excluding the gridpoints from within polygons...');

    for i = 1:ExPoly_num
        Poly = evalin('base',strcat('ExPoly_',num2str(i)));
        % Exclude grid points
        in = inpolygon(XG,YG,Poly(:,1),Poly(:,2));
        k = find(in == 1);
        XG(k) = [];
        YG(k) = [];
        
        % Exclude data records from within polygons
        in = inpolygon(datas_(:,1),datas_(:,2),Poly(:,1),Poly(:,2));
        k = find(in == 1);
        datas_(k,:) = [];
        plate(k,:) = [];        
    
        waitbar(i/ExPoly_num)
    end
    close(h);
    
    disp('Polygons excluded')
    disp(' ')
end

n_G = length(XG);
Plate_G = cell(n_G,1);

% Plate motion/strain preparation:
if isequal(compare_pm,'yes') || isequal(plate_affil,'yes')

% Load plate polygons
fid = fopen(plate_boundaries);
X = 0;
Y = 0;
previous = 'test';
name = 'Na';

tline = fgetl(fid);
while ischar(tline)
    if isletter(tline)
        previous = sscanf(tline,'%c');
        name = vertcat(name,previous);               
    elseif isequal(tline,'*** end of line segment ***');
        X(1) = [];
        Y(1) = [];
        assignin ('base',strcat('PN_',previous), [X Y]);
        X = 0;
        Y = 0;
    else
        L = sscanf(tline,' %f,%f');
        X = [ X; L(1) ];
        Y = [ Y; L(2) ];
    end
    
    tline = fgetl(fid);
end
fclose(fid);
name(1,:) = [];

% Assign plate names to gridpoints
h = waitbar(0,'Assigning plate names to gridpoints...');

for i = 1:length(name)
    Poly = evalin('base',strcat('PN_',name(i,:)));    
    in = inpolygon(XG,YG,Poly(:,1),Poly(:,2));
    k = find(in == 1);
    
    if isequal('NY',name(i,:))
        Plate_G(k,:) = cellstr('NA');
    elseif isequal('PB',name(i,:))
        Plate_G(k,:) = cellstr('PA');
    elseif isequal('PC',name(i,:))
        Plate_G(k,:) = cellstr('PA');
    elseif isequal('KF',name(i,:))
        Plate_G(k,:) = cellstr('KE');
    elseif isequal('BW',name(i,:))
        Plate_G(k,:) = cellstr('BR');
    elseif isequal('AV',name(i,:))
        Plate_G(k,:) = cellstr('AU');
    else
        Plate_G(k,:) = cellstr(name(i,:));
    end
    
    waitbar(i/length(name))
end
close(h);

disp('Tectonic plates initialised')
disp(' ')
end

clear previous fid X Y tline L PN_* Poly i j h name in k... 
    plate_boundaries ExPoly_* exclude_poly polygon_exclude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of Plate Motion

if isequal(compare_pm,'yes')
% Transform the coordinates of the gridpoints
blong = XG * deg2rad; 
blat = (90-YG) * deg2rad; 
xcord = (R_earth * sin(blat)) .* cos(blong);
ycord = (R_earth * sin(blat)) .* sin(blong);
zcord = R_earth * cos(blat);

% Get the Euler pole rotations
fid = fopen(euler_poles);
for i = 1:num_plates
    tline = fgetl(fid);
    ep(i,:) = textscan(tline,'%s %f %f %f %*500c','Delimiter',' ','MultipleDelimsAsOne',1);
end
fclose(fid);
euler_plate = [ep{:,1}]';

% Compute the Euler Pole omega vector   
b_lat  =  (90 - [ep{:,2}]') * deg2rad;
b_long =  [ep{:,3}]' * deg2rad;
b_rot  =  [ep{:,4}]' * deg2rad;
    
w1 = b_rot .* sin(b_lat) .* cos(b_long);
w2 = b_rot .* sin(b_lat) .* sin(b_long);
w3 = b_rot .* cos(b_lat);


% Derive absolute plate motion for each gridpoint
apm = zeros(1,n_G);

h = waitbar(0,'Compute absolute plate motion at gridpoints...');

for i = 1:n_G
k = find( strcmp(Plate_G(i),euler_plate));

if k > 0
    % Compute the rotation vector at the gridpoints
    vx = w2(k) * zcord(i) - w3(k) * ycord(i);
    vy = w3(k) * xcord(i) - w1(k) * zcord(i);
    vz = w1(k) * ycord(i) - w2(k) * xcord(i);

    dnord = 0.0001 * (-vx * cos(blat(i)) * cos(blong(i)) - vy * cos(blat(i)) * sin(blong(i)) + vz * sin(blat(i)));
    deast = 0.0001 * (-vx * sin(blong(i)) + vy * cos(blong(i)));
    dver  = 0.0001 * ( vx * sin(blat(i)) * cos(blong(i)) + vy * sin(blat(i)) * sin(blong(i)) + vz * cos(blat(i)));
    dver  = abs(dver);

    leng = sqrt(dnord^2 + deast^2);
    
    if (dnord^2 + deast^2) == 0
        apm(i) = 0;    
    elseif (dnord >= 0 && deast >= 0)
        apm(i) = (1/deg2rad) * asin(deast/leng);
    elseif (dnord <  0 && deast >= 0)
        apm(i) = 90 + abs((1/deg2rad) * asin(dnord/leng));
    elseif (dnord <  0 && deast < 0)
        apm(i) = abs((1/deg2rad) * asin(deast/leng));
    elseif (dnord >= 0 && deast < 0)
        apm(i) = 180.0 - abs((1/deg2rad) * asin(deast/leng));
    end

else
    apm(i) = -999;
end

waitbar(i/n_G)
end
close(h)

disp('Plate motion computed')
disp(' ')

end

clear blat blong *cord b_* w1 w2 w3 dver leng dnord deast vx vy...
    vz k i fid ep euler_poles R_earth euler_plate tline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the mean SHmax orientation
% Initialise the required variables
SH = cell(n_G,(num_r*6));
gmt_coord = cell(n_G,4);
gmt_mazi = cell(n_G,3);
gmt_shpm = cell(n_G,1);

% This is the computation of the mean values for different search radii
h = waitbar(0,'Computing Mean SH...');
for i = 1:n_G
    distij = ddistance(YG(i),XG(i),datas_(:,2),datas_(:,1));
    
    if arte_thres > 0 && min(distij) > arte_thres
                continue;
    end    
    
    for k = 1:length(R_range)
        R_search = R_range(k);        
        ids_R = find(distij < R_search);
        
        % Check if the datapoints are on the same plate as the gridpoint
        if isequal(plate_affil,'yes')
            cp = find( strcmp(Plate_G(i),plate(ids_R)) - 1 );
            ids_R(cp) = [];
        end
        
        N_in_R = length(ids_R);        
        sq_R = sum(datas_(ids_R,7));

        if N_in_R < min_data        % not enough data within search radius
            continue;
        elseif N_in_R == 1
            s0 = 0;
            meanSH = datas_(ids_R,6);
            mdr = distij(ids_R)/R_search;
        else
            mdr = mean(distij(ids_R))/R_search;
            
            dist_threshold_scal = R_search*dist_threshold;            
            if isequal(dist_weight,'inverse')
                wd_R = 1./max(dist_threshold_scal,distij(ids_R));
            elseif isequal(dist_weight,'linear')
                wd_R = R_search + 1 - max(dist_threshold_scal,distij(ids_R));
            else
                wd_R = ones(length(ids_R),1);
            end
        
            sw_R = sum(wd_R.*datas_(ids_R,7).*datas_(ids_R,5));
            array_sin2 = wd_R.*datas_(ids_R,7).*datas_(ids_R,5).*datas_(ids_R,3);
            array_cos2 = wd_R.*datas_(ids_R,7).*datas_(ids_R,5).*datas_(ids_R,4);
        
            % mean value
            sumsin2 = sum(array_sin2);
            sumcos2 = sum(array_cos2);
            meansin2 = sumsin2/sw_R;
            meancos2 = sumcos2/sw_R;
            meanR = sqrt(meansin2^2+meancos2^2);
                    
            s0 = sqrt(-2*log(meanR))/(2 * deg2rad);
            meanSH_rad = atan2(meansin2,meancos2)/2;
            meanSH = meanSH_rad/deg2rad;
        
            if meanSH < 0
                meanSH = meanSH + 180;
            end
        end
        
        % Excel output:
        ml = ((k-1) * 6);
        if isequal(compare_pm,'yes')
            if apm(i) ~= -999 && s0 < threshold
                spm = meanSH - apm(i);
                SH(i,ml+2) = num2cell(spm);
            end
        end
        
        SH(i,(ml+1)) = num2cell(meanSH);
        SH(i,(ml+3)) = num2cell(s0);
        SH(i,(ml+4)) = num2cell(R_search);
        SH(i,(ml+5)) = num2cell(mdr);
        SH(i,(ml+6)) = num2cell(N_in_R);
        
        % GMT Output
        if s0 < threshold % get homogeneous stress
            gmt_coord(i,:) = num2cell([ XG(i), YG(i), s0, R_search ]);
            gmt_mazi(i,:) = num2cell([ XG(i), YG(i), meanSH ]);
            if isequal(compare_pm,'yes')
                if apm(i) ~= -999
                    gmt_shpm(i) = num2cell(meanSH - apm(i));
                end
            end

        elseif R_search == min(R_range)
            gmt_coord(i,:) = num2cell([ XG(i), YG(i), s0, R_search ]);
        end
    end
    
    waitbar(i/n_G)
end
close(h)
disp('Mean SHmax computed')
disp(' ')
clear dist_threshold_scal i sum* array_* meansin* meancos*...
    N_in_R ids_R sq_R mdr cos_thetaij sum* ml k meanR distij h...
    wd_R mn_data dist_threshold num_plates cp datas_ meanSH meanSH_rad...
    n_G pbe s0 spm sw_R R_search gridsize no_artis arti_thres...
    threshold plate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot output data
if isequal(plot_output,'yes')

load coastlines;
size = 50; % Size of the coloured areas

YA = cell2mat(gmt_mazi(:,2));
XA = cell2mat(gmt_mazi(:,1));
azi = cell2mat(gmt_mazi(:,3));

Y = cell2mat(gmt_coord(:,2));
X = cell2mat(gmt_coord(:,1));
rad = cell2mat(gmt_coord(:,4));
var = cell2mat(gmt_coord(:,3));

if isequal(compare_pm,'yes')
    apmn = apm;
    XP = XG;
    YP = YG;
    k = find(apmn == -999);
    apmn(k) = [];
    XP(k) = [];
    YP(k) = [];
    
    k = find(cellfun('isempty',gmt_shpm));
    XS = XG;
    YS = YG;
    XS(k) = [];
    YS(k) = [];
    shpm = cell2mat(gmt_shpm);
    
end


subplot(2,3,1)
    plot(coastlon,coastlat,'k');
    hold on
    quiver(init_1,init_2,sin(init_orient*deg2rad),cos(init_orient*deg2rad),'k','ShowArrowHead','off');
    axis([west east south north]);
    title('Input data')
    hold off

subplot(2,3,2)
    hold on
    quiver(XA,YA,sin(azi*deg2rad),cos(azi*deg2rad),'k','ShowArrowHead','off');
    plot(coastlon,coastlat,'k');
    axis([west east south north]);
    title('Mean S_{Hmax} Azimuth')
    hold off

if isequal(compare_pm,'yes')
subplot(2,3,3)
    hold on
    quiver(XP,YP,sin(apmn'*deg2rad),cos(apmn'*deg2rad),'k','ShowArrowHead','off');
    plot(coastlon,coastlat,'k');
    axis([west east south north]);
    title('Absolut plate motion')
    hold off
end

subplot(2,3,4);
    hold on
    scatter(X,Y,size,rad,'s','filled');
    cb1 = colorbar('northoutside');
    xlabel(cb1,'Search radius [km]');
    plot(coastlon,coastlat,'k');
    axis([west east south north]);
    hold off

subplot(2,3,5);
    hold on
    scatter(X,Y,size,var,'s','filled');
    cb2 = colorbar('northoutside');
    xlabel(cb2,'Standard Deviation [Degree]');
    plot(coastlon,coastlat,'k');
    axis([west east south north]);
    hold off

if isequal(compare_pm,'yes')
subplot(2,3,6);
    hold on
    scatter(XS,YS,size,shpm,'s','filled');
    cb3 = colorbar('northoutside');
    xlabel(cb3,'mean S_{Hmax} - Plate motion [Degree]');
    plot(coastlon,coastlat,'k');
    axis([west east south north]);
    hold off
end

end

clear XA YA X Y XP YP XS YS azi var rad shpm cb...
    size k r apmn lat long deg2rad plot_* south west north east...
    init_*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export data
% This section is dedicated to the nameing, labeling, and export of the data.
% Excel Output
if isequal(output,'both') || isequal(output,'excel')
    outputfile = num2cell([XG,YG]);
    row1 = cat(2,cellstr(['LAT']),cellstr(['LON']));
    
    if isequal(compare_pm,'yes') || isequal(plate_affil,'yes')
        outputfile = cat(2,outputfile,Plate_G);
        row1 = cat(2,row1,cellstr(['Plate']));
    end
 
    if isequal(compare_pm,'yes')
        outputfile = cat(2,outputfile,num2cell(apm'));
        row1 = cat(2,row1,cellstr(['APM']));
        
        k = find(eq(cell2mat(outputfile(:,4)),-999));
        outputfile(k,4) = {[]};
    end

    outputfile = cat(2,outputfile,SH);
    
    ssnf = cat(2,cellstr(['meanSH ',num2str(R_range(1))]),cellstr(['SH-APM']),cellstr(['SD']),cellstr(['R']),cellstr(['mdr']),cellstr(['N']));
    for i = 1:(num_r-1)
        ssn = cat(2,cellstr(['meanSH ',num2str(R_range(i+1))]),cellstr(['SH-APM']),cellstr(['SD']),cellstr(['R']),cellstr(['mdr']),cellstr(['N']));
        ssnf = cat(2,ssnf,ssn);
    end
    row1 = cat(2,row1,ssnf);    

    outputfile = cat(1,row1,outputfile);
    
    xlswrite('wsm_script_output.xls',outputfile);
    
    clear outputfile row1 k ssnf ssn i
end

clear R_range num_r pb_affil

% Casmi/GMT Output
% Export the data for a GMT based plotting
if isequal(output,'both') || isequal(output,'gmt')
    % Mean Azimuth:
    len = arrowlength * ones(length(cell2mat(gmt_mazi(:,1))),1);
    mazi = [ cell2mat(gmt_mazi(:,1)), cell2mat(gmt_mazi(:,2)), cell2mat(gmt_mazi(:,3)), len ];
    
    fid = fopen('mean_azi.dat','w');
    for i = 1:length(mazi(:,1))
        fprintf(fid,'%2.3f %3.3f %3.2f %1.1f\n',mazi(i,:));
    end
    fclose(fid);
    
    
    % Wavelength/Grid results:
    grire = [ cell2mat(gmt_coord(:,1)), cell2mat(gmt_coord(:,2)), cell2mat(gmt_coord(:,4)) ];
    
    fid = fopen('wavelength.dat','w');
    for i = 1:length(grire(:,1))
        fprintf(fid,'%2.3f %3.3f %4d\n',grire(i,:));
    end
    fclose(fid);
    
    % Standard Deviation:
    std = [ cell2mat(gmt_coord(:,1)), cell2mat(gmt_coord(:,2)), cell2mat(gmt_coord(:,3)) ];
    
    fid = fopen('std.dat','w');
    for i = 1:length(std(:,1))
        fprintf(fid,'%2.3f %3.3f %4d\n',std(i,:));
    end
    fclose(fid);
    
    
    if isequal(compare_pm,'yes')
        % APM:
        Abspm = [ XG, YG, apm' ];
        plate_abs = char(Plate_G);
        k = find(Abspm(:,3) == -999);
        Abspm(k,:) = [];
        plate_abs(k,:) = [];
        len = arrowlength * ones(length(Abspm(:,1)),1);
        Abspm = [ Abspm, len ];
        
        fid = fopen('apm.dat','w');
        for i = 1:length(len)
            fprintf(fid,'%2.3f %3.3f %3.2f %1.1f %s\n',Abspm(i,:),plate_abs(i,:));
        end
        fclose(fid);
        
        % SHmax - APM:
        grid = [ XG, YG ];
        plate_shab = char(Plate_G);
        k = find(cellfun(@isempty,gmt_shpm));
        grid(k,:) = [];
        plate_shab(k,:) = [];
        len = arrowlength * ones(length(grid(:,1)),1);
        shab = [ grid, cell2mat(gmt_shpm(:,1)), len ];
        
        fid = fopen('SHmax-APM.dat','w');
        for i = 1:length(len)
            fprintf(fid,'%2.3f %3.3f %3.2f %1.1f %s\n',shab(i,:),plate_shab(i,:));
        end
        fclose(fid);
        
    end
      
end

disp('Output written to files')
disp(' ')

clear i k mazi grire Abspm plate_abs grid plate_shab shab len fid
clear all