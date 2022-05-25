function [QD_lat, QD_lon, Apex_lat, MLT, F1, F2] = qdipole2(t, r, theta, phi, varargin)
% [QD_lat, QD_lon, Apex_lat, MLT, F1, F2] = qdipole(t, r, theta, phi);
% [QD_lat, QD_lon, Apex_lat, MLT, F1, F2] = qdipole(t, r, theta, phi, filename_grid);
%
% Input:
%    t:         MD2000 (days)
%    r:         radius (in units of a = 6371.2 km)
%    theta:     geographic co-latitude (radians)
%    phi:       geographic longitude   (radians)
%    filename_grid: QD coefficient file [apexsh_1995-2015.txt]
% Output:
%    QD_lat:    Quasi-dipole latitude (in degrees)
%    QD_lon:    Quais-dipole longitude (in degrees)
%    Apex_lat:  (modified) Apex latitude (in degrees)
%    MLT:       Magnetic Local Time (in hours)
%    F1, F2:    Base vectors F1 and F2
%
% Nils Olsen, DTU Space, August 2010
%

rad = pi/180;
a = 6371.2;
R_E = 6371.009; % mean radius according to WGS84
h_R = 110; % reference altitude for (modified) Apex latitude

% determine size of input arrays
max_size = max([size(t); size(r); size(theta); size(phi)]);
max_length = max_size(1)*max_size(2);
% convert to matrix if input parameter is scalar
if length(t)     == 1; t = t*ones(max_size); end;
if length(r)     == 1; r = r*ones(max_size); end;
if length(theta) == 1; theta = theta*ones(max_size); end;
if length(phi)   == 1; phi = phi*ones(max_size); end;
% check for equal length of all input parameters
if size(t) ~= size(r) | size(t) ~= size(theta) | size(t) ~= size(phi);
    error('Variables must be of equal size (or scalars)');
    return
end
% convert to row vector
t = reshape(t, max_length, 1);
r = reshape(r, max_length, 1);
theta = reshape(theta, max_length, 1);
phi = reshape(phi, max_length, 1);

if nargin == 4
    filename =  'c:/nio/m/tools/apexsh_1995-2015.txt'; % File with Apex coefficients
else
    filename = varargin{1};
end

if exist(filename, 'file')
    
    % check if time is in prioper interval
    fid = fopen(filename, 'r');
    C = textscan(fid, '%d %d %d %d %d', 1);
    epochs_grid = textscan(fid, '%f', C{1});
    fclose(fid);
    
    if min(t)/365.25+2000 < min(epochs_grid{1}); error(['epoch = ' num2str(min(t)/365.25+2000) ' earlier than epochs in file ' filename]); end;
    if max(t)/365.25+2000 > max(epochs_grid{1}); error(['epoch = ' num2str(max(t)/365.25+2000) ' later than epochs in file ' filename]); end;
    
    gd_lon = phi/rad;
    
    [QD_lat, QD_lon, MLT, F11, F12, F21, F22] = qdipole2_m(t(:), r(:)*a, 90-theta(:)/rad, gd_lon(:), filename);
    
    if nargout > 2 % also compute Apex latitude
        % (modified) Apex latitude according to eq. (15) of Emmert et al 2010
        [~, h] = geoc2geod(r*a, theta);
        Apex_lat = acos(sqrt((R_E+h_R)./(R_E+h)).*cos(QD_lat*rad))/rad;
        % same sign as QD_lat
        Apex_lat = abs(Apex_lat).*sign(QD_lat);
        Apex_lat= reshape(Apex_lat, max_size);
        MLT= reshape(MLT, max_size);
    end
    
    % convert output variables to matrices of same size as input variables
    QD_lat   = reshape(QD_lat, max_size);
    QD_lon   = reshape(QD_lon, max_size);
    F1      = [reshape(F11, max_length, 1) reshape(F12, max_length, 1)];
    F2      = [reshape(F21, max_length, 1) reshape(F22, max_length, 1)];
else
    error(['File ', filename, ' not found']);
end

