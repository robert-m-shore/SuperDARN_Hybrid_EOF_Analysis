function [alpha, h, varargout] = geoc2geod(r, theta, varargin)
% [alpha, h]       = geoc2geod(r, theta);
% [alpha, h, X, Z] = geoc2geod(r, theta, B_r, B_theta);
%
% Input:   geographic co-latitude theta (rad)
%          geocentric radius r (km)
%          B_r, B_theta
% Output:  geodetic latitude alpha (rad)
%          geodetic altitude h [km]
%          X, Z

% Nils Olsen, DTU Space, June 2011

% Ellipsoid GRS 80 (identical for WGS84)
a = 6378.137;
b = 6356.752;
rad = pi/180;

RTOD = 57.2957795130823;
DTOR = 0.01745329251994330;

E2 = 1.-(b/a)^2;
E4 = E2*E2;
E6 = E4*E2;
E8 = E4*E4;
OME2REQ = (1.-E2)*a;
A21 =     (512.*E2 + 128.*E4 + 60.*E6 + 35.*E8)/1024.;
A22 =     (                        E6 +     E8)/  32.;
A23 = -3.*(                     4.*E6 +  3.*E8)/ 256.;
A41 =    -(           64.*E4 + 48.*E6 + 35.*E8)/1024.;
A42 =     (            4.*E4 +  2.*E6 +     E8)/  16.;
A43 =                                   15.*E8 / 256.;
A44 =                                      -E8 /  16.;
A61 =  3.*(                     4.*E6 +  5.*E8)/1024.;
A62 = -3.*(                        E6 +     E8)/  32.;
A63 = 35.*(                     4.*E6 +  3.*E8)/ 768.;
A81 =                                   -5.*E8 /2048.;
A82 =                                   64.*E8 /2048.;
A83 =                                 -252.*E8 /2048.;
A84 =                                  320.*E8 /2048.;

GCLAT = 90-theta/rad;
SCL = sin(GCLAT*DTOR);

RI = a./r;
A2 = RI.*(A21+RI.*(A22+RI.* A23));
A4 = RI.*(A41+RI.*(A42+RI.*(A43+RI.*A44)));
A6 = RI.*(A61+RI.*(A62+RI.* A63));
A8 = RI.*(A81+RI.*(A82+RI.*(A83+RI.*A84)));

CCL = sqrt(1-SCL.^2);
S2CL = 2.*SCL.*CCL;
C2CL = 2.*CCL.*CCL-1.;
S4CL = 2.*S2CL.*C2CL;
C4CL = 2.*C2CL.*C2CL-1.;
S8CL = 2.*S4CL.*C4CL;
S6CL = S2CL.*C4CL+C2CL.*S4CL;

DLTCL = S2CL.*A2+S4CL.*A4+S6CL.*A6+S8CL.*A8;
alpha = (DLTCL*RTOD+GCLAT)*rad;
h = r.*cos(DLTCL)-a.*sqrt(1-E2.*sin(alpha).^2);

% convert also magnetic components
psi = sin(alpha).*sin(theta) - cos(alpha).*cos(theta);
if nargin > 2
    psi = sin(alpha).*sin(theta) - cos(alpha).*cos(theta);
    B_r = varargin{1};
    B_theta = varargin{2};
    varargout{1}  = -cos(psi).*B_theta - sin(psi).*B_r; % X
    varargout{2}  = +sin(psi).*B_theta - cos(psi).*B_r; % Z
end

