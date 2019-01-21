function [ dist ] = ddistance(lat_g,lon_g,lat,lon)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%lat = [ 50; 30];
%lon = [ 40; 32];
%lon_g = 12;
%lat_g = 20;

r = 6371;               % earth radius
A = pi/180;

phi1 = A * lon_g;          % azimuthal angle
phi2 = A .* lon;
theta1 = A *(90-lat_g);     % polar angle / colatitude
theta2 = A .*(90-lat);

x1 = r .* sin(theta1) .* cos(phi1);
y1 = r .* sin(theta1) .* sin(phi1);
z1 = r .* cos(theta1);
x2 = r .* sin(theta2) .* cos(phi2);
y2 = r .* sin(theta2) .* sin(phi2);
z2 = r .* cos(theta2);

dist = sqrt( (x1-x2).^2 + (y1-y2).^2 + (z1-z2).^2 );

end

