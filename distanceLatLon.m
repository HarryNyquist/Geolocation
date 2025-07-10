% Haversine Distance Calculation in MATLAB
% This script computes the great-circle distance between two latitude/longitude points

clc; clear; close all;

% Earth's radius in meters
h = 500e3;
R = 6371e3; 

% Define latitude and longitude for two locations (in degrees)
lat1 = 48.0; % New York City
lon1 = 11.0;

lat2 = 49.0; % Los Angeles
lon2 = 11.0;

% Convert degrees to radians
lat1 = deg2rad(lat1);
lat2 = deg2rad(lat2);
lon1 = deg2rad(lon1);
lon2 = deg2rad(lon2);
del_lat = (lat2 - lat1);
del_lon = (lon2 - lon1);

% Haversine formula
a = sin(del_lat/2) * sin(del_lat/2) + cos(lat1) * cos(lat2) * sin(del_lon/2) * sin(del_lon/2);
c = 2 * atan2(sqrt(a), sqrt(1 - a));

% Compute the distance
d = (R+h) * c; % Distance in meters

% Display the result
fprintf('Distance between (%.4f, %.4f) and (%.4f, %.4f):\n', rad2deg(lat1), rad2deg(lon1), rad2deg(lat2), rad2deg(lon2));
fprintf('%.2f meters (%.2f km)\n', d, d/1000);


% lat1 = 48; lat2 = 49; lon1 = 11 ; lon2 = 12;
% h = 500e3 ;
% R = 6371e3; % Earth's radius in meters
% rad = R+h ;
% distance = dist_p(rad,lat1, lon1, lat2, lon2)
% 
% 
% function dist_p = distance_LatLon(rad,lat1, lon1, lat2, lon2)
% % Convert degrees to radians
%     lat1 = deg2rad(lat1);
%     lat2 = deg2rad(lat2);
%     delta_lat = (lat2 - lat1);
%     delta_lon= deg2rad(lon2 - lon1);
% 
%     % Haversine formula
%     a = sin(delta_lat/2) * sin(delta_lat/2) + cos(lat1) * cos(lat2) * sin(delta_lon/2) * sin(delta_lon/2);
%     c = 2 * atan2(sqrt(a), sqrt(1 - a));
% 
%     % Compute the distance
%     dist_p = rad * c; % Distance in meters
%    
% end
% 
% 
