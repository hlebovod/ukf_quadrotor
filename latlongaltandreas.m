clc;clear all;close all
%% 
%Convert GPS latitude longitude altitude data to cartesian XYZ coordinates

% ==========================
% Import GPS data
% ==========================
load('latlongaltandreas.mat')

gps_position = [slashmatrice210v2slashgpsposition1(:,1) slashmatrice210v2slashgpsposition1(:,2) slashmatrice210v2slashgpsposition1(:,3)];
length(gps_position);

n = size(gps_position,1); % Number of Measurement
freq_gps = 50; % Frequency of gps sensor
sample_gps = 1/freq_gps; % Sample time of gps sensor i.e 0.02sec
t = sample_gps : sample_gps : n*sample_gps;

%%% LLA2ECEF - convert latitude, longitude, and altitude to
%            earth-centered, earth-fixed (ECEF) cartesian
% 
% USAGE:
% [x,y,z] = lla2ecef(lat,lon,alt)
% 
% x = ECEF X-coordinate (m)
% y = ECEF Y-coordinate (m)
% z = ECEF Z-coordinate (m)
% lat = geodetic latitude (radians)
% lon = longitude (radians)
% alt = height above WGS84 ellipsoid (m)

% function [x,y,z]=lla2ecef(lat,lon,alt)
% % WGS84 ellipsoid constants:
% a = 6378137;
% e = 8.1819190842622e-2;
% % intermediate calculation
% % (prime vertical radius of curvature)
% N = a ./ sqrt(1 - e^2 .* sin(lat).^2);
% % results:
% x = (N+alt) .* cos(lat) .* cos(lon);
% y = (N+alt) .* cos(lat) .* sin(lon);
% z = ((1-e^2) .* N + alt) .* sin(lat);
% return

a = 6378137;
e = 8.1819190842622e-2;
% intermediate calculation
% (prime vertical radius of curvature)
N = a ./ sqrt(1 - e^2 .* sin(gps_position(:,2)).^2);
% results:
x = (N+gps_position(:,3)) .* cos(gps_position(:,1)) .* cos(gps_position(:,2));
y = (N+gps_position(:,3)) .* cos(gps_position(:,1)) .* sin(gps_position(:,2));
z = ((1-e^2) .* N + gps_position(:,3)) .* sin(gps_position(:,1));

%save('xyzAndreas','x','y','z')

figure(1)
plot(x,y)
title('Converted Cartesian World Coordinates')
xlabel('x distance (m)')
ylabel('y distance (m)') 
legend('Cartesian XY - GPS')
%hold on

%Way point calculations using lat lon GPS coordinates
x1 = linspace(2.179621925259402e+06,2.179075575899777e+06,50);
y1 = linspace(-4.845267808250103e+06,-4.845173441494693e+06,50);

x2 = linspace(2.179075575899777e+06,2.179505926420107e+06,50);
y2 = linspace(-4.845173441494693e+06,-4.844714701934233e+06,50);

x3 = linspace(2.179505926420107e+06,2.180040161138976e+06,50);
y3 = linspace(-4.844714701934233e+06,-4.844818365425290e+06,50);

x4 = linspace(2.180040161138976e+06,2.179621925259402e+06,50);
y4 = linspace(-4.844818365425290e+06,-4.845267808250103e+06,50);

%Waypoints corners ECEF converted as ground truth
figure(2)
plot(x1, y1,'b'); hold on; grid on;
xlabel('X coordinate(m)');ylabel('Y coordinate(m)');
plot(x2, y2, 'b'); hold on; grid on;
plot(x3, y3, 'b'); hold on; grid on;
plot(x4, y4, 'b'); hold off; grid on;
legend('Cartesian XY - GPS')


%GPS Position
figure(3)
subplot(3,1,1)
plot(t,gps_position(:,1));xlabel('Time (sec)');ylabel('Latitude (rad)');grid on
title('Geodetic Coordinates - Latitude, Longitude and Altitute')
legend('Lat - GPS')
subplot(3,1,2)
plot(t,gps_position(:,2));xlabel('Time (sec)');ylabel('Longitude (rad))');grid on
legend('Lon - GPS')
subplot(3,1,3)
plot(t,gps_position(:,3));xlabel('Time (sec)');ylabel('Altitude (m)');grid on
legend('Alt - GPS')

figure(4)
subplot(3,1,1)
plot(t,x);xlabel('Time (sec)');ylabel('Position X (m)');grid on
title('Converted to Cartesian XYZ World Coordinates')
legend('Cartesian X - GPS')
subplot(3,1,2)
plot(t,y);xlabel('Time (sec)');ylabel('Position Y (m)');grid on
legend('Cartesian Y - GPS')
subplot(3,1,3)
plot(t,z);xlabel('Time (sec)');ylabel('Position Z (m)');grid on
legend('Cartesian Z - GPS')

figure(5)
plot3(x,y,z,'r','linewidth',1.25)
title('Converted to Cartesian XYZ World Coordinates')
xlabel('X(m)');ylabel('Y(m)');zlabel('Z(m)')
legend('Cartesian XYZ - GPS')
grid on;

