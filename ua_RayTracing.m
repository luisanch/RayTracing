clear variables;
clc
close all;

%% Inputs
h = -3500; %height to sea bottom
grad_inv_h = 0; %height at which speed gradient will be inverted
source = -750;
theta = 89;
grad_c = 1.63e-2; %speed gradient
dt = 0.0001;
tTot = 30;
cBase = 1450; %base speed for c
h_Range = -(0:1:abs(h));
speedProfile = arrayfun(@(h) getC(h,grad_c,cBase, grad_inv_h),h_Range);
time = (0:1:(tTot/dt));

%% Script

points = getPoints(source, theta ,dt, grad_c, tTot, cBase, h, grad_inv_h);
% points2 = getPoints(source, theta + 7 ,dt, grad_c, tTot, cBase, h, grad_inv_h);
% points3 = getPoints(source, theta + 15 ,dt, grad_c, tTot, cBase, h, grad_inv_h);
subplot(1,4,2:4);

plot(points(:,1), points(:,2))
% hold on
% plot(points2(:,1), points2(:,2))
% hold on
% plot(points3(:,1), points3(:,2))
% hold off
title("Ray Trayectory x,z")
subplot(1,4,1);
plot(arrayfun(@(h) getC(h,grad_c,cBase, grad_inv_h),h_Range), h_Range) 
title("Speed Profile")


%% Functions
function tracingPoints = getPoints(source, thetaOrigin, dt, grad, totalT, cStart, h, grad_inv_h)
x = 0;
z = source;
t = 0;
tracingPoints = zeros(totalT/dt, 4);
theta =  thetaOrigin ;
C = getC(z,grad,cStart, grad_inv_h);
index = 1;

while t < totalT
    % compute step
    [dz, dx]= stepDt(theta, dt, C);

    % Log values
    tracingPoints(index,1) = x;
    tracingPoints(index,2) = z;
    tracingPoints(index,3) = C;
    tracingPoints(index,4) = theta;

    % Check reflection
    if z + dz > 0 || z + dz < h
        theta = 180 - theta;
        x = x + dx;
        t = t + dt;
        index = index + 1;
        continue
    end

    % Get new C and theta
    newC = getC((z + dz), grad, cStart, grad_inv_h);
    newTheta =  getTheta(newC, C, theta);

    % Update values for next iteration
    C = newC;
    theta = newTheta;
    t = t + dt;
    x = x + dx;
    z = z + dz;
    index = index + 1;
end
end

function newTheta = getTheta(newC, C, theta)

if abs(theta - 90) < 0.01 % check if ray is aproaching 90 from below or above
    theta = theta - sign(theta - 90)*0.02; % bump ray depending on the direction it was coming from
end

% asin((C_i / C_i-1) * sin(theta)) function is symetrical and mirrored at
% 90 degrees, angles beyond 90 degrees have to be translated to the
% opposite plane in order for the output to match the input
if theta > 90
    newTheta = 180 - asind((newC / C) * sind(theta));
else
    newTheta = asind((newC / C) * sind(theta));
end
end

function [dz, dx] = stepDt(theta, dt, c)
cdt = c * dt;
dz = -cdt*cosd(theta);
dx = cdt*sind(theta);
end

function c = getC(z, grad, base, h_inv)
if z >= h_inv %this keeps the gradient function  continuous
    c = base + (-grad * abs(z)) + (grad * abs(h_inv)) - (-grad * abs(h_inv));
else
    c = base + (grad * abs(z));
end
end