function [theta_a, theta_g, phi] = find_GPR_transmission_angles(er,h,xa,ya,xq,yq,zq)
% Samuel Wagner, UC Davis ECE MML, 2021

% find_GPR_transmission_angles
% is a speed > accuracy function to find the ray-based refraction
% angles for a GPR scenario. It is inaccurate with low depth.
% it uses an iterative Newton-Raphson method to solve for the angles
% which result in the light's "path of least time"

% er - relative permittivity of ground
% h  - height of antenna above ground
% d  - depth of buried target
% xA - antenna phase center x 
% yA - antenna phase center y
% xq - x position of target
% yq - y position of target
% zq - z position of target

% thetat - transmitted (in air) angle
% thetar - refracted (in soil) angle
% phi - azimuth angle

% note that theta = 0 is pointing downwards from antenna to ground
% and       phi   = 0 is pointing along the scan direction


% convert to a single, just for speed in case we do a lot of iterations.
xa = single(xa);
ya = single(ya);
xq = single(xq);
yq = single(yq);

max_n_steps = 100; % if there is no convergence after 100 steps, return

% move antenna to origin, shift target accordingly
yq=yq-ya;          
xq=xq-xa;

% calculate phi and re-center to [0,2*pi]
phi = atan2(yq, xq);
if(phi < 0)
    phi = 2*pi + phi;
end

% xi - intermediate x-point at which the "ray" intersects gnd
% initial guess of xi=xq generally has stable solutions.
xi = xq; 

% loop variables
tol     = 1e-4;    % tolerance (meters)
step    = 2*tol;   % dummy setting to enter while() loop
nsteps  = 0;     % keep track of number of steps

% --- explanation
% take F(xi) = time which it takes the light ray to travel from antenna (0,0,-h) to target (xq,yq,zq).
% while passing through point xi.
% Find xi such that F(x) is minimized. The corresponding angles (theta_a and theta_g) will be the true
% incident transmission/refraction angles
%
% To minimize, F(xi), use Newton iterative method (based on Taylor Series) to find the zero crossing of 
% F'(xi) - first derivative. This value of xi will give the minimum time in air and angles
% that will work with Snell's law.
% the equation of F(xi) is easily found from a geometric drawing, and F'(xi) further derived.
% 
% update each step as xi(n+1) = xi(n) - F'(xi(n))/(F''(xi(n)))
% find F''(xi(n)) by calculating F'(xi(n)) and F'(xi(n)+tol).
% I found that halving the step size makes the method more stable
% This is MUCH faster compared to calculating F(x) for many x and finding the minimum. Oftentimes
% only 3-5 iterations are needed!
while(abs(step)>tol)

    %calculate value of F'(xi)
    t_center = xi*(tan(phi)^2 + 1)/sqrt(xi^2 + (xi*tan(phi))^2 + h^2) - ...
        sqrt(er)*(tan(phi)*(yq-xi*tan(phi)) + (xq-xi))/sqrt((xq-xi)^2+(yq-xi*tan(phi))^2 +zq^2);

    %calculate value of F'(xi+delta)
    t_right = (xi+tol)*(tan(phi)^2 + 1)/sqrt((xi+tol)^2 + ((xi+tol)*tan(phi))^2 + h^2) - ...
        sqrt(er)*(tan(phi)*(yq-(xi+tol)*tan(phi)) + (xq-(xi+tol)))/sqrt((xq-(xi+tol))^2+(yq-(xi+tol)*tan(phi))^2 +zq^2);

    %calculate step size as F'(xi)/F''(xi). /2 for stability
    step = -t_center/((t_right-t_center)/tol)/2;

    %update xi with step size
    xi = xi + step;
    
    %keep track of how many times we've done this
    nsteps = nsteps+1;
    if(nsteps > max_n_steps)
        xi = xq/2;          % if we are error, take a guess of halfway.
        break;
    end
end

% based on the xi found above, calculate theta_a and theta_g (from problem geometry)
theta_a=acos(h/sqrt(xi^2 + (xi.*tan(phi)).^2 + h^2));
theta_g=acos(zq/sqrt((xq-xi)^2 + (yq - xi.*tan(phi))^2 + zq^2));