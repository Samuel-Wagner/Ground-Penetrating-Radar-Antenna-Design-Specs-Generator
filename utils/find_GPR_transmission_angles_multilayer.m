% find_GPR_transmission_angles_multilayer.m
% a function to calculate the transmission angles in a 
% multilayer dielectric system
% inputs __________________________________
% xa - antenna x position (m)
% ya - antenna y position (m)
% xt - target x position  (m)
% yt - target y position  (m)
% d  - heights and depths. d(1) = antenna height. d(end) = last target
%      depth in medium n. vector (n x 1)
% er - relative permittivity (n x 1)
% no_input_checks - boolean. true = disable all input checks for speed 
% outputs _________________________________
% theta_guess - output angles w.r.t. normal in each dielectric layer 
% phi         - output angle (total) w.r.t. azimuth. same for ea layer.
%
% Example
% xa = 0;
% ya = 0;
% xt = 1;
% yt = 0.3;
% d  = [0.25, 0.2, 0.3];
% er = [1,   7,  9];    
% no_input_checks = true;
%
% [theta_v, phi] = GPR_transmission_angles_multilayer(xa,ya,xt,yt,d,er,no_input_checks);

function [theta_guess ,phi] = find_GPR_transmission_angles_multilayer(xa, ya, xt, yt, d, er, no_input_checks)
    
    targ_illum_dist_thresh = 1e-10;
    limit_width_thresh = pi/10^6;
    if(~exist('no_input_checks','var'))
        no_input_checks = 0;
    end
    
    n = numel(d);
    
    if(~no_input_checks)
        % input checks 
        assert(is_double_single_scalar(xa),"xa must be a scalar double or single");
        assert(is_double_single_scalar(ya),"ya must be a scalar double or single");
        assert(is_double_single_scalar(xt),"xt must be a scalar double or single");
        assert(is_double_single_scalar(yt),"yt must be a scalar double or single");
        assert(is_double_single_vector(d), "d must be a vector double or single");
        assert(is_double_single_vector(er),"er must be a vector double or single");
        assert(numel(er)==n, "er and d must have the same number of elements");
    end


    cum_d    = cumsum(d);
    targ_xyz = single([xt-xa, yt-ya, cum_d(end)]);
    phi      = atan2(targ_xyz(2),targ_xyz(1));
    phi(phi<0)=pi-phi;
    
    theta_guess     = zeros(n,1,'single');
    if(targ_xyz(1)<0)
        is_neg_flag = true;
        targ_xyz(1)=-targ_xyz(1); 
    else
        is_neg_flag = false;
    end
    
    theta_0_ig      = atan(targ_xyz(1)/cum_d(n));
    theta_guess(1)  = theta_0_ig;

    % upper (uL) and lower (lL) limits for binary search 
    uL      = pi/2 - 0.01;
    lL      = 0.01;

    for ii = 2:n
        theta_guess(ii) = asin(sqrt(er(ii-1))/sqrt(er(ii))*sin(theta_guess(ii-1)));
    end

    PTS         = zeros(n+1,3,'single');
    PTS(1,:)    = [0,0,0];
    
    target_not_illuminated = true;
    limits_are_wide        = true;
    while(target_not_illuminated && limits_are_wide)
        % let's recalculate the path

        % find thetas according to snell's law
        for ii = 2:n
            theta_guess(ii) = asin(sqrt(er(ii-1))/sqrt(er(ii))*sin(theta_guess(ii-1)));
        end

        % find the intersection points based on 
        for ii = 2:n+1
            xp = PTS(ii-1,1) + d(ii-1)*tan(theta_guess(ii-1));
            yp = PTS(ii-1,2) + d(ii-1)*tan(theta_guess(ii-1))*tan(phi);
            zp = cum_d(ii-1);
            if(ii == n+1)
               xp = xp + d(ii-1)*tan(theta_guess(ii-1));
               yp = yp + d(ii-1)*tan(theta_guess(ii-1))*tan(phi);
               zp = zp + d(ii-1); 
            end
            PTS(ii,:) = [xp, yp, zp];
        end


        % check if target is illuminated
        dist_total = radial_distance(PTS(end-1,:),PTS(end,:),1);
        dist_seg1  = radial_distance(PTS(end-1,:),targ_xyz,1);
        dist_seg2  = radial_distance(targ_xyz,PTS(end,:),1);

        extra_distance         = abs(dist_total - dist_seg1 - dist_seg2);
        target_not_illuminated = extra_distance > targ_illum_dist_thresh;
        
        
    % find out whether it is left or right, to update (in a binary search method).
        if(target_not_illuminated)
            angle_pts  = atan2(PTS(end,3) - PTS(end-1,3), PTS(end,1) - PTS(end-1,1));
            angle_targ = atan2(targ_xyz(3)- PTS(end-1,3), targ_xyz(1)- PTS(end-1,1));

            if(angle_pts < 0)
                angle_pts = pi-angle_pts; 
            end 
            if(angle_targ < 0)
                angle_targ = pi-angle_targ;
            end       

            is_left  = angle_pts > angle_targ;
            is_right = angle_pts <= angle_targ;

            if(is_left) % should increase theta
                lL = theta_guess(1);
            elseif(is_right) % should decrease theta
                uL = theta_guess(1);
            end
            theta_guess(1) = (lL + uL) / 2;
            limits_are_wide = abs(lL-uL) > limit_width_thresh;
        end
    end
    if(is_neg_flag)
        theta_guess = -theta_guess;
    end

end
