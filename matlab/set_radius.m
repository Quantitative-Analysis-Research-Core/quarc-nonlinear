function [rec, diag_hist, vertical_hist, rad_final,A] = set_radius(x,a,radius_start,radius_end,threshold,type,iter)
        % Find the radius to provide target percent recurrence
        % If radius_start is too small
% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

        % Bracket the target from below. Bounded, because %REC is a step
        % function of radius for x taking few distinct values -- binary or
        % categorical series, for instance -- so a target between two steps is
        % unreachable and the search would otherwise never terminate.
        maxAdjust = 200;
        [rec, ~, ~, ~] = line_hist(x,a,radius_start,type);
        k = 0;
        while rec == 0 || rec > threshold
            k = k + 1;
            if k > maxAdjust
                error('set_radius:targetUnreachable', ...
                    ['Could not bracket %%REC = %g from below after %d ' ...
                     'adjustments (radius %g gives %%REC = %g). The distance ' ...
                     'matrix takes too few distinct values for this target to ' ...
                     'be reachable, which is usual for binary or categorical ' ...
                     'data. Set the radius directly with param="rad", or pick ' ...
                     'a target that lies on an achievable step.'], ...
                    threshold, maxAdjust, radius_start, rec);
            end
            if rec == 0
                radius_start = radius_start*2;
            elseif rec > threshold
                radius_start = radius_start / 1.5;
            end
            [rec, ~, ~, ~] = line_hist(x,a,radius_start,type);
        end

        % if radius_end is too large
        [rec, ~, ~, ~] = line_hist(x,a,radius_end,type);
        k = 0;
        while rec < threshold
            k = k + 1;
            if k > maxAdjust
                error('set_radius:targetUnreachable', ...
                    ['Could not bracket %%REC = %g from above after %d ' ...
                     'adjustments (radius %g gives %%REC = %g).'], ...
                    threshold, maxAdjust, radius_end, rec);
            end
            radius_end = radius_end*2;
            [rec, ~, ~, ~] = line_hist(x,a,radius_end,type);
        end

        % Search for radius with target percent recurrence
        lv = radius_start; % set low value
        hv = radius_end; % set high value
        target = threshold;  % designate what percent recurrence is wanted
        for  i1 = 1:iter
            mid(i1) = (lv(i1)+hv(i1))/2; % find midpoint between hv and lv
            rad(i1) = mid(i1); % new radius for this iteration

            % Compute recurrence matrix with new radius
            [rec, diag_hist, vertical_hist,A] = line_hist(x,a, rad(i1),type);
            rec_iter(i1) = rec;  % set percent recurrence
            if rec_iter(i1) < target
                % if percent recurrence is below target percent recurrence,
                % update low value
                hv(i1+1) = hv(i1);
                lv(i1+1) = mid(i1);
            else
                % if percent recurrence is above or equal to target percent
                % recurrence, update high value
                lv(i1+1) = lv(i1);
                hv(i1+1) = mid(i1);
            end
        end
        rec_final = rec_iter(end); % set final percent recurrence
        rad_final = rad(end); % set radius for final percent recurrence
        disp(['% recurrence = ',num2str(rec_final),', radius = ',num2str((rad_final))])
    end