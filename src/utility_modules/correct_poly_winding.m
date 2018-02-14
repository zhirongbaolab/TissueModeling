function [ uniformly_wound_faces ] = correct_poly_winding(f, v)
    tic;
    fprintf('validating mesh polygons, correcting CW/CCW differences\n');

    % adapted from: https://stackoverflow.com/questions/1165647/...
    % how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order

    uniformly_wound_faces = zeros(size(f));
    
    ccw_count = 0;
    cw_count = 0;
    % iterate over the faces
    for i=1:size(f, 1)
% -------------------------------------------------
        % calculate three axes
%         b_a = v(f(i, 2), :) - v(f(i, 1), :);
%         c_a = v(f(i, 3), :) - v(f(i, 1), :);
%         cross_ = cross(b_a, c_a);
%         
%         % determine the order
%         if cross_(3) > 0
%             ccw_count = ccw_count + 1;
% %           uniformly_wound_faces(i, :) = [f(i, 1) f(i, 3) f(i, 2)];
%             uniformly_wound_faces(i, :) = f(i, :);
%         else
%             cw_count = cw_count + 1;
%             % flip to ccw before adding to list
%             uniformly_wound_faces(i, :) = fliplr(f(i, :));
%         end       
        
%         (x2 - x1)(y2 + y1)
%         (x3 - x2)(y3 + y2)
%         (x1 - x3)(y1 + y3)
%         
%         x1 = v(f(i, 1), 1) 
%         y1 = v(f(i, 1), 2)
%         x2 = v(f(i, 2), 1)
%         y2 = v(f(i, 2), 2)
%         x3 = v(f(i, 3), 1)
%         y3 = v(f(i, 3), 2)
% -------------------------------------------------
%         sum = ((v(f(i, 2), 1) - v(f(i, 1), 1)) * v(f(i, 2), 2) + v(f(i, 1), 2)) + ...
%             ((v(f(i, 3), 1) - v(f(i, 2), 1)) * (v(f(i, 3), 2) + v(f(i, 2), 2))) + ...
%             ((v(f(i, 1), 1) - v(f(i, 3), 1)) * (v(f(i, 1), 2) + v(f(i, 3), 2)));
% -------------------------------------------------
        % determinate: (x2 - x1)(y3 - y1) - (x3 - x1)(y2 - y1)
%         det = ((v(f(i, 2), 1) - v(f(i, 1), 1)) * (v(f(i, 3), 2) - v(f(i, 1), 2))) -...
%             ((v(f(i, 3), 1) - v(f(i, 1), 1)) * (v(f(i, 2), 2) - v(f(i, 1), 2)));
% -------------------------------------------------
        
%         % determine the order
%         if det < 0
%             cw_count = cw_count + 1;
% %           uniformly_wound_faces(i, :) = [f(i, 1) f(i, 3) f(i, 2)];
%             uniformly_wound_faces(i, :) = fliplr(f(i, :));
%         else
%             ccw_count = ccw_count + 1;
%             % flip to ccw before adding to list
%             uniformly_wound_faces(i, :) = f(i, :);
%         end  
% -------------------------------------------------

    %ispolycw
    if ispolycw([v(f(i, 1), 1) v(f(i, 2), 1) v(f(i, 3), 1)], [v(f(i, 1), 2) v(f(i, 2), 2) v(f(i, 3), 2)])
       g = 0; 
    else
       h = 0;
    end

    end
    
    fprintf('mesh corrected and validated: %d ccw faces flipped\n', ccw_count);
    toc;
end

