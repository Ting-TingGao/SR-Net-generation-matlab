function [Line_list, Line_polygon, polygon_list, meta] = generate_sr_network(varargin)
% Scale-Rich (SR) 2D network generator 

% USAGE (defaults shown):
%   [Line_list, Line_polygon, polygon_list, meta] = generate_sr_network( ...
%       'side_len', 1, ...                 % domain side length (square)
%       'alpha', [], ...                  % decay exponent, range [0,1]
%       'v0', [], ...                    % initial thickness
%       'max_lines', 500, ...              % stop after this many lines
%       'max_iter_per_line', 1e2, ...      % retries for finding a valid split
%       'max_try_lines', 1e3, ...          % jamming retry limit
%       'initial_poly', polyshape([0 0 1 1],[1 0 0 1]), ... % custom domain
%       'first_point', [], ...             % [] => random interior (default)
%       'first_angle_deg', [], ...         % [] => random angle in (0,180)
%       'rng_seed', [], ...                % set to integer for reproducibility
%       'plot_result', true, ...           % draw resulting added polygons
%       'save_prefix', '' ...              % if nonempty, save .mat outputs
%   );
%
% OUTPUTS
%   Line_list    : [t x 12] numeric, columns are
%                  [mid_x1 mid_y1 mid_x2 mid_y2, left_x1 left_y1 left_x2 left_y2, right_x1 right_y1 right_x2 right_y2]
%   Line_polygon : {t x 1} cell, polyshape of added rectangular strip at step t
%   polygon_list : {t+1} cell, each entry is a cell array of void polygons at step t
%   meta         : struct with fields: t, density, angles, mid_points, v_hist, opts
%
% Author: Ting-Ting Gao

% -------------------- Parse options --------------------
p = inputParser;
addParameter(p,'side_len', 1, @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'alpha', 0.9, @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'v0', 0.05, @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'max_lines', 500, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'max_iter_per_line', 1e2, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'max_try_lines', 1e3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'initial_poly', polyshape([0 0 1 1],[1 0 0 1]));
addParameter(p,'first_point', [], @(x)(isnumeric(x)&&numel(x)==2)||isempty(x));
addParameter(p,'first_angle_deg', [], @(x)(isnumeric(x)&&isscalar(x))||isempty(x));
addParameter(p,'rng_seed', []);
addParameter(p,'plot_result', true, @(x)islogical(x)&&isscalar(x));
addParameter(p,'save_prefix', '', @(x)ischar(x)||isstring(x));
parse(p,varargin{:});
opts = p.Results;

if ~isempty(opts.rng_seed)
    rng(opts.rng_seed);
end

% -------------------- Initialization --------------------
side_len_initial = opts.side_len;
alpha            = opts.alpha;
v0               = opts.v0;
max_lines        = opts.max_lines;
max_iter_line    = opts.max_iter_per_line;
max_try_lines    = opts.max_try_lines;
initial_bound    = opts.initial_poly; 


mid_points   = [];
angles       = [];
Line_list    = [];
Line_polygon = {};
mid_angle    = [];

polygon_list       = {};
polygon_list{1}{1} = initial_bound;

density    = 0;
t          = 0;
try_time   = 1;
Notjammed  = true;


vt = side_len_initial^2/(zeta(2*alpha)*sqrt(2));

v_hist = [];

% -------------------- Main growth loop --------------------
while try_time < max_try_lines && Notjammed

    iter = 1;
    pointFound = false;

    while (iter <= max_iter_line) && ~pointFound
        % Check if any candidate polygon has a side long enough (> vt)
        max_candidate_side_length = zeros(t+1,1);
        for i = 1:t+1
            [xcentr, ycentr] = boundary(polygon_list{t+1}{i});
            side_length = sqrt((xcentr(2:end)-xcentr(1:end-1)).^2 + (ycentr(2:end)-ycentr(1:end-1)).^2);
            if ~isempty(side_length)
                max_candidate_side_length(i) = max(side_length);
            else
                max_candidate_side_length(i) = 0;
            end
        end

        if max(max_candidate_side_length) >= vt
            % Choose a polygon proportional to its area
            poly_area = zeros(numel(polygon_list{t+1}),1);
            for i = 1:numel(polygon_list{t+1})
                poly_area(i) = area(polygon_list{t+1}{i});
            end
            prob = poly_area ./ sum(poly_area);
            cumulativeProbabilities = cumsum(prob);
            r = rand();
            selectedPolygonIndex = find(cumulativeProbabilities >= r, 1);

            [tmp_xb, tmp_yb] = boundary(polygon_list{t+1}{selectedPolygonIndex});

            if t == 1 && ~isempty(opts.first_point) && isinterior(polygon_list{t+1}{selectedPolygonIndex}, opts.first_point(1), opts.first_point(2))
                xt = opts.first_point(1);
                yt = opts.first_point(2);
            else
                rdn_point = randomPointInPolygon([tmp_xb(1:end-1), tmp_yb(1:end-1)]);
                xt = rdn_point(1);
                yt = rdn_point(2);
            end

            is_inside = isinterior(polygon_list{t+1}{selectedPolygonIndex}, xt, yt);
            if ~is_inside
                iter = iter + 1; 
                continue;
            end
            chosed_polygon = selectedPolygonIndex;

            [xb, yb] = boundary(polygon_list{t+1}{chosed_polygon});
            mid_points = [mid_points; [xt, yt]]; 

            if t == 0 && ~isempty(opts.first_angle_deg)
                at = opts.first_angle_deg; 
            else
                at = 180 * rand(1);  
            end
            angles = [angles; at];

            rt = deg2rad(at);     
            lt = sqrt(2*side_len_initial);         
            x1t = xt + cos(rt)*lt; y1t = yt + sin(rt)*lt;
            x2t = xt - cos(rt)*lt; y2t = yt - sin(rt)*lt;

            
            Intersect_list = zeros(2,3); % [edgeIndex, xi, yi] for two hits
            I = 0;
            for i = 1:length(xb)-1
                [~, ints] = checkLineIntersection([x1t,y1t],[x2t,y2t],[xb(i),yb(i)],[xb(i+1),yb(i+1)]);
                if ~isempty(ints)
                    I = I+1;
                    if I<=2
                        Intersect_list(I,:) = [i, ints(1), ints(2)];
                    end
                end
            end

            vt = v0 * (t+1)^(-alpha);
            
            v_hist(end+1,1) = vt; 

            if at ~= 180 && at ~= 0
                % general case: oblique angle
                mid_left_x  = xt - (vt/2)/sin(rt);
                mid_right_x = xt + (vt/2)/sin(rt);
                line_left_t_start  = [mid_left_x +  cos(rt)*sqrt(2), y1t];
                line_left_t_end    = [mid_left_x -  cos(rt)*sqrt(2), y2t];
                line_right_t_start = [mid_right_x + cos(rt)*sqrt(2), y1t];
                line_right_t_end   = [mid_right_x - cos(rt)*sqrt(2), y2t];

                [ok, tmp_polygon, new_xb, new_yb, inters_left, inters_right] = ...
                    attempt_split(xb, yb, line_left_t_start, line_left_t_end, line_right_t_start, line_right_t_end, Intersect_list);

            else
                % horizontal (0 or 180 deg): offset in y
                line_left_t_start  = [x1t, yt - (vt/2)];
                line_left_t_end    = [x2t, yt - (vt/2)];
                line_right_t_start = [x1t, yt + (vt/2)];
                line_right_t_end   = [x2t, yt + (vt/2)];

                [ok, tmp_polygon, new_xb, new_yb, inters_left, inters_right] = ...
                    attempt_split(xb, yb, line_left_t_start, line_left_t_end, line_right_t_start, line_right_t_end, Intersect_list);
            end

            if ok
                
                signedArea = sum((new_xb(2:end) - new_xb(1:end-1)) .* (new_yb(2:end) + new_yb(1:end-1))) / 2;
                if signedArea > 0
                    New_bx = new_xb(1:end-1); New_by = new_yb(1:end-1);
                elseif signedArea < 0
                    New_bx = flip(new_xb(1:end-1)); New_by = flip(new_yb(1:end-1));
                else
                    iter = iter + 1; continue; 
                end

                ref_node = inters_left(2,2:3);
                new_start = find(all(ref_node == [New_bx, New_by],2),1);
                if isempty(new_start)
                    iter = iter + 1; continue;
                end
                front_x = New_bx(new_start:end); back_x = New_bx(1:new_start-1);
                front_y = New_by(new_start:end); back_y = New_by(1:new_start-1);
                tmp_x = [front_x; back_x]; tmp_y = [front_y; back_y];

                p1 = find(all(inters_left(1,2:3) == [tmp_x,tmp_y],2),1);
                p2 = find(all(inters_right(1,2:3) == [tmp_x,tmp_y],2),1);
                p3 = find(all(inters_right(2,2:3) == [tmp_x,tmp_y],2),1);
                if isempty(p1) || isempty(p2) || isempty(p3)
                    iter = iter + 1; continue;
                end

                new_polygon1 = polyshape([tmp_x(1:p1)], [tmp_y(1:p1)]);
                new_polygon2 = polyshape([tmp_x(p2:p3)], [tmp_y(p2:p3)]);

                % area conservation check
                if abs((area(new_polygon1)+area(new_polygon2)+area(tmp_polygon)) - area(polygon_list{t+1}{chosed_polygon})) < 1e-8
                    % record line and polygon strip
                    Line_list = [Line_list; [ ...
                        Intersect_list(1,2:3), Intersect_list(2,2:3), ... % centerline endpoints
                        inters_left(1,2:3),  inters_left(2,2:3), ...       % left strip edge
                        inters_right(1,2:3), inters_right(2,2:3)  ...      % right strip edge
                    ]]; 
                    Line_polygon{t+1,1} = tmp_polygon; 
                    mid_angle = [mid_angle; [xt, yt, at]]; 

                    % update density
                    if area(tmp_polygon) ~= 0
                        density = density + area(tmp_polygon);
                    else
                        density = density + vt * sqrt(sum((Intersect_list(1,2:3)-Intersect_list(2,2:3)).^2));
                    end

                    % update void list
                    polygon_list{t+2} = polygon_list{t+1}; 
                    polygon_list{t+2}(chosed_polygon) = [];
                    tmp_length = numel(polygon_list{t+2});
                    polygon_list{t+2}{tmp_length+1} = new_polygon1;
                    polygon_list{t+2}{tmp_length+2} = new_polygon2;

                    pointFound = true;
                else
                    iter = iter + 1; % failed area check
                end
            else
                iter = iter + 1; % failed to create valid strip
            end
        else
            % no polygon can host the current vt -> jammed
            break;
        end
    end

    if pointFound
        Notjammed = true;
        t = t + 1;
        if t >= max_lines
            break;
        end
    else
        try_time = try_time + 1;
        if try_time >= max_try_lines
            Notjammed = false;
            break;
        end
    end

    % Early stop
    if t >= max_lines
        break;
    end
end

% -------------------- Save & Plot --------------------
if ~isempty(opts.save_prefix)
    save([char(opts.save_prefix) '_line_polygon.mat'], 'Line_polygon');
    save([char(opts.save_prefix) '_line_list.mat'], 'Line_list');
end

if opts.plot_result
    figure; hold on;
    for i = 1:t
        plot(Line_polygon{i});
    end
    hold off;
    xlim([0 side_len_initial]); ylim([0 side_len_initial]);
    axis equal; axis off;
end

% -------------------- Meta outputs --------------------
meta = struct();
meta.t          = t;
meta.density    = density;
meta.angles_deg = angles;
meta.mid_points = mid_points;
meta.v_hist     = v_hist;
meta.opts       = opts;

end 

function [ok, tmp_polygon, new_xb, new_yb, inters_left, inters_right] = ...
    attempt_split(xb, yb, line_left_t_start, line_left_t_end, line_right_t_start, line_right_t_end, Intersect_list)

Il = 1; Ir = 1;
inters_left  = zeros(2,3);
inters_right = zeros(2,3);

for ii = 1:length(xb)-1
    [~, intsl] = checkLineIntersection(line_left_t_start,  line_left_t_end,  [xb(ii),yb(ii)], [xb(ii+1),yb(ii+1)]);
    [~, intsr] = checkLineIntersection(line_right_t_start, line_right_t_end, [xb(ii),yb(ii)], [xb(ii+1),yb(ii+1)]);
    if ~isempty(intsl)
        inters_left(Il,:) = [ii, intsl]; Il = Il+1;
    end
    if ~isempty(intsr)
        inters_right(Ir,:) = [ii, intsr]; Ir = Ir+1;
    end
end

ok = false; tmp_polygon = polyshape(); new_xb = []; new_yb = [];
if any(isnan(Intersect_list(:))) || any(Intersect_list(:,1)==0)
    return;
end

if ~isempty(inters_left) && ~isempty(inters_right) && Il == 3 && Ir == 3
    if sum(inters_left(:,1) - inters_right(:,1)) == 0
        tmp_polygon = polyshape([inters_left(1,2)  inters_right(1,2)  inters_right(2,2)  inters_left(2,2)], ...
                                [inters_left(1,3)  inters_right(1,3)  inters_right(2,3)  inters_left(2,3)]);
    else
        if inters_left(1,1) ~= inters_right(1,1) && inters_left(2,1) ~= inters_right(2,1)
            tmp_polygon = polyshape([inters_left(1,2)  xb(inters_left(1,1)+1:inters_right(1,1))'  inters_right(1,2)  inters_right(2,2)  xb(inters_right(2,1)+1:inters_left(2,1))'  inters_left(2,2)], ...
                                    [inters_left(1,3)  yb(inters_left(1,1)+1:inters_right(1,1))'  inters_right(1,3)  inters_right(2,3)  yb(inters_right(2,1)+1:inters_left(2,1))'  inters_left(2,3)]);
        elseif inters_left(1,1) ~= inters_right(1,1) && inters_left(2,1) == inters_right(2,1)
            tmp_polygon = polyshape([inters_left(1,2)  xb(inters_left(1,1)+1:inters_right(1,1))'  inters_right(1,2)  inters_right(2,2)  inters_left(2,2)], ...
                                    [inters_left(1,3)  yb(inters_left(1,1)+1:inters_right(1,1))'  inters_right(1,3)  inters_right(2,3)  inters_left(2,3)]);
        else
            tmp_polygon = polyshape([inters_left(1,2)  inters_right(1,2)  inters_right(2,2)  xb(inters_right(2,1)+1:inters_left(2,1))'  inters_left(2,2)], ...
                                    [inters_left(1,3)  inters_right(1,3)  inters_right(2,3)  yb(inters_right(2,1)+1:inters_left(2,1))'  inters_left(2,3)]);
        end
    end

    % Stitch new boundary polygon after inserting strip edges
    sorted_line_v = sortrows([inters_left;inters_right],1);
    new_xb = []; new_yb = [];

    if sorted_line_v(1,1) == sorted_line_v(2,1)
        if hypot(sorted_line_v(1,2)-xb(sorted_line_v(1,1)), sorted_line_v(1,3)-yb(sorted_line_v(1,1))) < ...
           hypot(sorted_line_v(2,2)-xb(sorted_line_v(1,1)), sorted_line_v(2,3)-yb(sorted_line_v(1,1)))
            new_xb = [xb(1:sorted_line_v(1,1)); sorted_line_v(1,2); sorted_line_v(2,2)];
            new_yb = [yb(1:sorted_line_v(1,1)); sorted_line_v(1,3); sorted_line_v(2,3)];
        else
            new_xb = [xb(1:sorted_line_v(1,1)); sorted_line_v(2,2); sorted_line_v(1,2)];
            new_yb = [yb(1:sorted_line_v(1,1)); sorted_line_v(2,3); sorted_line_v(1,3)];
        end
        if sorted_line_v(3,1) == sorted_line_v(4,1)
            if hypot(sorted_line_v(3,2)-xb(sorted_line_v(3,1)), sorted_line_v(3,3)-yb(sorted_line_v(3,1))) < ...
               hypot(sorted_line_v(4,2)-xb(sorted_line_v(3,1)), sorted_line_v(4,3)-yb(sorted_line_v(3,1)))
                new_xb = [new_xb; xb(sorted_line_v(2,1)+1:sorted_line_v(3,1)); sorted_line_v(3,2); sorted_line_v(4,2); xb(sorted_line_v(4,1)+1:end)];
                new_yb = [new_yb; yb(sorted_line_v(2,1)+1:sorted_line_v(3,1)); sorted_line_v(3,3); sorted_line_v(4,3); yb(sorted_line_v(4,1)+1:end)];
            else
                new_xb = [new_xb; xb(sorted_line_v(2,1)+1:sorted_line_v(3,1)); sorted_line_v(4,2); sorted_line_v(3,2); xb(sorted_line_v(4,1)+1:end)];
                new_yb = [new_yb; yb(sorted_line_v(2,1)+1:sorted_line_v(3,1)); sorted_line_v(4,3); sorted_line_v(3,3); yb(sorted_line_v(4,1)+1:end)];
            end
        end
    else
        new_xb = [xb(1:sorted_line_v(1,1)); sorted_line_v(1,2); xb(sorted_line_v(2,1)); sorted_line_v(2,2)];
        new_yb = [yb(1:sorted_line_v(1,1)); sorted_line_v(1,3); yb(sorted_line_v(2,1)); sorted_line_v(2,3)];
        if sorted_line_v(3,1) == sorted_line_v(4,1)
            if hypot(sorted_line_v(3,2)-xb(sorted_line_v(3,1)), sorted_line_v(3,3)-yb(sorted_line_v(3,1))) < ...
               hypot(sorted_line_v(4,2)-xb(sorted_line_v(3,1)), sorted_line_v(4,3)-yb(sorted_line_v(3,1)))
                new_xb = [new_xb; xb(sorted_line_v(2,1)+1:sorted_line_v(3,1)); sorted_line_v(3,2); sorted_line_v(4,2); xb(sorted_line_v(4,1)+1:end)];
                new_yb = [new_yb; yb(sorted_line_v(2,1)+1:sorted_line_v(3,1)); sorted_line_v(3,3); sorted_line_v(4,3); yb(sorted_line_v(4,1)+1:end)];
            else
                new_xb = [new_xb; xb(sorted_line_v(2,1)+1:sorted_line_v(3,1)); sorted_line_v(4,2); sorted_line_v(3,2); xb(sorted_line_v(4,1)+1:end)];
                new_yb = [new_yb; yb(sorted_line_v(2,1)+1:sorted_line_v(3,1)); sorted_line_v(4,3); sorted_line_v(3,3); yb(sorted_line_v(4,1)+1:end)];
            end
        else
            new_xb = [new_xb; xb(sorted_line_v(2,1)+1:sorted_line_v(3,1)); sorted_line_v(3,2); xb(sorted_line_v(4,1)); sorted_line_v(4,2); xb(sorted_line_v(4,1)+1:end)];
            new_yb = [new_yb; yb(sorted_line_v(2,1)+1:sorted_line_v(3,1)); sorted_line_v(3,3); yb(sorted_line_v(4,1)); sorted_line_v(4,3); yb(sorted_line_v(4,1)+1:end)];
        end
    end

    ok = true;
end

end 


% ======================================================
% =========== Helper: randomPointInPolygon =============
% ======================================================
function randomPoint = randomPointInPolygon(vertices)
    minX = min(vertices(:, 1));
    maxX = max(vertices(:, 1));
    minY = min(vertices(:, 2));
    maxY = max(vertices(:, 2));

    pointInPolygon = false;
    randomPoint = [0, 0];

    while ~pointInPolygon
        randomPoint(1) = minX + (maxX - minX) * rand();
        randomPoint(2) = minY + (maxY - minY) * rand();
        pointInPolygon = inpolygon(randomPoint(1), randomPoint(2), vertices(:, 1), vertices(:, 2));
    end
end

% ======================================================
% ========= Helper: checkLineIntersection ==============
% ======================================================
function [hasIntersection, intersectionPoint] = checkLineIntersection(p1, p2, p3, p4)

    v1 = p2 - p1;
    v2 = p4 - p3;

    cross1 = v1(1) * v2(2) - v1(2) * v2(1);
    cross2 = v2(1) * (p1(2) - p3(2)) - v2(2) * (p1(1) - p3(1));
    cross3 = v1(1) * (p1(2) - p3(2)) - v1(2) * (p1(1) - p3(1));

    threshold = 1e-16;

    if abs(cross1) < threshold %&& abs(cross2) < threshold && abs(cross3) < threshold
            hasIntersection = false;
            intersectionPoint = [];
    else
        t1 = cross2 / cross1;
        t2 = cross3 / cross1;

        if t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1
            intersectionPoint = p1 + t1 * v1;
            hasIntersection = true;
        else
            hasIntersection = false;
            intersectionPoint = [];
        end
    end
end

