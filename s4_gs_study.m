function [gs_params, varargout] = s4_gs_study(shield_params, varargin)
% function [gs_params, varargout] = s4_gs_study(shield_params, varargin)
% This function calculates the groundhield dimensions a small aperture
% telescope (SAT) for a wide range of input parameters based on the
% Double Diffraction Criteria listed below:
%   C1: No portion of the absorptive forebaffle may protrude above the height of the reflective ground shield.
%   C2: No ray originating from any receiver window may couple directly with the reflective ground shield.
%   C3: No ray within the field of view (FOV) cone may couple with the forebaffle or groundshield.
%
% [Input]
%   shield_params:  Struct containing experimental parameters. All float types.
%       .fb_h       Forebaffle height in meters
%       .fov        Field of view in degrees
%       .win_d      window diameter in meters
%       .dk_off     Deck offset in meters. Vector: [x,y]
%       .el_off     Elevation offset in meters. Vector: [x,y]
%       .az_off     Azimuth offset in meters. Vector: [x,y]
%       .min_el     Minimum observing elevation in degrees.
%       .n_rx       Number of receivers.
%       [.gs_dim]   Forces groundshield dimensions. Useful for comparing to
%                   pre-existing experiments. Optional.
%
%   shield_params:  Optionally, input a string corresponding to an
%                   experiment in the function get_shield_params to
%                   automatically retrieve parameters. Useful for
%                   experiments with known input parameters.
%                   Example: s4_gs_study('BICEP2')
%
%   varargin:       {Input type | default} Miscellaneous arguments the changes behavior
%       PLOT        {Bool | 0} Turning plotting on (1) or off (0).
%       anim        {Bool | 0} Outputs an invisible figure to varargout for gif-making.
%       singlestat  {Bool | 0} If true, will pack rx's together under a single forebaffle.
%       threeshield {Float/Int | 0} Inserts a tertiary shield (called a scoop)
%                   with the input length in meters. Exclusion ray will be
%                   defined by the length of the forebaffle + length of
%                   the scoop.
%       ts_dim      {bool | 1} Determines whether to define the scoop by
%                   'radius' mor 'height'. Default true for radius. select
%                   false for height.
%       spacing     {Float/Int} Additional window-to-window spacing (m) for
%                   singlestat option (only for 3-rx for now)
%       fixwindist  {Float | 0} If value given, don't close-pack rxs and 
%                   instead fix the enclosed radius
%       minscoop    {Float} Automatically determines the scoop length based on by
%                   finding the intersection between the bottom FOV ray of
%                   the RX at the bottom of the drum and the Exclusion ray
%                   of the RX rotated to the top of the drum.
%       showalt     {Bool | 1} When plotting, also plots RX's at the peak
%                   forebaffle height.
%       
%
% [Output]
%   gs_params       {Units} Struct containing derived groundshield data.
%       .excl_ang   {degrees} Exclusion angle defined by the forebaffle or
%                   forebaffle + scoop heights.
%       .el_peak    {degrees} Elevation at which the forebaffles, in their
%                   highest configuration, are at maximum height.
%       .fb_h       {meters} Input forebaffle height. Do we need this?
%       .fb_r       {meters} Forebaffle radius required for input
%                   forebaffle height to graze the instrument FOV.
%       .gs_dim     {meters} Minimum [x,y] groundshield dimensions that satisfy all criterion.
%
%   varargout       When 'anim' is toggled, this will output the figure
%                   for gif-making.
%
% [Examples]
%   gsp = s4_gs_study(shield_params);
%   gsp = s4_gs_study('keck');
%   gsp = s4_gs_study('BA','PLOT',true,'threeshield',2);
%   [gsp, fig] = s4_gs_study('BA','anim',true);
%
%
% function [gs_params, varargout] = s4_gs_study(shield_params, varargin)


% Initialize variables.
% 3 GS panels. Determine panel locations.
g3 = 0.25;
g1 = (1-g3)/2;
g2 = g1;
gd = 25;

% Extract variables from the input struct.
if isempty(shield_params)
    shield_params = get_shield_params('BA');
elseif ischar(shield_params)
    expt = shield_params;
    shield_params = get_shield_params(shield_params);
elseif isfield(shield_params,'expt')
    expt = shield_params.expt;
end
% Otherwise, struct should already contain requisite parameters

vars = fieldnames(shield_params);
for i = 1:length(vars)
    eval([vars{i} '= shield_params.(''' vars{i} ''');']);
end

opts = {'PLOT','anim','singlestat','threeshield','smargin',...
    'axis_window','INTEXT','OUTTEXT','LEGEND','TITLE',...
    'spacing','fixwindist','ts_dim', 'minscoop','showalt'};
defs = {false, false, false, false, 2,...
    10, false, false, true,true,...
    0, false,true, false, true};

for i = 1:length(varargin)
    for j = 1:length(opts)
        %keyboard()
        if ischar(varargin{i}) & strcmp(varargin{i},opts{j})
            if ischar(varargin{i+1})
                s = ([varargin{i} '=' varargin{i+1} ';']);
            else
                s = ([varargin{i} '=' num2str(varargin{i+1}) ';']);
            end
            eval(char(s));
        end
    end
end

for i = 1:length(opts)
    if ~exist(opts{i},'var')
        eval([opts{i} '=' num2str(defs{i}) ';']);
    end
end

if singlestat
    switch n_rx
        case 1
        case 2
            win_d = 2*win_d;
        case 3
            %win_d = (1+2/sqrt(3))*win_d;
            % Take spacing into account
            
            win_r = win_d / 2.0;
            R = win_r + (win_r + (spacing/2))*(2*sqrt(3)/3);
            win_d = 2*R;
        case 4
            win_d = (1+sqrt(2))*win_d;
        case 5
            win_d = (1+sqrt(2*(1+1/sqrt(5))))*win_d;
    end
    n_rx = 1;
end

if PLOT
    %clr = {[1,1,1]*0.6,[0,0,0]};
    clr = {[0,0,0],[1,1,1]*0.6};
    if anim
        fig = figure('Visible','off');
    else
        fig = figure(1);
    end
    clf(fig)
    %set(fig,'Position',[450,150,625,590])
    hold off
end

gs_params = [];


% Do this twice. Once at min el for exclusion ray.
% And again at peak forebaffle height to find inclusion ray.
for i = 1:2
    % Position and pointing vectors
    pos_bottom = zeros(1, 2);   % Main RX bottom of the 'drum'
    pos_top = zeros(1,2);       % Main RX top of the 'drum'
    org = zeros(1, 2);
    bs = pos_bottom; % Boresight position
    pnt = [0, 1];
    pnt2 = [1,0];
    
    % Alternate RX for display only.
    altpos_bottom = [0,0];
    altpos_top = [0,0];
    altpnt = [0,1];
    altpnt2 = [1,0];
    
    % Az, El, Dk axis vectors
    az_ax = org + az_off;
    el_ax = org + el_off;
    dk_ax = org + dk_off;
    
    % start by finding location at which the fov ray is fb_height above window.
    
    fb_r = (win_d/2 + fb_h*tand(fov/2));
    fb_l = pnt2*fb_r;
    
    % Translate window from boresight if multiple rx's.
    if n_rx==2
        
        pos_bottom = pos_bottom + [1 0]*fb_r;
        altpos_bottom = altpos_bottom - [1 0]*fb_r;
        pos_top = pos_top - [1 0]*fb_r;
        altpos_top = altpos_top + [1 0]*fb_r;
        
    elseif n_rx == 3
        if fixwindist
            encl_r = fixwindist;
        else
            encl_r = fb_r*(1+2/sqrt(3))-fb_r; % 3 Circle Packing radius
        end
        
        pos_bottom = pos_bottom + [1 0]*encl_r;
        altpos_bottom = altpos_bottom - [1 0]*encl_r/2;
        pos_top = pos_top - [1 0]*encl_r;
        altpos_top = altpos_top + [1 0]*encl_r/2;
        
    elseif n_rx == 4
        if fixwindist
            encl_r = fixwindist;
        else
            encl_r = fb_r*(1+sqrt(2))-fb_r; % 4 Circle Packing radius
        end
        pos_bottom = pos_bottom + [1 0]*encl_r;
        altpos_bottom = altpos_bottom - [1 0]*encl_r;
        pos_top = pos_top - [1 0]*encl_r;
        altpos_top = altpos_top + [1 0]*encl_r;
        
    elseif n_rx == 5
        if fixwindist
            encl_r = fixwindist;
        else
            encl_r = fb_r*(1+sqrt(2*(1+1/sqrt(5))))-fb_r; % 5 Circle Packing radius
        end
        pos_bottom = pos_bottom + [1 0]*encl_r;
        altpos_bottom = altpos_bottom - [1 0]*encl_r;
        pos_top = pos_top - [1 0]*encl_r;
        altpos_top = altpos_top + [1 0]*encl_r;
        
    end
    
    % Translate by dk_off
    pos_bottom = pos_bottom + dk_off;
    altpos_bottom = altpos_bottom + dk_off;
    pos_top = pos_top + dk_off;
    altpos_top = altpos_top + dk_off;
    bs = bs + dk_off;
    
    % Rotate in elevation.
    if i == 2
        fb_point = -fb_l+pos_top+fb_h*pnt;
        el_peak = acosd(dot(fb_point,[0,1])/(norm(fb_point)));
        el_ang = el_peak;
    else
        el_ang = 90-min_el;
    end
    
    pnt = rotate_2d(pnt,el_ang);
    pnt2 = rotate_2d(pnt2,el_ang);
    
    pos_bottom = rotate_2d(pos_bottom,el_ang);
    altpos_bottom = rotate_2d(altpos_bottom,el_ang);
    pos_top = rotate_2d(pos_top,el_ang);
    altpos_top = rotate_2d(altpos_top,el_ang);
    
    bs = rotate_2d(bs,el_ang);
    dk_ax = rotate_2d(dk_ax,el_ang);
    
    % Translate by any el offset.
    pos_bottom = pos_bottom + el_off;
    altpos_bottom = altpos_bottom + el_off;
    pos_top = pos_top + el_off;
    altpos_top = altpos_top + el_off;
    
    bs = bs + el_off;
    %dk_ax = dk_ax + el_off
    
    % Translate by Az offset
    pos_bottom = pos_bottom + az_off;
    altpos_bottom = altpos_bottom + az_off;
    pos_top = pos_top + az_off;
    altpos_top = altpos_top + az_off;
    
    bs = bs + az_off;
    %dk_ax = dk_ax + az_off;
    %el_ax = el_ax + az_off;
    
    
    % Make window
    winv = pnt2*win_d/2;
    
    % Make FOV cone
    % This is the +/- rotation from normal to window edges
    fov1 = rotate_2d(pnt,-fov/2)*1000;
    fov2 = rotate_2d(pnt,fov/2)*1000;
    
    % Remake forebaffle.
    fb_l = pnt2*(win_d/2 + fb_h*tand(fov/2));
    
    % Diffraction Safety Margin
    m2 = tand(smargin); % GS lip diffraction safety margin
    
    if i == 2
        % Make inclusion ray
        fb_point = -fb_l+pos_top+fb_h*pnt;
        
        % Find intersection of inc. and excl. rays.
        
        %px = (b2-b1)/(m1-m2);
        %P = [px m2*px+b2]; % <--- Min groundshield dimensions
        if minscoop
            b2 = P_int(2)+m2*P_int(1);
            P = get_intersection(excl_line,[m2 b2]);
        else
            b2 = fb_point(2);
            P = get_intersection(line1,[m2 b2]);
            P_alt = get_intersection(line2, [m2,b2]);
            Pline = [line1];
            Ps = [P];
            
            if ~singlestat
                Pline = [Pline; line2];
                Ps = [Ps; P_alt];
            end
            
            if threeshield & ~singlestat
                P_graze = get_intersection(grazeline, [m2,b2]);
                Ps = [Ps; P_graze];
                Pline = [Pline; grazeline];
            end
            
            [Pmax, Imax] = max(Ps(:,1));
            P = Ps(Imax,:);
            Pline = Pline(Imax,:);
        end
        
    end
    %else
    % Make exclusion ray
    % If we're including a tertiary shield (or scoop),
    % the new exclusion ray is defined by the ray that is traced
    % from the top-most part of the bottom-most window
    % to the tip of the scoop.
    if minscoop & i==1
        % Find the exclusion ray of the RX at the top of the drum.
        [excl_alt, excl2_alt, excl_line, excl_ang_alt] = make_excl_rays(pos_top,pnt,pnt2,fb_l,fb_r,fb_h,win_d);
        
        % Make the FOV rays
        [f1, f2, fov_line, fov_ang] = make_fov_rays(pos_bottom,pnt,pnt2,fb_l,fb_r,fb_h,win_d);
        
        % Find the intersection of that ray and the fov ray
        P_int = get_intersection(fov_line,excl_line);
        
        % Determine the scoop height and radius
        Q = P_int-pos_bottom;
        th = acos(dot(Q/norm(Q),pnt/norm(pnt)));
        %h_tot = norm(Q)*cos(th);
        threeshield = norm(Q)*sin(th);
        ts_dim = true;
        %gs_dim = [0 0]; %P_int+[1,0];
        
    end
    
    
    
    if threeshield
        % start by finding location at which the fov ray is fb_height+scoop above window.
        if ts_dim
            ts_r = threeshield;
            ts_h = (ts_r-win_d/2)/tand(fov/2);
        else
            ts_h = fb_h+threeshield;
            ts_r = (win_d/2 + (ts_h)*tand(fov/2));
        end
        
        ts_l = pnt2*ts_r;
        
        
        if i==1
            [excl, excl2, line1, excl_ang] = make_excl_rays(pos_bottom,pnt,pnt2,ts_l,ts_r,ts_h,win_d);
            [excl_alt, excl2_alt, line2, excl_ang_alt] = make_excl_rays(pos_top,pnt,pnt2,fb_l,fb_r,fb_h,win_d);
            
            exposure_angle = excl_ang-excl_ang_alt;
            % Make the exclusion ray that just grazes the scoop
            ts_pos = ts_l+pos_bottom+(ts_h*pnt);
            b_graze = ts_pos(2)-line2(1)*ts_pos(1);
            grazeline = [line2(1), b_graze];
            th_line = atand(line2(1));
            grazeline_start = [ts_pos(1)-ts_h/cosd(90-excl_ang_alt)*cosd(th_line), 0];
            grazeline_start(2) = grazeline_start(1)*line2(1)+b_graze;
            
        else
            [excl_pk, excl2_pk, line1_pk, excl_ang_pk] = make_excl_rays(pos_top,pnt,pnt2,fb_l,fb_r,fb_h,win_d);
            [excl_alt_pk, excl2_alt_pk, line2_pk, excl_ang_alt_pk] = make_excl_rays(altpos_bottom,pnt,pnt2,fb_l,fb_r,fb_h,win_d);
        end
        
    else
        
        if i==1
            [excl, excl2, line1, excl_ang] = make_excl_rays(pos_bottom,pnt,pnt2,fb_l,fb_r,fb_h,win_d);
            [excl_alt, excl2_alt, line2, excl_ang_alt] = make_excl_rays(pos_top,pnt,pnt2,fb_l,fb_r,fb_h,win_d);
            
            
        else
            [excl_pk, excl2_pk, line1_pk, excl_ang_pk] = make_excl_rays(pos_top,pnt,pnt2,fb_l,fb_r,fb_h,win_d);
            [excl_alt_pk, excl2_alt_pk, line2_pk, excl_ang_alt_pk] = make_excl_rays(altpos_bottom,pnt,pnt2,fb_l,fb_r,fb_h,win_d);
        end
        m1 = line1(1);
        b1 = line1(2);
    end
    
    
    %end
    
    % Grab pertinent GS parameters for return variable
    if i == 1
        gs_params.excl_ang = excl_ang;
        if threeshield
            gs_params.ts_dim = ts_l+pos_bottom+(ts_h*pnt);
            gs_params.exp_angle = exposure_angle;
        end
    else
        gs_params.el_peak = 90-el_peak;
        gs_params.fb_h = fb_h;
        gs_params.fb_r = fb_r;
        if P(1)>0
           gs_params.gs_dim = P;
        else
           P = [1000,1000];
           gs_params.gs_dim = [NaN,NaN];
        end
        
        
    end
    
    % Plotting
    if PLOT
        if showalt | i==1
        % These plots are purely for getting the legend right.
        % Please ignore.
        plot([0 0],[0 0],'g')
        hold on
        plot([0 0],[0 0],'b')
        plot([0 0],[0 0],'m')
        plot([0 0],[0 0],'k')
        plot([0 0],[0 0],'r')
        
        % Az, El, Dk axis vecs
        quiver(org(1),org(2),az_ax(1),az_ax(2),0,'k') % To Azimuth
        quiver(az_ax(1),az_ax(2),el_ax(1),el_ax(2),0,'k') % Az to El
        x = pnt(1)*dk_off(2);
        y = pnt(2)*dk_off(2);
        quiver(el_ax(1),el_ax(2),x,y,0,'k') % El to Dk y
        xy2 = pnt2 * dk_off(1);
        quiver(el_ax(1)+x,el_ax(2)+y,xy2(1),xy2(2),0,'color',[0.7 0 0]) % El to Dk y
        
        
        % window Vecs
        quiver(pos_bottom(1),pos_bottom(2),winv(1),winv(2),0,'Color',clr{i},'ShowArrowHead','off','LineWidth',4)
        quiver(pos_bottom(1),pos_bottom(2),-winv(1),-winv(2),0,'Color',clr{i},'ShowArrowHead','off','LineWidth',4)
        
        
        % Forebaffle Vecs
        quiver(pos_bottom(1),pos_bottom(2),fb_l(1),fb_l(2),0,'Color',clr{i},'ShowArrowHead','off')
        quiver(pos_bottom(1),pos_bottom(2),-fb_l(1),-fb_l(2),0,'Color',clr{i},'ShowArrowHead','off')
        quiver(fb_l(1)+pos_bottom(1),fb_l(2)+pos_bottom(2),fb_h*pnt(1),fb_h*pnt(2),0,'Color',clr{i},'ShowArrowHead','off')
        quiver((-fb_l(1)+pos_bottom(1)),(-fb_l(2)+pos_bottom(2)),fb_h*pnt(1),fb_h*pnt(2),0,'Color',clr{i},'ShowArrowHead','off')
        %end
        if i == 1
            % FOV Vecs
            quiver(winv(1)+pos_bottom(1),winv(2)+pos_bottom(2),fov2(1),fov2(2),0,'g','ShowArrowHead','off')
            quiver((-winv(1)+pos_bottom(1)),(-winv(2)+pos_bottom(2)),fov1(1),fov1(2),0,'g','ShowArrowHead','off')
            
            % Exclusion Ray Vecs
            
            quiver((pos_bottom(1)-winv(1)),(pos_bottom(2)-winv(2)),excl2(1)*1000,excl2(2)*1000,0,'m','ShowArrowHead','off')
            
            if ~singlestat
                %quiver(altpos(1)+winv(1),altpos(2)+winv(2),excl_alt(1)*1000,excl_alt(2)*1000,0,'Color',mag_clr*0.5,'ShowArrowHead','off')
                quiver((pos_top(1)-winv(1)),(pos_top(2)-winv(2)),excl2_alt(1)*1000,excl2_alt(2)*1000,0,'m--','ShowArrowHead','off')
            end
            % Draw the scoop if we want to.
            if threeshield
                quiver((ts_l(1)+pos_bottom(1)),(ts_l(2)+pos_bottom(2)),ts_h*pnt(1),ts_h*pnt(2),0,'Color',[0 0 0.3]+0.2,'ShowArrowHead','off','LineWidth',3)
                quiver((pos_bottom(1)),(pos_bottom(2)),ts_l(1),ts_l(2),0,'Color',[0 0 0.3]+0.2,'ShowArrowHead','off','LineWidth',3)
            end
            
            % Draw the graze line
            if threeshield & ~singlestat
                plot([grazeline_start(1),1000],grazeline_start(2)+[0,1000*line2(1)],'m-.')
            end
        end
        
        % Draw extra window / forebaffles
        if n_rx > 1
            clrx = 0.1;

            quiver(pos_top(1),pos_top(2),winv(1),winv(2),0,'Color',clr{i}+clrx,'ShowArrowHead','off','LineWidth',4)
            quiver(pos_top(1),pos_top(2),-winv(1),-winv(2),0,'Color',clr{i}+clrx,'ShowArrowHead','off','LineWidth',4)
            quiver(pos_top(1),pos_top(2),fb_l(1),fb_l(2),0,'Color',clr{i}+clrx,'ShowArrowHead','off')
            quiver(pos_top(1),pos_top(2),-fb_l(1),-fb_l(2),0,'Color',clr{i}+clrx,'ShowArrowHead','off')

            
            quiver(fb_l(1)+pos_top(1),fb_l(2)+pos_top(2),fb_h*pnt(1),fb_h*pnt(2),0,'Color',clr{i}+clrx,'ShowArrowHead','off')
            quiver((-fb_l(1)+pos_top(1)),(-fb_l(2)+pos_top(2)),fb_h*pnt(1),fb_h*pnt(2),0,'Color',clr{i}+clrx,'ShowArrowHead','off')
        end
        end
        % Draw Ground Shield
        if i == 2
            % Constraints Valid area
            % Inclusion Ray Vecs
            % Min GS tip location
            if minscoop
                fill([P(1) 100 100],[P(2) m2*100+P(2) excl_line(1)*100+excl_line(2)],[0 0.8 0]+0.2,'EdgeColor','m','LineStyle','--')
                plot([-1000,P_int(1)],[1000*tand(smargin)+b2,P_int(2)],'b--')
            else
                fill([P(1) 100 100],[P(2) m2*100+P(2) Pline(1)*100+Pline(2)],[0 0.8 0]+0.2,'EdgeColor','m','LineStyle','--')
            end
            plot([-1000,1000],b2+[-1000,1000]*tand(smargin),'b')
            plot(P(1),P(2),'rx')
            plot(fb_point(1),fb_point(2),'bx')
            
            % Force MAPO groundshield dimensions:
            if exist('gs_dim','var')
                P = gs_dim;
            end
            
            quiver(0,0,P(1)*g1,0,0,'Color','k','ShowArrowHead','off','LineWidth',1)
            quiver(P(1)*g1,0,P(1)*g2,P(1)*g2*tand(gd),0,'Color','k','ShowArrowHead','off','LineWidth',1)
            quiver(P(1)*(g1+g2),P(1)*g2*tand(gd),P(1)*g3,(P(2)-P(1)*g2*tand(gd)),0,'Color','k','ShowArrowHead','off','LineWidth',1)
            
            quiver(0,0,-P(1)*g1,0,0,'Color','k','ShowArrowHead','off','LineWidth',1)
            quiver(-P(1)*g1,0,-P(1)*g2,P(1)*g2*tand(gd),0,'Color','k','ShowArrowHead','off','LineWidth',1)
            quiver(-P(1)*(g1+g2),P(1)*g2*tand(gd),-P(1)*g3,(P(2)-P(1)*g2*tand(gd)),0,'Color','k','ShowArrowHead','off','LineWidth',1)
            
            plot([P(1),P(1)],[-1000,P(2)],'r')
            
        end
        
        
    end
end

if PLOT
    axis equal
    xlim([-1 1]*axis_window)
    ylim([-.20 1.80]*axis_window/(2.5-(INTEXT | OUTTEXT)))
    txt_x = .9;
    txt_y = .85;
    
    if INTEXT
        txt_str = {sprintf('Inputs:'),...
            sprintf('Forebaffle Height (m): %2.1f',fb_h),...
            sprintf('Window Diameter (cm): %2.1f',win_d*100),...
            sprintf('FOV (deg): %2.1f',fov),...
            sprintf('Min Obs El (deg): %2.1f',min_el),...
            sprintf('Diffraction Safety Margin (deg): %2.1f',smargin)};
        %sprintf('Elevation Offset (x,y in m): %2.1f, %2.1f',dk_off(1),dk_off(2)),...
        %sprintf('Number of RX''s: %2.0f',n_rx),...
        if threeshield
            if ts_dim
                ts_dim_txt = 'Radius';
            else
                ts_dim_txt = 'Height';
            end
            txt_str{end+1} = sprintf(['Tertiary Shield ', ts_dim_txt,' (m): %2.1f'],threeshield);
        end
    else
        txt_str = {''};
    end
    if OUTTEXT
        
        t = {sprintf(''),...
            sprintf('Outputs:'),...
            sprintf('El @ FB peak (deg): %2.1f',90-el_peak),...
            sprintf('Forebaffle Radius (m): %2.1f',fb_r),...
            sprintf('Excl. Ray Angle (deg): %2.1f',excl_ang),...
            sprintf('Min GS dim [x,y] (m): [ %2.1f, %2.1f ]',gs_params.gs_dim(1),gs_params.gs_dim(2))...
            };
        
        txt_str = {txt_str{:} t{:}};
        
        if threeshield
            txt_str{end+1} = sprintf('Tertiary Shield Tip [x,y] (m): [ %2.1f, %2.1f ]',gs_params.ts_dim(1),gs_params.ts_dim(2));
            txt_str{end+1} = sprintf('TS Exposure Angle (deg): %2.1f',exposure_angle);
        end
        
    end
    
    
    
    if PLOT
        text(-axis_window*txt_x,axis_window*txt_y,txt_str)
    end
    if LEGEND
        legend('FOV ray','Incl. ray','Excl. Ray','Forebaffle','GS min dist.','Location','southwest')
        
        
    end
    ylabel('Height (m)')
    xlabel('Distance (m)')
    
    if TITLE
        
        if exist('expt','var')
            tname = expt;
        else
            rxtitle = {'Single','Double','Triple','Quad','Quint'};
            tname = [rxtitle{n_rx} ' RX'];
        end
        
        if singlestat ~= 0
            tname = [tname ', shared cryostat'];
        end
        if threeshield ~= 0
            tname = [tname ', with scoop'];
        end
        title(tname)
        grid on
        
        if anim
            varargout{1} = fig;
        end
    end
end

function [pntp] = rotate_2d(pnt,th)
pnt2 = [pnt(2) -pnt(1)];
pntp = pnt*cosd(th)+pnt2*sind(th);

function P = get_intersection(line1,line2)
% Gives the intersection point of two lines
% [input]
%   line1 - [slope, intercept]
%   line2 - [slope, intercept]

px = (line2(2)-line1(2))/(line1(1)-line2(1));
P = [px line2(1)*px+line2(2)];

function [fov1, fov2, line, fov_ang] = make_fov_rays(pos,pnt,pnt2,fb_l,fb_r,fb_h,win_d)
fov1 = (fb_l)+pos+fb_h*pnt-(pnt2*win_d/2+pos);
fov2 = -(fb_l)+pos+fb_h*pnt-(-pnt2*win_d/2+pos);
% Find slope / intercept of ray
m1 = fov1(2)/fov1(1);
b1 = (pnt2(2)*win_d/2+pos(2))-m1*(pnt2(1)*win_d/2+pos(1));

line = [m1, b1];

% Excl. ray angle
fov_ang = atand(fb_h/(fb_r+win_d/2));

function [excl, excl2, line1, excl_ang] = make_excl_rays(pos,pnt,pnt2,fb_l,fb_r,fb_h,win_d)
excl = -(fb_l)+pos+fb_h*pnt-(pnt2*win_d/2+pos);
excl2 = (fb_l)+pos+fb_h*pnt-(-pnt2*win_d/2+pos);
% Find slope / intercept of ray
m1 = excl2(2)/excl2(1);
b1 = (-pnt2(2)*win_d/2+pos(2))-m1*(-pnt2(1)*win_d/2+pos(1));

line1 = [m1, b1];

% Excl. ray angle
excl_ang = atand(fb_h/(fb_r+win_d/2));


