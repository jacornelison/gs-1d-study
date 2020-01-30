function [gs_params, varargout] = s4_gs_study(shield_params,varargin)%PLOT,anim,singlestat,threeshield,varargin)
%function gs_params = s4_gs_study(fb_h,fov,win_d,el_off,az_off,min_el,nRx)
%
% fb_h - forebaffle height in meters
% fov - field of view in degrees
% win_d - window diameter in meters
% el_off - elevation offset in meters. Either a vector or scalar (assumes
%       y-axis if scalar
% az_off - Azimuth offset in meters. Scalar
% min_el - Minimum observing elevation in degrees.
% n_rx - number of receivers.
% PLOT - turning plotting on (1) or off (0). Off by default.
% anim - outputs an invisible figure to varargout for gif-making.
% singlestat - If 1, will pack rx's together under a single forebaffle. 0
%       by default.
% threeshield - If 1, uses bottom FOV ray as exclusion ray to simulate
%       3-shield configuration. Default 0.

%%
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
end

vars = fieldnames(shield_params);
for i = 1:length(vars)
eval([vars{i} '= shield_params.(''' vars{i} ''');']);
end

opts = {'PLOT','anim','singlestat','threeshield','smargin',...
    'axis_window','INTEXT','OUTTEXT','LEGEND','TITLE'};
defs = {false, false, false, false, 2,...
    10, false, false, true,true};

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
            win_d = (1+2/sqrt(3))*win_d;
        case 4
            win_d = (1+sqrt(2))*win_d;
        case 5
            win_d = (1+sqrt(2*(1+1/sqrt(5))))*win_d;
    end
    n_rx = 1;
end

if PLOT
    clr = {[1,1,1]*0.6,[0,0,0]};
    if anim
        fig = figure('Visible','off');
    else
        fig = figure(1);
    end
    clf(fig)
    set(fig,'Position',[450,150,625,590])
    hold off
end

gs_params = [];

% Do this twice. Once at min el for exclusion ray.
% And again at peak forebaffle height to find inclusion ray.
for i = 1:2
    % Position and pointing vectors
    pos = zeros(1, 2);
    org = zeros(1, 2);
    bs = pos; % Boresight position
    pnt = [0, 1];
    pnt2 = [1,0];
    
    % Alternate RX for display only.
    altpos = [0,0];
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
        if i ==1
            pos = pos + [1 0]*fb_r;
            altpos = altpos - [1 0]*fb_r;
        else
            pos = pos - [1 0]*fb_r;
            altpos = altpos + [1 0]*fb_r;
        end
    elseif n_rx == 3
        encl_r = fb_r*(1+2/sqrt(3))-fb_r; % 3 Circle Packing radius
        if i ==1
            altpos = altpos - [1 0]*encl_r/2;
            pos = pos + [1 0]*encl_r;
        else
            pos = pos - [1 0]*encl_r;
            altpos = altpos + [1 0]*encl_r/2;
        end
    elseif n_rx == 4
        encl_r = fb_r*(1+sqrt(2))-fb_r; % 4 Circle Packing radius
        if i ==1
            pos = pos + [1 0]*encl_r;
            altpos = altpos - [1 0]*encl_r;
        else
            pos = pos - [1 0]*encl_r;
            altpos = altpos + [1 0]*encl_r;
        end
     elseif n_rx == 5
        encl_r = fb_r*(1+sqrt(2*(1+1/sqrt(5))))-fb_r; % 5 Circle Packing radius
        if i ==1
            pos = pos + [1 0]*encl_r;
            altpos = altpos - [1 0]*encl_r;
        else
            pos = pos - [1 0]*encl_r;
            altpos = altpos + [1 0]*encl_r;
        end
    end
    
    % Translate by dk_off
    pos = pos + dk_off;
    altpos = altpos + dk_off;
    bs = bs + dk_off;
    
    % Rotate in elevation.
    if i == 2
        fb_point = -fb_l+pos+fb_h*pnt;
        el_peak = acosd(dot(fb_point,[0,1])/(norm(fb_point)));
        el_ang = el_peak;
    else
        el_ang = 90-min_el;
    end
    
    pnt = rotate_2d(pnt,el_ang);
    pnt2 = rotate_2d(pnt2,el_ang);
    pos = rotate_2d(pos,el_ang);
    altpos = rotate_2d(altpos,el_ang);
    bs = rotate_2d(bs,el_ang);
    dk_ax = rotate_2d(dk_ax,el_ang);
    
    % Translate by any el offset.
    pos = pos + el_off;
    altpos = altpos + el_off;
    bs = bs + el_off;
    %dk_ax = dk_ax + el_off
    
    % Translate by Az offset
    pos = pos + az_off;
    altpos = altpos + az_off;
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
    
    if i == 2
        % Make inclusion ray
        fb_point = -fb_l+pos+fb_h*pnt;
        
        % Find intersection of inc. and excl. rays.
        m2 = tand(smargin); % GS lip diffraction safety margin
        b2 = fb_point(2);
        px = (b2-b1)/(m1-m2);
        P = [px m2*px+b2]; % <--- Min forebaffle dimensions
        
    else
        % Make exclusion ray
        % If we're including a scoop, the exclusion ray becomes the
        % bottom-most FOV ray.
        if threeshield
            excl = fov1;
            excl2 = fov2;
            
            % Find slope / intercept of ray
            m1 = excl2(2)/excl2(1);
            b1 = (pnt2(2)*win_d/2+pos(2))-m1*(pnt2(1)*win_d/2+pos(1));
            
            % Excl. ray angle
            excl_ang = atand(fb_h/(fb_r+win_d/2));
        else
            
            excl = -(fb_l)+pos+fb_h*pnt-(pnt2*win_d/2+pos);
            excl2 = (fb_l)+pos+fb_h*pnt-(-pnt2*win_d/2+pos);
            % Find slope / intercept of ray
            m1 = excl2(2)/excl2(1);
            b1 = (-pnt2(2)*win_d/2+pos(2))-m1*(-pnt2(1)*win_d/2+pos(1));
            
            % Excl. ray angle
            excl_ang = atand(fb_h/(fb_r+win_d/2));
        end
        
        
    end
    
    % Grab pertinent GS parameters for return variable
    if i == 1
        gs_params.excl_ang = excl_ang;
    else
        gs_params.el_peak = 90-el_peak;
        gs_params.fb_h = fb_h;
        gs_params.fb_r = fb_r;
        gs_params.gs_dim = P;
        
    end
    
    % Plotting
    if PLOT
        % These plots are purely for getting the legend right.
        % Please ignore.
        plot([-10 11],[-1000 -1000],'g')
        hold on
        plot([-10 11],[-1000 -1000],'b')
        plot([-10 11],[-1000 -1000],'m')
        plot([-10 11],[-1000 -1000],'k')
        plot([-10 11],[-1000 -1000],'r')
        
        % Az, El, Dk axis vecs
        quiver(org(1),org(2),az_ax(1),az_ax(2),0,'k') % To Azimuth
        quiver(az_ax(1),az_ax(2),el_ax(1),el_ax(2),0,'k') % Az to El
        x = pnt(1)*dk_off(2);
        y = pnt(2)*dk_off(2);
        quiver(el_ax(1),el_ax(2),x,y,0,'k') % El to Dk y
        xy2 = pnt2 * dk_off(1);
        quiver(el_ax(1)+x,el_ax(2)+y,xy2(1),xy2(2),0,'color',[0.7 0 0]) % El to Dk y
        %quiver(dk_ax(1)-el_off(1)+el_ax(1),dk_ax(2)-el_off(2)+el_ax(2),dk_ax(1),dk_ax(2),0,'k') % El to Dk
        % pnting vecs
        %quiver(pos(1),pos(2),pnt(1)*100,pnt(2)*100,0)
        %quiver(pos(1),pos(2),pnt2(1)*100,pnt2(2)*100,0)
        
        
        % window Vecs
        quiver(pos(1),pos(2),winv(1),winv(2),0,'Color',clr{i},'ShowArrowHead','off','LineWidth',4)
        quiver(pos(1),pos(2),-winv(1),-winv(2),0,'Color',clr{i},'ShowArrowHead','off','LineWidth',4)
        
        
        % fov Vec
        if i ==2
            quiver(winv(1)+pos(1),winv(2)+pos(2),fov2(1),fov2(2),0,'g','ShowArrowHead','off')
            quiver((-winv(1)+pos(1)),(-winv(2)+pos(2)),fov1(1),fov1(2),0,'g','ShowArrowHead','off')
        end
        % Forebaffle Vecs
        quiver(pos(1),pos(2),fb_l(1),fb_l(2),0,'Color',clr{i},'ShowArrowHead','off')
        quiver(pos(1),pos(2),-fb_l(1),-fb_l(2),0,'Color',clr{i},'ShowArrowHead','off')
        quiver(fb_l(1)+pos(1),fb_l(2)+pos(2),fb_h*pnt(1),fb_h*pnt(2),0,'Color',clr{i},'ShowArrowHead','off')
        quiver((-fb_l(1)+pos(1)),(-fb_l(2)+pos(2)),fb_h*pnt(1),fb_h*pnt(2),0,'Color',clr{i},'ShowArrowHead','off')
        
        if i == 1
            % Exclusion Ray Vecs
            if threeshield
                quiver(pos(1)-winv(1),pos(2)-winv(2),excl(1)*1000,excl(2)*1000,0,'m','ShowArrowHead','off')
                quiver((pos(1)+winv(1)),(pos(2)+winv(2)),excl2(1)*1000,excl2(2)*1000,0,'m','ShowArrowHead','off')
            else
                quiver(pos(1)+winv(1),pos(2)+winv(2),excl(1)*1000,excl(2)*1000,0,'m','ShowArrowHead','off')
                quiver((pos(1)-winv(1)),(pos(2)-winv(2)),excl2(1)*1000,excl2(2)*1000,0,'m','ShowArrowHead','off')
            end
            
            
        else
            % Constraints Valid area
            fill([P(1) 100 100],[P(2) m2*100+b2 m1*100+b1],[0 0.8 0]+0.2,'EdgeColor','m')
            
            % Inclusion Ray Vecs
            plot([-1000,1000],fb_point(2)+[-1000,1000]*tand(2),'b')
            

            % Min GS tip location
            plot([P(1) P(1)],[-1000 1000],'r')
            plot(P(1),P(2),'rx')
            plot(fb_point(1),fb_point(2),'bx')
        end
        
        % Draw extra window / forebaffles
        if n_rx > 1
            clrx = 0.3;
            quiver(altpos(1),altpos(2),winv(1),winv(2),0,'Color',clr{i}+clrx,'ShowArrowHead','off','LineWidth',4)
            quiver(altpos(1),altpos(2),-winv(1),-winv(2),0,'Color',clr{i}+clrx,'ShowArrowHead','off','LineWidth',4)
            quiver(altpos(1),altpos(2),fb_l(1),fb_l(2),0,'Color',clr{i}+clrx,'ShowArrowHead','off')
            quiver(altpos(1),altpos(2),-fb_l(1),-fb_l(2),0,'Color',clr{i}+clrx,'ShowArrowHead','off')
            quiver(fb_l(1)+altpos(1),fb_l(2)+altpos(2),fb_h*pnt(1),fb_h*pnt(2),0,'Color',clr{i}+clrx,'ShowArrowHead','off')
            quiver((-fb_l(1)+altpos(1)),(-fb_l(2)+altpos(2)),fb_h*pnt(1),fb_h*pnt(2),0,'Color',clr{i}+clrx,'ShowArrowHead','off')
        end
        
        % Draw Ground Shield and maybe scoop
        if i == 2% & 1
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
            
        if threeshield == 1
            st = 0.4;
            quiver(P(1)*st,0,P(1)*(1-st),P(2),0,'Color',[0 0 0.3]+0.4,'ShowArrowHead','off','LineWidth',2)
        end
        
        end
        
        
    end
end

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
        
else
    txt_str = {''};
end
if OUTTEXT
    
    t = {sprintf(''),...
        sprintf('Outputs:'),...
        sprintf('El @ FB peak (deg): %2.1f',90-el_peak),...
        sprintf('Forebaffle Radius (m): %2.1f',fb_r),...
        sprintf('Min GS Dist. (m): %2.1f',P(1)),...
        sprintf('Min GS Ht. (m): %2.1f',P(2)),...
        sprintf('Excl. Ray Angle (deg): %2.1f',excl_ang)        
        };
    txt_str = {txt_str{:} t{:}};
end



if PLOT
    text(-axis_window*txt_x,axis_window*txt_y,txt_str)
end 
if LEGEND
    legend('FOV ray','Incl. ray','Excl. Ray','Forebaffle','GS min dist.','Location','northwest')

    
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

function [pntp] = rotate_2d(pnt,th)
pnt2 = [pnt(2) -pnt(1)];
pntp = pnt*cosd(th)+pnt2*sind(th);
