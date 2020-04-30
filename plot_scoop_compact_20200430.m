
figdir = './figs/';

%% Compact Configuration
sp.fov = 29.0;
sp.win_d = 0.72;
sp.dk_off = [0., 0.9271];
sp.el_off = [0., 2.3];
sp.az_off = [0., 0.];
sp.min_el = 50.0;
sp.n_rx = 3;
%sp.gs_dim = [0,0];

scoop_h = 0:0.01:3;
fb_h = [1.75, 2.5, 3.5];
gs_r = nan(length(fb_h), length(scoop_h));
gs_h = nan(length(fb_h), length(scoop_h));
fb_ang = nan(length(fb_h), length(scoop_h));
exp_angle = nan(length(fb_h), length(scoop_h));

for fbs = 1:length(fb_h)
    sp.fb_h = fb_h(fbs);
    parm = s4_gs_study(sp,'PLOT',false,'OUTTEXT',true,'axis_window',15,'spacing',0.1,'singlestat',true,'ts_dim',false);
    %printname = [figdir sprintf('SAT_3RX_compact_fb_%i_noscoop',ceil(fbs))];
    %print(1,printname,'-dpng')
    
    for scoops = 1:length(scoop_h)
        parm = s4_gs_study(sp,'PLOT',false,'OUTTEXT',true,'axis_window',15,'spacing',0.1,'singlestat',true,'ts_dim',false,'threeshield',scoop_h(scoops));
        %printname = [figdir sprintf('SAT_3RX_compact_fb_%i_scoop_%i',ceil(fbs),scoops)];
        %print(1,printname,'-dpng')
	gs_r(fbs, scoops) = parm.gs_dim(1);
	gs_h(fbs, scoops) = parm.gs_dim(2);
	if ~isnan(parm.gs_dim(1));
	  fb_ang(fbs, scoops) = parm.excl_ang;
	  if isfield(parm, 'exp_angle')
	    exp_angle(fbs, scoops) = parm.exp_angle;
	    fb_ang(fbs, scoops) = parm.excl_ang - parm.exp_angle;
	  end
	end
    end
end

figure(1); clf;
subplot(3,1,[1 2])
plot(scoop_h, gs_r(1,:), 'r', 'LineWidth',2);
hold on;
plot(scoop_h, gs_h(1,:), '--r', 'LineWidth',2);
plot(scoop_h, gs_r(2,:), 'g', 'LineWidth',2);
plot(scoop_h, gs_h(2,:), '--g', 'LineWidth',2);
plot(scoop_h, gs_r(3,:), 'b', 'LineWidth',2);
plot(scoop_h, gs_h(3,:), '--b', 'LineWidth',2);
ylabel('Ground Screen Size [m]')
ylim([0, 30]);
grid on;
legend('1.75 m FB, Radius','1.75 m FB, Height',...
    '2.5 m FB, Radius','2.5 m FB, Height',...
    '3.5 m FB, Radius','3.5 m FB, Height');
title('Three Shields, 50 deg min el, Compact')

subplot(3,1,3)
plot(scoop_h, exp_angle(1,:), 'r', 'LineWidth',2);
hold on;
plot(scoop_h, fb_ang(1,:), '--r', 'LineWidth',2);
plot(scoop_h, exp_angle(2,:), 'g', 'LineWidth',2);
plot(scoop_h, fb_ang(2,:), '--g', 'LineWidth',2);
plot(scoop_h, exp_angle(3,:), 'b', 'LineWidth',2);
plot(scoop_h, fb_ang(3,:), '--b', 'LineWidth',2);
ylabel('Exposure angle [deg]');
xlabel('Scoop Height [m]');
ylim([0, 60])
grid on;
legend('1.75 m FB, scoop ang','1.75 m FB, FB ang',...
    '2.5 m FB, scoop ang', '2.5 m FB, FB ang',...
    '3.5 m FB, scoop ang', '3.5 m FB, FB ang',...
    'Location','West','Numcolumns',2)
f = gcf;
f.Position = [100 100 600 600];
print(1, [figdir ...
      'compact_scoop_50deg'],'-dpng')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compact Configuration
sp.fov = 29.0;
sp.win_d = 0.72;
sp.dk_off = [0., 0.9271];
sp.el_off = [0., 2.3];
sp.az_off = [0., 0.];
sp.min_el = 45.0;
sp.n_rx = 3;
%sp.gs_dim = [0,0];

scoop_h = 0:0.01:3;
fb_h = [1.75, 2.5, 3.5];
gs_r = nan(length(fb_h), length(scoop_h));
gs_h = nan(length(fb_h), length(scoop_h));
fb_ang = nan(length(fb_h), length(scoop_h));
exp_angle = nan(length(fb_h), length(scoop_h));

for fbs = 1:length(fb_h)
    sp.fb_h = fb_h(fbs);
    parm = s4_gs_study(sp,'PLOT',false,'OUTTEXT',true,'axis_window',15,'spacing',0.1,'singlestat',true,'ts_dim',false);
    %printname = [figdir sprintf('SAT_3RX_compact_fb_%i_noscoop',ceil(fbs))];
    %print(1,printname,'-dpng')
    
    for scoops = 1:length(scoop_h)
        parm = s4_gs_study(sp,'PLOT',false,'OUTTEXT',true,'axis_window',15,'spacing',0.1,'singlestat',true,'ts_dim',false,'threeshield',scoop_h(scoops));
        %printname = [figdir sprintf('SAT_3RX_compact_fb_%i_scoop_%i',ceil(fbs),scoops)];
        %print(1,printname,'-dpng')
	gs_r(fbs, scoops) = parm.gs_dim(1);
	gs_h(fbs, scoops) = parm.gs_dim(2);
	if ~isnan(parm.gs_dim(1));
	  fb_ang(fbs, scoops) = parm.excl_ang;
	  if isfield(parm, 'exp_angle')
	    exp_angle(fbs, scoops) = parm.exp_angle;
	    fb_ang(fbs, scoops) = parm.excl_ang - parm.exp_angle;
	  end
	end
    end
end

figure(2); clf;
subplot(3,1,[1 2])
plot(scoop_h, gs_r(1,:), 'r', 'LineWidth',2);
hold on;
plot(scoop_h, gs_h(1,:), '--r', 'LineWidth',2);
plot(scoop_h, gs_r(2,:), 'g', 'LineWidth',2);
plot(scoop_h, gs_h(2,:), '--g', 'LineWidth',2);
plot(scoop_h, gs_r(3,:), 'b', 'LineWidth',2);
plot(scoop_h, gs_h(3,:), '--b', 'LineWidth',2);
ylabel('Ground Screen Size [m]')
ylim([0, 30]);
grid on;
legend('1.75 m FB, Radius','1.75 m FB, Height',...
    '2.5 m FB, Radius','2.5 m FB, Height',...
    '3.5 m FB, Radius','3.5 m FB, Height');
title('Three Shields, 45 deg min el, Compact')

subplot(3,1,3)
plot(scoop_h, exp_angle(1,:), 'r', 'LineWidth',2);
hold on;
plot(scoop_h, fb_ang(1,:), '--r', 'LineWidth',2);
plot(scoop_h, exp_angle(2,:), 'g', 'LineWidth',2);
plot(scoop_h, fb_ang(2,:), '--g', 'LineWidth',2);
plot(scoop_h, exp_angle(3,:), 'b', 'LineWidth',2);
plot(scoop_h, fb_ang(3,:), '--b', 'LineWidth',2);
ylabel('Exposure angle [deg]');
xlabel('Scoop Height [m]');
ylim([0, 60])
grid on;
legend('1.75 m FB, scoop ang','1.75 m FB, FB ang',...
    '2.5 m FB, scoop ang', '2.5 m FB, FB ang',...
    '3.5 m FB, scoop ang', '3.5 m FB, FB ang',...
    'Location','West','Numcolumns',2)
f = gcf;
f.Position = [100 100 600 600];
print(2, [figdir 'compact_scoop_45deg'],'-dpng')
