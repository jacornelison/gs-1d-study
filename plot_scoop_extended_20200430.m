
figdir = './figs/';

%% Extended Configuration
sp.fov = 29.0;
sp.win_d = 0.72;
sp.dk_off = [0., 0.9271];
sp.el_off = [0., 2.3];
sp.az_off = [0., 0.];
sp.min_el = 50.0;
sp.n_rx = 3;
%sp.gs_dim = [0,0];

scoop_h = 0:0.01:3;
fb_h = [1.25, 1.5, 1.75]%, 2.5, 3.5];
gs_r = nan(length(fb_h), length(scoop_h));
gs_h = nan(length(fb_h), length(scoop_h));
fb_ang = nan(length(fb_h), length(scoop_h));
exp_angle = nan(length(fb_h), length(scoop_h));

for fbs = 1:length(fb_h)
    sp.fb_h = fb_h(fbs);
    
    for scoops = 1:length(scoop_h)
        parm = s4_gs_study(sp,'PLOT',false,'OUTTEXT',true,'axis_window',15,'fixwindist',0.9388,'ts_dim',false,'threeshield',scoop_h(scoops));

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
plot(scoop_h, gs_r(1,:), 'm', 'LineWidth',2);
hold on;
plot(scoop_h, gs_h(1,:), '--m', 'LineWidth',2);
plot(scoop_h, gs_r(2,:), 'c', 'LineWidth',2);
plot(scoop_h, gs_h(2,:), '--c', 'LineWidth',2);
plot(scoop_h, gs_r(3,:), 'r', 'LineWidth',2);
plot(scoop_h, gs_h(3,:), '--r', 'LineWidth',2);
ylabel('Ground Screen Size [m]')
ylim([0, 30]);
grid on;
legend('1.25 m FB, Radius','1.25 m FB, Height',...
    '1.5 m FB, Radius','1.5 m FB, Height',...
    '1.75 m FB, Radius','1.75 m FB, Height',...
    'Location','NorthEast');
title('Three Shields, 50 deg min el, Extended')

subplot(3,1,3)
plot(scoop_h, exp_angle(1,:), 'm', 'LineWidth',2);
hold on;
plot(scoop_h, fb_ang(1,:), '--m', 'LineWidth',2);
plot(scoop_h, exp_angle(2,:), 'c', 'LineWidth',2);
plot(scoop_h, fb_ang(2,:), '--c', 'LineWidth',2);
plot(scoop_h, exp_angle(3,:), 'r', 'LineWidth',2);
plot(scoop_h, fb_ang(3,:), '--r', 'LineWidth',2);
ylabel('Exposure angle [deg]');
xlabel('Scoop Height [m]');
ylim([0, 60]);
grid on;
legend('1.25 m FB, scoop ang','1.25 m FB, FB ang',...
    '1.5 m FB, scoop ang', '1.5 m FB, FB ang',...
    '1.75 m FB, scoop ang', '1.75 m FB, FB ang',...
    'Location','West','Numcolumns',2)


f = gcf;
f.Position = [100 100 600 600];
print(1, [figdir 'extended_scoop_50deg'],'-dpng')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extended Configuration
sp.fov = 29.0;
sp.win_d = 0.72;
sp.dk_off = [0., 0.9271];
sp.el_off = [0., 2.3];
sp.az_off = [0., 0.];
sp.min_el = 45.0;
sp.n_rx = 3;
%sp.gs_dim = [0,0];

scoop_h = 0:0.01:3;
fb_h = [1.25, 1.5, 1.75]%, 2.5, 3.5];
gs_r = nan(length(fb_h), length(scoop_h));
gs_h = nan(length(fb_h), length(scoop_h));
fb_ang = nan(length(fb_h), length(scoop_h));
exp_angle = nan(length(fb_h), length(scoop_h));

for fbs = 1:length(fb_h)
    sp.fb_h = fb_h(fbs);
    
    for scoops = 1:length(scoop_h)
        parm = s4_gs_study(sp,'PLOT',false,'OUTTEXT',true,'axis_window',15,'fixwindist',0.9388,'ts_dim',false,'threeshield',scoop_h(scoops));

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
plot(scoop_h, gs_r(1,:), 'm', 'LineWidth',2);
hold on;
plot(scoop_h, gs_h(1,:), '--m', 'LineWidth',2);
plot(scoop_h, gs_r(2,:), 'c', 'LineWidth',2);
plot(scoop_h, gs_h(2,:), '--c', 'LineWidth',2);
plot(scoop_h, gs_r(3,:), 'r', 'LineWidth',2);
plot(scoop_h, gs_h(3,:), '--r', 'LineWidth',2);
ylabel('Ground Screen Size [m]')
ylim([0, 30]);
grid on;
legend('1.25 m FB, Radius','1.25 m FB, Height',...
    '1.5 m FB, Radius','1.5 m FB, Height',...
    '1.75 m FB, Radius','1.75 m FB, Height',...
    'Location','NorthEast');
title('Three Shields, 45 deg min el, Extended')

subplot(3,1,3)
plot(scoop_h, exp_angle(1,:), 'm', 'LineWidth',2);
hold on;
plot(scoop_h, fb_ang(1,:), '--m', 'LineWidth',2);
plot(scoop_h, exp_angle(2,:), 'c', 'LineWidth',2);
plot(scoop_h, fb_ang(2,:), '--c', 'LineWidth',2);
plot(scoop_h, exp_angle(3,:), 'r', 'LineWidth',2);
plot(scoop_h, fb_ang(3,:), '--r', 'LineWidth',2);
ylabel('Exposure angle [deg]');
xlabel('Scoop Height [m]');
ylim([0, 60]);
grid on;
legend('1.25 m FB, scoop ang','1.25 m FB, FB ang',...
    '1.5 m FB, scoop ang', '1.5 m FB, FB ang',...
    '1.75 m FB, scoop ang', '1.75 m FB, FB ang',...
    'Location','West','Numcolumns',2)


f = gcf;
f.Position = [100 100 600 600];
print(2, [figdir 'extended_scoop_45deg'],'-dpng')
