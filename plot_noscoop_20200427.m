
figdir = './figs/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 45

sp.fov = 29.0;
sp.win_d = 0.72;
sp.dk_off = [0., 0.9271];
sp.el_off = [0., 2.3];
sp.az_off = [0., 0.];
sp.min_el = 45.0;
sp.n_rx = 3;

% Look at a range of forebaffle heights
fb_h1 = 0.5:0.01:5;

gs_r1 = nan(1,length(fb_h1));
gs_h1 = nan(1,length(fb_h1));
fb_ang1 = nan(1,length(fb_h1));

for ff = 1:length(fb_h1)
  sp.fb_h = fb_h1(ff);
  parm = s4_gs_study(sp,'PLOT',false,'fixwindist',0.9388);
  gs_r1(ff) = parm.gs_dim(1);
  gs_h1(ff) = parm.gs_dim(2);
  if ~isnan(parm.gs_dim(1));
    fb_ang1(ff) = parm.excl_ang;
  end
end

% Now 3-shooter ('singlestat')
fb_h3 = 1.5:0.01:5;
gs_r3 = nan(1,length(fb_h3));
gs_h3 = nan(1,length(fb_h3));
fb_ang3 = nan(1,length(fb_h3));

for ff = 1:length(fb_h3)
  sp.fb_h = fb_h3(ff);
  parm = s4_gs_study(sp,'PLOT',false,'singlestat',true,'spacing',0.1);
  gs_r3(ff) = parm.gs_dim(1);
  gs_h3(ff) = parm.gs_dim(2);
  if ~isnan(parm.gs_dim(1));
    fb_ang3(ff) = parm.excl_ang;
  end
end

figure(1); clf;
subplot(3,1,[1 2])
plot(fb_h1, gs_r1,'r','LineWidth',2);
hold on;
plot(fb_h1, gs_h1,'--r','LineWidth',2);
plot(fb_h3, gs_r3,'-b','LineWidth',2);
plot(fb_h3, gs_h3,'--b','LineWidth',2);
% Mark points of comparison
plot([1.75],[6.18],'rx','LineWidth',2,'MarkerSize',10)
plot([1.75],[19.45],'rx','LineWidth',2,'MarkerSize',10)
plot([4.0],[8.46],'bx','LineWidth',2,'MarkerSize',10)
plot([4.0],[25.41],'bx','LineWidth',2,'MarkerSize',10)
xlim([0.75, 5.0])
ylim([0, 30]);
ylabel('Ground Screen Size [m]')
title('Two Shields, 45 deg min el')
grid on;
legend('Extended GS Radius','Extended GS Height', ...
      'Compact GS Radius','Compact GS height',...
      'Location','SouthEast','Numcolumns',2);

subplot(3,1,3)
plot(fb_h1,fb_ang1, 'r','LineWidth',2);
hold on;
plot(fb_h3,fb_ang3, 'b','LineWidth',2);
plot([1.75],[56.1761],'rx','LineWidth',2,'MarkerSize',10);
plot([4.0],[55.97],'bx','LineWidth',2,'MarkerSize',10);
xlim([0.75, 5.0]);
ylim([40, 70]);
xlabel('Forebaffle Height [m]');
ylabel('FB Exposure Angle [deg]')
grid on;
legend('Extended','Compact','Location','SouthEast');
f = gcf;
f.Position = [100 100 600 600];

print(1,[figdir 'noscoop_45deg'], '-dpng')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 50

sp.fov = 29.0;
sp.win_d = 0.72;
sp.dk_off = [0., 0.9271];
sp.el_off = [0., 2.3];
sp.az_off = [0., 0.];
sp.min_el = 50.0;
sp.n_rx = 3;

% Look at a range of forebaffle heights
fb_h1 = 0.5:0.01:5;

gs_r1 = nan(1,length(fb_h1));
gs_h1 = nan(1,length(fb_h1));
fb_ang1 = nan(1,length(fb_h1));

for ff = 1:length(fb_h1)
  sp.fb_h = fb_h1(ff);
  parm = s4_gs_study(sp,'PLOT',false,'fixwindist',0.9388);
  gs_r1(ff) = parm.gs_dim(1);
  gs_h1(ff) = parm.gs_dim(2);
  if ~isnan(parm.gs_dim(1));
    fb_ang1(ff) = parm.excl_ang;
  end
end

% Now 3-shooter ('singlestat')
fb_h3 = 1.5:0.01:5;
gs_r3 = nan(1,length(fb_h3));
gs_h3 = nan(1,length(fb_h3));
fb_ang3 = nan(1,length(fb_h3));

for ff = 1:length(fb_h3)
  sp.fb_h = fb_h3(ff);
  parm = s4_gs_study(sp,'PLOT',false,'singlestat',true,'spacing',0.1);
  gs_r3(ff) = parm.gs_dim(1);
  gs_h3(ff) = parm.gs_dim(2);
  if ~isnan(parm.gs_dim(1));
    fb_ang3(ff) = parm.excl_ang;
  end
end

figure(2); clf;
subplot(3,1,[1 2])
plot(fb_h1, gs_r1,'r','LineWidth',2);
hold on;
plot(fb_h1, gs_h1,'--r','LineWidth',2);
plot(fb_h3, gs_r3,'-b','LineWidth',2);
plot(fb_h3, gs_h3,'--b','LineWidth',2);
% Mark points of comparison
plot([1.75],[5.9],'rx','LineWidth',2,'MarkerSize',10)
plot([1.75],[12.39],'rx','LineWidth',2,'MarkerSize',10)
plot([4.0],[8.1],'bx','LineWidth',2,'MarkerSize',10)
plot([4.0],[16.0],'bx','LineWidth',2,'MarkerSize',10)
xlim([0.75, 5.0])
ylim([0, 30]);
title('Two Shields, 50 deg min el')
ylabel('Ground Screen Size [m]')
grid on;
legend('Extended GS Radius','Extended GS Height', ...
      'Compact GS Radius','Compact GS height',...
      'Location','SouthEast','Numcolumns',2);

subplot(3,1,3)
plot(fb_h1,fb_ang1, 'r','LineWidth',2);
hold on;
plot(fb_h3,fb_ang3, 'b','LineWidth',2);
plot([1.75],[56.1761],'rx','LineWidth',2,'MarkerSize',10);
plot([4.0],[55.97],'bx','LineWidth',2,'MarkerSize',10);
xlim([0.75, 5.0]);
ylim([40, 70]);
xlabel('Forebaffle Height [m]');
ylabel('FB Exposure Angle [deg]')
grid on;
legend('Extended','Compact','Location','SouthEast');
f = gcf;
f.Position = [100 100 600 600];

print(2,[figdir 'noscoop_50deg'], '-dpng')
