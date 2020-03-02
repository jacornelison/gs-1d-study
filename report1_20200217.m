% Generate plots for S4 shielding study presented to SAT group 
% on 20200217

figdir = './';

% This is a 1-shooter (3 separate cryostats) that are pushed out as
% far as possible.  
% On 20200120 telecon, decided on 162.6 cm window spacing
% This implies a maximum fb radius of 81.3 cm
% Given window diameter and FOV, this implies a max height of 1.7516 m
sp.fb_h = 1.7516;
sp.fov = 29.0;
sp.win_d = 0.72;
sp.dk_off = [0., 0.9271];
sp.el_off = [0., 2.3];
sp.az_off = [0., 0.];
sp.min_el = 50.0;
sp.n_rx = 3;

parm = s4_gs_study(sp,'PLOT',true,'axis_window',20,'OUTTEXT',true)
print(1,[figdir '1_shooter_nominal'],'-dpng')

% Another point of comparison to Fred
sp.fb_h = 1.0
parm = s4_gs_study(sp,'PLOT',true,'axis_window',35,'OUTTEXT',true)
print(1,[figdir '1_shooter_1m'],'-dpng')

% Look at a range of forebaffle heights
fb_h1 = 0.5:0.01:1.75;
gs_r1 = zeros(1,length(fb_h1));
gs_h1 = zeros(1,length(fb_h1));

for ff = 1:length(fb_h1)
  sp.fb_h = fb_h1(ff);
  parm = s4_gs_study(sp,'PLOT',false);
  gs_r1(ff) = parm.gs_dim(1);
  gs_h1(ff) = parm.gs_dim(2);
end

% Do it again holding the window distance fixed 
gs_r1b = zeros(1,length(fb_h1));
gs_h1b = zeros(1,length(fb_h1));

for ff = 1:length(fb_h1)
  sp.fb_h = fb_h1(ff);
  parm = s4_gs_study(sp,'PLOT',false,'fixwindist',0.9388);
  gs_r1b(ff) = parm.gs_dim(1);
  gs_h1b(ff) = parm.gs_dim(2);
end

figure(1); clf;
plot(fb_h1, gs_r1);
hold on;
plot(fb_h1, gs_h1);
xlim([0.9, 1.75])
ylim([0, 40]);
xlabel('Forebaffle Height [m]')
ylabel('Ground Screen Size [m]')
grid on;
legend('GS Radius','GS Height');
print(1,[figdir '1_shooter_varyfb'], '-dpng')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now 3-shooter ('singlestat')
fb_h3 = 1.97:0.01:5;
gs_r3 = zeros(1,length(fb_h3));
gs_h3 = zeros(1,length(fb_h3));

for ff = 1:length(fb_h3)
  sp.fb_h = fb_h3(ff);
  parm = s4_gs_study(sp,'PLOT',false,'singlestat',true,'spacing',0.1);
  gs_r3(ff) = parm.gs_dim(1);
  gs_h3(ff) = parm.gs_dim(2);
  %print(1,[figdir '3_shooter_' num2str(ff)], '-dpng')
end

figure(1); clf;
plot(fb_h3, gs_r3);
hold on;
plot(fb_h3, gs_h3);
xlim([2.0, 5.0])
ylim([0, 40]);
xlabel('Forebaffle Height [m]')
ylabel('Ground Screen Size [m]')
grid on;
legend('GS Radius','GS Height');
print(1,[figdir '3_shooter_varyfb'], '-dpng')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare the two situations

figure(2); clf;
plot(fb_h1, gs_r1b,'r','LineWidth',2);
hold on;
plot(fb_h1, gs_h1b,'--r','LineWidth',2);
%plot(fb_h1, gs_r1b,'m','LineWidth',2);
%plot(fb_h1, gs_h1b,'--m','LineWIdth',2);
plot(fb_h3, gs_r3,'-b','LineWidth',2);
plot(fb_h3, gs_h3,'--b','LineWidth',2);
% Mark points of comparison
% (1.75, 5.93) (1.75, 12.38)
% (3.5, 7.64) (3.5, 16.77)
plot([1.75],[5.93],'rx','LineWidth',2,'MarkerSize',10)
plot([1.75],[12.38],'rx','LineWidth',2,'MarkerSize',10)
plot([3.5],[7.64],'bx','LineWidth',2,'MarkerSize',10)
plot([3.5],[16.77],'bx','LineWidth',2,'MarkerSize',10)
xlim([0.9, 5.0])
ylim([0, 30]);
xlabel('Forebaffle Height [m]')
ylabel('Ground Screen Size [m]')
grid on;
legend('Extended GS Radius','Extended GS Height', ...
      'Close-packed GS Radius','Close-packed GS height');

print(2,[figdir 'comparison'], '-dpng')
