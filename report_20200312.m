
figdir = './figs/';

%% Case 1: Extended configuration

sp.fb_h = 1.7516;
sp.fov = 29.0;
sp.win_d = 0.72;
sp.dk_off = [0., 0.9271];
sp.el_off = [0., 2.3];
sp.az_off = [0., 0.];
sp.min_el = 50.0;
sp.n_rx = 3;
%sp.gs_dim = [0,0];
%sp = rmfield(sp,'gs_dim');

parm = s4_gs_study(sp,'PLOT',true,'OUTTEXT',true,'axis_window',15,'fixwindist',0.9388,'ts_dim',false);
print(1,[figdir 'SAT_3RX_extended_noscoop'],'-dpng')

for scoops = 1:3
    parm = s4_gs_study(sp,'PLOT',true,'OUTTEXT',true,'axis_window',15,'fixwindist',0.9388,'ts_dim',false,'threeshield',scoops);
    printname = [figdir sprintf('SAT_3RX_extended_scoop_%i',scoops)];
    print(1,printname,'-dpng')
end

%% Case 2: Compact Configuration
sp.fov = 29.0;
sp.win_d = 0.72;
sp.dk_off = [0., 0.9271];
sp.el_off = [0., 2.3];
sp.az_off = [0., 0.];
sp.min_el = 50.0;
sp.n_rx = 3;
%sp.gs_dim = [0,0];

for fbs = [1.75, 3.0, 3.5]
    sp.fb_h = fbs;
    parm = s4_gs_study(sp,'PLOT',true,'OUTTEXT',true,'axis_window',15,'spacing',0.1,'singlestat',true,'ts_dim',false);
    printname = [figdir sprintf('SAT_3RX_compact_fb_%i_noscoop',ceil(fbs))];
    print(1,printname,'-dpng')
    for scoops = 1:3
        parm = s4_gs_study(sp,'PLOT',true,'OUTTEXT',true,'axis_window',15,'spacing',0.1,'singlestat',true,'ts_dim',false,'threeshield',scoops);
        printname = [figdir sprintf('SAT_3RX_compact_fb_%i_scoop_%i',ceil(fbs),scoops)];
        print(1,printname,'-dpng')
    end
end


