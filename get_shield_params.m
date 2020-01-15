function sp = get_shield_params(expt_name)

% Put everything in meters. or degrees.
switch expt_name
    case "Custom"
        fb_h = 0.7366;
        fov = 20;
        win_d = 0.31;
        dk_off = [-1.1, 1.5964];
        el_off = [0., 1.1750];
        az_off = [0., 0.];
        min_el = 50;
        n_rx = 5;
        
    case "keck"
        % Forebaffle dimensions from MSG's posting
        % Optical params from Lorenzo's instrument parameter spreadsheet
        % Mount values from Colin's pointing model defaults
        fb_h = 0.7366+0.14; %plus dist to aperture
        fov = 9.6*2;
        win_d = 0.264; % aperture diameter
        dk_off = [-1.1, 1.5964];
        el_off = [0., 1.1750];
        az_off = [0., 0.];
        min_el = 46;
        n_rx = 5;
        gs_dim = [7.62 5.1]; % DASI groundshield
        
    case "BA"
        % Taken from SOLIDWorks "AsBuild" models in the repo.
        % FB dimensions from Fig. 3.7 of NWP's 20190207_BA_Forebaffle_Concepts posting in BA logbook.
        % Optical params from Lorenzo's instrument parameter spreadsheet
        fb_h = 1.086;
        fov = 29.6;
        win_d = 0.34*2;
        dk_off = [0., 0.9271];
        el_off = [0., 2.3];
        az_off = [0., 0.];
        min_el = 45;
        n_rx = 4;
        gs_dim = [7.62 5.1]; % DASI groundshield
        
    
    case "BICEP3"
        % Forebaffle dimensions from Mike Crumrin's designs on the TWiki:
        % http://polar-array.stanford.edu/twiki/bin/view/BICEP3/Logbook140825SplitForebaffle
        % Optical params from Lorenzo's instrument parameter spreadsheet
        % Mount values from Solidworks drawings
        fb_h = 1.5;
        fov = 14.1*2;
        win_d = 0.5776;
        dk_off = [0., 1.335-0.1915];
        el_off = [0., 0.0];
        az_off = [0., 0.];
        min_el = 54;
        n_rx = 1;
        gs_dim = [5.3 3];
        g3 = 0.5;
        g1 = 0.49;
        g2 = 0.01;
        gd = 0;
        
    case "BICEP1"
        % mount, fb, and gs dimensions pulled from paper vertex drawings...
        fb_h = 0.9689;
        fov = 9.68*2;
        win_d = 0.264;
        dk_off = [0., .9772];
        el_off = [0., 0.];
        az_off = [0., 0.];
        min_el = 50;
        n_rx = 1;
        gs_dim = [3.9878 2.1749];
        g3 = 0.5;
        g1 = 0.49;
        g2 = 0.01;
        gd = 0;
        
    case "S4SAT3SHOT"
        % Update Params
        fb_h = 0.7366;
        fov = 20;
        win_d = 0.31;
        dk_off = [-1.1, 1.5964];
        el_off = [0., 1.1750];
        az_off = [0., 0.];
        min_el = 50;
        n_rx = 3;
    otherwise
        
end

% Add all of the variable to a struct and return that.
vars = whos;
sp = [];
for i = 1:length(vars)
eval(['sp.(''' vars(i).name ''')=' vars(i).name ';'])
end
