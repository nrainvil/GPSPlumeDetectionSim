function [ snr_m ] = load_single_snr( rcvr, doy, year, varargin)
%LOAD_SINGLE_SNR 
%       snr_m   - [Time, SV, EL, AZ, L1_SNR, L2_SNR, Voltage, Temp, TS]
del_file = 0;
print_status = 1;
force_load = 0;
%Check if runnining on Windows or Linux
if(ispc)
    mat_dir = 'C:\snr_files\mat\';
else
    mat_dir = '/home/nrainvil/mat/';
end

%Hide output if requested
if nargin>=4
    print_status = varargin{1};
end
if nargin>=5
    force_load = varargin{2};
end

mat_file = strcat(num2str(rcvr),'_',num2str(year,'%02.0f'),...
                    '_',num2str(doy,'%03.0f'),'.mat');
if(print_status)
    fprintf('Loading %s\n',mat_file);
end

%% Load SNR .mat file
try
    mat_fname = strcat(mat_dir,mat_file);
    mat_dir = dir(mat_fname);
    if (~exist(mat_fname) || mat_dir.bytes < 1000 || force_load)
            %Copy mat file from Dirt
            fprintf('Loading %s from dirt.colorado.edu\n',...
                        mat_fname);
            if (ispc)
                bash_cmd = ...
                  strcat('bash -c "scp nira6794@dirt.colorado.edu:mat_',...
                     num2str(year,'%02.0f'),'/',mat_file,' /mnt/c/snr_files/mat/"');
                ps_cmd = 'powershell Get-Content .\nl.txt | ';
                system(strcat(ps_cmd,bash_cmd));
            else
                sh_cmd = ...
                   strcat('scp nira6794@dirt.colorado.edu:mat_',...
                     num2str(year,'%02.0f'),'/',mat_file,' /home/nrainvil/mat/');
                system(sh_cmd);
            end
    end
    snr_s = load(mat_fname);
    snr_m = snr_s.snr_m;
    if del_file
        if(ispc)
            system(strcat('del'," ",mat_file));
        else
            system(strcat('rm'," ",mat_file));
        end
    end
catch
    snr_m = zeros(1,9);
end

end