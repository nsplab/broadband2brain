function [data_filename, info_filename] = data_file(data_set, elec, seg, cond)
% Gives filename of data .mat file (and info .mat file)

% Inputs:
% data_set - data set number (1 or 2)
% elec - electrode number (or 0 if N/A, i.e. you want info_filename)
% seg - segment (file) number (of 0 if N/A)
% cond - data conditioning type (folder name)

% Outputs:
% data_filename - filename for neural data
% info_filename - filename for info such as arm movements

if nargin < 2 || elec == 0
    chan = channels_to_use(data_set);
    elec = chan(1);
end
if nargin < 3 || seg == 0
    seg = 1;
end
if nargin < 4
    cond = 'raw';
end

data_filename = [root_dir() 'data/' cond '/data_' int2str(data_set) '_' int2str(elec) '_' int2str(seg) '.mat'];

if data_set == 1
    info_filename = [root_dir() 'data/matthew20080721.mat'];
elseif data_set == 2
    info_filename = [root_dir() 'data/nemo20080511.mat'];
elseif data_set == 3
    info_filename = [root_dir() 'data/matthew20080801.mat'];
elseif data_set == 4
    info_filename = [root_dir() 'data/nemo20080604.mat'];
elseif data_set == 5
    info_filename = [root_dir() 'data/nemo20080502.mat'];
elseif data_set == 6
    info_filename = [root_dir() 'data/nemo20080512.mat'];
elseif data_set == 7
    info_filename = [root_dir() 'data/nemo20080521.mat'];
elseif data_set == 8
    info_filename = [root_dir() 'data/nemo20080522.mat'];
elseif data_set == 9
    info_filename = [root_dir() 'data/nemo20080523.mat'];
elseif data_set == 10
    info_filename = [root_dir() 'data/nemo20080524.mat'];
elseif data_set == 11
    info_filename = [root_dir() 'data/nemo20080503.mat'];
elseif data_set == 12
    info_filename = [root_dir() 'data/nemo20080504.mat'];
elseif data_set == 13
    info_filename = [root_dir() 'data/nemo20080505.mat'];
elseif data_set == 14
    info_filename = [root_dir() 'data/nemo20080507.mat'];
elseif data_set == 15
    info_filename = [root_dir() 'data/nemo20080508.mat'];
elseif data_set == 16
    info_filename = [root_dir() 'data/nemo20080509.mat'];
elseif data_set == 17
    info_filename = [root_dir() 'data/nemo20080513.mat'];
elseif data_set == 18
    info_filename = [root_dir() 'data/nemo20080518.mat'];
elseif data_set == 19
    info_filename = [root_dir() 'data/nemo20080519.mat'];
end

end
