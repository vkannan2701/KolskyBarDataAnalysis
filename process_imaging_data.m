function [filename_struc]=process_imaging_data()
%% Read datafile
path=uigetdir();
filename_struc=dir(path);
