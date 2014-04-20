function [filenames,glandnames]=findFileName(dirpath)
    %This file search for all the file names of the glands and the raw data
    %Author: Tan H. Nguyen
    %University of Illinois at Urbana-Champaign
    %Inputs
    %   dirpath: path to the searching folder
    %Outputs:
    %   filenames, glandnames: 2 cells containing the file names of the
    %   data
    %------------------------------------------------------------
    filenames = cell(4,1);
    glandnames = cell(4,1);
    %Search for all files of class 3+3
    glandlistname = dir(strcat(dirpath,'33*Glands.tif'));
    curfilenames = cell(0,1);
    curglandnames = cell(0,1);
    nfiles = length(glandlistname);
    for filenameidx=1:nfiles
        curglandnames{end+1,1} =strcat(dirpath,glandlistname(filenameidx,1).name);
        tempname = glandlistname(filenameidx,1).name; %Crop out the filename only
        curfilenames{end+1,1}=strcat(dirpath,tempname(1:end-10),tempname(end-3:end));
    end
    filenames{1,1}=curfilenames;
    glandnames{1,1}=curglandnames;
    
    %3+4 class
    glandlistname = dir(strcat(dirpath,'34*Glands.tif'));
    curfilenames = cell(0,1);
    curglandnames = cell(0,1);
    nfiles = length(glandlistname);
    for filenameidx=1:nfiles
        curglandnames{end+1,1} =strcat(dirpath,glandlistname(filenameidx,1).name);
        tempname = glandlistname(filenameidx,1).name; %Crop out the filename only
        curfilenames{end+1,1}=strcat(dirpath,tempname(1:end-10),tempname(end-3:end));
    end
    filenames{2,1}=curfilenames;
    glandnames{2,1}=curglandnames;
    
    %4+3 class
    glandlistname = dir(strcat(dirpath,'43*Glands.tif'));
    curfilenames = cell(0,1);
    curglandnames = cell(0,1);
    nfiles = length(glandlistname);
    for filenameidx=1:nfiles
        curglandnames{end+1,1} =strcat(dirpath,glandlistname(filenameidx,1).name);
        tempname = glandlistname(filenameidx,1).name; %Crop out the filename only
        curfilenames{end+1,1}=strcat(dirpath,tempname(1:end-10),tempname(end-3:end));
    end
    filenames{3,1}=curfilenames;
    glandnames{3,1}=curglandnames;
    
    %4+4 class
    glandlistname = dir(strcat(dirpath,'44*Glands.tif'));
    curfilenames = cell(0,1);
    curglandnames = cell(0,1);
    nfiles = length(glandlistname);
    for filenameidx=1:nfiles
        curglandnames{end+1,1} =strcat(dirpath,glandlistname(filenameidx,1).name);
        tempname = glandlistname(filenameidx,1).name; %Crop out the filename only
        curfilenames{end+1,1}=strcat(dirpath,tempname(1:end-10),tempname(end-3:end));
    end
    filenames{4,1}=curfilenames;
    glandnames{4,1}=curglandnames;

end

