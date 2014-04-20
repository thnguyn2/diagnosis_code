function [label_im]=findLabelImg(org_im,gland_im)
    %Find the label image give the name of the feature image together with
    %the original, final size. Here, img is the downsampled image
    %Inputs:
    %  gland_im: an input image of the glands only
    %  org_im: an original image of the phase values
    %Outputs:
    %   label_im: the gland image
    %   ds_gland_im: downsampled gland image
    %   inner_label: label of the inner image
    %% Find the human-marking image
    org_im = imresize(org_im,[2048 2048]);
    gland_im = imresize(gland_im,[2048 2048]);
    
    gland_im=double(gland_im);
    org_im = double(org_im);
    ds_mean_gland=mean(mean(gland_im(1:30,1:30))); %Find the phase of the background region
    gland_idx=find((gland_im>1.08*ds_mean_gland)); %Make sure what we have is good for training
     %Calculate the stroma image
    mask = zeros(size(org_im));
    mask(gland_idx)=1;
    se = strel('disk',6);
    newmask = imclose(mask,se);
    gland_idx = find(newmask==1);
    
    label_im=zeros(size(gland_im));
    label_im(gland_idx)=1;
    %Next, eliminates all the pixels that has been assigned a label
    org_im(gland_idx)=mean(mean(org_im(1:10,1:10)));
    %Find the index of all glands image
    bg_mean=mean(mean(org_im(1:10,1:10))); %Find the phase of the background region
    stroma_idx=find((org_im>1.2*bg_mean)); %Make sure what we have is good for training
    mask = zeros(size(org_im));
    mask(stroma_idx)=1;
    se = strel('disk',8);
    newmask = imclose(mask,se);
    
    stroma_idx = find(newmask==1);
    label_im(stroma_idx)=2;
end
