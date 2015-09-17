function H=convert2hsv(phasemap,newmap,maxproj,minproj)
%Overlaid the heat map and the phase image. The staturation will be set to
%1. The heat map will be used to choose the hue channel. The phasemap will
%be used to choose the bright ness
       %Create an hsv image
       H = zeros(size(newmap,1),size(newmap,2),3);
       H(:,:,1)=(newmap-minproj)/(maxproj - minproj);
       H(:,:,2)=1; %Value
       H(:,:,3)=phasemap; %Value
       H =im2uint16(hsv2rgb(H));

end