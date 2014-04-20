function imrad = phase_im_convert(imdouble)
%This function perform a linear conversion where a double image (0 to 1 is converted into a phase image from -0.5*pi to 3.5*pi)
%Input: imdouble: an inpute image with grayscale values range from 0 to 1
%imrad: output image with phase values from -0.5 to 3.5
    min_val = -0.5;
    max_val = 3.5;
    imrad = (imdouble)*(max_val-min_val)+min_val;
end