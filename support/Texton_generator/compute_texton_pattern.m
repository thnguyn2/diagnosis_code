function [texton_image]=compute_texton_pattern(filter,text_coord)
%This function computes the texton correspond to a set of filter response
%f_res and the set of coefficients text_coord
%Author: Tan H. Nguyen
%University of Illinois at Urbana-Champaign
%Inputs
%   f_res: filter responses. This is a 3D datablock where each image is a
%   PSF
%   text_coord: coordinates of the texton in the coordinate system of
%   filter responses. The purpose of this program is constructing the
%   filter that will give you necessary response
%Output(s):
%   texton_image:  3D block of texton patterns
%--------------------------------------------------------------------
    nfilters= size(filter,3);    %Number of filters responses
    nrows=size(filter,1);
    ncols=size(filter,2);
    npixels=nrows*ncols;
    H=zeros(nfilters,npixels);
    for filterIdx=1:nfilters
        curf=filter(:,:,filterIdx);
        hi =reshape(curf,[1 npixels]);
        H(filterIdx,:)=hi;
    end
    ntexton=size(text_coord,1); %Number of texton
    Hinv = H'*inv(H*H'+1e-3); %Look for the minimum norm solution
    texton_image=zeros(nrows,ncols,ntexton);
    for textonIdx=1:ntexton
        cur_vect=text_coord(textonIdx,:)';
        pattern_vect = Hinv*cur_vect;
        errorval = norm(cur_vect-H*pattern_vect,'fro');
        disp(['Error norm ' num2str(errorval)]);
        texton_image(:,:,textonIdx)=reshape(pattern_vect,nrows,ncols);
        
    end
    
end