function [fe,fo,emag,edgeloc,edir]=computeOrientedEnergyLM(filterRes,nscales,ninvscales,ndir)
%This function computes the oriented engery at different scales
%fe and fo are 3D data block where the third dimension corresponds to the
%scales. The first nscales * ndir are filter responses to the odd symmetric
%filter while the next nscales * ndirs are filter responses to the even
%symmetric filters 
%fe = filter responses to anisotropic, even filters
%fo = filter responses to anisotropic, odd filters
%edgeloc = edge map 
%gmag = oriented energy at diffent scale

    computing_platform = 0;
    nrows = size(filterRes,1);
    ncols = size(filterRes,2);
    fe = zeros(nrows,ncols,nscales);
    fo = zeros(nrows,ncols,nscales);
    edgeloc = zeros(nrows,ncols,nscales); %Edge index
    edir = zeros(nrows,ncols,nscales);
    base_size = nrows*ncols;
    id = [1:base_size]';
    tic;
    for scaleIdx=1:nscales
       curOddBlock = filterRes(:,:,(scaleIdx-1)*ndir+1:scaleIdx*ndir);
       curEvenBlock = filterRes(:,:,ndir*nscales+(scaleIdx-1)*ndir+1:ndir*nscales+scaleIdx*ndir);
       if (computing_platform==0)
           curEvenBlocksqr = curEvenBlock.^2;
           curOddBlocksqr = curOddBlock.^2;
           curenergy = curEvenBlocksqr+curOddBlocksqr;
           [max_val,maxidx]=max(curenergy,[],3);
           edir(:,:,scaleIdx)=maxidx;
           gmag = sum(curenergy,3);
           %Pick out the response coorespond to max energy accross direction of max direction
           fe(:,:,scaleIdx)=reshape(curEvenBlock(id+(maxidx(:)-1)*base_size),[nrows ncols]);
           fo(:,:,scaleIdx)=reshape(curOddBlock(id+(maxidx(:)-1)*base_size),[nrows ncols]);
           eg = (gmag>0.005*max(gmag(:))); %Find the index of all pixels that has good gradient magnitude
           temp_map = (fe(:,:,scaleIdx)>0)-(fe(:,:,scaleIdx)<0); %Convert into the binary map
           h = eg & [temp_map(2:nrows,:)~=temp_map(1:nrows-1,:);zeros(1,ncols)];
           v = eg & [temp_map(:,2:ncols)~=temp_map(:,1:ncols-1),zeros(nrows,1)];
           clear curenergy;
           clear curEvenBlocksqr;
           clear curOddBlocksqr;
           clear curOddBlock;
           clear curEvenBlocl;
           edgeloc(:,:,scaleIdx)=(h|v);
           gmag = gmag/max(gmag(:));
           emag(:,:,scaleIdx) = gmag.*edgeloc(:,:,scaleIdx); %Oriented energ
       else
           %To be implemented with opencv-mex
       end
     end
    timeelap = toc
    
   
end