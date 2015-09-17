function classnum=cell2number(curcell)
    %Convert the cell to classnumber '+1'-> 1, '-1'>-1
    nsamples=size(curcell,1);
    classnum = zeros(size(curcell));
    for sampleidx=1:nsamples
        classnum(sampleidx,1)=str2double(cell2mat(curcell(sampleidx,1)));
    end
end
