function [Di,Fi]=compute_feature_statistics_cancer_class_for_em(featvect,dataclass,donorm,indivglandfeat)
    %Generate the data for em algorithm where Di is the label of all cancer
    %cases 1 for (3+3), 2 for (3+4), 4 for 4+3 and 5 for 4+4
    Di = zeros(0,1);
    Fi = cell(0,1);
    idxclass33=[find(strcmp(dataclass,'33')==1);
                find(strcmp(dataclass,'23')==1)]; 
    featurearr = [2;5];
     if (~isempty(idxclass33))
        Ni = length(idxclass33);
        Di(end+1:end+Ni,1)=1;    
        for sampleidx=1:Ni
            Fi{end++1,1}=indivglandfeat{idxclass33(sampleidx)}(:,featurearr);
        end
     end
    
    idxclass34=[find(strcmp(dataclass,'34')==1)]; 
    if (~isempty(idxclass34))
        Ni = length(idxclass34);
        Di(end+1:end+Ni,1)=2;    
        for sampleidx=1:Ni
            Fi{end+1,1}=indivglandfeat{idxclass34(sampleidx)}(:,featurearr);
        end
    end

    idxclass43=[find(strcmp(dataclass,'43')==1)]; 
    if (~isempty(idxclass43))
        Ni = length(idxclass43);
        Di(end+1:end+Ni,1)=3;    
        for sampleidx=1:Ni
            Fi{end+1,1}=indivglandfeat{idxclass43(sampleidx)}(:,featurearr);
        end
    end
    
    idxclass44=[find(strcmp(dataclass,'44')==1)]; 
    if (~isempty(idxclass44))
        Ni = length(idxclass44);
        Di(end+1:end+Ni,1)=4;    
        for sampleidx=1:Ni
            Fi{end+1,1}=indivglandfeat{idxclass44(sampleidx)}(:,featurearr);
        end
    end
end