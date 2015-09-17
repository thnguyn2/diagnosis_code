clc;
clear all;
close all;
load texton_data.mat;

ntexton_per_im=100;
nfiles = 67;
prior = zeros(ntexton_per_im*nfiles,3);
for fileidx=1:nfiles
    prior((fileidx-1)*ntexton_per_im+1:fileidx*ntexton_per_im,:)=...
        repmat(prior_over_images(fileidx,:),[ntexton_per_im 1]);
end
figure(1);
imagesc(prior);
colorbar;
title('Prior');

%Compute the probability of texton;
texton_prob = sum(hist_over_images.*prior,2);

%texton_set = unique(texton_map_over_images,'rows');

%Now, get rid of all textons that has the using frequency less than the
%average
prob_thresh = 1/ntexton_per_im;

goodtextonidx = find(texton_prob>prob_thresh);
hist_over_images = hist_over_images(goodtextonidx,:);
texton_map_over_images = texton_map_over_images(goodtextonidx,:);
%Normalization w.r.t energy
texton_norm = 1./sqrt(sum(texton_map_over_images.^2,2));

%Normalization w.r.t a feature
%texton_norm = 1./texton_map_over_images(:,6);

length_norm = repmat(texton_norm,[1 20]);
texton_map_over_images =texton_map_over_images.*length_norm; 



prior = prior(goodtextonidx,:);
texton_prob = texton_prob(goodtextonidx,:);

%Next, compute the joint distribution
%p(texton = i| class = j)*p(class=j) = p(texton = i, class = j)
joint_prob = hist_over_images.*prior;
%By Bayes law,p(texton = i, class = j) = p(class= j|texton = i)*p(texton = i)
%Normalize each row and choose the class one that gives p(class=j|texton =
%i)>0.7*(sum over k {p(class = k|texton = i)})
normalizing_fact = sum(joint_prob,2);
joint_prob = joint_prob./repmat(normalizing_fact,[1 3]);

all_texton_map_over_images = texton_map_over_images;
all_prior = prior;
all_joint_prob = joint_prob;
all_texton_prob = texton_prob;
all_hist_over_images = hist_over_images;

%Second reduction by choosing the textons that have highest discrimination power
maxval=max(joint_prob,[],2);
retainedidx = find(maxval>0.75);
clear maxval;
clear normalizing_fact;
clear goodtextonidx;

joint_prob = joint_prob(retainedidx,:);
hist_over_images = hist_over_images(retainedidx,:);
texton_map_over_images = texton_map_over_images(retainedidx,:);
prior = prior(retainedidx,:);
texton_prob = texton_prob(retainedidx,:);
[maxval,classidx]=max(joint_prob,[],2);
%Now pick out all texton of a specific glass
lumenidx = find(classidx==1);
glandidx = find(classidx==2);
stromaidx = find(classidx==3);

figure(1);
imagesc(joint_prob);
title('Joint probability (normalized)'); 
colorbar
%Find all textons in common
gland_textons = texton_map_over_images(glandidx,:);
stroma_textons = texton_map_over_images(stromaidx,:);

figure(2);
for idx=1:size(stroma_textons,1)
    plot(stroma_textons(idx,:),'g');
    hold on;
end
for idx=1:size(gland_textons,2)
    plot(gland_textons(idx,:),'r');
    hold on;
end
title('green: stroma, r: glands')
%Now, perform FLD on our data set and see the result


gland_textons = gland_textons';
stroma_textons = stroma_textons';
%Now, doing the FLD to find the best low-dimensionality reduction
nglands = size(gland_textons,2);
nstromas = size(stroma_textons,2);
mugland = sum(gland_textons,2)/nglands;
mustroma = sum(stroma_textons,2)/nstromas;
mu = sum([gland_textons stroma_textons],2)/(nglands+nstromas);
figure(3)
plot(mugland)
hold on
plot(mustroma,'r')
Sglands = (gland_textons - repmat(mugland,[1 nglands]))*(gland_textons - repmat(mugland,[1 nglands]))';
Sstromas = (stroma_textons - repmat(mustroma,[1 nstromas]))*(stroma_textons - repmat(mustroma,[1 nstromas]))';
Sw=Sglands+Sstromas;
Sb = nglands*(mugland-mu)*(mugland-mu)'+nstromas*(mustroma-mu)*(mustroma-mu)';
[V,D]=eigs(Sb,Sw,2);

%Now, take only the first data and do the projections
nfeatures = 2;
W=real(V(:,1:nfeatures));
glandproj = W'*gland_textons;
stromaproj = W'*stroma_textons;

if (nfeatures==1)
    figure(4);
    plot(glandproj');
    hold on;
    plot(stromaproj','r');
end

%Plot the distribution of two projections
if (nfeatures==2)
    figure(4);
    for glandidx=1:nglands
        plot(glandproj(1,glandidx),glandproj(2,glandidx),'+r');
        hold on;
    end
    for stromaidx=1:nstromas
        plot(stromaproj(1,stromaidx),stromaproj(2,stromaidx),'ob');
        hold on;
    end
    
end
threshval = -0.3;

%Now, verify the stability of our finding by finding all the texton
[maxval,maxidx]=max(all_joint_prob,[],2);
testidx = find(maxidx~=1);
new_joint_prob = all_joint_prob(testidx,:);
figure(5);
imagesc(new_joint_prob);
title('Probabilty of 2 classes');
colorbar;
new_texton_map_over_images = all_texton_map_over_images(testidx,:);
new_texton_map_over_images = new_texton_map_over_images';
proj_data=W'*new_texton_map_over_images;
glandclassidx = find(proj_data<threshval);
groundtruthgland = find(maxidx==2); %Find all the texton that support gland
disp(['Matching ratio: ' num2str(length(intersect(groundtruthgland,glandclassidx))/length(glandclassidx))]);

