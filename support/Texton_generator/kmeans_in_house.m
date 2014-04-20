function [idx,c]=kmeans_in_house(X,k)
%k-means clustering the data into k-partitions
%Inputs:
%   X: a data set to be partitioned where each row is a data vector
%   k: number of clusters
%Outputs:
%   idx: an array telling the label of each data vector. This is an
%   integer number between 1 and k
%   c: coordinate of the centroid
%Note that the algorithm is not robust to noise and outliers.
%----------------------------------------------------------------------
% Examples:
% %Generate testing data
% x=[24;34;-4;5;9];
% y=[2;4;6;9;2];
% figure
% plot(x,y,'or');
% %Generate data point around each point
% ndata_point=20;
% var=0.1;
% data=zeros(2,0);
% for i=1:length(x)
%     new_data_x=x(i)+var*randn(1,ndata_point);
%     new_data_y=y(i)+var*randn(1,ndata_point);
%     new_data=[new_data_x;new_data_y];
%     data(:,end+1:end+ndata_point)=new_data;
% end
% hold on
% plot(data(1,:),data(2,:),'+b');
% [idx,c]=kmeans(data,5);
% hold on;
% plot(c(1,:),c(2,:),'*g');
% color_arr='rgbmky';
% for cluster_idx=1:5
%     cur_idx=find(idx==cluster_idx);
%     plot(data(1,cur_idx),data(2,cur_idx),strcat('+',color_arr(mod(cluster_idx,6)+1)));
% end

%% ----------------------------------------------------------------------
disp('Applying k-means algorithm....');
%Step 1: initialize k cluster centroids randomly
    nsamples=size(X,1); %Number of samples
    idx_arr=randperm(nsamples);
    c=X(idx_arr(1:k),:); %Initialize the centroids
%Step 2: Assign every point to the nearest centroids and updating the
%centroids
    cur_rmse=1e+14;
    iter = 1;
    while (1)
        distance_arr=zeros(k,nsamples);
        for cluster_idx=1:k
             error_arr=X-repmat(c(cluster_idx,:),nsamples,1);%Compute the distance from every point to its centroid
             distance_arr(cluster_idx,:)=sqrt(sum(error_arr.^2,2));         
        end
        [min_dist,idx]=min(distance_arr,[],1);
        new_rmse=sqrt(sum(min_dist.^2,2)/nsamples);
        disp(['#Iter: ' num2str(iter), ' RMSE= ', num2str(new_rmse)]);
        %Check the stopping condition
        if ((abs((cur_rmse-new_rmse))/cur_rmse)<0.01)
            break;
        else
            cur_rmse=new_rmse;
        end
        %Update the centroid
        for cluster_idx=1:k
            cur_idx=find(idx==cluster_idx); %Find all point that is currently assigned to cluster_idx
            c(cluster_idx,:)=mean(X(cur_idx,:),1);
        end
        iter=iter+1;
    end
end

