function [idx,c]=kmedoids(X,k)
%k-medoid clustering the data into k-partitions that minimize the pairwise
%dissimilarity
%Inputs:
%   X: a data set to be partitioned where each column is a data vector
%   k: number of clusters
%Outputs:
%   idx: an array telling the label of each data vector. This is an
%   integer number between 1 and k
%   c: coordinate of the centroid
%Note that the algorithm is robust to noise and outliers. Furthermore, it
%is more beneficial than k-mean when we want to characterize each single
%texture
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
% [idx,c]=kmedoids(data,5);
% hold on;
% plot(c(1,:),c(2,:),'*g');
% color_arr='rgbmky';
% for cluster_idx=1:5
%     cur_idx=find(idx==cluster_idx);
%     plot(data(1,cur_idx),data(2,cur_idx),strcat('+',color_arr(mod(cluster_idx,6)+1)));
% end

%% ----------------------------------------------------------------------
disp('Applying k-medoid algorithm....');
%Step 1: initialize k cluster centroids randomly
    nsamples=size(X,2); %Number of samples
    idx_arr=randperm(nsamples);
    c=X(:,idx_arr(1:k)); %Initialize the centroids

%Step 2: Assign every point to the nearest centroids and updating the
%centroids
    cur_rmse=1e+14;
    iter = 1;
    while (1)
        distance_arr=zeros(k,nsamples);
        for cluster_idx=1:k
             error_arr=X-repmat(c(:,cluster_idx),1,nsamples);%Compute the distance from every point to its centroid
             distance_arr(cluster_idx,:)=sqrt(sum(error_arr.^2,1));         
        end
        [min_dist,idx]=min(distance_arr,[],1);
        new_rmse=sqrt(sum(min_dist.^2,2)/nsamples);
        disp(['#Iter: ' num2str(iter), ' RMSE= ', num2str(new_rmse)]);
        %Check the stopping condition
        if ((abs((cur_rmse-new_rmse))/(cur_rmse+1e-6))<0.01)
            break;
        else
            cur_rmse=new_rmse;
        end
        
        %Update the centroid by searching over the neiboring points of each
        %cluster center.
        for cluster_idx=1:k
            best_inner_rmse=cur_rmse;
            best_candidate=-1; %%This is the index for the best candidate to swap the mediod
            for candidate_idx=1:nsamples %Go through every other candidate and pick the one that gives smallest error
                new_cluster_center_candidate=X(:,candidate_idx);
                %Compute the distance of every point to the new candidate
                temp_c=c;
                temp_c(:,cluster_idx)=new_cluster_center_candidate;
                %Calculate the error for the new setup
                for inner_cluster_idx=1:k
                    temp_error_arr=X-repmat(temp_c(:,inner_cluster_idx),1,nsamples);%Compute the distance from every point to its centroid
                    temp_distance_arr(inner_cluster_idx,:)=sqrt(sum(temp_error_arr.^2,1));         
                end
                [temp_min_dist,new_idx]=min(temp_distance_arr,[],1);
                temp_rmse=sqrt(sum(temp_min_dist.^2,2)/nsamples);
                if (temp_rmse<best_inner_rmse)
                    best_inner_rmse=temp_rmse;
                    best_candidate=candidate_idx;
                end
            end
            if (best_candidate~=-1)
                c(:,cluster_idx)=X(:,best_candidate);%if a better candidate found, then swap
                disp(['Best temp RMSE = ', num2str(temp_rmse)]);
            end
        end
        iter=iter+1;
    end
end

