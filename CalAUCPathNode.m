function [ path_mat ] = CalAUCPathNode( single_pathway )
    Ids=unique(single_pathway(:,1)); 
    [m,n]=size(single_pathway);
    path_mat=zeros(length(Ids),n-3);
    for i=1:length(Ids)
        tid=Ids(i);
        tmpmat=single_pathway(single_pathway(:,1)==tid,:);
        disp(i);
            for j=4:n
                 [~,~,~,AUC] = perfcurve(tmpmat(:,2),tmpmat(:,j),1);
                 path_mat(i,j-3)=AUC;
            end
    end

end

