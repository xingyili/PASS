function [ path_mat ] = CalAUCPath( single_pathway )
    Ids=unique(single_pathway(:,1)); 
    [m,n]=size(single_pathway);
    path_mat=zeros(length(Ids),n-4);
    for i=1:length(Ids)
        tid=Ids(i);
        tmpmat=single_pathway(single_pathway(:,1)==tid,:);
        Tmp_label = single_pathway(single_pathway(:,1)==tid,2);
        value_of_Tmplabel = sum(Tmp_label);
        [x,y]=size(tmpmat);
        disp(i);
        if x==1                                      
            path_mat(i,:)=tmpmat(1,5:n);
        elseif  value_of_Tmplabel == x               
            for k = 5:n
            path_mat(i,k-4)=mean(tmpmat(:,k));
            end
        else
            for j=5:n
                 [~,~,~,AUC] = perfcurve(tmpmat(:,2),tmpmat(:,j),1);
                 path_mat(i,j-4)=AUC;
            end
        end
    end

end

