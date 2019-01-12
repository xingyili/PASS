function [ Ref_network_matrix,Pathway_in_Normal_location,Pathway_Fullconnected_edges,Label_pathway_fullconnected_final,...
    Label_pathway_network_final] = Ref_network_construction( Normal_sample_name,Normal_sample_value,...
    network_name,Label_pathway_fullconnected,Label_pathway_network)
[Normal_gene_number,Normal_sample_number] = size(Normal_sample_value); 
name_A = network_name(:,1);
name_B = network_name(:,2);
[PPI_row,~] = size(name_A); 
PPI_Normal_value = [];                      
Pathway_Fullconnected_edges = [];
PPI_delete = [];
Pathway_in_Normal_location = [];

for i = 1:PPI_row
    [~,PPIA_in_Normal_location,~] = intersect(Normal_sample_name,name_A(i),'rows');
    [~,PPIB_in_Normal_location,~] = intersect(Normal_sample_name,name_B(i),'rows');
    if isempty(PPIA_in_Normal_location) || isempty(PPIB_in_Normal_location)     
        PPI_delete = [PPI_delete;i];        
    else
       Pathway_in_Normal_location = [Pathway_in_Normal_location;PPIA_in_Normal_location,PPIB_in_Normal_location];
       PPIA_Normal_value = Normal_sample_value(PPIA_in_Normal_location,:);     
       PPIB_Normal_value = Normal_sample_value(PPIB_in_Normal_location,:);
       PPI_Normal_value = [PPI_Normal_value;PPIA_Normal_value,PPIB_Normal_value];
       Pathway_Fullconnected_edges = [Pathway_Fullconnected_edges;name_A(i),name_B(i)];   
    end 
end 
Label_pathway_fullconnected(PPI_delete,:) = [];
Label_pathway_fullconnected_final = Label_pathway_fullconnected;
Label_pathway_network(PPI_delete,:) = [];
Label_pathway_network_final = Label_pathway_network;
Ref_network_matrix = zeros(size(PPI_Normal_value,1),1);
for i = 1:size(PPI_Normal_value,1)
    [Pearson_C_temp,Pearson_P_temp] = corrcoef(PPI_Normal_value(i,1:Normal_sample_number),PPI_Normal_value(i,(Normal_sample_number+1):end));
    Ref_network_matrix(i) = Pearson_C_temp(2,1);
end
end

