function [Dif_network_matrix] = Perturbed_network_construction(Sample_value,Pathway_in_Normal_location)
Sample_number = size(Sample_value,2);
PPI_Dif_value = [Sample_value(Pathway_in_Normal_location(:,1),1:Sample_number),Sample_value(Pathway_in_Normal_location(:,2),1:Sample_number)];

Dif_network_matrix = zeros(size(PPI_Dif_value,1),1);
for i = 1:size(PPI_Dif_value,1)
    [Pearson_C_temp,~] = corrcoef(PPI_Dif_value(i,1:Sample_number),PPI_Dif_value(i,(Sample_number+1):end));
    Dif_network_matrix(i) = Pearson_C_temp(2,1);
end
end


