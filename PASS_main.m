clear all
clc

dataPath = './data/';
load([dataPath,'/Genes_expression_data.mat'])
Genes_ID = Genes_expression_data(2:end,1);
Label = Genes_expression_data(1,2:end);  
[~,n1] = find(Label==0);
[~,n2] = find(Label==-1);
[~,n3] = find(Label==1);
label1 = Label(Label == -1);
label2 = Label(Label == 1);
NormalLabel = Label(Label == 0);
DiseaseLabel = [label1,label2];
Normal_sample_value = Genes_expression_data(2:end,n1+1);  
Normal_sample_name = Genes_ID;
Normal_number = size(Normal_sample_value,2);
Disease_sample_valueA = Genes_expression_data(2:end,n2+1);
Disease_sample_valueB = Genes_expression_data(2:end,n3+1);
Disease_sample_value = [Disease_sample_valueA,Disease_sample_valueB];
Disease_number = size(Disease_sample_value,2);


%% *********************Constructing pathway full connected network***********

AllFile=dir([dataPath,'/Pathway_data/*.txt']);                                            
AllFileNum = length(AllFile);                                                    
FullConnectedNetwork = cell(2,1);                                               

for i = 1:AllFileNum
    Temp = textread([dataPath,'/Pathway_data/',AllFile(i).name]);                     
    TempNode = [Temp(:,1);Temp(:,2)];                                            
    UniqueNode = unique(TempNode);                                             
    UniqueNode = intersect(UniqueNode,Genes_ID);         
    NodeSize = length(UniqueNode);                                             
    FullConnectedNetworkTemp = [];                                            
  if NodeSize ==1 || NodeSize ==0
      i
      continue;
  end
    for j = 1:NodeSize-1
        NodeID = UniqueNode(j);                                                
        NetworkA = repmat(NodeID,(NodeSize-j),1);                               
        NetworkB = UniqueNode((j+1):end) ;                                       
        FullConnectedNetworkTemp = [FullConnectedNetworkTemp;NetworkA,NetworkB]; 
    end
    [~,~,ifullconnected_1] = intersect(Temp,FullConnectedNetworkTemp,'rows');  
    Temp2 = [Temp(:,2),Temp(:,1)];                                              
    [~,~,ifullconnected_2] = intersect(Temp2,FullConnectedNetworkTemp,'rows');
    c = zeros(size(FullConnectedNetworkTemp,1),1);
    c(ifullconnected_1) = 1;
    c(ifullconnected_2) = 1;                                                 
    ILabel = ones(size(c,1),1)*i;                                            
    if sum(c) == 0
        i
        continue;
    end      
    FullConnectedNetworkTemp = [FullConnectedNetworkTemp,c,ILabel];
    FullConnectedNetwork = {FullConnectedNetwork{1,:},FullConnectedNetworkTemp; ... 
                               FullConnectedNetwork{2,:},AllFile(i).name(1:8)};
                                             
  
end
 FullConnectedNetwork(:,1) = [];  
%% **************************AUCpath for nodes****************************************
NumberOfNetwork = size(FullConnectedNetwork,2);

AllNodeInPathway = [];
All_node = [];
NodeAndLabelTemp = [];
for n = 1:NumberOfNetwork
    network_name_temp = FullConnectedNetwork{1,n};
    network_name = network_name_temp(:,1:2);
    Node_temp = [network_name(:,1);network_name(:,2)];
    Unique_Node_temp = unique(Node_temp);
    AllNodeInPathway = [AllNodeInPathway;Unique_Node_temp];  
end
AllNodeInPathway = unique(AllNodeInPathway);

for m = 1:NumberOfNetwork
    network_name_temp = FullConnectedNetwork{1,m};
    network_name = network_name_temp(:,1:2);
    Node_temp = [network_name(:,1);network_name(:,2)];
    Unique_Node_temp = unique(Node_temp);
    [~,~,itemp1] = intersect(AllNodeInPathway,Genes_ID,'rows');
    [itemp1_number,~] = size(itemp1);
    NodeExpressionValue = zeros(itemp1_number,Disease_number);
   for j = 1:Disease_number
       for k = 1:itemp1_number
       NodeExpressionValue(k,j) = Disease_sample_value(itemp1(k,1),j)./mean(Normal_sample_value(itemp1(k,1),:));
       end
   end
    [~,~,itemp2] = intersect(Unique_Node_temp,AllNodeInPathway,'rows');
    C = zeros(size(AllNodeInPathway,1),1);
    C(itemp2) = 1;
    NodeLabel = ones(size(AllNodeInPathway,1),1)*m;
    NodeAndLabelTemp = [NodeLabel,C,AllNodeInPathway,NodeExpressionValue]; 
    All_node = [All_node;NodeAndLabelTemp];      
end
path_mat1 = CalAUCPathNode(All_node);
 
 %% *************************************AUCpath for edges****************************************  
All_network = [];
NumberOfNetwork = size(FullConnectedNetwork,2);
for j = 1: NumberOfNetwork        
    j
    network_name_temp = FullConnectedNetwork{1,j};
    network_name = network_name_temp(:,1:2);
    Label_pathway_fullconnected = network_name_temp(:,3);     
    Label_pathway_network = network_name_temp(:,4);          
[Ref_network_matrix,Pathway_in_Normal_location,Pathway_Fullconnected_edges,Label_pathway_fullconnected_final,...
    Label_pathway_network_final] = Ref_network_construction(Normal_sample_name,Normal_sample_value,network_name,...
    Label_pathway_fullconnected,Label_pathway_network);
Every_network = zeros(size(Pathway_in_Normal_location,1),Disease_number);

for i = 1:Disease_number
    Sample_value = [Normal_sample_value,Disease_sample_value(:,i)];
    Perturbed_network_matrix = Perturbed_network_construction(Sample_value,Pathway_in_Normal_location);
    Dif_network_matrix = Perturbed_network_matrix - Ref_network_matrix;
    Dif_network_zscore = Dif_network_matrix./((1-Ref_network_matrix.^2)/(Normal_number-1));
    
    Every_network(:,i) = Dif_network_zscore; 
    
end

All_network = [All_network;Label_pathway_network_final,Label_pathway_fullconnected_final,Pathway_Fullconnected_edges,Every_network];

end

path_mat2 = CalAUCPath(All_network);
path_mat = [path_mat1;path_mat2];
PASSresult = [DiseaseLabel;path_mat];

save PASSresult.mat PASSresult
