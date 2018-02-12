function [Strain]=HalfBridgeCkt_process(data,Gf,Vi)
N_data=size(data);
Strain=zeros(N_data(1),1);
Strain(:,1)=2*data(:,1)/(Gf*Vi);
end

