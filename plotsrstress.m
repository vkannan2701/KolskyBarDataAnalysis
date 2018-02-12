function []=plotsrstress(stress_time,strain_time,image_time,sr_time)
clc;
close all;
%% 
N_imgarr=size(image_time);
N_img=N_imgarr(1);
N_stressarr=size(stress_time);
N_stress=N_stressarr(1);
N_strainarr=size(strain_time);
N_strain=N_strainarr(1);
N_strainrarr=size(sr_time);
N_strainr=N_strainrarr(1);
%%
imgnum=0;
timetemp=image_time(1,1);
for i=2:N_img
    if image_time(i,2)>=1.5
        if image_time(i,1)-timetemp>=150e-9
            imgnum=imgnum+1;
            timetemp=image_time(i,1);
        end
    end
end
%%
SCMg=zeros(imgnum,5);
timetemp=image_time(1,1);
imnum=0;
for i=2:N_img
    if image_time(i,2)>=1.5
        if image_time(i,1)-timetemp>=150e-9
            imnum=imnum+1;
            SCMg(imnum,1)=imnum;
            SCMg(imnum,2)=image_time(i,1);
            for istrain=1:N_strain
                if roundn(stress_time(istrain,1),-7)==roundn(image_time(i,1),-7)
                    SCMg(imnum,3)=strain_time(istrain,2);
                    break;
                end
            end
            for istress=1:N_stress
                if roundn(stress_time(istress,1),-7)==roundn(image_time(i,1),-7)
                    SCMg(imnum,4)=stress_time(istress,2);
                    break;
                end
            end
            for isr=1:N_strainr
                if roundn(sr_time(isr,1),-7)==roundn(image_time(i,1),-7)
                    SCMg(imnum,5)=sr_time(isr,2);
                    break;
                end
            end
            timetemp=image_time(i,1);
        end
    end
end
offset=180-SCMg(imnum,1);
SCMg(:,1)=SCMg(:,1)+offset;
%% Plots
figure;
hold on;
% plot(image_time(:,1),image_time(:,2));
plot(SCMg(:,1),SCMg(:,3),'o','MarkerSize',10);
%%
figure;
hold on;
plot(strain_time(:,1),strain_time(:,2)*10^2);
plot(SCMg(:,2),SCMg(:,3)*10^2,'o','MarkerSize',10);
%%
figure;
hold on;
plot(stress_time(:,1),stress_time(:,2)*10^-6);
plot(SCMg(:,2),SCMg(:,4)*10^-6,'o','MarkerSize',10);
%%
figure;
hold on;
plot(sr_time(:,1),sr_time(:,2));
plot(SCMg(:,2),SCMg(:,5),'o','MarkerSize',10);
end
        
    
    
