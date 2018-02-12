%FEB 20, 2015_CKB Data Processing- Vignesh Kannan, Ramesh Lab (Latrobe 026)/
% 9/16" MARAGING STEEL BARS
%-----------INPUT PARAMETERS---------------------------------------/
%datafile- Perception data file with all of the signals in one file
%GFi - strain gauge factor for incident gauge
%GFt - strain gauge factor for transmitted gauge
%Vii- excitation voltage for the incident gauge
%Vit- excitation voltage for the transmitted gauge
%Ab - area of cross-section of the bar
%Aso - area of cross-section of the specimen
%Eb - young's modulus of the bar
%cb - bar wave speed
%Ls - Length of the specimen
%ti- initial time
%tf - final time
function []=KolskyBar_dataprocessing()
clc;
close all;
%% INPUT PARAMETERS
Ab=(pi*12.7e-3*12.7e-3)/4;
Aso=(3.575e-3)*(4.025e-3);
Eb=200e9;
cb=5000;
Ls=2.997e-3;
%% Read raw data file & Calculate Bar Strains from function
[barstrain_incident,barstrain_transmitted,ti,tf,camout]=BridgeCktAnalysis();   %compute barstrains
%% Apply a manual Filter to signals (Filter written by Vignesh Kannan, Ramesh Lab, Latrobe 026)
incident=Filter_Kannan_realsignal(barstrain_incident(:,1),barstrain_incident(:,2));
reflected=Filter_Kannan_realsignal(barstrain_incident(:,1),barstrain_incident(:,2));
transmitted=Filter_Kannan_realsignal(barstrain_incident(:,1),barstrain_transmitted(:,2)); 
cam=camout(:,2);
%% NO FILTER
% incident=barstrain_incident(:,2);
% reflected=barstrain_incident(:,2);
% transmitted=barstrain_transmitted(:,2); 
%% Matchup the 3 signals to the specimen
%Signal matchup
timeincident=barstrain_incident(:,1)+(1045/5)*(10^-6)+4.8e-6;
timereflected=barstrain_incident(:,1)-(1045/5)*(10^-6)-4.8e-6;
timetransmitted=barstrain_transmitted(:,1)-(535/5)*(10^-6)-8.5e-6+1.2e-6;
timecam=barstrain_incident(:,1);
%% Plot matched-up signals
figure;
plot(timeincident*(10^6),incident*(10^6));
hold on;
plot(timereflected*(10^6),-reflected*(10^6),'color','red');
plot(timetransmitted*(10^6),transmitted*(10^6),'color','green');
plot(timecam*(10^6),cam*(100),'color',[0 0.5 0]);
title('Bar Strain Data- Matched up to specimen','FontSize',18,'FontName','Arial Narrow')
xlabel('Time(\mus)','FontSize',14,'FontName','Arial');
ylabel('Strain x 10^{-6} ','FontSize',14,'FontName','Arial');
legend('Incident pulse','Reflected pulse','Transmitted pulse');
grid on;
hold off;
pause;
%% Identify indices in the three traces for a common start point
buffer1=0;
for i=1:size(timeincident)
    for l=1:size(timecam)
        if round(timeincident(i),7)==round(timecam(l),7)
            for j=1:size(timereflected)
                if (round(timecam(l),7)==round(timereflected(j),7))
                    for k=1:size(timetransmitted)
                        if(round(timetransmitted(k),7)==round(timereflected(j),7))
                            indexincident=i;
                            indexreflected=j;
                            indextransmitted=k;
                            indexcam=l;
                            buffer1=1; %to determine the common start time of the three shifted pulses
                            break;
                        end
                    end
                else
                    continue;
                end
                if buffer1==1
                    break; %means we have found a common starting point for the three shifted signals
                end
            end
            if buffer1==1
                break;
            end
        end
    end
    if buffer1==1
        break;
    end
end
%% Identify image number
cam

figure;
plot(timecam,cam,'b');
hold on;
pause;

N_image=size(cam);
N_im=N_image(1);

N_im

imnum=zeros(180,3);
indcam=1;
timetemp=timecam(1);
for i=1:N_im
    sprintf('InLoop');
    if cam(i)>=2.125 && cam(i)<=2.25
        sprintf('Run');
        if timecam(i)-timetemp>=200e-9
            imnum(indcam,1)=timecam(i);
            imnum(indcam,2)=cam(i);
            imnum(indcam,3)=indcam;
            timetemp=imnum(indcam,1);
            indcam=indcam+1;
        end
    end
end
indcam
imnum
plot(imnum(:,1),imnum(:,2),'-ro');
hold off;
pause;
%%
% for l=1:size(camera)
%     if(roundn(barstrain_incident(l,1),-7)==roundn(timetransmitted(indextransmitted),-7))   
%         indexcam=l;
%         break;
%     end
% end 
%% Identify start and end points of interest in the signal
for i=1:size(timeincident)
    if(barstrain_incident(i,1)>=ti)
        breakpoint1=i;
        break;
    end
end

for i=1:size(timeincident)
    if(barstrain_incident(i,1)>tf)
        breakpoint2=i;
        break;
    end
end
%% Truncate data
sizedata=breakpoint2-indexincident; %indexincident is used because the breakpoint is measured with respect to the incident data
for k=1:sizedata
    datamaterial{1,1}(k)=timeincident(indexincident);
    datamaterial{1,2}(k)=incident(indexincident);
    datamaterial{1,3}(k)=reflected(indexreflected);
    datamaterial{1,4}(k)=transmitted(indextransmitted);
    datamaterial{1,5}(k)=cam(indexcam);
    indexincident=indexincident+1;
    indexreflected=indexreflected+1;
    indextransmitted=indextransmitted+1;
    indexcam=indexcam+1;
end
%% Truncate data further between startpoint and end point
%record data for just the region of interest between breakpoint1 and
%breakpoint2-datamat hold the data we require for the first loading pulse
for l=1:sizedata
    if (datamaterial{1,1}(l)==timeincident(breakpoint1))
        startpoint=l;
        break;
    end
end
datamaterial{1,1}(startpoint)
m=1;
offset1=0;
offset2=0;
offset3=0;
%offset1=datamaterial{1,2}(startpoint);
%offset2=datamaterial{1,3}(startpoint);
%offset3=datamaterial{1,4}(startpoint);
for k=startpoint:sizedata
    datamat{1,1}(m)=datamaterial{1,1}(k); %time
    datamat{1,2}(m)=datamaterial{1,2}(k)-offset1; %incident 
    datamat{1,3}(m)=datamaterial{1,3}(k)-offset2; %reflected
    datamat{1,4}(m)=datamaterial{1,4}(k)-offset3; %transmitted
    datamat{1,5}(m)=datamaterial{1,5}(k); %camera signal
    m=m+1;
end 
%% Plot Truncated Data before dispersion correction
figure;
plot(datamat{1,1}*(10^6),datamat{1,2}*(10^6));
hold on;
plot(datamat{1,1}*(10^6),-datamat{1,3}*(10^6),'r');
plot(datamat{1,1}*(10^6),datamat{1,4}*(10^6),'g');
plot(datamat{1,1}*(10^6),datamat{1,5},'g');
plot(imnum(:,1)*10^6,imnum(:,3)*10,'ob','MarkerSize',10);
title('Truncated Bar Strain Data- Matched up to specimen','FontSize',18,'FontName','Arial Narrow')
xlabel('Time(\mus)','FontSize',14,'FontName','Arial');
ylabel('Strain x 10^{-6} ','FontSize',14,'FontName','Arial');
legend('Incident Gauge','Reflected Gauge','Transmitted Gauge','Camera Signal');
grid on;
hold off;
%% Dispersion Correction
datamat{1,1}=transpose(datamat{1,1});
datamat{1,2}=transpose(datamat{1,2});
datamat{1,3}=transpose(datamat{1,3});
datamat{1,4}=transpose(datamat{1,4});
datamat{1,5}=transpose(datamat{1,5});
[DispCorr_datamat2,~]=Disp_corr_latest_Bancroft(datamat{1,2},datamat{1,1},622.5e-3);
[DispCorr_datamat3,~]=Disp_corr_latest_Bancroft(datamat{1,3},datamat{1,1},-622.5e-3);
% [DispCorr_datamat4,~]=Disp_corr_latest_Bancroft(datamat{1,4},datamat{1,1},-727e-3);
% datamat{1,1}=T1;
datamat{1,2}=DispCorr_datamat2;
datamat{1,3}=DispCorr_datamat3;
% datamat{1,4}=DispCorr_datamat4;
%% Plot Truncated Data after dispersion correction
figure;
plot(datamat{1,1}*(10^6),datamat{1,2}*(10^6));
hold on;
plot(datamat{1,1}*(10^6),-datamat{1,3}*(10^6),'r');
plot(datamat{1,1}*(10^6),datamat{1,4}*(10^6),'g');
plot(datamat{1,1}*(10^6),datamat{1,5},'g');
title('Truncated Bar Strain Data- Matched up to specimen','FontSize',18,'FontName','Arial Narrow')
xlabel('Time(\mus)','FontSize',14,'FontName','Arial');
ylabel('Strain x 10^{-6} ','FontSize',14,'FontName','Arial');
legend('Incident Gauge','Reflected Gauge','Transmitted Gauge','Camera Signal');
grid on;
hold off;
%% Concatenate Data
ti_final=input('Enter start time for final processed data');
N_data=size(datamat{1,1});
for i=1:N_data(1)
    if roundn(datamat{1,1}(i),-7)==roundn(ti_final,-7)
        n=i;
        break;
    end
end
size(datamat{1,1})
size(datamat{1,5})
for i=n:N_data(1)
    data{1,1}(i-n+1,1)=datamat{1,1}(i,1);
    data{1,2}(i-n+1,1)=datamat{1,2}(i,1);
    data{1,3}(i-n+1,1)=datamat{1,3}(i,1);
    data{1,4}(i-n+1,1)=datamat{1,4}(i,1);
    data{1,5}(i-n+1,1)=datamat{1,5}(i,1);
end
%% Plot Truncated Data
figure;
plot(data{1,1}*(10^6),data{1,2}*(10^6));
hold on;
plot(data{1,1}*(10^6),-data{1,3}*(10^6),'r');
plot(data{1,1}*(10^6),data{1,4}*(10^6),'g');
plot(datamat{1,1}*(10^6),datamat{1,5}*100,'g');
title('Truncated Bar Strain Data- Matched up to specimen','FontSize',18,'FontName','Arial Narrow')
xlabel('Time(\mus)','FontSize',14,'FontName','Arial');
ylabel('Strain x 10^{-6} ','FontSize',14,'FontName','Arial');
legend('Incident Gauge','Reflected Gauge','Transmitted Gauge','Camera Signal');
grid on;
hold off;  
%% Specimen stress, strainrate, strain computation
stress=Eb*data{1,4}*(Ab/Aso);
strainrate= -(2*cb/Ls)*data{1,3};  %strain rate from reflected signal data
strain=cumtrapz(data{1,1},strainrate); %integrate strainrate with respect to time to get strain
N_strain=size(strain);
true_strain=zeros(N_strain(1),1);
true_stress=zeros(N_strain(1),1);
true_strainrate=zeros(N_strain(1),1);
for i=1:N_strain(1)
true_strain(i)=-log(1-strain(i));
true_stress(i)=stress(i)*(1-strain(i));
true_strainrate(i)=strainrate(i)/(1-strain(i));
end
%% SAVE ALL DATA INTO ONE VARIABLE
DATAFULL{1,1}=data{1,1};
DATAFULL{1,2}=true_stress;
DATAFULL{1,3}=true_strain;
DATAFULL{1,4}=true_strainrate;
%% Plot FINAL
figure;
Kolsky(1)=subplot(2,1,1); 
hold on;
plot(data{1,1}*(10^6),data{1,2}*(10^6),'linewidth',0.5,'color','blue');
plot(data{1,1}*(10^6),data{1,3}*(10^6),'linewidth',0.5,'color','red');
plot(data{1,1}*(10^6),data{1,4}*(10^6),'linewidth',0.5,'color','green');
title('Raw strain signals matched up in time to the specimen','FontSize',18,'FontName','Arial Narrow')
xlabel('Time(\mus)','FontSize',14,'FontName','Arial');
ylabel('Strain x 10^{-6} ','FontSize',14,'FontName','Arial');
legend('Incident Gauge','Reflected Gauge','Transmitted Gauge');
grid on;

Kolsky(2)=subplot(2,1,2); 
hold on;
[X_1,Y1_1,Y2_1]=plotyy(data{1,1}*(10^6),true_stress*(10^-6),data{1,1}*(10^6),true_strain);
plot(data{1,1}*(10^6),data{1,5}*10);
plot(imnum(:,1)*(10^6),imnum(:,3));
title('Stress-Time and Strain-Time plot','FontSize',18,'FontName','Arial Narrow')
xlabel('Time(\mus)','FontSize',14,'FontName','Arial');
ylabel(X_1(1),'Stress(MPa) ','FontSize',14,'FontName','Arial');
ylabel(X_1(2),'Strain','FontSize',14,'FontName','Arial');
set(Y1_1,'linewidth',0.5,'color','blue');
set(Y2_1,'linewidth',0.5,'color','red');
grid on;
hold off;

figure;  
hold on;
[X_2,Y1_2,Y2_2]=plotyy(true_strain*(100),true_stress*(10^-6),true_strain*100,true_strainrate);
plot(true_strain*100,data{1,5}*10);
% T1={datafile;'Stress-Strain Plot/Strain rate-Strain plot'};
title('Stress-Strain Plot/Strain rate-Strain plot','FontSize',18,'FontName','Arial Narrow')
xlabel('Strain(%)','FontSize',14,'FontName','Arial');
ylabel(X_2(1),'True Stress(MPa) ','FontSize',14,'FontName','Arial');
ylabel(X_2(2),'Strainrate (s^{-1})','FontSize',14,'FontName','Arial');
set(Y1_2,'linewidth',0.5,'color','blue');
set(Y2_2,'linewidth',0.5,'color','green');
hold off;

figure;
title('Force Balance across the specimen');
hold on;
Force_INC= (data{1,2}+data{1,3})*Eb*Aso;
Force_TR=data{1,4}*Eb*Aso;
plot(data{1,1}*(10^6),Force_INC,'linewidth',0.5,'color','blue');
plot(data{1,1}*(10^6),Force_TR,'linewidth',0.5,'color','red');
plot(data{1,1}*(10^6),Force_TR-Force_INC,'linewidth',1.5,'color','green');
title('Force Balance across specimen','FontSize',18,'FontName','Arial Narrow')
xlabel('Time(\mus)','FontSize',14,'FontName','Arial');
ylabel('Force (N) ','FontSize',14,'FontName','Arial');
legend('Force at Incident-specimen interface','Force at specimen-transmitted interface','Force Balance');
hold off;

linkaxes(Kolsky,'x');
%% Plot bulk stress-strain data with images
N_bulk=size(true_strain);
N_b=N_bulk(1);
buf=0;
figure;
plot(imnum(:,1));
hold on;
plot(DATAFULL{1,1});
hold off;
pause;
for j=1:180
    for i=1:N_b
        if roundn(imnum(j,1),-7)==roundn(DATAFULL{1,1}(i),-7)
            startind=j;
            imnum(startind,4)=DATAFULL{1,3}(i);
            imnum(startind,5)=DATAFULL{1,2}(i);
            buf=1;
            break;
        end
        if buf==1
            break;
        end
    end
end
for i=startind+1:180
   for j=1:N_b
      if roundn(imnum(i,1),-7)==roundn(DATAFULL{1,1}(j),-7)
          imnum(i,4)=DATAFULL{1,3}(j);
          imnum(i,5)=DATAFULL{1,2}(j);
          break;
      end
   end
end
figure;
plot(DATAFULL{1,3},DATAFULL{1,2});
hold on;
plot(imnum(:,4),imnum(:,5),'o');
%% Images
IMG(180).image=zeros(768,924);
% Get folder path
folder=uigetdir('D:\Kannan_Research\Lab\KANNAN_PROCESSING DATA\2015_10(Oct)_Single Crystal Mg\2015_OCT_SINGLE CRYSTAL\SCMg_a_011\SCMg011\jpeg');
for i=1:180
    filename=strcat(sprintf('SCMg011 %03d',i-1),'.jpg');
    fullpath=fullfile(folder,filename);
    img=imread(fullpath);
%     IMG(i).image=imcrop(img,[300,142,390,464]);
    IMG(i).image=img;
    imshow(IMG(i).image);
end
%% VideoWriter
mkdir(folder,'procIm');
cfolder=fullfile(folder,'procIm');
movfilename=fullfile(cfolder,'processedtimemov.avi');
outputVid=VideoWriter(movfilename,'Uncompressed AVI');
outputVid.FrameRate=2;
open(outputVid);
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:170
    subplot(1,2,1);
    plot((DATAFULL{1,3}(1:107150))*10^2,DATAFULL{1,2}(1:107150)*10^-6,'-o','LineWidth',6,'Color','black');
    hold on;
    plot(imnum(i,4)*10^2,imnum(i,5)*10^-6,'o','MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','g');
%     plot((DATAFULL{1,1}(1:6000)-DATAFULL{1,1}(1))*10^6,DATAFULL{1,2}(1:6000)*10^-6,'-o','LineWidth',6,'Color','black');
%     hold on;
%     plot((imnum(i,1)-imnum(1,1))*10^6,imnum(i,5)*10^-6,'o','MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','g');
    xlabel('True Strain (%)','FontWeight','bold','FontSize',30,'FontName','Times New Roman');
%     xlabel('Time (\mu s)','FontWeight','bold','FontSize',30,'FontName','Times New Roman');
    ylabel('True Stress (MPa)','FontWeight','bold','FontSize',30,'FontName','Times New Roman');
%     set(gca,'LineWidth',2,'Position',[0.08,0.11,0.441,0.815],'FontSize',30,'FontName','Arial Narrow','FontWeight','bold','XLim',[0 60],'YLim',[0 600],'XMinorTick','on','YMinorTick','on');
    set(gca,'LineWidth',2,'Position',[0.08,0.11,0.441,0.815],'FontSize',30,'FontName','Arial Narrow','FontWeight','bold','XLim',[0 25],'YLim',[0 600],'XMinorTick','on','YMinorTick','on');
    hold off;
    subplot(1,2,2);
    imshow(IMG(imnum(i,3)).image);
    set(gcf,'PaperPositionMode','auto');
    %% save as tiff
    file=fullfile(cfolder,sprintf('img%03',i));
    matname=strcat(file,'.fig');
    savefig(matname);
    print(file,'-dpng','-r100');
    readfile=strcat(file,'.png');
    %% Write to vid
    fig=imread(readfile);
    writeVideo(outputVid,fig);
end
close(outputVid);
end
