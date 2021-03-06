%% demo 
clc
clear all
close all
addpath('./minFunc');
N=5000;
Fs=10000;  % sampling frequency
fts=[0:1:N-1]*Fs/N;
t=[0:1:N-1]/Fs;
load('sim_fault.mat');
load('sim_outlier.mat');
load('sim_noise.mat');
sx1=sim_fault+sim_outlier+sim_noise;
%%
close all
aa=300;bb=150;
figure
plot(t,sim_fault,'LineWidth',1)
ylabel('Amplitude','fontsize',12)
xlabel('Time (s)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);
ylim([-2,2])

figure
plot(t,sim_outlier,'LineWidth',1)
ylabel('Amplitude','fontsize',12)
xlabel('Time (s)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);
ylim([-5,7])
yticks([-5,0,7])

figure
plot(t,sim_noise,'LineWidth',1)
ylabel('Amplitude','fontsize',12)
xlabel('Time (s)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);

figure
plot(t,sx1,'LineWidth',1)
ylabel('Amplitude','fontsize',12)
xlabel('Time (s)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);
ylim([-5,7])
yticks([-5,0,7])
%%
figure
plot(fts,abs(fft(sx1))*2/N,'LineWidth',1)
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,300,150]);
% ylim([0,0.3])
xlim([0,5000])
figure
plot(fts,abs(fft(abs(hilbert(sx1))))*2/N,'LineWidth',1)
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
%set(gcf,'position',[200,300,300,150]);
ylim([0,0.1])
xlim([0,800])
set(gcf,'position',[200,300,600,200]);
%%
[~,rec3]=min_lplq(sx1,40,0,1,2);
fts1=(0:length(rec3)-1)*Fs/length(rec3);
NN=length(rec3);
figure
plot(t(1:NN),rec3,'LineWidth',1)
ylabel('Amplitude','fontsize',12)
xlabel('Time (s)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);
ylim([-6,6])
yticks([-6:3:6])
figure
plot(fts1,abs(fft(abs(hilbert(rec3))))*2/length(rec3),'LineWidth',1)
ylabel('Amplitude','fontsize',12)
xlabel('Frequency (Hz)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);
xlim([0,400]);
ylim([0,0.4])