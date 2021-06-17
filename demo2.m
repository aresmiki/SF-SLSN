%% demo 
%% Code By He Liu (aresmiki@163.com)
clc
clear all
close all
addpath('./minFunc');
N=5000;
Fs=10000;  % sampling frequency
fts=[0:1:N-1]*Fs/N;
t=[0:1:N-1]/Fs;
load('sim_fault.mat');
load('sim_noise.mat');
sx1=sim_fault+0.4*sim_noise;
%%
close all
[optWQ,recQ,obj,optWQ2,recQ2,obj2]=min_lplq_optimation(sx1,100,5,1,2);
fts1=(0:length(recQ)-1)*Fs/length(recQ);
NN=length(recQ);
figure
plot(t(1:NN),recQ,'LineWidth',1)
ylabel('Amplitude','fontsize',12)
xlabel('Time (s)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);
ylim([-6,6])
yticks([-6:3:6])

figure
plot(fts1,abs(fft(((recQ))))*2/length(recQ),'LineWidth',1)
ylabel('Amplitude','fontsize',12)
xlabel('Frequency (Hz)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);


figure
plot(fts1,abs(fft(abs(hilbert(recQ))))*2/length(recQ),'LineWidth',1)
ylabel('Amplitude','fontsize',12)
xlabel('Frequency (Hz)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);
xlim([0,400]);
ylim([0,0.4])

fts1=(0:length(recQ2)-1)*Fs/length(recQ2);
NN=length(recQ2);
figure
plot(t(1:NN),recQ2,'LineWidth',1)
ylabel('Amplitude','fontsize',12)
xlabel('Time (s)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);
ylim([-6,6])
yticks([-6:3:6])

figure
plot(fts1,abs(fft(((recQ2))))*2/length(recQ2),'LineWidth',1)
ylabel('Amplitude','fontsize',12)
xlabel('Frequency (Hz)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);

figure
plot(fts1,abs(fft(abs(hilbert(recQ2))))*2/length(recQ2),'LineWidth',1)
ylabel('Amplitude','fontsize',12)
xlabel('Frequency (Hz)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);
xlim([0,400]);
ylim([0,0.4])

figure
plot(obj)
hold on
plot(obj2)
