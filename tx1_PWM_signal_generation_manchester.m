% Author: LIU Liqiong
% Email: liuliqiongrabby@gmail.com
% Institution: Department of Information Engineering, 
% The Chinese University of Hong Kong.
% Created: March 2021
% Modified: November 2021
% tx1_PWM_signal_generation_manchester is for generating PAM4 signals input to AWG
% with Manchester-like flicker mitigation scheme . 
% The PAM4 signals is further modulated by PWM scheme. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;
%% parameter setting
symbol_rate =32;% Symbol rate of PAM4 signals (kbaud/s)
frame_rate = 30;

header_length = 20;% Number of symbols for header
pack_length = floor(1/frame_rate/(1/symbol_rate/1000)/3); % Number of symbols in a packet
net_pack_length = pack_length-header_length ;% Number of symbols in a packet excluding header
net_pack_num = 10;% Number of transmitted packets 
rep_num = 4;
pam_n3 = repmat([0 0 0 ],1,rep_num);% PWM signals for negative 3
pam_n1 = repmat([ 0 1 0 ],1,rep_num);% PWM signals for negative 1
pam_p1 = repmat([ 1  0 1 ],1,rep_num);% PWM signals for positive 1
pam_p3 = repmat([ 1  1  1 ],1,rep_num);% PWM signals for positive 3

net_data_rate = net_pack_length/pack_length*symbol_rate/3*0.5*2;
%% Packetization and Manchester coding
load('sequence_pam4.mat');% randomly generated binary sequence
signal = sequence_pam4(1:ceil(net_pack_length/2)*2*net_pack_num);
signal = reshape(signal,2,[]);
signal_pam = bi2de(signal','left-msb');
signal_temp = 3-signal_pam;
signal = [signal_pam';signal_temp'];% Manchester coding
signal = signal(:);
signal = reshape(signal,[],net_pack_num);% Packetization
signal = signal(1:net_pack_length,:);
signal = repmat(signal,3,1);% Packet repetation
signal = signal(:);

%% PWM signal generation
signal_pwm = zeros(length(signal),length(pam_n3));
for i = 1:length(signal)
    switch signal(i)
        case 0
            signal_pwm(i,:) = pam_n3;
        case 1
            signal_pwm(i,:) = pam_n1;
        case 2
            signal_pwm(i,:) = pam_p1;
        case 3
            signal_pwm(i,:) = pam_p3;
    end
end
signal_pwm = signal_pwm';
signal_pwm = signal_pwm(:);

%% Header insertion

header = [ 1 1 1 1, 1 0 1 0, 1 0 1 0 ];% header for synchronization and timing recovery
header = repmat(header,rep_num*3,1);
header = header(:);
header = repmat(header,1,net_pack_num);

No = [0 0 1 1, 0 1 1 0, 1 0 0 1, 1 1 0 0];% header to indicate the order in which the packets appear
No = repmat(No,1,ceil(net_pack_num/4));% These four headers appear in cycles
No = [1 0 1 0 , No];% header 1010 represents the first packet
No = No(1:net_pack_num*4);% PWM modulation
No = repmat(No,rep_num*3,1);
No = reshape(No,4*rep_num*3,[]);

RPT_1 = [1 0 1 0];% header to indicate the order in which the repeatitive packets appear
RPT_1 = repmat(RPT_1,rep_num*3,1);
RPT_1 = RPT_1(:);
RPT_1 = repmat(RPT_1,1,net_pack_num);
RPT_2 = [0 1 0 1];
RPT_2 = repmat(RPT_2,rep_num*3,1);
RPT_2 = RPT_2(:);
RPT_2 = repmat(RPT_2,1,net_pack_num);
RPT_3 = [1 0 0 1];
RPT_3 = repmat(RPT_3,rep_num*3,1);
RPT_3 = RPT_3(:);
RPT_3 = repmat(RPT_3,1,net_pack_num);

signal_pwm = reshape(signal_pwm,[],net_pack_num);
signal_pwm = [header;RPT_1;No;signal_pwm(1:net_pack_length*rep_num*3,:);...
              header;RPT_2;No;signal_pwm(net_pack_length*rep_num*3+1:2*net_pack_length*rep_num*3,:);...
              header;RPT_3;No;signal_pwm(2*net_pack_length*rep_num*3+1:end,:)];
signal_pwm = signal_pwm(:)';
signal_pwm = signal_pwm';

%% insert addtional symbols to fit the data length (2^19-1) in AWG
temp = ceil((2^19-1)/length(signal_pwm));
signal_pwm = repmat(signal_pwm,temp,1);
signal_pwm = signal_pwm(1:(floor((2^19-1)/(pack_length*rep_num*3*3))-1) * (pack_length*rep_num*3*3));

temp = (2^19-1)-length(signal_pwm);
comp = repmat([zeros(12,1);ones(12,1)],ceil(temp/12),1);% insert zeros and ones to the signals
comp = comp(1:temp);
signal_pwm = [signal_pwm;comp];

modulation_fre = symbol_rate*rep_num*3 % modualtion frequency for PWM signals (kHz)
modulation_fre_awg = modulation_fre*1000/(2^19-1) % modulation frequency set in AWG (Hz)


%% save
signal_pwm(signal_pwm==0) = -1;
csvwrite('seq.csv',signal_pwm);