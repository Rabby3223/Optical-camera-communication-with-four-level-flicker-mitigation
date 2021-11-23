% Author: LIU Liqiong
% Email: liuliqiongrabby@gmail.com
% Institution: Department of Information Engineering, 
% The Chinese University of Hong Kong.
% Created: March 2021
% Modified: November 2021
% tx1_PWM_signal_generation_manchester is for generating PAM4 signals input to AWG
% with proposed flicker mitigation scheme . 
% The PAM4 signals is further modulated by PWM scheme. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;
%% paramter setting
symbol_rate =20;% Symbol rate of PAM4 signals (kbaud/s)
frame_rate = 30;

header_length = 20;% Number of symbols for header
pack_length = floor(1/frame_rate/(1/symbol_rate/1000)/3); % Number of symbols in a packet
net_pack_length = pack_length-header_length;% Number of symbols in a packet excluding header
net_pack_num = 10;% Number of transmitted packets 
rep_num = 4;
pam_n3 = repmat([0 0 0 ],1,rep_num);% PWM signals for negative 3
pam_n1 = repmat([ 0 1 0 ],1,rep_num);% PWM signals for negative 1
pam_p1 = repmat([ 1  0 1 ],1,rep_num);% PWM signals for positive 1
pam_p3 = repmat([ 1  1  1 ],1,rep_num);% PWM signals for positive 3
    
block_bit_length = 30;% number of bits in a regular block
r_bit_pre = [1;1];% redundant bits added at the beginning of the block
r_bit_post = [0;0];% redundant bits added at the end of the block
r_bit_length = length(r_bit_pre)+length(r_bit_post);% number of redundant bits in a block
block_num = floor(net_pack_length/(block_bit_length/2));% number of regular blocks in a packet
block_bit_length_temp = (net_pack_length - block_num*block_bit_length/2)*2;% number of bits in the last block

%% Four level flicker mitigation scheme
load('sequence_pam4.mat');% Randomly generated binary sequence
signal = sequence_pam4(1:(block_bit_length-r_bit_length)*block_num*net_pack_num+(block_bit_length_temp-r_bit_length)*net_pack_num);
signal = reshape(signal,[],net_pack_num);% Packetization
signal = repmat(signal,3,1);
signal = reshape(signal,[],net_pack_num*3);% Packet repetation
signal_1 = signal(1:(block_bit_length-r_bit_length)*block_num,:);% Binary bits in regular blocks
signal_2 = signal((block_bit_length-r_bit_length)*block_num+1:end,:);% Binary bits in the last block of each packet

signal_1 = reshape(signal_1,block_bit_length-r_bit_length,[]);
signal_1 = signal_1';
[LENGTH,~] = size(signal_1);
out_1 = zeros(LENGTH,block_bit_length/2);% Regular blocks after flicker mitigation 
signal_2 = reshape(signal_2,block_bit_length_temp-r_bit_length,[]);
signal_2 = signal_2';
[LENGTH,~] = size(signal_2);
out_2 = zeros(LENGTH,block_bit_length_temp/2);% The last block of each packet after flicker mitigation 

RD = -1;% To indicate accumulated AC intensity is smaller than 0 or not.
I = 0;% Accumulated AC intensity.
index_1 = 1;
row_1 = 1;
row_2 = 1;
for i = 1:block_num*net_pack_num*3+net_pack_num*3
    if index_1<block_num+1
        sig = signal_1(row_1,:)';
        sig = [r_bit_pre;sig;r_bit_post];% Insert redundant bits
        sig = reshape(sig,2,[]);
        sig = bi2de(sig','left-msb');% Convert binary bits to PAM4 symbols
        sig1 = sig';% Original PAM4 symbols
        sig2 = 3-sig1;% Inverted PAM4 symbols
        if sum(sig1)>sum(sig2)
            sig = sig2;
            sig2 = sig1;
            sig1 = sig;
        end
        % Select the block based on the accumulated AC intensity
        if(RD==-1)
            out_1(row_1,1:block_bit_length/2) = sig1;
        else
            out_1(row_1,1:block_bit_length/2) = sig2;
        end
        % Calcualte the accumulated AC intensity
        I = I+sum(out_1(row_1,1:block_bit_length/2))-block_bit_length/2/2*3;
        if(I<0)
            RD = 1;
        else
            RD = -1;
        end
        index_1 = index_1+1;
        row_1 = row_1+1;
    else
        index_1 = 1;
        sig = signal_2(row_2,:)';
        sig = [r_bit_pre;sig;r_bit_post];% Insert redundant bits
        sig = reshape(sig,2,[]);
        sig = bi2de(sig','left-msb');% Convert binary bits to PAM4 symbols
        sig1 = sig';% Original PAM4 symbols
        sig2 = 3-sig1;% Inverted PAM4 symbols
        if sum(sig1)>sum(sig2)
            sig = sig2;
            sig2 = sig1;
            sig1 = sig;
        end
        % Select the block based on the accumulated AC intensity
        if(RD==-1)
            out_2(row_2,1:block_bit_length_temp/2) = sig1;
        else
            out_2(row_2,1:block_bit_length_temp/2) = sig2;
        end
        % Calcualte the accumulated AC intensity
        I = I+sum(out_2(row_2,1:block_bit_length_temp/2))-block_bit_length_temp/2/2*3;
        if(I<0)
            RD = 1;
        else
            RD = -1;
        end
        row_2 = row_2+1;
    end
end
% Combine the blocks
out_1 = out_1';
out_1 = reshape(out_1,[],net_pack_num*3);
out_2 = out_2';
out = [out_1;out_2];
signal = out;
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
I = repmat([zeros(12,1);ones(12,1)],ceil(temp/12),1);
I = I(1:temp);
signal_pwm = [signal_pwm;I];

modulation_fre = symbol_rate*rep_num*3 % modualtion frequency for PWM signals (kHz)
modulation_fre_awg = modulation_fre*1000/(2^19-1) % modulation frequency set in AWG (Hz)


%% save
signal_pwm(signal_pwm==0) = -1;
csvwrite('seq.csv',signal_pwm);