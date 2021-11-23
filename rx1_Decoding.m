% Author: LIU Liqiong
% Email: liuliqiongrabby@gmail.com
% Institution: Department of Information Engineering, 
% The Chinese University of Hong Kong.
% Created: March 2021
% Modified: November 2021
% rx1_Decoding is for decoding PAM4 signal processed by flicker mitigation scheme. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
%% Parameter and data loading

% Set parameters based on the folder name (pam4_symbolRate_redundantBitLength_blockLength).

% The data length of AWG is fixed, so some additional data should be
% inserted to fullfill the requirement. To avoid large number of additional
% data, the symbol rate may not be an integral number.

disp('Parameter setting: start.');
symbol_rate =23;% Symbol rate of PAM4 signals (kbaud/s)
flicker_mitigation = 1;% 1 for Manchester-like scheme, and 2 for proposed scheme.
block_bit_length = 20;
r_bit = [1;1;0;0];

disp('Parameter setting: generate transmitted PAM4 symbols');
if flicker_mitigation == 2
    [sequence_R,pack_length] = PAM4_generation(symbol_rate,block_bit_length,r_bit);
    filename = strcat('pam4_',num2str(symbol_rate),'k_rbit',num2str(length(r_bit)),'_',...
        num2str(block_bit_length));
end
if flicker_mitigation == 1
    [sequence_R,pack_length] = PAM4_generation(symbol_rate);
    filename = strcat('pam4_',num2str(symbol_rate),'k_manchester');
end
disp('Parameter setting: calcualte scaling factors');
img_num =200;% number of images used for averaging
[scaling_factor] = scaling_function(filename,img_num);
% load('sequence_save_v2.mat');% load the transmitted PAM4 symbols, symbol rate, and packet length.
% load('scaling_factor');% load the scaling factors.
sampling_rate = 70.4260;% unit: kSa/s, the rolling shutter sampling rate for Oneplus 5T with resolution 960x1280.
header_index = 0;% header index (0,1,2,3,4) indicates the No in the header of current packet
% LMS equalizer initial coefficients 
beta =[-0.0161506051240003,-0.0100264324566792,0.0724609484713563,-0.141985195139301,-0.0289347870496333,0.314372438070751,1.35875809047925,1.28199915685150,-0.526149373747433,-0.0899502573590000];
seq = [sequence_R(end-1000:end)';sequence_R'];
seq = seq';

BER = [];LOSS = [];ERROR = [];
header = [1 1 1 1,1 0 1 0,1 0 1 0];
signal_frame_head = [];
signal_frame = [];
disp('Parameter setting: finished.');

disp('Decding: start.');
image_index_total = 299;
for image_index =0:image_index_total % process images
    
    filename_img = strcat(filename,'/',num2str(image_index),'.csv');
    [sig] = textread(filename_img,'%s');% Read the csv file
    [row_num,~] = size(sig);
    if row_num~=0
        
        % decode signals in base64 format
        
        % The function base64decode is from poltly graphing library for
        % MATLAB, which can be found from the following link.
        % https://github.com/plotly/plotly_matlab/blob/master/plotly/plotly_aux/base64decode.m
        
        sig_temp = char(sig(1));
        sig_temp_start = base64decode(sig_temp);
        [~,col_num] = size(sig_temp_start);
        sig_rx = uint8(zeros(row_num-1,col_num));
        for i = 1:row_num-1
            sig_temp = char(sig(i));
            sig_rx(i,:) = base64decode(sig_temp);
        end
        sig_temp = char(sig(row_num));
        sig_temp_end = base64decode(sig_temp);

        sig_rx = sig_rx';
        sig_rx = sig_rx(:);
        sig_rx = [sig_rx;sig_temp_end'];
        sig_rx = double(sig_rx);
        
        sig_rx = sig_rx./scaling_factor;% Scaling
        sig_scaled = (sig_rx-mean(sig_rx))*2;% Normalization and clipping
        sig_scaled(sig_scaled>1.5) = 1.5;
        sig_scaled(sig_scaled<-1.5) = -1.5;

        % Synchronization by correlating header with the received signal
        
        % The function corrSyn was initially developed by Dr. Shuang Gao,
        % and further optimized by memebers from Lightwave Lab.
        [DeWaveform,P]=corrSyn(sig_scaled,header-0.5,symbol_rate,sampling_rate);
        P = P+1;
        % Process signals from the next image if the synchronization is not successful. 
        if length(find(DeWaveform>1))>length(DeWaveform)*0.4
            continue;
        end
        
        % oversampling the synchronized signals by two times.
        ratio  = 0.5;
        q = 1:length(DeWaveform);
        xi = 1:ratio:length(DeWaveform);
        sig_resamp=interp1(q,DeWaveform,xi,'linear') ;
        
        % LMS equalization
        upsample = 2;
        learning_rate = 1e-3;% learning rate in LMS equalizer
        sig_decoded = zeros(1,floor(length(sig_resamp)/2)-4);
        for i = 3:floor(length(sig_resamp)/2)-2
            sig_resamp_temp = [sig_resamp((i-3)*upsample+1:(i+2)*upsample)];
            Test = sum(beta.*sig_resamp_temp);
            if(Test>2)
                sig_decoded(i-2) = 3;
            else
                if(Test>0 && Test<2)
                sig_decoded(i-2) = 1;
                
                else
                    if(Test>-2 && Test<0)
                        sig_decoded(i-2) = -1;
                    else
                        sig_decoded(i-2) = -3;
                    end
                end
            end 
            beta = beta + learning_rate*((sig_decoded(i-2))-Test)*sig_resamp_temp;% Update coefficients
        end
        sig_decoded = (sig_decoded+3)./2;

        P = P-3;
        if P<1
            P = P+pack_length;
        end
        % Packet reconstruction
        [signal_frame_head,header_index,signal_frame] = packet_reconstruction_PAM4(sig_decoded,P,...
        signal_frame_head,pack_length,...
            header_index,image_index,flicker_mitigation,r_bit,block_bit_length,signal_frame);
    end
    if (image_index+1)/100==fix((image_index+1)/100)
        disp(strcat('Decoding:',num2str((image_index+1)),'/',num2str(image_index_total+1)));
    end
    close all;
end
disp('Decding: finished.');
%% BER calculation
disp('BER calculation: start.');
seq = sequence_R;seq = [seq,seq];corr = zeros(1,length(seq));
if flicker_mitigation == 2
    block_num = floor((pack_length-20)/(block_bit_length/2));% number of regular blocks in a packet
    delete_index = zeros( 1,(block_num+1)*length(r_bit) );% index for bits to delete, such as redundant bits
    % Calculate the index for redundant bits
    kkk = 1;
    for i = 1:length(r_bit)/2
        delete_index(kkk:kkk+block_num) = i:block_bit_length:(pack_length-20)*2;
        kkk = kkk+block_num+1;
    end
    for i = 1:length(r_bit)/2
        delete_index(kkk:kkk+block_num-1) = block_bit_length-length(r_bit)/2+i:block_bit_length:(pack_length-20)*2;
        delete_index(kkk+block_num) = (pack_length-20)*2-length(r_bit)/2+i;
        kkk = kkk+block_num+1;
    end

    [row,col] = size(signal_frame);
    signal_frame(:,1:20) = [];% Remove header
    for i = 1:row
        temp = signal_frame(i,:);
        temp(temp==-1) = 1.5;
        temp = temp-1.5;
        % find the transmitted sequence by correlation method
        for j = 21:pack_length:length(seq)-length(temp)
                sequence_temp = seq(j:j+length(temp)-1)-1.5;
                corr(j) = sum(temp.*sequence_temp);
        end
        hh = find(corr==max(corr));
        hh = hh(1);
        ss = seq(hh:hh+length(temp)-1);% corresponding transmitted PAM4 symbols
        % Convert the PAM4 symbols to binary data and calcualte BER
        index_temp_1 = find(temp==0);% Find block where the redundant bits have errors                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  p==0);
        if ~isempty(index_temp_1)
            index_temp_1 = [(index_temp_1-1)*2+1;index_temp_1*2];% Calculate the bit indexes
            index_temp = union(index_temp_1,delete_index);% Combine the index_temp_1 and delete_index with no repetitions
        else
            index_temp = delete_index;
        end
        % Convert the PAM4 symbols to binary data and calcualte BER
        temp = temp+1.5;temp(temp==1.5) = 0;
        temp_bit = de2bi(temp,'left-msb',2);
        temp_bit = temp_bit';temp_bit = temp_bit(:);
        ss_bit = de2bi(ss,'left-msb',2);
        ss_bit = ss_bit';ss_bit = ss_bit(:);
        % Delete the redundant bits and blocks where the redundant bits have errors
        temp_bit(index_temp) = [];
        ss_bit(index_temp) = [];
        if ~ isempty(ss_bit)
            BER = [BER;sum(abs(temp_bit-ss_bit))/length(ss_bit)];
            LOSS = [LOSS;length(index_temp_1)/(block_bit_length)];
        end
    end
    % Burst error occurs when image sensor captures complementary data 
    % inserted to the transmitted seuqence to fulfill the data length required
    % in AWG 2^19-1.Burst errors are not included when calculating the total BER.

    % This burst error will not occur if data with arbitray length can be
    % uploaded to the waveform generator
    ii = find(LOSS>15);
    LOSS(ii) = [];
%     LOSS = floor(LOSS);
    BER(ii) = [];
    disp(strcat('BER:',num2str(mean(BER))));
    disp(strcat('Block loss rate:',num2str(sum(LOSS)/row/(pack_length-20))));
end
if flicker_mitigation == 1
    [row,col] = size(signal_frame);
    signal_frame(:,1:20) = [];% Remove header
    % calculate BER packet by packet
    for i = 1:row
        temp = signal_frame(i,:);
        temp = temp-1.5;
        % find the transmitted sequence by correlation method
        for j = 21:pack_length:length(seq)-length(temp)
                sequence_temp = seq(j:j+length(temp)-1)-1.5;
                corr(j) = sum(temp.*sequence_temp);
        end
        hh = find(corr==max(corr));
        hh = hh(1);
        ss = seq(hh:hh+length(temp)-1);% corresponding transmitted PAM4 symbols
        % Convert the PAM4 symbols to binary data and calcualte BER
        temp = temp+1.5;
        temp_bit = de2bi(temp,'left-msb',2);
        temp_bit = temp_bit';temp_bit = temp_bit(:);
        ss_bit = de2bi(ss,'left-msb',2);
        ss_bit = ss_bit';ss_bit = ss_bit(:);
        temp_bit(2:2:end) = [];
        ss_bit(2:2:end) = [];
        BER = [BER;sum(abs(temp_bit-ss_bit))/length(ss_bit)];
    end
    % Burst error occurs when image sensor captures complementary data 
    % inserted to the transmitted seuqence to fulfill the data length required
    % in AWG 2^19-1.Burst errors are not included when calculating the total BER.

    % This burst error will not occur if data with arbitray length can be
    % uploaded to the waveform generator
    BER(BER>0.1) = [];
    disp(strcat('BER:',num2str(mean(BER))));
end
disp('BER calculation: finished.');
