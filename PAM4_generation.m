function [sequence_R,pack_length] = PAM4_generation(symbol_rate,block_bit_length,r_bit)
% Author: LIU Liqiong
% Email: liuliqiongrabby@gmail.com
% Institution: Department of Information Engineering, 
% The Chinese University of Hong Kong.
% Created: March 2021
% Modified: November 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin == 3
        %% parameter setting
        frame_rate = 30;
        header_length = 20;% Number of symbols for header
        pack_length = floor(1/frame_rate/(1/symbol_rate/1000)/3);% Number of symbols in a packet
        net_pack_length = pack_length-header_length ;% Number of symbols in a packet excluding header
        net_pack_num =10;% Number of transmitted packets 
        if symbol_rate==28 && block_bit_length==30
            net_pack_num = 8;% Limited by the input size of AWG
        end

        r_bit_length = length(r_bit);% number of redundant bits in a block
        r_bit_pre = r_bit(1:r_bit_length/2);% redundant bits added at the beginning of the block
        r_bit_post = r_bit(r_bit_length/2+1:end);% redundant bits added at the end of the block

        block_num = floor(net_pack_length/(block_bit_length/2));% number of regular blocks in a packet
        block_bit_length_temp = (net_pack_length - block_num*block_bit_length/2)*2;% number of bits in the last block

        net_data_rate = symbol_rate*2/3*(net_pack_length-(block_num+1)*r_bit_length/2)/pack_length;
        disp(strcat('net data rate:',num2str(net_data_rate),'kbit/s'));
        %% Packetization and flicker mitigation
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
                sig = [r_bit_pre;sig;r_bit_post]; 
                sig = reshape(sig,2,[]);
                sig = bi2de(sig','left-msb');
                sig1 = sig';
                sig2 = 3-sig1;
                if(RD==-1)
                    out_1(row_1,1:block_bit_length/2) = sig1;
                else
                    out_1(row_1,1:block_bit_length/2) = sig2;
                end
                I = I+sum(out_1(row_1,1:block_bit_length/2))-block_bit_length/2/2*3;
                if(I<0)
                    RD = -1;
                else
                    RD = -1;
                end
                index_1 = index_1+1;
                row_1 = row_1+1;
            else
                index_1 = 1;
                sig = signal_2(row_2,:)';
                sig = [r_bit_pre;sig;r_bit_post];
                sig = reshape(sig,2,[]);
                sig = bi2de(sig','left-msb');
                sig1 = sig';
                sig2 = 3-sig1;
                if(RD==-1)
                    out_2(row_2,1:block_bit_length_temp/2) = sig1;
                else
                    out_2(row_2,1:block_bit_length_temp/2) = sig2;
                end
                I = I+sum(out_2(row_2,1:block_bit_length_temp/2))-block_bit_length_temp/2/2*3;
                if(I<0)
                    RD = -1;
                else
                    RD = -1;
                end
                row_2 = row_2+1;
            end
        end
        out_1 = out_1';
        out_1 = reshape(out_1,[],net_pack_num*3);
        out_2 = out_2';
        out = [out_1;out_2];
    end
    if nargin == 1
       %% parameter setting
        frame_rate = 30;

        header_length = 20;% Number of symbols for header
        pack_length = floor(1/frame_rate/(1/symbol_rate/1000)/3);% Number of symbols in a packet
        net_pack_length = pack_length-header_length ;% Number of symbols in a packet excluding header
        net_pack_num =10;% Number of transmitted packets 
        net_data_rate = symbol_rate*2/3*(net_pack_length/2)/pack_length;
        disp(strcat('net data rate:',num2str(net_data_rate),'kbit/s'));
       %% Packetization and Manchester coding
        load('sequence_pam4.mat');% randomly generated binary sequence
        signal = sequence_pam4(1:ceil(net_pack_length/2)*2*net_pack_num);% 
        signal = reshape(signal,2,[]);
        signal_pam = bi2de(signal','left-msb');
        signal_temp = 3-signal_pam;
        signal = [signal_pam';signal_temp'];% Manchester coding
        signal = signal(:);
        signal = reshape(signal,[],net_pack_num);% Packetization
        signal = signal(1:net_pack_length,:);
        signal = repmat(signal,3,1);% Packet repetation
        sequence_R = signal(:);
        out = reshape(sequence_R,[],net_pack_num);
    end

    %% Header insertion
    header = [ 3 3 3 3, 3 0 3 0, 3 0 3 0 ];% header for synchronization and timing recovery
    header = header(:);
    header = repmat(header,1,net_pack_num);

    No = [0 0 3 3, 0 3 3 0, 3 0 0 3, 3 3 0 0];% header to indicate the order in which the packets appear
    No = repmat(No,1,ceil(net_pack_num/4));
    No = [3 0 3 0 , No];
    No = No(1:net_pack_num*4);
    No = reshape(No,4,[]);

    RPT_1 = [3 0 3 0];% header to indicate the order in which the repeatitive packets appear
    RPT_1 = RPT_1(:);
    RPT_1 = repmat(RPT_1,1,net_pack_num);
    RPT_2 = [0 3 0 3];
    RPT_2 = RPT_2(:);
    RPT_2 = repmat(RPT_2,1,net_pack_num);
    RPT_3 = [3 0 0 3];
    RPT_3 = RPT_3(:);
    RPT_3 = repmat(RPT_3,1,net_pack_num);

    sequence_R = reshape(out,[],net_pack_num);
    sequence_R = [header;RPT_1;No;sequence_R(1:net_pack_length,:);...
                  header;RPT_2;No;sequence_R(net_pack_length+1:2*net_pack_length,:);...
                  header;RPT_3;No;sequence_R(2*net_pack_length+1:end,:)];
    sequence_R = sequence_R(:)';
    sequence_R = [sequence_R];
    %% save
%     r_bit = [r_bit_pre;r_bit_post];

%     save('sequence_save_v2.mat','sequence_R',...
%         'pack_length','symbol_rate','r_bit','block_bit_length');
end