function [ signal_output ] = flicker_mitigation_decoding_PAM4( signal,r_bit,block_bit_length,pack_length,CASE )
%% parameter setting
    r_bit_length = length(r_bit);% number of redundant bits in a block
    r_pam_length = ceil(r_bit_length/2/2)*2;% number of redundant symbols in a block
    r_bit_pre = r_bit(1:r_bit_length/2);% redundant bits added at the beginning of the block
    r_bit_post = r_bit(r_bit_length/2+1:end);% redundant bits added at the end of the block
    net_pack_length = pack_length-20;% Number of symbols in a packet excluding header
    block_num = floor(net_pack_length/(block_bit_length/2));% number of regular blocks in a packet
    block_bit_length_temp = (net_pack_length - block_num*block_bit_length/2)*2;% number of bits in the last block
    
    signal_output = signal;
%% Decoding
    switch CASE % Three cases for decoding:'head','tail','packet'
        case 'head' % If the incomplete packet is the head of a packet
            if length(signal)>20
                signal = signal(21:end);% Remove header.
            else
                signal_output = signal;% Stop decdoing if only header is included.
                return;
            end
            % If the incomplete packet include less than a block.
            if length(signal)+1>r_pam_length/2 && length(signal)<block_bit_length/2
                temp = signal;% PAM4 symbols to be decoded
                % Estimate whether the packet is inverted or not based on
                % the redundant bits r_bit_pre 
                temp_pam = [temp(1:r_pam_length/2)];% Extract redundant symbols
                temp_bit = de2bi(temp_pam','left-msb',2);% Convert PAM4 symbols to binary bits
                temp_bit = temp_bit';temp_bit = temp_bit(:);
                temp_bit = temp_bit(1:r_bit_length/2);
                BER = sum(abs(temp_bit-r_bit_pre))/(r_bit_length/2);
                if(BER==1)% The inverted redundant bits are received
                     signal = 3-temp;% Invert the symbols
                end
                if(BER>0 && BER<1)% Error occurs in the redundant bits
                     signal = ones(1,length(signal))*(-1);% Set all the signals as -1
                end
            % If the incomplete pacekt inlude one or more than one block,    
            else
                % Decode for the complete blocks.The r_bit_pre and
                % r_bit_post can both be used for inversion estimation
                for i = 1:floor(length(signal)/(block_bit_length/2))
                    temp = signal((i-1)*block_bit_length/2+1:i*block_bit_length/2);% PAM4 symbols to be decoded
                    temp_pam = [temp(1:r_pam_length/2),temp(end-r_pam_length/2+1:end)];% Extract redundant symbols
                    temp_bit = de2bi(temp_pam','left-msb',2);% Convert PAM4 symbols to binary bits
                    temp_bit = temp_bit';temp_bit = temp_bit(:);
                    temp_bit = [temp_bit(1:r_bit_length/2);temp_bit(end-r_bit_length/2+1:end)];
                    BER = sum(abs(temp_bit-r_bit))/r_bit_length;
                    if(BER==1)% The inverted redundant bits are received
                        signal((i-1)*block_bit_length/2+1:i*block_bit_length/2) = 3-temp;% Invert the symbols
                    end
                    if(BER>0 && BER<1)% Error occurs in the redundant bits
                        signal((i-1)*block_bit_length/2+1:i*block_bit_length/2) = -1;% Set all the signals as -1
                    end
                end
                i = floor(length(signal)/(block_bit_length/2))+1;
                % Decode for the last incomplete block based on r_bit_pre.
                if (i-1)*block_bit_length/2+r_pam_length/2 < length(signal)+1
                    temp = signal((i-1)*block_bit_length/2+1:end);% PAM4 symbols to be decoded
                    temp_pam = [temp(1:r_pam_length/2)];% Extract redundant symbols
                    temp_bit = de2bi(temp_pam','left-msb',2);% Convert PAM4 symbols to binary bits
                    temp_bit = temp_bit';temp_bit = temp_bit(:);
                    temp_bit = [temp_bit(1:r_bit_length/2)];
                    BER = sum(abs(temp_bit-r_bit_pre))/(r_bit_length/2);
                    if(BER==1)% The inverted redundant bits are received
                         signal((i-1)*block_bit_length/2+1:end) = 3-temp;% Invert the symbols
                    end
                    if(BER>0 && BER<1)% Error occurs in the redundant bits
                         signal((i-1)*block_bit_length/2+1:end) = -1;% Set all the signals as -1
                    end
                end
            end
            signal_output = [signal_output(1:20),signal];
        case 'tail'% If the incomplete packet is the tail of a packet
            % Estimate whether the last block with length block_bit_length_temp is inverted or not based on
            % the redundant bits r_bit_post
            if length(signal)+1>r_pam_length/2 && length(signal)<block_bit_length_temp/2
                temp = signal;% PAM4 symbols to be decoded
                temp_pam = [temp(end-block_bit_length_temp/2+1:end)];% Extract redundant symbols
                temp_bit = de2bi(temp_pam','left-msb',2);% Convert PAM4 symbols to binary bits
                temp_bit = temp_bit';temp_bit = temp_bit(:);
                temp_bit = temp_bit(end-r_bit_length/2+1:end);
                BER = sum(abs(temp_bit-r_bit_post))/(r_bit_length/2);
                if(BER==1)% The inverted redundant bits are received
                     signal = 3-temp;% Invert the symbols
                end
                if(BER>0 && BER<1)% Error occurs in the redundant bits
                     signal = ones(1,length(signal))*(-1);% Set all the signals as -1
                end
            % If more blocks are captured,    
            else
                end_index = length(signal)-block_bit_length_temp/2+1;
                start_index = end_index-floor((end_index-1)/(block_bit_length/2))*block_bit_length/2;
                index = start_index:block_bit_length/2:end_index;
                % Decode for the complete regular blocks.The r_bit_pre and
                % r_bit_post can both be used for inversion estimation
                for i = 1:length(index)-1
                        index_temp = index(i);
                        temp = signal(index_temp:index_temp+block_bit_length/2-1);% PAM4 symbols to be decoded
                        temp_pam = [temp(1:r_pam_length/2),temp(end-r_pam_length/2+1:end)];% Extract redundant symbols
                        temp_bit = de2bi(temp_pam','left-msb',2);% Convert PAM4 symbols to binary bits
                        temp_bit = temp_bit';temp_bit = temp_bit(:);
                        temp_bit = [temp_bit(1:r_bit_length/2);temp_bit(end-r_bit_length/2+1:end)];
                        BER = sum(abs(temp_bit-r_bit))/r_bit_length;
                        if(BER==1)% The inverted redundant bits are received
                            signal(index_temp:index_temp+block_bit_length/2-1) = 3-temp;% Invert the symbols
                        end
                        if(BER>0 && BER<1)% Error occurs in the redundant bits
                            signal(index_temp:index_temp+block_bit_length/2-1) = -1;% Set all the signals as -1
                        end
                end
                % Decode for the last block.The r_bit_pre and
                % r_bit_post can both be used for inversion estimation
                temp = signal(index(end):end);% PAM4 symbols to be decoded
                temp_pam = [temp(1:r_pam_length/2),signal(end-r_pam_length/2+1:end)];% Extract redundant symbols
                temp_bit = de2bi(temp_pam','left-msb',2);% Convert PAM4 symbols to binary bits
                temp_bit = temp_bit';temp_bit = temp_bit(:);
                temp_bit = [temp_bit(1:r_bit_length/2);temp_bit(end-r_bit_length/2+1:end)];
                BER = sum(abs(temp_bit-r_bit))/r_bit_length;
                if(BER==1)% The inverted redundant bits are received
                    signal(index(end):end) = 3-temp;% Invert the symbols
                end
                if(BER>0 && BER<1)% Error occurs in the redundant bits
                    signal(index(end):end) = -1;% Set all the signals as -1
                end
                % Decode for the fist imcomplete block.The r_bit_post can 
                % both be used for inversion estimation
                if index(1)>r_pam_length/2
                    temp = signal(1:index(1)-1);% PAM4 symbols to be decoded
                    temp_pam = [temp(end-r_pam_length/2+1:end)];% Extract redundant symbols
                    temp_bit = de2bi(temp_pam','left-msb',2);% Convert PAM4 symbols to binary bits
                    temp_bit = temp_bit';temp_bit = temp_bit(:);
                    temp_bit = [temp_bit(end-r_bit_length/2+1:end)];
                    BER = sum(abs(temp_bit-r_bit_post))/(r_bit_length/2);
                    if(BER==1)% The inverted redundant bits are received
                        signal(1:index(1)-1) = 3-temp;% Invert the symbols
                    end
                    if(BER>0 && BER<1)% Error occurs in the redundant bits
                        signal(1:index(1)-1) = -1;% Set all the signals as -1
                    end
                end
            end
            signal_output = signal;
            case 'packet'
                signal = signal(21:end);
                index = 1:block_bit_length/2:length(signal);
                % Decode for the complete blocks.The r_bit_pre and
                % r_bit_post can both be used for inversion estimation
                for i = 1:length(index)-1
                    index_temp = index(i);
                    temp = signal(index_temp:index_temp+block_bit_length/2-1);% PAM4 symbols to be decoded
                    temp_pam = [temp(1:r_pam_length/2),temp(end-r_pam_length/2+1:end)];% Extract redundant symbols
                    temp_bit = de2bi(temp_pam','left-msb',2);% Convert PAM4 symbols to binary bits
                    temp_bit = temp_bit';temp_bit = temp_bit(:);
                    temp_bit = [temp_bit(1:r_bit_length/2);temp_bit(end-r_bit_length/2+1:end)];
                    BER = sum(abs(temp_bit-r_bit))/r_bit_length;
                    if(BER==1)% The inverted redundant bits are received
                        signal(index_temp:index_temp+block_bit_length/2-1) = 3-temp;% Invert the symbols
                    end
                    if(BER>0 && BER<1)% Error occurs in the redundant bits
                        signal(index_temp:index_temp+block_bit_length/2-1) = -1;% Set all the signals as -1
                    end
                end
                temp = signal(index(end):end);% PAM4 symbols to be decoded
                temp_pam = [temp(1:r_pam_length/2),temp(end-r_pam_length/2+1:end)];% Extract redundant symbols
                temp_bit = de2bi(temp_pam','left-msb',2);% Convert PAM4 symbols to binary bits
                temp_bit = temp_bit';temp_bit = temp_bit(:);
                temp_bit = [temp_bit(1:r_bit_length/2);temp_bit(end-r_bit_length/2+1:end)];
                BER = sum(abs(temp_bit-r_bit))/r_bit_length;
                if(BER==1)% The inverted redundant bits are received
                    signal(index(end):end) = 3-temp;% Invert the symbols
                end
                if(BER>0 && BER<1)% Error occurs in the redundant bits
                    signal(index(end):end) = -1;% Set all the signals as -1
                end
                signal_output = [signal_output(1:20),signal];
    end

    
end

