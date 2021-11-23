function [signal_frame_head,No_index,signal_frame] = packet_reconstruction_PAM4(sig_decoded,P,...
        signal_frame_head,pack_length,...
            No_index,image_index,flicker_mitigation,r_bit,block_bit_length,signal_frame)
% This function is for packet reconstruction. 
% This function will first check the number No and RPT.A complete packet will be generated from one or two images. 
% Input:
%       sig_decoded: Decoded PAM4 symbols.
%       P: the index for header in the sig_decoded.
%       signal_frame_head: Incomplete packet from the last image, which may
%       be used to generate a complete packet with signals in current image. 
%       pack_length: number of symbols in a packet
%       No_index: No of signal_frame_head if not empty
%       image_index: Image index
%       flicker_mitigation: 0 for no flicker mitigation, 1 for Manchester,
%       2 for proposed scheme. 
%       re_bit: redundant bits
%       block_bit_length: number of binary bits in a block
%       signal_frame: complete packets
% Output: 
%       signal_frame_head: Incomplete packet in current image.
%       header_index: The No of signal_frame_head if not empty.
%       signal_frame: complete packets.
% Created: Mar 2021
% Modified: Nov 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Header detection
    % Find all all headers
    header_length = 16;% header length excluding the length of rpt
    net_pack_length = pack_length-16;
    index_header1 = P;% indexes of headers
    if P+pack_length+header_length < length(sig_decoded)
        index_header1 = [P,P+pack_length];
    end
    if P-pack_length>0
        index_header1 = [P-pack_length,P];
    end

    % Find the RPT index for all the headers
    P_temp = index_header1';
    if isempty(P_temp)
        % disp(strcat(num2str(image_index),': warning: no header is found'));
        signal_frame_head = [];
        return;% If no header has been found, process signals from the next image
    end
    RPT_index_temp = zeros(length(P_temp),1);% PRT index
    for i = 1:length(P_temp)
        if P_temp(i)+20<length(sig_decoded)
            
            RPT = sig_decoded(P_temp(i)+12)*1000+sig_decoded(P_temp(i)+13)*100+sig_decoded(P_temp(i)+14)*10+sig_decoded(P_temp(i)+15)*1;
            switch RPT
                case 3030
                    RPT_index_temp(i) = 0;
                    continue;
                case 303
                    RPT_index_temp(i) = 1;
                    continue;
                case 3003
                    RPT_index_temp(i) = 2;
                    continue;
                otherwise
                    continue;
            end
        end
    end
    % Sort the header index and PRT index
    a = [P_temp,RPT_index_temp];
    a = sortrows(a,1);
    P_temp = a(:,1);
    RPT_index_temp = a(:,2);
    % The distance between the two headers is equal to the packet length
    % Wrong header indexes can be removed by checking the distance between
    % the two headers.
    P = [];
    RPT_index = [];
    if numel(P_temp)==1
        P=P_temp;
        RPT_index = RPT_index_temp;
    else
        for i = 1:length(P_temp)-1
            for j = i+1:length(P_temp)
                if P_temp(j) - P_temp(i)==pack_length
                    P = [P_temp(i);P_temp(j)];
                    RPT_index = [RPT_index_temp(i),RPT_index_temp(j)];
                else

                end
            end
        end
    end
    temp = [];
    for i = 1:length(P)-1
        if P(i+1)==P(i)
            temp = [temp,i];
        end
    end
    P(temp) = [];
    RPT_index(temp) = [];
    % If no header is effective, process signals from the next image
    if isempty(P)
        % disp(strcat(num2str(image_index),': warning: no header is found'));
        signal_frame_head = [];
        return;
    end

%% Packet reconstruction

    reuse_index=0;
    % Case 1: If one packet is captured without complete header and its RPT index
    % is 3, the payload will be saved in signal_frame.
    if P(1)>net_pack_length && RPT_index(1)==0
            p_temp = P(1)-net_pack_length;
            header = [3 3 3 3,3 0 3 0,3 0 3 0,3 0 0 3];
            current_No = sig_decoded(p_temp)*1000+sig_decoded(p_temp+1)*100+sig_decoded(p_temp+2)*10+sig_decoded(p_temp+3)*1;
            switch current_No
                case 33
                    current_No = 0;
                case 330
                    current_No = 1;
                case 3003
                    current_No = 2;
                case 3300
                    current_No = 3;
                case 3030
                    current_No = 5;
                otherwise
                    current_No = 4;
            end
            signal_current = [header,sig_decoded(P(1)-net_pack_length:P(1)-1)]; 
            if flicker_mitigation==2
                signal_current = flicker_mitigation_decoding_PAM4( signal_current,r_bit,block_bit_length,pack_length,'packet');
            end
            signal_frame = [signal_frame;signal_current];
            signal_frame_head = [];
    end
    
    % The incomplete packets from the last image and current image may form a complete packet. 
    if ~isempty(signal_frame_head)
        overlap_length = P(1)-1+numel(signal_frame_head)-pack_length;% Calculate the number of overlapped symbols in the two incomplete packets.
        % The packet cannot be formed if the two incomplete packets do not
        % have overlapped symbols
        if RPT_index(1)==0 && overlap_length<0
            signal_frame_head = [];
            %disp(strcat(num2str(image_index),': cannot form a packet'));
            return;
        else
            if RPT_index(1)==2 && overlap_length<0
%                 disp('danger: packet loss possible ',num2str(image_index));
                signal_frame_head = [];
            else
                % If the length of symbols before the first header is
                % larger than the packet length, the existing complete
                % packet should be saved later.
                
                % The index of current incomplet packet and overlapping
                % length should be recalculated.
                % The No index should be checked to make sure the two
                % incomplete packets from the same transmitted packet.
                if RPT_index(1)==1 && overlap_length>0 
                    if P(1)>pack_length
                        
                        P_temp = P(1) - pack_length;
                        overlap_length = P_temp-1+numel(signal_frame_head)-pack_length;
                        current_No = sig_decoded(P_temp+header_length)*1000+sig_decoded(P_temp+header_length+1)*100+sig_decoded(P_temp+header_length+2)*10+sig_decoded(P_temp+header_length+3)*1;
                        switch current_No
                            case 33
                                current_No = 0;
                            case 330
                                current_No = 1;  
                            case 3003
                                current_No = 2;
                            case 3300
                                current_No = 3; 
                            case 3030
                                current_No = 5; % Indicate the first transmitted packet.
                            otherwise
                                current_No = 4; % No effective index can be detected. 
                        end
                        if current_No~=4 && current_No~=5
                            current_No = current_No-1;
                            curent_index = mod(current_No,4);
                        end
                        if overlap_length>0 && (current_No==No_index||No_index==5)
                            % Combine the two incomplete packets. Two cases
                            % depend on whether the overlapping length is
                            % odd or even.
                            if mod(overlap_length,2)==1
                                diff = length(signal_frame_head) - floor(overlap_length/2) - floor( (length(signal_frame_head) - floor(overlap_length/2))/4 )*4;
                                seq_head = signal_frame_head(1:end-floor(overlap_length/2)-diff);
                                if flicker_mitigation==2
                                    seq_head = flicker_mitigation_decoding_PAM4( seq_head,r_bit,block_bit_length,pack_length,'head');
                                end
                                seq_tail = sig_decoded(floor(overlap_length/2)+2-diff:P_temp-1);
                                if flicker_mitigation==2
                                    seq_tail = flicker_mitigation_decoding_PAM4( seq_tail,r_bit,block_bit_length,pack_length,'tail');
                                end
                                signal_current = [seq_head,seq_tail];                      
                                signal_frame_head = [];% Clear the signal_frame_head
                                if numel(signal_current)==pack_length
                                      signal_frame = [signal_frame;signal_current];% Add the new pacekt to signal_frame
                                end
                            else
                                diff = length(signal_frame_head) - overlap_length/2 - floor( (length(signal_frame_head) - overlap_length/2)/4 )*4;
                                seq_head = signal_frame_head(1:end-overlap_length/2-diff);
                                if flicker_mitigation==2
                                    seq_head = flicker_mitigation_decoding_PAM4( seq_head,r_bit,block_bit_length,pack_length,'head');
                                end
                                seq_tail = sig_decoded(overlap_length/2+1-diff:P_temp-1);
                                if flicker_mitigation==2
                                    seq_tail = flicker_mitigation_decoding_PAM4( seq_tail,r_bit,block_bit_length,pack_length,'tail');
                                end
                                signal_current  = [seq_head,seq_tail];                                                        
                                signal_frame_head = [];% Clear the signal_frame_head
                                if numel(signal_current)==pack_length
                                      signal_frame = [signal_frame;signal_current];% Add the new pacekt to signal_frame
                                end
                            end
                        end
                        % Save the existing complete packet
                        header = [3 3 3 3,3 0 3 0,3 0 3 0,0 0 3 3];
                        signal_current = [header,Decoded_R1(P(1)-net_pack_length:P(1)-1)];
                        signal_frame = [signal_frame;signal_current];
                        reuse_index = 1;% Indicate that the symbols berfore the first header has been used for pacekt reconstruction.
                    end
                else 
                    % if there is no complete packet before the first
                    % detected header
                     current_No = sig_decoded(P(1)+header_length)*1000+sig_decoded(P(1)+header_length+1)*100+sig_decoded(P(1)+header_length+2)*10+sig_decoded(P(1)+header_length+3)*1;
                     switch current_No
                            case 33
                                current_No = 0;
                            case 330
                                current_No = 1;
                            case 3003
                                current_No = 2;
                            case 3300
                                current_No = 3;
                            case 3030
                                current_No = 5;% Indicate the first transmitted packet.
                            otherwise
                                current_No = 4;% No effective index can be detected. 
                     end
                     % Check whether the two incomplete packets from the same
                     % transmitted pacekt
                     b1 = (current_No~=4) &&  ((RPT_index(1) == 0 && (mod((current_No - 1 + 4),4) == No_index)) || (RPT_index(1) == 2 && mod(current_No,4)== No_index));
                     b2 = (No_index==5);
                     if overlap_length>0 && (b1||b2)
                        % Combine the two incomplete packets. Two cases
                        % depend on whether the overlapping length is
                        % odd or even.
                        if mod(overlap_length,2)==1
                            diff = length(signal_frame_head) - floor(overlap_length/2) - floor( (length(signal_frame_head) - floor(overlap_length/2))/4 )*4;
                            seq_head = signal_frame_head(1:end-floor(overlap_length/2)-diff);
                            if flicker_mitigation==2
                                seq_head = flicker_mitigation_decoding_PAM4( seq_head,r_bit,block_bit_length,pack_length,'head');
                            end
                            seq_tail = sig_decoded(floor(overlap_length/2)+2-diff:P(1)-1);
                            if flicker_mitigation==2
                                seq_tail = flicker_mitigation_decoding_PAM4( seq_tail,r_bit,block_bit_length,pack_length,'tail');
                            end
                            signal_current = [seq_head,seq_tail];
                            if numel(signal_current)==pack_length
                                signal_frame = [signal_frame;signal_current];% Add the new pacekt to signal_frame
                            end
                        else
                            diff = length(signal_frame_head) - overlap_length/2 - floor( (length(signal_frame_head) - overlap_length/2)/4 )*4;
                            seq_head = signal_frame_head(1:end-overlap_length/2-diff);
                            if flicker_mitigation==2
                                seq_head = flicker_mitigation_decoding_PAM4( seq_head,r_bit,block_bit_length,pack_length,'head');
                            end
                            seq_tail = sig_decoded(overlap_length/2+1-diff:P(1)-1);
                            if flicker_mitigation==2
                                seq_tail = flicker_mitigation_decoding_PAM4( seq_tail,r_bit,block_bit_length,pack_length,'tail');
                            end
                            signal_current  = [seq_head,seq_tail];
                            if numel(signal_current)==pack_length
                                signal_frame = [signal_frame;signal_current];% Add the new pacekt to signal_frame
                            end
                        end
                        signal_frame_head = [];% Clear the signal_frame_head 
                        reuse_index = 1;% Indicate that the symbols berfore the first header has been used for pacekt reconstruction.
                    else
                        signal_frame_head = [];% Clear the signal_frame_head
                    end
                end
            end
        end
    end
    signal_frame_head = [];% Clear the signal_frame_head
    
    % Save the complete pacekt received in the current image or combine two
    % incomplete packet from the current image to form a complete one. 
    if (P(1) + pack_length -1 < numel(sig_decoded))
             % Save the complete pacekt received in the current image
             current_No = sig_decoded(P(1)+header_length)*1000+sig_decoded(P(1)+header_length+1)*100+sig_decoded(P(1)+header_length+2)*10+sig_decoded(P(1)+header_length+3)*1;
             switch current_No
                 case 33
                    current_No = 0;
                case 330
                    current_No = 1;
                case 3003
                    current_No = 2;
                case 3300
                    current_No = 3;
                case 3030
                    current_No = 5;% Indicate the first transmitted packet.
                otherwise
                    current_No = 4;% No effective index can be detected. 
             end
             signal_current = sig_decoded(P(1):P(1)+pack_length-1);
             if flicker_mitigation==2
                 signal_current = flicker_mitigation_decoding_PAM4( signal_current,r_bit,block_bit_length,pack_length,'packet');
             end
             signal_frame = [signal_frame;signal_current];% Add the new pacekt to signal_frame
    else 
         if( (P(1)+20)<length(sig_decoded) )
             % Combine two incomplete packet from the current image to form a complete one
             signal_frame_head = sig_decoded(P(1):end);
             No_index =  sig_decoded(P(1)+header_length)*1000+sig_decoded(P(1)+header_length+1)*100+sig_decoded(P(1)+header_length+2)*10+sig_decoded(P(1)+header_length+3)*1;
             switch No_index
                case 33
                    No_index = 0;
                case 330
                    No_index = 1;
                case 3003
                    No_index = 2;
                case 3300
                    No_index = 3;
                case 3030
                    No_index = 5;% Indicate the first transmitted packet.
                otherwise
                    No_index = 4;% No effective index can be detected. 
             end
             overlap_length = P(1)-1+numel(signal_frame_head)-pack_length;
             % Check whether the two incomplete packets from the same
             % transmitted pacekt
             if ( RPT_index(1)==1 || RPT_index(1)==2 ) ...
                     && overlap_length>0 && reuse_index==0 
                  % Combine the two incomplete packets. Two cases
                  % depend on whether the overlapping length is
                  % odd or even.
                 if mod(overlap_length,2)==1 && length(signal_frame_head) - floor(overlap_length/2)>0
                     diff = length(signal_frame_head) - floor(overlap_length/2) - floor( (length(signal_frame_head) - floor(overlap_length/2))/4 )*4;   
                     signal_current = [signal_frame_head(1:end-floor(overlap_length/2)-diff),sig_decoded(floor(overlap_length/2)+2-diff:P(1)-1)];
                     if numel(signal_current)==pack_length
                         if flicker_mitigation==2
                             signal_current = flicker_mitigation_decoding_PAM4( signal_current,r_bit,block_bit_length,pack_length,'packet');
                         end
                         signal_frame = [signal_frame;signal_current];
                     end
                 else if mod(overlap_length,2)==0 && length(signal_frame_head) - overlap_length/2>0
                         diff = length(signal_frame_head) - overlap_length/2 - floor( (length(signal_frame_head) - overlap_length/2)/4 )*4;
                         seq_head = signal_frame_head(1:end-overlap_length/2-diff);
                         if flicker_mitigation==2
                             seq_head = flicker_mitigation_decoding_PAM4(seq_head,r_bit,block_bit_length,pack_length,'head');
                         end
                         seq_tail = sig_decoded(overlap_length/2+1-diff:P(1)-1);
                         if flicker_mitigation==2
                             seq_tail = flicker_mitigation_decoding_PAM4(seq_tail,r_bit,block_bit_length,pack_length,'tail');
                         end
                         signal_current  = [seq_head,seq_tail];
                         if numel(signal_current)==pack_length
                            signal_frame = [signal_frame;signal_current];
                         end                  
                     end
                 end
                 signal_frame_head = [];% Clear the signal_frame_head
             else
                 if  (RPT_index(1)==1||RPT_index(1)==2) && reuse_index==1
                     signal_frame_head = [];% Clear the signal_frame_head
                 end
             end
         end
    end
    %  Check whether there is incomplete packet that has not be used to form a complete packet           
    if numel(P)==2 && numel(sig_decoded)-P(2)+1>20      
         signal_frame_head = sig_decoded(P(2):end); 
         No_index =  sig_decoded(P(2)+header_length)*1000+sig_decoded(P(2)+header_length+1)*100+sig_decoded(P(2)+header_length+2)*10+sig_decoded(P(2)+header_length+3)*1;
         switch No_index
             case 33
                No_index = 0;
            case 330
                No_index = 1;
            case 3003
                No_index = 2;
            case 3300
                No_index = 3;
            case 3030
                No_index = 5;
            otherwise
                No_index = 4; 
                signal_frame_head = [];% Clear the signal_frame_head
         end 
         if RPT_index(2)==1 || RPT_index(2)==2
            signal_frame_head = [];
         end
    end
    %  Clear signal_frame_head if it only includes header,
    if length(signal_frame_head)==header_length
        signal_frame_head = [];% Clear the signal_frame_head
    end
    %  Clear signal_frame_head if no effective No_index can be detected,
    if No_index==4
        signal_frame_head = [];% Clear the signal_frame_head
    end
end

