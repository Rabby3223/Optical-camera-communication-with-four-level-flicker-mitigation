function [scaling_factor] = scaling_function(filename,img_num)
% Author: LIU Liqiong
% Email: liuliqiongrabby@gmail.com
% Institution: Department of Information Engineering, 
% The Chinese University of Hong Kong.
% Created: March 2021
% Modified: November 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    scaling_factor = zeros(960,1); % each row signal has a scaling factor
    for index =0:img_num-1 % scaling factor is derived by averaging signals from 300 frames
        filename_img = strcat(filename,'/',num2str(index),'.csv');% read csv file 
        [data] = textread(filename_img,'%s');
        [row_num,~] = size(data);
        if row_num~=0 

            data_temp = char(data(1));
            data_decoded_temp = base64decode(data_temp);% decode signals in base64 format
            [~,col_num] = size(data_decoded_temp);
            data_decoded = uint8(zeros(row_num-1,col_num));
            for i = 1:row_num-1
                data_temp = char(data(i));
                data_decoded(i,:) = base64decode(data_temp);
            end
            data_temp = char(data(row_num));
            data_decoded_end = base64decode(data_temp);

            data_decoded = data_decoded';
            data_decoded = data_decoded(:);
            data_decoded = [data_decoded;data_decoded_end'];
            data_decoded = double(data_decoded);
            scaling_factor = scaling_factor+data_decoded./img_num;
        end
    end
% save('scaling_factor.mat','scaling_factor');
end