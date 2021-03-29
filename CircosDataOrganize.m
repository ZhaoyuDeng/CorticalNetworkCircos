function [] = CircosDataOrganize(working_dir,RawDataCircos,P_threshold,link_mode)
% FORMAT [] = CircosDataOrganize(working_dir,RawDataCircos,P_threshold,link_mode)
% Data organization for Circos plot via the format of txt
% Input:
%   working_dir - working directory that generate scripts and Circos figure
%   RawDataCircos - raw data according to format that manual defined
%   P_threshold - threshold of P value, to filter
%   colorbar_roundn_sens - sensitivity of colorbar range
%   link_mode - select link mode: 1-even link width, 2-ratio link width
%__________________________________________________________________________
% Written by DENG Zhao-Yu 210322 for DPARBI.
% Institute of Psychology, Chinese Academy of Sciences
% dengzy@psych.ac.cn
%__________________________________________________________________________
%%
% change working directory
cd(working_dir);

% define variables
colorbar_roundn_sens = -2;
link_transparency = 0.5;

% define the amount of networks and regions
region_num = length(RawDataCircos.ElementLabel);
network_num = length(RawDataCircos.HigherOrderNetworkLabel(:,1));
% width of a region band
if link_mode == 1
    spacing = region_num; % even link width
elseif link_mode == 2
    spacing = 100; % ratio link width, default 100
elseif link_mode == 3
    spacing = []; %ratio link width varies with sum of correlation
end

% networks and regions
tab = tabulate(RawDataCircos.HigherOrderNetworkIndex);
network_region_num = tab(:,2);
networks = RawDataCircos.HigherOrderNetworkLabel(:,2);
regions = RawDataCircos.ElementLabel(:,2);

% generate correlation matrix for links, filter threshold
corr_mat = RawDataCircos.P_Corrected < P_threshold; % filter the correlation matrix of regions that P less then 0.01
[cor_row,cor_col] = find(triu(corr_mat)); % withdraw upper triangle matrix
cor_num = length(cor_row); % number of correlation
plot_mat = zeros(length(cor_row),1); % initialize matrix that store values that plots
for k = 1:length(cor_row)
    plot_mat(k) = RawDataCircos.Matrix(cor_row(k),cor_col(k)); % store values in matrix
end
plot_mat_max = max(abs(min(plot_mat)),abs(max(plot_mat)));
ratio_mat = plot_mat/roundn(plot_mat_max,colorbar_roundn_sens); % normalization, sensitivity 0.01

% calculate correlation for plot
link_index_logical = false(cor_num,cor_num);
link_num = zeros(cor_num,1); 
link_width_ratio = zeros(region_num,region_num); % initialize ratio of regions' link width
for k = 1:region_num
    link_index_logical(:,k) = cor_row==k|cor_col==k;
    link_num(k) = sum(link_index_logical(:,k)); % regions' links amount
    link_width_ratio(k,1:link_num(k)) = abs(ratio_mat(link_index_logical(:,k))');
end
link_width_ratio_sum = sum(link_width_ratio,2); % calculate sum of ratios

% calculate the color of links
color_mat = zeros(cor_num,3); % initialize matrix that store color of links
for k = 1:length(cor_row)
    if ratio_mat(k) > 0
        color_mat(k,:) = fix([255,255-ratio_mat(k)*255,255-ratio_mat(k)*255]); % the bigger positive number is, the more red
    else
        color_mat(k,:) = fix([255+ratio_mat(k)*255,255+ratio_mat(k)*255,255]); % the smaller negative number is, the more blue
    end
end


% write data of networks and regions
fid = fopen('CircosInput1_band.txt','w');
% describe external networks, FORMAT: chr - ID label start end attribute
if link_mode==1 || link_mode==2
    % isometry band
    for k = 1:network_num
        fprintf(fid,'chr - %s %s ',['net',num2str(k)],cell2mat(networks(k)));
        fprintf(fid,'%u %u %s',0,network_region_num(k)*spacing,['chr',num2str(k)]); %color may not needed
        fprintf(fid,'\n');
    end
elseif link_mode==3
    % not isometry band
    link_width_ratio100 = floor(link_width_ratio*100);
    band_width = sum(link_width_ratio100,2);
    band_width(band_width==0) = floor(min(band_width(band_width~=0)));
    network_width = zeros(3,1); % initialize
    for k = 1:network_num
        network_width(k) = sum(band_width(RawDataCircos.HigherOrderNetworkIndex==k));
    end
    for k = 1:network_num
        fprintf(fid,'chr - %s %s ',['net',num2str(k)],cell2mat(networks(k)));
        fprintf(fid,'%u %u %s',0,network_width(k),['chr',num2str(k)]); %color may not needed
        fprintf(fid,'\n');
    end
end
% describe internal bands, FORMAT: band ID label label start end attribute
if link_mode==1 || link_mode==2
    index = 1;
    for k = 1:network_num
        for l = 1:network_region_num(k)
            fprintf(fid,'band %s %s %s ',['net',num2str(k)],cell2mat(regions(index)),cell2mat(regions(index)));
            fprintf(fid,'%u %u %s',(l-1)*spacing,l*spacing,['chr',num2str(index)]);
            fprintf(fid,'\n');
            index = index + 1;
        end
    end
elseif link_mode==3
    index = 1;
    for k = 1:network_num
        former_sum = 0;
        for l = 1:network_region_num(k)
            fprintf(fid,'band %s %s %s ',['net',num2str(k)],cell2mat(regions(index)),cell2mat(regions(index)));
            fprintf(fid,'%u %u %s',former_sum,former_sum+band_width(index),['chr',num2str(index)]);
            fprintf(fid,'\n');
            former_sum = former_sum + band_width(index);
            index = index + 1;
        end
    end
end
fclose(fid);

% write data of band labels
fid = fopen('CircosInput2_label.txt','w');
% label karyotype band, FORMAT: ID start end label
if link_mode==1 || link_mode==2
    index = 1;
    for k = 1:network_num
        for l = 1:network_region_num(k)
            fprintf(fid,'%s %u %u %s',['net',num2str(k)],(l-1)*spacing,l*spacing,cell2mat(regions(index)));
            fprintf(fid,'\n');
            index = index + 1;
        end
    end
elseif link_mode==3
    index = 1;
    for k = 1:network_num
        former_sum = 0;
        for l = 1:network_region_num(k)
            fprintf(fid,'%s %u %u %s',['net',num2str(k)],former_sum,former_sum+band_width(index),cell2mat(regions(index)));
            fprintf(fid,'\n');
            former_sum = former_sum + band_width(index);
            index = index + 1;
        end
    end
end
fclose(fid);

% write data of links
fid = fopen('CircosInput3_link.txt','w'); 
% describe links, FORMAT: Chromosome1 Start1 End1 Chromosome2 Start2 End2 Attributes
if link_mode == 1 % even link width mode
    for k = 1:cor_num
        % calculate chromsome1's network and region
        cor_row_net = RawDataCircos.HigherOrderNetworkIndex(cor_row(k));
        cor_row_reg = cor_row(k);
        if cor_row_net > 1
            for l = 1:cor_row_net-1
                cor_row_reg = cor_row_reg - network_region_num(l);
            end
        end
        % calculate chromsome2's network and region
        cor_col_net = RawDataCircos.HigherOrderNetworkIndex(cor_col(k));
        cor_col_reg = cor_col(k);
        if cor_col_net > 1
            for l = 1:cor_col_net-1
                cor_col_reg = cor_col_reg - network_region_num(l);
            end
        end
        % calculate start and end
        cor_row_start = (cor_row_reg-1)*spacing+(cor_col_reg-1);
        cor_row_end = cor_row_start+1;
        cor_col_start = (cor_col_reg-1)*spacing+(cor_row_reg-1);
        cor_col_end = cor_col_start+1;
        % print on txt according to format
        fprintf(fid,'net%u %u %u ',cor_row_net,cor_row_start,cor_row_end);
        fprintf(fid,'net%u %u %u ',cor_col_net,cor_col_start,cor_col_end);
        fprintf(fid,'color=%u,%u,%u',color_mat(k,1),color_mat(k,2),color_mat(k,3));
        fprintf(fid,'\n');
    end
elseif link_mode == 2 % ratio link width mode
    link_width_percent = zeros(region_num,region_num);
    % normalize link width percent
    for k = 1:region_num
        if link_width_ratio_sum(k) ~= 0
            link_width_percent(k,:) = floor(link_width_ratio(k,:)/link_width_ratio_sum(k)*spacing);
            % use 'floor', sum is less than spacing(100), compensate the rest to first order
            if sum(link_width_percent(k,:)) ~= spacing
                link_width_percent(k,1) = link_width_percent(k,1) + (spacing - sum(link_width_percent(k,:)));
            end
        end
    end
    link_order_ploted = zeros(region_num,1); % initialize, record ploted order
    for k = 1:cor_num
        % calculate chromsome1's network and region
        cor_row_net = RawDataCircos.HigherOrderNetworkIndex(cor_row(k));
        cor_row_reg = cor_row(k);
        if cor_row_net > 1
            for l = 1:cor_row_net-1
                cor_row_reg = cor_row_reg - network_region_num(l);
            end
        end
        % calculate chromsome2's network and region
        cor_col_net = RawDataCircos.HigherOrderNetworkIndex(cor_col(k));
        cor_col_reg = cor_col(k);
        if cor_col_net > 1
            for l = 1:cor_col_net-1
                cor_col_reg = cor_col_reg - network_region_num(l);
            end
        end
        % calculate chromsome1's start and end
        if link_order_ploted(cor_row(k)) == 0
            cor_row_start = (cor_row_reg-1)*spacing;
        else
            cor_row_start = (cor_row_reg-1)*spacing + sum(link_width_percent(cor_row(k),1:link_order_ploted(cor_row(k))));
        end
        link_order_ploted(cor_row(k)) = link_order_ploted(cor_row(k)) + 1;
        cor_row_end = cor_row_start + link_width_percent(cor_row(k),link_order_ploted(cor_row(k))) - 1;
        % calculate chromsome2's start and end
        if link_order_ploted(cor_col(k)) == 0
            cor_col_start = (cor_col_reg-1)*spacing;
        else
            cor_col_start = (cor_col_reg-1)*spacing + sum(link_width_percent(cor_col(k),1:link_order_ploted(cor_col(k))));
        end
        link_order_ploted(cor_col(k)) = link_order_ploted(cor_col(k)) + 1;
        cor_col_end = cor_col_start + link_width_percent(cor_col(k),link_order_ploted(cor_col(k))) - 1;
        % print on txt according to format
        fprintf(fid,'net%u %u %u ',cor_row_net,cor_row_start,cor_row_end);
        fprintf(fid,'net%u %u %u ',cor_col_net,cor_col_start,cor_col_end);
        fprintf(fid,'color=%u,%u,%u,%.1f',color_mat(k,1),color_mat(k,2),color_mat(k,3),link_transparency);
        fprintf(fid,'\n');
    end
elseif link_mode==3
    link_order_ploted = zeros(region_num,1); % initialize, record ploted order
    % calculate start point of regions
    reg_start = zeros(3,max(network_region_num)); % initialize
    index = 1;
    for k = 1:network_num
        former_sum = 0;
        for l = 1:network_region_num(k)
            reg_start(k,l) = former_sum;
            former_sum = former_sum + band_width(index);
            index = index + 1;
        end
    end
    for k = 1:cor_num
        % calculate chromsome1's network and region
        cor_row_net = RawDataCircos.HigherOrderNetworkIndex(cor_row(k));
        cor_row_reg = cor_row(k);
        if cor_row_net > 1
            for l = 1:cor_row_net-1
                cor_row_reg = cor_row_reg - network_region_num(l);
            end
        end
        % calculate chromsome2's network and region
        cor_col_net = RawDataCircos.HigherOrderNetworkIndex(cor_col(k));
        cor_col_reg = cor_col(k);
        if cor_col_net > 1
            for l = 1:cor_col_net-1
                cor_col_reg = cor_col_reg - network_region_num(l);
            end
        end
        % calculate chromsome1's start and end
        if link_order_ploted(cor_row(k)) == 0
            cor_row_start = reg_start(cor_row_net,cor_row_reg);
        else
            cor_row_start = reg_start(cor_row_net,cor_row_reg) + sum(link_width_ratio100(cor_row(k),1:link_order_ploted(cor_row(k))));
        end
        link_order_ploted(cor_row(k)) = link_order_ploted(cor_row(k)) + 1;
        cor_row_end = cor_row_start + link_width_ratio100(cor_row(k),link_order_ploted(cor_row(k))) - 1;
        % calculate chromsome2's start and end
        if link_order_ploted(cor_col(k)) == 0
            cor_col_start = reg_start(cor_col_net,cor_col_reg);
        else
            cor_col_start = reg_start(cor_col_net,cor_col_reg) + sum(link_width_ratio100(cor_col(k),1:link_order_ploted(cor_col(k))));
        end
        link_order_ploted(cor_col(k)) = link_order_ploted(cor_col(k)) + 1;
        cor_col_end = cor_col_start + link_width_ratio100(cor_col(k),link_order_ploted(cor_col(k))) - 1;
        % print on txt according to format
        fprintf(fid,'net%u %u %u ',cor_row_net,cor_row_start,cor_row_end);
        fprintf(fid,'net%u %u %u ',cor_col_net,cor_col_start,cor_col_end);
        fprintf(fid,'color=%u,%u,%u,%.1f',color_mat(k,1),color_mat(k,2),color_mat(k,3),link_transparency);
        fprintf(fid,'\n');
    end
end
fclose(fid);


