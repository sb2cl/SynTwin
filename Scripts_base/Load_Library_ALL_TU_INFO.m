%This script Load_Library_ALL_TU_INFO.m 
% creates the cell data type  Cell_TU_info_ALL with all the LIBRARY general information concerning the names of Oris,
% promoters and RBSs and the names of the combinatorial TUs. 
% 
% These names are used to get the experimental data: 
%
% ORI  pGreen (high copy):
%                       J23106                  J23102              J23101                                                   J23100                          J23116                        J23117
% B0030    'pIAKA1_11'=1       'pIAKA1_41'=6      'pIAKA1_51'=11                            'pSMKA1_81'=31
% B0032    'pIAKA1_12'=2       'pIAKA1_42'=7      'pIAKA1_52'=12                            'pSMKA1_82'=32
% B0034    'pIAKA1_13'=3       'pIAKA1_43'=8      'pIAKA1_53'=13                            'pSMKA1_83'=33           'pSMKA1_93'=38        'pSMKA1_103'=43
% J61100   'pIAKA1_14'=4       'pIAKA1_44'=9      'pIAKA1_54'=14                            'pSMKA1_84'=34
% J61101   'pIAKA1_15'=5       'pIAKA1_45'=10    'pIAKA1_55'=15                            'pSMKA1_85'=35
% ORI  pSC101 (low copy):
%                       J23106                            J23102                                J23101                          J23100                             J23116                                            J23117
%  B0030    'pIAKA1S271_11'=16    'pIAKA1S271_41'=21  'pIAKA1S271_51'=26
%  B0032    'pIAKA1S271_12'=17    'pIAKA1S271_42'=22  'pIAKA1S271_52'=27
%  B0034    'pIAKA1S271_13'=18    'pIAKA1S271_43'=23  'pIAKA1S271_53'=28                        
%  J61100   'pIAKA1S271_14'=19    'pIAKA1S271_44'=24  'pIAKA1S271_54'=29
%  J61101   'pIAKA1S271_15'=20    'pIAKA1S271_45'=25  'pIAKA1S271_55'=30

% RBS COLOR CODES:
RBS_colors=['A8DADC';'F4A261';'B5E48C';'CDB4DB';'FFE066'];
Promoter_colors = ['#264653';'#E76F51';'#2A9D8F';'#6A0572'];

% IMPORTANT: the cell data type Cell_TU_info must follow the
% ordering logical structure {Ori, promoter, RBS}, and contain:
%   Cell_TU_info{idx_o,idx_p,dx_r}.name_plasmid
%   Cell_TU_info{idx_o,idx_p,dx_r}.name_promoter
%   Cell_TU_info{idx_o,idx_p,dx_r}.name_rbs
%   Cell_TU_info{idx_o,idx_p,dx_r}.color_rbs hexadecimal code for RBS color
%   Cell_TU_info{idx_o,idx_p,dx_r}.name_TU

%% HIGH COPY:

Cell_TU_info_ALL{1,1,1}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,1,1}.name_promoter = 'J23106';
Cell_TU_info_ALL{1,1,1}.name_rbs = 'B0030';
Cell_TU_info_ALL{1,1,1}.color_rbs = RBS_colors(1,:);
Cell_TU_info_ALL{1,1,1}.name_TU = 'pIAKA1_11';

Cell_TU_info_ALL{1,1,2}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,1,2}.name_promoter = 'J23106';
Cell_TU_info_ALL{1,1,2}.name_rbs = 'B0032';
Cell_TU_info_ALL{1,1,2}.color_rbs = RBS_colors(2,:);
Cell_TU_info_ALL{1,1,2}.name_TU = 'pIAKA1_12';

Cell_TU_info_ALL{1,1,3}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,1,3}.name_promoter = 'J23106';
Cell_TU_info_ALL{1,1,3}.name_rbs = 'B0034';
Cell_TU_info_ALL{1,1,3}.color_rbs = RBS_colors(3,:);
Cell_TU_info_ALL{1,1,3}.name_TU = 'pIAKA1_13';

Cell_TU_info_ALL{1,1,4}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,1,4}.name_promoter = 'J23106';
Cell_TU_info_ALL{1,1,4}.name_rbs = 'J61100';
Cell_TU_info_ALL{1,1,4}.color_rbs = RBS_colors(4,:);
Cell_TU_info_ALL{1,1,4}.name_TU = 'pIAKA1_14';

Cell_TU_info_ALL{1,1,5}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,1,5}.name_promoter = 'J23106';
Cell_TU_info_ALL{1,1,5}.name_rbs = 'J61101';
Cell_TU_info_ALL{1,1,5}.color_rbs = RBS_colors(5,:);
Cell_TU_info_ALL{1,1,5}.name_TU = 'pIAKA1_15';

%%
Cell_TU_info_ALL{1,2,1}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,2,1}.name_promoter = 'J23102';
Cell_TU_info_ALL{1,2,1}.name_rbs = 'B0030';
Cell_TU_info_ALL{1,2,1}.color_rbs = RBS_colors(1,:);
Cell_TU_info_ALL{1,2,1}.name_TU = 'pIAKA1_41';

Cell_TU_info_ALL{1,2,2}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,2,2}.name_promoter = 'J23102';
Cell_TU_info_ALL{1,2,2}.name_rbs = 'B0032';
Cell_TU_info_ALL{1,2,2}.color_rbs = RBS_colors(2,:);
Cell_TU_info_ALL{1,2,2}.name_TU = 'pIAKA1_42';

Cell_TU_info_ALL{1,2,3}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,2,3}.name_promoter = 'J23102';
Cell_TU_info_ALL{1,2,3}.name_rbs = 'B0034';
Cell_TU_info_ALL{1,2,3}.color_rbs = RBS_colors(3,:);
Cell_TU_info_ALL{1,2,3}.name_TU = 'pIAKA1_43';

Cell_TU_info_ALL{1,2,4}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,2,4}.name_promoter = 'J23102';
Cell_TU_info_ALL{1,2,4}.name_rbs = 'J61100';
Cell_TU_info_ALL{1,2,4}.color_rbs = RBS_colors(4,:);
Cell_TU_info_ALL{1,2,4}.name_TU = 'pIAKA1_44';

Cell_TU_info_ALL{1,2,5}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,2,5}.name_promoter = 'J23102';
Cell_TU_info_ALL{1,2,5}.name_rbs = 'J61101';
Cell_TU_info_ALL{1,2,5}.color_rbs = RBS_colors(5,:);
Cell_TU_info_ALL{1,2,5}.name_TU = 'pIAKA1_45';

%%

Cell_TU_info_ALL{1,3,1}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,3,1}.name_promoter = 'J23101';
Cell_TU_info_ALL{1,3,1}.name_rbs = 'B0030';
Cell_TU_info_ALL{1,3,1}.color_rbs = RBS_colors(1,:);
Cell_TU_info_ALL{1,3,1}.name_TU = 'pIAKA1_51';

Cell_TU_info_ALL{1,3,2}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,3,2}.name_promoter = 'J23101';
Cell_TU_info_ALL{1,3,2}.name_rbs = 'B0032';
Cell_TU_info_ALL{1,3,2}.color_rbs = RBS_colors(2,:);
Cell_TU_info_ALL{1,3,2}.name_TU = 'pIAKA1_52';

Cell_TU_info_ALL{1,3,3}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,3,3}.name_promoter = 'J23101';
Cell_TU_info_ALL{1,3,3}.name_rbs = 'B0034';
Cell_TU_info_ALL{1,3,3}.color_rbs = RBS_colors(3,:);
Cell_TU_info_ALL{1,3,3}.name_TU = 'pIAKA1_53';

Cell_TU_info_ALL{1,3,4}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,3,4}.name_promoter = 'J23101';
Cell_TU_info_ALL{1,3,4}.name_rbs = 'J61100';
Cell_TU_info_ALL{1,3,4}.color_rbs = RBS_colors(4,:);
Cell_TU_info_ALL{1,3,4}.name_TU = 'pIAKA1_54';

Cell_TU_info_ALL{1,3,5}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,3,5}.name_promoter = 'J23101';
Cell_TU_info_ALL{1,3,5}.name_rbs = 'J61101';
Cell_TU_info_ALL{1,3,5}.color_rbs = RBS_colors(5,:);
Cell_TU_info_ALL{1,3,5}.name_TU = 'pIAKA1_55';

%%

Cell_TU_info_ALL{1,4,1}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,4,1}.name_promoter = 'J23100';
Cell_TU_info_ALL{1,4,1}.name_rbs = 'B0030';
Cell_TU_info_ALL{1,4,1}.color_rbs = RBS_colors(1,:);
Cell_TU_info_ALL{1,4,1}.name_TU = 'pSMKA1_81';

Cell_TU_info_ALL{1,4,2}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,4,2}.name_promoter = 'J23100';
Cell_TU_info_ALL{1,4,2}.name_rbs = 'B0032';
Cell_TU_info_ALL{1,4,2}.color_rbs = RBS_colors(2,:);
Cell_TU_info_ALL{1,4,2}.name_TU = 'pSMKA1_82';

Cell_TU_info_ALL{1,4,3}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,4,3}.name_promoter = 'J23100';
Cell_TU_info_ALL{1,4,3}.name_rbs = 'B0034';
Cell_TU_info_ALL{1,4,3}.color_rbs = RBS_colors(3,:);
Cell_TU_info_ALL{1,4,3}.name_TU = 'pSMKA1_83';

Cell_TU_info_ALL{1,4,4}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,4,4}.name_promoter = 'J23100';
Cell_TU_info_ALL{1,4,4}.name_rbs = 'J61100';
Cell_TU_info_ALL{1,4,4}.color_rbs = RBS_colors(4,:);
Cell_TU_info_ALL{1,4,4}.name_TU = 'pSMKA1_84';

Cell_TU_info_ALL{1,4,5}.name_plasmid =  'pGreen';
Cell_TU_info_ALL{1,4,5}.name_promoter = 'J23100';
Cell_TU_info_ALL{1,4,5}.name_rbs = 'J61101';
Cell_TU_info_ALL{1,4,5}.color_rbs = RBS_colors(5,:);
Cell_TU_info_ALL{1,4,5}.name_TU = 'pSMKA1_85';

%%

Cell_TU_info_ALL{1,5,3}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,5,3}.name_promoter = 'J23116';
Cell_TU_info_ALL{1,5,3}.name_rbs = 'B0034';
Cell_TU_info_ALL{1,5,3}.color_rbs = RBS_colors(3,:);
Cell_TU_info_ALL{1,5,3}.name_TU = 'pSMKA1_93';

%%

Cell_TU_info_ALL{1,6,3}.name_plasmid = 'pGreen';
Cell_TU_info_ALL{1,6,3}.name_promoter = 'J23117';
Cell_TU_info_ALL{1,6,3}.name_rbs = 'B0034';
Cell_TU_info_ALL{1,6,3}.color_rbs = RBS_colors(3,:);
Cell_TU_info_ALL{1,6,3}.name_TU = 'pSMKA1_103';

%% LOW COPY:

Cell_TU_info_ALL{2,1,1}.name_plasmid = 'pSC101';
Cell_TU_info_ALL{2,1,1}.name_promoter = 'J23106';
Cell_TU_info_ALL{2,1,1}.name_rbs = 'B0030';
Cell_TU_info_ALL{2,1,1}.color_rbs = RBS_colors(1,:);
Cell_TU_info_ALL{2,1,1}.name_TU = 'pIAKA1S271_11';

Cell_TU_info_ALL{2,1,2}.name_plasmid = 'pSC101';
Cell_TU_info_ALL{2,1,2}.name_promoter = 'J23106';
Cell_TU_info_ALL{2,1,2}.name_rbs = 'B0032';
Cell_TU_info_ALL{2,1,2}.color_rbs = RBS_colors(2,:);
Cell_TU_info_ALL{2,1,2}.name_TU = 'pIAKA1S271_12';

Cell_TU_info_ALL{2,1,3}.name_plasmid = 'pSC101';
Cell_TU_info_ALL{2,1,3}.name_promoter = 'J23106';
Cell_TU_info_ALL{2,1,3}.name_rbs = 'B0034';
Cell_TU_info_ALL{2,1,3}.color_rbs = RBS_colors(3,:);
Cell_TU_info_ALL{2,1,3}.name_TU = 'pIAKA1S271_13';

Cell_TU_info_ALL{2,1,4}.name_plasmid = 'pSC101';
Cell_TU_info_ALL{2,1,4}.name_promoter = 'J23106';
Cell_TU_info_ALL{2,1,4}.name_rbs = 'J61100';
Cell_TU_info_ALL{2,1,4}.color_rbs = RBS_colors(4,:);
Cell_TU_info_ALL{2,1,4}.name_TU = 'pIAKA1S271_14';

Cell_TU_info_ALL{2,1,5}.name_plasmid = 'pSC101';
Cell_TU_info_ALL{2,1,5}.name_promoter = 'J23106';
Cell_TU_info_ALL{2,1,5}.name_rbs = 'J61101';
Cell_TU_info_ALL{2,1,5}.color_rbs = RBS_colors(5,:);
Cell_TU_info_ALL{2,1,5}.name_TU = 'pIAKA1S271_15';

%%
Cell_TU_info_ALL{2,2,1}.name_plasmid = 'pSC101';
Cell_TU_info_ALL{2,2,1}.name_promoter = 'J23102';
Cell_TU_info_ALL{2,2,1}.name_rbs = 'B0030';
Cell_TU_info_ALL{2,2,1}.color_rbs = RBS_colors(1,:);
Cell_TU_info_ALL{2,2,1}.name_TU = 'pIAKA1S271_41';

Cell_TU_info_ALL{2,2,2}.name_plasmid = 'pSC101';
Cell_TU_info_ALL{2,2,2}.name_promoter = 'J23102';
Cell_TU_info_ALL{2,2,2}.name_rbs = 'B0032';
Cell_TU_info_ALL{2,2,2}.color_rbs = RBS_colors(2,:);
Cell_TU_info_ALL{2,2,2}.name_TU = 'pIAKA1S271_42';

Cell_TU_info_ALL{2,2,3}.name_plasmid = 'pSC101';
Cell_TU_info_ALL{2,2,3}.name_promoter = 'J23102';
Cell_TU_info_ALL{2,2,3}.name_rbs = 'B0034';
Cell_TU_info_ALL{2,2,3}.color_rbs = RBS_colors(3,:);
Cell_TU_info_ALL{2,2,3}.name_TU = 'pIAKA1S271_43';

Cell_TU_info_ALL{2,2,4}.name_plasmid = 'pSC101';
Cell_TU_info_ALL{2,2,4}.name_promoter = 'J23102';
Cell_TU_info_ALL{2,2,4}.name_rbs = 'J61100';
Cell_TU_info_ALL{2,2,4}.color_rbs = RBS_colors(4,:);
Cell_TU_info_ALL{2,2,4}.name_TU = 'pIAKA1S271_44';

Cell_TU_info_ALL{2,2,5}.name_plasmid = 'pSC101';
Cell_TU_info_ALL{2,2,5}.name_promoter = 'J23102';
Cell_TU_info_ALL{2,2,5}.name_rbs = 'J61101';
Cell_TU_info_ALL{2,2,5}.color_rbs = RBS_colors(5,:);
Cell_TU_info_ALL{2,2,5}.name_TU = 'pIAKA1S271_45';

%%

Cell_TU_info_ALL{2,3,1}.name_plasmid = 'pSC101';
Cell_TU_info_ALL{2,3,1}.name_promoter = 'J23101';
Cell_TU_info_ALL{2,3,1}.name_rbs = 'B0030';
Cell_TU_info_ALL{2,3,1}.color_rbs = RBS_colors(1,:);
Cell_TU_info_ALL{2,3,1}.name_TU = 'pIAKA1S271_51';

Cell_TU_info_ALL{2,3,2}.name_plasmid = 'pSC101';
Cell_TU_info_ALL{2,3,2}.name_promoter = 'J23101';
Cell_TU_info_ALL{2,3,2}.name_rbs = 'B0032';
Cell_TU_info_ALL{2,3,2}.color_rbs = RBS_colors(2,:);
Cell_TU_info_ALL{2,3,2}.name_TU = 'pIAKA1S271_52';

Cell_TU_info_ALL{2,3,3}.name_plasmid = 'pSC101';
Cell_TU_info_ALL{2,3,3}.name_promoter = 'J23101';
Cell_TU_info_ALL{2,3,3}.name_rbs = 'B0034';
Cell_TU_info_ALL{2,3,3}.color_rbs = RBS_colors(3,:);
Cell_TU_info_ALL{2,3,3}.name_TU = 'pIAKA1S271_53';

Cell_TU_info_ALL{2,3,4}.name_plasmid = 'pSC101';
Cell_TU_info_ALL{2,3,4}.name_promoter = 'J23101';
Cell_TU_info_ALL{2,3,4}.name_rbs = 'J61100';
Cell_TU_info_ALL{2,3,4}.color_rbs = RBS_colors(4,:);
Cell_TU_info_ALL{2,3,4}.name_TU = 'pIAKA1S271_54';

Cell_TU_info_ALL{2,3,5}.name_plasmid = 'pSC101';
Cell_TU_info_ALL{2,3,5}.name_promoter = 'J23101';
Cell_TU_info_ALL{2,3,5}.name_rbs = 'J61101';
Cell_TU_info_ALL{2,3,5}.color_rbs = RBS_colors(5,:);
Cell_TU_info_ALL{2,3,5}.name_TU = 'pIAKA1S271_55';





