
function [Mu_exp_nut, Nu_t_exp, Mu_exp_PhiR, PhiR_exp, color_Cremer_phiR] = Load_Data_Nut_PhiR()

opts = spreadsheetImportOptions("NumVariables", 2);
opts.VariableTypes = ["double", "double"];
opts.Sheet = "Data_nut";
%opts.DataRange = 'A3:B73';
opts.DataRange = 'A3:B67'; %without Bremer's last data
path_Cremer_data = 'ecoli_peptide_elongation_rates_updated';
Cremer_nut = readtable(path_Cremer_data,opts, "UseExcel", false);
data_Cremer_nut = table2array(Cremer_nut);
Mu_exp_nut = data_Cremer_nut(:,1);
Nu_t_exp = data_Cremer_nut(:,2);

opts = spreadsheetImportOptions("NumVariables", 6);
opts.VariableTypes = ["double", "double","char","char","char","double"];
opts.Sheet = "Data_phiR";
%opts.DataRange = 'A2:F283';
opts.DataRange = 'A2:F272'; %without Bremer's last data
path_Cremer_data = 'ecoli_ribosomal_mass_fractions_updated';
Cremer_phiR = readtable(path_Cremer_data,opts, "UseExcel", false);
Mu_exp_PhiR = table2array(Cremer_phiR(:,1));
PhiR_exp = table2array(Cremer_phiR(:,2));
color_Cremer_phiR = table2array(Cremer_phiR(:,6));

end