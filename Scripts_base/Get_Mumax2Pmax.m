function Data_ConstructMP_Mumax2Pmax = Get_Mumax2Pmax(Data_struct_construct,slack_time_mumax,slack_time_pmax)
%
% UPDATED 26/12/2024 to get also the syntesis rate Pi_p and the growth rate Mu
%
% This function gets the MEFL w.r.t. Particles data from the time instant
% at which mu_max is reached until the one in which the maximum number of
% Particles is reached. The inputs lack_time_mumax,slack_time_pmax allow
% for a time slack. Introduce them as positive to delay and negative for
% advance
%
arguments
    Data_struct_construct struct
    slack_time_mumax double
    slack_time_pmax double
end

Data_ConstructMP_Mumax2Pmax={};
 shortest_index = 100000;
 num_instances = length(Data_struct_construct.Stats_Mu.List_instances);
 for i=1:num_instances
     t_mumax_mean = Data_struct_construct.Stats_Mu.List_instances{i}.Index_Time_Mumax_mean+slack_time_mumax;
     t_particlesmax_mean = Data_struct_construct.Stats_Particles.List_instances{i}.Index_Time_Particlesmax_mean+slack_time_pmax;
     Data_ConstructMP_Mumax2Pmax.List_Instances{i}.MEFL_mean = Data_struct_construct.Stats_MEFL.List_instances{i}.Data_mean(t_mumax_mean:t_particlesmax_mean);
     Data_ConstructMP_Mumax2Pmax.List_Instances{i}.Particles_mean = Data_struct_construct.Stats_Particles.List_instances{i}.Data_mean(t_mumax_mean:t_particlesmax_mean);
     Data_ConstructMP_Mumax2Pmax.List_Instances{i}.Pi_p_mean = Data_struct_construct.Stats_Pi_p.List_instances{i}.Data_mean(t_mumax_mean:t_particlesmax_mean);
     Data_ConstructMP_Mumax2Pmax.List_Instances{i}.Mu_mean = Data_struct_construct.Stats_Mu.List_instances{i}.Data_mean(t_mumax_mean:t_particlesmax_mean);
     Data_ConstructMP_Mumax2Pmax.List_Instances{i}.MEFL_std = Data_struct_construct.Stats_MEFL.List_instances{i}.Data_std(t_mumax_mean:t_particlesmax_mean);
     Data_ConstructMP_Mumax2Pmax.List_Instances{i}.Particles_std = Data_struct_construct.Stats_Particles.List_instances{i}.Data_std(t_mumax_mean:t_particlesmax_mean);
     Data_ConstructMP_Mumax2Pmax.List_Instances{i}.Pi_p_std = Data_struct_construct.Stats_Pi_p.List_instances{i}.Data_std(t_mumax_mean:t_particlesmax_mean);
     Data_ConstructMP_Mumax2Pmax.List_Instances{i}.Mu_std = Data_struct_construct.Stats_Mu.List_instances{i}.Data_std(t_mumax_mean:t_particlesmax_mean);
     Data_ConstructMP_Mumax2Pmax.List_Instances{i}.Exp_IDnum = Data_struct_construct.Stats_Particles.List_instances{i}.Experiment_IDnum;
     Data_ConstructMP_Mumax2Pmax.List_Instances{i}.Exp_date = Data_struct_construct.Stats_Particles.List_instances{i}.Experiment_date;
     Data_ConstructMP_Mumax2Pmax.List_Instances{i}.Mu_max = Data_struct_construct.Stats_Mu.List_instances{i}.Value_Mumax_mean;

     no_wells = length(Data_struct_construct.Stats_Mu.List_instances{i}.Indices_Times_Mumax_wells);
     for j=1:no_wells
          t_mumax_well = Data_struct_construct.Stats_Mu.List_instances{i}.Indices_Times_Mumax_wells(j)+slack_time_mumax;
          t_particlesmax_well = Data_struct_construct.Stats_Particles.List_instances{i}.Index_Time_Particlesmax_wells(j)+slack_time_pmax;
          Data_ConstructMP_Mumax2Pmax.List_Instances{i}.Data_well{j}.MEFL = Data_struct_construct.Stats_MEFL.List_instances{i}.Data_raw(t_mumax_well:t_particlesmax_well,j);
          Data_ConstructMP_Mumax2Pmax.List_Instances{i}.Data_well{j}.Particles = Data_struct_construct.Stats_Particles.List_instances{i}.Data_raw(t_mumax_well:t_particlesmax_well,j);
          Data_ConstructMP_Mumax2Pmax.List_Instances{i}.Data_well{j}.Pi_p = Data_struct_construct.Stats_Pi_p.List_instances{i}.Data_raw(t_mumax_well:t_particlesmax_well,j);
          Data_ConstructMP_Mumax2Pmax.List_Instances{i}.Data_well{j}.Mu = Data_struct_construct.Stats_Mu.List_instances{i}.Data_raw(t_mumax_well:t_particlesmax_well,j);
           if size(Data_ConstructMP_Mumax2Pmax.List_Instances{i}.Data_well{j}.Particles,1) < shortest_index 
               shortest_index = size(Data_ConstructMP_Mumax2Pmax.List_Instances{i}.Data_well{j}.Particles,1);
           end
     end% for no_wells
 end % for num_instances_BK

 Data_ConstructMP_Mumax2Pmax.Shortest_index = shortest_index;
 t_mumax_global = Data_struct_construct.Stats_Mu.Global_index_Time_Mumax_mean+slack_time_mumax;
 t_particlesmax_global =  Data_struct_construct.Stats_Particles.Global_index_Time_Particlesmax_mean+slack_time_pmax;
 Data_ConstructMP_Mumax2Pmax.MEFL_global_mean =  Data_struct_construct.Stats_MEFL.Global_data_mean(t_mumax_global:t_particlesmax_global);
 Data_ConstructMP_Mumax2Pmax.Particles_global_mean =  Data_struct_construct.Stats_Particles.Global_data_mean(t_mumax_global:t_particlesmax_global);
 Data_ConstructMP_Mumax2Pmax.Pi_p_global_mean =  Data_struct_construct.Stats_Pi_p.Global_data_mean(t_mumax_global:t_particlesmax_global);
 Data_ConstructMP_Mumax2Pmax.Mu_global_mean =  Data_struct_construct.Stats_Mu.Global_data_mean(t_mumax_global:t_particlesmax_global);
 Data_ConstructMP_Mumax2Pmax.MEFL_global_std =  Data_struct_construct.Stats_MEFL.Global_data_std(t_mumax_global:t_particlesmax_global);
 Data_ConstructMP_Mumax2Pmax.Particles_global_std =  Data_struct_construct.Stats_Particles.Global_data_std(t_mumax_global:t_particlesmax_global);
 Data_ConstructMP_Mumax2Pmax.Pi_p_global_std =  Data_struct_construct.Stats_Pi_p.Global_data_std(t_mumax_global:t_particlesmax_global);
 Data_ConstructMP_Mumax2Pmax.Mu_global_std =  Data_struct_construct.Stats_Mu.Global_data_std(t_mumax_global:t_particlesmax_global);
 end % function