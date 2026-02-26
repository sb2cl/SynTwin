function Plot_Results_Lib6_reduced_model(Results_Lib,Lower_Bounds,Upper_Bounds,indices_plasmids,indices_promoters,indices_rbss, title_text)

color_blue = hex2rgb('#00008B');
color_grey = hex2rgb('#4b4b4b');
color_grey_boira = hex2rgb('#F2F2F2');
color_grey_neutre = hex2rgb('#CCCCCC');
RBS_colors= ['#A8DADC';'#F4A261';'#B5E48C';'#CDB4DB';'#FFE066'];
Promoter_colors = ['#264653';'#E76F51';'#2A9D8F';'#6A0572'];

num_plasmids = length(indices_plasmids);
num_promoters=length(indices_promoters);
num_rbss = length(indices_rbss);

x_L = Lower_Bounds;
x_U = Upper_Bounds;

%% ESTIMATED PARAMETERS HISTOGRAMS



% HISTOGRAMS OF RBSs k0/sigma0
figure();
tiledlayout(1,num_rbss);
% First the k/sigma's
for k=1:num_rbss
    nexttile;
    Data_rbs_k0_sigma0_mean = Results_Lib{1,1,k}.RBS_k0_sigma0_mean;
    Data_rbs_k0_sigma0_raw = Results_Lib{1,1,k}.RBS_k0_sigma0_raw;
    pd1 = histfit(Data_rbs_k0_sigma0_raw,[],'normal');%nbins set to default round(sqrt(num_data))
    %pd1 = histfit(Data_Data_rbs_k_mean_raw,[],'kernel','Kernel','epanechnikov');%nbins set to default round(sqrt(num_data))
    grid on;
    pd1(1).FaceColor = RBS_colors(k,:);
    pd1(2).Color = RBS_colors(k,:);
    pd1(1).FaceAlpha = 1;
    hold on
    xline( Data_rbs_k0_sigma0_mean,'Color','k','LineWidth',1 ); 
    xlabel('$\frac{k_0}{\sigma_0}\, (\mathrm{molec}^{-1})$','FontSize',20,'Interpreter','latex')
    if k==1
      t= title([title_text, Results_Lib{1,1,k}.TU_RBS]);
    else
       t= title([Results_Lib{1,1,k}.TU_RBS]);
    end
    t.FontSize = 18;
    ax=gca;
    ax.XAxis.FontSize = 18;
    ax.YAxis.FontSize = 18;
    ax.XLabel.FontSize = 18;
    ax.YLabel.FontSize = 18;
end


% END OF HISTOGRAMS OF ESTIMATED PARAMETERS

%% PREDICTED SYNTHESIS RATES VERSUS EXPERIMENTAL ONES

% TILED PLOT OF PI_PRED, PI_EXP VERSUS MU (95% CONFIDENCE INTERVAL)

figure();
tiledlayout(num_plasmids,num_promoters);
for i=1:num_plasmids %Plasmids
   for j=1:num_promoters %Promotorers
      for k=1:num_rbss  %RBS
          if not(isempty(Results_Lib{i,j,k}))
            tile_num = j + (num_plasmids+1)*(i-1); 
            nexttile(tile_num)
            hold on
            Bioparts_TU= Results_Lib{i,j,k}.TU_Bioparts;
            Color_TU= Results_Lib{i,j,k}.TU_color_code;
            num_slices = length(Results_Lib{i,j,k}.MC_mu_slices);
            for q=1:num_slices
              xmu(q) = Results_Lib{i,j,k}.MC_mu_slices(q).Mu_slice;
              yl(q) = Results_Lib{i,j,k}.MC_mu_slices(q).Pi_pred_q2p5;
              yu(q) = Results_Lib{i,j,k}.MC_mu_slices(q).Pi_pred_q97p5;
              ymean(q) = Results_Lib{i,j,k}.MC_mu_slices(q).Pi_pred_q50;
            end
            plot(xmu,ymean,'linestyle','-','Color',color_blue,'LineWidth',2,'HandleVisibility','off')
            plot(xmu,yl,'linestyle','-','Color',color_blue,'LineWidth',1,'HandleVisibility','off')
            plot(xmu, yu,'linestyle','-','Color',color_blue,'LineWidth',1,'HandleVisibility','off')
            p = fill([xmu,xmu(end:-1:1)],[yl,yu(end:-1:1)],Color_TU,'FaceAlpha',0.2);
            % Full set of experimental data 
            num_instances = length(Results_Lib{i,j,k}.Instances);
            for q=1:num_instances
              % For the average of the wells in each instance of each experiment use this below
               y = Results_Lib{i,j,k}.Instances{q}.Mu_mumax_pmax_instance_mean;
               x = Results_Lib{i,j,k}.Instances{q}.Pi_mumax_pmax_instance_mean;
               s= plot(x,y,'Color',color_grey_neutre,'LineStyle','-','LineWidth',0.5,'HandleVisibility','off');
            end %instances
            grid on
            ax=gca;
            ax.XLim=[xmu(1),xmu(end)];
            ax.XAxis.FontSize = 18;
            ax.YAxis.FontSize = 18;
            if i==1 && j==1 && k==1
               t= title([title_text, Bioparts_TU]);
            else
               t= title(Bioparts_TU);
            end
            t.FontSize = 18;
            xlabel('$\mu$','FontSize',18,'Interpreter','latex')
            ylabel('$\Pi_p$','FontSize',18,'Interpreter','latex')
          end %not empty TU
      end
   end
end

% TILED PLOT OF PI_PRED, PI_EXP VERSUS MU (95% CONFIDENCE INTERVAL IN BOTH)
figure();
tiledlayout(num_plasmids,num_promoters);
for i=1:num_plasmids %Plasmids
   for j=1:num_promoters %Promotorers
      for k=1:num_rbss  %RBS
          if not(isempty(Results_Lib{i,j,k}))
            tile_num = j +  (num_plasmids+1)*(i-1); 
            nexttile(tile_num)
            hold on
            Bioparts_TU= Results_Lib{i,j,k}.TU_Bioparts;
            Color_TU= Results_Lib{i,j,k}.TU_color_code;
            num_slices = length(Results_Lib{i,j,k}.MC_mu_slices);
            for q=1:num_slices
              xmu(q) = Results_Lib{i,j,k}.MC_mu_slices(q).Mu_slice;
              yl(q) = Results_Lib{i,j,k}.MC_mu_slices(q).Pi_pred_q2p5;
              yu(q) = Results_Lib{i,j,k}.MC_mu_slices(q).Pi_pred_q97p5;
              ymean(q) = Results_Lib{i,j,k}.MC_mu_slices(q).Pi_pred_q50;
            end
            plot(xmu,ymean,'linestyle','-','Color',Color_TU,'LineWidth',2,'HandleVisibility','off')
            plot(xmu,yl,'linestyle','-','Color',Color_TU,'LineWidth',1,'HandleVisibility','off')
            plot(xmu, yu,'linestyle','-','Color',Color_TU,'LineWidth',1,'HandleVisibility','off')
            p = fill([xmu,xmu(end:-1:1)],[yl,yu(end:-1:1)],Color_TU,'FaceAlpha',0.99);

             y_g = Results_Lib{i,j,k}.Pi_mumax_pmax_global_mean;
             x_g = Results_Lib{i,j,k}.Mu_mumax_pmax_global_mean;
             y_g_std = Results_Lib{i,j,k}.Pi_mumax_pmax_global_std;
             x_gn = x_g(end:-1:1);
             y_gn = y_g(end:-1:1);
             y_gn_std = y_g_std(end:-1:1);
             yln = y_gn - 2*y_gn_std;
             yun = y_gn + 2*y_gn_std;
             plot(x_gn,yln,'linestyle','-','Color',color_grey,'LineWidth',0.5,'HandleVisibility','off')
             plot(x_gn, yun,'linestyle','-','Color',color_grey,'LineWidth',0.5,'HandleVisibility','off')
             p = fill([x_gn;x_gn(end:-1:1)],[yln;yun(end:-1:1)],color_grey,'FaceAlpha',0.25);
             plot(x_g,y_g,'Color',color_grey_neutre,'LineWidth',2,'HandleVisibility','on')
             grid on
             ax=gca;
             ax.XLim=[x_gn(1),x_gn(end)];
             ax.YLim=[max(0.8*yln(1),1e-4),1.2*yun(end)];
            ax.XAxis.FontSize = 18;
            ax.YAxis.FontSize = 18;
             if i==1 && j==1 && k==1
               t= title([title_text, Bioparts_TU]);
             else
               t= title(Bioparts_TU);
             end
             t.FontSize = 18;
             if j==1
                  ylabel('$\Pi_p$','FontSize',20,'Interpreter','latex')
             end
             if k==num_rbss
                  xlabel('$\mu$','FontSize',20,'Interpreter','latex')
             end
            %  yscale("log")
          end %not empty
      end %rbs
   end %promoter
end %plasmid
