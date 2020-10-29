clear; close all; clc
load('Pa.mat')
pa = par;
clear('par');
for exp_Fip = 1:1
    for i_VPLC = 1:1
        for i_KRyR = 1
        %% Parameters
            param.VPLC =0.007;
            % Parameters of diffusion
            Dc=5.0;     % Calcium
            Dp=283;%300.0;   % Ip3
            param.Dp = Dp;
            De=5.0;     % Calcium in the ER

            % Parameters of transmission
            Fc = 0;     % Transmission parameter Calcium
            Fip = 0;%10^(exp_Fip);% Transmission parameter Ip3
            param.Fc = Fc;
            param.Fip = Fip;

            % Parameters for time scheme
            delt = 0.05;  % Time step
            tend = 50; % 302;   % Total time
            t_keep = tend;

        %% Number of cells
            % needs to be at max the number of cells

            
            cells_to_simulate = 1:1:7;
            param.cells_to_simulate = cells_to_simulate;
            n_cells = max(size(cells_to_simulate));
            
            %% Parameters of the cells

            % RyR
            
            param.V_RyR=0.15;
           
            param.K_RyR=0.2;
            param.K_RyR2 = 17;
            param.d_RyR = 1;
            param.n_RyR = 4;
            param.m_RyR = 0;
            
            
            % IPR
            param.k_beta =0.4;
            param.K_p = 0.2;
            param.K_c = 0.2;
            param.K_h = 0.08;
            param.kIPR = 45;
            param.kIPR_min = 0;
            
            % Serca
            param.V_p = 0.9;
            param.k_p = 0.2;
            param.K_bar = 0.00001957;
            
            % Calcium
            param.gamma=0.185;

            % IP3
            param.V_3K=0.05;
            param.V_5K=0.05;
            param.K_PLC = 0.07;
            param.K3K=0.4;

            % h
            param.tau_max = 100;
            param.K_tau = 0.1;
            
            % h RyR
            param.K_hRyR = 0.15;
            param.tau = 1;

            param.IPRdn = 0.8;
            param.IPRdf = 0.11;
            param.PLCds = 0.8;
            param.PLCdl = 1.2;


        %% Initial conditions

            c0_init = [0.0717];
            h_init = param.K_h^4./(param.K_h^4+c0_init.^4);
            param.ct = 3;
            
            c0 = [c0_init,c0_init*ones(1,n_cells-1)]; % Initial condition Calcium in each cell
            ip_init = 0;%.5*(2e-4/0.1)*((param.K3K^2+c0_init^2)/(c0_init^2));
            ip0 = [ip_init,ip_init*ones(1,n_cells-1)]; % Initial condition Ip3 in each cell
            ce0 = ((param.ct-c0)/param.gamma);            % Initial condition ER Calcium in each cell
            h0 = 0*h_init*ones(1,n_cells);
            hRyR0 = 0*c0.^2./(c0.^2+param.K_hRyR^2);
            
            
%             for j = 1:7
%                 clust_cell(j) = Cell(j,pa);
%             end
%             IC = [932.091078856757,25.1677793230965,120.759224637064,50.1419328613779,10.0648588190545,0.000124703178788534,1031.52379977483,25.9199705639459,120.076778386228,49.5479003504709,10.4005361880529,0.000133957613013055,457.709196680091,31.5538862321531,114.366863796841,43.1408259960324,9.90875357389654,0.000134280939623646,1104.25588640659,26.5948909518354,119.317710739548,49.0061846131914,10.6373931600690,0.000141201494527687,1087.50081493686,25.8451584156504,120.018985336922,49.7408791720075,10.4667201122196,0.000135294785319313,990.562160228820,27.8649098374267,118.000347983568,47.5990966133658,10.6858640709528,0.000145521997031992,996.363670964924,26.3668476652166,119.514897933036,49.0632145460712,10.4342437607514,0.000135828436796344,120.814355490194,4.72530872743208,120.937177931981,4.73819477704778,120.938557647822,4.73852671961441,120.891998539914,4.73398451239475,122.211931255056,4.76622321003168,120.802037516969,4.71865539565386,120.760207696534,4.72068310147413,120.929836250664,4.73696333401751,120.838186346924,4.73091055123133,120.841727256752,4.73326374568338,120.865185372808,4.73416653620952,120.808594212337,4.73161812865412,120.843710692305,4.73316123185455,120.863524165038,4.73327557086435,120.822165471082,4.73006703127439,120.752315339664,4.72614675411527,120.750902986183,4.72079562990252,120.820175715921,4.73052553008102,120.813939545093,4.73260191014275];
            IC = pa.IC_new;

        %% Name of the meshes
            file_mesh = {'out_N4_p3-p2-p4-7tet.msh'};
            file_mesh = {'out_N4_p3-p2-p4-1tet.msh','out_N4_p3-p2-p4-2tet.msh','out_N4_p3-p2-p4-3tet.msh',...
                'out_N4_p3-p2-p4-4tet.msh','out_N4_p3-p2-p4-5tet.msh',...
                'out_N4_p3-p2-p4-6tet.msh','out_N4_p3-p2-p4-7tet.msh'};

            show = 0;
            mod_basal = 1;


        %% save_file

        str = [];
        
        fields = fieldnames(param);

        for i = 1:numel(fields)
          str = [str,fields{i},num2str(param.(fields{i}))];
        end
        str(ismember(str,' ,:;!')) = [];
        str(ismember(str,'.')) = 'p';
        str = [str,file_mesh];
        str = ['Fc',num2str(Fc),'Fip',num2str(Fip),'t',num2str(tend),...
            'nc',num2str(n_cells),DataHash(str)]

        save(str)
        end
    end
end







