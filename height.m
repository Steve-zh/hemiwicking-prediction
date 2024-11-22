clear;
clc;
h_fig = figure('Name', [ 'Experiments_l_change_height']);
c=[0 0 0;1 0 0]; %color
data=xlsread('experiments.xlsx');
errorbar(data(:,1),data(:,2)*1000,data(:,4)*1000/2.0,data(:,4)*1000/2.0,data(:,3)/2.0,data(:,3)/2.0,'o','color',c(2,:),'MarkerSize', 9,'linewidth',0.5,'CapSize',10);
hold on;
T_cond = 20;
g = 9.81;

% Properties for Water (REFPROP curve fitting)
h_fg = -1.376e-05*T_cond^4 - 0.009548*T_cond^3 + 0.361*T_cond^2 - 2372*T_cond + 2.501e+06;
rho_l = -8.7425E-08*T_cond^4 + 3.27336E-05*T_cond^3 - 0.006930398*T_cond^2 + 0.037679074*T_cond + 999.9096389;
rho_v = 4.777e-09*T_cond^4 - 1.312e-07*T_cond^3 + 2.359e-05*T_cond^2 + 0.0001006*T_cond + 0.005891;
mu_l = 1.87739E-11*T_cond^4 - 6.19172E-09*T_cond^3 + 7.88134E-07*T_cond^2 - 5.03822E-05*T_cond + 0.001747344;
k_l = - 2.452e-10*T_cond^4 + 8.725e-08*T_cond^3 - 1.837e-05*T_cond^2 + 0.002421*T_cond + 0.5561;
gamma_l = 2.5414E-13*T_cond^4 + 1.75636E-10*T_cond^3 - 2.96437E-07*T_cond^2 - 0.000139726*T_cond + 0.075647774;
L_plate=0.017;
load("Lookup_Table.mat");
for geometry = 1:1

    if geometry ==1 %dimension of micropillar
        D_wick = 0.000007;
        l_pitch = 0.000020;
        t_wick = 0.000008;
    end
    eps = 1 - pi*D_wick^2/4/l_pitch^2;   % Porosity 
    s = 1 - eps ;             % Solid Fraction

    %% Effective Height
    d = D_wick;
    l = l_pitch;
    h = t_wick;
    
    d_l_ratio = floor(l/d);

    % Locations in Lookup Table: 1st column:d/l ratio,2nd column:theta
    %5th column: effective height reduction,6th column:meniscus interface
    %projected area
    if l == 0.000010 && d == 0.000007
        aaa = 1;
    elseif l == 0.000020 && d == 0.000007
        aaa = 18;
    elseif l == 0.000040 && d == 0.000007
        aaa = 35;
    elseif l == 0.0000081 && d == 0.00000282
        aaa = 52;
    elseif l == 0.0000055 && d == 0.0000029
        aaa = 69;    
    elseif l == 0.00001 && d == 0.000005
        aaa = 86;
    elseif l == 0.00004 && d == 0.00002
        aaa = 86;
    elseif l == 0.00002 && d == 0.000005
        aaa = 103;    
    end   

    %% Initial Constant Permeability 

    % Sangani & Acrivos
    k_2D =   l^2 * ( log(s^-0.5) - 0.7105 + 0.6648*s + 0.7086*s^2 - 1.097*s^3)/4/pi ;    
    beta = sqrt(eps/k_2D);

    % Kappa for last cell (receding CA = 15)
    kappa = k_2D ;
    %% Plot Figures
    theta = 15;     
    dP_cap_max_YL = gamma_l*pi*D_wick*cosd(theta)/(l_pitch*l_pitch-pi*(D_wick/2)^2);  

    t_3d = 0.0;
    N = floor (L_plate / l_pitch);
    num_ele = 1;
    t_3Dx = zeros(1,N);
    length_list = zeros(1,N);
    length = num_ele*l_pitch;

    t_2d = 0.0;
    t_2Dx = zeros(1,N);
    while length <= L_plate

        %% Initial Values

        P_org = @(y) -dP_cap_max_YL/length *y ;

        for i = 1:num_ele+1
            P_r2(i) = P_org(l_pitch * (i-1));  
        end

        for i = 1:1:num_ele             % i stands for i_th element on the plate         
            P_local(i) = ( P_r2(i) + P_r2(i+1) ) / 2 ;             
            K_local(i) = -P_local(i)/2/gamma_l *d/2;               
            cos_theta(i) = K_local(i) * ( ( l/(d/2) )^2 - pi*1^2 ) / pi/ 1 ;   
            theta(i) = acosd( cos_theta(i) );        
        end
            % Refer to Lookup Table     
        for i = 1:1:num_ele

            for iii = aaa:1:aaa+16       
                dtheta(i) = theta(i) - Lookup_Table(iii,2) ;   

                if dtheta(i) < 5
                   dH(i) = d/2 * (Lookup_Table(iii,5) + 0.2*dtheta(i)*( Lookup_Table(iii+1,5) - Lookup_Table(iii,5)) );
                   H_eff(i) = h + dH(i);
                   AR(i) = Lookup_Table(iii,6) + 0.2*dtheta(i)*( Lookup_Table(iii+1,6) - Lookup_Table(iii,6) ) ;
                   dH2(i) = d/2 * (Lookup_Table(iii,4) + 0.2*dtheta(i)*( Lookup_Table(iii+1,4) - Lookup_Table(iii,4)) );
                   H_eff2(i) = h + dH2(i);
                   dH3(i) = d/2 * (Lookup_Table(iii,7) + 0.2*dtheta(i)*( Lookup_Table(iii+1,7) - Lookup_Table(iii,7)) );
                   H_eff3(i) = h + dH3(i);
                  
                   break
                end
            end

        end
        u_total3d = 0.0;
        u_total2d = 0.0;
        for i=1:num_ele
            kappa_3D = kappa/mu_l/eps*( 1 - ( exp( 2*H_eff3(i)*beta ) - 1 )/( H_eff3(i)*beta * (exp( 2*H_eff3(i)*beta ) + 1)));
            u_total3d = u_total3d+kappa_3D* (dP_cap_max_YL/length);
        end
        u_3d = u_total3d/num_ele;
        %
        t_3Dx(num_ele) = t_3d;
        del_t_3d = l_pitch / u_3d;
        t_3d = t_3d + del_t_3d;

        length_list(num_ele) = length*1000;
        num_ele =num_ele+1;
        length = num_ele*l_pitch;

    end
    plot(t_3Dx,length_list,'linestyle','-','linewidth',2.5,'Color',c(geometry,:));      %velocity_average of each unit-15o receding angle
    hold on;
    t_3Dx=t_3Dx(2:end);
    length_list=length_list(2:end);
end
set(gcf,'unit','normalized','position',[0.2,0.2,0.4,0.5]);
set(gca,'FontName','Times New Roman','FontSize',15);
xlabel('Time (s)','FontName','Times New Roman','FontSize',22,'FontWeight','bold','Color','k')
ylabel('Hemiwicking length (mm)','FontName','Times New Roman','FontSize',22,'FontWeight','bold','Color','k')
ylim([0 17]);
xlim([0 8]);
saveas(h_fig, h_fig.Name, 'svg')