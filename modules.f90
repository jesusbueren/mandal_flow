module dimensions
    implicit none
    integer,parameter::P_max=8 ! Set the maximum number of plots in an adjacency
    integer,parameter::K=5,par=5,M=2,COV=4,types_a=4 !K: points of support of flow; M:types of moonzoons; type_a: types of areas
    integer,parameter::sims=20
end
    
module cadastral_maps
    use dimensions
    implicit none
    double precision:: rho_sc=0.9d0
    integer,parameter::plots_in_map=1909,villages=14,unobs_types=3,wealth_quantiles=4
    integer,parameter,dimension(villages)::plots_v=(/1794,302,912,517,292,535,939,637,405,837,973,1844,443,1909/) !plots in each village
    double precision,dimension(villages):: mean_area
    integer,dimension(plots_in_map,plots_in_map,villages)::neighbors_map
    integer,dimension(plots_in_map,2,villages,sims)::PA_type !number of neighbors, a type for each plot in the map
    integer,dimension(plots_in_map,villages,sims)::unobs_types_i
    integer,dimension(plots_in_map,P_max,villages,sims)::neighbors !identify neighbors for each plot in the map
    !Zombie pr
    double precision,dimension(types_a):: pr_non_zombie=0.50d0
    integer,dimension(plots_in_map,villages,sims)::active_plots
    character(len=92)::file_map="C:\Users\jbueren\Google Drive\overdrilling\fortran\mandal_flow\cadastral_maps\fortran_files\"
    end
    
module primitives
    use dimensions; use cadastral_maps
    implicit none
    character(len=74)::path_primitives="C:\Users\jbueren\Google Drive\overdrilling\fortran\mandal_flow\primitives\"
    !q: Discharge distribution (flow for each point of the support)
    double precision,dimension(K,1)::q
    !PI_k: Discharge pr. (first position indicates 1 well) depending on n0 of wells
    double precision,dimension(2*P_max,K,M,villages,unobs_types)::PI_k
    !PI_m: unconditional pr of high and low monzoon
    double precision,dimension(M,villages)::PI_m 
    !PI_f: Pr. of failure (first position indicates 1 well; last position indicates all plots have 2 wells)
    double precision,dimension(2*P_max,M,villages,unobs_types)::PI_fm
    double precision,dimension(2*P_max-1,3,P_max,villages,unobs_types)::PI_f_v
    !PI_s: Pr. of success (first position indicates no wells; last position indicates all plots 2 wells but one with one well)
    double precision,dimension(2*P_max,villages)::PI_s
    double precision,dimension(2*P_max-1,3,P_max,villages)::PI_s_v
    !c_d: fixed cost of failing to drill;c_s: fixed cost of succeeding to drill; c_e: cost of electricity by well
    double precision::c_s=43.0d0,beta=0.95d0,c_d=22.0d0,c_e=8.5d0 
    !extreme value distribution shocks
    double precision,parameter::gamma=0.577215664901533d0
    double precision::v_nod=0.0d0
    double precision,dimension(2)::rho=3.0d0
    !area of plots
    double precision,dimension(types_a)::area=(/1.0d0,2.0d0,3.0d0,5.1d0/) 
    double precision,dimension(types_a-1)::area_lims=(/1.3d0,2.3d0,4.0d0/) 
    !pr of unobserved heterogeneity type
    double precision,dimension(unobs_types)::pr_unobs_t=1.0d0/dble(unobs_types)
    double precision::pr_z_type2_to_pr_z=1.0d0
    !Taxation parameters
    double precision::T_g=0.0d0,tau=0.0d0
    double precision,dimension(2*P_max-1,P_max,types_a,villages,unobs_types)::smthg=0.06d0
    integer::social
    
end
     
module simulation
use cadastral_maps
    implicit none
    character(len=68)::path_estimation="C:\Users\jbueren\Google Drive\overdrilling\fortran\mandal_flow\data\"
    character(len=71)::path_results="C:\Users\jbueren\Google Drive\overdrilling\fortran\mandal_flow\Results\"
    !maximum number of functioning wells seen in the data drill_export_.xls
    integer,parameter::max_NFW=10,simulations=1
    !Parameters of the simulated panel-> I: number of indv; T: number of periods
    integer,parameter::T_sim=5,plots_i=1052
    ! dec_it: drilling decision
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,villages,unobs_types)::F_est
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP_est
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP_final=0.0d0
    integer::bootstrap=0
    
    !Data
    integer,dimension(plots_i)::V_type,P_type,A_type,impute_i,can_be_zombie_i
    double precision,dimension(plots_i)::wealth_i
    double precision,dimension(unobs_types,plots_i)::UHE_type
    double precision,dimension(unobs_types,plots_i)::UHE_type_model
    integer,dimension(plots_i)::modal_UHE_type
    integer,dimension(T_sim,plots_i,simulations)::drilling_it
    integer,dimension(T_sim,plots_i)::n_data,wealth_q !number of wells in reference plot
    double precision,dimension(max_NFW+1,T_sim,plots_i)::Pr_N_data !pr of number of functioning wells in the adjacency
    integer,dimension(T_sim,plots_i)::MNF
    integer,dimension(T_sim,plots_i)::modal_N !pr of number of functioning wells in the adjacency
    double precision,dimension(plots_i)::N_bar
    double precision,dimension(types_a,2)::moment_own_nxa_data
    double precision,dimension(wealth_quantiles)::moment_w_data
    double precision,dimension(types_a,3,villages)::shares_n_a_v
    double precision,dimension(types_a,villages,wealth_quantiles)::shares_w_v
    double precision,dimension(villages)::shares_v_n
    double precision,dimension(2*P_max-1)::pr_N_n_data
    
    !Unobsverded heterogeneity from beliefs
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::Pr_u_X
    
    double precision:: max_mle=99999999.0d0
    integer::it_est=1
    
end module    
    