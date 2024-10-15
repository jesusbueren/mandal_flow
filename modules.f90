module dimensions
    implicit none
    integer,parameter::P_max=8 ! Set the maximum number of plots in an adjacency
    integer,parameter::K=5,par=5,M=2,COV=4,types_a=5 !K: points of support of flow; M:types of moonzoons; type_a: types of areas
    integer,parameter::sims=5,sims_tr=250 !10 !250  !sims needs to be smaller than sims_tr
end
    
module cadastral_maps
    use dimensions
    implicit none
    double precision:: rho_sc=0.9d0,shrinkage_p
    double precision,dimension(2):: logit_constrain_p
    integer,parameter::plots_in_map=1909,villages=14,unobs_types=2,wealth_quantiles=3
    !plots in each village
    integer,parameter,dimension(villages)::plots_v=(/1794,302,912,517,292,535,939,637,405,837,973,1844,443,1909/) 
    double precision,dimension(villages):: total_area
    integer,dimension(plots_in_map,plots_in_map,villages)::neighbors_map
    integer,dimension(plots_in_map,2,villages,sims_tr)::PA_type !number of neighbors, a type for each plot in the map
    integer,dimension(plots_in_map,villages,sims_tr)::unobs_types_i
    integer,dimension(plots_in_map,P_max,villages,sims_tr)::neighbors !identify neighbors for each plot in the map
    
    integer,dimension(plots_in_map,villages,sims_tr)::active_plots,wealth_plots
    double precision,dimension(villages,plots_in_map)::areas
    character(len=92+10)::file_map="C:\Users\jbueren\Google Drive\overdrilling\fortran\representative_sample\cadastral_maps\fortran_files\"
    integer::Island_economy=0
    end
    
module primitives
    use dimensions; use cadastral_maps
    implicit none
    character(len=74+10)::path_primitives="C:\Users\jbueren\Google Drive\overdrilling\fortran\representative_sample\primitives\"
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
    double precision::c_s= 59.8d0,beta=0.95d0,c_d=28.8d0,c_e=8.5d0
    !extreme value distribution shocks
    double precision,parameter::gamma=0.577215664901533d0
    double precision::v_nod=0.0d0
    double precision,dimension(2)::rho=3.0d0
    !area of plots
    double precision,dimension(types_a)::area=(/0.5d0,1.0d0,2.0d0,3.0d0,5.3d0/) !(/1.3d0,4.0d0/) 
    double precision,dimension(types_a-1)::area_lims=(/0.9d0,1.5d0,2.3d0,4.0d0/)!(/2.3d0/) 
    !wealth quantiles
    double precision,dimension(wealth_quantiles-1)::w_lims=(/12.59d0,13.33d0/) 
    !pr of unobserved heterogeneity type
    double precision,dimension(unobs_types)::pr_unobs_t=(/0.364d0, 1.0d0-0.364d0/)!1.0d0/dble(unobs_types)
    double precision::pr_developed=0.3d0
    !Taxation parameters
    double precision::T_g=0.0d0,tau=0.0d0,tau_per_N=0.0d0
    double precision,dimension(2*P_max-1,P_max,types_a,villages,unobs_types)::smthg=0.06d0
    integer::social
    double precision,dimension(villages,types_a,3)::pr_v_na
    double precision,dimension(villages,types_a)::pr_v_a
    double precision,dimension(villages)::pr_v
    double precision,dimension(types_a)::Pr_D_a_data,Pr_D_a_var
    
end
     
module simulation
use cadastral_maps
    implicit none
    character(len=68+10)::path_estimation="C:\Users\jbueren\Google Drive\overdrilling\fortran\representative_sample\data\"
    character(len=71+10)::path_results="C:\Users\jbueren\Google Drive\overdrilling\fortran\representative_sample\Results\"
    !maximum number of functioning wells seen in the data drill_export_.xls
    integer,parameter::max_NFW=10,simulations=1
    !Parameters of the simulated panel-> I: number of indv; T: number of periods
    integer,parameter::T_sim=5,plots_i=1052
    ! dec_it: drilling decision
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,villages,unobs_types)::F_est
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP_est
    !double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP_final=0.0d0
    integer::bootstrap=0,bs_l
    character(len=3)::n2s
    integer,parameter::bs_samples=100
    double precision,dimension(types_a,2,bs_samples)::moment_own_nxa_bs
    double precision,dimension(types_a,bs_samples)::moments_D_bs
    
    !Data
    integer,dimension(plots_i)::V_type,P_type,A_type,impute_i,can_be_zombie_i,wealth_q
    double precision,dimension(plots_i)::wealth_i
    double precision,dimension(unobs_types,plots_i)::UHE_type
    double precision,dimension(unobs_types,plots_i)::UHE_type_model
    integer,dimension(plots_i)::modal_UHE_type
    integer,dimension(T_sim,plots_i,simulations)::drilling_it
    integer,dimension(T_sim,plots_i)::n_data !number of wells in reference plot
    double precision,dimension(max_NFW+1,T_sim,plots_i)::Pr_N_data !pr of number of functioning wells in the adjacency
    integer,dimension(T_sim,plots_i)::MNF
    integer,dimension(T_sim,plots_i)::modal_N !pr of number of functioning wells in the adjacency
    double precision,dimension(plots_i)::N_bar
    double precision,dimension(types_a,2)::moment_own_nxa_data,var_own_nxa
    double precision,dimension(wealth_quantiles,types_a)::moment_wa_data,var_wa_data
    double precision,dimension(villages)::shares_v
    double precision,dimension(types_a,3,villages)::shares_n_a_v
    double precision,dimension(types_a,villages)::shares_a_v
    double precision,dimension(types_a,villages,wealth_quantiles)::shares_w_v
    double precision,dimension(villages,3)::shares_v_n
    double precision,dimension(villages,wealth_quantiles,3,types_a)::share_v_wn
    double precision,dimension(2*P_max-1,3)::pr_N_n_data
    double precision,dimension(2*P_max-1)::var_N
    double precision,dimension(3)::pr_little_n_data,var_little_n
    double precision,dimension(plots_i,types_a)::p_w_a
    integer,dimension(types_a)::counter_w_a
    double precision,dimension(3,villages,types_a)::pr_D_data
    
    !Unobsverded heterogeneity from beliefs
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::Pr_u_X
    
    double precision:: max_mle=99999999.0d0
    !generate artificial data
    integer::it_est=1,generate_data=0,montecarlo=0
    
end module    
    
