
include "nrtype.f90"
include "nr.f90"
include "nrutile.f90"
include "modules.f90"    
include "amoeba.f90"
include "compute_eq_F_CCP.f90"
include "counterfactual_2.f90"
include "estimation2.f90"
include "expected_productivity.f90"
include "compute_moments.f90"
include "generate_beliefs.f90"
include "generate_C.f90"
include "generate_transition_beliefs.f90"
include "input_primitives.f90"
include "load_cadastral_maps.f90"
include "policy_fct_it.f90"
include "transitional_dynamics.f90"
include "update_CCP.f90"
include "valuation.f90"
include "value_fct_it.f90"
    
    
program main
use dimensions; use cadastral_maps; use simulation
implicit none
integer,allocatable::seed(:)
double precision,dimension(par)::params_true,params_MLE
double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::empty=-9.0d0
double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,villages,unobs_types)::F_true
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP_true
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types)::V_fct,V_social
    integer,dimension(plots_in_map,villages)::n_dist
    double precision,dimension(villages)::mean_N,social_output,private_output
double precision::log_likeli
integer::v_l,p_l,ns
character::end_key



!Call seed number
call random_seed(size = v_l)
allocate(seed(v_l))
seed=321
call random_seed(PUT=seed)

!Define the primitives (flow, fail, success)
call input_primitives()
    
!Compute moments in the data
call compute_moments()

!Create a map with plots, neighbors and types
call load_cadastral_maps()

!Simulate data from a DGP and check estimation code recovers the DGP in a Montecarlo
call monte_carlo()

print*,'Start estimation'
!Generate a random CCP for computing initial beliefs
CCP_est=-9.0d0
do P_l=1,P_max
    CCP_est(1:2*P_l-1,1:2,P_l,:,:,1)=0.07d0
end do
!(set open mp & set mem:50000000)
bootstrap=0

!call estimation2(params_MLE,log_likeli)
open(unit=12, file=path_results//"parameters.txt")
    read(12,'(<par>f20.12)'),params_MLE
close(12)

!Bootstrap se
do bs_l=46,bs_samples
    bootstrap=1
    write (n2s, "(I3.3)") bs_l    
    call estimation2(params_MLE,log_likeli)
end do


print*,'estimated parameters',params_MLE

!call counterfactual_2(params_MLE)

!(take out open mp: set mem:800000000)
call transitional_dynamics(params_MLE)
!call quantifiying_strategic_interactions(params_MLE)

read*,end_key

end program main