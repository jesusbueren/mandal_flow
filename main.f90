program main
use dimensions; use cadastral_maps; use simulation
implicit none
integer(8),dimension(1)::seed=321
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
call random_seed(PUT=seed)

!Input primitives of the model
call input_primitives()

!Load panel data of drilling decisions
call load_estimation_data()

!Create a map with plots, neighbors and types
!call load_cadastral_maps()




!Uncomment if you want to simulate and check recovery of the simulated parameters.
    !Choose true parameter vector and compute associated equilibrium beliefs and a final distribution of wells in all plots
    !params_true=(/8.15d0,0.84d0,6.18d0/)
    !
    !do v_l=1,1!villages
    !    print*,'village,',v_l
    !    if (v_l==1) then
    !        ccp_true(:,:,:,:,v_l,:)=0.07d0
    !        n_dist(:,v_l)=1
    !        v_fct=0.0d0
    !    else
    !        ccp_true(:,:,:,:,v_l,:)=ccp_true(:,:,:,:,1,:)
    !        n_dist(:,v_l)=n_dist(:,1)
    !    end if        
    !    call compute_eq_f_ccp(params_true,f_true(:,:,:,:,:,v_l,:),ccp_true(:,:,:,:,v_l,:),v_fct,v_social,n_dist(:,v_l),v_l,mean_n(v_l),social_output(v_l),private_output(v_l),pr_u_x(:,:,:,v_l,:))
    !end do
    !
    !!Using CCPs simulate panel data with drilling decision.
    !OPEN(UNIT=12, FILE="CCP_true.txt")
    !    write(12,*),CCP_true,n_dist,F_true,Pr_u_X
    !close(12)
    !!do ns=1,100
    !call random_seed(PUT=seed)
    !call simulate_panel(CCP_true,n_dist,seed)
    !call random_seed(get=seed)
    
    !print*,'Loading panel'
!
!
!print*,'end simulation'
!Compute moments in the data
call compute_moments(dble(drilling_it(:,:,1)),"data",empty,moment_own_nxa_data,moment_w_data)
OPEN(UNIT=12, FILE=path_results//"data"//"_own_nxa.txt")
    write(12,*),moment_own_nxa_data
close(12)
OPEN(UNIT=12, FILE=path_results//"data"//"_wealth.txt")
    write(12,*),moment_w_data
close(12)

print*,'Start estimation'
!Generate a random CCP for computing initial beliefs
CCP_est=sqrt(-1.0d0)
do P_l=1,P_max
    CCP_est(1:2*P_l-1,1:2,P_l,:,:,1)=0.07d0
end do

!call estimation2(params_MLE,log_likeli)


open(unit=12, file=path_results//"parameters.txt")
    read(12,'(<par>f20.12)'),params_MLE
close(12)


!call bootstrap_se()
!open(unit=12, file=path_results//"bootstrapped_parameters_85_nem.txt")
!    read(12,*),params_mle
!close(12)
print*,'estimated parameters',params_MLE
!end do
!call generate_panel_sample(params_MLE)
call counterfactual_2(params_MLE)
!call counterfactual_1(params_MLE)

!call transitional_dynamics(params_MLE)
!call quantifiying_strategic_interactions(params_MLE)

read*,end_key

end program main