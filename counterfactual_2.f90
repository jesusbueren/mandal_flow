subroutine counterfactual_2(params_MLE)
    use dimensions; use cadastral_maps; use simulation; use primitives
    implicit none
    double precision,dimension(par)::params_true,params_MLE,params
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,villages,types_a)::F_true
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP_true
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types)::V_fct,V_social
    integer,dimension(plots_in_map,villages)::n_dist
    double precision,dimension(villages)::mean_N,social_output,private_output
    double precision::mle,d_rate,delta_PV_model
    double precision,dimension(types_a)::pr_D_model
    integer::v_l,p_l,it,i_l,j_l,ind
    character::end_key
    integer,parameter::nkk=20
    double precision,dimension(nkk)::tau_grid
    double precision,dimension(plots_i,unobs_types)::posterior_type
    double precision,dimension(types_a,2)::pr_d_a_n
    double precision,dimension(2*P_max-1,3)::Pr_N_n
    double precision,dimension(3,types_a)::pr_na
    double precision,dimension(plots_i):: pr_non_zombie_i
    double precision,dimension(types_a)::pr_non_zombie_a
    integer,dimension(types_a)::counter_a
     interface 
        function log_likelihood2(params_MLE)
            use dimensions
            double precision,dimension(par),intent(in)::params_MLE
            double precision::log_likelihood2
        end function log_likelihood2
     end interface     
    
    tau_grid(1)=0.0d0
    if (nkk>1) then
        do p_l=2,nkk
            tau_grid(p_l)=tau_grid(p_l-1)+1 !1.0d0
        end do
    end if
    
    
        
    !I want to compute the optimal tax of production giving a subsidy as a lumpsum.
    !Set a tax to production, find the lumpsum transfer that makes government transfer to be in eq.
    !Look for the optimal tax that maximizes average NPV
    tau=0.0d0
    do v_l=1,villages;do p_l=1,nkk
        print*,'exp',p_l
        print*,'village,',v_l 
        tau=tau_grid(p_l)
        tau_per_N=tau_grid(p_l)
        if (p_l==1) then
            V_fct=0.0d0
            V_social=0.0d0
            !T_g=0.0d0
        end if
        
        !open(unit=12, file=path_results//"initial_beliefs.txt")
        !    read(12,*),n_dist,CCP_true
        !close(12)
        n_dist=1
        CCP_true=0.07d0

        call compute_eq_F_CCP(params_MLE,F_true(:,:,:,:,:,v_l,:),CCP_true(:,:,:,:,v_l,:),V_fct,V_social,n_dist(:,v_l),v_l,&
            mean_N(v_l),social_output(v_l),private_output(v_l),pr_d_a_n,pr_N_n,pr_na,d_rate,pr_D_model,delta_PV_model)
        
        if (p_l==1 .and. v_l==1) then
            OPEN(UNIT=12, FILE=path_results//"counterfactuals_new.txt")
            write(12,'(F10.3,I4,F10.3,F10.3,F10.3,F10.3,F10.3,F10.3)'),tau,v_l,mean_N(v_l),social_output(v_l),private_output(v_l),d_rate
            close(12)
        else
            OPEN(UNIT=12, FILE=path_results//"counterfactuals_new.txt",access='append')
            write(12,'(F10.3,I4,F10.3,F10.3,F10.3,F10.3,F10.3,F10.3)'),tau,v_l,mean_N(v_l),social_output(v_l),private_output(v_l),d_rate
            close(12)
        end if        
    end do;end do


    
end subroutine