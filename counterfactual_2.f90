subroutine counterfactual_2(params_MLE)
    use dimensions; use cadastral_maps; use simulation; use primitives
    implicit none
    double precision,dimension(par)::params_true,params_MLE,params
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,villages,unobs_types)::F_true
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP_true
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types)::V_fct,V_social
    integer,dimension(plots_in_map,villages)::n_dist
    double precision,dimension(villages)::mean_N,social_output,private_output
    double precision::mle
    integer::v_l,p_l,it,i_l,j_l,ind
    character::end_key
    integer,parameter::nkk=25
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
    
    ! rho_sc=params_MLE(6)
     shrinkage_p=params_MLE(5)
     logit_constrain_p=(/params_MLE(4),params_MLE(6)/)
     

     call load_cadastral_maps()
    
    tau_grid(1)=-5.0d0
    do p_l=2,nkk
        tau_grid(p_l)=tau_grid(p_l-1)+1.0d0
    end do
    
    
        
    !I want to compute the optimal tax of production giving a subsidy as a lumpsum.
    !Set a tax to production, find the lumpsum transfer that makes government transfer to be in eq.
    !Look for the optimal tax that maximizes average NPV
    tau=0.0d0
    do v_l=1,villages;do p_l=1,nkk
        print*,'exp',p_l
        print*,'village,',v_l 
        tau=tau_grid(p_l)
        tau_per_N=0.0d0!tau_grid(p_l)
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

        call compute_eq_F_CCP(params_MLE,F_true(:,:,:,:,:,v_l,:),CCP_true(:,:,:,:,v_l,:),V_fct,V_social,n_dist(:,v_l),v_l,mean_N(v_l),social_output(v_l),private_output(v_l),Pr_u_X(:,:,:,:,v_l,:),pr_d_a_n,pr_N_n,pr_na)
        
        if (p_l==1 .and. v_l==1) then
            OPEN(UNIT=12, FILE=path_results//"counterfactuals_new.txt")
            write(12,'(F10.3,I4,F10.3,F10.3,F10.3,F10.3,F10.3)'),tau,v_l,mean_N(v_l),social_output(v_l),private_output(v_l)
            close(12)
        else
            OPEN(UNIT=12, FILE=path_results//"counterfactuals_new.txt",access='append')
            write(12,'(F10.3,I4,F10.3,F10.3,F10.3,F10.3,F10.3)'),tau,v_l,mean_N(v_l),social_output(v_l),private_output(v_l)
            close(12)
        end if        
    end do;end do


    
    end subroutine
    

     
     
