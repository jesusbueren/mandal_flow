subroutine transitional_dynamics(params_MLE)
    use dimensions; use cadastral_maps; use simulation; use primitives
    implicit none
    integer,parameter::transition_time=700
    double precision,dimension(par)::params_true,params_MLE,params
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,types_a)::F_ini,F_final,slope,intercept
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,types_a,transition_time)::F_tr
    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types)::CCP_ini,CCP_final
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types)::V_fct_final,V_fct_initial,V_social
    integer,dimension(plots_in_map,villages)::n_dist
    double precision,dimension(villages)::mean_N,social_output,private_output
    integer::v_l,p_l,it,i_l,j_l,ind,t_l
    character::end_key
    double precision,dimension(types_a,2)::pr_d_a_n
    double precision,dimension(2*P_max-1,3)::Pr_N_n
    double precision,dimension(3,types_a)::pr_na
    double precision,dimension(2)::stats_out
    integer,parameter::nkk=20
    double precision,dimension(nkk)::tau_grid
     interface 
        function log_likelihood2(params_MLE)
            use dimensions
            double precision,dimension(par),intent(in)::params_MLE
            double precision::log_likelihood2
        end function log_likelihood2
        subroutine  solve_path(params,V_fct_final,V_fct_initial,F_tr,transition_time,v_l,CCP_ini,CCP_final,stats_out,n_tr,av_drilling_rate)  
            use cadastral_maps; use dimensions; use primitives; use simulation
            implicit none
            integer,intent(in)::transition_time,v_l
            double precision,dimension(par),intent(in)::params
            double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types),intent(in)::V_fct_final,V_fct_initial
            double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,types_a,transition_time),intent(inout)::F_tr 
            double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types),intent(in)::CCP_ini,CCP_final
            double precision,allocatable,dimension(:),intent(out)::n_tr,av_drilling_rate
            double precision,dimension(2),intent(out)::stats_out
        end subroutine solve_path
     end interface
     double precision,allocatable,dimension(:)::n_tr,n_tr_max,av_drilling_rate,av_drilling_rate_max
     
     double precision:: max_val
    
    !Assign contrained plots in map
     shrinkage_p=params_MLE(5)
     logit_constrain_p=(/params_MLE(4),params_MLE(6)/)
     call load_cadastral_maps()
      
    !Select village number
    do v_l=1,villages
        max_val=-1.0/0.0
        print*,'village number:',v_l
        !Compute initial CCP function
        tau=0.0d0
        V_fct_initial=0.0d0
        V_social=0.0d0  
        n_dist=1
        CCP_ini=0.07d0
        print*,'Initial steady-state'
        call compute_eq_F_CCP(params_MLE,F_ini,CCP_ini,V_fct_initial,V_social,n_dist(:,v_l),v_l,mean_N(v_l),social_output(v_l),private_output(v_l),Pr_u_X(:,:,:,:,v_l,:),pr_d_a_n,pr_N_n,pr_na)
  
        !Compute final CCP function for different taxes
        tau_grid(1)=0.0d0
        do p_l=2,nkk
            tau_grid(p_l)=tau_grid(p_l-1)+1.0d0
        end do
    
        do p_l=1,nkk
            print*,'policy: ',p_l
            tau=tau_grid(p_l)
            V_fct_final=0.0d0
            V_social=0.0d0  
            n_dist=1
            CCP_final=0.07d0
            print*,'Final steady-state'
            call compute_eq_F_CCP(params_MLE,F_final,CCP_final,V_fct_final,V_social,n_dist(:,v_l),v_l,mean_N(v_l),social_output(v_l),private_output(v_l),Pr_u_X(:,:,:,:,v_l,:),pr_d_a_n,pr_N_n,pr_na)
    
            !Compute initial guess of beliefs along the transitions
            slope=(F_final-F_ini)/(dble(transition_time)-1.0d0)
            intercept=F_ini-slope
            do t_l=1,transition_time
                F_tr(:,:,:,:,:,:,t_l)=slope*dble(t_l)+intercept
            end do

            !Solve for the transition
            print*,'Transition'
            call solve_path(params_MLE,V_fct_final,V_fct_initial,F_tr,transition_time,v_l,CCP_ini,CCP_final,stats_out,n_tr,av_drilling_rate)  
    
            if (p_l==1 .and. v_l==1) then
                OPEN(UNIT=12, FILE=path_results//"counterfactuals_dynamics_no_dismantle.txt")
                write(12,'(F10.3,I4,F10.3,F10.3)'),tau,v_l,stats_out(1),stats_out(2)
                close(12)
            else
                OPEN(UNIT=12, FILE=path_results//"counterfactuals_dynamics_no_dismantle.txt",access='append')
                write(12,'(F10.3,I4,F10.3,F10.3)'),tau,v_l,stats_out(1),stats_out(2)
                close(12)    
            end if 
            if (max_val<stats_out(1)) then
                max_val=stats_out(1)
                n_tr_max=n_tr
                av_drilling_rate_max=av_drilling_rate
            end if           
        end do
        if (v_l==1) then
            open(unit=12, file=path_results//"transition_N_no_dismantle.txt")
        else
            open(unit=12, file=path_results//"transition_N_no_dismantle.txt",access='append')
        end if
        do t_l=1,transition_time+150
            write(12,*),v_l,n_tr_max(t_l),av_drilling_rate_max(t_l)
        end do
        close(12)
    end do
    
    

    
end subroutine
   
   
subroutine solve_path(params,V_fct_final,V_fct_initial,F_tr,transition_time,v_l,CCP_ini,CCP_final,stats_out,n_tr,av_drilling_rate)  
    use cadastral_maps; use dimensions; use primitives; use simulation
    implicit none
    integer,intent(in)::transition_time,v_l
    double precision,dimension(par),intent(in)::params
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types),intent(in)::V_fct_final,V_fct_initial
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,types_a,transition_time),intent(inout)::F_tr 
    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types),intent(in)::CCP_ini,CCP_final
    double precision,allocatable,dimension(:),intent(out)::n_tr,av_drilling_rate
    double precision,dimension(2),intent(out)::stats_out
    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types,transition_time)::CCP_tr
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types,transition_time)::V_fct_tr
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::Ef_v
    integer::u_l,a_l,t_l,p_l,it
    integer::total_time
    double precision,allocatable,dimension(:)::n_tr_old
    
    double precision::dist
    character::pause_k
    
    
    it=0
    total_time=transition_time+150
    allocate (n_tr(total_time))
    allocate (av_drilling_rate(total_time))
    allocate (n_tr_old(total_time))
    
    n_tr_old=-9.0d0
    
    rho=params(3)
    
    !Compute expected productivity 
    do a_l=1,types_a;do u_l=1,unobs_types
        call expected_productivity((/params(2),params(1),params(7)/),area(a_l),Ef_v(:,:,:,a_l,v_l,u_l),v_l,u_l) 
    end do;end do
    
    V_fct_tr(:,:,:,:,:,transition_time)=V_fct_final
    CCP_tr=0.0d0
    CCP_tr(:,:,:,:,:,transition_time)=CCP_final
    !Taking as given the value function of the next period compute the value function and CCP of the next period
1   do t_l=transition_time,2,-1;do P_l=1,P_max; do a_l=1,types_a; do u_l=1,unobs_types
        call one_step_value_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,v_l,u_l)&
                            ,F_tr(1:2*P_l-1,1:2*P_l-1,:,:,P_l,a_l,t_l) &
                            ,P_l &
                            ,CCP_tr(1:2*P_l-1,:,P_l,a_l,u_l,t_l) &
                            ,CCP_tr(1:2*P_l-1,:,P_l,a_l,u_l,t_l-1),v_l,u_l & !CCP_tr(1:2*P_l-1,:,P_l,a_l,u_l,500)
                            ,V_fct_tr(1:2*P_l-1,:,P_l,a_l,u_l,t_l) &
                            ,V_fct_tr(1:2*P_l-1,:,P_l,a_l,u_l,t_l-1))
     end do; end do;end do;end do
   
    call generate_transition_beliefs(CCP_ini,CCP_tr,F_tr,n_tr,transition_time,total_time,v_l,V_fct_tr,V_fct_initial,Ef_v(:,:,:,:,v_l,:),stats_out,av_drilling_rate)
    
    
    
    
    dist=maxval(abs(n_tr(total_time-transition_time:total_time)-n_tr_old(total_time-transition_time:total_time)),1)
    print*,dist
    
    !print*,'press key to continue'
    !read*,pause_k
    
    if (dist>0.01d0 .and. it<5) then
        n_tr_old=n_tr
        it=it+1
        go to 1
    end if
    
    if (it==5)then
        print*,'reached max number of iterations'
    end if
    
    


end subroutine
   
