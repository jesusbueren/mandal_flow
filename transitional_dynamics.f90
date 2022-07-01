subroutine transitional_dynamics(params_MLE)
    use dimensions; use cadastral_maps; use simulation; use primitives
    implicit none
    integer,parameter::T=500,Sims2=1000
    double precision,dimension(par)::params_true,params_MLE
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,unobs_types,T)::F_in,F_out
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,unobs_types)::slope,intercept
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types,T)::CCP_true
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types,T)::V_fct,V_fct2
    integer,dimension(plots_in_map,T)::n_dist
    integer,dimension(plots_in_map,Sims2)::n_ini
    double precision,dimension(T)::mean_N,social_output,private_output,ccp_mean,diss_N
    integer::v_l,p_l,it,t_l,s_l
    character::end_key
    
    !rho=params_MLE(3)
    !Compute the value functions, beliefs and ditribution of wells for the two steady states
    v_l=1
    tau=0.0d0
    t_l=1
    CCP_true(:,:,:,:,v_l,:,t_l)=0.07d0
    call compute_eq_F_CCP(params_MLE,F_in(:,:,:,:,:,:,t_l),CCP_true(:,:,:,:,v_l,:,t_l),V_fct(:,:,:,:,:,t_l),V_fct2(:,:,:,:,:,t_l),n_dist(:,t_l),v_l,&
                                                mean_N(t_l),social_output(t_l),private_output(t_l),Pr_u_X(:,:,:,:,v_l,:))
    print*,'NPV at begining',social_output(t_l)
    
    tau=c_e
    t_l=T
    CCP_true(:,:,:,:,v_l,:,t_l)=0.07d0
    call compute_eq_F_CCP(params_MLE,F_in(:,:,:,:,:,:,t_l),CCP_true(:,:,:,:,v_l,:,t_l),V_fct(:,:,:,:,:,t_l),V_fct2(:,:,:,:,:,t_l),n_dist(:,t_l),v_l,&
                                                mean_N(t_l),social_output(t_l),private_output(t_l),Pr_u_X(:,:,:,:,v_l,:))
    print*,'NPV at end',social_output(t_l)
    
    !Initial guess of beliefs
    slope=(F_in(:,:,:,:,:,:,T)-F_in(:,:,:,:,:,:,1))/(dble(T)-1.0d0)
    intercept=F_in(:,:,:,:,:,:,1)-slope
    do t_l=2,T-1
        F_in(:,:,:,:,:,:,t_l)=slope*dble(t_l)+intercept
    end do
    do s_l=1,Sims2
        n_ini(:,s_l)=n_dist(:,1)
    end do
    call solve_path(params_MLE,T,Sims2,n_ini,F_in,v_l,V_fct,mean_N,social_output,ccp_mean,diss_N)    

end subroutine
    
subroutine solve_path(params,T_path,Sims2,n_ini,F_in,v_l,V_in,mean_N,social_output,ccp_mean,diss_N)
    use cadastral_maps; use dimensions; use primitives
    implicit none
    integer,intent(in)::T_path,Sims2
    double precision,dimension(par),intent(in)::params
    integer,dimension(plots_in_map,Sims2),intent(inout)::n_ini
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,unobs_types,T_path),intent(inout)::F_in
    integer,intent(in)::v_l  
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types,T_path),intent(in)::V_in
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types,T_path)::V_fct
    double precision,dimension(T_path),intent(out)::mean_N,social_output,ccp_mean,diss_N
    
    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types,T_path)::CCP_old,CCP,CCP_mid
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::Ef_v !Ef_v: expected productivity
    double precision::dist
    integer::p_l,a_l,n_l,P_l2,ind,counter_all,counter_bad,u_l,t_l
    integer(8),dimension(2*P_max-1,3,3,P_max,T_path)::iterations
    
    print*,'got into solve path'
    !Set scale parameter Gumbel distribution of shocks
    !rho=params(4)

    !Compute expected productivity 
    do u_l=1,unobs_types;do a_l=1,types_a
        call expected_productivity(params(1:3),area(a_l),Ef_v(:,:,:,a_l,v_l,u_l),v_l,u_l)
    end do;end do
    
    !Taking as given the value function of the next period compute the value function and CCP of the next period
    CCP_old=0.0d0
1   do P_l=2,P_max; do a_l=1,types_a; do u_l=1,unobs_types;do t_l=T_path,1,-1
        if (t_l>20) then
            tau=c_e
        else
            tau=0.0d0
        end if        
        if (t_l==T_path) then
            call one_step_value_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,v_l,u_l)&
                            ,F_in(1:2*P_l-1,1:2*P_l-1,:,:,P_l,u_l,t_l) &
                            ,P_l &
                            ,CCP(1:2*P_l-1,:,P_l,a_l,u_l,t_l),v_l,u_l &
                            ,V_in(1:2*P_l-1,:,P_l,a_l,u_l,t_l) &
                            ,V_fct(1:2*P_l-1,:,P_l,a_l,u_l,t_l))
        elseif (t_l<=20) then
            call one_step_value_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,v_l,u_l)&
                                ,F_in(1:2*P_l-1,1:2*P_l-1,:,:,P_l,u_l,t_l) &
                                ,P_l &
                                ,CCP(1:2*P_l-1,:,P_l,a_l,u_l,t_l),v_l,u_l &
                                ,V_in(1:2*P_l-1,:,P_l,a_l,u_l,1) &
                                ,V_fct(1:2*P_l-1,:,P_l,a_l,u_l,t_l))
        else           
            call one_step_value_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,v_l,u_l)&
                                ,F_in(1:2*P_l-1,1:2*P_l-1,:,:,P_l,u_l,t_l) &
                                ,P_l &
                                ,CCP(1:2*P_l-1,:,P_l,a_l,u_l,t_l),v_l,u_l &
                                ,V_fct(1:2*P_l-1,:,P_l,a_l,u_l,t_l+1) &
                                ,V_fct(1:2*P_l-1,:,P_l,a_l,u_l,t_l))
        end if
     end do; end do;end do;end do
    
    call generate_transition_beliefs(T_path,Sims2,CCP,Ef_v(:,:,:,:,v_l,:),n_ini,F_in,v_l,V_fct,iterations,mean_N,social_output,ccp_mean,diss_N)
    F_in(:,:,:,:,:,:,20)=F_in(:,:,:,:,:,:,1)
    print*,'NPV at begining',social_output(1)
    print*,'NPV at end',social_output(T_path)
    dist=0.0
    t_l=1
    do P_l=2,P_max; do n_l=1,2;do ind=1,2*P_l-1; do t_l=1,T_path;
        dist=dist+dble(sum(iterations(ind,n_l,1:3,P_l,t_l)))/dble(sum(iterations(:,1:2,1:3,:,:)))*sum(abs(CCP_old(ind,n_l,P_l,:,:,t_l)-CCP(ind,n_l,P_l,:,:,t_l)))/dble(types_a)/dble(unobs_types)
    end do;end do; end do;end do
    
    OPEN(UNIT=12, FILE="transitional_dynamics.txt")
    do t_l=1,T_path
        write(12,'(F10.4,I4,F10.4,F10.4,F10.4,F10.4)'),tau,v_l,mean_N(t_l),social_output(t_l),ccp_mean(t_l),diss_N(t_l)
    end do
    close(12)
    
    print*,'dist CCP',dist
    !Compute beliefs over the transition path given ccps
    if (dist>1.0d-4) then
        CCP_old=CCP
        go to 1
    end if

end subroutine
    
