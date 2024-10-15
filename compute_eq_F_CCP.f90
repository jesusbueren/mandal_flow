subroutine compute_eq_F_CCP(params,F,CCP,V_fct_ts,V_social,n_initial,v_l,mean_N,social_output,private_output,pr_d_a_n,pr_N_n,pr_na,d_rate,pr_D_model,delta_PV_model)
    use cadastral_maps; use dimensions; use primitives
    implicit none
    double precision,dimension(par),intent(in)::params
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,types_a),intent(out)::F
    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types),intent(inout)::CCP
    integer,dimension(plots_in_map,1),intent(inout)::n_initial
    integer,dimension(plots_in_map,1)::n_start
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types),intent(inout)::V_fct_ts,V_social
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types)::V_fct
    integer,intent(in)::v_l    
    double precision,intent(out)::mean_N,social_output,private_output,delta_PV_model
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types)::joint_pr
    double precision,dimension(3,unobs_types)::pr_n_u,pr_n_u_old
    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types)::CCP_old,CCP2
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::Ef_v !Ef_v: expected productivity Ef_v(1,3,1,4,2,:)
    double precision::dist,Mean_N_old,u
    integer::p_l,a_l,n_l,P_l2,ind,counter_all,counter_bad,u_l,it,i,s_l
    integer(8),dimension(2*P_max-1,3,3,P_max,types_a)::iterations
    double precision,dimension(types_a,2),intent(out)::pr_d_a_n
    double precision,dimension(2*P_max-1,3),intent(out)::pr_N_n
    double precision,dimension(3,types_a),intent(out)::pr_na
    double precision,intent(out)::d_rate
    double precision,dimension(types_a),intent(out)::pr_D_model
    character::pause_k
    integer,allocatable::seed(:)
    interface
    subroutine random_value ( seed, r )
    implicit none
    double precision,intent(out):: r
    integer,dimension(:),intent(inout):: seed
    end
    end interface
    
        
    
    !Call seed number
    call random_seed(size = p_l)
    allocate(seed(p_l))
    seed=321321
    call random_seed(PUT=seed)

    !Define developed/undeveloped plots
    do s_l=1,sims_tr;do i=1,plots_v(v_l)
        call random_value( seed, u )
        if (u<params(5)) then
            active_plots(i,v_l,s_l)=1
        else
            active_plots(i,v_l,s_l)=0 
        end if
    end do;end do
    
    
    
    Mean_N_old=0.0d0
    rho=params(3)
    
    !Compute expected productivity 
    do a_l=1,types_a;do u_l=1,unobs_types
        call expected_productivity((/params(2),params(1),params(4)/),area(a_l),Ef_v(:,:,:,a_l,v_l,u_l),v_l,u_l)
    end do;end do
    
    if (Island_economy==1) then
        P_l=1
        ind=1
        F=0.0d0
        do n_l=1,3
            F(1,1,n_l,1:min(n_l+1,3),P_l,:)=1.0d0
        end do
        do P_l=1,1; do a_l=1,types_a; do u_l=1,unobs_types
            call value_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,v_l,u_l)& 
                                ,F(1:2*P_l-1,1:2*P_l-1,:,:,P_l,a_l) &
                                ,P_l &
                                ,CCP(1:2*P_l-1,:,P_l,a_l,u_l),v_l,u_l & 
                                ,V_fct_ts(1:2*P_l-1,:,P_l,a_l,u_l)) 
        end do;end do;end do
        return
    end if

    !Generate beliefs consitent with CCP
!    F=1.0d0

    it=0
    pr_n_u_old=0.0d0
    n_start=n_initial
    V_fct=0.0d0
1   n_initial=1
    
    call generate_beliefs(CCP,V_fct,V_social,Ef_v(:,:,:,:,v_l,:),n_initial,F,v_l,iterations,mean_N,social_output,private_output,&
        joint_pr,pr_d_a_n,pr_N_n,pr_na,d_rate,pr_D_model,delta_PV_model) 
    do u_l=1,unobs_types;do n_l=1,3
        pr_n_u(n_l,u_l)=sum(joint_pr(:,n_l,:,:,u_l))/sum(joint_pr(:,:,:,:,u_l))
    end do;end do


    
    !For each plot type obtain a new CCP given beliefs
    !print*,'policy step'
    CCP_old=CCP
    dist=0.0d0
    counter_bad=0
    counter_all=0
    
    do P_l=1,P_max; do a_l=1,types_a; do u_l=1,unobs_types
        call value_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,v_l,u_l)& 
                            ,F(1:2*P_l-1,1:2*P_l-1,:,:,P_l,a_l) &
                            ,P_l &
                            ,CCP(1:2*P_l-1,:,P_l,a_l,u_l),v_l,u_l & 
                            ,V_fct_ts(1:2*P_l-1,:,P_l,a_l,u_l)) 
        social=0
        call policy_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,v_l,u_l)&
                            ,F(1:2*P_l-1,1:2*P_l-1,:,:,P_l,a_l) &
                            ,P_l &
                            ,CCP(1:2*P_l-1,:,P_l,a_l,u_l),CCP2(1:2*P_l-1,:,P_l,a_l,u_l),v_l,u_l &
                            ,V_fct(1:2*P_l-1,:,P_l,a_l,u_l),a_l)
        social=1
        call policy_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,v_l,u_l)&
                            ,F(1:2*P_l-1,1:2*P_l-1,:,:,P_l,a_l) & 
                            ,P_l &
                            ,CCP(1:2*P_l-1,:,P_l,a_l,u_l),CCP2(1:2*P_l-1,:,P_l,a_l,u_l),v_l,u_l & 
                            ,V_social(1:2*P_l-1,:,P_l,a_l,u_l),a_l) 
    end do;end do;end do
    

    
    P_l2=P_max
    dist=0.0
    dist=sum((pr_n_u-pr_n_u_old)**2.0d0)
    pr_n_u_old=pr_n_u
    mean_N_old=mean_N

    if (v_l==1) then
        print*,dist,mean_N
    end if

    
    !print*,'press any key to continue'
    !read*,pause_k
    it=it+1
    if ((dist>1.0d-3.and. it<10) .or. it==1) then !1.0d-4 
        go to 1 
    end if
    
    if (it==10 .and. dist>1.0d-2) then
        print*,'equilibrium not reached in village',v_l,dist
    end if
    

    end subroutine
    