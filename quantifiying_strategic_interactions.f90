!subroutine quantifiying_strategic_interactions(params)
!    use cadastral_maps; use dimensions; use primitives
!    implicit none
!    double precision,dimension(par),intent(in)::params
!    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,unobs_types)::F
!    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types)::CCP_mid
!    integer,dimension(plots_in_map,1)::n_initial
!    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types)::V_fct,V_social
!    integer::v_l    
!    double precision::mean_N,social_output,private_output
!    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types)::Pr_u_X
!    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types)::CCP_old,CCP,CCP2
!    
!    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::Ef_v !Ef_v: expected productivity
!    double precision::dist
!    integer::p_l,a_l,n_l,n_l2,P_l2,ind,counter_all,counter_bad,u_l
!    integer(8),dimension(2*P_max-1,3,3,P_max,unobs_types)::iterations
!    character::pause_k
!    rho=params(3)
!    v_l=1
!    !Compute expected productivity 
!    do a_l=1,types_a;do u_l=1,unobs_types
!        call expected_productivity(params(1:3),area(a_l),Ef_v(:,:,:,a_l,v_l,u_l),v_l,u_l)
!        print*,'Type',u_l,a_l
!        print*, 'private return',(Ef_v(1,2,1,a_l,v_l,u_l))/(1.0d0-beta*(1.0d0-PI_f_v(1,2,1,v_l,u_l)))-c_s
!        print*, 'social return',(Ef_v(1,2,1,a_l,v_l,u_l)-c_e),(Ef_v(1,2,1,a_l,v_l,u_l)-c_e)/(1.0d0-beta*(1.0d0-PI_f_v(1,2,1,v_l,u_l)))-c_s
!    end do;end do
!
!    CCP_mid=0.07d0
!    n_initial(:,1)=1
!    V_fct=0.0d0 
!    V_social=0.0d0
!    
!    !Generate beliefs consitent with CCP
!    F=1.0d0
!    
!
!!   print*,'generating beliefs'
!    n_initial=1
!    F=0.0d0
!    do P_l=1,P_max; do ind=1,2*P_l-1; do n_l=1,3; do n_l2=1,min(n_l+1,3); do u_l=1,unobs_types
!        F(ind,1,n_l,n_l2,P_l,u_l)=1.0d0
!    end do;end do;end do;end do;end do
!    tau=0.0d0
!    do P_l=1,P_max; do u_l=1,unobs_types; do a_l=1,types_a !do P_l=1,P_max; do u_l=1,unobs_types; do a_l=1,types_a
!        call value_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,v_l,u_l)&
!                            ,F(1:2*P_l-1,1:2*P_l-1,:,:,P_l,u_l) &
!                            ,P_l &
!                            ,CCP(1:2*P_l-1,:,P_l,a_l,u_l),v_l,u_l &
!                            ,V_fct(1:2*P_l-1,:,P_l,a_l,u_l))
!        call policy_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,v_l,u_l)&
!                            ,F(1:2*P_l-1,1:2*P_l-1,:,:,P_l,u_l) &
!                            ,P_l &
!                            ,CCP(1:2*P_l-1,:,P_l,a_l,u_l),CCP2(1:2*P_l-1,:,P_l,a_l,u_l),v_l,u_l &
!                            ,V_social(1:2*P_l-1,:,P_l,a_l,u_l),a_l)
!    end do; end do;end do
!    
!    !print*,'EDP value',V_social(1,1,1,2,3)
!    CCP_mid=CCP
!
!1    call generate_beliefs(CCP_mid,V_fct,V_social,Ef_v(:,:,:,:,v_l,:),n_initial,F,v_l,iterations,mean_N,social_output,private_output,Pr_u_X)
!    read*,pause_k
!    !For each plot type obtain a new CCP given beliefs
!    !print*,'policy step'
!    CCP_old=CCP
!    dist=0.0d0
!    counter_bad=0
!    counter_all=0
!    do P_l=1,P_max; do u_l=1,unobs_types; do a_l=1,types_a
!        call value_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,v_l,u_l)&
!                            ,F(1:2*P_l-1,1:2*P_l-1,:,:,P_l,u_l) &
!                            ,P_l &
!                            ,CCP(1:2*P_l-1,:,P_l,a_l,u_l),v_l,u_l &
!                            ,V_fct(1:2*P_l-1,:,P_l,a_l,u_l))
!    end do; end do;end do
!
!    !V_fct(1,3,2:7,1,1)
!   
!
!    P_l2=P_max
!    dist=0.0
!    do P_l=2,P_max; do n_l=1,2;do ind=1,2*P_l-1; 
!        dist=dist+dble(sum(iterations(ind,n_l,1:3,P_l,:)))/dble(sum(iterations(:,1:2,1:3,:,:)))*sum(abs(CCP_old(ind,n_l,P_l,:,:)-CCP(ind,n_l,P_l,:,:))/CCP_old(ind,n_l,P_l,:,:))/dble(types_a)/dble(unobs_types)
!    end do;end do; end do
!    !print*,'village',v_l
!    print*,'dist CCP',dist,'social_output',social_output
!    
!    !New guess of the ccp is half way through
!    CCP_mid=CCP*0.5d0+CCP_old*0.5d0
!    
!    !print*,'npv',mean_NPV,'mean wells',mean_N
!    
!    !print*,'press any key to continue'
!    !read*,pause_k
!    if (dist>6.0d-3) then !1.0d-4 
!        go to 1 
!    end if
!    
!    !call generate_beliefs(CCP_mid,V_fct,Ef_v(:,:,:,:,v_l,:),n_initial,F,v_l,iterations,mean_N,social_output,private_output,Pr_u_X)
!    !print*,'dist CCP',dist,'social_output',social_output
!end subroutine