subroutine compute_eq_F_CCP(params,F,CCP,V_fct,V_social,n_initial,v_l,mean_N,social_output,private_output,joint_pr,pr_d_a_n,pr_N_n,pr_n)
    use cadastral_maps; use dimensions; use primitives
    implicit none
    double precision,dimension(par),intent(in)::params
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,unobs_types),intent(out)::F
    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types),intent(inout)::CCP
    integer,dimension(plots_in_map,1),intent(inout)::n_initial
    integer,dimension(plots_in_map,1)::n_start
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types),intent(inout)::V_fct,V_social
    integer,intent(in)::v_l    
    double precision,intent(out)::mean_N,social_output,private_output
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types),intent(inout)::joint_pr
    double precision,dimension(3,unobs_types)::pr_n_u,pr_n_u_old
    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types)::CCP_old,CCP2
    double precision,dimension(villages)::village_fe
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::Ef_v !Ef_v: expected productivity Ef_v(1,3,1,4,2,:)
    double precision::dist,Mean_N_old
    integer::p_l,a_l,n_l,P_l2,ind,counter_all,counter_bad,u_l,it
    integer(8),dimension(2*P_max-1,3,3,P_max,unobs_types)::iterations
    double precision,dimension(types_a,2),intent(out)::pr_d_a_n
    double precision,dimension(2*P_max-1),intent(out)::pr_N_n
    double precision,dimension(3),intent(out)::pr_n
    character::pause_k
    
    Mean_N_old=0.0d0
    rho=params(3)
    village_fe(1:villages)=-9.0d0
    
    !Compute expected productivity 
    do a_l=1,types_a;do u_l=1,unobs_types
        call expected_productivity((/params(2),params(1)/),area(a_l),Ef_v(:,:,:,a_l,v_l,u_l),v_l,u_l) !Ef_v(1,:,1,1,1,:)
        !if ( a_l==2 .and. u_l==3) then
        !    print*,'Type',u_l,a_l
        !    print*, 'private return',(ef_v(4,2,8,a_l,v_l,u_l))/(1.0d0-beta*(1.0d0-pi_f_v(1,2,1,v_l,u_l)))-c_s-c_d/(1-PI_s(1,v_l))
        !    print*, 'social return',(ef_v(1,2,1,a_l,v_l,u_l)-c_e),(ef_v(1,2,1,a_l,v_l,u_l)-c_e)/(1.0d0-beta*(1.0d0-pi_f_v(1,2,1,v_l,u_l)))-c_s-c_d/(1-PI_s(1,v_l))
        !end if
    end do;end do

    !Generate beliefs consitent with CCP
!    F=1.0d0
    

!   print*,'generating beliefs'
    it=0
    pr_n_u_old=0.0d0
    n_start=n_initial
    1 n_initial=1
    call generate_beliefs(CCP,V_fct,V_social,Ef_v(:,:,:,:,v_l,:),n_initial,F,v_l,iterations,mean_N,social_output,private_output,joint_pr,pr_d_a_n,pr_N_n) 
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
        call value_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,v_l,u_l)& !Ef_v(1:2*P_l-1,2,P_l,a_l,v_l,:)
                            ,F(1:2*P_l-1,1:2*P_l-1,:,:,P_l,u_l) &
                            ,P_l &
                            ,CCP(1:2*P_l-1,:,P_l,a_l,u_l),v_l,u_l & !CCP(1,:,1,2,:)
                            ,V_fct(1:2*P_l-1,:,P_l,a_l,u_l)) !V_fct(1:2*P_l-1,:,P_l,a_l,:)
        social=0
        call policy_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,v_l,u_l)&
                            ,F(1:2*P_l-1,1:2*P_l-1,:,:,P_l,u_l) &
                            ,P_l &
                            ,CCP(1:2*P_l-1,:,P_l,a_l,u_l),CCP2(1:2*P_l-1,:,P_l,a_l,u_l),v_l,u_l &
                            ,V_fct(1:2*P_l-1,:,P_l,a_l,u_l),a_l)
        social=1
        call policy_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,v_l,u_l)&
                            ,F(1:2*P_l-1,1:2*P_l-1,:,:,P_l,u_l) &
                            ,P_l &
                            ,CCP(1:2*P_l-1,:,P_l,a_l,u_l),CCP2(1:2*P_l-1,:,P_l,a_l,u_l),v_l,u_l & 
                            ,V_social(1:2*P_l-1,:,P_l,a_l,u_l),a_l) !V_social(1,1,1,4,:) CCP(1,:,1,4,3) V_social(1,1,1,4,:)
    end do;
    !if (minval(CCP(1:2*P_l-1,:,P_l,2,u_l)-CCP(1:2*P_l-1,:,P_l,1,u_l))<0.0d0) then
    !    print*,'ufff',minval(CCP(1:2*P_l-1,:,P_l,2,u_l)-CCP(1:2*P_l-1,:,P_l,1,u_l))
    !end if
    end do;end do
    
    
!CCP(1:2*6-1,:,6,2,:)
    !CCP(1,1,8,:,1)
   
    !V_fct(1,3,1,:,4)
    !P_l2=P_max
    !do u_l=1,unobs_types;do a_l=1,types_a; do n_l=1,1;do ind=1,2*P_l2-1; 
    !    if (abs(CCP_old(ind,n_l,P_l2,a_l,u_l)-CCP(ind,n_l,P_l2,a_l,u_l))>0.02) then
    !        print*,'U_T',u_l,'N',n_l-2+ind,' ;A',a_l,'n_l',n_l,'ind',ind
    !        print '(F8.4,F8.4,I12,I12)',real(CCP(ind,n_l,P_l2,a_l,u_l)),real(CCP_old(ind,n_l,P_l2,a_l,u_l)),iterations(ind,n_l,1,P_l2),iterations(ind,n_l,2,P_l2)
    !        print*,''
    !    end if
    !end do; end do;end do;end do
    
    P_l2=P_max
    dist=0.0
    dist=sum((pr_n_u-pr_n_u_old)**2.0d0)
    pr_n_u_old=pr_n_u
    mean_N_old=mean_N

    if (v_l==1) then
        print*,dist
    end if
    !do P_l=2,P_max; do n_l=1,2;do ind=1,2*P_l-1; 
    !    dist=dist+dble(sum(iterations(ind,n_l,1:3,P_l,:)))/dble(sum(iterations(:,1:2,1:3,:,:)))*sum(abs(CCP_old(ind,n_l,P_l,:,:)-CCP(ind,n_l,P_l,:,:)))/dble(types_a)/dble(unobs_types)
    !end do;end do; end do
    !print*,'village',v_l
    !print*,'dist CCP',dist
    
    
    !print*,'npv',mean_NPV,'mean wells',mean_N
    
    !print*,'press any key to continue'
    !read*,pause_k
    it=it+1
    if ((dist>1.0d-3.and. it<10) .or. it==1) then !1.0d-4 
        go to 1 
    end if
    
    if (it==10 .and. dist>1.0d-3) then
        print*,'equilibrium not reached in village',v_l,dist
    end if
    
    do n_l=1,3
        pr_n(n_l)=sum(pr_n_u(n_l,:)*pr_unobs_t)
    end do
    
    
    !call generate_beliefs(CCP_mid,V_fct,Ef_v(:,:,:,:,v_l,:),n_initial,F,v_l,iterations,mean_N,social_output,private_output,Pr_u_X)
    !print*,'dist CCP',dist,'social_output',social_output
    end subroutine
    