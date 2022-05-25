!subroutine counterfactual_1(params_MLE)
!    use dimensions; use cadastral_maps; use simulation; use primitives
!    implicit none
!    double precision,dimension(par)::params_true,params_MLE
!    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,villages)::F_true
!    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP_true
!    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types)::V_fct
!    integer,dimension(plots_in_map,villages)::n_dist
!    double precision,dimension(villages)::mean_budget
!    integer,parameter::samples=10
!    double precision,dimension(samples)::mean_N,mean_NPV
!    integer::v_l,p_l,it,ns
!    character::end_key
!    integer,parameter::nkk=8
!    double precision,dimension(nkk)::fraction_grid
!
!    
!    tau=0.0d0
!    
!    fraction_grid(1)=0.95
!    do p_l=2,nkk
!        fraction_grid(p_l)=fraction_grid(p_l-1)-0.1d0
!    end do
!        
!    !I want to look for the optimal NPV by selecting a given number of farmers
!    do v_l=1,1!villages
!        print*,'village,',v_l 
!        !Run benchmark and compute NPV
!        CCP_true(:,:,:,:,v_l,:)=0.07d0
!        n_dist(:,v_l)=1
!        V_fct=0.0d0
!        T_g=0.0d0
!        call compute_eq_F_CCP(params_MLE,F_true(:,:,:,:,:,v_l),CCP_true(:,:,:,:,v_l,:),V_fct,n_dist(:,v_l),v_l,mean_N(1),mean_NPV(1),mean_budget(v_l),Pr_u_X(:,:,:,:,v_l,:))
!        print*,'mean_NPV',mean_NPV(v_l)
!        OPEN(UNIT=12, FILE=path_results//"counterfactuals_1.txt")
!        write(12,'(F20.3,F20.3,F20.3)'),0.0,mean_N(1),mean_NPV(1)
!        close(12)
!        
!        do p_l=1,nkk
!            do ns=1,samples
!                print*,'p=',p_l,' ns=',ns
!                call select_active_plots(fraction_grid(p_l),v_l)
!                call compute_eq_F_CCP(params_MLE,F_true(:,:,:,:,:,v_l),CCP_true(:,:,:,:,v_l,:),V_fct,n_dist(:,v_l),v_l,mean_N(ns),mean_NPV(ns),mean_budget(v_l),Pr_u_X(:,:,:,:,v_l,:))
!                print*,'mean_NPV: ',mean_NPV(ns)
!            end do
!            OPEN(UNIT=12, FILE=path_results//"counterfactuals_1.txt",access='append')
!            write(12,'(F20.3,F20.3,F20.3)'),fraction_grid(p_l),mean_N(maxloc(mean_NPV)),maxval(mean_NPV)
!            close(12)
!        end do
!       
!    end do
!    
!    !Compute NPV and mean number of wells in counterfactual
!
!    
!end subroutine
!    
!
!
!
!subroutine select_active_plots(fraction,v_l)
!use cadastral_maps
!    implicit none
!    double precision,intent(in):: fraction
!    integer,intent(in)::v_l
!    integer::plots_active,plot_count,ind_p,i,j,ind
!    double precision::u
!    
!    !Reload neighbor info of cadastral maps
!    call load_cadastral_maps()
!    
!    active_plots(1:plots_v(v_l),v_l)=0
!    plots_active=fraction*plots_v(v_l)
!    plot_count=0
!    do while (plot_count<plots_active)
!        !call compute_new_number of neighbors()
!        !call sort_most_productive_plots()
!        call RANDOM_NUMBER(u)
!        ind_p=u*(plots_v(v_l)-1)+1
!        if (active_plots(ind_p,v_l)==0) then
!            active_plots(ind_p,v_l)=1
!            plot_count=plot_count+1
!        end if
!    end do
!        
!    !Store which are my new neighbors
!    do i=1,plots_in_map;
!        neighbors_map(i,i,v_l)=1
!        if (active_plots(i,v_l)==0) then
!            neighbors_map(1:i-1,i,v_l)=0
!            neighbors_map(i+1:plots_v(v_l),i,v_l)=0
!            neighbors_map(i,1:i-1,v_l)=0
!            neighbors_map(i,i+1:plots_v(v_l),v_l)=0
!        end if
!        ind=0
!        do j=1,plots_in_map
!            if (neighbors_map(i,j,v_l)==1 .and. ind<P_max .and. active_plots(ind_p,v_l)==1) then !the number of neighbors cannot be greater than P_max ortherwise I select the firt P_max neighbors
!                ind=ind+1
!                neighbors(i,ind,v_l)=j
!            end if
!        end do
!        !I need to make sure that the reference plot is a neighbor plot for plots whose
!        !number of neighboring plots is larger than P_max in the data so in case it is not, I force the last neighbor to be the reference plot
!        if (ind==P_max .and. ind<i) then
!            neighbors(i,P_max,v_l)=i
!        end if     
!    end do
!    !number of neighbors
!    PA_type(1:plots_v(v_l),1,v_l)=min(sum(neighbors_map(1:plots_v(v_l),1:plots_v(v_l),v_l),2),P_max)   
!    
!    
!    
!end subroutine
!    
!subroutine compute_NPV_SP(params,F,CCP_mid,V_fct,n_initial,v_l,mean_N,mean_NPV,mean_budget)
!    use cadastral_maps; use dimensions; use primitives
!    implicit none
!    double precision,dimension(par),intent(in)::params
!    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max),intent(out)::F
!    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types),intent(inout)::CCP_mid
!    integer,dimension(plots_in_map,1),intent(inout)::n_initial
!    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types),intent(inout)::V_fct
!    integer,intent(in)::v_l    
!    double precision,intent(out)::mean_N,mean_NPV,mean_budget
!    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types)::CCP_old,CCP
!    
!    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::Ef_v !Ef_v: expected productivity
!    double precision::dist
!    integer::a_l,u_l
!    integer(8),dimension(2*P_max-1,3,3,P_max)::iterations
!    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::Pr_u_X
!    character::pause_k
!    
!    !Set scale parameter Gumbel distribution of shocks
!    !rho=params(4)
!
!    !Compute expected productivity 
!    do u_l=1,unobs_types;do a_l=1,types_a
!        call expected_productivity(params(1:2),area(a_l),Ef_v(:,:,:,a_l,v_l,u_l),v_l,u_l)
!    end do;end do
!
!    !Generate beliefs consitent with CCP
!    F=1.0d0
!    CCP=CCP_mid
!!   print*,'generating beliefs'
!    n_initial=1
!    call generate_beliefs(CCP_mid,V_fct,Ef_v(:,:,:,:,v_l,:),n_initial,F,v_l,iterations,mean_N,mean_NPV,mean_budget,Pr_u_X(:,:,:,:,v_l,:))
!end subroutine
!    