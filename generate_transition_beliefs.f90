!subroutine generate_transition_beliefs(T_path,Sims2,CCP,Ef_v,n_ini,F_new,v_l,V_fct,iterations,mean_N,social_output,ccp_mean,diss_N)
!    use cadastral_maps; use primitives
!    implicit none
!    integer,intent(in)::T_path,Sims2
!    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types,T_path),intent(in)::CCP
!    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types),intent(in)::Ef_v 
!    integer,dimension(plots_in_map,Sims2),intent(inout)::n_ini
!    integer,dimension(plots_in_map,T_path+1)::n_initial
!    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,T_path),intent(out)::F_new
!    integer,intent(in)::v_l
!    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types,T_path),intent(in)::V_fct
!    integer(8),dimension(2*P_max-1,3,3,P_max,T_path),intent(out)::iterations
!    integer,dimension(plots_in_map,3)::state,state_old
!    integer::i_l,j_l,t_l,ind,N_all,n_l,P,A,P_l,n_l2,it,m_l,it_min,s_l
!    double precision::u_d,u_s,u_f,u_m,it2
!    double precision,dimension(Sims2,T_path)::NPV,total_N,CCP_av,no_N
!    double precision,dimension(T_path),intent(out)::mean_N,social_output,ccp_mean,diss_N
!    integer(8),dimension(2*P_max-1,2*P_max-1,3,3,P_max)::beliefs_c
!    integer,dimension(1)::seed=123,seed2
!    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max)::F
!    double precision,dimension(2*P_max-1)::CCP_aux
!    character::continue_k
!    
!    
!    
!    CCP_aux=1.0d0/(1.0d0+exp(-(-PI_s_v(1:2*P_max-1,2,P_max,v_l)*c_s-(1.0d0-PI_s_v(1:2*P_max-1,2,P_max,v_l))*c_d)/rho(2)))
!    
!    print*,'got into trans beliefs'
!    !Call seed number
!    call random_seed(GET=seed2)
!    call random_seed(PUT=seed)
!    
!    beliefs_c=0
!    iterations=0
!    F_new=-9.0d0
!    it=0
!    
!    !Store the state for each plot and simulate decision to drill
!    NPV=0.0d0
!    CCP_av=0.0d0
!    no_N=0.0d0
!    
!    !CCP_av(1,:)
!    do s_l=1,Sims2
!    do t_l=1,T_path+1;  
!        it2=0.0d0
!        if (t_l==1) then
!            n_initial(:,t_l)=n_ini(:,s_l)
!        end if
!        !print*,t_l
!        !simulate monsoon next period
!        call RANDOM_NUMBER(u_m)
!        if (u_m<PI_m(1,v_l))then
!            m_l=1
!        else
!            m_l=2
!        end if
!        beliefs_c=0
!
!        state(:,1)=n_initial(:,t_l)
!        do i_l=1,plots_v(v_l)
!            if (active_plots(i_l,v_l,s_l)==1) then
!                N_all=1 !Indicates the number of wells in the adjacency
!                !Loop over all neighbors
!                do j_l=1,PA_type(i_l,1,v_l) !PA_type(i_l,1) stores the number of plots in the adjacency
!                    if (state(neighbors(i_l,j_l,v_l),1)==2)  then
!                        N_all=N_all+1 !number of wells (there is one well)
!                    elseif (state(neighbors(i_l,j_l,v_l),1)==3)  then
!                        N_all=N_all+2 !number of wells (there is two wells)
!                    end if
!                end do
!                state(i_l,2)=N_all !second column in state: number of plots with one well
!                n_l=state(i_l,1) !number of well in reference plot
!
!                P=PA_type(i_l,1,v_l) !number of plots in the adjacency
!                A=PA_type(i_l,2,v_l) !area of the reference plot
!            
!                !Locate position in the state space wrt to the CCP, PI_s_v and ,PI_f_v
!                if (n_l==1) then
!                    ind=N_all 
!                elseif (n_l==2) then
!                    ind=N_all-1
!                elseif (n_l==3) then
!                    ind=N_all-2
!                else
!                    print*,'error generating beliefs'
!                end if         
!                state(i_l,3)=ind
!                !Count transitions (in the first iteration state_old is undefined: no problem)
!                if (t_l>1) then
!                    beliefs_c(state_old(i_l,3),state(i_l,3),state_old(i_l,1),state(i_l,1),P)=&
!                    beliefs_c(state_old(i_l,3),state(i_l,3),state_old(i_l,1),state(i_l,1),P)+1
!                end if
!                
!                if (t_l<=T_path) then
!                    !Compute NPV
!                    if (n_l==1 .or. n_l==2)then  
!                        it2=it2+1.0d0
!                        CCP_av(s_l,t_l)=(it2-1.0d0)/it2*CCP_av(s_l,t_l)+1.0d0/it2*CCP(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l),t_l)
!                        NPV(s_l,t_l)=dble(i_l-1)/dble(i_l)*NPV(s_l,t_l)+1.0d0/dble(i_l)*(Ef_v(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l))- &
!                                                CCP(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l),t_l)*(PI_s_v(ind,n_l,P,v_l)*c_s+(1.0d0-PI_s_v(ind,n_l,P,v_l))*c_d)-c_e*dble(n_l-1))
!                    else
!                        NPV(s_l,t_l)=dble(i_l-1)/dble(i_l)*NPV(s_l,t_l)+1.0d0/dble(i_l)*(Ef_v(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l))-c_e*dble(n_l-1))
!                    end if
!                    
!                    !!Decision to dismantle well
!
!                    if (n_l==2 .and. V_fct(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l),t_l)<V_fct(ind,1,P,A,unobs_types_i(i_l,v_l,s_l),t_l) .and. t_l==21) then
!                        n_initial(i_l,t_l+1)=1
!                        if (t_l<=T_path) then
!                            no_N(s_l,t_l)=no_N(s_l,t_l)+1.0d0
!                        end if
!                    elseif (n_l==3 .and. V_fct(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l),t_l)<V_fct(ind,2,P,A,unobs_types_i(i_l,v_l,s_l),t_l) .and. t_l==21) then
!                        n_initial(i_l,t_l+1)=2
!                        if (t_l<=T_path) then
!                            no_N(s_l,t_l)=no_N(s_l,t_l)+1.0d0
!                        end if
!                        if (V_fct(ind,2,P,A,unobs_types_i(i_l,v_l,1),t_l)<V_fct(ind,1,P,A,unobs_types_i(i_l,v_l,1),t_l)) then
!                            n_initial(i_l,t_l+1)=1
!                            if (t_l<=T_path) then
!                                no_N(s_l,t_l)=no_N(s_l,t_l)+1.0d0
!                            end if
!                        end if
!                        
!                    !Well drilling decision and failures/successes
!                    elseif (n_l==1) then !no well
!                        call RANDOM_NUMBER(u_d)
!                        if (u_d<CCP(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l),t_l)) then !decides to drill
!                            call RANDOM_NUMBER(u_s)
!                            if (u_s<PI_s_v(ind,n_l,P,v_l)) then !successful attempt
!                                n_initial(i_l,t_l+1)=n_l+1
!                            else !unsuccessful attempt
!                                n_initial(i_l,t_l+1)=n_l
!                            end if
!                        else !decides not to drill
!                            n_initial(i_l,t_l+1)=n_l
!                        end if
!                    elseif (n_l==2) then !one well
!                        call RANDOM_NUMBER(u_d)
!                        if (u_d<CCP(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l),t_l)) then !decides to drill
!                            call RANDOM_NUMBER(u_s)
!                            if (u_s<PI_s_v(ind,n_l,P,v_l)) then !successful attempt
!                                call RANDOM_NUMBER(u_f)
!                                if (u_f<PI_fm(N_all-1,m_l,unobs_types_i(i_l,v_l,s_l),v_l)) then !failure of the previous well
!                                    n_initial(i_l,t_l+1)=n_l
!                                else
!                                    n_initial(i_l,t_l+1)=n_l+1
!                                end if
!                            else !unsuccessful attempt
!                                call RANDOM_NUMBER(u_f)
!                                if (u_f<PI_fm(N_all-1,m_l,unobs_types_i(i_l,v_l,s_l),v_l)) then !failure of the previous well
!                                    n_initial(i_l,t_l+1)=n_l-1
!                                else
!                                    n_initial(i_l,t_l+1)=n_l
!                                end if
!                            end if
!                        else !decides not to drill
!                            call RANDOM_NUMBER(u_f)
!                            if (u_f<PI_fm(N_all-1,m_l,unobs_types_i(i_l,v_l,s_l),v_l)) then !failure of the previous well
!                                n_initial(i_l,t_l+1)=n_l-1
!                            else
!                                n_initial(i_l,t_l+1)=n_l
!                            end if 
!                        end if 
!                    elseif(n_l==3) then !two wells
!                        call RANDOM_NUMBER(u_f)
!                        if (u_f<PI_fm(N_all-1,m_l,unobs_types_i(i_l,v_l,s_l),v_l)**2) then !failure of the two wells
!                            n_initial(i_l,t_l+1)=n_l-2
!                        elseif (u_f<PI_fm(N_all-1,m_l,unobs_types_i(i_l,v_l,s_l),v_l)**2+(1.0d0-PI_fm(N_all-1,m_l,unobs_types_i(i_l,v_l,s_l),v_l))**2) then !failure of none
!                            n_initial(i_l,t_l+1)=n_l
!                        else !failure of one
!                            n_initial(i_l,t_l+1)=n_l-1
!                        end if 
!                    else
!                        print*,'error in gen beliefs 2'
!                    end if
!                end if
!            else
!                NPV(s_l,t_l)=dble(i_l-1)/dble(i_l)*NPV(s_l,t_l)+1.0d0/dble(i_l)*0.0d0
!            end if
!        end do
!        !Store current state
!        state_old=state
!        !Compute beliefs
!        if (t_l>1) then
!            F=0.0d0
!            do P_l=2,P_max; do ind=1,2*P_l-1; do n_l=1,3; do n_l2=1,3
!                    if (n_l==1 .and. n_l2==3) then
!                        F(ind,:,n_l,n_l2,P_l)=-9.0d0
!                    elseif ((sum(beliefs_c(ind,:,n_l,n_l2,P_l)))==0) then
!                        F(ind,1:2*P_l-1,n_l,n_l2,P_l)=-9.0d0
!                    else
!                        F(ind,:,n_l,n_l2,P_l)=dble(beliefs_c(ind,:,n_l,n_l2,P_l))/dble(sum(beliefs_c(ind,:,n_l,n_l2,P_l)))
!                        iterations(ind,n_l,n_l2,P_l,t_l-1)=iterations(ind,n_l,n_l2,P_l,t_l-1)+sum(beliefs_c(ind,:,n_l,n_l2,P_l))
!                        F_new(ind,:,n_l,n_l2,P_l,t_l-1)=dble(sum(beliefs_c(ind,:,n_l,n_l2,P_l)))/dble(iterations(ind,n_l,n_l2,P_l,t_l-1))*F(ind,:,n_l,n_l2,P_l) +&
!                                                    dble(iterations(ind,n_l,n_l2,P_l,t_l-1)-sum(beliefs_c(ind,:,n_l,n_l2,P_l)))/dble(iterations(ind,n_l,n_l2,P_l,t_l-1))*F_new(ind,:,n_l,n_l2,P_l,t_l-1) 
!                        if (minval(F_new(ind,:,n_l,n_l2,P_l,t_l-1))<0)then
!                            print*,''
!                        end if
!                    end if
!            end do;end do;end do;end do
!        end if
!        if (t_l<=T_path) then
!            total_N(s_l,t_l)=dble(sum(n_initial(1:plots_v(v_l),t_l))-plots_v(v_l))
!        end if
!        it=0
!        it=it+1
!    end do
!    n_ini(:,s_l)=n_initial(:,19)
!    end do
!
!    !close(13)
!    
!
!    ! In case I don't have observations for a given state, I consider that the transition pr 
!    ! is the same for all possible future states
!    do t_l=1,T_path;do P_l=2,P_max; do ind=1,2*P_l-1; do n_l=1,3; do n_l2=1,min(n_l+1,3)
!        it_min=300000
!        if (iterations(ind,n_l,n_l2,P_l,t_l)==0) then
!            F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,t_l)=1.0d0/(2.0d0*dble(P_l)-1.0d0)
!        end if
!        if (isnan(F_new(ind,1,n_l,n_l2,P_l,t_l))) then
!            print*,'error in generate beliefs'
!        end if
!        F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,t_l)=F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,t_l)/sum(F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,t_l))
!        !print*,sum(F_new(ind,:,n_l,n_l2,P_l))
!    end do;end do;end do;end do;end do
!    
!    !Beliefs for plots with no neighbors are degenerate: they know what the future will look like cond on their own outcomes
!    P_l=1
!    ind=1
!    do n_l=1,3; do n_l2=1,min(n_l+1,3)
!        F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,:)=1.0d0
!    end do;end do
!    
!    social_output=sum(NPV,1)/dble(Sims2)/mean_area(v_l)
!    mean_N=sum(total_N,1)/dble(Sims2)/dble(plots_v(v_l))
!    diss_N=sum(no_N,1)/sum(total_N,1)
!    ccp_mean=sum(CCP_av,1)/dble(Sims2)
!    !print*,'av drilling',ccp_mean
!    
!    call random_seed(PUT=seed2)
!    
!
!!print*,'press key to continue'    
!!read*,continue_k
!end subroutine
!    