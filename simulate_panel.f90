subroutine simulate_panel(CCP,n_initial,seed_new)
    use cadastral_maps; use primitives; use simulation
    implicit none
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types),intent(in)::CCP
    integer,dimension(plots_in_map,villages),intent(inout)::n_initial
    integer(8),dimension(plots_in_map,3)::state,state_old
    integer(8),dimension(plots_in_map)::drill_old
    integer(8)::i_l,j_l,t_l,ind,N_all,n_l,P,A,P_l,n_l2,it,m_l,it_min,a_l,u_l,ind2,N_all2,i,v_l
    double precision::u_d,u_s,u_f,u_m,it2
    integer(8),dimension(1)::seed=321,seed2
    integer(8),parameter::burn_t=1000

    character::continue_k
    integer(8),dimension(1),intent(inout):: seed_new
    
    !seed_new=987654321
    
    it=0

    v_l=1
    V_type=1
    
    do t_l=1,T_sim;   
        !print*,t_l
        !simulate monsoon next period
        call random_value( seed_new, u_m )
        if (u_m<PI_m(1,v_l))then
            m_l=1
        else
            m_l=2
        end if

        state(:,1)=n_initial(:,v_l)
        !print*,'t_l',t_l,'av number of wells per plot',real(sum(n_initial(1:plots_v(v_l),1))-plots_v(v_l))/real(plots_v(v_l))
        do i_l=1,plots_v(v_l)
            if (active_plots(i_l,v_l,1)==1) then
                N_all=1 !Indicates the number of wells in the adjacency
                !Loop over all neighbors
                do j_l=1,PA_type(i_l,1,v_l) !PA_type(i_l,1) stores the number of plots in the adjacency
                    if (state(neighbors(i_l,j_l,v_l),1)==2)  then !neighbors(42,:,v_l)
                        N_all=N_all+1 !number of wells (there is one well)
                    elseif (state(neighbors(i_l,j_l,v_l),1)==3)  then
                        N_all=N_all+2 !number of wells (there is two wells)
                    end if
                end do
                state(i_l,2)=N_all !second column in state: number of plots with one well
                n_l=state(i_l,1) !number of well in reference plot

                P=PA_type(i_l,1,v_l) !number of plots in the adjacency
                A=PA_type(i_l,2,v_l) !area of the reference plot
                
                if (i_l<=plots_i) then
                    P_type(i_l)=P
                    A_type(i_l)=A
                    n_data(t_l,i_l)=n_initial(i_l,v_l)
                    Pr_N_data(:,t_l,i_l)=0
                    Pr_N_data(min(N_all,max_NFW+1),t_l,i_l)=1.0d0
                end if
            
                !Locate position in the state space wrt to the CCP, PI_s_v and ,PI_f_v
                if (n_l==1) then
                    ind=N_all 
                elseif (n_l==2) then
                    ind=N_all-1
                elseif (n_l==3) then
                    ind=N_all-2
                else
                    print*,'error generating beliefs'
                end if 
                if (ind==0) then
                    print*,'paused'
                    print'(<P_max+2>I6)',i_l,PA_type(i_l,1,v_l),neighbors(i_l,:,v_l)
                    read*,continue_k
                end if
                state(i_l,3)=ind

                !Well drilling decision and failures/successes
                drill_old(i_l)=0
                if (n_l==1) then !no well
                    call random_value( seed_new, u_d )
                    if (u_d<CCP(ind,n_l,P,A,v_l,unobs_types_i(i_l,v_l,1))) then !decides to drill unobs_types_i(1:1052,1)
                        drill_old(i_l)=1
                        call random_value( seed_new, u_s )
                        if (u_s<PI_s_v(ind,n_l,P,v_l)) then !successful attempt
                            n_initial(i_l,v_l)=n_l+1
                        else !unsuccessful attempt
                            n_initial(i_l,v_l)=n_l
                        end if
                    else !decides not to drill
                        n_initial(i_l,v_l)=n_l
                    end if
                elseif (n_l==2) then !one well
                    call random_value( seed_new, u_d )
                    if (u_d<CCP(ind,n_l,P,A,v_l,unobs_types_i(i_l,v_l,1))) then !decides to drill
                        drill_old(i_l)=1
                        call random_value( seed_new, u_s )
                        if (u_s<PI_s_v(ind,n_l,P,v_l)) then !successful attempt
                            call random_value( seed_new, u_f )
                            if (u_f<PI_fm(N_all-1,m_l,unobs_types_i(i_l,v_l,1),v_l)) then !failure of the previous well
                                n_initial(i_l,v_l)=n_l
                            else
                                n_initial(i_l,v_l)=n_l+1
                            end if
                        else !unsuccessful attempt
                            call random_value( seed_new, u_f )
                            if (u_f<PI_fm(N_all-1,m_l,unobs_types_i(i_l,v_l,1),v_l)) then !failure of the previous well PI_fm(:)
                                n_initial(i_l,v_l)=n_l-1
                            else
                                n_initial(i_l,v_l)=n_l
                            end if
                        end if
                    else !decides not to drill
                        call random_value( seed_new, u_f )
                        if (u_f<PI_fm(N_all-1,m_l,unobs_types_i(i_l,v_l,1),v_l)) then !failure of the previous well
                            n_initial(i_l,v_l)=n_l-1
                        else
                            n_initial(i_l,v_l)=n_l
                        end if 
                    end if 
                elseif(n_l==3) then !two wells
                    call random_value( seed_new, u_f )
                    if (u_f<PI_fm(N_all-1,m_l,unobs_types_i(i_l,v_l,1),v_l)**2) then !failure of the two wells
                        n_initial(i_l,v_l)=n_l-2
                    elseif (u_f<PI_fm(N_all-1,m_l,unobs_types_i(i_l,v_l,1),v_l)**2+(1.0d0-PI_fm(N_all-1,m_l,unobs_types_i(i_l,v_l,1),v_l))**2) then !failure of none
                        n_initial(i_l,v_l)=n_l
                    else !failure of one
                        n_initial(i_l,v_l)=n_l-1
                    end if 
                else
                    print*,'error in gen beliefs 2'
                end if
                
                if (i_l<=plots_i) then
                    drilling_it(t_l,i_l,1)=drill_old(i_l)
                end if
            end if  
            
        end do
        !Store current state
        state_old=state
    end do
    
    !Modal number of functionning wells
    do i_l=1,plots_i;do t_l=1,T_sim
        modal_N(t_l,i_l)=maxloc(Pr_N_data(:,t_l,i_l),1)-(n_data(t_l,i_l)-1) 
    end do;end do
    N_bar=dble(sum(modal_N,1))/dble(T_sim)

end subroutine

