subroutine generate_transition_beliefs(CCP_ini,CCP_tr,F_new,n_tr,transition_time,total_time,v_l,V_fct_tr,V_fct_initial,Ef_v,stats_out,av_drilling_rate)
    use cadastral_maps; use primitives; use simulation
    implicit none
    integer,intent(in)::transition_time,total_time,v_l
    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types),intent(in)::CCP_ini
    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types,transition_time),intent(in)::CCP_tr
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types,transition_time),intent(in)::V_fct_tr
    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types)::CCP
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types)::V_fct_t
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,types_a,transition_time),intent(out)::F_new
    double precision,dimension(total_time),intent(out)::n_tr,av_drilling_rate
    double precision,dimension(2),intent(out)::stats_out
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types),intent(in)::V_fct_initial
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types),intent(in)::Ef_v
    
    integer,dimension(2*P_max-1,3,3,P_max,types_a,transition_time)::iterations
    integer(8),dimension(1):: seed_new
    integer::P_l,ind,n_l,n_l2,a_l,m_l,t_l,r_l,i_l,s_l,P,A,j_l,N_all,ss_time,count
    integer(8),dimension(plots_in_map,3)::state,state_old
    integer,dimension(plots_in_map,1)::n_initial
    double precision::u_m,u_s,u_d,u_f,drilling
    double precision,dimension(8)::draw
    integer(8),dimension(2*P_max-1,2*P_max-1,3,3,P_max,types_a,transition_time)::beliefs_c
    integer,dimension(total_time,sims_tr)::total_N
    double precision,dimension(transition_time,sims_tr,2)::av_output_per_acre
    double precision,dimension(total_time,sims_tr)::av_drilling_rate_s
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,types_a)::F
    double precision::dc,y
    character::pause_k
    
    ss_time=total_time-transition_time

    !prior of beliefs
    F_new=-9.0
    iterations=0
    do P_l=1,P_max;do ind=1,2*P_l-1; do n_l=1,3; do n_l2=1,1; do a_l=1,types_a
        F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l,:)=1.0d0/dble(2*P_l-1)
        iterations(ind,n_l,n_l2,P_l,a_l,:)=1
    end do;end do;end do;end do;end do
    beliefs_c=0
    
    !set seed
    seed_new=4321
    
    av_output_per_acre=0.0d0
    av_drilling_rate_s=0.0d0
    
    do s_l=1,sims_tr
        state_old=1
        do t_l=1,total_time
            if (t_l<=ss_time) then
                CCP=CCP_ini    
                V_fct_t=V_fct_initial
            else
                CCP=CCP_tr(:,:,:,:,:,t_l-ss_time)
                V_fct_t=V_fct_tr(:,:,:,:,:,t_l-ss_time)
            end if
            
            !simulate monsoon next period
            call random_value( seed_new, u_m )
            if (u_m<PI_m(1,v_l))then
                m_l=1
            else
                m_l=2
            end if
            if (t_l>1) then
                state(:,1)=n_initial(:,1)
            else
                n_initial(:,1)=1
                state(:,1)=1
            end if
            count=0
            do i_l=1,plots_v(v_l)
                do r_l=1,8
                    call random_value(seed_new,draw(r_l))
                end do
                
                if (active_plots(i_l,v_l,s_l)==1) then  
                    count=count+1
                    N_all=1 !Indicates the number of wells in the adjacency
                    !Loop over all neighbors
                    do j_l=1,PA_type(i_l,1,v_l,s_l) !PA_type(i_l,1) stores the number of plots in the adjacency
                        if (state(neighbors(i_l,j_l,v_l,s_l),1)==2)  then  !neighbors(i_l,:,v_l,s_l) 
                            N_all=N_all+1 !number of wells (there is one well)
                        elseif (state(neighbors(i_l,j_l,v_l,s_l),1)==3)  then !neighbors(1,:,v_l,s_l)
                            N_all=N_all+2 !number of wells (there is two wells)
                        end if
                    end do
                    state(i_l,2)=N_all !second column in state: number of plots with one well
                    n_l=state(i_l,1) !number of well in reference plot

                    P=PA_type(i_l,1,v_l,s_l) !number of plots in the adjacency
                    A=PA_type(i_l,2,v_l,s_l) !area of the reference plot
            
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
                    
                    state(i_l,3)=ind
                    !Count transitions (in the first iteration state_old is undefined: no problem)
                    if (t_l>ss_time) then
                        beliefs_c(state_old(i_l,3),state(i_l,3),state_old(i_l,1),1,P,A,t_l-ss_time)=&
                        beliefs_c(state_old(i_l,3),state(i_l,3),state_old(i_l,1),1,P,A,t_l-ss_time)+1
                    end if 
                    
                    state_old(i_l,:)=state(i_l,:)
                    
                    !I allow well-dismantling
                    !if (t_l==ss_time+1) then
                    !    if (n_l==3 .and. V_fct_tr(ind,3,P,A,unobs_types_i(i_l,v_l,s_l),1)<V_fct_tr(ind,2,P,A,unobs_types_i(i_l,v_l,s_l),1)) then
                    !        n_l=n_l-1
                    !        N_all=N_all-1
                    !        !print*,'well dismantled'
                    !    end if
                    !    if (n_l==2 .and. V_fct_tr(ind,2,P,A,unobs_types_i(i_l,v_l,s_l),1)<V_fct_tr(ind,1,P,A,unobs_types_i(i_l,v_l,s_l),1)) then
                    !        n_l=n_l-1
                    !        N_all=N_all-1
                    !        !print*,'well dismantled'
                    !    end if
                    !end if
                        
                    
                    !Well drilling decision and failures/successes
                    y=Ef_v(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l))
                    dc=0.0d0
                    drilling=0.0d0
                    if (n_l==1) then !no well
                        u_d=draw(1) 
                        if (u_d<CCP(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l))) then !decides to drill
                            drilling=1.0d0
                            u_s=draw(2) 
                            if (u_s<PI_s_v(ind,n_l,P,v_l)) then !successful attempt
                                n_initial(i_l,1)=n_l+1
                                dc=c_s
                            else !unsuccessful attempt
                                n_initial(i_l,1)=n_l
                                dc=c_d
                            end if
                        else !decides not to drill
                            n_initial(i_l,1)=n_l
                        end if
                    elseif (n_l==2) then !one well
                        u_d=draw(3) 
                        if (u_d<CCP(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l))) then !decides to drill
                            drilling=1.0d0
                            u_s=draw(4) 
                            if (u_s<PI_s_v(ind,n_l,P,v_l)) then !successful attempt
                                u_f=draw(5) 
                                dc=c_s
                                if (u_f<PI_fm(N_all-1,m_l,v_l,unobs_types_i(i_l,v_l,s_l))) then !failure of the previous well
                                    n_initial(i_l,1)=n_l
                                else
                                    n_initial(i_l,1)=n_l+1
                                end if
                            else !unsuccessful attempt
                                u_f=draw(6) 
                                dc=c_d
                                if (u_f<PI_fm(N_all-1,m_l,v_l,unobs_types_i(i_l,v_l,s_l))) then !failure of the previous well PI_fm(:)
                                    n_initial(i_l,1)=n_l-1

                                else
                                    n_initial(i_l,1)=n_l
                                end if
                            end if
                        else !decides not to drill
                            u_f=draw(7) 
                            if (u_f<PI_fm(N_all-1,m_l,v_l,unobs_types_i(i_l,v_l,s_l))) then !failure of the previous well
                                n_initial(i_l,1)=n_l-1
                            else
                                n_initial(i_l,1)=n_l
                            end if 
                        end if 
                    elseif(n_l==3) then !two wells
                        u_f=draw(8) 
                        if (u_f<PI_fm(N_all-1,m_l,v_l,unobs_types_i(i_l,v_l,s_l))**2) then !failure of the two wells !PI_fm(:,1,3,1) !PI_fm(:,m_l,v_l,unobs_types_i(i_l,v_l,s_l))
                            n_initial(i_l,1)=n_l-2
                        elseif (u_f<(1.0d0-PI_fm(N_all-1,m_l,v_l,unobs_types_i(i_l,v_l,s_l)))**2+PI_fm(N_all-1,m_l,v_l,unobs_types_i(i_l,v_l,s_l))**2) then !failure of none (this is correct)
                            n_initial(i_l,1)=n_l
                        else !failure of one
                            n_initial(i_l,1)=n_l-1
                        end if 
                    else
                        print*,'error in gen beliefs 2'
                    end if 
                    if (t_l>ss_time) then
                        av_output_per_acre(t_l-ss_time,s_l,1)=dble(count-1)/dble(count)*av_output_per_acre(t_l-ss_time,s_l,1)+1.0d0/dble(count)*beta**dble(t_l-ss_time-1)*(y-dc-dble(n_l-1)*c_e)/area(A) 
                        av_output_per_acre(t_l-ss_time,s_l,2)=dble(count-1)/dble(count)*av_output_per_acre(t_l-ss_time,s_l,2)+1.0d0/dble(count)*beta**dble(t_l-ss_time-1)*(y-dc)/area(A) 
                    end if
                    av_drilling_rate_s(t_l,s_l)=dble(count-1)/dble(count)*av_drilling_rate_s(t_l,s_l)+1.0d0/dble(count)*drilling
                end if
            end do
            
            total_N(t_l,s_l)=sum(n_initial(:,1)-1)
            
            
            !Compute beliefs
            if (t_l>ss_time) then
                F=0.0d0
                do P_l=2,P_max; do ind=1,2*P_l-1; do n_l=1,3; do n_l2=1,1;do a_l=1,types_a
                        if (n_l==1 .and. n_l2==3) then
                            F(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)=-9.0d0
                        elseif(sum(beliefs_c(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l,t_l-ss_time))>0) then !beliefs_c(1,1,2,1,5,4,:)
                            F(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)=dble(beliefs_c(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l,t_l-ss_time))/dble(sum(beliefs_c(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l,t_l-ss_time)))
                            iterations(ind,n_l,n_l2,P_l,a_l,t_l-ss_time)=iterations(ind,n_l,n_l2,P_l,a_l,t_l-ss_time)+sum(beliefs_c(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l,t_l-ss_time))
                            F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l,t_l-ss_time)=dble(sum(beliefs_c(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l,t_l-ss_time)))/dble(iterations(ind,n_l,n_l2,P_l,a_l,t_l-ss_time))*F(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l) +&
                                                        dble(iterations(ind,n_l,n_l2,P_l,a_l,t_l-ss_time)-sum(beliefs_c(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l,t_l-ss_time)))/dble(iterations(ind,n_l,n_l2,P_l,a_l,t_l-ss_time))*F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l,t_l-ss_time) 
                            if (minval(F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l,t_l-ss_time))<0)then 
                                print*,'error in beliefs'
                                print*,F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l,t_l-ss_time)
                                read*,pause_k
                            end if
                        end if
                end do;end do;end do;end do;end do
            end if
        end do
    end do
!av
        ! In case I don't have observations for a given state, I consider that the transition pr 
        ! is the same for all possible future states
        
        do t_l=1,transition_time;do P_l=1,P_max; do ind=1,2*P_l-1; do n_l=1,3; do n_l2=1,1; do a_l=1,types_a
            if (iterations(ind,n_l,n_l2,P_l,a_l,t_l)==0) then !iterations(1,2,1,5,4,:)
                F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l,t_l)=1.0d0/dble(2*P_l-1)
            end if
            if (isnan(F_new(ind,1,n_l,n_l2,P_l,a_l,t_l))) then
                print*,'error in generate beliefs'
            end if
            F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l,t_l)=F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l,t_l)/sum(F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l,t_l))
        end do;end do;end do;end do; end do;end do
        
        !Beliefs don't depend on n'
        F_new(:,:,:,2,:,:,:)=F_new(:,:,:,1,:,:,:)
        F_new(:,:,:,3,:,:,:)=F_new(:,:,:,1,:,:,:)
        F_new(:,:,1,3,:,:,:)=0.0d0
        
        n_tr=sum(total_N,2)/(count*dble(sims_tr))
        do P_l=1,2
            stats_out(P_l)=sum(av_output_per_acre(:,:,P_l))/(dble(sims_tr))
        end do
        
        av_drilling_rate=sum(av_drilling_rate_s,2)/(dble(sims_tr))
    
    
end subroutine
    
    
    
    
    
    
