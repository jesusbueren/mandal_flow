subroutine generate_beliefs(CCP,V_fct,V_social,Ef_v,n_initial,F_new,v_l,iterations,mean_N,social_output,private_output,joint_pr, &
        pr_d_a_n,pr_N_n,pr_na,d_rate,pr_D_model,delta_PV_model)
    use cadastral_maps; use primitives; use simulation
    implicit none
    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types),intent(in)::CCP
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types),intent(in)::V_fct,V_social
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types),intent(in)::Ef_v 
    integer,dimension(plots_in_map,1),intent(inout)::n_initial
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,types_a),intent(out)::F_new
    integer,intent(in)::v_l
    integer(8),dimension(2*P_max-1,3,3,P_max,types_a),intent(out)::iterations
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types),intent(out)::joint_pr
    double precision,dimension(types_a,2),intent(out)::pr_d_a_n
    integer(8),dimension(2*P_max-1,3,P_max,types_a,unobs_types)::counter_u
    integer(8)::T,d_it
    integer(8),dimension(plots_in_map,3)::state,state_old
    integer(8)::i_l,j_l,t_l,ind,N_all,n_l,P,A,P_l,n_l2,it,m_l,it_min,a_l,u_l,ind2,N_all2,i,r_l,it2,t_index,ind_t5
    integer::s_l
    double precision::u_d,u_s,u_f,u_m,it3,y,dc,drilling
    integer(8):: its
    double precision,allocatable,dimension(:,:)::NPV,total_N,NPV_PV
    integer,allocatable,dimension(:,:,:)::pr_N_n_it_c
    double precision,intent(out)::mean_N,social_output,private_output,d_rate,delta_PV_model
    double precision,dimension(types_a),intent(out)::pr_D_model
    integer(8),dimension(2*P_max-1,2*P_max-1,3,3,P_max,types_a)::beliefs_c
    integer(8),parameter::burn_t=99
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,types_a)::F
    double precision,dimension(P_max)::dist
    character::continue_k
    integer,allocatable:: seed_new(:)
    character(LEN=1)::s_c1,v_s1
    character(LEN=2)::s_c2,v_s2
    double precision,dimension(8)::draw
    integer(8),dimension(3,types_a)::count_a_n
    integer(8),dimension(types_a)::count_a
    integer(8),dimension(types_a)::counter_nd
    double precision,dimension(2*P_max-1,3),intent(out)::pr_N_n
    double precision,dimension(3,types_a),intent(out)::pr_na
    integer(8),dimension(plots_in_map,types_a)::aux_i
    double precision,dimension(plots_in_map,types_a)::PV_i
    double precision,dimension(sims)::total_a
    interface
    subroutine random_value ( seed, r )
    implicit none
    double precision,intent(out):: r
    integer,dimension(:),intent(inout):: seed
    end
    end interface

    !active_plots(415,1,1)
    
    !Set number of iterations
    T=10000000/plots_v(v_l)  !
    its=T-burn_t-1
    allocate ( NPV(its,sims))
    allocate ( total_N(its,sims))
    allocate ( NPV_PV(its,sims))
    allocate ( pr_N_n_it_c(2*P_max-1,3,sims))

    
    !prior of beliefs
    F_new=-9.0
    iterations=0
    do P_l=1,P_max;do ind=1,2*P_l-1; do n_l=1,3; do n_l2=1,1; do a_l=1,types_a
        F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)=1.0d0/dble(2*P_l-1)
        iterations(ind,n_l,n_l2,P_l,a_l)=1
    end do;end do;end do;end do;end do

    
    !Prior initial conditions
    counter_u=0
    do P_l=1,P_max;do n_l=1,3; do u_l=1,unobs_types;do a_l=1,types_a
        counter_u(1:2*P_l-1,n_l,P_l,a_l,u_l)=1
    end do;end do;end do;end do
    
    
    !print*,'smthg'
    !Call seed number
    call random_seed(size = s_l)
    allocate(seed_new(s_l))
    seed_new=4321
    call random_seed(PUT=seed_new)
    
    
    NPV=0.0d0
    NPV_PV=0.0d0
    pr_N_n_it_c=0
    pr_d_a_n=0.0d0
    total_N=0.0d0

    d_rate=0.0d0
    pr_D_model=0.0d0
    beliefs_c=0
    ind_t5=0
    it=0

    count_a_n=0
    counter_nd=0
    total_a=0.0d0
    
    !count number of plots with a given area
    count_a=0
    do i_l=1,plots_v(v_l)
        count_a(PA_type(i_l,2,v_l,1))=count_a(PA_type(i_l,2,v_l,1))+1
    end do

    
    if (generate_data==1)then
        if (v_l>9) then
            write (v_s2, "(I2)") v_l
            OPEN(UNIT=12, FILE=path_results//"artificial_data_v"//v_s2//".txt")
        else
            write (v_s1, "(I1)") v_l
            OPEN(UNIT=12, FILE=path_results//"artificial_data_v"//v_s1//".txt")
        end if
    end if
    
    if (montecarlo==1)then
        OPEN(UNIT=12, FILE=path_results//"artificial_data_montecarlo.txt")
    end if

    
    do s_l=1,sims;
        t_index=0
        state_old=1
        do t_l=1,T-1;
        t_index=t_index+1
        !all plots are started indexed as undeveloped as they have wells or try to drill they will be marked as developed
        if (t_index==1) then
            aux_i=0
            PV_i=0.0d0
        end if
        !print*,t_l
        !simulate monsoon next period
        call random_value( seed_new, u_m )
        if (u_m<PI_m(1,v_l))then
            m_l=1
        else
            m_l=2
        end if

        if (t_l>burn_t) then
            it2=0
            it3=0.0d0            
        end if
        if (t_l>1) then
            state(:,1)=n_initial(:,1)
        else
            n_initial(:,1)=1
            state(:,1)=1
        end if
        !print*,'t_l',t_l,'av number of wells per plot',real(sum(n_initial(1:plots_v(v_l),1))-plots_v(v_l))/real(plots_v(v_l))
        do i_l=1,plots_v(v_l)
            do r_l=1,8
                call random_value(seed_new,draw(r_l))
            end do  
            d_it=0
            if (active_plots(i_l,v_l,s_l)==1) then  
                if (t_l==1) then
                    total_a(s_l)=total_a(s_l)+area(PA_type(i_l,2,v_l,s_l))
                end if
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
                
                !Commment/uncomment island model

                !print*,'Island model!'
                !N_all=n_l
                !n_l=state(i_l,1)
                !P=1
            
                !Locate position in the state space wrt to the CCP, PI_s_v and ,PI_f_v
                if (n_l==1) then
                    ind=N_all 
                elseif (n_l==2) then
                    aux_i(i_l,A)=1
                    ind=N_all-1
                elseif (n_l==3) then
                    aux_i(i_l,A)=1
                    ind=N_all-2
                else
                    print*,'error generating beliefs'
                end if 
                    
                if (ind==0) then
                    print*,'paused'
                    !print'(<P_max+2>I6)',i_l,PA_type(i_l,1,v_l,s_l),neighbors(i_l,:,v_l,s_l)
                    read*,continue_k
                end if
                state(i_l,3)=ind
                !Count transitions (in the first iteration state_old is undefined: no problem)
                if (t_l>=burn_t) then
                    beliefs_c(state_old(i_l,3),state(i_l,3),state_old(i_l,1),1,P,A)=&
                    beliefs_c(state_old(i_l,3),state(i_l,3),state_old(i_l,1),1,P,A)+1
                    !Compute joint distribution state variables and unobserved heterogeneity type
                    
                    counter_u(state(i_l,3),state(i_l,1),P,A,unobs_types_i(i_l,v_l,s_l))=counter_u(state(i_l,3),state(i_l,1),P,A,&
                        unobs_types_i(i_l,v_l,s_l))+1
                end if
                !Compute NPV
                if (t_l>burn_t) then  
                    it2=it2+1
                    pr_N_n_it_c(ind,n_l,s_l)=pr_N_n_it_c(ind,n_l,s_l)+1.0d0                    
                    total_N(t_l-burn_t,s_l)=total_N(t_l-burn_t,s_l)+dble(n_l-1)
                    if (n_l>1) then
                        it3=it3+1
                    end if
                    count_a_n(n_l,A)=count_a_n(n_l,A)+1
                    if (n_l<3) then
                        pr_d_a_n(A,n_l)=dble(count_a_n(n_l,A)-1)/dble(count_a_n(n_l,A))*pr_d_a_n(A,n_l)+&
                            1.0d0/dble(count_a_n(n_l,A))*CCP(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l))
                    end if
                end if
                
                !Well drilling decision and failures/successes
                y=Ef_v(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l))
                dc=0.0d0
                drilling=0.0d0
                if (n_l==1) then !no well
                    u_d=draw(1) 
                    if (u_d<CCP(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l))) then !decides to drill
                        aux_i(i_l,A)=1
                        drilling=1.0d0
                        u_s=draw(2) 
                        d_it=1
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
                        d_it=1
                        u_s=draw(4) 
                        if (u_s<PI_s_v(ind,n_l,P,v_l)) then !successful attempt
                            dc=c_d
                            u_f=draw(5) 
                            if (u_f<PI_fm(N_all-1,m_l,v_l,unobs_types_i(i_l,v_l,s_l))) then !failure of the previous well
                                n_initial(i_l,1)=n_l
                            else
                                n_initial(i_l,1)=n_l+1
                            end if
                        else !unsuccessful attempt
                            dc=c_s
                            u_f=draw(6) 
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
                    if (u_f<PI_fm(N_all-1,m_l,v_l,unobs_types_i(i_l,v_l,s_l))**2) then !failure of the two wells 
                        n_initial(i_l,1)=n_l-2
                    elseif (u_f<(1.0d0-PI_fm(N_all-1,m_l,v_l,unobs_types_i(i_l,v_l,s_l)))**2+&
                    PI_fm(N_all-1,m_l,v_l,unobs_types_i(i_l,v_l,s_l))**2) then !failure of none (this is correct)
                        n_initial(i_l,1)=n_l
                    else !failure of one
                        n_initial(i_l,1)=n_l-1
                    end if 
                else
                    print*,'error in gen beliefs 2'
                end if 
                PV_i(i_l,A)=V_fct(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l))/area(A)
                if (t_l>burn_t) then
                    if (isnan((dble(n_l)-1.0d0)*((1.0d0-PI_s_v(ind,1,P,v_l))/PI_s_v(ind,1,P,v_l)*c_d+c_s))) then
                        print*,'error in computing expected drilling cost'
                    end if
                    NPV(t_l-burn_t,s_l)=NPV(t_l-burn_t,s_l)+V_social(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l))
                    NPV_PV(t_l-burn_t,s_l)=NPV_PV(t_l-burn_t,s_l)+V_fct(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l))
                    d_rate=dble(sum(count_a_n)-1)/dble(sum(count_a_n))*d_rate+1.0d0/dble(sum(count_a_n))*drilling
                end if
            else
                counter_nd(PA_type(i_l,2,v_l,s_l))=counter_nd(PA_type(i_l,2,v_l,s_l))+1
            end if   
            
            
            if ((generate_data==1 .or. montecarlo==1) .and. t_l>T-6) then
                if (active_plots(i_l,v_l,s_l)==0) then  
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
                    n_l=state(i_l,1)
                end if
                write(12,'(I10,I10,I10,I10,F10.2,I10,I10,I10,I10)'),v_l,s_l,i_l,active_plots(i_l,v_l,s_l), &
                                                areas(v_l,i_l),t_l-(T-6),n_l-1,N_all-(n_l-1)-1,d_it
            end if
        end do
        if (t_index==5) then
            ind_t5=ind_t5+1
            pr_D_model=dble(ind_t5-1)/dble(ind_t5)*pr_D_model+1.0d0/dble(ind_t5)*dble(sum(aux_i(1:plots_v(v_l),:),1))/dble(count_a)
            delta_PV_model=dble(ind_t5-1)/dble(ind_t5)*delta_PV_model+1.0d0/dble(ind_t5)*(sum(aux_i*PV_i)/dble(sum(aux_i))-sum((1-aux_i)*PV_i)/dble(plots_v(v_l)-sum(aux_i)))
            !print '(<types_a> F6.3)',pr_D_model
            t_index=0
            aux_i=0
        end if
        if (t_l>burn_t) then
            total_N(t_l-burn_t,s_l)=total_N(t_l-burn_t,s_l)/total_a(s_l)
            NPV(t_l-burn_t,s_l)=NPV(t_l-burn_t,s_l)/total_a(s_l)
            NPV_PV(t_l-burn_t,s_l)=NPV_PV(t_l-burn_t,s_l)/total_a(s_l)
        end if
            
        !pr_N_n_it_c
        !Store current state
        state_old=state
        !Compute beliefs
        if (t_l>burn_t) then
            F=0.0d0
            do P_l=2,P_max; do ind=1,2*P_l-1; do n_l=1,3; do n_l2=1,1;do a_l=1,types_a
                    if (n_l==1 .and. n_l2==3) then
                        F(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)=-9.0d0
                    elseif(sum(beliefs_c(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l))>0) then
                        F(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)=&
                            dble(beliefs_c(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l))/dble(sum(beliefs_c(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)))
                        iterations(ind,n_l,n_l2,P_l,a_l)=&
                            iterations(ind,n_l,n_l2,P_l,a_l)+sum(beliefs_c(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l))
                        F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)=&
                            dble(sum(beliefs_c(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)))/dble(iterations(ind,n_l,n_l2,P_l,a_l))&
                            *F(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l) +&
                            dble(iterations(ind,n_l,n_l2,P_l,a_l)-sum(beliefs_c(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)))&
                            /dble(iterations(ind,n_l,n_l2,P_l,a_l))*F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l) 
                        if (minval(F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l))<0)then 
                            print*,'error in beliefs'
                        end if
                    end if
            end do;end do;end do;end do;end do
        end if
        it=it+1
    end do !total_N(1:200,s_l)
    end do 
    
    if (generate_data==1 .or. montecarlo==1)then
        close(12)
    end if
    
    ! In case I don't have observations for a given state, I consider that the transition pr 
    ! is the same for all possible future states
    it=0
    it2=0
    do P_l=1,P_max; do ind=1,2*P_l-1; do n_l=1,3; do n_l2=1,1; do a_l=1,types_a
        it2=it2+1
        if (iterations(ind,n_l,n_l2,P_l,a_l)<20) then
            it=it+1
            F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)=0.0d0!1.0d0/dble(2*P_l-1)
            F_new(ind,ind,n_l,n_l2,P_l,a_l)=1.0d0
        end if
        if (isnan(F_new(ind,1,n_l,n_l2,P_l,a_l))) then
            print*,'error in generate beliefs'
        end if
        F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)=F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)/sum(F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)) 
    end do;end do;end do;end do; end do
    !Beliefs don't depend on n'
    F_new(:,:,:,2,:,:)=F_new(:,:,:,1,:,:)
    F_new(:,:,:,3,:,:)=F_new(:,:,:,1,:,:)
    F_new(:,:,1,3,:,:)=0.0d0
    
    !Beliefs for plots with no neighbors are degenerate: they know what the future will look like cond on their own outcomes
    P_l=1
    ind=1
    do n_l=1,3; do n_l2=1,min(n_l+1,3);do a_l=1,types_a
        F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)=1.0d0
    end do;end do;end do
    
    social_output=sum(NPV)/dble(sims*its)
    private_output=sum(NPV_PV)/dble(sims*its)
    mean_N=sum(total_N)/dble(its*sims)    
    
    do n_l=1,3;do ind=1,2*P_max-1
        pr_N_n(ind,n_l)=dble(sum(pr_N_n_it_c(ind,n_l,:)))/dble(sum(pr_N_n_it_c(:,n_l,:))) 
    end do;end do
    
    do a_l=1,types_a
        pr_na(:,a_l)=dble(count_a_n(:,a_l))/dble(sum(count_a_n(:,a_l)))
    end do
    
    joint_pr=0.0d0
     do P_l=1,P_max; do u_l=1,unobs_types;do a_l=1,types_a
         if (sum(counter_u(:,:,P_l,a_l,u_l))==0.0d0)then
             !print*,'got here'
             joint_pr(:,:,P_l,a_l,u_l)=0.0d0
         else
            joint_pr(1:2*P_l-1,:,P_l,a_l,u_l)=dble(counter_u(1:2*P_l-1,:,P_l,a_l,u_l))&
                /dble(sum(counter_u(1:2*P_l-1,:,P_l,a_l,u_l)))
        end if    
     end do;end do;end do
     
     !Compute the unconditional probability of drilling (taking into account undeveloped plots)
     do a_l=1,types_a
        pr_d_a_n(a_l,1)=pr_d_a_n(a_l,1)*dble(count_a_n(1,a_l))/dble(count_a_n(1,a_l)+counter_nd(a_l))
     end do
     
          
 
end subroutine

    
    subroutine random_value ( seed, r )

!*****************************************************************************80
!
!! RANDOM_VALUE generates a random value R.
!
!  Discussion:
!
!    This is not a good random number generator.  It is a SIMPLE one.
!    It illustrates a model which works by accepting an integer seed value
!    as input, performing some simple operation on the seed, and then
!    producing a "random" real value using some simple transformation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) R, the random value.
!
  implicit none

  double precision,intent(out):: r
  integer,dimension(:),intent(inout):: seed

  seed = mod ( seed, 65536 )
  seed = mod ( ( 3125 * seed ), 65536 )
  r = dble ( seed(1) ) / 65536.0D+00

  return
end