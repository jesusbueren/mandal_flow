subroutine generate_beliefs(CCP,V_fct,V_social,Ef_v,n_initial,F_new,v_l,iterations,mean_N,social_output,private_output,joint_pr,pr_d_a_n,pr_N_n,pr_na)
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
    integer(8)::i_l,j_l,t_l,ind,N_all,n_l,P,A,P_l,n_l2,it,m_l,it_min,a_l,u_l,ind2,N_all2,i,r_l,s_l
    double precision::u_d,u_s,u_f,u_m,it2,it3,y,dc
    integer(8):: its
    double precision,allocatable,dimension(:,:)::NPV,total_N,NPV_PV,total_f,cultivated_area,high_type_fraction
    integer,allocatable,dimension(:,:,:)::pr_N_n_it_c
    double precision,allocatable,dimension(:,:,:,:)::NPV_by_a_n
    double precision,intent(out)::mean_N,social_output,private_output
    integer(8),dimension(2*P_max-1,2*P_max-1,3,3,P_max,types_a)::beliefs_c
    integer(8),dimension(1)::seed=321,seed2
    integer(8),parameter::burn_t=99
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,types_a)::F
    double precision,dimension(P_max)::dist
    character::continue_k
    integer(8),dimension(1):: seed_new
    character(LEN=1)::s_c1
    character(LEN=2)::s_c2,v_s
    double precision,dimension(8)::draw
    integer,dimension(3,types_a)::count_a_n
    double precision,dimension(2*P_max-1,3),intent(out)::pr_N_n
    double precision,dimension(3,types_a),intent(out)::pr_na
    double precision,dimension(plots_in_map)::aux_i
    
    !active_plots(415,1,1)
    
    
    !Set number of iterations
    T=8000000/plots_v(v_l) 
    its=T-burn_t-1
    allocate ( NPV(its,sims))
    allocate ( cultivated_area(its,sims))
    allocate ( high_type_fraction(its,sims))
    allocate ( total_N(its,sims))
    allocate ( total_f(its,sims))
    allocate ( NPV_PV(its,sims))
    allocate ( pr_N_n_it_c(2*P_max-1,3,sims))
    allocate ( NPV_by_a_n(its,sims,types_a,3))

    
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
    
    seed_new=4321
    NPV=0.0d0
    NPV_by_a_n=0.0d0
    pr_N_n_it_c=0
    pr_d_a_n=0.0d0

    total_f=0
    beliefs_c=0
    it=0

    count_a_n=0
    
    write (v_s, "(I2)") v_l
    if (generate_data==1)then
        OPEN(UNIT=12, FILE=path_results//"artificial_data_v"//v_s//".txt")
    end if

    
    do s_l=1,sims;
        state_old=1
        do t_l=1,T-1;
        !print*,t_l
        !simulate monsoon next period
        call random_value( seed_new, u_m )
        if (u_m<PI_m(1,v_l))then
            m_l=1
        else
            m_l=2
        end if

        if (t_l>burn_t) then
            NPV(t_l-burn_t,s_l)=0.0d0
            cultivated_area(t_l-burn_t,s_l)=0.0d0
            high_type_fraction(t_l-burn_t,s_l)=0.0d0
            NPV_PV(t_l-burn_t,s_l)=0.0d0
            it2=0.0d0
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
                    ind=N_all-1
                elseif (n_l==3) then
                    ind=N_all-2
                else
                    print*,'error generating beliefs'
                end if 
                    
                if (ind==0) then
                    print*,'paused'
                    print'(<P_max+2>I6)',i_l,PA_type(i_l,1,v_l,s_l),neighbors(i_l,:,v_l,s_l)
                    read*,continue_k
                end if
                state(i_l,3)=ind
                !Count transitions (in the first iteration state_old is undefined: no problem)
                if (t_l>=burn_t) then
                    beliefs_c(state_old(i_l,3),state(i_l,3),state_old(i_l,1),1,P,A)=&
                    beliefs_c(state_old(i_l,3),state(i_l,3),state_old(i_l,1),1,P,A)+1
                    !Compute joint distribution state variables and unobserved heterogeneity type
                    counter_u(state(i_l,3),state(i_l,1),P,A,unobs_types_i(i_l,v_l,s_l))=counter_u(state(i_l,3),state(i_l,1),P,A,unobs_types_i(i_l,v_l,s_l))+1
                end if
                !Compute NPV
                if (t_l>burn_t) then  
                    it2=it2+1.0d0
                    pr_N_n_it_c(ind,n_l,s_l)=pr_N_n_it_c(ind,n_l,s_l)+1.0d0                    
                    total_N(t_l-burn_t,s_l)=dble(it2-1)/dble(it2)*total_N(t_l-burn_t,s_l)+1.0d0/dble(it2)*(n_l-1)
                    if (n_l>1) then
                        it3=it3+1
                        cultivated_area(t_l-burn_t,s_l)=cultivated_area(t_l-burn_t,s_l)+area(A)
                        if (unobs_types_i(i_l,v_l,s_l)==unobs_types) then
                            high_type_fraction(t_l-burn_t,s_l)=dble(it3-1)/dble(it3)*high_type_fraction(t_l-burn_t,s_l)+1.0d0/dble(it3)*1.0d0
                        else
                            high_type_fraction(t_l-burn_t,s_l)=dble(it3-1)/dble(it3)*high_type_fraction(t_l-burn_t,s_l)+1.0d0/dble(it3)*0.0d0
                        end if
                    end if
                    count_a_n(n_l,A)=count_a_n(n_l,A)+1
                    NPV_by_a_n(t_l-burn_t,s_l,A,n_l)=dble(count_a_n(n_l,A)-1)/dble(count_a_n(n_l,A))*NPV_by_a_n(t_l-burn_t,s_l,A,n_l)+1.0d0/dble(count_a_n(n_l,A))*(V_fct(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l)))
                    if (n_l<3) then
                        pr_d_a_n(A,n_l)=dble(count_a_n(n_l,A)-1)/dble(count_a_n(n_l,A))*pr_d_a_n(A,n_l)+1.0d0/dble(count_a_n(n_l,A))*CCP(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l))
                    end if
                end if
                
                !Well drilling decision and failures/successes
                y=Ef_v(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l))
                dc=0.0d0
                if (n_l==1) then !no well
                    u_d=draw(1) 
                    if (u_d<CCP(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l))) then !decides to drill
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
                if (t_l>burn_t) then
                    aux_i(it2)=V_social(ind,1,P,A,unobs_types_i(i_l,v_l,s_l))/area(A) !aux_i(1:648)
                    if (isnan((dble(n_l)-1.0d0)*((1.0d0-PI_s_v(ind,1,P,v_l))/PI_s_v(ind,1,P,v_l)*c_d+c_s))) then
                    print*,'error in computing expected drilling cost'
                    end if
                    NPV(t_l-burn_t,s_l)=dble(it2-1)/dble(it2)*NPV(t_l-burn_t,s_l)+1.0d0/dble(it2)*V_social(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l))/area(A) !(V_social(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l))-(dble(n_l)-1.0d0)*((1.0d0-PI_s_v(ind,1,P,v_l))/PI_s_v(ind,1,P,v_l)*c_d+c_s))/area(A) 
                    NPV_PV(t_l-burn_t,s_l)=dble(it2-1)/dble(it2)*NPV_PV(t_l-burn_t,s_l)+1.0d0/dble(it2)*V_fct(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l))/area(A) !(V_fct(ind,n_l,P,A,unobs_types_i(i_l,v_l,s_l))-(dble(n_l)-1.0d0)*((1.0d0-PI_s_v(ind,1,P,v_l))/PI_s_v(ind,1,P,v_l)*c_d+c_s))/area(A) 
                end if
            end if   
            
            
            if (generate_data==1 .and. t_l>T-6) then
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
                write(12,'(<8>I10)'),v_l,s_l,i_l,active_plots(i_l,v_l,s_l),t_l-(T-6),n_l-1,N_all-(n_l-1)-1,d_it
            end if
        end do
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
                        F(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)=dble(beliefs_c(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l))/dble(sum(beliefs_c(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)))
                        iterations(ind,n_l,n_l2,P_l,a_l)=iterations(ind,n_l,n_l2,P_l,a_l)+sum(beliefs_c(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l))
                        F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)=dble(sum(beliefs_c(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)))/dble(iterations(ind,n_l,n_l2,P_l,a_l))*F(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l) +&
                                                    dble(iterations(ind,n_l,n_l2,P_l,a_l)-sum(beliefs_c(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)))/dble(iterations(ind,n_l,n_l2,P_l,a_l))*F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l) 
                        if (minval(F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l))<0)then 
                            print*,'error in beliefs'
                        end if
                    end if
            end do;end do;end do;end do;end do
        end if
        it=it+1
    end do 
    end do 
    
    if (generate_data==1)then
        close(12)
    end if
    
    ! In case I don't have observations for a given state, I consider that the transition pr 
    ! is the same for all possible future states
    it=0
    it2=0
    do P_l=1,P_max; do ind=1,2*P_l-1; do n_l=1,3; do n_l2=1,1; do a_l=1,types_a
        it2=it2+1
        if (iterations(ind,n_l,n_l2,P_l,a_l)==0) then
            it=it+1
            F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)=1.0d0/dble(2*P_l-1)
        end if
        if (isnan(F_new(ind,1,n_l,n_l2,P_l,a_l))) then
            print*,'error in generate beliefs'
        end if
        F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)=F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)/sum(F_new(ind,1:2*P_l-1,n_l,n_l2,P_l,a_l)) !F_new(2,2,1,1,4,:)
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

    !print*,'cultivated area',sum(cultivated_area)/dble(its*sims)
    !print*,'Fraction of active plots from a high type ',sum(high_type_fraction)/dble(its*sims)
    
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
            joint_pr(1:2*P_l-1,:,P_l,a_l,u_l)=dble(counter_u(1:2*P_l-1,:,P_l,a_l,u_l))/dble(sum(counter_u(1:2*P_l-1,:,P_l,a_l,u_l)))
        end if    
     end do;end do;end do

     
 
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
  integer(8),dimension(1),intent(inout):: seed

  seed = mod ( seed, 65536 )
  seed = mod ( ( 3125 * seed ), 65536 )
  r = dble ( seed(1) ) / 65536.0D+00

  return
end