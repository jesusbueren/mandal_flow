module data_moments
    use dimensions
    implicit none
    integer,parameter::mc_size=100
    double precision,dimension(types_a,2,mc_size)::pr_d_a_n_0
    double precision,dimension(types_a,mc_size)::pr_D_0
    double precision:: best_L
    character(len=3)::mc2s    
    integer::mc_num
end module
    
subroutine monte_carlo
    use nrtype; use dimensions; use simulation; use data_moments
    implicit none
    real(DP),dimension(par)::p_0
    double precision,dimension(types_a,2)::pr_d_a_n
    double precision,dimension(types_a)::pr_D_model
    double precision,dimension(9,plots_v(1),5,sims)::sim_data
    double precision::log_L,ftol,fret
    double precision,dimension(par+1,par)::p_g
    double precision,dimension(par+1)::y
    integer::iter,p_l,mc_l
    double precision,dimension(par,par)::xi
    interface 
        function log_likelihood_mc(p)
            use dimensions
            double precision,dimension(:),intent(in)::p
            double precision::log_likelihood_mc
        end function log_likelihood_mc
        SUBROUTINE amoeba(p,y,ftol,func,iter)
	        USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,iminloc,nrerror,swap
	        IMPLICIT NONE
	        INTEGER(I4B), INTENT(OUT) :: iter
	        REAL(DP), INTENT(IN) :: ftol
	        REAL(DP), DIMENSION(:), INTENT(INOUT) :: y
	        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: p
	        INTERFACE
		        FUNCTION func(x)
		        USE nrtype
		        IMPLICIT NONE
		        REAL(DP), DIMENSION(:), INTENT(IN) :: x
		        double precision :: func
		        END FUNCTION func
            END INTERFACE
        end subroutine
        SUBROUTINE powell(p,xi,ftol,iter,fret)
	        USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	        USE nr, ONLY : linmin
	        IMPLICIT NONE
	        REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
	        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: xi
	        INTEGER(I4B), INTENT(OUT) :: iter
	        REAL(DP), INTENT(IN) :: ftol
	        REAL(DP), INTENT(OUT) :: fret
        end subroutine
    end interface
    
    !Set a DGP
    !share water (CES)
    p_0(1)=0.5d0 
    !Productivity p
    p_0(2)=20.0d0 
    !Var taste shock
    p_0(3)=2.0d0 
    !flow p
    p_0(4)=0.4d0 
    !share non-zombie
    p_0(5)=0.6d0 
    
    open(unit=12, file=path_results//"\montecarlo\parameters_nod.txt",status='replace')
            write(12,'(<5>f10.5)'),p_0
        close(12)
    
    !Solve model for the DGP & generate simulated data
    montecarlo=1
    call compute_eq_village_1(p_0,pr_d_a_n,pr_D_model)
    montecarlo=0
    
    print*,'map moments D',pr_D_model
    print*,'map moments pr d',pr_d_a_n
    
    !Read simulated data
    OPEN(UNIT=12, FILE=path_results//"artificial_data_montecarlo.txt")
        read(12,*),sim_data
    close(12)

    !Create a sample with 2862 observations and compute the moments
    do mc_l=1,100
        call create_artificial_moments(sim_data,pr_D_0(:,mc_l),pr_d_a_n_0(:,:,mc_l))
    end do
    do mc_l=1,100
        mc_num=mc_l
        write (mc2s, "(I3.3)") mc_l     
    
        best_L=1000.0d0
        
        !Set a guess of parameters
        p_g(1,:)=(/0.48972d0,20.43374d0,2.08024d0,0.41467d0,0.54988d0/)

 
        do p_l=2,par+1
            p_g(p_l,:)=p_g(1,:)
            p_g(p_l,p_l-1)=p_g(1,p_l-1)*0.95d0
        end do
        
        !Change parameters to the (-Inf;Inf) real line
        do p_l=1,par+1
            p_g(p_l,1)=log(p_g(p_l,1)/(1.0d0-p_g(p_l,1)))
            p_g(p_l,2:3)=log(p_g(p_l,2:3))
            p_g(p_l,4)=log(p_g(p_l,4))
            p_g(p_l,5)=log(p_g(p_l,5)/(1.0d0-p_g(p_l,5)))
            y(p_l)=log_likelihood_mc(p_g(p_l,:))
            if (p_l<par+1)then
                xi(p_l,p_l)=1.0d0
            end if

        end do
        
        ftol=1.0d-2
        call amoeba(p_g,y,ftol,log_likelihood_mc,iter)
        
        p_g(1,1)=1.0d0/(1.0d0+exp(-p_g(1,1))) !p_g(1,:)
        p_g(1,2:3)=exp(p_g(1,2:3))
        p_g(1,4)=exp(p_g(1,4))
        p_g(1,5)=1.0d0/(1.0d0+exp(-p_g(1,5)))

    end do
    print*,'End of MonteCarlo'
    pause
    

end subroutine
    
subroutine create_artificial_moments(sim_data,pr_D_model,pr_d_a_n)
use nrtype; use dimensions;use cadastral_maps; use simulation
implicit none
double precision,dimension(9,plots_v(1),5,sims),intent(in)::sim_data
double precision,dimension(types_a),intent(out)::pr_D_model
double precision,dimension(types_a,2),intent(out)::pr_d_a_n
double precision,dimension(types_a)::counter_D
double precision,dimension(types_a,2)::counter_a_n
integer::s_l,s_ind,i_l,i,i_ind,A,n_l,t_l
double precision::u

pr_D_model=0.0d0
    pr_d_a_n=0
    counter_D=0
    counter_a_n=0
    do i=1,2862
        !sample a simulation
        call RANDOM_NUMBER(u)
        s_l=-9
        s_ind=1
        do while (s_l==-9)
            if (u<dble(s_ind)*1.0d0/dble(sims)) then
                s_l=s_ind
            else
                s_ind=s_ind+1
            end if
        end do
                  
        !sample a plot
        call RANDOM_NUMBER(u)
        i_l=-9
        i_ind=1
        do while (i_l==-9)
            if (u<dble(i_ind)*1.0d0/dble(plots_v(1))) then
                i_l=i_ind
            else
                i_ind=i_ind+1
            end if
        end do
        
        !some drilling in 5 periods
        A=PA_type(i_l,2,1,s_l)
        if (sum(sim_data(9,i_l,:,s_l))>0 .or. sum(sim_data(7,i_l,:,s_l))>0) then
            pr_D_model(A)=pr_D_model(A)+1
        end if
        counter_D(A)=counter_D(A)+1
        
        !drilling by area and number of wells
        do t_l=1,5
            n_l=sim_data(7,i_l,t_l,s_l)+1
            if (n_l<3) then
                if (sim_data(9,i_l,t_l,s_l)==1) then
                    pr_d_a_n(A,n_l)=pr_d_a_n(A,n_l)+1
                end if
                counter_a_n(A,n_l)=counter_a_n(A,n_l)+1
            end if                
        end do
        
    end do
    
    pr_D_model=pr_D_model/counter_D
    pr_d_a_n=pr_d_a_n/counter_a_n
    
end subroutine
    
subroutine compute_eq_village_1(p,pr_d_a_n,pr_D_model)
use dimensions; use cadastral_maps
implicit none
double precision,dimension(par),intent(in)::p
double precision,dimension(types_a),intent(out)::pr_D_model
double precision,dimension(types_a,2),intent(out)::pr_d_a_n
integer::v_l    
double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,types_a)::F
double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types)::CCP
double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types)::V_fct,V_social
integer,dimension(plots_in_map)::n_dist
double precision::mean_N,social_output,private_output,d_rate,delta_PV_model
double precision,dimension(2*P_max-1,3)::pr_N_n_v
double precision,dimension(3,types_a,villages)::pr_na_v

    CCP=0.07d0
    do v_l=1,1 !only one village
        print*,'village',v_l
        call compute_eq_f_ccp(p,F,CCP,v_fct,v_social,&
            n_dist,v_l,mean_n,social_output,private_output,&
            pr_d_a_n,pr_N_n_v,pr_na_v,d_rate,pr_D_model,delta_PV_model)
    end do

    end subroutine
    
    
function log_likelihood_mc(p)
    use dimensions; use data_moments; use simulation
    implicit none
    double precision,dimension(:),intent(in)::p
    double precision,dimension(par)::params
    double precision::log_likelihood_mc
    double precision,dimension(types_a)::pr_D_model
    double precision,dimension(types_a,2)::pr_d_a_n
    
    params(1)=1.0d0/(1.0d0 + exp(-p(1))) 
    params(2:3)=exp(p(2:3))
    params(4)=exp(p(4))
    params(5)=1.0d0/(1.0d0 + exp(-p(5))) 

    print*,'parameters',params
    call compute_eq_village_1(params,pr_d_a_n,pr_D_model)
    
    
    log_likelihood_mc=sum((pr_d_a_n-pr_d_a_n_0(:,:,mc_num))**2.0d0)+sum((pr_D_model-pr_D_0(:,mc_num))**2.0d0)
    print*,'map moments D',pr_D_model
    print*,'map moments pr d',pr_d_a_n
    
    print*,'log_likelihood_mc',log_likelihood_mc
    
    if (log_likelihood_mc<best_L) then
        best_L=log_likelihood_mc
        open(unit=12, file=path_results//"\montecarlo\parameters_far"//trim(mc2s)//".txt",status='replace')
            write(12,'(<6>f10.5)'),params,log_likelihood_mc
        close(12)
        open(unit=12, file=path_results//"\montecarlo\moment_D_far"//trim(mc2s)//".txt",status='replace')
            write(12,'(<5>f10.5)'),pr_D_model
        close(12)
        open(unit=12, file=path_results//"\montecarlo\moment_pr_d_far"//trim(mc2s)//".txt",status='replace')
            write(12,'(<10>f10.5)'),pr_d_a_n
        close(12)
    end if
end function log_likelihood_mc





    
  