!Estimate parameters
subroutine estimation2(params_MLE,log_likeli)
    use dimensions; use primitives; use simulation
    implicit none
    double precision,dimension(par),intent(inout)::params_MLE
    double precision,intent(out)::log_likeli
    double precision::log_L,ftol,fret
    double precision,dimension(par+1,par)::p_g
    double precision,dimension(par+1)::y
    integer::iter,p_l,a_l,P_l2,v_l,ind,n_l,u_l
    interface 
        function log_likelihood2(p_MLE)
            use dimensions
            double precision,dimension(:),intent(in)::p_MLE
            double precision::log_likelihood2
        end function log_likelihood2
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
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP_old,CCP_mid
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::V_fct,V_social
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::Ef_v 
    double precision::dist
    integer::it
    integer(8),dimension(2*P_max-1,3,3,P_max,villages,unobs_types)::iterations_all
    double precision,dimension(par,par)::xi
    double precision, dimension(villages)::mean_N,mean_NPV,mean_budget
    double precision,dimension(villages)::village_fe
    character::pause_k
    
    
    !Fixing beliefs, estimate parameter
    !print*,'Initial Conditions'
    
    max_mle=99999999.0d0
    
    p_g(1,:)=params_MLE
    !!share water (CES)
    !p_g(1,1)=0.795d0 !0.5d0 
    !!Productivity p
    !p_g(1,2)=15.153d0 !20.0d0
    !!Var taste shock
    !p_g(1,3)=1.144d0 !2.1d0 !2.0d0
    !!flow p
    !p_g(1,4)=0.402d0 !0.5d0
    !!share non-zombie
    !p_g(1,5)=0.719d0 !0.7d0



    do p_l=2,par+1
        p_g(p_l,:)=p_g(1,:)
        p_g(p_l,p_l-1)=p_g(1,p_l-1)*0.9d0
    end do

    xi=0.0d0    
    !Change parameters to the (-Inf;Inf) real line
    do p_l=1,par+1
        p_g(p_l,1)=log(p_g(p_l,1)/(1.0d0-p_g(p_l,1)))
        p_g(p_l,2:3)=log(p_g(p_l,2:3))
        p_g(p_l,4)=log(p_g(p_l,4))
        p_g(p_l,5)=log(p_g(p_l,5)/(1.0d0-p_g(p_l,5)))
        y(p_l)=log_likelihood2(p_g(p_l,:))  

        if (p_l<par+1)then
            xi(p_l,p_l)=1.0d0
        end if
        !read*,pause_k
    end do 

    !print*,'likelihood_ini',y(1)
        
    ftol=1.0d-3
    call amoeba(p_g,y,ftol,log_likelihood2,iter)
     print*,'likelihood amoeba',y(1)
    !call powell(p_g(1,:),xi,ftol,iter,fret)
    
    log_likeli=y(1)
    p_g(1,1)=1.0d0/(1.0d0+exp(-p_g(1,1)))
    p_g(1,2:3)=exp(p_g(1,2:3))
    p_g(1,4)=exp(p_g(1,4))
    p_g(1,5)=1.0d0/(1.0d0+exp(-p_g(1,5)))


    
    params_MLE=p_g(1,:)
    print*,'likelihood amoeba',y(1)
    
    
    if (dist<0.0) then
        print*,'error in precision of iterations'   
    end if

   
end subroutine
    
function log_likelihood2(params_MLE)
    use dimensions; use simulation; use cadastral_maps; use primitives
    implicit none
    double precision,dimension(:),intent(in)::params_MLE
    double precision,dimension(par)::params
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP,CCP_opt   !CCP(1,2,1,1,1,:)
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::V_fct,V_social
    integer::i_l,t_l,type_l,a_l,p_l,v_l,ind,u_l,j_l,s_l,t,missing_x1,j_l2,ind2,n_l,r_l,q_l
    double precision::log_likelihood2,pr_non_zombie_II
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::Ef_v 
    double precision,dimension(unobs_types)::likelihood_i,likelihood_it,P_N2_N1,P_BigN2_BigN1
    character::end_k
    double precision,dimension(T_sim,plots_i,unobs_types)::av_CCP_uhe
    double precision,dimension(T_sim,plots_i)::av_CCP_it
    double precision,dimension(plots_i,unobs_types)::posterior_type
    character::pause_k
    double precision,dimension(types_a,2)::moment_own_nxa_model
    double precision,dimension(COV,1)::X
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,villages,types_a)::F
    integer,dimension(plots_in_map,villages)::n_dist,n_dist_opt
    double precision,dimension(villages)::mean_N,social_output,private_output,delta_PV_model

    double precision,dimension(unobs_types)::expost_types
    double precision,dimension(types_a,2,villages)::pr_d_a_n
    double precision,dimension(types_a,2)::pr_d_a_n_av
    double precision,dimension(types_a,2)::pr_d_a_av
    double precision,dimension(villages)::pr_drill_v
    double precision::fraction_constrained=0.7d0
    double precision,dimension(2*P_max-1,3,villages)::pr_N_n_v
    double precision,dimension(3,villages)::pr_n_v
    double precision,dimension(3,types_a,villages)::pr_na_v,pr_non_zombie_anv
    double precision,dimension(3,types_a,villages,2)::pr_na_vD
    double precision,dimension(2*P_max-1,3)::pr_N_n_av
    double precision,dimension(3)::pr_n_av
    double precision::u
    double precision,dimension(plots_i):: pr_non_zombie_i
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::empty=-9.0d0
    double precision,dimension(types_a)::pr_non_zombie_a_s
    double precision,dimension(wealth_quantiles,types_a)::pr_non_zombie_aw_s
    integer,dimension(types_a)::counter_a
    integer,dimension(wealth_quantiles,types_a)::counter_aw
    double precision,dimension(wealth_quantiles,types_a)::pr_d_wa
    double precision,dimension(wealth_quantiles,types_a,villages)::pr_non_zombie_wav
    double precision,dimension(villages)::d_rate
    double precision,dimension(types_a,villages)::pr_D_model
    double precision,dimension(types_a)::Pr_D_a_model
    
    
    params(1)=1.0d0/(1.0d0 + exp(-params_MLE(1))) 
    params(2:3)=exp(params_MLE(2:3))
    params(4)=exp(params_MLE(4))
    params(5)=1.0d0/(1.0d0 + exp(-params_MLE(5))) 

    print*,' parameters',params
         
    n_dist=1
    CCP=0.07d0
    
    
    
    !$omp parallel  default(private) &
    !$omp&  shared(params,f,ccp,v_fct,v_social,n_dist,v_l,mean_n,social_output,private_output,pr_d_a_n,pr_N_n_v,pr_na_v,d_rate,pr_D_model,delta_PV_model)  
    !$omp  do
    do v_l=1,villages
        print*,'village',v_l
        call compute_eq_f_ccp(params,f(:,:,:,:,:,v_l,:),ccp(:,:,:,:,v_l,:),v_fct(:,:,:,:,v_l,:),v_social(:,:,:,:,v_l,:),&
            n_dist(:,v_l),v_l,mean_n(v_l),social_output(v_l),private_output(v_l),&
            pr_d_a_n(:,:,v_l),pr_N_n_v(:,:,v_l),pr_na_v(:,:,v_l),d_rate(v_l),pr_D_model(:,v_l),delta_PV_model(v_l))
    end do
    !$omp end do  nowait
    !$omp end parallel 
    
    print*,'value of land per acre',sum(delta_PV_model*pr_v)
    
    !weight by the number of observations in each village
    do n_l=1,2;do a_l=1,types_a
        pr_d_a_n_av(a_l,n_l)=sum(pr_d_a_n(a_l,n_l,:)*pr_v_na(:,a_l,n_l))
    end do;end do    
    
    Pr_D_a_model=0.0d0
    do a_l=1,types_a;do v_l=1,villages
        Pr_D_a_model(a_l)=Pr_D_a_model(a_l)+pr_D_model(a_l,v_l)*pr_v_a(v_l,a_l)
    end do;end do
    !do v_l=1,villages
    !    print*,'village',v_l
    !    print*,'M1',pr_d_a_n(:,:,v_l)
    !    print*,'M2',pr_D_model(:,v_l)
    !end do
    if (bootstrap==0) then
        log_likelihood2=sum((pr_d_a_n_av-moment_own_nxa_data)**2.0d0/var_own_nxa)+sum((Pr_D_a_model-Pr_D_a_data)**2.0d0/Pr_D_a_var)
    else
        log_likelihood2=sum((pr_d_a_n_av-moment_own_nxa_bs(:,:,bs_l))**2.0d0/var_own_nxa)+sum((Pr_D_a_model-moments_D_bs(:,bs_l))**2.0d0/Pr_D_a_var)   
    end if
    
    
    if (isnan(log_likelihood2)) then
        print*,'paused in estimation2'
1        log_likelihood2=1000.0d0
    end if
     

    
    !GMM
    !log_likelihood=sum(((moment_own_nxa_model-moment_own_nxa_data))**2)
    if (log_likelihood2<max_mle) then
        max_mle=log_likelihood2
        if (bootstrap==0)  then
            open(unit=12, file=path_results//"parameters.txt",status='replace')
                write(12,'(f20.12,f20.12,f20.12,f20.12,f20.12)'),params,log_likelihood2
                !write(12,'(<par>f20.12)'),params,log_likelihood2
            close(12)
            OPEN(UNIT=12, FILE=path_results//"modl"//"_own_nxa_v.txt")
            write(12,*),pr_d_a_n
            close(12)  
            OPEN(UNIT=12, FILE=path_results//"modl"//"_D.txt")
            write(12,*),Pr_D_a_model
            close(12)  
            OPEN(UNIT=12, FILE=path_results//"land_values.txt")
                write(12,*),private_output,pr_v
            close(12) 
        else
            open(unit=12, file=path_results//"\bootstrap_samples\parameters"//trim(n2s)//".txt",status='replace')
                write(12,'(f20.12,f20.12,f20.12,f20.12,f20.12)'),params,log_likelihood2
                !write(12,'(<par>f20.12)'),params,log_likelihood2
            close(12)
        end if            
    end if
    

    it_est=it_est+1
    print*,'likelihood',log_likelihood2
    print*,'missing_x1',missing_x1
    
    

end function