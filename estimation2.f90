!Estimate parameters
subroutine estimation2(params_MLE,log_likeli)
    use dimensions; use primitives; use simulation
    implicit none
    double precision,dimension(par),intent(out)::params_MLE
    double precision,intent(out)::log_likeli
    double precision::log_L,ftol,fret
    double precision,dimension(par+1,par)::p_g
    double precision,dimension(par+1)::y
    integer::iter,p_l,a_l,P_l2,v_l,ind,n_l,u_l
    interface 
        function log_likelihood2(params_MLE)
            use dimensions
            double precision,dimension(par),intent(in)::params_MLE
            double precision::log_likelihood2
        end function log_likelihood2
    end interface
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP_old,CCP_mid
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::V_fct,V_social
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::Ef_v 
    double precision::dist
    integer::it
    integer(8),dimension(2*P_max-1,3,3,P_max,villages,unobs_types)::iterations_all
    double precision,dimension(par,par)::xi
    integer,dimension(1)::seed_c
    double precision, dimension(villages)::mean_N,mean_NPV,mean_budget
    double precision,dimension(villages)::village_fe
    character::pause_k
    
      
    call random_seed(PUT=seed_c)
    !Fixing beliefs, estimate parameter
    !print*,'Initial Conditions'
    
    max_mle=99999999.0d0
    
    !share water (CES)
    p_g(1,1)=0.727d0 
    !Productivity p
    p_g(1,2)=18.878d0
    !Var taste shock
    p_g(1,3)=1.298d0
    !Intercept logit contrained
    p_g(1,4)=-22.338d0
    !Shrinkage parameter 
    p_g(1,5)=2.521d0  
    !Spatial correlation
    !p_g(1,6)=0.95d0 
    !wealth effect logit constrained
    p_g(1,6)=1.947d0
    !flow p
    p_g(1,7)=0.18d0



    do p_l=2,par+1
        p_g(p_l,:)=p_g(1,:)
        p_g(p_l,p_l-1)=p_g(1,p_l-1)*0.9d0
    end do

    xi=0.0d0    
    !Change parameters to the (-Inf;Inf) real line
    do p_l=1,par+1
        p_g(p_l,1)=log(p_g(p_l,1)/(1.0d0-p_g(p_l,1)))
        p_g(p_l,2:3)=log(p_g(p_l,2:3))
        !p_g(p_l,6)=log(p_g(p_l,6)/(1.0d0-p_g(p_l,6)))
        p_g(p_l,7)=log(p_g(p_l,7))
        y(p_l)=log_likelihood2(p_g(p_l,:))  
        if (p_l<par+1)then
            xi(p_l,p_l)=1.0d0
        end if
        !read*,pause_k
    end do 

    !print*,'likelihood_ini',y(1)
        
    ftol=1.0d-3
    call amoeba(p_g,y,ftol,log_likelihood2,iter)
    call powell(p_g(1,:),xi,ftol,iter,fret)
    
    log_likeli=y(1)
    p_g(1,1)=1.0d0/(1.0d0+exp(-p_g(1,1)))
    p_g(1,2:3)=exp(p_g(1,2:3))
    !p_g(1,6)=1.0d0/(1.0d0+exp(-p_g(1,6)))
    p_g(1,7)=exp(p_g(1,7))



    
    params_MLE=p_g(1,:)
    print*,'likelihood amoeba',y(1)
    
    
    if (dist<0.0) then
        print*,'error in precision of iterations'   
    end if

   
end subroutine
    
function log_likelihood2(params_MLE)
    use dimensions; use simulation; use cadastral_maps; use primitives
    implicit none
    double precision,dimension(par),intent(in)::params_MLE
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
    double precision,dimension(villages)::mean_N,social_output,private_output
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::joint_pr
    double precision,dimension(unobs_types)::expost_types
    double precision,dimension(types_a,2,villages)::pr_d_a_n
    double precision,dimension(types_a,2)::pr_d_a_n_av
    double precision,dimension(types_a,2)::pr_d_a_av
    double precision::fraction_constrained=0.7d0
    double precision,dimension(2*P_max-1,3,villages)::pr_N_n_v
    double precision,dimension(3,villages)::pr_n_v
    double precision,dimension(3,types_a,villages)::pr_na_v,pr_non_zombie_anv
    double precision,dimension(2*P_max-1,3)::pr_N_n_av
    double precision,dimension(3)::pr_n_av
    integer(8),dimension(1)::seed=321
    double precision::u
    double precision,dimension(plots_i):: pr_non_zombie_i
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::empty=-9.0d0
    double precision,dimension(types_a)::pr_non_zombie_a_s
    double precision,dimension(wealth_quantiles,types_a)::pr_non_zombie_aw_s
    integer,dimension(types_a)::counter_a
    integer,dimension(wealth_quantiles,types_a)::counter_aw
    double precision,dimension(wealth_quantiles,types_a)::pr_d_wa
    double precision,dimension(wealth_quantiles,types_a,villages)::pr_non_zombie_wav
    
    
    params(1)=1.0d0/(1.0d0 + exp(-params_MLE(1))) 
    params(2:3)=exp(params_MLE(2:3))
    params(4)=params_MLE(4)
    params(5)=params_MLE(5)
    !params(6)=1.0d0/(1.0d0 + exp(-params_MLE(6))) 
    params(6)=params_MLE(6) 
    params(7)=exp(params_MLE(7))

    !rho_sc=params(6)
    shrinkage_p=params(5)
    logit_constrain_p=(/params(4),params(6)/)
    
    counter_a=0
    pr_non_zombie_a_s=0.0d0
    counter_aw=0
    pr_non_zombie_aw_s=0.0d0
    do i_l=1,plots_i
        pr_non_zombie_i(i_l)=1.0d0/(1.0d0+exp(-(logit_constrain_p(1)+logit_constrain_p(2)*wealth_i(i_l))))
        counter_a(A_type(i_l))=counter_a(A_type(i_l))+1
        pr_non_zombie_a_s(A_type(i_l))=pr_non_zombie_a_s(A_type(i_l))+pr_non_zombie_i(i_l)
        counter_aw(wealth_q(i_l),A_type(i_l))=counter_aw(wealth_q(i_l),A_type(i_l))+1
        pr_non_zombie_aw_s(wealth_q(i_l),A_type(i_l))=pr_non_zombie_aw_s(wealth_q(i_l),A_type(i_l))+pr_non_zombie_i(i_l)
    end do
    pr_non_zombie_a_s=pr_non_zombie_a_s/dble(counter_a)
    pr_non_zombie_aw_s=pr_non_zombie_aw_s/dble(counter_aw) 
    
    if (1.0d0/(1.0d0+exp(-(logit_constrain_p(1)+logit_constrain_p(2)*(12.97-shrinkage_p))))<0.1d0) then
        go to 1
    end if
     
    

    print*,'Non-constrained pr sample',pr_non_zombie_a_s

    
    impute_i=0
     
    
    call load_cadastral_maps()
    
    !Call seed number
    call random_seed(PUT=seed)
    

    print*,' parameters',params
    
    log_likelihood2=0.0d0
    missing_x1=0
    
    !open(unit=12, file=path_results//"initial_beliefs.txt")
    !    read(12,*),n_dist,CCP
    !close(12)   
    n_dist=1
    CCP=0.07d0
    joint_pr=0.0d0   
    
    !$omp parallel default(private) shared(params,f,ccp,v_fct,v_social,n_dist,v_l,mean_n,social_output,private_output,joint_pr,pr_d_a_n,pr_N_n_v,pr_na_v) 
    !$omp  do
    do v_l=1,villages
        print*,'village',v_l
        call compute_eq_f_ccp(params,f(:,:,:,:,:,v_l,:),ccp(:,:,:,:,v_l,:),v_fct(:,:,:,:,v_l,:),v_social(:,:,:,:,v_l,:),n_dist(:,v_l), &
            v_l,mean_n(v_l),social_output(v_l),private_output(v_l),joint_pr(:,:,:,:,v_l,:),pr_d_a_n(:,:,v_l),pr_N_n_v(:,:,v_l),pr_na_v(:,:,v_l))
    end do
    !$omp end do  nowait
    !$omp end parallel 
    
    open(unit=12, file=path_results//"checks.txt")
        write(12,*),pr_d_a_n,pr_N_n_v,pr_na_v !pr_d_a_n(:,1,7) pr_na_v(1,1,:) pr_non_zombie_anv(1,:,1) 
    close(12)  
    
    
    do a_l=1,types_a;do v_l=1,villages    
        pr_non_zombie_anv(1,a_l,v_l)=pr_na_v(1,a_l,v_l)*pr_non_zombie_a_s(a_l)/(pr_na_v(1,a_l,v_l)*pr_non_zombie_a_s(a_l)+(1.0d0-pr_non_zombie_a_s(a_l)))
    end do ;end do
    pr_non_zombie_anv(2:3,:,:)=1.0d0
    
    do a_l=1,types_a;do v_l=1,villages; do q_l=1,wealth_quantiles  
        pr_non_zombie_wav(q_l,a_l,v_l)=pr_na_v(1,a_l,v_l)*pr_non_zombie_aw_s(q_l,a_l)/(pr_na_v(1,a_l,v_l)*pr_non_zombie_aw_s(q_l,a_l)+(1.0d0-pr_non_zombie_aw_s(q_l,a_l)))       
    end do; end do; end do
     
    open(unit=12, file=path_results//"results.txt")
        write(12,*),F,CCP,joint_pr,pr_d_a_n
    close(12)
    !do v_l=1,villages
    !    f(:,:,:,:,:,v_l,:)=f(:,:,:,:,:,9,:)
    !    ccp(:,:,:,:,v_l,:)=ccp(:,:,:,:,9,:)
    !    joint_pr(:,:,:,:,v_l,:)=joint_pr(:,:,:,:,9,:)
    !    pr_d_a_n(:,:,v_l)=pr_d_a_n(:,:,9)
    !    pr_n_v(:,v_l)=pr_n_v(:,9) !pr_n_v(:,1)
    !end do
    pr_non_zombie_II=1.0d0


    
!3   do s_l=1,simulations
!    do i_l=1,plots_i;
!        X(:,1)=(/1.0d0,N_bar(i_l),dble(A_type(i_l)-1),dble(P_type(i_l))/)
!        if (impute_i(i_l)==0) then
!            UHE_type_model(:,i_l)=0.0d0
!            if (P_type(i_l)>1) then !more than one neighbor
!                likelihood_i=1.0d0
!                do t_l=1,T_sim
!                    likelihood_it=0.0d0
!                    av_CCP_uhe(t_l,i_l,:)=0.0d0
!                    !apanyo  
!                    !change this
!                    j_l=maxloc(Pr_N_data(:,t_l,i_l),1)!MNF(t_l,i_l)
!                    !position in the state space wrt to the CCP, PI_s_v and ,PI_f_v
!                    if (n_data(t_l,i_l)==1) then
!                        ind=j_l 
!                    elseif (n_data(t_l,i_l)==2) then
!                        ind=j_l-1
!                    elseif (n_data(t_l,i_l)==3) then
!                        ind=j_l-2
!                    else
!                        print*,'error in estimation'
!                    end if
! 
!                    if (drilling_it(t_l,i_l,s_l)==1 .and. n_data(t_l,i_l)<3) then 
!                        likelihood_it=likelihood_it+CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:)
!                        av_CCP_uhe(t_l,i_l,:)=av_CCP_uhe(t_l,i_l,:)+CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:)*joint_pr(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:)
!                    elseif (drilling_it(t_l,i_l,s_l)==0 .and. n_data(t_l,i_l)<3) then
!                        likelihood_it=likelihood_it+(1.0d0-CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:))
!                        av_CCP_uhe(t_l,i_l,:)=av_CCP_uhe(t_l,i_l,:)+CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:)
!                    else
!                        likelihood_it=1.0d0 
!                    end if
!                        
!                    if (t_l<T_sim) then
!                    
!                        if (n_data(t_l,i_l)==1.and. drilling_it(t_l,i_l,s_l)==0 .and. n_data(t_l+1,i_l)==1 ) then
!                                    
!                            P_N2_N1=1.0d0
!                                        
!                        elseif (n_data(t_l,i_l)==1.and. drilling_it(t_l,i_l,s_l)==1 .and. n_data(t_l+1,i_l)==1 ) then
!                                    
!                            P_N2_N1=1.0d0-PI_s_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l))
!                                    
!                        elseif (n_data(t_l,i_l)==1.and. drilling_it(t_l,i_l,s_l)==1 .and. n_data(t_l+1,i_l)==2 ) then
!                                    
!                            P_N2_N1=PI_s_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l))
!                                    
!                        elseif (n_data(t_l,i_l)==2.and. drilling_it(t_l,i_l,s_l)==0 .and. n_data(t_l+1,i_l)==1 ) then
!                                    
!                            P_N2_N1=PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:)
!                                    
!                        elseif (n_data(t_l,i_l)==2.and. drilling_it(t_l,i_l,s_l)==0 .and. n_data(t_l+1,i_l)==2 ) then
!                                    
!                            P_N2_N1=1.0d0-PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:)
!                                    
!                        elseif (n_data(t_l,i_l)==2.and. drilling_it(t_l,i_l,s_l)==1 .and. n_data(t_l+1,i_l)==1 ) then
!                                    
!                            P_N2_N1=(1.d0-PI_s_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l)))*PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:)
!                                    
!                        elseif (n_data(t_l,i_l)==2.and. drilling_it(t_l,i_l,s_l)==1 .and. n_data(t_l+1,i_l)==2 ) then 
!                                    
!                            P_N2_N1=(1.d0-PI_s_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l)))*(1.0d0-PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:))+&
!                                                PI_s_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l))*PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:)
!                                
!                        elseif (n_data(t_l,i_l)==2.and. drilling_it(t_l,i_l,s_l)==1 .and. n_data(t_l+1,i_l)==3 ) then
!                                    
!                            P_N2_N1=PI_s_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l))*(1.0d0-PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:))
!                                    
!                        elseif (n_data(t_l,i_l)==3 .and. n_data(t_l+1,i_l)==1 ) then
!                                    
!                            P_N2_N1=PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:)**2.0d0
!                                
!                        elseif (n_data(t_l,i_l)==3  .and. n_data(t_l+1,i_l)==2 ) then !n_data(:,i_l) 
!                                    
!                            P_N2_N1=2.0d0*PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:)*(1.0d0-PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:))
!                                    
!                        elseif (n_data(t_l,i_l)==3 .and. n_data(t_l+1,i_l)==3 ) then
!                                    
!                            P_N2_N1=(1.0d0-PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:))**2.0d0
!                                    
!                        end if
!                            
!                        !apanyo
!                        !change this
!                        j_l2=maxloc(Pr_N_data(:,t_l+1,i_l),1)!MNF(t_l+1,i_l)
!                        if (n_data(t_l+1,i_l)==1) then
!                            ind2=j_l2 !position in the state space wrt to the CCP, PI_s_v and ,PI_f_v
!                        elseif (n_data(t_l+1,i_l)==2) then
!                            ind2=j_l2-1
!                        elseif (n_data(t_l+1,i_l)==3) then
!                            ind2=j_l2-2
!                        else
!                            print*,'error in estimation'
!                        end if
!                     
!                    
!                        P_BigN2_BigN1=F(ind,ind2,n_data(t_l,i_l),n_data(t_l+1,i_l),P_type(i_l),V_type(i_l),:) 
!                        !if (sum(P_BigN2_BigN1)==0)then
!                        !    print*,''
!                        !end if
!                        !if (sum(P_N2_N1)==0)then
!                        !    print*,''
!                        !end if
!                        !if (sum(likelihood_it)==0)then
!                        !    print*,''
!                        !end if
!                            
!                        likelihood_it=likelihood_it*P_N2_N1*P_BigN2_BigN1
!
!                    end if
!                    
!
!
!                    if (t_l==1) then
!                        !change this
!                        !Full likelihood
!                        likelihood_it=likelihood_it*joint_pr(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:) 
!                        !Cond. likelihood 
!                        !print*,joint_pr(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:)/sum(joint_pr(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:)) !joint_pr(2,2,3,1,:) joint_pr(9,2,5,1,:) 
!                        !likelihood_it=likelihood_it*joint_pr(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:)/sum(joint_pr(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:))
!                        if (isnan(sum(likelihood_it))) then
!                            likelihood_it=0.0d0
!                        end if
!                    end if  
!                    
!
!                    if (sum(likelihood_it)==0)then
!                        print*,'zero likelihood',i_l,t_l
!                    end if
!
!                    likelihood_i=likelihood_i*likelihood_it
!
!
!                    !if (isnan(sum(likelihood_i))) then
!                    !    print*,'pb in likelihood',i_l,ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),2),drilling_it(t_l,i_l,s_l)
!                    !    read*,end_k
!                    !end if
!                    !if (sum(likelihood_i)==0.0d0) then
!                    !    print*,'pb in likelihood',i_l,ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),2),drilling_it(t_l,i_l,s_l)
!                    !    read*,end_k
!                    !end if
!                        
!                end do;
!                
!                if (sum(likelihood_i)==0)then
!                        print*,'zero likelihood',i_l
!                    end if
!                
!                posterior_type(i_l,:)=likelihood_i*pr_unobs_t/sum(likelihood_i*pr_unobs_t)  !posterior_type(:,1)
!                
!
!                !Model 5
!                
!                log_likelihood2=log_likelihood2+log(sum(likelihood_i*pr_unobs_t)*pr_non_zombie_II) 
!                !log_likelihood2=log_likelihood2+log(sum(likelihood_i*pr_unobs_t)*pr_non_zombie_II) 
!                
!
!
!                do t_l=1,T_sim
!                    if ((drilling_it(t_l,i_l,s_l)==1 .or. drilling_it(t_l,i_l,s_l)==0) .and. sum(likelihood_i*pr_unobs_t)/=0.0d0 .and. n_data(t_l,i_l)<3) then 
!                        if (can_be_zombie_i(i_l)==0) then
!                            av_CCP_it(t_l,i_l)=sum(av_CCP_uhe(t_l,i_l,:)*(likelihood_i*pr_unobs_t/sum(likelihood_i*pr_unobs_t))) 
!                        else
!                            av_CCP_it(t_l,i_l)=sum(av_CCP_uhe(t_l,i_l,:)*(likelihood_i*pr_unobs_t*pr_non_zombie_II/sum(likelihood_i*pr_unobs_t*pr_non_zombie_II+(1-pr_non_zombie_II))))
!                        end if    
!                    else
!                        av_CCP_it(t_l,i_l)=-9.0d0
!                    end if
!                end do
!                        
!            end if
!        end if
!    end do
!    end do
!        !close(12)
!
!    log_likelihood2=-log_likelihood2
    
    
    do n_l=1,2;do a_l=1,types_a
        pr_d_a_n_av(a_l,n_l)=sum(pr_d_a_n(a_l,n_l,:)*pr_non_zombie_anv(n_l,a_l,:)*shares_n_a_v(a_l,n_l,:)) 
    end do;end do
    
    do q_l=1,wealth_quantiles;do a_l=1,types_a
        pr_d_wa(q_l,a_l)=sum(pr_d_a_n(a_l,1,:)*pr_non_zombie_wav(q_l,a_l,:)*share_v_wn(:,q_l,1,a_l))
    end do; end do
    
    do n_l=1,3;do ind=1,2*P_max-1
        pr_N_n_av(ind,n_l)=sum(pr_N_n_v(ind,n_l,:)*shares_v_n(:,n_l)) !pr_N_n_v(:,2,4)
    end do;end do
    
    do n_l=1,3
        pr_n_av(n_l)=sum(pr_na_v(n_l,:,:)*shares_a_v) !pr_na_v(:,4,1)
        if (isnan(sum(pr_N_n_av(1:max_NFW+1,n_l)))) then
        pr_N_n_av(1:max_NFW+1,n_l)=1.0d0/dble(max_NFW+1)
        print*,'some n not visited',n_l
    end if
    end do
    
    pr_n_av=pr_n_av*sum(pr_non_zombie_i(:))/dble(plots_i)
    pr_n_av(1)=pr_n_av(1)+(1.0d0-sum(pr_non_zombie_i(:))/dble(plots_i))
    
    
    print*,sum(pr_N_n_av(1,:)*pr_n_av),sum(pr_N_n_data(1,:)*pr_little_n_data)

    log_likelihood2=sum((pr_d_a_n_av-moment_own_nxa_data)**2.0d0/var_own_nxa)+(pr_n_av(1)-pr_little_n_data(1))**2.0d0/var_little_n(1) + &
                     sum((pr_d_wa-moment_wa_data)**2.0d0/var_wa_data) + (sum(pr_N_n_av(1,:)*pr_n_av)-sum(pr_N_n_data(1,:)*pr_little_n_data))**2.0d0/var_N(1)
    
    if (isnan(log_likelihood2)) then
        print*,'paused in estimation2'
        print*,pr_N_n_v(1,:,:)
1        log_likelihood2=1000.0d0
        !read*, end_k
    end if
    


    
    OPEN(UNIT=12, FILE=path_results//"modl"//"_own_nxa.txt")
        write(12,*),pr_d_a_n_av
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//"modl"//"_wealth.txt")
        write(12,*),pr_d_wa
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//"modl"//"_BigN_n2.txt")
        write(12,*),pr_N_n_av
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//"modl"//"_n.txt")
        write(12,*),pr_n_av
    close(12)

    
    !GMM
    !log_likelihood=sum(((moment_own_nxa_model-moment_own_nxa_data))**2)
    if (bootstrap==0 .and. log_likelihood2<max_mle) then
        open(unit=12, file=path_results//"parameters.txt",status='replace')
            write(12,'(<par>f20.12,f20.12)'),params,log_likelihood2
        close(12)
        open(unit=12, file=path_results//"posterior_type.txt",status='replace')
            write(12,*),posterior_type
        close(12)
        !open(unit=12, file=path_results//"initial_beliefs.txt",status='replace')
        !write(12,*),n_dist,CCP
        !close(12)

        max_mle=log_likelihood2
        !n_dist_opt=n_dist
        !CCP_opt=CCP
    end if
    

    it_est=it_est+1
    print*,'likelihood',log_likelihood2
    print*,'missing_x1',missing_x1
    
    !open(unit=12, file=path_results//"welfare_in_sample.txt",status='replace')
    !    write(12,'(<par>f20.12,f20.12)'),params,log_likelihood2
    !close(12)
    

end function