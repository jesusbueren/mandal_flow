!Estimate parameters
subroutine estimation2(params_MLE,log_likeli)
    use dimensions; use primitives; use simulation
    implicit none
    double precision,dimension(par),intent(out)::params_MLE
    double precision,intent(out)::log_likeli
    double precision::log_L,ftol
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
    
    p_g(1,1)=0.4d0 
    p_g(1,2:par-1)=10.68d0
    p_g(1,par)=3.86d0


    do p_l=2,par+1
        p_g(p_l,:)=p_g(1,:)
        p_g(p_l,p_l-1)=p_g(1,p_l-1)*0.9d0
    end do

        
    !Change parameters to the (-Inf;Inf) real line
    do p_l=1,par+1
        p_g(p_l,1)=log(p_g(p_l,1)/(1.0d0-p_g(p_l,1)))
        p_g(p_l,2:par)=log(p_g(p_l,2:par))
        y(p_l)=log_likelihood2(p_g(p_l,:))  
        !read*,pause_k
    end do 

    !print*,'likelihood_ini',y(1)
        
    ftol=1.0d-7
    
    call amoeba(p_g,y,ftol,log_likelihood2,iter)
    
    log_likeli=y(1)
    p_g(:,1)=1.0d0/(1.0d0+exp(-p_g(:,1)))
    p_g(:,2:par)=exp(p_g(:,2:par))

    
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
    integer::i_l,t_l,type_l,a_l,p_l,v_l,ind,u_l,j_l,s_l,t,missing_x1,j_l2,ind2
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
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,villages,unobs_types)::F
    integer,dimension(plots_in_map,villages)::n_dist,n_dist_opt
    double precision,dimension(villages)::mean_N,social_output,private_output
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::joint_pr
    double precision,dimension(unobs_types)::expost_types
    
    params(1)=1.0d0/(1.0d0 + exp(-params_MLE(1))) 
    params(2:par)=exp(params_MLE(2:par))
        

    print*,' parameters',params
    
    log_likelihood2=0.0d0
    missing_x1=0
    
    !open(unit=12, file=path_results//"initial_beliefs.txt")
    !    read(12,*),n_dist,CCP
    !close(12)   
    n_dist=1
    CCP=0.07d0
    joint_pr=0.0d0   
    
    !$omp parallel default(private) shared(params,f,ccp,v_fct,v_social,n_dist,v_l,mean_n,social_output,private_output,joint_pr)
    !$omp  do
    do v_l=1,villages
        print*,'village',v_l
        call compute_eq_f_ccp(params,f(:,:,:,:,:,v_l,:),ccp(:,:,:,:,v_l,:),v_fct(:,:,:,:,v_l,:),v_social(:,:,:,:,v_l,:),n_dist(:,v_l),v_l,mean_n(v_l),social_output(v_l),private_output(v_l),joint_pr(:,:,:,:,v_l,:))
    end do
    !$omp end do  
    !$omp end parallel 
    

    !open(unit=12, file=path_results//"results.txt")
    !    read(12,*),F,CCP,joint_pr
    !close(12)
    !do v_l=1,villages
    !    f(:,:,:,:,:,v_l,:)=f(:,:,:,:,:,1,:)
    !    ccp(:,:,:,:,v_l,:)=ccp(:,:,:,:,1,:)
    !    joint_pr(:,:,:,:,v_l,:)=joint_pr(:,:,:,:,1,:)
    !end do
    pr_non_zombie_II=1.0d0


    
3   do s_l=1,simulations
    do i_l=1,plots_i;
        X(:,1)=(/1.0d0,N_bar(i_l),dble(A_type(i_l)-1),dble(P_type(i_l))/)
        if (impute_i(i_l)==0) then
            UHE_type_model(:,i_l)=0.0d0
            if (P_type(i_l)>1) then !more than one neighbor
                likelihood_i=1.0d0
                do t_l=1,T_sim
                    likelihood_it=0.0d0
                    av_CCP_uhe(t_l,i_l,:)=0.0d0
                    !apanyo  
                    !change this
                    j_l=maxloc(Pr_N_data(:,t_l,i_l),1)!MNF(t_l,i_l)
                    !position in the state space wrt to the CCP, PI_s_v and ,PI_f_v
                    if (n_data(t_l,i_l)==1) then
                        ind=j_l 
                    elseif (n_data(t_l,i_l)==2) then
                        ind=j_l-1
                    elseif (n_data(t_l,i_l)==3) then
                        ind=j_l-2
                    else
                        print*,'error in estimation'
                    end if
 
                    if (drilling_it(t_l,i_l,s_l)==1 .and. n_data(t_l,i_l)<3) then 
                        likelihood_it=likelihood_it+CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:)
                        av_CCP_uhe(t_l,i_l,:)=av_CCP_uhe(t_l,i_l,:)+CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:)*joint_pr(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:)
                    elseif (drilling_it(t_l,i_l,s_l)==0 .and. n_data(t_l,i_l)<3) then
                        likelihood_it=likelihood_it+(1.0d0-CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:))
                        av_CCP_uhe(t_l,i_l,:)=av_CCP_uhe(t_l,i_l,:)+CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:)
                    else
                        likelihood_it=1.0d0 
                    end if
                        
                    if (t_l<T_sim) then
                    
                        if (n_data(t_l,i_l)==1.and. drilling_it(t_l,i_l,s_l)==0 .and. n_data(t_l+1,i_l)==1 ) then
                                    
                            P_N2_N1=1.0d0
                                        
                        elseif (n_data(t_l,i_l)==1.and. drilling_it(t_l,i_l,s_l)==1 .and. n_data(t_l+1,i_l)==1 ) then
                                    
                            P_N2_N1=1.0d0-PI_s_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l))
                                    
                        elseif (n_data(t_l,i_l)==1.and. drilling_it(t_l,i_l,s_l)==1 .and. n_data(t_l+1,i_l)==2 ) then
                                    
                            P_N2_N1=PI_s_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l))
                                    
                        elseif (n_data(t_l,i_l)==2.and. drilling_it(t_l,i_l,s_l)==0 .and. n_data(t_l+1,i_l)==1 ) then
                                    
                            P_N2_N1=PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:)
                                    
                        elseif (n_data(t_l,i_l)==2.and. drilling_it(t_l,i_l,s_l)==0 .and. n_data(t_l+1,i_l)==2 ) then
                                    
                            P_N2_N1=1.0d0-PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:)
                                    
                        elseif (n_data(t_l,i_l)==2.and. drilling_it(t_l,i_l,s_l)==1 .and. n_data(t_l+1,i_l)==1 ) then
                                    
                            P_N2_N1=(1.d0-PI_s_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l)))*PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:)
                                    
                        elseif (n_data(t_l,i_l)==2.and. drilling_it(t_l,i_l,s_l)==1 .and. n_data(t_l+1,i_l)==2 ) then 
                                    
                            P_N2_N1=(1.d0-PI_s_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l)))*(1.0d0-PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:))+&
                                                PI_s_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l))*PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:)
                                
                        elseif (n_data(t_l,i_l)==2.and. drilling_it(t_l,i_l,s_l)==1 .and. n_data(t_l+1,i_l)==3 ) then
                                    
                            P_N2_N1=PI_s_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l))*(1.0d0-PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:))
                                    
                        elseif (n_data(t_l,i_l)==3 .and. n_data(t_l+1,i_l)==1 ) then
                                    
                            P_N2_N1=PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:)**2.0d0
                                
                        elseif (n_data(t_l,i_l)==3  .and. n_data(t_l+1,i_l)==2 ) then !n_data(:,i_l) 
                                    
                            P_N2_N1=2.0d0*PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:)*(1.0d0-PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:))
                                    
                        elseif (n_data(t_l,i_l)==3 .and. n_data(t_l+1,i_l)==3 ) then
                                    
                            P_N2_N1=(1.0d0-PI_f_v(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:))**2.0d0
                                    
                        end if
                            
                        !apanyo
                        !change this
                        j_l2=maxloc(Pr_N_data(:,t_l+1,i_l),1)!MNF(t_l+1,i_l)
                        if (n_data(t_l+1,i_l)==1) then
                            ind2=j_l2 !position in the state space wrt to the CCP, PI_s_v and ,PI_f_v
                        elseif (n_data(t_l+1,i_l)==2) then
                            ind2=j_l2-1
                        elseif (n_data(t_l+1,i_l)==3) then
                            ind2=j_l2-2
                        else
                            print*,'error in estimation'
                        end if
                     
                    
                        P_BigN2_BigN1=F(ind,ind2,n_data(t_l,i_l),n_data(t_l+1,i_l),P_type(i_l),V_type(i_l),:) 
                        !if (sum(P_BigN2_BigN1)==0)then
                        !    print*,''
                        !end if
                        !if (sum(P_N2_N1)==0)then
                        !    print*,''
                        !end if
                        !if (sum(likelihood_it)==0)then
                        !    print*,''
                        !end if
                            
                        likelihood_it=likelihood_it*P_N2_N1*P_BigN2_BigN1

                    end if
                    


                    if (t_l==1) then
                        !change this
                        !Full likelihood
                        likelihood_it=likelihood_it*joint_pr(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:) 
                        !Cond. likelihood 
                        !print*,joint_pr(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:)/sum(joint_pr(ind,n_data(t_l,i_l),P_type(i_l),V_type(i_l),:)) !joint_pr(2,2,3,1,:) joint_pr(9,2,5,1,:) 
                        !likelihood_it=likelihood_it*joint_pr(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:)/sum(joint_pr(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:))
                        if (isnan(sum(likelihood_it))) then
                            likelihood_it=0.0d0
                        end if
                    end if  
                    

                    if (sum(likelihood_it)==0)then
                        print*,'zero likelihood',i_l,t_l
                    end if

                    likelihood_i=likelihood_i*likelihood_it


                    !if (isnan(sum(likelihood_i))) then
                    !    print*,'pb in likelihood',i_l,ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),2),drilling_it(t_l,i_l,s_l)
                    !    read*,end_k
                    !end if
                    !if (sum(likelihood_i)==0.0d0) then
                    !    print*,'pb in likelihood',i_l,ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),2),drilling_it(t_l,i_l,s_l)
                    !    read*,end_k
                    !end if
                        
                end do;
                
                if (sum(likelihood_i)==0)then
                        print*,'zero likelihood',i_l
                    end if
                
                posterior_type(i_l,:)=likelihood_i*pr_unobs_t/sum(likelihood_i*pr_unobs_t)  !posterior_type(:,1)
                

                !Model 5
                
                log_likelihood2=log_likelihood2+log(sum(likelihood_i*pr_unobs_t)*pr_non_zombie_II) 
                !log_likelihood2=log_likelihood2+log(sum(likelihood_i*pr_unobs_t)*pr_non_zombie_II) 
                


                do t_l=1,T_sim
                    if ((drilling_it(t_l,i_l,s_l)==1 .or. drilling_it(t_l,i_l,s_l)==0) .and. sum(likelihood_i*pr_unobs_t)/=0.0d0 .and. n_data(t_l,i_l)<3) then 
                        if (can_be_zombie_i(i_l)==0) then
                            av_CCP_it(t_l,i_l)=sum(av_CCP_uhe(t_l,i_l,:)*(likelihood_i*pr_unobs_t/sum(likelihood_i*pr_unobs_t))) 
                        else
                            av_CCP_it(t_l,i_l)=sum(av_CCP_uhe(t_l,i_l,:)*(likelihood_i*pr_unobs_t*pr_non_zombie_II/sum(likelihood_i*pr_unobs_t*pr_non_zombie_II+(1-pr_non_zombie_II))))
                        end if    
                    else
                        av_CCP_it(t_l,i_l)=-9.0d0
                    end if
                end do
                        
            end if
        end if
    end do
    end do
        !close(12)

    log_likelihood2=-log_likelihood2
    
    call compute_moments(av_CCP_it,"modl",joint_pr,moment_own_nxa_model)
    log_likelihood2=sum((moment_own_nxa_model-moment_own_nxa_data)**2.0d0)

    
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
        
        call compute_moments(av_CCP_it,"modl",joint_pr,moment_own_nxa_model)
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