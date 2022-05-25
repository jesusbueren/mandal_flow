!subroutine bootstrap_se()
!    use simulation
!    implicit none
!    double precision,dimension(20,T_sim,plots_i)::data_csv
!    integer,dimension(plots_i)::V_type_or,P_type_or,A_type_or
!    double precision,dimension(unobs_types,plots_i)::UHE_type_or
!    integer,dimension(T_sim,plots_i,simulations)::drilling_it_or
!    integer,dimension(T_sim,plots_i)::n_data_or
!    double precision,dimension(max_NFW+1,T_sim,plots_i)::Pr_N_data_or
!    integer::i_l,t_l,i_sim,bs_l 
!    double precision::u,likeli
!    integer,parameter::bs_samples=500
!    double precision,dimension(par,500)::params_bs
!    
!    !make it a one for not saving the model implied moments
!    bootstrap=1
!    
!    OPEN(UNIT=12, FILE=path_estimation//"drill_export_r.csv")
!    read(12,*),data_csv
!    close(12)
!    
!    V_type_or=data_csv(1,1,:)
!    P_type_or=data_csv(2,1,:)
!    A_type_or=data_csv(3,1,:)
!    n_data_or=data_csv(4,:,:)+1
!    Pr_N_data_or=data_csv(5:15,:,:)
!    UHE_type_or=0.25d0!data_csv(16:18,1,:)
!    drilling_it_or(:,:,1)=data_csv(19,:,:)
!    
!    do bs_l=1,bs_samples
!        print*,'Sample ',bs_l,' out of ',bs_samples
!        do i_l=1,plots_i
!            call RANDOM_NUMBER(u)
!            i_sim=ceiling(u*dble(plots_i))
!            V_type(i_l)=V_type_or(i_sim)
!            P_type(i_l)=P_type_or(i_sim)
!            A_type(i_l)=A_type_or(i_sim)
!            n_data(:,i_l)=n_data_or(:,i_sim)
!            Pr_N_data(:,:,i_l)=Pr_N_data_or(:,:,i_sim)
!            UHE_type(:,i_l)=UHE_type_or(:,i_sim)
!            drilling_it(:,i_l,1)=drilling_it_or(:,i_sim,1)
!        end do
!    
!        !Modal value of unobserved heterogeneity
!        do i_l=1,plots_i
!            modal_UHE_type(i_l)=maxloc(UHE_type(:,i_l),1)
!        end do
!    
!        !Modal number of functionning wells
!        do i_l=1,plots_i;do t_l=1,T_sim
!            modal_N(t_l,i_l)=maxloc(Pr_N_data(:,t_l,i_l),1)
!        end do;end do
!    
!        call estimation(params_bs(:,bs_l),likeli)
!        OPEN(UNIT=12, FILE=path_results//"bootstrapped_parameters.txt",access='append')
!        write(12,'(F20.12,F20.12,F20.12,F20.12)'),params_bs(1,bs_l),params_bs(2,bs_l),likeli
!        close(12)
!    
!    end do
!    
!    
!    
!end subroutine
