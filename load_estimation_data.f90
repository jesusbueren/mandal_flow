subroutine load_estimation_data
    use simulation; use primitives
    implicit none
    double precision,dimension(23,T_sim,plots_i)::data_csv
    integer::i_l,t_l,n_l,a_l,q_l
    double precision,dimension(types_a)::moment_a
    double precision,dimension(max_NFW+1)::moment_N
    double precision,dimension(2)::moment_own_n
    double precision,dimension(unobs_types)::moment_uhe
    double precision,dimension(P_max)::moment_P

    
        
    
    OPEN(UNIT=12, FILE=path_estimation//"drill_export_r.csv")
    read(12,*),data_csv
    close(12)
    
    
    V_type=data_csv(1,1,:)
    
    P_type=data_csv(2,1,:)
    wealth_i=data_csv(23,1,:)

    

    A_type=data_csv(3,1,:)

    n_data=data_csv(4,:,:)+1

    Pr_N_data=data_csv(5:15,:,:)

    UHE_type=1.0d0/dble(unobs_types)
    drilling_it(:,:,1)=data_csv(19,:,:) 
    impute_i=data_csv(20,1,:)
    can_be_zombie_i=data_csv(21,1,:)
    MNF=data_csv(22,:,:)+n_data
    print*,'no zombie in estimation data'
    impute_i=can_be_zombie_i
    
    !can_be_zombie_i=0
    impute_i=0
    


    !UHE_type(selected_type,:)=1.0d0
    
    !Modal value of unobserved heterogeneity
    do i_l=1,plots_i
        modal_UHE_type(i_l)=maxloc(UHE_type(:,i_l),1)
    end do
    
    !Modal number of functionning wells
    do i_l=1,plots_i;do t_l=1,T_sim
        modal_N(t_l,i_l)=maxloc(Pr_N_data(:,t_l,i_l),1)-(n_data(t_l,i_l)-1) !MNF(t_l,i_l)-(n_data(t_l,i_l)-1)!
    end do;end do
    N_bar=dble(sum(modal_N,1))/dble(T_sim)
    
    
    shares_n_a_v=0.0d0
    do i_l=1,plots_i;do t_l=1,T_sim
        shares_n_a_v(A_type(i_l),n_data(t_l,i_l),V_type(i_l))=shares_n_a_v(A_type(i_l),n_data(t_l,i_l),V_type(i_l))+1.0d0
    end do;end do

    shares_v=sum(sum(shares_n_a_v(:,:,:),2),1)/sum(shares_n_a_v)
    
    do n_l=1,3
        shares_v_n(:,n_l)=sum(shares_n_a_v(:,n_l,:),1)/sum(shares_n_a_v(:,n_l,:))
    end do
    shares_a_v=sum(shares_n_a_v(:,:,:),2)/sum(shares_n_a_v)

    do n_l=1,3;do a_l=1,types_a
        shares_n_a_v(a_l,n_l,:)=shares_n_a_v(a_l,n_l,:)/sum(shares_n_a_v(a_l,n_l,:))   
    end do;end do
    
    !Pdf of wealth across areas
    counter_w_a=0
    p_w_a=-9
    do i_l=1,plots_i
        counter_w_a(A_type(i_l))=counter_w_a(A_type(i_l))+1
        p_w_a(counter_w_a(A_type(i_l)),A_type(i_l))=wealth_i(i_l)
    end do
    
    share_v_wn=0.0d0
    do i_l=1,plots_i
        if (wealth_i(i_l)<w_lims(1))then 
            wealth_q(i_l)=1
        elseif (wealth_i(i_l)<w_lims(2))then 
            wealth_q(i_l)=2
        else 
            wealth_q(i_l)=3
        end if
        do t_l=1,T_sim
            share_v_wn(V_type(i_l),wealth_q(i_l),n_data(t_l,i_l),A_type(i_l))=share_v_wn(V_type(i_l),wealth_q(i_l),n_data(t_l,i_l),A_type(i_l))+1.0d0
         end do
    end do
    
    do q_l=1,wealth_quantiles; do n_l=1,3; do a_l=1,types_a
        share_v_wn(:,q_l,n_l,a_l)=share_v_wn(:,q_l,n_l,a_l)/sum(share_v_wn(:,q_l,n_l,a_l)) 
    end do;end do;end do
    

end subroutine
