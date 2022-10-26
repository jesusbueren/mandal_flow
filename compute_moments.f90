subroutine compute_moments(data_in,string_name)
    use simulation; use simulation
    implicit none
    double precision,dimension(T_sim,plots_i),intent(in)::data_in
    CHARACTER (LEN=4),intent(in) ::string_name
    double precision,dimension(types_a,2)::moment_own_nxa
    double precision,dimension(wealth_quantiles,types_a)::moment_wa
    double precision,dimension(max_NFW+1)::counter_N,moment_N
    double precision,dimension(types_a,2)::counter_own_nxa
    double precision,dimension(wealth_quantiles,types_a)::counter_wa
    double precision,dimension(unobs_types,2)::counter_uhe,moment_uhe
    double precision,dimension(P_max)::counter_P,moment_P
    double precision,dimension(villages)::counter_v,moment_v
    double precision,dimension(max_NFW+1)::counter_Nbar,moment_Nbar
    double precision,dimension(max_NFW+1)::counter_N2
    double precision,dimension(3)::counter_own_n
    double precision,dimension(max_NFW+1,3)::counter_BigN_n,moment_BigN_n
    double precision,dimension(3)::counter_little_n,moment_little_n
    integer::i_l,t_l,u_l,n_l
    
    !Initialize to zero
    counter_N=0.0d0
    moment_N=0.0d0
    counter_own_nxa=0.0d0
    counter_wa=0.0d0
    moment_wa=0.0d0
    moment_own_nxa=0.0d0
    counter_uhe=0.0d0
    moment_uhe=0.0d0
    counter_P=0.0d0
    moment_P=0.0d0
    counter_v=0.0d0
    moment_v=0.0d0
    counter_Nbar=0.0d0
    moment_Nbar=0.0d0
    counter_BigN_n=0.0d0
    counter_little_n=0
    
    do i_l=1,plots_i
        if (impute_i(i_l)==0) then
            do t_l=1,T_sim
                counter_N2(modal_N(t_l,i_l))=counter_N2(modal_N(t_l,i_l))+1.0d0
                counter_own_n(n_data(t_l,i_l))=counter_own_n(n_data(t_l,i_l))+1
                if (data_in(t_l,i_l)/=-9.0d0 .and. n_data(t_l,i_l)<3 ) then

                    !Moments across number of functioning wells in adjacency
                    counter_N(modal_N(t_l,i_l))=counter_N(modal_N(t_l,i_l))+1.0d0
                    moment_N(modal_N(t_l,i_l))=(counter_N(modal_N(t_l,i_l))-1.0)/counter_N(modal_N(t_l,i_l))*moment_N(modal_N(t_l,i_l))&
                                                +1.0d0/counter_N(modal_N(t_l,i_l))*data_in(t_l,i_l)

                    !Moments across number of owned functioning wells and owned wells
                    counter_own_nxa(A_type(i_l),n_data(t_l,i_l))=counter_own_nxa(A_type(i_l),n_data(t_l,i_l))+1.0d0
                    moment_own_nxa(A_type(i_l),n_data(t_l,i_l))=(counter_own_nxa(A_type(i_l),n_data(t_l,i_l))-1.0)/counter_own_nxa(A_type(i_l),n_data(t_l,i_l))*moment_own_nxa(A_type(i_l),n_data(t_l,i_l))&
                                                    +1.0d0/counter_own_nxa(A_type(i_l),n_data(t_l,i_l))*data_in(t_l,i_l)
                    
                    !!Moments across wealth
                    if (n_data(t_l,i_l)==1) then
                        counter_wa(wealth_q(i_l),A_type(i_l))=counter_wa(wealth_q(i_l),A_type(i_l))+1.0d0
                        moment_wa(wealth_q(i_l),A_type(i_l))=(counter_wa(wealth_q(i_l),A_type(i_l))-1.0)/counter_wa(wealth_q(i_l),A_type(i_l))*moment_wa(wealth_q(i_l),A_type(i_l)) &
                                                        +1.0d0/counter_wa(wealth_q(i_l),A_type(i_l))*data_in(t_l,i_l)
                    end if

                    !Moments across number of unobserved heterogeneity types
                    do u_l=1,unobs_types
                        counter_uhe(u_l,n_data(t_l,i_l))=counter_uhe(u_l,n_data(t_l,i_l))+UHE_type(u_l,i_l)
                        moment_uhe(u_l,n_data(t_l,i_l))=moment_uhe(u_l,n_data(t_l,i_l))&
                                                        +data_in(t_l,i_l)*UHE_type(u_l,i_l)
                    end do

                        !Moments across number of plots
                        counter_P(P_type(i_l))=counter_P(P_type(i_l))+1.0d0
                        moment_P(P_type(i_l))=(counter_P(P_type(i_l))-1.0)/counter_P(P_type(i_l))*moment_P(P_type(i_l))&
                                                +1.0d0/counter_P(P_type(i_l))*data_in(t_l,i_l) 


                        !Moments across villages
                        counter_v(V_type(i_l))=counter_v(V_type(i_l))+1.0d0
                        moment_v(V_type(i_l))=(counter_v(V_type(i_l))-1.0)/counter_v(V_type(i_l))*moment_v(V_type(i_l))&
                                                +1.0d0/counter_v(V_type(i_l))*data_in(t_l,i_l)  
                        
                        !Moments across N_bar
                        counter_Nbar(NINT(N_bar(i_l)))=counter_Nbar(NINT(N_bar(i_l)))+1.0d0
                        moment_Nbar(NINT(N_bar(i_l)))=(counter_Nbar(NINT(N_bar(i_l)))-1.0)/counter_Nbar(NINT(N_bar(i_l)))*moment_Nbar(NINT(N_bar(i_l)))&
                                                +1.0d0/counter_Nbar(NINT(N_bar(i_l)))*data_in(t_l,i_l)  
                        
                end if 
                if (t_l==1) then

                    counter_BigN_n(modal_N(t_l,i_l),n_data(t_l,i_l))=counter_BigN_n(modal_N(t_l,i_l),n_data(t_l,i_l))+1
                    counter_little_n(n_data(t_l,i_l))=counter_little_n(n_data(t_l,i_l))+1

                end if
                
            end do
        end if
    end do
    
    
    
    OPEN(UNIT=12, FILE=path_results//string_name//"_N.txt")
        write(12,*),moment_N
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//"counter_N.txt")
        write(12,*),counter_N
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//string_name//"_own_nxa.txt")
        write(12,*),moment_own_nxa
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//"counter_own_nxa.txt")
        write(12,*),counter_own_nxa
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//"counter_wealth.txt")
        write(12,*),counter_wa
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//string_name//"_wealth.txt")
        write(12,*),moment_wa
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//string_name//"_uhe.txt")
        write(12,*),moment_uhe/counter_uhe
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//"counter_uhe.txt")
        write(12,*),counter_uhe
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//string_name//"_P.txt")
        write(12,*),moment_P
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//"counter_P.txt")
        write(12,*),counter_P
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//string_name//"_V.txt")
        write(12,*),moment_v
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//"counter_V.txt")
        write(12,*),counter_v
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//string_name//"_Nbar.txt")
        write(12,*),moment_Nbar
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//"counter_Nbar.txt")
        write(12,*),counter_Nbar
    close(12)
    

    
    OPEN(UNIT=12, FILE=path_results//"counter_own_n.txt")
        write(12,*),counter_own_n/sum(counter_own_n)
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//"counter_all_N.txt")
        write(12,*),counter_N2/sum(counter_N2)
    close(12)

        
    OPEN(UNIT=12, FILE=path_results//string_name//"_BigN_n.txt")
        pr_N_n_data=0.0d0
        do n_l=1,3
            pr_N_n_data(1:max_NFW+1,n_l)=counter_BigN_n(:,n_l)/sum(counter_BigN_n(:,n_l))
            write(12,*),pr_N_n_data(1:max_NFW+1,n_l) !pr_N_n_data(:,3)
        end do
    close(12)
        
    pr_little_n_data=dble(counter_little_n)/dble(sum(counter_little_n))
    OPEN(UNIT=12, FILE=path_results//string_name//"_little_n.txt")
        write(12,*),counter_little_n/sum(counter_little_n)
    close(12)
        

    moment_own_nxa_data=moment_own_nxa 
    moment_wa_data=moment_wa 
    
    
end subroutine
    
    