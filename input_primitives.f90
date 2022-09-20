subroutine input_primitives()
    use primitives;use dimensions
    implicit none
    integer::P,n_l,N,v_l,u_l
    double precision,dimension(2*P_max)::PI_f
    double precision,dimension(3,villages)::rain_success_csv
    double precision,dimension(10,villages,2,unobs_types,18)::flow_fail_prob_csv 
    double precision,dimension(1)::prueba

     OPEN(UNIT=12, FILE=path_primitives//"rain_success_pr.csv")
        read(12,*),rain_success_csv
     close(12)
     OPEN(UNIT=12, FILE=path_primitives//"flow_fail_prob_r.csv")
        read(12,*),flow_fail_prob_csv 
     close(12)

     
     
    !print*,'kill heterogeneity across villages in pr of good moonzoon and success pr'
    !!pr of good moonzoon from drill_export.xls (hanan's data)
    PI_m(2,:)=rain_success_csv(2,:) 
    !PI_s: Pr. of success (varies by village)
    PI_s(1,:)=rain_success_csv(3,:) 
    
    !Discharge distribution
    q(:,1)=(/0.10d0,0.25d0,0.50d0,0.75d0,1.0d0/)

    do v_l=1,villages;do u_l=1,unobs_types
        
        !Unconditional probabilities of high and low monzoon
        PI_m(1,v_l)=1.0d0-PI_m(2,v_l)
        
        !Probability of failure (first position indicates one well, last position indicates all plots with 2 wells)

        PI_fm(1:2*P_max,1,v_l,u_l)=flow_fail_prob_csv(5,v_l,1,u_l,1:2*P_max) 
        PI_fm(1:2*P_max,2,v_l,u_l)=flow_fail_prob_csv(5,v_l,2,u_l,1:2*P_max) 

        PI_f=PI_fm(:,1,v_l,u_l)*PI_m(1,v_l)+PI_fm(:,2,v_l,u_l)*PI_m(2,v_l)
        
        
        !PI_s: Pr. of success doesn't vary with number of wells around
        PI_s(:,v_l)=PI_s(1,v_l)       
    
        !Discharge pr. depends of moonsoon and in number of wells in the adjacency
        PI_k(:,:,1,v_l,u_l)=transpose(flow_fail_prob_csv(6:10,v_l,1,u_l,1:2*P_max))
        PI_k(:,:,2,v_l,u_l)=transpose(flow_fail_prob_csv(6:10,v_l,2,u_l,1:2*P_max))

        !PI_k(:,1,1,v_l)

        !Adjust to sum exactly one
        do N=1,2*P_max
            PI_k(N,:,1,v_l,u_l)=PI_k(N,:,1,v_l,u_l)/sum(PI_k(N,:,1,v_l,u_l))
            PI_k(N,:,2,v_l,u_l)=PI_k(N,:,2,v_l,u_l)/sum(PI_k(N,:,2,v_l,u_l))
        end do

        !Fill in the probability of success/failure in vector form of dimension 2*P-1
        PI_s_v(:,:,:,v_l)=sqrt(-1.0d0)
        PI_f_v(:,:,:,v_l,u_l)=sqrt(-1.0d0)
        do P=1,P_max       
            do n_l=1,3
                !No own wells
                if (n_l==1) then
                    PI_s_v(1:P*2-1,n_l,P,v_l)=PI_s(1:P*2-1,v_l)
                end if
                !One own well
                if (n_l==2) then 
                    PI_s_v(1:P*2-1,n_l,P,v_l)=PI_s(2:P*2,v_l)
                    PI_f_v(1:P*2-1,n_l,P,v_l,u_l)=PI_f(1:2*P-1)
                end if
                !Two own wells
                if (n_l==3) then
                    PI_f_v(1:P*2-1,n_l,P,v_l,u_l)=PI_f(2:2*P)
                end if
            end do
        end do
    end do;end do

end subroutine