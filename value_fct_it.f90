subroutine value_fct_it(Ef_v,F,P,CCP,v_l,u_l,V_new)
    use dimensions;use primitives
    implicit none
    integer,intent(in)::P
    double precision,dimension(2*P-1,3),intent(in)::Ef_v
    double precision,dimension(2*P-1,2*P-1,3,3),intent(in)::F
    integer,intent(in)::v_l,u_l
    double precision,dimension(2*P-1,3)::V_old
    double precision,dimension((2*P-1)*3)::Vec_old,Vec_new
    double precision::crit,dist
    double precision,dimension(2*P-1,2),intent(out)::CCP
    double precision,dimension(2*P-1,3),intent(out)::V_new
    character::pause_k
    integer::it,it_max=8000

    !Define intial guess of value function
    V_old=0.0d0
    dist=1.0d0
    crit=1.0d-12
    it=0
    CCP=0.0d0
        
    do while (dist>crit .and. it<it_max)
        call one_step_value_fct_it(Ef_v,F,P,CCP,v_l,u_l,V_old,V_new)
        it=it+1       
        !Check contraction
        !!!!!!!!!!!!!!!!!!
            Vec_old=reshape(V_old,(/(2*P-1)*3/))
            Vec_new=reshape(V_new,(/(2*P-1)*3/))
            dist=maxval(abs(Vec_old-Vec_new))
            if (dist==1.0d0/0.0d0)then
                print*,Vec_old
                read*,pause_k
            end if
            if (it==it_max)then
                print*,'no convergence in vfi',dist
            end if
            if (dist>crit)then
                V_old=V_new
            end if
    end do  
    if (it==it_max)then
        print*,'no convergence in vfi',dist
    end if
    
    !do ind=2,2*P-1
    !    if (V_new(ind-1,1)<V_new(ind,1)) then
    !        print*,'hi there 1'
    !    end if
    !    if (V_new(ind-1,2)<V_new(ind,2)) then
    !        print*,'hi there 2'
    !    end if
    !    if (V_new(ind-1,3)<V_new(ind,3)) then
    !        print*,'hi there 3'
    !    end if
    !    !F(3,:,3,1)
    !end do
end subroutine
    
        
subroutine one_step_value_fct_it(Ef_v,F,P,CCP,v_l,u_l,V_old,V_new)
    use dimensions;use primitives
    implicit none
    integer,intent(in)::P
    double precision,dimension(2*P-1,3),intent(in)::Ef_v
    double precision,dimension(2*P-1,2*P-1,3,3),intent(in)::F
    integer,intent(in)::v_l,u_l
    double precision,dimension(2*P-1,3),intent(in)::V_old
    double precision,dimension(2*P-1,3),intent(out)::V_new
    double precision,dimension((2*P-1)*3)::Vec_old,Vec_new
    double precision,dimension(2*P-1)::v_00,v_0I,v_10,v_1I,vec_N
    double precision,dimension(2*P-1,2),intent(inout)::CCP
    integer::ind
    character::pause_k
    
    CCP=max(CCP,1d-18)
    
    V_new=0.0d0
    
    do ind=1,2*P-1
        vec_N(ind)=dble(ind)-1.0d0
    end do

        
            !No well (n=1)
            !!!!!!!!!!!!!!
                !No attempt
                v_00(1:2*P-1)=T_g+beta*matmul(F(1:2*P-1,1:2*P-1,1,1),V_old(1:2*P-1,1))
                !Attempt               
                v_0I(1:2*P-1)=T_g+PI_s_v(1:2*P-1,1,P,v_l)*(-c_s+beta*(matmul(F(1:2*P-1,1:2*P-1,1,2),V_old(1:2*P-1,2)))) & !success
                           +(1.0d0-PI_s_v(1:2*P-1,1,P,v_l))*(-c_d+beta*(matmul(F(1:2*P-1,1:2*P-1,1,1),V_old(1:2*P-1,1)))) !Failure
                !Value function
                V_new(1:2*P-1,1)=rho(1)*log(1.0d0+exp((v_0I(1:2*P-1)-v_00(1:2*P-1))/rho(1)))+v_00(1:2*P-1)+rho(1)*gamma

            !!One well (n=2)
            !!!!!!!!!!!!!!!!!
                !No attempt
                v_10(1:2*P-1)=T_g+Ef_v(1:2*P-1,2)-tau-vec_N*tau_per_N &
                    +beta*((1.0d0-PI_f_v(1:2*P-1,2,P,v_l,u_l))*matmul(F(1:2*P-1,1:2*P-1,2,2),V_old(1:2*P-1,2))) & !No failure
                    +beta*(PI_f_v(1:2*P-1,2,P,v_l,u_l)*matmul(F(1:2*P-1,1:2*P-1,2,1),V_old(1:2*P-1,1)))  !Failure

    
                !Attempt
                v_1I(1:2*P-1)=T_g+Ef_v(1:2*P-1,2)-tau-vec_N*tau_per_N &
                        -PI_s_v(1:2*P-1,2,P,v_l)*c_s-(1.0d0-PI_s_v(1:2*P-1,2,P,v_l))*c_d &
                        +beta*(PI_s_v(1:2*P-1,2,P,v_l)*(1.0d0-PI_f_v(1:2*P-1,2,P,v_l,u_l))*matmul(F(1:2*P-1,1:2*P-1,2,3),V_old(1:2*P-1,3))) & !Success in the new and no failure of the old
                        +beta*(PI_s_v(1:2*P-1,2,P,v_l)*PI_f_v(1:2*P-1,2,P,v_l,u_l)+(1.0d0-PI_s_v(1:2*P-1,2,P,v_l))*(1.0d0-PI_f_v(1:2*P-1,2,P,v_l,u_l)))*matmul(F(1:2*P-1,1:2*P-1,2,2),V_old(1:2*P-1,2)) & ! Success and failure of the old or failure of the new but no faile of the old
                        +beta*((1.0d0-PI_s_v(1:2*P-1,2,P,v_l))*PI_f_v(1:2*P-1,2,P,v_l,u_l)*matmul(F(1:2*P-1,1:2*P-1,2,1),V_old(1:2*P-1,1)))  !failure in both

                
                !Value function
                V_new(1:2*P-1,2)=rho(2)*log(1.0d0+exp((v_1I(1:2*P-1)-v_10(1:2*P-1))/rho(2)))+v_10(1:2*P-1)+rho(2)*gamma
                
                
            
            !!Two wells (n=3)
            !!!!!!!!!!!!!!!!!!
                do ind=1,2*P-1
                    if (CCP(ind,2)/=0.0d0 .and. CCP(ind,2)/=1.0d0) then
                        V_new(ind,3)=T_g+Ef_v(ind,3)-2.0d0*tau-vec_N(ind)*tau_per_N & 
                                        + CCP(ind,2)*(rho(2)*gamma-rho(2)*log(CCP(ind,2))) & 
                                        + (1.0d0-CCP(ind,2))*(rho(2)*gamma-rho(2)*log(1.0d0-CCP(ind,2))) &
                                        + beta*(1.0d0-PI_f_v(ind,3,P,v_l,u_l))**2.0d0*sum(F(ind,1:2*P-1,3,3)*V_old(1:2*P-1,3)) & !none fails
                                        + 2.0d0*beta*(1.0d0-PI_f_v(ind,3,P,v_l,u_l))*PI_f_v(ind,3,P,v_l,u_l)*sum(F(ind,1:2*P-1,3,2)*V_old(1:2*P-1,2)) & !one fails
                                        + beta*PI_f_v(ind,3,P,v_l,u_l)**2.0d0*sum(F(ind,1:2*P-1,3,1)*V_old(1:2*P-1,1))  !both fail
                    else
                        V_new(ind,3)=T_g+Ef_v(ind,3)-2.0d0*tau-vec_N(ind)*tau_per_N & 
                                        +rho(2)*gamma &                        
                                        + beta*(1.0d0-PI_f_v(ind,3,P,v_l,u_l))**2.0d0*sum(F(ind,1:2*P-1,3,3)*V_old(1:2*P-1,3)) & !none fails
                                        + 2.0d0*beta*(1.0d0-PI_f_v(ind,3,P,v_l,u_l))*PI_f_v(ind,3,P,v_l,u_l)*sum(F(ind,1:2*P-1,3,2)*V_old(1:2*P-1,2)) & !one fails
                                        + beta*PI_f_v(ind,3,P,v_l,u_l)**2.0d0*sum(F(ind,1:2*P-1,3,1)*V_old(1:2*P-1,1))  !both fail
                    end if
                end do
              
                                
               CCP(1:2*P-1,1)=1.0d0/(1.0d0+exp(v_00(1:2*P-1)/rho(2)-v_0I(1:2*P-1)/rho(2)))
               CCP(1:2*P-1,2)=1.0d0/(1.0d0+exp(v_10(1:2*P-1)/rho(2)-v_1I(1:2*P-1)/rho(2))) 
                
    
               
    
    if (isnan(sum(CCP)) .or. isnan(sum(V_new))) then
        print*,'error in vfi'
        print*,'v_00',v_00(1:2*P-1)
        print*,'v_0I',v_0I(1:2*P-1)
        print*,'v_10',v_10(1:2*P-1)
        print*,'v_1I',v_1I(1:2*P-1)
        print*,'CCP',CCP(1:2*P-1,:)
        read*,pause_k
    end if
    
    

end subroutine