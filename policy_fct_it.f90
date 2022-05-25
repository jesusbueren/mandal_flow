subroutine policy_fct_it(Ef_v,F,P,CCP_in,CCP_out,v_l,u_l,V_new,a_l)
    use dimensions;use primitives
    implicit none
    integer,intent(in)::P
    double precision,dimension(2*P-1,3),intent(in)::Ef_v
    double precision,dimension(2*P-1,2*P-1,3,3),intent(in)::F
    double precision,dimension(2*P-1,2),intent(in)::CCP_in
    integer,intent(in)::v_l,u_l,a_l
    double precision,dimension(2*P-1,2),intent(out)::CCP_out
    double precision,dimension((2*P-1),3),intent(out)::V_new
    double precision,dimension(2*P-1)::CCP_aux
    double precision,dimension(2*P-1,2)::CCP_old
    double precision,dimension((2*P-1)*3,(2*P-1)*3)::C
    double precision,dimension((2*P-1)*3,1)::V
    double precision::dist,crit
    integer::it
    character::pause_k
    
    it=0
    CCP_old=CCP_in
    crit=1.0d-9
    
    do it=1,2*P-1
        if (sum(F(it,:,1,1))<0.98d0) then
            print*,'error in pol it',F(it,:,1,1)
            print*,'sum',sum(F(it,:,1,1))
            read*,pause_k
        end if
    end do
    
    CCP_aux=1.0d0/(1.0d0+exp(-(-PI_s_v(1:2*P-1,2,P,v_l)*c_s-(1.0d0-PI_s_v(1:2*P-1,2,P,v_l))*c_d)/rho(2)))
    
        
    call generate_C(CCP_old,F,P,C,v_l,u_l) 
    call valuation(CCP_old,C,Ef_v,P,V,v_l,u_l,a_l,CCP_aux)
    call update_CCP(V,F,Ef_v,P,CCP_out,v_l,u_l,a_l)
    
    dist=maxval(abs((reshape(CCP_old,(/(2*P-1)*2/))-reshape(CCP_out,(/(2*P-1)*2/)))))
    it=it+1
    
    if (isnan(sum(CCP_out)))then
        print*,'problem in policy fct it: ',P,v_l,u_l
        dist=1000
        CCP_out=CCP_old
        read*,pause_k
    end if
    
    if (it>1000) then
        print*,'P',P
        print*,CCP_old(1:2*P-1,1)
        print*,CCP_out(1:2*P-1,1)
        print*,'policy fct it not converged'
        return
    end if

    !if (dist>crit)then
    !    CCP_old=CCP_out
    !    go to 1
    !end if
    
    V_new=reshape(V,(/(2*P-1),3/))

end subroutine