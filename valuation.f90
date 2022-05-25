subroutine valuation(CCP,C,Ef_v,P,V,v_l,u_l,a_l,CCP_aux)
    use dimensions; use primitives
    integer,intent(in)::P
    double precision,dimension((2*P-1)*3,(2*P-1)*3),intent(in)::C
    double precision,dimension(2*P-1,2),intent(in)::CCP
    double precision,dimension(2*P-1),intent(in)::CCP_aux
    double precision,dimension(2*P-1,3),intent(in)::Ef_v
    integer,intent(in)::v_l,u_l,a_l
    double precision,dimension((2*P-1)*3,1),intent(out)::V
    double precision,dimension((2*P-1)*3,1)::U
    double precision,dimension(2*P-1,3)::U_small
    integer::n_l,m_l,ind
    

    
    U_small=0.0d0
    do ind=1,2*P-1; do n_l=1,3
        if (n_l==1) then
            if (CCP(ind,n_l)==1.0d0)then
                U_small(ind,n_l)=-c_s*PI_s_v(ind,n_l,P,v_l)-c_d*(1.0d0-PI_s_v(ind,n_l,P,v_l))              
            elseif (CCP(ind,n_l)==0.0d0)then
                U_small(ind,n_l)=v_nod
            else
                U_small(ind,n_l)=CCP(ind,n_l)*(-c_s*PI_s_v(ind,n_l,P,v_l)-c_d*(1.0d0-PI_s_v(ind,n_l,P,v_l)))
            end if
        elseif (n_l==2) then
            if (CCP(ind,n_l)==1.0d0)then
                U_small(ind,n_l)=-c_s*PI_s_v(ind,n_l,P,v_l)-c_d*(1.0d0-PI_s_v(ind,n_l,P,v_l))              
            elseif (CCP(ind,n_l)==0.0d0)then
                U_small(ind,n_l)=0.0d0
            else
                U_small(ind,n_l)=CCP(ind,n_l)*(-c_s*PI_s_v(ind,n_l,P,v_l)-c_d*(1.0d0-PI_s_v(ind,n_l,P,v_l)))
            end if
        else     
            U_small(ind,n_l)=0.0d0
        end if
        !if (n_l<3) then
        !    U_small(ind,n_l)=U_small(ind,n_l)+CCP(ind,n_l)*(rho(1)*gamma-rho(1)*log(CCP(ind,n_l)))+&
        !                            (1.0d0-CCP(ind,n_l))*(rho(1)*gamma-rho(1)*log(1.0d0-CCP(ind,n_l)))-rho(1)*gamma
        !else
        !    U_small(ind,n_l)=U_small(ind,n_l)+CCP(ind,2)*(rho(1)*gamma-rho(1)*log(CCP(ind,2)))+&
        !                            (1.0d0-CCP(ind,2))*(rho(1)*gamma-rho(1)*log(1.0d0-CCP(ind,2)))-rho(1)*gamma
        !end if
    end do; end do
    
    U_small=U_small+Ef_v
    !Tax on having a well
    if (social==1) then
        U_small(:,2)=U_small(:,2)-c_e
        U_small(:,3)=U_small(:,3)-2.0d0*c_e
    end if
    
    
    do n_l=1,3
        ind=(2*P-1)*(n_l-1)+1
        U(ind:ind+2*P-1-1,1)=U_small(:,n_l)
    end do
    V=matmul(C,U)
   
    
end subroutine