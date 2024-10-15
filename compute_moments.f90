subroutine compute_moments()
    use simulation; use primitives
    implicit none

    double precision,dimension(types_a*2,types_a*2)::var_cov
    double precision,dimension(types_a*2)::diag
    double precision,dimension(15,bs_samples)::moments_all
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    double precision,dimension(3,villages,types_a)::pr_D_va
    integer::i,n_l,a_l,v_l
    double precision,dimension(types_a,types_a)::var_cov2
    
    OPEN(UNIT=12, FILE=path_primitives//"pr_drilling.csv")
        read(12,*),moment_own_nxa_data 
    close(12)
    
    OPEN(UNIT=12, FILE=path_primitives//"pr_drilling_var.csv")
        read(12,*),var_cov 
    close(12)
    
    OPEN(UNIT=12, FILE=path_primitives//"pr_v_na.csv")
        read(12,*),pr_v_na 
    close(12)
    
    OPEN(UNIT=12, FILE=path_primitives//"pr_Da.csv")
        read(12,*),Pr_D_a_data
    close(12)
    
    OPEN(UNIT=12, FILE=path_primitives//"pr_Da_var.csv")
        read(12,*),var_cov2
    close(12)
    
    do v_l=1,villages;do a_l=1,types_a
        pr_v_a(v_l,a_l)=sum(pr_v_na(v_l,a_l,:))/sum(pr_v_na(:,a_l,:)) 
    end do;end do
    
    do n_l=1,3;do a_l=1,types_a
        pr_v_na(:,a_l,n_l)=pr_v_na(:,a_l,n_l)/sum(pr_v_na(:,a_l,n_l))
    end do; end do
    

    do i=1,types_a*2
        diag(i)=var_cov(i,i)
    end do
    var_own_nxa=reshape(diag,(/types_a,2/))
    
    
    do i=1,types_a
        Pr_D_a_var(i)=var_cov2(i,i)
    end do
    

    do bs_l=1,bs_samples
        write (n2s, "(I3.3)") bs_l 
        OPEN(UNIT=12, FILE=path_results//"bootstrap_samples\"//"bs_mom"//trim(n2s)//".csv")
            read(12,*),moment_own_nxa_bs(:,:,bs_l)
        close(12)
        OPEN(UNIT=12, FILE=path_results//"bootstrap_samples\"//"bs_mom"//trim(n2s)//".csv")
            read(12,*),moments_all(:,bs_l)  !moments_all(:,1) 
        close(12)
        moments_D_bs(:,bs_l)=moments_all(11:15,bs_l)
    end do
    
    
    
    end subroutine
    
    
double precision function c4_normal_01 (  )
!--------------------------------------------------------------------------------------------------------------------------------
  implicit none
  double precision, parameter :: r4_pi=3.14159265358979323846264338327950288419716939937510
  double precision:: v1
  double precision:: v2
  double precision:: x_c
  double precision:: x_r 
  call random_number(v1)
  call random_number(v2)
  x_r =sqrt(-2.0d0*log(v1))*cos(2.0d0*r4_pi*v2)
  x_c =sqrt(-2.0d0*log(v1))*sin(2.0d0*r4_pi*v2)
  c4_normal_01=x_r
  return
end
    
    