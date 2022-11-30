!subroutine simulate_spatial_correlation()
!    use cadastral_maps; use primitives;use simulation
!    implicit none
!    character(LEN=1)::s_c1
!    character(LEN=2)::s_c2
!    integer::v_l,i,j,ind,a_l,t_l,count_unique_mode,count_mode_equal_type,s_l
!    double precision::u
!    double precision,dimension(villages,plots_in_map)::areas
!    integer,dimension(villages,plots_in_map)::area_type
!    integer,dimension(P_max,2)::PA_stat
!    integer,dimension(1)::seed=321
!    integer,dimension(plots_in_map,villages,sims)::unobs_types_i_new
!    double precision,dimension(unobs_types)::pr_type
!    integer,dimension(unobs_types)::counter_type,counter_all
!
!    do s_l=1,sims; do v_l=1,villages
!        if (s_l==1) then
!            !if (v_l<10) then
!            !    Write( s_c1, '(I1)' )  v_l
!            !    OPEN(UNIT=12, FILE=file_map//"map_types_iid_sc_"//s_c1//".txt")
!            !else
!            !    Write( s_c2, '(I2)' )  v_l
!            !    OPEN(UNIT=12, FILE=file_map//"map_types_iid_sc_"//s_c2//".txt")
!            !end if        
!            !do i=1,plots_v(v_l)
!            !    write(12,'(I3)'),unobs_types_i(i,v_l,s_l)
!            !end do
!            !close(12)
!                
!        end if
!    
!        do t_l=1,500
!        counter_all=0
!        count_unique_mode=0
!        count_mode_equal_type=0
!        do i=1,plots_v(v_l)
!            counter_type=0.0d0
!            do j=1,PA_type(i,1,v_l,s_l) 
!                if (neighbors(i,j,v_l,s_l)/=i) then !neighbors(i,:,v_l)
!                    counter_type(unobs_types_i(neighbors(i,j,v_l,s_l),v_l,s_l))=counter_type(unobs_types_i(neighbors(i,j,v_l,s_l),v_l,s_l))+1
!                end if
!            end do
!            if (counter_type(1)==counter_type(2) .and. counter_type(3)==counter_type(2)) then
!                pr_type=1.0d0/3.0d0
!            elseif (counter_type(1)==counter_type(2) .and. counter_type(1)>counter_type(3)) then
!                pr_type(1:2)=(1.0d0+rho_sc)/4.0d0
!                pr_type(3)=(1.0d0-rho_sc)/2.0d0
!            elseif (counter_type(1)==counter_type(3) .and. counter_type(1)>counter_type(2)) then
!                pr_type(1)=(1.0d0+rho_sc)/4.0d0
!                pr_type(3)=(1.0d0+rho_sc)/4.0d0
!                pr_type(2)=(1.0d0-rho_sc)/2.0d0
!            elseif (counter_type(2)==counter_type(3) .and. counter_type(2)>counter_type(1)) then
!                pr_type(2:3)=(1.0d0+rho_sc)/4.0d0
!                pr_type(1)=(1.0d0-rho_sc)/2.0d0
!            else
!                pr_type=(1.0d0-rho_sc)/2.0d0
!                pr_type(maxloc(counter_type))=rho_sc
!                count_unique_mode=count_unique_mode+1
!                if (unobs_types_i(i,v_l,s_l)==maxloc(counter_type,1)) then
!                    count_mode_equal_type=count_mode_equal_type+1
!                end if
!            end if
!            call RANDOM_NUMBER(u)
!            if (u<pr_type(1)) then
!                unobs_types_i_new(i,v_l,s_l)=1
!            elseif (u<sum(pr_type(1:2))) then
!                unobs_types_i_new(i,v_l,s_l)=2
!            else
!                unobs_types_i_new(i,v_l,s_l)=3
!            end if
!            counter_all(unobs_types_i_new(i,v_l,s_l))=counter_all(unobs_types_i_new(i,v_l,s_l))+1
!            unobs_types_i(i,v_l,s_l)=unobs_types_i_new(i,v_l,s_l)
!        end do
!        !print '(I5,<4>F20.3)',t_l,dble(counter_all)/dble(plots_v(v_l)),dble(count_mode_equal_type)/dble(count_unique_mode)
!        unobs_types_i(:,v_l,s_l)=unobs_types_i_new(:,v_l,s_l)
!        end do
!        
!        if (s_l==1) then
!            !if (v_l<10) then
!            !    Write( s_c1, '(I1)' )  v_l
!            !    OPEN(UNIT=12, FILE=file_map//"map_types_99_sc_"//s_c1//".txt")
!            !else
!            !    Write( s_c2, '(I2)' )  v_l
!            !    OPEN(UNIT=12, FILE=file_map//"map_types_99_sc_"//s_c2//".txt")
!            !end if        
!            !do i=1,plots_v(v_l)
!            !    write(12,'(I3)'),unobs_types_i(i,v_l,s_l)
!            !end do
!            !close(12)
!        end if
!        
!    end do;end do
!    
!    
!    
!    
!    
!    
!    
!end subroutine