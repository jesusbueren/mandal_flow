subroutine load_cadastral_maps()
    use cadastral_maps; use primitives;use simulation
    implicit none
    character(LEN=1)::s_c1
    character(LEN=2)::s_c2
    integer::v_l,i,j,ind,a_l,j2,ind2,s_l
    double precision::u,u_m,aux_w
    double precision,dimension(villages)::area_v
    integer,dimension(villages,plots_in_map)::area_type
    integer,dimension(P_max,types_a,villages,unobs_types)::counter_type
    integer,allocatable::seed(:)
    
    

    PA_type=-9
    
    OPEN(UNIT=12, FILE=file_map//"area_type.csv")
        read(12,*),areas
    close(12)
    
    OPEN(UNIT=12, FILE=path_primitives//"area_v.csv")
        read(12,*),area_v
    close(12)
    pr_v=area_v/sum(area_v)
    
    !Area Type
    total_area=0.0d0
    do v_l=1,villages
        do i=1,plots_v(v_l)
            do a_l=1,types_a-1
                if (areas(v_l,i)<=area_lims(a_l)) then
                    area_type(v_l,i)=a_l
                    exit
                else
                    area_type(v_l,i)=a_l+1
                end if
            end do
            total_area(v_l)=total_area(v_l)+area(area_type(v_l,i))
            PA_type(i,2,v_l,:)=area_type(v_l,i)
        end do 
    end do
       
    !Load map (who is connected to who)
    neighbors_map=0
    do v_l=1,villages
        if (v_l<10) then
            Write( s_c1, '(I1)' )  v_l
            OPEN(UNIT=12, FILE=file_map//"map_"//s_c1//".csv")
        else
            Write( s_c2, '(I2)' )  v_l
            OPEN(UNIT=12, FILE=file_map//"map_"//s_c2//".csv")
        end if
            read(12,*),neighbors_map(1:plots_v(v_l),1:plots_v(v_l),v_l)
        close(12)
        
        !Store which are my neighbors
        do i=1,plots_in_map
            neighbors_map(i,i,v_l)=1
            ind=1
            neighbors(i,ind,v_l,:)=i
            ind=ind+1
            do j=1,plots_in_map
                !the number of neighbors cannot be greater than P_max ortherwise I select the firt P_max neighbors 
                if (neighbors_map(i,j,v_l)==1 .and. j/=i .and. ind<=P_max) then 
                    neighbors(i,ind,v_l,:)=j 
                    ind=ind+1
                end if                  
            end do
        end do
        !number of neighbors
        do s_l=1,sims_tr
            PA_type(1:plots_v(v_l),1,v_l,s_l)=min(sum(neighbors_map(1:plots_v(v_l),1:plots_v(v_l),v_l),2),P_max)
        end do
        !if (v_l<10) then
        !    Write( s_c1, '(I1)' )  v_l
        !    OPEN(UNIT=12, FILE=file_map//"new_map_"//s_c1//".txt")
        !else
        !    Write( s_c2, '(I2)' )  v_l
        !    OPEN(UNIT=12, FILE=file_map//"new_map_"//s_c2//".txt")
        !end if
        !    write(12,*),neighbors_map(1:plots_v(v_l),1:plots_v(v_l),v_l)
        !close(12)
    end do



    !Generate permanent unobserved heterogeneity type
    
    !Call seed number
    call random_seed(size = v_l)
    allocate(seed(v_l))
    seed=321
    
    counter_type=0
    do s_l=1,sims_tr
        do v_l=1,villages;do i=1,plots_v(v_l)
            call RANDOM_NUMBER(u_m)
            if (u_m<pr_unobs_t(1)) then
                unobs_types_i(i,v_l,s_l)=1 !unobs_types_i(1,1,:)
            else
                unobs_types_i(i,v_l,s_l)=2
            end if
            counter_type(PA_type(i,1,v_l,s_l),PA_type(i,2,v_l,s_l),v_l,unobs_types_i(i,v_l,s_l))=counter_type(PA_type(i,1,v_l,s_l),&
                PA_type(i,2,v_l,s_l),v_l,unobs_types_i(i,v_l,s_l))+1
        end do;end do
    end do
 
    
    print*,''
    print*,dble(counter_type(3,1,1,:))/dble(sum(counter_type(3,1,1,:))) 
    

    end subroutine
    