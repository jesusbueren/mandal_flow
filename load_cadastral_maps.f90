subroutine load_cadastral_maps()
    use cadastral_maps; use primitives;use simulation
    implicit none
    character(LEN=1)::s_c1
    character(LEN=2)::s_c2
    integer::v_l,i,j,ind,a_l,j2,ind2,s_l
    double precision::u,u_m
    double precision,dimension(villages,plots_in_map)::areas
    integer,dimension(villages,plots_in_map)::area_type
    integer,dimension(P_max,types_a,villages,unobs_types)::counter_type
    integer,dimension(P_max,2)::PA_stat
    integer,dimension(1)::seed=321
    
    PA_type=-9
    
    OPEN(UNIT=12, FILE=file_map//"area_type.csv")
        read(12,*),areas
    close(12)
    
    !Area Type
    mean_area=0.0d0
    do v_l=1,villages;do s_l=1,sims
        do i=1,plots_v(v_l)
            do a_l=1,types_a-1
                if (areas(v_l,i)<=area_lims(a_l)) then
                    area_type(v_l,i)=a_l
                    exit
                else
                    area_type(v_l,i)=a_l+1
                end if
            end do
            mean_area(v_l)=dble(i-1)/dble(i)*mean_area(v_l)+1.0d0/dble(i)*area(area_type(v_l,i))
            PA_type(i,2,v_l,s_l)=area_type(v_l,i)
        end do 
    end do; end do 
       
    !Load map (who is connected to who)
    neighbors_map=0
    do v_l=1,villages;do s_l=1,sims
        if (v_l<10) then
            Write( s_c1, '(I1)' )  v_l
            OPEN(UNIT=12, FILE=file_map//"map_"//s_c1//".csv")
        else
            Write( s_c2, '(I2)' )  v_l
            OPEN(UNIT=12, FILE=file_map//"map_"//s_c2//".csv")
        end if
            read(12,*),neighbors_map(1:plots_v(v_l),1:plots_v(v_l),v_l)
        close(12)
        
        !Set zombie plots
        do i=1,plots_v(v_l)
            call RANDOM_NUMBER(u)
            if (u<pr_non_zombie(area_type(v_l,i))) then
                active_plots(i,v_l,s_l)=1
            else
                active_plots(i,v_l,s_l)=0 
            end if
        end do
        

        !Store which are my neighbors
        do i=1,plots_in_map
            neighbors_map(i,i,v_l)=1
            ind=1
            neighbors(i,ind,v_l,s_l)=i
            ind=ind+1
            do j=1,plots_in_map
                if (neighbors_map(i,j,v_l)==1 .and. j/=i .and. ind<=P_max) then !the number of neighbors cannot be greater than P_max ortherwise I select the firt P_max neighbors 
                    neighbors(i,ind,v_l,s_l)=j 
                    ind=ind+1
                end if                  
            end do
        end do
        !number of neighbors
        PA_type(1:plots_v(v_l),1,v_l,s_l)=min(sum(neighbors_map(1:plots_v(v_l),1:plots_v(v_l),v_l),2),P_max)
        !if (v_l<10) then
        !    Write( s_c1, '(I1)' )  v_l
        !    OPEN(UNIT=12, FILE=file_map//"new_map_"//s_c1//".txt")
        !else
        !    Write( s_c2, '(I2)' )  v_l
        !    OPEN(UNIT=12, FILE=file_map//"new_map_"//s_c2//".txt")
        !end if
        !    write(12,*),neighbors_map(1:plots_v(v_l),1:plots_v(v_l),v_l)
        !close(12)
    end do;end do
    !PA_type(:,1,1)
    
    !neighbors(398,:,1,1)=0 PA_type(957,1,1,1)
    
    v_l=1
    PA_stat=0
    do i=1,plots_v(v_l);do s_l=1,sims
        PA_stat(PA_type(i,1,v_l,s_l),1)=PA_stat(PA_type(i,1,v_l,s_l),1)+1
        PA_stat(PA_type(i,2,v_l,s_l),2)=PA_stat(PA_type(i,2,v_l,s_l),2)+1
    end do;end do
    !print*,'distribution of number of neighbors'
    !print*,dble(PA_stat(1:P_max,1))/dble(plots_v(v_l))
    !print*,'distribution of area types'
    !print*,dble(PA_stat(1:types_a,2))/dble(plots_v(v_l))

    !Generate permanent unobserved heterogeneity type
    counter_type=0
    do s_l=1,sims;do v_l=1,villages;do i=1,plots_v(v_l)
        call RANDOM_NUMBER(u_m)
        if (u_m<pr_unobs_t(1)) then
            unobs_types_i(i,v_l,s_l)=1
        elseif (u_m<sum(pr_unobs_t(1:2))) then
            unobs_types_i(i,v_l,s_l)=2
        else
            unobs_types_i(i,v_l,s_l)=3
        end if
        counter_type(PA_type(i,1,v_l,s_l),PA_type(i,2,v_l,s_l),v_l,unobs_types_i(i,v_l,s_l))=counter_type(PA_type(i,1,v_l,s_l),PA_type(i,2,v_l,s_l),v_l,unobs_types_i(i,v_l,s_l))+1
    end do;end do;end do
 
    !print*,dble(counter_type(6,:,11,:))/dble(sum(counter_type(3,1,1,:)))

    !call simulate_spatial_correlation()

    
    counter_type=0
    do s_l=1,sims;do v_l=1,villages;do i=1,plots_v(v_l)
        counter_type(PA_type(i,1,v_l,s_l),PA_type(i,2,v_l,s_l),v_l,unobs_types_i(i,v_l,s_l))=counter_type(PA_type(i,1,v_l,s_l),PA_type(i,2,v_l,s_l),v_l,unobs_types_i(i,v_l,s_l))+1
    end do;end do;end do
    
    print*,''
    print*,dble(counter_type(3,1,1,:))/dble(sum(counter_type(3,1,1,:))) 
    

    end subroutine
    