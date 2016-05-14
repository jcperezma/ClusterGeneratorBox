!LOCATES FIBERS ACCORDING TO A FILE WITH ORIENTATION

!*********************************************
!*********************************************

program fiber_position
use distances
implicit none
real(8), parameter                   :: pi=3.14159
real(8), dimension(3)                :: box_size
real(8), dimension(3)                :: box_location
real(8), dimension(3)                :: ini_pos
real(8), dimension(3)                :: vec_unit
real(8)                              :: fiber_length            
real(8)                              :: direction
real(8)                              :: max_theta
real(8)                              :: phi
real(8)                              :: theta
real(8)                              :: total_length
real(8)                              :: total_length_simulation
real(8)                              :: volumetric_fraction
real(8)                              :: fiber_diameter
real(8)                              :: vol_frac
real(8)                              :: r
real(8)                              :: segment_length
real(8),dimension(:,:), allocatable  :: coords_ini
real(8),dimension(:,:), allocatable  :: coords_end
real(8),dimension(3)                 :: coords
integer                              :: i, j, l, k, nbr_bins, kk, jj, ll
integer                              :: nbr_fibers
integer                              :: nbr_segments
integer                              :: orientation
!*********************************************
real(8), dimension(3)                :: ra, ra_end, rb, rb_end, pa, pb, Gab
real(8)                              :: Sab, Sba, la, lb, Gab_norm
logical                              :: flag, fibs_collide
real(8), dimension(:),  allocatable  :: bin_centers
real(8), dimension(:),  allocatable  :: bin_heights
real(8), dimension(:,:),allocatable  :: bin_limits
real(8)                              :: del_bin
!*********************************************

namelist /input/ fiber_diameter,&         !Average fiber diameter
                 fiber_length,&           !Average fiber length
                 nbr_segments,&           !number of segments per fiber
                 volumetric_fraction,&    !Volumetric Fraction
                 max_theta,&              !Maximum fiber inclination
                 box_size,&               !Cluster Size
                 box_location             !Cluster location

                         
open(3,file='../INPUT/input.in', status='old')

open(30,file='../INPUT/initial_histogram.in', status='old')

read(30,*),nbr_bins
print *,"hola", nbr_bins
allocate (bin_limits(nbr_bins,2))
allocate (bin_heights(nbr_bins))
print *,"hola22", nbr_bins
do i=1,nbr_bins
    read(30,*), bin_heights(i)
end do
print *,"hola2222", nbr_bins
del_bin=pi/nbr_bins

do i=1,nbr_bins
    bin_limits(i,1)=(i-1)*del_bin
    bin_limits(i,2)=(i)  *del_bin
end do
print *,"hola33333", nbr_bins
read(3,nml = input)
print *,"hola44444", nbr_bins
 close(3)
print *, "Box Location", box_location
print *, "Fiber Lenght", fiber_length
print *, "Fiber Diameter", fiber_diameter

!stop
!*********************************************
!*********************************************
open(1, file="../OUTPUT/coords.txt")
!open(2, file="../OUTPUT/coords_matlab.txt")

  
!print *, "FIB DIAM", fiber_diameter

segment_length=fiber_length/nbr_segments
print *,"Segment Length", segment_length

if (segment_length<=(2.5d0*fiber_diameter)) then
	print *, "Segment Length is too small when compared to the radius"
	stop
end if

if (fiber_length>(0.95*box_size(1))) then
	print *, "Fiber Length is close to or exceeds Box Size"
	stop
end if

if (fiber_length>(0.95*box_size(2))) then
	print *, "Fiber Length is close to or exceeds Box Size"
	stop
end if

!if (fiber_length>(0.95*box_size(3))) then
!	print *, "Fiber Length is close to or exceeds Box Size"
!	stop
!end if

total_length_simulation=((box_size(1)*box_size(2)*box_size(3))*volumetric_fraction)/(pi*(fiber_diameter/2d0)**2d0)

total_length=0

nbr_fibers=floor(total_length_simulation/fiber_length)
bin_heights=nint(nbr_fibers*bin_heights/sum(bin_heights))
nbr_fibers=sum(bin_heights)

do i=1, ubound(bin_heights,1)
    print *, bin_heights(i)
end do


print *, "nbr_fibers", nbr_fibers

write(1,*),nbr_fibers
!write(2,*),nbr_fibers
allocate(coords_ini(nbr_fibers,3))
allocate(coords_end(nbr_fibers,3))


vol_frac=0
total_length=0

kk=1

do i=1,nbr_bins
    print *,"BIN",i
    !print *,"KK",kk
    do ll=1,bin_heights(i)
        flag=.false.
        do while (flag .eqv. .false.)
	        flag=.true.
	        do j=1,3
		        call random_number(r)
			    ini_pos(j)=box_location(j)+(box_size(j)/2)*(2*(r-0.5))
	        end do
	        !print *, flag 
	        direction=1
	
	        call random_number(r)
	        phi=pi/2+r*(bin_limits(i,2)-bin_limits(i,1))+bin_limits(i,1)

            call random_number(r)
	        theta=r*max_theta
			        
            vec_unit(1)=cos(theta)*cos(phi)
            vec_unit(2)=cos(theta)*sin(phi)
            vec_unit(3)=sin(theta)
            !print *, "vec_Length", sqrt(dot_product(vec_unit, vec_unit))
	        do j=1,3
		        coords_ini(kk,j)=ini_pos(j)
		        coords_end(kk,j)=ini_pos(j)+fiber_length*direction*vec_unit(j)
            end do				
	        do j=1, 3
 		        if 	((coords_ini(kk,j) .gt. box_location(j)+box_size(j)/2).or.&
	    	   	      coords_ini(kk,j) .lt. box_location(j)-box_size(j)/2)then
			          flag=.false.
		        end if
	        end do
            !print *,'ddd',flag
	        do j=1, 3
 		        if 	((coords_end(kk,j) .gt. box_location(j)+box_size(j)/2).or.&
	    	   	      coords_end(kk,j) .lt. box_location(j)-box_size(j)/2)then
			        flag=.false.
		        end if
	        end do

	        do j=1,3
		        rb(j)=coords_ini(kk,j)
		        rb_end(j)=coords_end(kk,j)		 
	        end do
            !print *, "t",flag
            if (kk.ne.1) then
		        do k=1, kk-1
			        do j=1,3
				        ra(j)    =coords_ini(k,j)
				        ra_end(j)=coords_end(k,j)	 
			        end do
			        fibs_collide=.false.
			        call collide   (ra,&
                 				    ra_end,&
                 				    rb,&
                 				    rb_end,&
                 				    fiber_diameter,&
                 				    fibs_collide)
	
			        if (fibs_collide.eqv. .true.) then
				        flag=.false.
			        end if                        
           !print *,"ch4ttttttttttttttttttttttttttttttttttttttttt" 
		        end do
	        end if
        end do
	
        do j=1,3
	        coords(j)=coords_end(kk,j)-coords_ini(kk,j)
        end do
        total_length=total_length+fiber_length
        !print *, "Number of Segments", nbr_segments
        write (1,*), nbr_segments+1
        write (2,*), nbr_segments+1
        do l=1,nbr_segments+1
	        do j=1,3
		        coords(j)=ini_pos(j)+(l-1)*segment_length*direction*vec_unit(j)
	        end do
	        write (1,*),0, real(coords(1),4), real(coords(2),4), real(coords(3),4)
	        !write (2,*),   coords(1), coords(2), coords(3)
        end do
    	
        vol_frac=(total_length*pi*(fiber_diameter/2)**2)/(box_size(1)*box_size(2)*box_size(3))
        print *,"Fiber", kk,"of",nbr_fibers,"   ","Volumetric Fraction",vol_frac
        kk=kk+1
        
        if ((kk).gt. nbr_fibers)  stop
    end do
end do

end program fiber_position



