module distances
implicit none
contains
subroutine dist_segs(ra,&
                     ra_end,&
                     rb,&
                     rb_end,&
                     pa,&
                     pb,&
                     Gab,&
                     Gab_norm,&
                     Sba,&
                     coll_course)
implicit none
real(8), dimension(3):: ra, ra_end, rb, rb_end, pa, pb, Gab
real(8)              :: Sab, Sba, la, lb, Gab_norm
logical              :: coll_course

pa=ra_end-ra
pb=rb_end-rb

la=sqrt(dot_product(pa,pa))
lb=sqrt(dot_product(pb,pb))

pa=pa/la
pb=pb/lb

Sab=(dot_product(ra-rb, pb)*dot_product(pa,pb)-dot_product(ra-rb,pa))&
    /(1d0-dot_product(pa, pb)**2d0)

Sba=(dot_product(rb-ra, pa)*dot_product(pb,pa)-dot_product(rb-ra,pb))&
    /(1d0-dot_product(pb, pa)**2d0)

Gab=ra+Sab*pa-rb-Sba*pb
Gab_norm=sqrt(dot_product(Gab,Gab))

coll_course=.false.
if (Sab.le.la .and. Sab .ge. 0.0) then
	if (Sba.le.lb .and. Sba .ge. 0.0) then
		coll_course= .true.
	end if
end if
end subroutine dist_segs

!************************************************************************

subroutine dist_pt_seg(ra,&
                       rb,&
                       rb_end,&
		       pb,&
                       Gab,&
                       Gab_norm,&
                       Sba,&
                       coll_course)

real(8), dimension(3):: ra, rb, rb_end, pb, Gab
real(8)              :: Sba, Gab_norm, lb
logical              :: coll_course


pb=rb_end-rb

lb=sqrt(dot_product(pb,pb))

pb=pb/lb

Sba=dot_product(ra-rb, pb)

Gab=ra-rb-Sba*pb
Gab_norm=sqrt(dot_product(Gab,Gab))


coll_course=.false.

if (Sba.le.lb .and. Sba .ge. 0.0) then
	coll_course= .true.
end if

end subroutine dist_pt_seg

!************************************************************************

subroutine dist_pt_pt(ra,&
                      rb,&
                      Gab,&
                      Gab_norm)
real(8), dimension(3):: ra, rb, Gab
real(8)              :: Gab_norm


Gab=ra-rb
Gab_norm=sqrt(dot_product (Gab, Gab))

end subroutine dist_pt_pt

!------------------------------------------------------
subroutine collide   (ra,&
                     ra_end,&
                     rb,&
                     rb_end,&
                     fiber_diameter,&
                     fibs_collide)
                     

real(8), dimension(3):: ra, ra_end, rb, rb_end, pa, pb, Gab
real(8)              :: Sab, Sba, la, lb, Gab_norm, fiber_diameter
logical              :: coll_course
logical              :: fibs_collide 

fibs_collide=.false.

!********************************************

call dist_segs(rb_end,&
               rb,&
               ra_end,&
               ra,&
               pa,&
               pb,&
               Gab,&
               Gab_norm,&
               Sba,&
               coll_course)

if (Gab_norm.lt.fiber_diameter .and.&
	coll_course.eqv..true.) then		
	fibs_collide=.true.
end if

!********************************************

call dist_segs(rb,&
               rb_end,&
               ra_end,&
               ra,&
               pa,&
               pb,&
               Gab,&
               Gab_norm,&
               Sba,&
               coll_course)

if (Gab_norm.lt.fiber_diameter .and.&
	coll_course.eqv..true.) then		
	fibs_collide=.true.
end if


call dist_segs(rb_end,&
               rb,&
               ra,&
               ra_end,&
               pa,&
               pb,&
               Gab,&
               Gab_norm,&
               Sba,&
               coll_course)

if (Gab_norm.lt.fiber_diameter .and.&
	coll_course.eqv..true.) then		
	fibs_collide=.true.
end if
!*********************************************
call dist_segs(rb,&
               rb_end,&
               ra,&
               ra_end,&
               pa,&
               pb,&
               Gab,&
               Gab_norm,&
               Sba,&
               coll_course)

if (Gab_norm.lt.fiber_diameter .and.&
	coll_course.eqv..true.) then		
	fibs_collide=.true.
end if
!*********************************************
!*********************************************

call dist_segs(ra,&
               ra_end,&
               rb,&
               rb_end,&
               pa,&
               pb,&
               Gab,&
               Gab_norm,&
               Sba,&
               coll_course)

if (Gab_norm.lt.fiber_diameter .and.&
	coll_course.eqv..true.) then		
	fibs_collide=.true.
end if



!********************************************

call dist_segs(ra_end,&
               ra,&
               rb_end,&
               rb,&
               pa,&
               pb,&
               Gab,&
               Gab_norm,&
               Sba,&
               coll_course)

if (Gab_norm.lt.fiber_diameter .and.&
	coll_course.eqv..true.) then		
	fibs_collide=.true.
end if

!********************************************

call dist_segs(ra_end,&
               ra,&
               rb_end,&
               rb,&
               pa,&
               pb,&
               Gab,&
               Gab_norm,&
               Sba,&
               coll_course)

if (Gab_norm.lt.fiber_diameter .and.&
	coll_course.eqv..true.) then		
	fibs_collide=.true.
end if


call dist_segs(ra_end,&
               ra,&
               rb,&
               rb_end,&
               pa,&
               pb,&
               Gab,&
               Gab_norm,&
               Sba,&
               coll_course)

if (Gab_norm.lt.fiber_diameter .and.&
	coll_course.eqv..true.) then		
	fibs_collide=.true.
end if
!*********************************************
call dist_segs(ra,&
               ra_end,&
               rb_end,&
               rb,&
               pa,&
               pb,&
               Gab,&
               Gab_norm,&
               Sba,&
               coll_course)

if (Gab_norm.lt.fiber_diameter .and.&
	coll_course.eqv..true.) then		
	fibs_collide=.true.
end if
!*********************************************

call dist_segs(ra,&
               ra_end,&
               rb,&
               rb_end,&
               pa,&
               pb,&
               Gab,&
               Gab_norm,&
               Sba,&
               coll_course)

if (Gab_norm.lt.fiber_diameter .and.&
	coll_course.eqv..true.) then		
	fibs_collide=.true.
end if
!*********************************************
!*********************************************
call dist_pt_seg(ra,&
                 rb,&
                 rb_end,&
		 pb,&
                 Gab,&
                 Gab_norm,&
                 Sba,&
                 coll_course)

if (Gab_norm.lt.fiber_diameter .and.&
	coll_course.eqv..true.) then		
	fibs_collide=.true.
end if

!------------------------------------------------------
call dist_pt_seg(ra_end,&
                 rb,&
                 rb_end,&
		 pb,&
                 Gab,&
                 Gab_norm,&
                 Sba,&
                 coll_course)

if (Gab_norm.lt.fiber_diameter .and.&
	coll_course.eqv..true.) then		
	fibs_collide=.true.
end if
!------------------------------------------------------
call dist_pt_seg(rb,&
                 ra,&
                 ra_end,&
		 pb,&
                 Gab,&
                 Gab_norm,&
                 Sba,&
                 coll_course)

if (Gab_norm.lt.fiber_diameter .and.&
	coll_course.eqv..true.) then		
	fibs_collide=.true.
end if
!------------------------------------------------------
call dist_pt_seg(rb_end,&
                 ra,&
                 ra_end,&
		 pb,&
                 Gab,&
                 Gab_norm,&
                 Sba,&
                 coll_course)

if (Gab_norm.lt.fiber_diameter .and.&
	coll_course.eqv..true.) then		
	fibs_collide=.true.
end if
!*********************************************
!*********************************************
call dist_pt_pt(ra,&
                rb,&
                Gab,&
                Gab_norm)
                        
if (Gab_norm.lt.fiber_diameter) then		
	fibs_collide=.true.
end if
!------------------------------------------------------   
call dist_pt_pt(ra,&
                rb_end,&
                Gab,&
                Gab_norm)
                        
if (Gab_norm.lt.fiber_diameter) then		
	fibs_collide=.true.
end if
!------------------------------------------------------   
call dist_pt_pt(ra_end,&
                rb,&
                Gab,&
                Gab_norm)
                        
if (Gab_norm.lt.fiber_diameter) then		
	fibs_collide=.true.
end if
!------------------------------------------------------
call dist_pt_pt(ra_end,&
                rb_end,&
                Gab,&
                Gab_norm)
                        
if (Gab_norm.lt.fiber_diameter) then		
	fibs_collide=.true.
end if
!------------------------------------------------------     
 
end subroutine collide

!*********************************************
!*********************************************
SUBROUTINE init_random_seed()
INTEGER :: i, n, clock
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)
DEALLOCATE(seed)
END SUBROUTINE init_random_seed
!*********************************************
 
end module distances

