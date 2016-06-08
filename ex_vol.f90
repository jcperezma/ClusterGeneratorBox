module distances
implicit none

type bound_box
real(8), dimension(3):: X
real(8), dimension(3):: W
end type bound_box

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

                      
                      
function clamp( a, b, c)
real(8)              ::a, b, c, clamp

if ( a < b ) then
clamp = b
return
end if

if ( a > c ) then
clamp = c
return
end if

clamp = a
end function clamp
 
subroutine dist_segments(p1,&
                     q1,&
                     p2,&
                     q2,&
                     s,&
                     t,&
                     Gab,&
                     Gab_norm)


real(8), dimension(3):: p1, q1, p2, q2, d1, d2, Gab, r, c1, c2
real(8)              ::Gab_norm, epsilon, s, t
real(8)              ::a, e, f, c, b, denom
logical              :: coll_course

epsilon  =1D-40

d1=q1-p1
d2=q2-p2

r = p1 - p2

a = dot_product(d1,d1)
e = dot_product(d2,d2)
f = dot_product(d2,r)

if ((a<=epsilon) .and. (e<=epsilon)  ) then
    s = 0.0D0
    t = 0.0D0
    
    c1 = p1
    c2 = p2
    
    Gab_norm = sqrt(dot_product(c1 - c2, c1 - c2))
    Gab =  (c1 - c2)
    return
end if
    
if ( a<= epsilon ) then
    s = 0D0
    t = f/e
    t = clamp(t, 0.0D0, 1.0D0)    
else
    c = dot_product(d1,r)
    
    if(e <= epsilon)then
        t= 0.0D0
        s = clamp(-c /a, 0.0D0, 1.0D0)
    else
        b = dot_product(d1, d2)
        denom = a*e-b*b
        
        if(denom /=	0.0D0 ) then
            s = clamp((b*f -c*e) / denom  , 0.0D0, 1.0D0) 
        else
         s =0.0D0
        end if
        
        t =(b*s +f) / e
        
        if (t < 0D0) then
            t = 0.0D0
            s =clamp(-c / a, 0.0D0 , 1.0D0)
        else if ( t> 1.0) then
            
            t= 1.0D0
            s =clamp((b-c)/a, 0.0D0, 1.0D0)
            
       end if
        
        
   end if


    
end if

    c1 = p1 + d1*s
    c2 = p2 + d2*t
    !print *, "c1 " , c1
    !print *, "c2 " , c2
     !print *, "s " , s
     !print *, "t " , t 
     !print *, "p1 " , p1 
     !print *, "p2 " , p2 
     !print *, "d1 " , d1 
     !print *, "d2 " , d2 
     !print *, " "
    Gab_norm = sqrt(dot_product(c1 - c2, c1 - c2))
    Gab =  (c1 - c2)


end subroutine dist_segments

!------------------------------------------------------
function are_boxes_intersect(hinge_pt1, hinge_pt2,segment_pt1, segment_pt2, fiber_diameter)

logical               ::are_boxes_intersect
real(8), dimension(3) ::hinge_pt1, hinge_pt2,segment_pt1, segment_pt2 
type(bound_box)       ::a, b
real(8)               :: fiber_diameter
are_boxes_intersect=.false.
a.X=(hinge_pt1+hinge_pt2)/2
a.W=abs(hinge_pt1-hinge_pt2)+fiber_diameter


b.X=(segment_pt1+segment_pt2)/2
b.W=abs(segment_pt1-segment_pt2)+fiber_diameter


are_boxes_intersect= (abs(a.X(1) - b.X(1)) * 2 < (a.W(1) + b.W(1))) .and.&
                     (abs(a.X(2) - b.X(2)) * 2 < (a.W(2) + b.W(2))) .and.&
                     (abs(a.X(3) - b.X(3)) * 2 < (a.W(3) + b.W(3)))  

end function are_boxes_intersect

                     
subroutine collide2   (ra,&
                     ra_end,&
                     rb,&
                     rb_end,&
                     fiber_diameter,&
                     fibs_collide)
                     

real(8), dimension(3):: ra, ra_end, rb, rb_end, Gab
real(8)              :: s, t, Gab_min, fiber_diameter
logical              :: fibs_collide 


call dist_segments(ra,&
               ra_end,&
               rb,&
               rb_end,&
                     s,&
                     t,&
                     Gab,&
                     Gab_min)

if (Gab_min .LT. fiber_diameter) fibs_collide = .true.


end subroutine collide2

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

