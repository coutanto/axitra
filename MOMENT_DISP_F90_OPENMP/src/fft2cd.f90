	subroutine  fft2cd (a,m,iwk)                                       
!                                  specifications for arguments         
	implicit real*8 (a-h,o-z)
	integer            m,iwk(*)                                       
	complex*16      a(*)                                           
!                                  specifications for local variables   
	integer     i,isp,j,jj,jsp,k,k0,k1,k2,k3,kb, &
                    kn,mk,mm,mp,n,n4,n8,n2,lm,nn,jk 
	real*8      rad,c1,c2,c3,s1,s2,s3,ck,sk,sq,a0,a1,a2,a3,     &
                    b0,b1,b2,b3,twopi,temp,                        &
                    zero,one,z0(2),z1(2),z2(2),z3(2)               
	complex*16  za0,za1,za2,za3,ak2                            
	equivalence  (za0,z0(1)),(za1,z1(1)),(za2,z2(1)),           &
                    (za3,z3(1)),(a0,z0(1)),(b0,z0(2)),(a1,z1(1)), & 
                    (b1,z1(2)),(a2,z2(1)),(b2,z2(2)),(a3,z3(1)),   &
                    (b3,z3(2))                                     
	data        sq,sk,ck,twopi/.7071068,.3826834,   &
                    .9238795,6.283185/                                
	data        zero/0.0/,one/1.0/                             
!                   sq=sqrt2/2,sk=sin(pi/8),ck=cos(pi/8) 
!                   twopi=2*pi                           
!                                  first executable statement           
	mp = m+1                                                          
	n = 2**m                                                          
	iwk(1) = 1                                                        
	mm = (m/2)*2                                                      
	kn = n+1                                                          
!                                  initialize work vector               
	do 5  i=2,mp                                                      
	   iwk(i) = iwk(i-1)+iwk(i-1)                                     
5	continue                                                          
	rad = twopi/n                                                     
	mk = m - 4                                                        
	kb = 1                                                            
	if (mm .eq. m) go to 15                                           
	k2 = kn                                                           
	k0 = iwk(mm+1) + kb                                               
10	k2 = k2 - 1                                                       
	k0 = k0 - 1                                                       
	ak2 = a(k2)                                                       
	a(k2) = a(k0) - ak2                                               
	a(k0) = a(k0) + ak2                                               
	if (k0 .gt. kb) go to 10                                          
15	c1 = one                                                          
	s1 = zero                                                         
	jj = 0                                                            
	k = mm - 1                                                        
	j = 4                                                             
	if (k .ge. 1) go to 30                                            
	go to 70                                                          
20	if (iwk(j) .gt. jj) go to 25                                      
	jj = jj - iwk(j)                                                  
	j = j-1                                                           
	if (iwk(j) .gt. jj) go to 25                                      
	jj = jj - iwk(j)                                                  
	j = j - 1                                                         
	k = k + 2                                                         
	go to 20                                                          
25	jj = iwk(j) + jj                                                  
	j = 4                                                             
30	isp = iwk(k)                                                      
	if (jj .eq. 0) go to 40                                           
!                                  reset trigonometric parameter(s       )
	c2 = jj * isp * rad                                               
	c1 = cos(c2)                                                      
	s1 = sin(c2)                                                      
35	c2 = c1 * c1 - s1 * s1                                            
	s2 = c1 * (s1 + s1)                                               
	c3 = c2 * c1 - s2 * s1                                            
	s3 = c2 * s1 + s2 * c1                                            
40	jsp = isp + kb                                                    
!                                  determine fourier coefficients       
!                                    in groups of 4                     
	do 50 i=1,isp                                                     
	   k0 = jsp - i                                                   
	   k1 = k0 + isp                                                  
	   k2 = k1 + isp                                                  
	   k3 = k2 + isp                                                  
	   za0 = a(k0)                                                    
	   za1 = a(k1)                                                    
	   za2 = a(k2)                                                    
	   za3 = a(k3)                                                    
	   if (s1 .eq. zero) go to 45                                     
	   temp = a1                                                      
	   a1 = a1 * c1 - b1 * s1                                         
	   b1 = temp * s1 + b1 * c1                                       
	   temp = a2                                                      
	   a2 = a2 * c2 - b2 * s2                                         
	   b2 = temp * s2 + b2 * c2                                       
	   temp = a3                                                      
	   a3 = a3 * c3 - b3 * s3                                         
	   b3 = temp * s3 + b3 * c3                                       
45	   temp = a0 + a2                                                 
	   a2 = a0 - a2                                                   
	   a0 = temp                                                      
	   temp = a1 + a3                                                 
	   a3 = a1 - a3                                                   
	   a1 = temp                                                      
	   temp = b0 + b2                                                 
	   b2 = b0 - b2                                                   
	   b0 = temp                                                      
	   temp = b1 + b3                                                 
	   b3 = b1 - b3                                                   
	   b1 = temp                                                      
	   a(k0) = cmplx(a0+a1,b0+b1)                                     
	   a(k1) = cmplx(a0-a1,b0-b1)                                     
	   a(k2) = cmplx(a2-b3,b2+a3)                                     
	   a(k3) = cmplx(a2+b3,b2-a3)                                     
50	continue                                                          
	if (k .le. 1) go to 55                                            
	k = k - 2                                                         
	go to 30                                                          
55	kb = k3 + isp                                                     
!                                  check for completion of final        
!                                    iteration                          
	if (kn .le. kb) go to 70                                          
	if (j .ne. 1) go to 60                                            
	k = 3                                                             
	j = mk                                                            
	go to 20                                                          
60	j = j - 1                                                         
	c2 = c1                                                           
	if (j .ne. 2) go to 65                                            
	c1 = c1 * ck + s1 * sk                                            
	s1 = s1 * ck - c2 * sk                                            
	go to 35                                                          
65	c1 = (c1 - s1) * sq                                               
	s1 = (c2 + s1) * sq                                               
	go to 35                                                          
70	continue                                                          
!                                  permute the complex vector in        
!                                    reverse binary order to normal     
!                                    order                              
	if(m .le. 1) go to 9005                                           
	mp = m+1                                                          
	jj = 1                                                            
!                                  initialize work vector               
	iwk(1) = 1                                                        
	do 75  i = 2,mp                                                   
	   iwk(i) = iwk(i-1) * 2                                          
75	continue                                                          
	n4 = iwk(mp-2)                                                    
	if (m .gt. 2) n8 = iwk(mp-3)                                      
	n2 = iwk(mp-1)                                                    
	lm = n2                                                           
	nn = iwk(mp)+1                                                    
	mp = mp-4                                                         
!                                  determine indices and switch a       
	j = 2                                                             
80	jk = jj + n2                                                      
	ak2 = a(j)                                                        
	a(j) = a(jk)                                                      
	a(jk) = ak2                                                       
	j = j+1                                                           
	if (jj .gt. n4) go to 85                                          
	jj = jj + n4                                                      
	go to 105                                                         
85	jj = jj - n4                                                      
	if (jj .gt. n8) go to 90                                          
	jj = jj + n8                                                      
	go to 105                                                         
90	jj = jj - n8                                                      
	k = mp                                                            
95	if (iwk(k) .ge. jj) go to 100                                     
	jj = jj - iwk(k)                                                  
	k = k - 1                                                         
	go to 95                                                          
100	jj = iwk(k) + jj                                                  
105	if (jj .le. j) go to 110                                          
	k = nn - j                                                        
	jk = nn - jj                                                      
	ak2 = a(j)                                                        
	a(j) = a(jj)                                                      
	a(jj) = ak2                                                       
	ak2 = a(k)                                                        
	a(k) = a(jk)                                                      
	a(jk) = ak2                                                       
110	j = j + 1                                                         
!                                  cycle repeated until limiting number 
!                                    of changes is achieved             
	if (j .le. lm) go to 80                                           
!                                                                       
9005	return                                                            
	end                                                               
