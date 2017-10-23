#define KIND_REAL 8
#define PI (3.141592653589793)
#define EPSILON_0  (8.854187818 * 1d-12)
  attributes(device) real(kind = KIND_REAL) function face_integral_cuda(xs, ys, zs, x, y, z)
    real(kind = KIND_REAL), intent(in) :: xs(3), ys(3), zs(3)
    real(kind = KIND_REAL), intent(in) :: x, y, z
    
    real(kind = KIND_REAL) :: r(3)
    real(kind = KIND_REAL) :: xi, xj, yi, dx, dy, t, l, m, d, ti, tj
    real(kind = KIND_REAL) :: theta, omega, q, g, zp, zpabs
    
    integer :: i, j
    real(kind = KIND_REAL) :: u(3), v(3), w(3) 
    real(kind = KIND_REAL) :: ox, oy, oz 
    
    r(:) = sqrt( (xs(:) - x)**2 + (ys(:) - y)**2 + (zs(:) - z)**2 )
    
    u(1) = xs(2) - xs(1);  v(1) = xs(3) - xs(2)
    u(2) = ys(2) - ys(1);  v(2) = ys(3) - ys(2)
    u(3) = zs(2) - zs(1);  v(3) = zs(3) - zs(2)
    call cross_product_g(u, v, w)
    w(:) = w(:) / sqrt( dot_product(w, w) )
    
    u(1) = x - xs(1);  u(2) = y - ys(1);  u(3) = z - zs(1)
    zp = dot_product(u, w)
    ox = x - zp * w(1);  oy = y - zp * w(2);  oz = z - zp * w(3)
    zpabs = abs(zp)
    
    face_integral_cuda = 0d0
    do i = 1, 3
       j = mod(i, 3) + 1
       u(1) = xs(j) - ox;  u(2) = ys(j) - oy; u(3) = zs(j) - oz;
       xj = sqrt( dot_product(u, u) )
       u(:) = u(:) / xj
       call cross_product_g(w, u, v) 
       xi = (xs(i) - ox) * u(1) + (ys(i) - oy) * u(2) + (zs(i) - oz) * u(3)
       yi = (xs(i) - ox) * v(1) + (ys(i) - oy) * v(2) + (zs(i) - oz) * v(3)
       
       dx = xj - xi;  dy = - yi ! (yj = 0)
       t = sqrt( dx**2 + dy**2 )
       l = dx/t;  m = dy/t
       d = l * yi - m * xi
       ti = l * xi + m * yi;  tj = l * xj ! (yj = 0)
       
       !            theta = sign(1d0, yi) * acos( xi / sqrt( xi**2 + yi**2 ) ) ! (xj > 0, yj = 0)
       theta = atan2(yi, xi)
       omega = theta - atan2( r(i) * d, zpabs * ti ) + atan2( r(j) * d, zpabs * tj )
       q = log( (r(j) + tj) / ( r(i) + ti ) )
       g = d * q - zpabs * omega 
       face_integral_cuda = face_integral_cuda + g
    enddo
    
    face_integral_cuda = abs(face_integral_cuda) / (4d0 * PI * EPSILON_0)       
    
  end function face_integral_cuda

  attributes(device) real(kind = KIND_REAL) function face_integral_cuda2(xs, ys, zs, x, y, z)
    real(kind = KIND_REAL), intent(in) :: xs(3), ys(3), zs(3)
    real(kind = KIND_REAL), intent(in) :: x, y, z
    
    real(kind = KIND_REAL) :: r(3)
    real(kind = KIND_REAL) :: xi, xj, yi, dx, dy, t, l, m, d, ti, tj
    real(kind = KIND_REAL) :: theta, omega, q, g, zp, zpabs, dtmp
    
    integer :: i, j
    real(kind = KIND_REAL) :: u(3), v(3), w(3) 
    real(kind = KIND_REAL) :: ox, oy, oz 
    
    r(:) = sqrt( (xs(:) - x)**2 + (ys(:) - y)**2 + (zs(:) - z)**2 )
    
    u(1) = xs(2) - xs(1);  v(1) = xs(3) - xs(2)
    u(2) = ys(2) - ys(1);  v(2) = ys(3) - ys(2)
    u(3) = zs(2) - zs(1);  v(3) = zs(3) - zs(2)
    call cross_product_g(u, v, w)
    ! w(1) = u(2) * v(3) - u(3) * v(2)
    ! w(2) = u(3) * v(1) - u(1) * v(3)
    ! w(3) = u(1) * v(2) - u(2) * v(1)

    dtmp = 1 / sqrt( dot_product(w, w) )
    w(:) = w(:) * dtmp 
    
    u(1) = x - xs(1);  u(2) = y - ys(1);  u(3) = z - zs(1)
    zp = dot_product(u, w)
    ox = x - zp * w(1);  oy = y - zp * w(2);  oz = z - zp * w(3)
    zpabs = abs(zp)
    
    face_integral_cuda2 = 0d0
    do i = 1, 3
       j = mod(i, 3) + 1
       u(1) = xs(j) - ox;  u(2) = ys(j) - oy; u(3) = zs(j) - oz;
       xj = sqrt( dot_product(u, u) )
       dtmp = 1 / xj
       u(:) = u(:) * dtmp
       call cross_product_g(w, u, v) 
       ! v(1) = w(2) * u(3) - w(3) * u(2)
       ! v(2) = w(3) * u(1) - w(1) * u(3)
       ! v(3) = w(1) * u(2) - w(2) * u(1)
       xi = (xs(i) - ox) * u(1) + (ys(i) - oy) * u(2) + (zs(i) - oz) * u(3)
       yi = (xs(i) - ox) * v(1) + (ys(i) - oy) * v(2) + (zs(i) - oz) * v(3)
       
       dx = xj - xi;  dy = - yi ! (yj = 0)
       t = 1/sqrt( dx**2 + dy**2 )
       l = dx*t;  m = dy*t
       d = l * yi - m * xi
       ti = l * xi + m * yi;  tj = l * xj ! (yj = 0)
       
       !            theta = sign(1d0, yi) * acos( xi / sqrt( xi**2 + yi**2 ) ) ! (xj > 0, yj = 0)
       theta = atan2(yi, xi)
       omega = theta - atan2( r(i) * d, zpabs * ti ) + atan2( r(j) * d, zpabs * tj )
       q = log( (r(j) + tj) / ( r(i) + ti ) )
       g = d * q - zpabs * omega 
       face_integral_cuda2 = face_integral_cuda2 + g
    enddo
    
    face_integral_cuda2 = abs(face_integral_cuda2) / (4d0 * PI * EPSILON_0)       
    
  end function face_integral_cuda2

!***********************************************************************
  attributes(device) subroutine cross_product_g(u, v, w)
    real(kind = KIND_REAL) :: u(3), v(3), w(3)
    
    w(1) = u(2) * v(3) - u(3) * v(2)
    w(2) = u(3) * v(1) - u(1) * v(3)
    w(3) = u(1) * v(2) - u(2) * v(1)
    
  end subroutine cross_product_g
#undef KIND_REAL 
#undef PI
#undef EPSILON_0
