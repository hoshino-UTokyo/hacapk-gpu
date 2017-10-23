#include "m_ppohBEM_user_func_cuda.inc"

  attributes(device) real(8) function ppohBEM_matrix_element_ij_cuda(i, j, nond, nofc, nond_on_fc, np,int_para_fc, nint_para_fc, dble_para_fc, ndble_para_fc, face2node)
  integer ,intent(in) :: i, j, nond, nofc, nond_on_fc, nint_para_fc, ndble_para_fc
  type(coordinate), intent(in),pointer :: np(:)
  integer, intent(in), pointer :: face2node(:, :), int_para_fc(:,:)
  real(8),pointer, intent(in) :: dble_para_fc(:, :)
  
  integer :: n(3)
  real(8) :: xf(3), yf(3), zf(3)
  real(8) :: xp, yp, zp

  n(1:3) = face2node(1:3, i) + 1
  xf(1:3) = np( n(1:3) )%x
  yf(1:3) = np( n(1:3) )%y
  zf(1:3) = np( n(1:3) )%z

  xp = (xf(1) + xf(2) + xf(3) ) / 3d0
  yp = (yf(1) + yf(2) + yf(3) ) / 3d0
  zp = (zf(1) + zf(2) + zf(3) ) / 3d0

  n(1:3) = face2node(1:3, j) + 1
  xf(1:3) = np( n(1:3) )%x
  yf(1:3) = np( n(1:3) )%y
  zf(1:3) = np( n(1:3) )%z
  
  ppohBEM_matrix_element_ij_cuda = face_integral_cuda(xf, yf, zf, xp, yp, zp)
end function ppohBEM_matrix_element_ij_cuda

  attributes(device) real(8) function ppohBEM_matrix_element_ij_cuda2(i, j, nond, nofc, nond_on_fc, np,int_para_fc, nint_para_fc, dble_para_fc, ndble_para_fc, face2node)
  integer ,intent(in) :: i, j, nond, nofc, nond_on_fc, nint_para_fc, ndble_para_fc
  type(coordinate), intent(in),pointer :: np(:)
  integer, intent(in), pointer :: face2node(:, :), int_para_fc(:,:)
  real(8),pointer, intent(in) :: dble_para_fc(:, :)
  
  integer :: n(3)
  real(8) :: xf(3), yf(3), zf(3)
  real(8) :: xp, yp, zp

  n(1:3) = face2node(1:3, i) + 1
  xf(1:3) = np( n(1:3) )%x
  yf(1:3) = np( n(1:3) )%y
  zf(1:3) = np( n(1:3) )%z

  xp = (xf(1) + xf(2) + xf(3) ) / 3d0
  yp = (yf(1) + yf(2) + yf(3) ) / 3d0
  zp = (zf(1) + zf(2) + zf(3) ) / 3d0

  n(1:3) = face2node(1:3, j) + 1
  xf(1:3) = np( n(1:3) )%x
  yf(1:3) = np( n(1:3) )%y
  zf(1:3) = np( n(1:3) )%z
  
  ppohBEM_matrix_element_ij_cuda2 = face_integral_cuda2(xf, yf, zf, xp, yp, zp)

end function ppohBEM_matrix_element_ij_cuda2

