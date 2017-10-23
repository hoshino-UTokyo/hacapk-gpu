#include "m_ppohBEM_matrix_element_ij_cuda.inc"

  attributes(device) real*8 function HACApK_entry_ij_cuda(i, j, zbemv)
  integer*4 :: i,j
  type(st_HACApK_calc_entry) :: zbemv	
  
  HACApK_entry_ij_cuda=ppohBEM_matrix_element_ij_cuda(i, j, &
                                            zbemv%nond, zbemv%nofc, zbemv%nond_on_face, zbemv%np, &
                                            zbemv%int_para_fc, zbemv%nint_para_fc, zbemv%dble_para_fc, &
                                            zbemv%ndble_para_fc, zbemv%face2node)
  
  if(zbemv%lp61==3) HACApK_entry_ij_cuda=HACApK_entry_ij_cuda*zbemv%ao(i)*zbemv%ao(j)

end function HACApK_entry_ij_cuda

  attributes(device) real*8 function HACApK_entry_ij_cuda2(i, j, zbemv)
  integer*4 :: i,j
  type(st_HACApK_calc_entry) :: zbemv	
  
  HACApK_entry_ij_cuda2=ppohBEM_matrix_element_ij_cuda2(i, j, &
                                            zbemv%nond, zbemv%nofc, zbemv%nond_on_face, zbemv%np, &
                                            zbemv%int_para_fc, zbemv%nint_para_fc, zbemv%dble_para_fc, &
                                            zbemv%ndble_para_fc, zbemv%face2node)
  
  if(zbemv%lp61==3) HACApK_entry_ij_cuda2=HACApK_entry_ij_cuda2*zbemv%ao(i)*zbemv%ao(j)

end function HACApK_entry_ij_cuda2
