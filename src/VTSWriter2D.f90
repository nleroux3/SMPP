! Writes 2D outputs files that can be read by Paraview
subroutine VTSWriter(TimeStep,nx,ny,x,y,Quantity,folder,opt)
!  opt      : variable of type character :
!              - 'ini' for the first call of VTSWriter                                
!              - 'int' for a standard call of  VTSWriter                              
!              - 'end' for the last call of  VTSWriter                                
use Declarations
implicit none
integer, intent(in) :: TimeStep    ! Number of iteration
integer, intent(in) :: nx, &       ! Number of horizontal nodes
                    &  ny          ! Number of vertical nodes
real(dp), intent(in) :: Quantity(nx-1,ny-1), & ! Variable that will be written
        &   x(nx,ny), &                        ! x coordinate of each node
        &   y(nx,ny)                           ! y coordinate of each node
character(100) :: num2char
character(7), intent(in):: folder             ! Folder in which outputs are created
character(200) :: FileName, formatperso,FileName2
character(3), intent(in) :: opt
 

!========== Writting a temporal file with the format Paraview =========
if (TimeStep >= 0) then
   write(num2char,'(i6.6)') TimeStep
   FileName = 'sol_'//trim(num2char)//'.vts'
   FileName2 = './outputs/'//folder//'/sol_'//trim(num2char)//'.vts'
   open(8,file=FileName2)
else
   open(8,file='./outputs/'//folder//'/sol_exacte.vts')
end if
write(num2char,*) 3*nx*ny
formatperso = '('//trim(num2char)//'(F15.6,1x))'
write(8,'(a)') '<?xml version="1.0"?>'
write(8,'(a)') '<VTKFile type="StructuredGrid">'
write(8,'(a,6i6,a)') '<StructuredGrid WholeExtent="', 0,nx-1,0,ny-1,0,0,'">'
write(8,'(a,6i6,a)') '<Piece Extent="',0,nx-1,0,ny-1,0,0,'">'
write(8,'(a)') '<Points>'
write(8,'(a)') '<DataArray type="Float32" NumberOfComponents="3"/>'

!================ Writting the coordinates ============================
DO j=1,ny
   write(8,formatperso) (x(i,j),y(i,j),0.,i=1,nx)
END DO
write(8,'(a)') '</Points>'
write(8,'(a)') '<CellData Scalars="'//folder//'">'
write(8,'(a)') '<DataArray type="Float32" Name="'//folder//'"/>'
write(num2char,*) nx*ny
formatperso = '('//trim(num2char)//'(F15.6,1x))'

!================ Writting the values of Quantity ====================
DO j=1,ny-1
   write(8,formatperso) (Quantity(i,j),i=1,nx-1)
END DO
write(8,'(a)') '</CellData>'
write(8,'(a)') '</Piece>'
write(8,'(a)') '</StructuredGrid>'
write(8,'(a)') '</VTKFile>'
close(8)

!============= Filling the file "Collection" ===========================
if (opt == 'ini' ) then
  open(10,file='./outputs/'//folder//'/sol.pvd')
  write(10,'(a)') '<?xml version="1.0"?>'
  write(10,*) '<VTKFile type="Collection">'
  write(10,*) '<Collection>'
else
  open(10,file='./outputs/'//folder//'/sol.pvd',position='append')
end if
if (TimeStep >= 0) write(10,*) '<DataSet timestep="',TimeStep,'" group="" part="0" file="',trim(FileName),'"/>'
if ( opt == 'end') then
   write(10,*) '</Collection>'
   write(10,*) '</VTKFile>'
end if
close(10)

end subroutine VTSWriter

