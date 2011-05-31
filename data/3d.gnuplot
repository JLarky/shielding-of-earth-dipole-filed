set output '3d.png'
set terminal png size 800,600
set view 80,105,1,1
splot './points.dat' w d t '3D'

set output '2d.png'
plot [20:-100][:] './points.dat' u 1:3 w d t 'X-Z plane'

set output 'XZ.png'
plot [20:-100][:] './XZ.dat' u 1:3 w d t 'X-Z plane, Y=0'