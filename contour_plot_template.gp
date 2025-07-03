#set term postscript eps enhanced color "Helvetica" 18
#set out 'IntMap.ps'

#for png
set term png size 850, 600
set out 'SlicedData.png'

set pm3d
unset surface
set view map

#set yrange[0:1]
#set xrange[0.5:1]
#set cbran[0:]

set palette defined (-8 "black", -4 "blue", -1 "red", 1 "orange", 5 "yellow", 10 "white")
#set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000') #jet
#set palette defined ( -1 '#000030', 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')  #jet2

#set encoding iso_8859_1
#set xlabel 'Q ({/Helvetica \305}^-^1)' font "Helvetica,22"    # (A^-1)
#set ylabel 'Normalized intensity (arb units)' font "Helvetica,22"

set xlabel 'Qx (A^{-1})'
set ylabel 'Qy (A^{-1})'
set cblabel 'Int'

splot "SlicedData.txt" u 1:2:($3-1)