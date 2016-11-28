set term png size 1920,1080 crop
#set term pdf
#set term postscript landscape
set yrange[0:1]
set xrange[0:3000]
set xlabel "dist, dist in pc"
set ylabel "<relative distance error>"
set output "res.png"
plot "testout.csv" u 1:2 w lp
