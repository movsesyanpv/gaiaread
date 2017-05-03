set term png size 1920,1080 crop
#set term pdf
set term postscript eps enhanced size 8in, 3in
#set yrange[0:1]
#set output "healpmeddist.узы"
#plot "meddist.dat" u 1 w lp
set xrange[0:5000]
set xlabel "dist, dist in pc"
set ylabel "<relative distance error>"
set output "distvserr.eps"
plot "testout.csv" u 1:2 w lp
set ylabel "number of objects"
set output "distvsN.eps"
plot "testout.csv" u 1:3 w lp
