set term png size 800,600 crop
#set term pdf
#set term postscript eps enhanced size 8in, 3in
#set yrange[0:1]
#set output "healpmeddist.узы"
#plot "meddist.dat" u 1 w lp
#set xrange[0:5000]
#set xlabel "dist, dist in pc"
#set ylabel "<relative distance error>"
#set output "distvserr.eps"
#plot "testout.csv" u 1:2 w lp
#set ylabel "number of objects"
#set output "distvsN.eps"
#plot "testout.csv" u 1:3 w lp
#set yrange[5:30]
set output "./graphs/UVW500.png"
plot "500.dat" u 1:2:3 w yerr t "U" ,\
	 "500.dat" u 1:4:5 w yerr t "V" ,\
	 "500.dat" u 1:6:7 w yerr t "W"

#set yrange[-7:4]
set output "./graphs/Omega500.png"
plot "500.dat" u 1:8:9 w yerr t "Wx" ,\
	 "500.dat" u 1:10:11 w yerr t "Wy" ,\
	 "500.dat" u 1:12:13 w yerr t "Wz"

#set yrange[-1:10]
set output "./graphs/M500.png"
plot "500.dat" u 1:14:15 w yerr t "M13+" ,\
	 "500.dat" u 1:16:17 w yerr t "M23+" ,\
	 "500.dat" u 1:18:19 w yerr t "M12+" ,\
	 "500.dat" u 1:20:21 w yerr t "M11*" ,\
	 "500.dat" u 1:22:23 w yerr t "M33"