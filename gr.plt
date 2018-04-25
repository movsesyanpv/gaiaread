#set term png size 800,600 crop
#set term pdf
set term postscript eps enhanced size 8in, 6in
#set yrange[0:1]
#set output "healpmeddist.узы"
#plot "meddist.dat" u 1 w lp
#set xrange[0:1000]
#set xlabel "dist, dist in pc"
#set ylabel "<relative distance error>"
#set output "distvserr.eps"
#plot "testout.csv" u 1:2 w lp
#set ylabel "number of objects"
#set output "distvsN.eps"
#plot "testout.csv" u 1:3 w lp
#set yrange[5:30]
#set output "./graphs/UVW100.png"
set output "./graphs/UVW100.eps"
plot "100.dat" u 1:2:3 w yerr t "U" ,\
	 '' using 1:2 w lines ls 1 notitle,\
	 "100.dat" u 1:4:5 w yerr t "V" ,\
	 '' using 1:4 w lines ls 1 notitle,\
	 "100.dat" u 1:6:7 w yerr t "W",\
	 '' using 1:6 w lines ls 1 notitle

#set yrange[-7:4]
#set output "./graphs/Omega100.png"
set output "./graphs/Omega100.eps"
plot "100.dat" u 1:8:9 w yerr t "Wx" ,\
	 '' using 1:8 w lines ls 1 notitle,\
	 "100.dat" u 1:10:11 w yerr t "Wy" ,\
	 '' using 1:10 w lines ls 1 notitle,\
	 "100.dat" u 1:12:13 w yerr t "Wz",\
	 '' using 1:12 w lines ls 1 notitle

#set yrange[-1:10]
#set output "./graphs/M100.png"
set output "./graphs/M100.eps"
plot "100.dat" u 1:14:15 w yerr t "M13+" ,\
	 '' using 1:14 w lines ls 1 notitle,\
	 "100.dat" u 1:16:17 w yerr t "M23+" ,\
	 '' using 1:16 w lines ls 1 notitle,\
	 "100.dat" u 1:18:19 w yerr t "M12+" ,\
	 '' using 1:18 w lines ls 1 notitle,\
	 "100.dat" u 1:20:21 w yerr t "M11*" ,\
	 '' using 1:20 w lines ls 1 notitle,\
	 "100.dat" u 1:22:23 w yerr t "M33",\
	 '' using 1:22 w lines ls 1 notitle
