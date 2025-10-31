#!/bin/bash
shopt -s extglob
echo "script for evaluation of the dendrite tip radius"
echo "as input the dendrite contour for every timestep is needed"
echo "the script ceates a new dendriteproperties folder with the actual date-and timestamp"
echo "after the run of the script the are two residual folders"
echo "--> raw data files (result files of OpenPhase)"
echo "--> result files..."
echo "the format of the plot can be adjusted"
echo "--> x-axis minimum and maximum and x-axis tics"
echo "--> y-axis maximum and y-axis tics (y-axis minimum is set to 0 per default)"
echo "in the result files folder there are subfolders for contour plots,"
echo "fit-data, and a list with the tipradius for every timestep (combined_....txt)"
echo "as fitrange the intersections of the fit of the dendrite tip and trunk are used"
echo "(for the first fit step a fitrange of 1/4 the length of the x-axis is used)"
echo "a user defined fit range can be entered later on"
echo "an automated gif of all plots for ongoing timesteps can be enabled later on"
echo "-----------------------------------------------------------------------------"
echo "name of the simulation (e.g. eta5dx1)"
read simname
#echo "intervallbreite fit trunk"
#read intertrunk
#negintertrunk=$($intertrunk*-1)
echo "x-axis minimum"
read xmin
echo "x-axis maximum"
read xmax
echo "x-axis tics"
read teilstrichex
echo "y-axis maximun"
read ymax
echo "y-axis tics"
read teilstrichey
echo "automated fit range (a) or user defined fit-range (u)"
read intertipgen
###########################
if [[ $intertipgen = a ]] ; then
##########################
  echo "automated contour of contour plots (y=yes;n=o)?"
  read video
  ############################################
  if [[ $video = y ]] ; then
  ############################################

  stamptime=`date +%d-%m-%Y_%H:%M`
  cp -rv DendriteProperties DendriteProperties_"$stamptime"
  #>>>>>>>>>>>>>>>>>
  cd DendriteProperties_"$stamptime"
  mkdir results_$simname
  mkdir results_$simname/rawdata_dendrcontour-files
  find . -type f -name "DendrContour_*.txt" -exec cp {} ./results_$simname/rawdata_dendrcontour-files/ \;
  #>>>>>>>>>>>>>>>>>
  cd results_$simname/rawdata_dendrcontour-files/
   for b in $(ls)
   do  
   mv $b $(echo $b | awk -F '[_]' '{printf "%07d" ,$2}')_halbe_dendritenkontur.txt
   done
   #ls -I to ignore special pattern
   ls -I "liste*.txt" > liste_halbe_dendritenkontur.txt
   ##create complete dendrite contour
    for d in `cat liste_halbe_dendritenkontur.txt`
     do
      timestep=$(echo $d | awk -F '_' '{print $1}')
      contvalue2=$(cat $d | awk '{print $2*-1}')
      cat $d | awk '{print $1}' > liste_zcoord_half_"$timestep".txt
      cat $d | awk '{print $1}' > liste_zcoord_half2_"$timestep".txt
      cat $d | awk '{print $2}' > liste_contour_coord_pos_"$timestep".txt
      cat $d | awk '{print $2*-1}' > liste_contour_coord_neg_"$timestep".txt
      #add one point in -x direction to get correct tip geometry
      cat liste_contour_coord_neg_"$timestep".txt | sed -e '1s/0/-1/' > liste_contour_coord_neg2_"$timestep".txt
      mv liste_contour_coord_neg2_"$timestep".txt liste_contour_coord_neg_"$timestep".txt
      awk 'NF' liste_zcoord_half_"$timestep".txt liste_zcoord_half2_"$timestep".txt > z_coord_komp_"$timestep".txt
      awk 'NF' liste_contour_coord_pos_"$timestep".txt liste_contour_coord_neg_"$timestep".txt > contour_coord_komp_"$timestep".txt
    paste z_coord_komp_"$timestep".txt contour_coord_komp_"$timestep".txt > "$timestep"_komplette_dendriten_contour.txt
    ##finished complete dendrite contour
    mkdir list-files
    find . -type f -name "liste_zcoord_half_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "liste_zcoord_half2_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "liste_contour_coord_pos_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "liste_contour_coord_neg_*.txt" -exec mv {} ./list-files \; 
    find . -type f -name "liste_contour_coord_neg2_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "z_coord_komp_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "contour_coord_komp_*.txt" -exec mv {} ./list-files \;
   done
  find . -type f -name "*_komplette_dendriten_contour.txt" > list_komplette_dendritenkontur.txt
  cat list_komplette_dendritenkontur.txt | cut -c 3- | sort -n > list_komplette_dendritenkontur2.txt
  mv list_komplette_dendritenkontur2.txt list_komplette_dendritenkontur.txt
   for y in `cat list_komplette_dendritenkontur.txt`
   do
   z=$(echo $y | awk -F'[_]' '{print $1}')
   yposmax=$(cat $y | awk 'NR==1{print $1}')
   yposmin=$(cat $y | awk 'NR==1{print $1-20}')
   echo "$yposmin" >> list_yposmax.txt
    gnuplot -persist <<- EOF
     reset
     set terminal pngcairo size 1000,1000 enhanced font "font Arial,22"
     set key tmargin
     set key font ",14"
     #set style line 5 lw 8 lc rgb '#F0F8FF' pt 4
     set style line 3 lw 5 lc rgb '#00ff00'
     set style line 4 lw 5 lc rgb '#00008B'
     set xlabel "Dendrite radius [{\/Symbol m}m]"
     set ylabel "Z-coordinate [-]"
     set xrange ["${xmin}":"${xmax}"]
     set yrange [0:"${ymax}"]
     #set arrow 1 from -"${intertip}","${yposmin}" to -"${intertip}","${yposmax}" nohead lw 4 lc rgb "#FF0000" dt "--"
     #set arrow 2 from  "${intertip}","${yposmin}" to  "${intertip}","${yposmax}" nohead lw 4 lc rgb "#FF0000" dt "--"
  #set object 11 circle at axis 0,"${yposcircle}" size "${intertip}"
  #set output "contour_timestep_${z}.png"
  set fit logfile "fit_data_full_contour_timestep_${z}.txt"
  f(x) = a*x**2+b*x+c
  fit f(x) "${y}" using 2:1 via a, b, c
  f2(x)= a2*x**2+b2*x+c2
  fit ["-${intertip}":"${intertip}"] f2(x) "${y}" using 2:1 via a2, b2, c2
  set print "fit_parameters_${z}.txt"
  print a
  print a2
  print b
  print d
  print c
  print f
  set print
  #plot "${y}" using 2:1 with points ps 5 title "timestep${z}", f(x) ls 3 title "parabolic fit trunk", f2(x) ls 4 title "parabolic fit tip", 1/0 lw 4 lc rgb "#FF0000" dt "--" title "fit range" 

EOF

find . -type f -name "fit_data_full_contour_timestep_*.txt" -exec rm {} \;
fita=$(cat fit_parameters_$z | awk 'NR==1{print ($1*-1)}')
fitb=$(cat fit_parameters_$z | awk 'NR==2{print ($1*-1)}')
fitc=$(cat fit_parameters_$z | awk 'NR==3{print $1}')
fitd=$(cat fit_parameters_$z | awk 'NR==4{print $1}')
fite=$(cat fit_parameters_$z | awk 'NR==5{print $1}')
fitf=$(cat fit_parameters_$z | awk 'NR==6{print $1}')
fitintervalpos=$(echo "((($fitb-$fite)/($fita-$fitd)2)+ sqrt((($fitb-$fite)/($fita-$fitd))/2)*(($fitb-$fite)/($fita-$fitd))/2)+(($fitc-$fitf)/($fita-$fitd))" | bc)
labelfitrangey=$(echo "ymax" -5 |bc)
labelfitrangex=$(echo "xmax" -20|bc)
####start gnuplot-script
gnuplot -persist <<- EOF
      reset
      firstfitrange = "${xmax}" * 0.25
      set terminal pngcairo size 1000,1000 enhanced font "font Arial,22"
      set key tmargin
      set key font ",14"
      #set style line 5 lw 8 lc rgb '#F0F8FF' pt 4
      set style line 3 lw 5 lc rgb '#00ff00'
      set style line 4 lw 5 lc rgb '#00008B'
      set xlabel "Dendrite radius [{\/Symbol m}m]"
      set ylabel "Z-coordinate [-]"
      set xrange ["${xmin}":"${xmax}"]
      set yrange [0:"${ymax}"]
      set output "contour_timestep_${z}.png"
      #### first fit with a total range and 0.25% of xmax
      set fit logfile "fit_data_full_contour_timestep_1step_${z}.txt"
      f(x) = a*x**2+b*x+c
      fit f(x) "${y}" using 2:1 via a, b, c
      f2(x)= a2*x**2+b2*x+c2
      fit [-firstfitrange:firstfitrange] f2(x) "${y}" using 2:1 via a2, b2, c2
      #### end of fit and save of fitparameter in txt-file
      set print "fit_parameters_1step_${z}.txt"
      print a
      print a2
      print b
      print b2
      print c
      print c2
      set print
      #### variables for analytical solution of intersection of two parabolas
      fita = a*-1
      fitb = b
      fitc = c
      fitd = a2*-1
      fite = b2
      fitf = c2
      #### analytical solution of parabola intersection
      fitinterval= (((fitb-fite)/(fita-fitd))/2)+ sqrt((((fitb-fite)/(fita-fitd))/2)*(((fitb-fite)/(fita-fitd))/2)+((fitc-fitf)/(fita-fitd)))
      #### second fit in range of analytical solution pos and neg
      set fit logfile "fit_data_full_contour_timestep_final_${z}.txt"
      f3(x)= h*x**2+i*x+j
      fit [-firstfitrange:firstfitrange] f3(x) "${y}" using 2:1 via h, i, j
      #### save the fitparameters for fit for dendrite tip radius calculation
      set print "fit_parameters_final_${z}.txt"
      print h
      print i
      print j
      #### creation of the final plot
      set arrow 1 from -fitinterval,"${yposmin}" to -fitinterval,"${yposmax}" nohead lw 4 lc rgb "#FF0000" dt "--"
      set arrow 2 from  fitinterval,"${yposmin}" to  fitinterval,"${yposmax}" nohead lw 4 lc rgb "#FF0000" dt "--"
      plot "${y}" using 2:1 with points ps 5 title "timestep${z}", f(x) ls 3 title "parabolic fit trunk", f2(x) ls 4 title "parabolic fit tip", 1/0 lw 4 lc rgb "#FF0000" dt "--" title "fit range" 
EOF
done
   #creation of results subfolders and data migration
   mkdir fit_data
   mkdir fit_data/fitparameter
   mkdir contour_plots
   mkdir files_volle_dendritecontour
   mkdir files_halbe_dendritecontour
   find . -type f -name "fit_data_full_contour_timestep_final_*.txt" -exec mv {} ./fit_data/ \;
   find . -type f -name "fit_data_full_contour_timestep_1step_*.txt" -exec mv {} ./fit_data/ \;
   find . -type f -name "contour_timestep_*.png" -exec mv {} ./contour_plots/ \;
   find . -type f -name "*_komplette_dendriten_contour.txt" -exec mv {} ./files_volle_dendritecontour/ \;
   find . -type f -name "*_halbe_dendritenkontur.txt" -exec mv {} ./files_halbe_dendritecontour/ \;
   find . -type f -name "fit_parameters_final_*" -exec mv {} ./fit_data/fitparameter/ \;
   find . -type f -name "fit_parameters_1step_*" -exec mv {} ./fit_data/fitparameter/ \;
   find . -type f -name "DendrContour_*.txt" -exec mv {} ./rawdata_dendritencontour/ \;
   mkdir files_tipradius
   find . -type f -name "DendrTipRadius_*.txt" -exec mv {} ./files_tipradius/ \;
   find . -type f -name "list_komplette_dendritenkontur.txt" -exec mv {} ./list-files/ \;
   #<<<<<<<<<<<<<
   cd ..
   cp -r rawdata_dendrcontour-files/contour_plots .
   cp -r rawdata_dendrcontour-files/files_volle_dendritecontour .
   cp -r rawdata_dendrcontour-files/files_halbe_dendritecontour .
   cp -r rawdata_dendrcontour-files/rawdata_dendritencontour .
   cp -r rawdata_dendrcontour-files/fit_data .
   cp -r rawdata_dendrcontour-files/files_tipradius .
   cp -r rawdata_dendrcontour-files/list-files .
   cp -r rawdata_dendrcontour-files/list_yposmax.txt .
   rm -rf rawdata_dendrcontour-files
   #creation of gif of contour plots
   cd contour_plots
     for i in *
     do
     mv $i $(echo $i | awk -F '[_.]' '{printf "%09d" ,$3}')_contour.png
     done
   ls | convert -delay $delay -loop $numberloops *_contour.png animation_dendrite_$simname.gif
   mv animation_dendrite_$simname.gif ../
   cd ..
### evaluation of dendrite tip radius
   cd fit_data/fitparameter/
   mkdir fitparameter_sorted
   find . -type f -name "fit_parameters_final_*.txt" -exec cp {} ./fitparameter_sorted/ \;
   cd fitparameter_sorted
     for k in $(ls)
     do
     onlynumber=$(echo "$k"|awk -F'[_]' '{print $4}' | awk -F'[.]' '{printf "%06d\n", $1}')
     mv "$k" "$onlynumber"_gnuplotfitparameters.txt
     done
   find . -name "*_gnuplotfitparameters.txt" > listoffiles.txt
   cat listoffiles.txt | cut -c 3- | sort -n > listoffiles2.txt
   rm listoffiles.txtecho ""
   mv listoffiles2.txt listoffiles.txt
     for j in $(find . -name "*_gnuplotfitparameters.txt")
     do 
     echo $j | awk -F'[_]' '{print $1}' >> listtimesteps.txt
     done
   cat listtimesteps.txt | cut -c 3- | sort -n > listtimesteps2.txt
   rm listtimesteps.txt
   mv listtimesteps2.txt listtimesteps.txt
     for l in `cat listoffiles.txt`
     do
     cat $l | sed -n '1p' >> list_fitparameter_tip.txt
     done
   cat list_fitparameter_tip.txt | awk  '{print 1/(-1*2*$0)}' > list_radius_tip.txt
   paste listtimesteps.txt list_radius_tip.txt > list_combined_timesteps_radius.txt
   cd ../../../
   cp  fit_data/fitparameter/fitparameter_sorted/list_combined_timesteps_radius.txt .
   cd ..
   mkdir rawdata_output_OP
   find . -type f -name "DendrTipRadius_*.txt" -exec mv {} ./rawdata_output_OP/ \;
   find . -type f -name "DendrContour_*.txt" -exec mv {} ./rawdata_output_OP/ \;
   
   ####################################################
   elif [[ $video = n ]] ; then
   ####################################################
   
  stamptime=`date +%d-%m-%Y_%H:%M`
  cp -rv DendriteProperties DendriteProperties_"$stamptime"
  #>>>>>>>>>>>>>>>>>
  cd DendriteProperties_"$stamptime"
  mkdir results_$simname
  mkdir results_$simname/rawdata_dendrcontour-files
  find . -type f -name "DendrContour_*.txt" -exec cp {} ./results_$simname/rawdata_dendrcontour-files/ \;
  #>>>>>>>>>>>>>>>>>
  cd results_$simname/rawdata_dendrcontour-files/
   for b in $(ls)
   do  
   mv $b $(echo $b | awk -F '[_]' '{printf "%07d" ,$2}')_halbe_dendritenkontur.txt
   done
   #ls -I to ignore special pattern
   ls -I "liste*.txt" > liste_halbe_dendritenkontur.txt
   ##create complete dendrite contour
    for d in `cat liste_halbe_dendritenkontur.txt`
     do
      timestep=$(echo $d | awk -F '_' '{print $1}')
      contvalue2=$(cat $d | awk '{print $2*-1}')
      cat $d | awk '{print $1}' > liste_zcoord_half_"$timestep".txt
      cat $d | awk '{print $1}' > liste_zcoord_half2_"$timestep".txt
      cat $d | awk '{print $2}' > liste_contour_coord_pos_"$timestep".txt
      cat $d | awk '{print $2*-1}' > liste_contour_coord_neg_"$timestep".txt
      #add one point in -x direction to get correct tip geometry
      cat liste_contour_coord_neg_"$timestep".txt | sed -e '1s/0/-1/' > liste_contour_coord_neg2_"$timestep".txt
      mv liste_contour_coord_neg2_"$timestep".txt liste_contour_coord_neg_"$timestep".txt
      awk 'NF' liste_zcoord_half_"$timestep".txt liste_zcoord_half2_"$timestep".txt > z_coord_komp_"$timestep".txt
      awk 'NF' liste_contour_coord_pos_"$timestep".txt liste_contour_coord_neg_"$timestep".txt > contour_coord_komp_"$timestep".txt
    paste z_coord_komp_"$timestep".txt contour_coord_komp_"$timestep".txt > "$timestep"_komplette_dendriten_contour.txt
    ##finished complete dendrite contour
    mkdir list-files
    find . -type f -name "liste_zcoord_half_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "liste_zcoord_half2_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "liste_contour_coord_pos_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "liste_contour_coord_neg_*.txt" -exec mv {} ./list-files \; 
    find . -type f -name "liste_contour_coord_neg2_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "z_coord_komp_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "contour_coord_komp_*.txt" -exec mv {} ./list-files \;
   done
  find . -type f -name "*_komplette_dendriten_contour.txt" > list_komplette_dendritenkontur.txt
  cat list_komplette_dendritenkontur.txt | cut -c 3- | sort -n > list_komplette_dendritenkontur2.txt
  mv list_komplette_dendritenkontur2.txt list_komplette_dendritenkontur.txt
   for y in `cat list_komplette_dendritenkontur.txt`
   do
   z=$(echo $y | awk -F'[_]' '{print $1}')
   yposmax=$(cat $y | awk 'NR==1{print $1}')
   yposmin=$(cat $y | awk 'NR==1{print $1-20}')
   echo "$yposmin" >> list_yposmax.txt
    gnuplot -persist <<- EOF
     reset
     set terminal pngcairo size 1000,1000 enhanced font "font Arial,22"
     set key tmargin
     set key font ",14"
     #set style line 5 lw 8 lc rgb '#F0F8FF' pt 4
     set style line 3 lw 5 lc rgb '#00ff00'
     set style line 4 lw 5 lc rgb '#00008B'
     set xlabel "Dendrite radius [{\/Symbol m}m]"
     set ylabel "Z-coordinate [-]"
     set xrange ["${xmin}":"${xmax}"]
     set yrange [0:"${ymax}"]
     #set arrow 1 from -"${intertip}","${yposmin}" to -"${intertip}","${yposmax}" nohead lw 4 lc rgb "#FF0000" dt "--"
     #set arrow 2 from  "${intertip}","${yposmin}" to  "${intertip}","${yposmax}" nohead lw 4 lc rgb "#FF0000" dt "--"
  #set object 11 circle at axis 0,"${yposcircle}" size "${intertip}"
  #set output "contour_timestep_${z}.png"
  set fit logfile "fit_data_full_contour_timestep_${z}.txt"
  f(x) = a*x**2+b*x+c
  fit f(x) "${y}" using 2:1 via a, b, c
  f2(x)= a2*x**2+b2*x+c2
  fit ["-${intertip}":"${intertip}"] f2(x) "${y}" using 2:1 via a2, b2, c2
  set print "fit_parameters_${z}.txt"
  print a
  print a2
  print b
  print d
  print c
  print f
  set print
  #plot "${y}" using 2:1 with points ps 5 title "timestep${z}", f(x) ls 3 title "parabolic fit trunk", f2(x) ls 4 title "parabolic fit tip", 1/0 lw 4 lc rgb "#FF0000" dt "--" title "fit range" 

EOF

find . -type f -name "fit_data_full_contour_timestep_*.txt" -exec rm {} \;
fita=$(cat fit_parameters_$z | awk 'NR==1{print ($1*-1)}')
fitb=$(cat fit_parameters_$z | awk 'NR==2{print ($1*-1)}')
fitc=$(cat fit_parameters_$z | awk 'NR==3{print $1}')
fitd=$(cat fit_parameters_$z | awk 'NR==4{print $1}')
fite=$(cat fit_parameters_$z | awk 'NR==5{print $1}')
fitf=$(cat fit_parameters_$z | awk 'NR==6{print $1}')
fitintervalpos=$(echo "((($fitb-$fite)/($fita-$fitd)2)+ sqrt((($fitb-$fite)/($fita-$fitd))/2)*(($fitb-$fite)/($fita-$fitd))/2)+(($fitc-$fitf)/($fita-$fitd))" | bc)
labelfitrangey=$(echo "ymax" -5 |bc)
labelfitrangex=$(echo "xmax" -20|bc)
####start gnuplot-script
gnuplot -persist <<- EOF
      reset
      firstfitrange = "${xmax}" * 0.25
      set terminal pngcairo size 1000,1000 enhanced font "font Arial,22"
      set key tmargin
      set key font ",14"
      #set style line 5 lw 8 lc rgb '#F0F8FF' pt 4
      set style line 3 lw 5 lc rgb '#00ff00'
      set style line 4 lw 5 lc rgb '#00008B'
      set xlabel "Dendrite radius [{\/Symbol m}m]"
      set ylabel "Z-coordinate [-]"
      set xrange ["${xmin}":"${xmax}"]
      set yrange [0:"${ymax}"]
      set output "contour_timestep_${z}.png"
      #### first fit with a total range and 0.25% of xmax
      set fit logfile "fit_data_full_contour_timestep_1step_${z}.txt"
      f(x) = a*x**2+b*x+c
      fit f(x) "${y}" using 2:1 via a, b, c
      f2(x)= a2*x**2+b2*x+c2
      fit [-firstfitrange:firstfitrange] f2(x) "${y}" using 2:1 via a2, b2, c2
      #### end of fit and save of fitparameter in txt-file
      set print "fit_parameters_1step_${z}.txt"
      print a
      print a2
      print b
      print b2
      print c
      print c2
      set print
      #### variables for analytical solution of intersection of two parabolas
      fita = a*-1
      fitb = b
      fitc = c
      fitd = a2*-1
      fite = b2
      fitf = c2
      #### analytical solution of parabola intersection
      fitinterval= (((fitb-fite)/(fita-fitd))/2)+ sqrt((((fitb-fite)/(fita-fitd))/2)*(((fitb-fite)/(fita-fitd))/2)+((fitc-fitf)/(fita-fitd)))
      #### second fit in range of analytical solution pos and neg
      set fit logfile "fit_data_full_contour_timestep_final_${z}.txt"
      f3(x)= h*x**2+i*x+j
      fit [-firstfitrange:firstfitrange] f3(x) "${y}" using 2:1 via h, i, j
      #### save the fitparameters for fit for dendrite tip radius calculation
      set print "fit_parameters_final_${z}.txt"
      print h
      print i
      print j
      #### creation of the final plot
      set arrow 1 from -fitinterval,"${yposmin}" to -fitinterval,"${yposmax}" nohead lw 4 lc rgb "#FF0000" dt "--"
      set arrow 2 from  fitinterval,"${yposmin}" to  fitinterval,"${yposmax}" nohead lw 4 lc rgb "#FF0000" dt "--"
      plot "${y}" using 2:1 with points ps 5 title "timestep${z}", f(x) ls 3 title "parabolic fit trunk", f2(x) ls 4 title "parabolic fit tip", 1/0 lw 4 lc rgb "#FF0000" dt "--" title "fit range" 
EOF
done
   #creation of results subfolders and data migration
   mkdir fit_data
   mkdir fit_data/fitparameter
   mkdir contour_plots
   mkdir files_volle_dendritecontour
   mkdir files_halbe_dendritecontour
   find . -type f -name "fit_data_full_contour_timestep_final_*.txt" -exec mv {} ./fit_data/ \;
   find . -type f -name "fit_data_full_contour_timestep_1step_*.txt" -exec mv {} ./fit_data/ \;
   find . -type f -name "contour_timestep_*.png" -exec mv {} ./contour_plots/ \;
   find . -type f -name "*_komplette_dendriten_contour.txt" -exec mv {} ./files_volle_dendritecontour/ \;
   find . -type f -name "*_halbe_dendritenkontur.txt" -exec mv {} ./files_halbe_dendritecontour/ \;
   find . -type f -name "fit_parameters_final_*" -exec mv {} ./fit_data/fitparameter/ \;
   find . -type f -name "fit_parameters_1step_*" -exec mv {} ./fit_data/fitparameter/ \;
   find . -type f -name "DendrContour_*.txt" -exec mv {} ./rawdata_dendritencontour/ \;
   mkdir files_tipradius
   find . -type f -name "DendrTipRadius_*.txt" -exec mv {} ./files_tipradius/ \;
   find . -type f -name "list_komplette_dendritenkontur.txt" -exec mv {} ./list-files/ \;
   #<<<<<<<<<<<<<
   cd ..
   cp -r rawdata_dendrcontour-files/contour_plots .
   cp -r rawdata_dendrcontour-files/files_volle_dendritecontour .
   cp -r rawdata_dendrcontour-files/files_halbe_dendritecontour .
   cp -r rawdata_dendrcontour-files/rawdata_dendritencontour .
   cp -r rawdata_dendrcontour-files/fit_data .
   cp -r rawdata_dendrcontour-files/files_tipradius .
   cp -r rawdata_dendrcontour-files/list-files .
   cp -r rawdata_dendrcontour-files/list_yposmax.txt .
   rm -rf rawdata_dendrcontour-files
   ####evaluation of dendrite tip radius
   cd fit_data/fitparameter/
   mkdir fitparameter_sorted
   find . -type f -name "fit_parameters_final_*.txt" -exec cp {} ./fitparameter_sorted/ \;
   cd fitparameter_sorted
     for k in $(ls)
     do
     onlynumber=$(echo "$k"|awk -F'[_]' '{print $4}' | awk -F'[.]' '{printf "%06d\n", $1}')
     mv "$k" "$onlynumber"_gnuplotfitparameters.txt
     done
   find . -name "*_gnuplotfitparameters.txt" > listoffiles.txt
   cat listoffiles.txt | cut -c 3- | sort -n > listoffiles2.txt
   rm listoffiles.txtecho ""
   mv listoffiles2.txt listoffiles.txt
     for j in $(find . -name "*_gnuplotfitparameters.txt")
     do 
     echo $j | awk -F'[_]' '{print $1}' >> listtimesteps.txt
     done
   cat listtimesteps.txt | cut -c 3- | sort -n > listtimesteps2.txt
   rm listtimesteps.txt
   mv listtimesteps2.txt listtimesteps.txt
     for l in `cat listoffiles.txt`
     do
     cat $l | sed -n '1p' >> list_fitparameter_tip.txt
     done
   cat list_fitparameter_tip.txt | awk  '{print 1/(-1*2*$0)}' > list_radius_tip.txt
   paste listtimesteps.txt list_radius_tip.txt > list_combined_timesteps_radius.txt
   cd ../../../
   cp  fit_data/fitparameter/fitparameter_sorted/list_combined_timesteps_radius.txt .
   cd ..
   mkdir rawdata_output_OP
   find . -type f -name "DendrTipRadius_*.txt" -exec mv {} ./rawdata_output_OP/ \;
   find . -type f -name "DendrContour_*.txt" -exec mv {} ./rawdata_output_OP/ \;
   
   #########################################
   fi
   #########################################
   
############################################
elif [[ $intertipgen = u ]] ; then
###########################################
echo "a-xis minimum fit range"
read userdefxminfitpar
echo "a-xis maximum fit range"
read userdefxmaxfitpar   
echo "automated contour of contour plots (y=yes;n=o)?"
  read video
  echo "automated contour of contour plots (y=yes;n=o)?"
  read video
  ############################################
  if [[ $video = y ]] ; then
  ############################################

  stamptime=`date +%d-%m-%Y_%H:%M`
  cp -rv DendriteProperties DendriteProperties_"$stamptime"
  #>>>>>>>>>>>>>>>>>
  cd DendriteProperties_"$stamptime"
  mkdir results_$simname
  mkdir results_$simname/rawdata_dendrcontour-files
  find . -type f -name "DendrContour_*.txt" -exec cp {} ./results_$simname/rawdata_dendrcontour-files/ \;
  #>>>>>>>>>>>>>>>>>
  cd results_$simname/rawdata_dendrcontour-files/
   for b in $(ls)
   do  
   mv $b $(echo $b | awk -F '[_]' '{printf "%07d" ,$2}')_halbe_dendritenkontur.txt
   done
   #ls -I to ignore special pattern
   ls -I "liste*.txt" > liste_halbe_dendritenkontur.txt
   ##create complete dendrite contour
    for d in `cat liste_halbe_dendritenkontur.txt`
     do
      timestep=$(echo $d | awk -F '_' '{print $1}')
      contvalue2=$(cat $d | awk '{print $2*-1}')
      cat $d | awk '{print $1}' > liste_zcoord_half_"$timestep".txt
      cat $d | awk '{print $1}' > liste_zcoord_half2_"$timestep".txt
      cat $d | awk '{print $2}' > liste_contour_coord_pos_"$timestep".txt
      cat $d | awk '{print $2*-1}' > liste_contour_coord_neg_"$timestep".txt
      #add one point in -x direction to get correct tip geometry
      cat liste_contour_coord_neg_"$timestep".txt | sed -e '1s/0/-1/' > liste_contour_coord_neg2_"$timestep".txt
      mv liste_contour_coord_neg2_"$timestep".txt liste_contour_coord_neg_"$timestep".txt
      awk 'NF' liste_zcoord_half_"$timestep".txt liste_zcoord_half2_"$timestep".txt > z_coord_komp_"$timestep".txt
      awk 'NF' liste_contour_coord_pos_"$timestep".txt liste_contour_coord_neg_"$timestep".txt > contour_coord_komp_"$timestep".txt
    paste z_coord_komp_"$timestep".txt contour_coord_komp_"$timestep".txt > "$timestep"_komplette_dendriten_contour.txt
    ##finished complete dendrite contour
    mkdir list-files
    find . -type f -name "liste_zcoord_half_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "liste_zcoord_half2_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "liste_contour_coord_pos_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "liste_contour_coord_neg_*.txt" -exec mv {} ./list-files \; 
    find . -type f -name "liste_contour_coord_neg2_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "z_coord_komp_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "contour_coord_komp_*.txt" -exec mv {} ./list-files \;
   done
  find . -type f -name "*_komplette_dendriten_contour.txt" > list_komplette_dendritenkontur.txt
  cat list_komplette_dendritenkontur.txt | cut -c 3- | sort -n > list_komplette_dendritenkontur2.txt
  mv list_komplette_dendritenkontur2.txt list_komplette_dendritenkontur.txt
   for y in `cat list_komplette_dendritenkontur.txt`
   do
   z=$(echo $y | awk -F'[_]' '{print $1}')
   yposmax=$(cat $y | awk 'NR==1{print $1}')
   yposmin=$(cat $y | awk 'NR==1{print $1-20}')
   echo "$yposmin" >> list_yposmax.txt
    gnuplot -persist <<- EOF
     reset
     set terminal pngcairo size 1000,1000 enhanced font "font Arial,22"
     set key tmargin
     set key font ",14"
     #set style line 5 lw 8 lc rgb '#F0F8FF' pt 4
     set style line 3 lw 5 lc rgb '#00ff00'
     set style line 4 lw 5 lc rgb '#00008B'
     set xlabel "Dendrite radius [{\/Symbol m}m]"
     set ylabel "Z-coordinate [-]"
     set xrange ["${xmin}":"${xmax}"]
     set yrange [0:"${ymax}"]
     #set arrow 1 from -"${intertip}","${yposmin}" to -"${intertip}","${yposmax}" nohead lw 4 lc rgb "#FF0000" dt "--"
     #set arrow 2 from  "${intertip}","${yposmin}" to  "${intertip}","${yposmax}" nohead lw 4 lc rgb "#FF0000" dt "--"
  #set object 11 circle at axis 0,"${yposcircle}" size "${intertip}"
  #set output "contour_timestep_${z}.png"
  set fit logfile "fit_data_full_contour_timestep_${z}.txt"
  f(x) = a*x**2+b*x+c
  fit f(x) "${y}" using 2:1 via a, b, c
  f2(x)= a2*x**2+b2*x+c2
  fit ["-${userdefxminfitpar}":"${userdefxminfitpar}"] f2(x) "${y}" using 2:1 via a2, b2, c2
  set print "fit_parameters_${z}.txt"
  print a
  print a2
  #plot "${y}" using 2:1 with points ps 5 title "timestep${z}", f(x) ls 3 title "parabolic fit trunk", f2(x) ls 4 title "parabolic fit tip", 1/0 lw 4 lc rgb "#FF0000" dt "--" title "fit range" 


EOF
done
   #creation of results subfolders and data migration
   mkdir fit_data
   mkdir fit_data/fitparameter
   mkdir contour_plots
   mkdir files_volle_dendritecontour
   mkdir files_halbe_dendritecontour
   find . -type f -name "fit_data_full_contour_timestep_*.txt" -exec mv {} ./fit_data/ \;
   find . -type f -name "contour_timestep_*.png" -exec mv {} ./contour_plots/ \;
   find . -type f -name "*_komplette_dendriten_contour.txt" -exec mv {} ./files_volle_dendritecontour/ \;
   find . -type f -name "*_halbe_dendritenkontur.txt" -exec mv {} ./files_halbe_dendritecontour/ \;
   find . -type f -name "fit_parameters_*" -exec mv {} ./fit_data/fitparameter/ \;
   find . -type f -name "DendrContour_*.txt" -exec mv {} ./rawdata_dendritencontour/ \;
   mkdir files_tipradius
   find . -type f -name "DendrTipRadius_*.txt" -exec mv {} ./files_tipradius/ \;
   find . -type f -name "list_komplette_dendritenkontur.txt" -exec mv {} ./list-files/ \;
   #<<<<<<<<<<<<<
   cd ..
   cp -r rawdata_dendrcontour-files/contour_plots .
   cp -r rawdata_dendrcontour-files/files_volle_dendritecontour .
   cp -r rawdata_dendrcontour-files/files_halbe_dendritecontour .
   cp -r rawdata_dendrcontour-files/rawdata_dendritencontour .
   cp -r rawdata_dendrcontour-files/fit_data .
   cp -r rawdata_dendrcontour-files/files_tipradius .
   cp -r rawdata_dendrcontour-files/list-files .
   cp -r rawdata_dendrcontour-files/list_yposmax.txt .
   rm -rf rawdata_dendrcontour-files
   #creation of gif of contour plots
   cd contour_plots
     for i in *
     do
     mv $i $(echo $i | awk -F '[_.]' '{printf "%09d" ,$3}')_contour.png
     done
   ls | convert -delay $delay -loop $numberloops *_contour.png animation_dendrite_$simname.gif
   mv animation_dendrite_$simname.gif ../
   cd ..
### evaluation of dendrite tip radius
   cd fit_data/fitparameter/
   mkdir fitparameter_sorted
   find . -type f -name "fit_parameters_final_*.txt" -exec cp {} ./fitparameter_sorted/ \;
   cd fitparameter_sorted
     for k in $(ls)
     do
     onlynumber=$(echo "$k"|awk -F'[_]' '{print $4}' | awk -F'[.]' '{printf "%06d\n", $1}')
     mv "$k" "$onlynumber"_gnuplotfitparameters.txt
     done
   find . -name "*_gnuplotfitparameters.txt" > listoffiles.txt
   cat listoffiles.txt | cut -c 3- | sort -n > listoffiles2.txt
   rm listoffiles.txtecho ""
   mv listoffiles2.txt listoffiles.txt
     for j in $(find . -name "*_gnuplotfitparameters.txt")
     do 
     echo $j | awk -F'[_]' '{print $1}' >> listtimesteps.txt
     done
   cat listtimesteps.txt | cut -c 3- | sort -n > listtimesteps2.txt
   rm listtimesteps.txt
   mv listtimesteps2.txt listtimesteps.txt
     for l in `cat listoffiles.txt`
     do
     cat $l | sed -n '1p' >> list_fitparameter_tip.txt
     done
   cat list_fitparameter_tip.txt | awk  '{print 1/(-1*2*$0)}' > list_radius_tip.txt
   paste listtimesteps.txt list_radius_tip.txt > list_combined_timesteps_radius.txt
   cd ../../../
   cp  fit_data/fitparameter/fitparameter_sorted/list_combined_timesteps_radius.txt .
   cd ..
   mkdir rawdata_output_OP
   find . -type f -name "DendrTipRadius_*.txt" -exec mv {} ./rawdata_output_OP/ \;
   find . -type f -name "DendrContour_*.txt" -exec mv {} ./rawdata_output_OP/ \;
   
   ####################################################
   elif [[ $video = n ]] ; then
   ####################################################
  
  stamptime=`date +%d-%m-%Y_%H:%M`
  cp -rv DendriteProperties DendriteProperties_"$stamptime"
  #>>>>>>>>>>>>>>>>>
  cd DendriteProperties_"$stamptime"
  mkdir results_$simname
  mkdir results_$simname/rawdata_dendrcontour-files
  find . -type f -name "DendrContour_*.txt" -exec cp {} ./results_$simname/rawdata_dendrcontour-files/ \;
  #>>>>>>>>>>>>>>>>>
  cd results_$simname/rawdata_dendrcontour-files/
   for b in $(ls)
   do  
   mv $b $(echo $b | awk -F '[_]' '{printf "%07d" ,$2}')_halbe_dendritenkontur.txt
   done
   #ls -I to ignore special pattern
   ls -I "liste*.txt" > liste_halbe_dendritenkontur.txt
   ##create complete dendrite contour
    for d in `cat liste_halbe_dendritenkontur.txt`
     do
      timestep=$(echo $d | awk -F '_' '{print $1}')
      contvalue2=$(cat $d | awk '{print $2*-1}')
      cat $d | awk '{print $1}' > liste_zcoord_half_"$timestep".txt
      cat $d | awk '{print $1}' > liste_zcoord_half2_"$timestep".txt
      cat $d | awk '{print $2}' > liste_contour_coord_pos_"$timestep".txt
      cat $d | awk '{print $2*-1}' > liste_contour_coord_neg_"$timestep".txt
      #add one point in -x direction to get correct tip geometry
      cat liste_contour_coord_neg_"$timestep".txt | sed -e '1s/0/-1/' > liste_contour_coord_neg2_"$timestep".txt
      mv liste_contour_coord_neg2_"$timestep".txt liste_contour_coord_neg_"$timestep".txt
      awk 'NF' liste_zcoord_half_"$timestep".txt liste_zcoord_half2_"$timestep".txt > z_coord_komp_"$timestep".txt
      awk 'NF' liste_contour_coord_pos_"$timestep".txt liste_contour_coord_neg_"$timestep".txt > contour_coord_komp_"$timestep".txt
    paste z_coord_komp_"$timestep".txt contour_coord_komp_"$timestep".txt > "$timestep"_komplette_dendriten_contour.txt
    ##finished complete dendrite contour
    mkdir list-files
    find . -type f -name "liste_zcoord_half_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "liste_zcoord_half2_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "liste_contour_coord_pos_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "liste_contour_coord_neg_*.txt" -exec mv {} ./list-files \; 
    find . -type f -name "liste_contour_coord_neg2_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "z_coord_komp_*.txt" -exec mv {} ./list-files \;
    find . -type f -name "contour_coord_komp_*.txt" -exec mv {} ./list-files \;
   done
  find . -type f -name "*_komplette_dendriten_contour.txt" > list_komplette_dendritenkontur.txt
  cat list_komplette_dendritenkontur.txt | cut -c 3- | sort -n > list_komplette_dendritenkontur2.txt
  mv list_komplette_dendritenkontur2.txt list_komplette_dendritenkontur.txt
   for y in `cat list_komplette_dendritenkontur.txt`
   do
   z=$(echo $y | awk -F'[_]' '{print $1}')
   yposmax=$(cat $y | awk 'NR==1{print $1}')
   yposmin=$(cat $y | awk 'NR==1{print $1-20}')
   echo "$yposmin" >> list_yposmax.txt
    gnuplot -persist <<- EOF
     reset
     set terminal pngcairo size 1000,1000 enhanced font "font Arial,22"
     set key tmargin
     set key font ",14"
     #set style line 5 lw 8 lc rgb '#F0F8FF' pt 4
     set style line 3 lw 5 lc rgb '#00ff00'
     set style line 4 lw 5 lc rgb '#00008B'
     set xlabel "Dendrite radius [{\/Symbol m}m]"
     set ylabel "Z-coordinate [-]"
     set xrange ["${xmin}":"${xmax}"]
     set yrange [0:"${ymax}"]
     #set arrow 1 from -"${intertip}","${yposmin}" to -"${intertip}","${yposmax}" nohead lw 4 lc rgb "#FF0000" dt "--"
     #set arrow 2 from  "${intertip}","${yposmin}" to  "${intertip}","${yposmax}" nohead lw 4 lc rgb "#FF0000" dt "--"
  #set object 11 circle at axis 0,"${yposcircle}" size "${intertip}"
  #set output "contour_timestep_${z}.png"
  set fit logfile "fit_data_full_contour_timestep_${z}.txt"
  f(x) = a*x**2+b*x+c
  fit f(x) "${y}" using 2:1 via a, b, c
  f2(x)= a2*x**2+b2*x+c2
  fit ["-${userdefxminfitpar}":"${userdefxminfitpar}"] f2(x) "${y}" using 2:1 via a2, b2, c2
  set print "fit_parameters_${z}.txt"
  print a
  print a2
  #plot "${y}" using 2:1 with points ps 5 title "timestep${z}", f(x) ls 3 title "parabolic fit trunk", f2(x) ls 4 title "parabolic fit tip", 1/0 lw 4 lc rgb "#FF0000" dt "--" title "fit range" 


EOF
done
   #creation of results subfolders and data migration
   mkdir fit_data
   mkdir fit_data/fitparameter
   mkdir contour_plots
   mkdir files_volle_dendritecontour
   mkdir files_halbe_dendritecontour
   find . -type f -name "fit_data_full_contour_timestep_*.txt" -exec mv {} ./fit_data/ \;
   find . -type f -name "contour_timestep_*.png" -exec mv {} ./contour_plots/ \;
   find . -type f -name "*_komplette_dendriten_contour.txt" -exec mv {} ./files_volle_dendritecontour/ \;
   find . -type f -name "*_halbe_dendritenkontur.txt" -exec mv {} ./files_halbe_dendritecontour/ \;
   find . -type f -name "fit_parameters_*" -exec mv {} ./fit_data/fitparameter/ \;
   find . -type f -name "DendrContour_*.txt" -exec mv {} ./rawdata_dendritencontour/ \;
   mkdir files_tipradius
   find . -type f -name "DendrTipRadius_*.txt" -exec mv {} ./files_tipradius/ \;
   find . -type f -name "list_komplette_dendritenkontur.txt" -exec mv {} ./list-files/ \;
   #<<<<<<<<<<<<<
   cd ..
   cp -r rawdata_dendrcontour-files/contour_plots .
   cp -r rawdata_dendrcontour-files/files_volle_dendritecontour .
   cp -r rawdata_dendrcontour-files/files_halbe_dendritecontour .
   cp -r rawdata_dendrcontour-files/rawdata_dendritencontour .
   cp -r rawdata_dendrcontour-files/fit_data .
   cp -r rawdata_dendrcontour-files/files_tipradius .
   cp -r rawdata_dendrcontour-files/list-files .
   cp -r rawdata_dendrcontour-files/list_yposmax.txt .
   rm -rf rawdata_dendrcontour-files
   #### evaluation of dendrite tip radius
   cd fit_data/fitparameter/
   mkdir fitparameter_sorted
   find . -type f -name "fit_parameters_final_*.txt" -exec cp {} ./fitparameter_sorted/ \;
   cd fitparameter_sorted
     for k in $(ls)
     do
     onlynumber=$(echo "$k"|awk -F'[_]' '{print $4}' | awk -F'[.]' '{printf "%06d\n", $1}')
     mv "$k" "$onlynumber"_gnuplotfitparameters.txt
     done
   find . -name "*_gnuplotfitparameters.txt" > listoffiles.txt
   cat listoffiles.txt | cut -c 3- | sort -n > listoffiles2.txt
   rm listoffiles.txtecho ""
   mv listoffiles2.txt listoffiles.txt
     for j in $(find . -name "*_gnuplotfitparameters.txt")
     do 
     echo $j | awk -F'[_]' '{print $1}' >> listtimesteps.txt
     done
   cat listtimesteps.txt | cut -c 3- | sort -n > listtimesteps2.txt
   rm listtimesteps.txt
   mv listtimesteps2.txt listtimesteps.txt
     for l in `cat listoffiles.txt`
     do
     cat $l | sed -n '1p' >> list_fitparameter_tip.txt
     done
   cat list_fitparameter_tip.txt | awk  '{print 1/(-1*2*$0)}' > list_radius_tip.txt
   paste listtimesteps.txt list_radius_tip.txt > list_combined_timesteps_radius.txt
   cd ../../../
   cp  fit_data/fitparameter/fitparameter_sorted/list_combined_timesteps_radius.txt .
   cd ..
   mkdir rawdata_output_OP
   find . -type f -name "DendrTipRadius_*.txt" -exec mv {} ./rawdata_output_OP/ \;
   find . -type f -name "DendrContour_*.txt" -exec mv {} ./rawdata_output_OP/ \;
   ######################################### 
   fi
   #########################################
   
############################################
fi
############################################