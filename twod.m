mov 3 #
		pl $1 $2
		do i=$2+1,$3,1 {\
			rd $i
			pla $1
		}
		#
pl	2	#
		#
		rd $2
		pla $1
pla	1	#
		#
		set i=1,$nr*$nh
		#
		set r=(int((i-0.5)/$nh) + 0.5)*$dr + $rin
                set th=(int(i-0.5 - $nh*int((i-0.5)/$nh)) + 0.5)*$dh
		#
		image($nr,$nh) $rin $rout 0 3.14159
		#
		set ir = int((r - $rin - 0.4*$dr)/$dr)
		set ih = int((th - 0.4*$dh)/$dh)
		#
		set image[ir,ih] = $1[i-1]
		#
		limits $rin $rout 0 3.14159
		erase
		box
		minmax min max echo $min $max
		if($min*$max < 0.) {\
			#define delta (($max-$min)/10.)
			#set lev=$min,-$delta,$delta
			#levels lev
			#ltype 2
			#contour
			#
			#set lev=$delta,$max,$delta
			#levels lev
			#ltype 0
			#contour
	                #
                        define delta ((-$min)/10.)
                        set lev=$min,-$delta,$delta
                        levels lev
                        ltype 2
                        contour
                        #
                        define delta ($max/10.)
                        set lev=$delta,$max,$delta
			#set lev=0,0.5,0.01
                        levels lev
                        ltype 0
                        contour
                        #
                        #ctype cyan
                        #set lev = 0
                        #levels lev
                        #contour
                        #ctype white
                        #
		} \
		else {\
			set lev=$min,$max,($max-$min)/10.
			#set lev=0,0.5,0.01
			levels lev
			ltype 0
			contour
		}
		#
rd	1	#
		if($1 < 10) {define num <00$1>} \
                else {if($1 < 100) {define num <0$1>} \
                else {define num <$1>}}
                echo $num
		rdp dump$num
rdpold	1	#
		#
		da dumps/$1
		lines 1 1
		read {_nr 1 _nh 2 _rin 3 _rout 4 _t 5}
		define nr (_nr)
		define nh (_nh)
		define rin (_rin)
		define rout (_rout)
		define dr ((_rout - _rin)/_nr)
		define dh (pi/$nh)
		lines 2 1000000
		#
		#
		read {r 1 h 2 rho 3 u 4 ur 5 uh 6 up 7 br 8 bh 9 bp 10}
		read {divb 11}
                read {b2 12}
		read {tm00 13 te00 14 tm03 15 te03 16 tm10 17 te10 18 tm13 19 te13 20 tm20 21 te20 22 tm23 23 te23 24}
		set ibeta=0.5*b2/((5/3-1)*u)
		#read {frd 11 fru 12 frur 13 fruh 14}
		set g = 5./3.
		#
		# auxiliary info.
		set p = (g - 1.)*u
		set K = p*rho**(-g)
		set cs = sqrt(g*p/rho)
		set ma = ur/cs
		#
		set pr = rho*ur
		set pp = rho*up
		set ph = rho*uh
		#
		set lrho = lg(rho)
		set luh = lg(abs(uh) + 1.e-20)
		set lu = lg(u)
		set ldivb = lg(abs(divb) + 1.e-20)
		#
rdp2	1	#
		#
		da dumps/$1
		lines 1 1000000
		#
		#
		read {rho 1 u 2 ur 3 uh 4 up 5 br 6 bh 7 bp 8 r 9 h 10}
                #
		gammienew
                #
rdp3	1	#
		#
		da dumps/$1
		lines 1 10000000
		#
		#
		read {x1 1 x2 2 r 3 h 4 rho 5 u 6 ur 7 uh 8 up 9 br 10 bh 11 bp 12 divb 13 bsq 14 tm00 15\
			  te00 16 tm03 17 te03 18 tm10 19 te10 20 tm13 21 te13 22 tm20 23 te20 24 tm23 25 te23 26}
                #
		gammienew3
                #
		#
gammieener    0 #
                da ener.out
		lines 1 10000000
                read {t 1 rmed 2 pp 3 e 4 pmid 5 pmed2 6 mdot 7 edot 8 mx1dot 9 mx2dot 10 ldot 11 bx1dot 12 bx2dot 13 bx3dot 14 mfladd 15 efladd 16 mx1fladd 17 mx2fladd 18 mx3fladd 19 bx1fladd 20 bx2fladd 21 bx3fladd 22 nstep 23 divbmax 24 divbavg 25}
		set realedot=edot-mdot
		#set realedot=edot
		# for use of jon's macro
		set min1d=mdot
		set ein1d=realedot
		set mx3in1d=ldot
		set td=t
		define GAMMIE (1)
                #
gammieenerold    0 #
                da ener.out
		lines 1 10000000
                read {t 1 rmed 2 pp 3 e 4 pmid 5 pmed2 6 mdot 7 edot 8 mx1dot 9 mx2dot 10 ldot 11 bx1dot 12 bx2dot 13 bx3dot 14}
		set realedot=edot-mdot
		# for use of jon's macro
		set min1d=mdot
		set ein1d=realedot
		set mx3in1d=ldot
		set td=t
		define GAMMIE (1)
                #
gammieenerold2    0 #
                da ener.out
		lines 1 10000000
                read {t 1 rmed 2 pp 3 e 4 pmid 5 pmed2 6 mdot 7 edot 8 ldot 9}
		set realedot=edot-mdot
                #
stressdoallg             #
		     ctype default
		     #rdnumd
		     #set start=2*$NUMDUMPS/3
		     #set end=$NUMDUMPS-1
		     set start=83
		     set end=125
		     avgtimeg 'dump' start end
		     # or just 1 dump
		     #rdp dump125.dat
		     #set Flmtime=-br*bp
		     #set entime=u
		     #set b2time=br**2+bh**2+bp**2
		     echo doing flm
		     thetaphiavg PI/6 Flmtime avgflm newx1
		     echo doing en
		     thetaphiavg PI/6 entime avgen newx1
		     echo doing b2
		     thetaphiavg PI/6 b2time avgb2 newx1
		     echo doing plot
		     define totalgrid (0)
		     stressplotgreal	
stressplotgreal #
		     erase
                     device postencap plotallgrmhdslimtimeavg.eps
		     set lgnewx1=LG(newx1)
		     if($totalgrid){\
		            set innerr = newx1[2]*.9999
		            set outerr = newx1[$nx-3]*1.0001
		         }
		     if($totalgrid==0){\
		            set innerr = newx1[0]*.9999
		            set outerr = newx1[$nx-1]*1.0001
		         }		         
		     set max=avgflm if((newx1>=innerr)&&(newx1<=outerr))
		     set press=avgen*($gam-1) if((newx1>=innerr)&&(newx1<=outerr))
		     set magpress=0.5*avgb2 if((newx1>=innerr)&&(newx1<=outerr))
		     set realx1=newx1 if((newx1>=innerr)&&(newx1<=outerr))		     
		     #limits lgnewx1 -8 -3
                     limits realx1 -8 -2
                     #ticksize -1 0 -1 0
		     ticksize 0 0 -1 0
		     fdraft
		     box
		     xla R c^2/GM
		     ltype 0
		     # maxwell stress
		     #connect lgnewx1 (LG(ABS(max))) 
		     connect realx1 (LG(ABS(max)))
		     # gas prssure
		     ltype 3
		     #connect lgnewx1 (LG(ABS(press)))
		     connect realx1 (LG(ABS(press)))
		     # magnetic pressure
		     ltype 4
		     connect realx1 (LG(ABS(magpress)))
		     device X11
		     #ctype red connect lgnewx1 (-2.45-2*lgnewx1)		
avgtimeg   3	# avgtimeg (e.g. avgtimeg 'dump' start end)
		#
                set h1=$1
		set h3='.dat'
		#
		set rtime=1,$nx*$ny*$nz,1
		set entime=1,$nx*$ny*$nz,1
                set cs2time=1,$nx*$ny*$nz,1
		set b2time=1,$nx*$ny*$nz,1
                set vtimex=1,$nx*$ny*$nz,1
                set vtimey=1,$nx*$ny*$nz,1
                set vtimez=1,$nx*$ny*$nz,1
		set btimex=1,$nx*$ny*$nz,1
		set btimey=1,$nx*$ny*$nz,1
		set btimez=1,$nx*$ny*$nz,1
		set Fltime=1,$nx*$ny*$nz,1
		set Flrtime=1,$nx*$ny*$nz,1
		set Flmtime=1,$nx*$ny*$nz,1
		set Mftimex=1,$nx*$ny*$nz,1
		set Mftimey=1,$nx*$ny*$nz,1
		set Lftimex=1,$nx*$ny*$nz,1
		set Lftimey=1,$nx*$ny*$nz,1
		set Fvzdxtime=1,$nx*$ny*$nz,1
		set Fvzdytime=1,$nx*$ny*$nz,1            
		set Fmdxtime=1,$nx*$ny*$nz,1
		set Fmdytime=1,$nx*$ny*$nz,1		
                #
                set rtime=rtime*0
                set entime=entime*0
                set cs2time=cs2time*0
                set b2time=b2time*0
                set vtimex=vtimex*0
                set vtimey=vtimey*0
                set vtimez=vtimez*0
		set Fltime=Fltime*0
		set Flrtime=Flrtime*0
		set Flmtime=Flmtime*0
		set Mftimex=Mftimex*0
		set Mftimey=Mftimey*0
		set Lftimex=Lftimex*0
		set Lftimey=Lftimey*0
		set btimex=btimex*0
		set btimey=btimey*0
		set btimez=btimez*0
		set Fvzdxtime=Fvzdxtime*0
		set Fvzdytime=Fvzdytime*0
		set Fmdxtime=Fmdxtime*0
		set Fmdytime=Fmdytime*0		
                #
		set numstart=$2
		set numend=$3
                set numtotal=numend-numstart+1
                do ii=numstart,numend,1 {
                  set h2=sprintf('%03d',$ii)
                  set _fname=h1+h2+h3
                  define filename (_fname)
		  rdp $filename
		  set rtime=rtime+rho
                  set entime=entime+u
                  set cs2time=cs2time+cs*cs
		  set b2time=b2time+br*br+bh*bh+bp*bp
                  set vtimex=vtimex+ur
                  set vtimey=vtimey+uh
                  set vtimez=vtimez+up
		  set btimex=btimex+br
		  set btimey=btimey+bh
		  set btimez=btimez+bp
		  set Fmdx=rho*ur
		  set Fmdy=rho*uh
		  set Fmdz=rho*up
		  set Fmx=Fmdx*x12*x12*sin(x22)
		  set Fmy=Fmdy*x12*sin(x22)
		  set Fmz=Fmdz*x12
		  set Fvzdx=x12*sin(x22)*(rho*ur*up)
		  set Fvzdy=x12*sin(x22)*(rho*uh*up)
		  set Fvzdz=x12*sin(x22)*(rho*up*up)
		  set Fvzx=Fvzdx*x12*x12*sin(x22)
		  set Fvzy=Fvzdy*x12*sin(x22)
		  set Fvzz=Fvzdz*x12
		  set Mftimex=Mftimex+Fmx
		  set Mftimey=Mftimey+Fmy
		  set Lftimex=Lftimex+Fvzx
		  set Lftimey=Lftimey+Fvzy
		  set Flrtime=Flrtime+(r*ur*up)
		  set Fltime=Fltime+(r*ur*up-bp*br)
		  set Flmtime=Flmtime+(-bp*br)
		  set Fvzdxtime=Fvzdxtime+Fvzdx
		  set Fvzdytime=Fvzdytime+Fvzdy
		  set Fmdxtime=Fmdxtime+Fmdx
		  set Fmdytime=Fmdytime+Fmdy
		}
                set rtime=rtime/numtotal
                set entime=entime/numtotal
                set cs2time=cs2time/numtotal
		set b2time=b2time/numtotal
                set vtimex=vtimex/numtotal
                set vtimey=vtimey/numtotal
                set vtimez=vtimez/numtotal
		set btimex=btimex/numtotal
		set btimey=btimey/numtotal
		set btimez=btimez/numtotal
		set Fltime=Fltime/numtotal
		set Flrtime=Flrtime/numtotal
		set Flmtime=Flmtime/numtotal
		set Mftimex=Mftimex/numtotal
		set Mftimey=Mftimey/numtotal
		set Lftimex=Lftimex/numtotal
		set Lftimey=Lftimey/numtotal
		set Fvzdxtime=Fvzdxtime/numtotal
		set Fvzdytime=Fvzdytime/numtotal
		set Fmdxtime=Fmdxtime/numtotal
		set Fmdytime=Fmdytime/numtotal
                #
gammie1  #
	 define interp (0)
	 define PLANE (3)
	 define nx 128
	 define ny 128
         define nz 1
         define Lx 1
         define Ly 1         
	 define Lz 1
         define Sx 0
         define Sy 0         
	 define Sz 0
         define dx ($Lx/($nx))
         define dy ($Ly/($ny))
         define dz ($Lz/($nz))
         define ncpux1 1
         define ncpux2 1
         define ncpux3 1
gammie12  #
         # used for reversed ordering of radius and theta
         set x12=h # my x1 is his theta
         set x22=r # my x2 is his radius
	 define nx 256
	 define ny 256
         define nz 1
	 define Sx (x12[0])
	 define Sy (x22[0])
         define Sz (1)
	 define Lx (x12[$nx-1])
	 define Ly (x22[$ny-1])
	 define Lz (1)
         define dx ($Lx/($nx))
         define dy ($Ly/($ny))
         define dz ($Lz/($nz))
         define ncpux1 1
         define ncpux2 1
         define ncpux3 1
	 define interp (0)
	 define coord (1) # don't want any special treatement since uniform grid
	 define x1label "\theta"
	 define x2label "radius"
	 set k=1,$nx*$ny,1
	 set k=0*k
	 #
gammienew  #
	 define gam (5.0/3.0)
         # used for normal ordering of radius and theta
	 # set whether boundary zones
	 define bc (0) # 2*2 or 0
	 # x1
	 set x12=r
         set x1=r
         set temptemp=x1[1]-x1[0]
	 set dx1=h*0+temptemp
	 set dx12=dx1
	 # x2
         set x22=h
         set x2=h
	 set temptemp=x2[$nx]-x2[0]
	 set dx2=h*0+temptemp
	 set dx22=dx2
         # x3
	 set x3=r*0
	 set x32=x3		
         set dx3=r*0
         set dx3=dx3+2*PI
	 set dx32=dx3		
	 define nx (128+$bc)
	 define ny (128+$bc)
         define nz 1
	 define Sx (x12[0])
	 define Sy (x22[$nx-1])
         define Sz (1)
	 define Lx (x12[$nx-1])
	 define Ly (x22[$ny*$nx-1])
	 define Lz (1)
         define dx ($Lx/($nx))
         define dy ($Ly/($ny))
         define dz ($Lz/($nz))
         define ncpux1 1
         define ncpux2 1
         define ncpux3 1
	 define interp (0)
	 define coord (1) # don't want any special treatement since uniform grid
	 define x1label "radius"
	 define x2label "\theta"
	 set k=1,$nx*$ny,1
	 set i=1,$nx*$ny,1
	 set j=1,$nx*$ny,1
	 set i=k%$nx
	 set j=INT(k/$nx)
	 set k=0*k
	 
	 #
gammienew3  #
	 define gam (5.0/3.0)
         # used for normal ordering of radius and theta
	 # set whether boundary zones
	 define bc (0) # 2*2 or 0
	 # x1
	 set x12=x1
	 #set x1=x1
         set temptemp=x1[1]-x1[0]
	 set dx1=h*0+temptemp
	 set dx12=dx1
	 # x2
         set x22=x2
	 #set x2=h
	 set temptemp=x2[$nx]-x2[0]
	 set dx2=h*0+temptemp
	 set dx22=dx2
         # x3
	 set x3=r*0
	 set x32=x3		
         set dx3=r*0
         set dx3=dx3+2*PI
	 set dx32=dx3		
	 define nx (64+$bc)
	 define ny (64+$bc)
         define nz 1
	 define Sx (x12[0])
	 define Sy (x22[$nx-1])
         define Sz (1)
	 define Lx (x12[$nx-1])
	 define Ly (x22[$ny*$nx-1])
	 define Lz (1)
         define dx ($Lx/($nx))
         define dy ($Ly/($ny))
         define dz ($Lz/($nz))
         define ncpux1 1
         define ncpux2 1
         define ncpux3 1
	 define interp (0)
	 define coord (1) # don't want any special treatement since uniform grid
	 define x1label "radius"
	 define x2label "\theta"
	 set k=1,$nx*$ny,1
	 set i=1,$nx*$ny,1
	 set j=1,$nx*$ny,1
	 set i=k%$nx
	 set j=INT(k/$nx)
	 set k=0*k
	 
	 #
gammie2  1 #
	 da $1
         read{which 1 cpu 2 i 3 j 4 pr 5 val 6}
         #set x12=-2*$dx,$Lx+2*$dx,$dx
         #set x22=-2*$dy,$Ly+2*$dy,$dy
	 set x12=j if((pr==0)&&(which==0)&&(cpu==0))
	 set x22=i if((pr==0)&&(which==0)&&(cpu==0))
	 define interp (0)
	 define coord (1)
	 define Sx (x12[0])
	 define Sy (x22[0])
	 #define Sz (x32[0])      # 3D
	 define Sz (x3[0])      # 2D
	#
localdiff #
          set vpos=1,4,1
          set vpos=vpos*0
          set maxdiff=rho*0
         do ii=0,$nx*$ny-1,1 {
           set iindex=INT($ii%($nx))
           set jindex=INT($ii/($nx))
           if((iindex>0)&&(iindex<($nx-1))&&(jindex>0)&&(jindex<($ny-1))){\
             set v0=rho[$ii]
             set vpos[0]=rho[$ii+1]
             set vpos[1]=rho[$ii-1]
             set vpos[2]=rho[$ii+$nx]
             set vpos[3]=rho[$ii-$nx]
             #set olddiff=ABS(v0-vpos[0])
             set olddiff=ABS(vpos[0]/v0)
             do jj=1,3,1 {
              #set newdiff=ABS(v0-vpos[$jj])
              set newdiff=ABS(vpos[$jj]/v0)
              if(newdiff>olddiff){ set olddiff=newdiff}
             }
             set maxdiff[$ii]=olddiff
           }\
           else{\
             #set maxdiff[$ii]=0
             set maxdiff[$ii]=1
           }
         }
