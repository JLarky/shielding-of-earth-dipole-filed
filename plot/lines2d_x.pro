;  this is an analog of MONLINES - plots the lines in the noon-midnight plane

function xticks,axis,index,value
return, ''
end
function yticks,axis,index,value
return, ''
end
;   The above functions are for producing axes without ticks/labels




function xticks,axis,index,value
return, ''
end
function yticks,axis,index,value
return, ''
end

linname=''
numname=''
magname=''


print , '   Standard names ?  (lines.dat, numpnt.dat, magnpaus.dat)'
print , '      1 - yes,  0 - no'
;read, stand
stand=1            ;  temporarily

if stand eq 1 then begin

linname='lines.dat'
numname='numpnt.dat'
magname='magnpaus.dat'


goto, conti
endif

print, '   Enter the line coordinate filename'
read, linname
print, '   Enter the line point enumeration filename'
read, numname
print, '   Enter the magnetopause coordinate filename '
print , '    (enter none if not needed or stand if standard one)'
read, magname

if magname eq 'stand' then magname='mgnplong.dat'
conti:
set_plot, 'x'


Perpdist=30
Sundist=20

print,'  Enter tailward tracing limit (|X_min| in Re) '
;read, Taildist
Taildist=30

xticksno=Taildist/10+2
yticksno=6

if Taildist lt 20. then begin
Perpdist=15
Sundist=15
xticksno=6
yticksno=6
endif

;Taildist=30.                      ;  temporarily


print, '  Enter line thickness (2 is good, but in some cases 1 is better)'
;read, thickline
thickline=1                 ;  temporarily


nummag=fltarr(100)

openr,9,numname
readf,9,nummag0
print, ' Number of lines: ',nummag0
nummag=fltarr(nummag0)
readf,9,nummag
close,9

ipoints=0
for i=0,nummag0-1 do begin
ipoints=ipoints+nummag(i)
; print, '  number of pts:',nummag(i)
endfor

print, ' Number of points: ',ipoints

linmag=fltarr(3,ipoints)


openr,10,linname

for i=0,ipoints-1 do begin
readf,10, aa,bb,cc
linmag(0,i)=aa
linmag(1,i)=bb
linmag(2,i)=cc
; print, linmag(*,i)
endfor
close,10

;openr,9,linname
;openr,10,numname
;readf,9,linmag
;readf,10,nummag



print,'  Drawing device:  enter 0 for screen and 1 for PostScript'
;read,inp
inp=1                               ;  temporarily


tickleng=0.005
print,'  Draw a grid ?   1 - yes, 0 - no'
;read, igri
igri=1                                ;  temporarily

if igri eq 1 then tickleng=1

if inp eq 1 then begin

print , '   Enter scale factor for the entire plot ( ~0.3 - ~1)'
;read ,  scl
scl=1

endif

if inp eq 1 then begin
SET_PLOT, 'PS'
device, filename='lines2d.ps',/portrait,xsize=16.,xoffset=5,$
ysize=16.*(2.*Perpdist/(Sundist+Taildist)), yoffset=8, /times, scale_factor=scl
endif


icount=0


for i=1,nummag0 do begin
if i eq 1 then plot,linmag(0,icount:icount+nummag(i-1)-1),$
linmag(2,icount:icount+nummag(i-1)-1),XRANGE=[Sundist,-Taildist],YRANGE=[-Perpdist,Perpdist],$
thick=thickline,BACKGROUND=255,color=0,font=0,xthick=5,ythick=5,$
xtitle='XGSM, R!DE',ytitle='ZGSM, R!DE',xstyle=1,ystyle=1,/data,$
xcharsize=1.,ycharsize=1.,xticks=xticksno,yticks=yticksno

IF i NE 1 THEN oplot,linmag(0,icount:icount+nummag(i-1)-1),$
linmag(2,icount:icount+nummag(i-1)-1),thick=thickline,COLOR=0

icount=icount+nummag(i-1)
endfor

if magname ne 'none' then magn=fltarr(2,800)  ;  this is for  the magnetopause


if magname ne 'none' then openr,8,magname

if magname ne 'none' then begin
 icoun=-1
 while (not eof(8)) do begin
 icoun=icoun+1
 readf,8,xxm,zzm
 magn(0,icoun)=xxm
 magn(1,icoun)=zzm
 endwhile
endif
if magname ne 'none' then close,8


magnp=fltarr(2,icoun)
for ic=0,icoun-1 do begin
magnp(0,ic)=magn(0,ic)
magnp(1,ic)=magn(1,ic)
endfor


axis,xax=0,/data,xthick=0.3,COLOR=0,xticks=xticksno,ticklen=tickleng,xtickformat='xticks',xgridstyle=1
axis,yax=0,/data,ythick=0.3,COLOR=0,yticks=yticksno,ticklen=tickleng,ytickformat='yticks',ygridstyle=1

if magname ne 'none' then oplot,magnp(0,*),magnp(1,*),linestyle=5,thick=4,COLOR=0
;if magname ne 'none' then oplot,magnp(0,*),-magnp(1,*),linestyle=5,thick=4,COLOR=0


; Now draw Earth's contour:

   a=findgen(33)*(!pi*2/32)
   oplot,0.92*cos(a),0.92*sin(a),color=0,thick=2

;  Fill Earth's interior with white/black dashing

  a=findgen(181)*(!pi/180)    ;   dash the dayside by white
  c=0.89*sin(a)
  s=0.89*cos(a)
  polyfill, c,s,color=255

  c=-0.89*sin(a)            ;   dash the nightside by black
  s=0.89*cos(a)
  polyfill, c,s,color=10




if inp eq 1 then device,/close

END




