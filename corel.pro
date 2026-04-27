; pixel scale by cross-correlation to determine shift between rings
pro corel, file1


img1 = float(readfits(file1))

 nx = (size(img1))[1] 

img1 -= median(img1)  ; remove background

   ft1 = fft(img1)

   ccf = float(shift(fft(ft1*conj(ft1), /inverse),nx/2,nx/2))
   ccf -= median(ccf) ; subtract median level

  m = 256 < nx/2  ; half-size of the central zone
;  ccf1 = ccf[nx/2-m:nx/2+m-1,nx/2:nx/2+m-1] ; upper half
  ccf1 = ccf[nx/2-m:nx/2+m-1,nx/2-m:nx/2+m-1] ; center


; Arrays of X and Y-coordinates   
   x = (findgen(2*m) - m) # replicate(1.,2*m) 
   y = transpose(x)

  r = shift(dist(2*m,2*m), m, m)
; r = sqrt(x^2 + y^2) : alternative method

  rmask = 30. ; mask central zone 2*rring 
  ccf1 [where(r lt rmask)] = 0. 
   
   tvscl, ccf1  ; display

; Optionally, save the fits file of the correlation
;  writefits, file1+'_corel.fits', ccf


   cmax = max(ccf1)
   peak = where((ccf1 gt 0.7*cmax) and (y gt 0))  ; pixels near peak, upper half 

   tmp = ccf1[peak]
   xc = total(x[peak]*tmp)/total(tmp)
   yc = total(y[peak]*tmp)/total(tmp)
   d = sqrt(xc^2 + yc^2)
   print, 'Peak distance [pix]: ', d

 stop

   end
