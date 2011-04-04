/*
 * apety_utils.i
 *
 * This file is part of the APETY package, the Ao Psf Estimation Tool in Yorick
 *
 *
 * This program is free software; you can redistribute it and/or  modify it
 * under the terms of the GNU General Public License  as  published  by the
 * Free Software Foundation; either version 2 of the License,  or  (at your
 * option) any later version.
 *
 * This program is distributed in the hope  that  it  will  be  useful, but
 * WITHOUT  ANY   WARRANTY;   without   even   the   implied   warranty  of
 * MERCHANTABILITY or  FITNESS  FOR  A  PARTICULAR  PURPOSE.   See  the GNU
 * General Public License for more details (to receive a  copy  of  the GNU
 * General Public License, write to the Free Software Foundation, Inc., 675
 * Mass Ave, Cambridge, MA 02139, USA).
 *
*/

require,"iterkolmo.i";

func circavg(a,center=,middle=)
/* DOCUMENT circavg
 *  
 * average=circavg(array[,center=,middle=])
 *
 * This routine returns the circular mean of an array. The routine dist is
 * used to compute an array of coordinate of points to be averaged.
 *
 * KEYWORDS :
 * a      (input) : The array of which we want the circular mean. It can be
 *                  a long, a float or a double. Complex are not supported
 * center (input) : An array of the form [x,y] which give the position of
 *                  the origin for the circular mean calculation
 * middle (input) : A flag to indicate that the origin for the circular mean
 *                  calculation is the middle of the array a
 *
 * SEE ALSO: dist, circavg2
 */ 
{
  s=dimsof(a);

  if (!is_set(middle)) middle=0;
  if (s(1) != 2) write,"error - invalid dimensions";
  if (s(3) != s(2)) write,"error - invalid dimensions";

  dim=s(2);
  
  if (center!=[]) {
    s=dimsof(center);
    if ((s(1) != 1) | (s(1) != 2)) \
       write,"error - center has invalid dimensions";

    center=long(center);

    if (middle) {
      center=long([0,0]);
      write,"error - center and middle are not compatible keywords";
    }
  } else { 
    if (middle) center=long([0,0]);
    else center=[dim/2,dim/2];
  }
  
  r=long(roll(long(dist(dim)+.5)+1,[center(1),center(2)]));
  j=long(max(r));
  n=array(long,j);
  sx=array(double,j);
  dim2=long(dim)*long(dim);
  
  for (i=1;i<=dim2;i++) {
    j=r(i);
    sx(j)=sx(j)+a(i);
    n(j)=n(j)+1;
  }
  
  return sx/n;
}


func calc_uij(mode1,mode2,pup,den)
/* DOCUMENT calc_uij
 *  This routine computes the Uij functions using the standard fourier domain formulae
 * KEYWORDS :
 * mode1 (input) : The first mode (Mi)
 * mode2 (input) : The second mode (Mj)
 * pup   (input) : the pupil
 */ 
{
  npix = dimsof(mode1)(2);
  mi = mj = p = uij = array(float,[2,2*npix,2*npix]);
  mi(1:npix,1:npix) = mode1;
  mj(1:npix,1:npix) = mode2;
  p(1:npix,1:npix)  = pup;

  //den  = (fft(fft(p,1)*conj(fft(p,1)),-1)).re;
  mask = den > max(den)*1.e-7;
  pix  = where(mask);
  //uij(pix) = (fft(2.*(fft(mi*mj*p,1)*conj(fft(p,1))-fft(mi*p,1)*fft(mj*p,1)).re,-1)(pix)).re/den(pix);
  uij(pix) = ((-2*fft((fft(mi*p,1)*conj(fft(mj*p,1))).re,-1)+2*fft((fft(mi*mj*p,1)*conj(fft(p,1))).re,-1)).re)(pix)/den(pix);
                                                                                                        
  return uij;   
}

func calc_uijf(ftMode1,ftMode2,mode1,mode2,den,conjftpup)
/* DOCUMENT calc_uijf
 *  This routine computes the Uij functions using pre-computed inputs (fft of mode 1 and conjugate fft of mode 2)
 * this function is faster but implies some pre-computations 
 */ 
{
  mask = den > max(den)*1.e-7;
  pix  = where(mask);
  uij = mode1*0.;
  uij(pix) = ((-2.*fft((ftMode1*conj(ftMode2)).re,-1)+2.*fft((fft(mode1*mode2,1)*conjftpup).re,-1)).re)(pix)/den(pix);
  return uij;
}

func calc_dphi(phase,pup,den)
/* DOCUMENT calc_uij
 *  This routine computes the Uij functions using the standard fourier domain formulae
 * KEYWORDS :
 * mode 1 (input) : The first mode (Mi)
 * mode 2 (input) : The second mode (Mj)
 * pup    (input) : the pupil
 */ 
{
  npix = dimsof(phase)(2);
  mi = p = dphi = array(float,[2,2*npix,2*npix]);
  mi(1:npix,1:npix) = phase;
  p(1:npix,1:npix)  = pup;

  //den  = (fft(fft(p,1)*conj(fft(p,1)),-1)).re;
  mask = den > max(den)*1.e-7;
  pix  = where(mask);
  dphi(pix) = fft(2*((fft(mi^2*p,1)*conj(fft(p,1))).re - abs(fft(mi*p,1))^2),-1).re(pix)/den(pix);

  return dphi;   
}

func calc_dphif(phase,pup,den,conjftpup)
/* DOCUMENT calc_uijf
 *  This routine computes the Uij functions using the standard fourier domain formulae
 * KEYWORDS :
 */ 
{
  npix = dimsof(phase)(2);
  mi = p = dphi = array(float,[2,2*npix,2*npix]);
  mi(1:npix,1:npix) = phase;
  p(1:npix,1:npix)  = pup;

  //den  = (fft(fft(p,1)*conj(fft(p,1)),-1)).re;
  mask = den > max(den)*1.e-7;
  pix  = where(mask);
  dphi(pix) = fft(2*((fft(mi^2*p,1)*conjftpup).re - abs(fft(mi*p,1))^2),-1).re(pix)/den(pix);
  
  return dphi;   
}

/////////////////////////////////////////////////////////////////////////////////////////////
//////////////////                                       ////////////////////////////////////
///////////////// THE FOLLOWING ROUTINES ARE NOT WORKING ////////////////////////////////////
////////////////                                         ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


func get_aliasing(dr0at05mic,zenithangle,nit)
/* DOCUMENT get_aliasing
 *  This routine requires an initialized yao environment to run
 *  This routine computes nit aliasing measurements for a specified r0 and zenith
 *  this routine is for test purposes only ... the aliasing component should be computed
 *  once for all and re-scaled using the proper r0
 *
 * KEYWORDS :
 *  dr0at05mic  (input) : the D/r0 @ 0.5 mic
 *  zenithangle (input) : the zenith angle
 *  nit         (input) : the # of iterations (for statistics)
 */ 
{
  // some init
  mInf   = *dm(1)._def;
  sz     = dimsof(mInf);
  size   = sim._size;
  nModes = sz(4);
  pupd   = sim.pupildiam;
  n1     = dm(1)._n1;
  n2     = dm(1)._n2;
  pup    = ipupil(n1:n2,n1:n2);
  pix    = where(ipupil);
  pixb   = where(pup);

  matPass = (LUsolve(mInf(*,)(+,)*mInf(*,)(+,))(,+)*mInf(*,)(,+))(,pixb);

  smes = array(float,[2,wfs(1)._nsub,nit]);
  nici_p = array(float,[2,size,size]);
  
  for (i=1;i<=nit/2;i++) {
    // generate phase screen
    tmp = CreatePhaseScreens_silent(size,size);
    // this fucntion generates two phase screens each time so the following loop is splittd in 2
    write,format="\rPhase screen %d out of %d",i,nit/2;
    
    phase = tmp(,,1)*weight;
    // remove piston
    phase *= ipupil;
    phase -= phase(where(ipupil))(avg);
    // project the phase screen onto the modes (using matPass)
    vect  = phase(pix)(+)*matPass(,+);
    // compute the miror component of the phase using this vector and the influence functions
    nici_p(n1:n2,n1:n2) = vect(+)*mInf(,,+);
    // compute the measurement on the orthogonal component of the phase (i.e. phase - mirror component)
    smes(,2*i-1) = CurvWfs(pupil,float(phase-nici_p),1)-*wfs(1)._refmes;

    // do the same thing with the second phase screen
    phase = tmp(,,2)*weight;
    phase *= ipupil;
    phase -= phase(where(ipupil))(avg);
    vect  = phase(pix)(+)*matPass(,+);
    nici_p(n1:n2,n1:n2) = vect(+)*mInf(,,+);
    smes(,2*i) = CurvWfs(pupil,float(phase-nici_p),1)-*wfs(1)._refmes;
  }

  return smes;
}

func get_ortho(r0,nit)
/* DOCUMENT get_ortho
 *  This routine computes the orthogonal phase structure function i.e. DphiOrtho
 *  this routine is for test purposes only ... the orthogonal component should be computed
 *  once for all and re-scaled using the proper r0
 *  this function is very similar to get_aliasing (see comments)
 *  This routine requires an initialized yao environment to run
 *
 * KEYWORDS :
 */ 
{
  size   = sim._size;
  
  nmodes = dm._nact(sum);
  mInf = array(0.0f,size,size,nmodes);
  cpt = 0;
  for (nm=1;nm<=ndm;nm++) {
    n1 = dm(nm)._n1; n2 = dm(nm)._n2;
    for (cc=1;cc<=dm._nact(nm);cc++) {
      cpt ++;
      scommand = float(modToAct(indexDm(1,nm):indexDm(2,nm),cc));
      mInf(n1:n2,n1:n2,cpt) = comp_dm_shape(nm,&scommand);
    }
  }
  
  n1     = dm(1)._n1;
  n2     = dm(1)._n2;
  sz = n2-n1+1;
  pupd = sim.pupildiam;
  ind0 = (sz - pupd)/2;
  n1 += ind0;
  n2 -= ind0;

  mInf = mInf(n1:n2,n1:n2,);
  mypup = ipupil(n1:n2,n1:n2);
  pvalid = where(mypup);  
  // computing the projection matrix from phase screen to modes i.e. phase2modes
  matPass = (LUsolve(mInf(*,)(+,)*mInf(*,)(+,))(,+)*mInf(*,)(,+))(,pvalid);
  
  // some init
  AB, pupd, A, B, istencil;
  phase = array(0.0,pupd,pupd);
  phase = extrude(phase,r0,A,B,istencil);
  for (i=1;i<=pupd;i++) phase = extrude(phase,r0,A,B,istencil);

  dphiOrtho = 0.;
  
  for (i=1;i<=nit;i++) {
    phase = extrude(phase,r0,A,B,istencil);
    tmp      = float(phase-phase(pvalid)(avg));
    vect  = tmp(pvalid)(+)*matPass(,+);
    
    dm_p = vect(+)*mInf(,,+); 
    tmp_ortho = float(tmp - dm_p);
    
    dphiOrtho += calc_dphif(tmp_ortho,mypup,den,conjftpup);
    write,format="\rPhase screen %d out of %d",i,nit;
   }
  
  dphiOrtho /= nit;

  return dphiOrtho;
}


func generateVKspectrum_silent(dim,k0,nalias=)
/* DOCUMENT func generateVKspectrum(sdim,bdim,k0)
   generate correct VoKarman spectrum including aliasing.
     
   SEE ALSO:
 */
{
  if (is_void(nalias)) nalias = 0;
  
  res = array(float,[2,dim,dim]);
  for (i=-nalias;i<=nalias;i++) {
    for (j=-nalias;j<=nalias;j++) {
      tmp = sqrt(dist(dim,xc=i*dim+dim/2,yc=j*dim+dim/2)^2.f+k0^2.);
      if ((i==0) && (j==0)) tmp = clip(tmp,1.,);
      amp = (6.88*0.00969)*tmp^(-11.f/6.f);
      res += amp;
      //print,i,j;tv,res; pause,500;
    }
  }
  roll,res;
  res = float(res);
  return res;
}

func generatePhaseWithL0_silent(dim,l0,nalias=)
  /* DOCUMENT generate_phase(size,l0)
     Generate by Fourier an un-normalized 2D phase screen from the 
     -11/3 amplitude and a randomn phase component. Returns the real and 
     complex parts.
     Uses fftVE and cosf/sinf to keep floats for RAM use consideration
     (the previous version of this routine was using the yorick fft,
     thus double complex, which limits things on my machine to 4096 screens).

     dim: desired dimension of the result phase screens
     l0: outer scale IN PIXELS
     
     F.Rigaut, 2001/11/10.
     SEE ALSO: CreatePhaseScreens, PhaseStructFunc.
  */

{
  if (l0 == 0.) { k0=0.0f; } else { k0 = float(dim)/l0; }
  randomize;
  //  tmp = clip(float(sqrt(dist(dim)^2.f+k0^2.)),1e-8,);
  //  amp = 6.88*0.00969*dim*eclat(tmp^(-11.f/6.f));
  //  amp = eclat(generateVKspectrum(dim,clip(2*dim,,2048),k0));
  //  amp = dim*eclat(generateVKspectrum(dim,1024,k0));
  amp = dim*generateVKspectrum_silent(dim,k0,nalias=nalias);
  amp = float(amp);
  // normalized so that the structure function is correct for r0 = 1 pixel
  tmp = [];
  pha = reform(float(2*pi)*gaussdev(dim*dim),dimsof(amp));
  re  = amp*cosf(pha);
  im  = amp*sinf(pha);
  amp = pha = [];
  re(1,1) = 0.0f;
  im(1,1) = 0.0f;
  phaout = fftVE(re,im,1);
  re = im = [];
  return phaout;
}

func CreatePhaseScreens_silent(dimx,dimy,l0=,prefix=,nalias=)
  /* DOCUMENT CreatePhaseScreens(dimx,dimy,prefix=)
     Create phase screens and save them in fits files.
     The saved phase screens have a dimension dimx*dimy.
     Number of phase screens = dimx/dimy.
     The phase screens are normalized so that 1 pixel = 1 r0,
     i.e. the variance of the squared difference of the screen
     with itself at one pixel interval is 6.88 (rd^2).

     dimx = long dimension of result screens
     dimy = short dimension of result screens
     prefix = Prefix to filename. if prefix is not set, the
       screens are returned by not saved.

     Example:
     CreatePhaseScreens,2048,256,prefix="screen256"
     
     F.Rigaut, 2001/11/10.
     modify 2003 Feb 24 to add dimy (before dimy=256) and prefix
     SEE ALSO: generate_phase, PhaseStructFunc.
  */ 

{
  if (is_void(l0)) l0 = 0.;
  nscreen = dimx/dimy*2;
  pscreen = generatePhaseWithL0_silent(dimx,l0,nalias=nalias);

  off = [1,5]; // spatial offset for structure function normalization
  psfunc = array(float,off(2));
  for (i=off(1);i<=off(2);i++ ) {
    fsx     = avg((pscreen(1:-i,,)-pscreen(i+1:,,))^2.);
    fsy     = avg((pscreen(,1:-i,)-pscreen(,i+1:,))^2.);
    psfunc(i)      = sqrt((fsx+fsy)/2.);
  }
  theo = sqrt(6.88*indgen(off(2))^1.66);
  
  nfact=avg(psfunc(off(1):off(2))/theo(off(1):off(2)));
  pscreen = pscreen/nfact;

  pscreen = reform(pscreen,[3,dimx,dimy,nscreen]);

  return pscreen;
}
