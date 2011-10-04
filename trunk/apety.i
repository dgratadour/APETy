/*
 * apety.i
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

require,"yao.i";
require,"yao_funcs.i";
require,"apety_utils.i";
YAO_SAVEPATH = "/tmp/";

func get_psf(cbmes, dphiOrtho, covmes_alias, den, gamma_eps=, cov_basis=, noise_method=, disp=)
/* DOCUMENT get_psf(cbmes,noisevar,dphiOrtho,cov_mes_alias,den,disp=)
 
 cbmes = circular buffer mesure
 den # OTF_tel # correlate(pupil, pupil)
 
 */ 
{

	extern cMat; // command matrix => mes2volt
	// RECONSTRUCTEUR
    
	if (disp == []) disp = 0;
	
	// noise calculation
	if (is_void(noise_method)) {
		noise_method = "DSP";
		noisevar     = get_noisevar(cbmes, noise_method, cbcom, iMat);
	} else {
		noisevar     = get_noisevar(cbmes, noise_method, cbcom, iMat);
	}

	// erreur de sous-modelisation
	under_modeling_var = 0.021 * (atm.dr0at05mic)^(5./3.);

	
	// compute the (noisy) measurements covariance : Cww
	covmes = (cbmes(,+)*cbmes(,+)) / dimsof(cbmes)(3)^2;
	// remove the noise variance on the diagonal
	covNoise = diag(noisevar);
	covmes  -= covNoise;
	
	// compute the modes covariance from the measurements covariance (x comMat) eq. 41
	cov_eps =(cMat(,+)*covmes(+,))(,+) * cMat(,+); // ~ D+(Cww)T(D+) mais pour cov des modes
	// compute the modes covariance from the aliasing measurements component Crr of eq. 44
	cov_alias =(cMat(,+)*(covmes_alias)(+,))(,+) * cMat(,+);
	
	// error caused by servo lag
	temporal_var = trace(cov_eps); // verif car beaucoup trop grand	
	
	// some inits
	nmodes = dm._nact(sum); // # actuators
	mInf = array(0.0f,sim._size,sim._size,nmodes);
	cpt = 0;
	for (nm=1;nm<=ndm;nm++) {
		n1 = dm(nm)._n1; n2 = dm(nm)._n2; nxy = n2-n1+1;
		for (cc=1;cc<=dm._nact(nm);cc++) {
			cpt ++;
			scommand = float(modToAct(indexDm(1,nm):indexDm(2,nm),cc)); 
			// compute the influence matrix with computed DM shape from given commands
			// # modes miroirs
			mInf(n1:n2,n1:n2,cpt) = comp_dm_shape(nm,&scommand)*ipupil(n1:n2,n1:n2);
		}
	}
	
	
	n1     = dm(1)._n1;
	n2     = dm(1)._n2;
	sz = n2-n1+1;
	pupd = sim.pupildiam;
	ind0 = (sz - pupd)/2;
	n1 += ind0;
	n2 -= ind0;
	
	dphiPara = dphiAlias = 0.;
	
	
	// compture the parallel component of the residual phase structure function
	// eq. 36
	/*
	 // slow version
	 for (i=1;i<=nmodes;i++) {
	 for (j=1;j<=i;j++) {
	 write,format="\rComputing U%d%d",i,j;
	 tmp = calc_uij(mInf(n1:n2,n1:n2,i),mInf(n1:n2,n1:n2,j),ipupil(n1:n2,n1:n2),den);
	 // pure parallel
	 dphiPara += (cov_eps(i,j)+cov_eps(j,i))*tmp;
	 // aliaising component
	 dphiAlias += (cov_alias(i,j)+cov_alias(j,i))*tmp;
	 }
	 }
	 */
	
	
	if (is_void(cov_basis) || (cov_basis != "vij")) {
	 
		// fast version
		write,"modes fft computation\n";
		// FFT(M(i)) = ftmi   M(i)  mode i
		// M(i)
		
		ftmi = array(complex,[3,2*pupd,2*pupd,nmodes]);
		mis  = array(float,[3,2*pupd,2*pupd,nmodes]);
		for (i=1;i<=nmodes;i++) {
			mis(1:pupd,1:pupd,i) = (mInf*ipupil)(n1:n2,n1:n2,i)(:pupd, :pupd);
			ftmi(,,i) = fft(mis(,,i),1);
			write,format=" \rComputing fft of mode %d",i;
		}
		
		p                = array(float,[2,2*pupd,2*pupd]);
		p(1:pupd,1:pupd) = ipupil(n1:n2,n1:n2)(:pupd, :pupd);
		conjftp          = conj(fft(p,1));


		 write,"Uij computation\n";
		for (i=1;i<=nmodes;i++) {
			for (j=1;j<=i;j++) {
				write,format=" \rComputing U%d%d",i,j;
				tmp = calc_uijf(ftmi(,,i),ftmi(,,j),mis(,,i),mis(,,j),den,conjftp);
				if (i != j) {
					// pure parallel
					dphiPara += (cov_eps(i,j)+cov_eps(j,i))*tmp; //Eq.36
					// aliasing component
					dphiAlias += (cov_alias(i,j)+cov_alias(j,i))*tmp;
				} else {
					// pure parallel
					dphiPara += cov_eps(i,j)*tmp;
					// aliasing component
					dphiAlias += cov_alias(i,j)*tmp;
				}
			}
		 }
		
		para  = dphiPara;
		alias = dphiAlias;
		dphi = para + alias;
		//error,"test uij";
	} else { 
		
		write,"modes computation\n";
		mis = array(float,[3,2*pupd,2*pupd,nmodes]);
		for (i=1;i<=nmodes;i++) {
			mis(1:pupd,1:pupd,i) = (mInf*ipupil)(n1:n2,n1:n2,i)(:pupd, :pupd);
			write,format=" \rComputing mode %d",i;
		}
		
		ftmi = array(complex,[3,2*pupd,2*pupd,nmodes]);
		p = array(float,[2,2*pupd,2*pupd]);
		p(1:pupd,1:pupd) = ipupil(n1:n2,n1:n2)(:pupd, :pupd);
		conjftp = conj(fft(p,1));
		
		// very fast version
		vii = array(float, [3, dimsof(p)(2), dimsof(p)(3), nmodes]);
		
		// Some inits 
		dim  = dimsof(mis);
		dim1 = dim(2);
		dim2 = dim(3);
		dim3 = dim(4);
		
		// Calculation of the eigen values/vectors 
		l = SVdec(cov_eps + cov_alias, u, vt, full=1);
		
		// reconstruction: alias+mes = (u(,+) * diag(l)(+,))(,+) * vt(+,)
		
		write, "\nVii computation\n";
		vii = mis(, , +) * vt(,+);
		
		for (i=1;i<=nmodes;i++) {
			ftmi(, , i) = fft(vii(, , i), 1);
			write,format=" \rComputing fft of mode %d",i;
		}	
		write, "";
		dphi = dphiOrtho * 0.;
		
		
		for (i=1;i<=nmodes;i++) {
			write, format=" \rComputing V%d%d mode", i;
			modei = vii(..,i);
			
			tmp = calc_Viif(ftmi(, , i), ftmi(, , i), modei, modei, den, conjftp);
			
			dphi  += tmp * l(i);
		}
		write, "\n";
		
	}
	
	// here we create a mask that defines the support on which the phase structure functions or non-nil
	p = array(float,[2,2*pupd,2*pupd]);
	p(1:pupd,1:pupd)  = ipupil(n1:n2,n1:n2)(:pupd, :pupd);
	
	//den  = correlate(p,p).re;
	mask = den > max(den)*1.e-7;
	mask = 1-mask;
	den  = den(where(mask));
	
	// building the otf following : OTF = exp(-1/2 * Dphi)
	//tmp2 = exp(-0.5*para)*exp(-0.5*alias)*exp(-0.5*(dphiOrtho)*(2*pi/(*target.lambda)(1))^2);
	
	if (gamma_eps) {
		tmp = exp(-0.5*(dphi+2.*gamma_eps))*exp(-0.5*(dphiOrtho)*(2.*pi/(*target.lambda)(1))^2);
	} else {
		tmp = exp(-0.5*dphi)*exp(-0.5*(dphiOrtho)*(2.*pi/(*target.lambda)(1))^2);
	}

	
	//tmp = exp(-0.5*eclat(rebin_n(dphi_exact, 256)))*(2.*pi/(*target.lambda)(1))^2;
	
	
	// clean-up the otf beyond the cut-off frequency
	tmp(where(mask)) = 0.;
	//tmp2(where(mask)) = 0.;
	
	sz = dimsof(tmp)(2);
	
	// building a larger support before backward transform
	// this is where you take into account the sampling of the PSF
	// basically shannon sampling means an array of [4 x cut-off freq (in pixels)] x [4 x cut-off freq]
	// here we depend on yao and actually everything is fixed by sim._size
	fto_turb = array(float,[2,sim._size,sim._size]);
	fto_turb(1:sz,1:sz) = eclat(tmp);
	
	// recentering the otf on its final support
	fto_turb = eclat(roll(fto_turb,[sim._size/2-sz/2,sim._size/2-sz/2]));
	//fto_turb = eclat(fto_turb)
	
	// building the telescope otf
	// the pixel size is fixed by the simulation parameters
	fto_tel   = telfto(1.65,0.01,8.,0.1125,(*target.lambda)(1)/tel.diam/4.85/(float(sim._size)/sim.pupildiam),sim._size);
	
	// building the various PSFs
	psftel = eclat(abs(fft(fto_tel,-1)));
	psf    = eclat(abs(fft(fto_tel*fto_turb,-1)));
	
	//PSF from yao
	psftest= imav(,,1,1)/sum(imav(,,1,1));
	
	// re-normalization of the reconstructed PSF and the telescope PSF
	psfrec = psf/sum(psf);
	psftel = psftel/sum(psftel);
	
	difract = circavg(psftel,middle=1);
	
	if (disp) {
		fma; limits,,5;
		plg, difract/max(difract),color="red";
		plg, circavg(psftest,middle=1)/max(difract);
		plg, circavg(psfrec,middle=1)/max(difract),color="green";
		xytitles, "pixels","Strehl ratio";
		pltitle, "PSF circ avg";
	}
	
	//error;
	return psfrec;
}

func test_apety(void)
{
	// run this line after line
	// the yao simulation is long ...
	// use this batch only to generate new dphi_ortho and mesAlias
	extern dphi_tot, dphi_ortho, dphi_para, 
	       smes_alias, smes_nonoise, matPass, den, conjftpup, iMat;
	
	fma; plg,[1]; fma;
	
	// reading parfile
	aoread, "nici.par";
	//aoread, "sh6x6_svipc.par";
	//aoread, "sh16x16_svipc.par";
	
	aoinit,clean=1,disp=0;
	aoloop,savecb=1,disp=0;
	
	// r0 prop lambda^(6/5);
	// r0at05mic = tel.diam/atm.dr0at05mic;
	// r0im = ((*target.lambda)(1)/0.5)^(6/5.)*r0at05mic;
	// pixsize = tel.diam/sim.pupildiam;
	// r0pix = r0im/pixsize;
	// dphiO = get_ortho(r0pix,50)
	
	size   = sim._size;
	nmodes = dm._nact(sum);
	mInf   = array(0.0f,size,size,nmodes);
	cpt    = 0;
	
	for (nm=1 ; nm<=ndm ; nm++) {
		n1 = dm(nm)._n1; n2 = dm(nm)._n2;
		
		for (cc=1 ; cc<=dm._nact(nm) ; cc++) {
			cpt ++;
			scommand = float(modToAct(indexDm(1, nm):indexDm(2, nm), cc));
			mInf(n1:n2, n1:n2, cpt) = comp_dm_shape(nm, &scommand);
		}
	}
	
	// computing the projection matrix from phase screen to modes i.e. phase2modes
	// and init some variables (for use in yao-apety.i)
	valid_pix  = where(ipupil);
	tmp        = mInf(*, )(valid_pix, );
	matPass    = (LUsolve(tmp(+, ) * tmp(+, ))(+, ) * tmp(*, )(, +));
	smes_alias = smes_nonoise \
	           = array(0., [2, dimsof(iMat)(2), 1]);
	dphi_ortho = 0.;
	dphi_tot   = 0.;
	dphi_para  = 0.;
	gamma_eps  = 0.;
	
	n1   = dm(1)._n1;
	n2   = dm(1)._n2;
	sz   = n2-n1+1;
	pupd = sim.pupildiam;
	ind0 = (sz - pupd) / 2;
	n1  += int(ind0);
	n2  -= int(ind0);
	
	pup    = ipupil(n1:n2, n1:n2);
	pix    = where(pup);
	
	tmp = array(float, [2, 2*pupd, 2*pupd]);
	tmp(:pupd, :pupd) = pup(:pupd, :pupd);
	
	// compute the pupil ft conjugate (for fast dphi computation)
	conjftpup = conj(fft(tmp, 1));
	
	// compute the pupil autocorrelation in order the get the u_ij computation faster
	den  = correlate(tmp, tmp).re;

	// AO loop ... acquiring circular buffers
	// + computation dphi_ortho and the real dphi_tot for comparison
	//for (i=1;i<=loop.niter;i++) go,1;
	go, all=1;
	
	smes = smes_alias(, 2:);
	// compute the covariance of the aliasing measurments
	covmes_alias = (smes(,+)*smes(,+))/dimsof(smes)(3);
	
	
	// and now ... reconstruct psf with display ...
	// in black the real psf
	// in green the reconstructed psf
	// in red the telescope only psf
	
	
	error;
	psf = get_psf(cbmes, dphi_ortho/loop.niter, covmes_alias, den, gamma_eps=, cov_basis="vij", noise_method="DSP", disp=1);
	
	return psf;
}