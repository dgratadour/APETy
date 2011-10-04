/*
 * yao_funcs.i
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


func mult_wfs(iter,disp=)
/* DOCUMENT func mult_wfs(iter,disp=)
   Goes through all WFS and concatenate the resulting measurement vectors.
   SEE ALSO:
 */
{
  if (iter == 1) write,"Using the modified mult_wfs version";
  mes = [];
  for (ns=1;ns<=nwfs;ns++) {

    offsets = wfs(ns).gspos;
    phase  = get_phase2d_from_optics(ns,"wfs");
    phase  += get_turb_phase(iter,ns,"wfs");
    // only look at DMs if not running in open loop
    if (loop.method != "open-loop") {
        phase  += get_phase2d_from_dms(ns,"wfs"); 
		// enleve la phase parallele  ====> obtient la phase residuelle
      }

    if (wfs(ns).correctUpTT) {
      phase = correct_uplink_tt(phase,ns);
    }

    //APETY HACK
    // computes the total residual phase structure function from the phase screen directly
    // i.e. THE TRUTH !!! ... for test purposes only off-course !
    // don't forget to init dphi_ortho to 0 before starting the loop : dphi_ortho=0.
    // divide by loop.niter (e.g. 10000) when loop is over    
	extern dphi_tot, dphi_para, dphi_ortho, gamma_eps, \
	       smes_alias, smes_nonoise, matPass, den, conjftpup, wfs;
    
    pix = where(ipupil);
	  
    //tmp = float(phase);
    tmp       = float(phase-phase(pix)(avg)); // remove piston mode
    vect      = (tmp*ipupil)(pix)(+)*matPass(,+);  // projection phase sur modes miroir => vect: coefficient sur chaque mode
    dm_p      = vect(+)*mInf(,,+); // phase residuelle para
    tmp_ortho = float(tmp - dm_p); // phase ortho
 
    n1   = dm(1)._n1;
    n2   = dm(1)._n2;
    sz   = n2-n1+1;
    pupd = sim.pupildiam;
    ind0 = (sz - pupd)/2;
    n1  += ind0;
    n2  -= ind0;
	  
	  
	  //dphi_ortho += calc_dphi(tmp_ortho(n1:n2,n1:n2),ipupil(n1:n2,n1:n2),den);
    dphi_tot   += calc_dphif(tmp(n1:n2,n1:n2),ipupil(n1:n2,n1:n2),den,conjftpup);
    dphi_ortho += calc_dphif(tmp_ortho(n1:n2,n1:n2),ipupil(n1:n2,n1:n2),den,conjftpup);
    dphi_para  += calc_dphif(dm_p(n1:n2,n1:n2),ipupil(n1:n2,n1:n2),den,conjftpup);	  
	gamma_eps  += calc_gamma(dm_p(n1:n2,n1:n2), tmp_ortho(n1:n2,n1:n2),ipupil(n1:n2,n1:n2),den,conjftpup);
    if (is_void(dphi_exact)) {
		dphi_exact = calc_dphis(dm_p(n1:n2,n1:n2));
	} else {
		dphi_exact += calc_dphis(dm_p(n1:n2,n1:n2));
	}

	  
	nsave = wfs(ns).noise;
    // switch off noise
    wfs(ns).noise = 0;
    // get aliasing & noise free measurements
    if (wfs(ns).type == "hartmann" ) {
      if (wfs(ns).disjointpup) {
        tmpmes = sh_wfs(disjointpup(,,ns),tmp_ortho,ns);
        //tmpmes2 = sh_wfs(disjointpup(,,ns),tmp,ns);
      } else {
        tmpmes = sh_wfs(ipupil,tmp_ortho,ns);
        //tmpmes2 = sh_wfs(ipupil,tmp,ns);
      }
    } else if (wfs(ns).type == "curvature") {
      tmpmes = curv_wfs(pupil,tmp_ortho,ns);
      tmpmes2 = curv_wfs(pupil,tmp,ns);
    } else if (wfs(ns).type == "pyramid") {
      tmpmes = pyramid_wfs(pupil,tmp_ortho,ns);
      tmpmes2 = pyramid_wfs(pupil,tmp,ns);
    } else if (wfs(ns).type == "zernike") {
      tmpmes = zernike_wfs(ipupil,tmp_ortho,ns);
      tmpmes2 = zernike_wfs(ipupil,tmp,ns);
    } else {
      // assign user_wfs to requested function/type:
      cmd = swrite(format="user_wfs = %s",wfs(ns).type);
      include,[cmd],1;
      tmpmes = user_wfs(ipupil,tmp_ortho,ns);
      tmpmes2 = user_wfs(ipupil,tmp,ns);
    }
    // switch back to original noise value
    wfs(ns).noise = nsave;
    grow,smes_alias,tmpmes;
    //grow,smes_nonoise,tmpmes2;
    
    // END OF APETY HACK

    
    // get the measurements:
    if (wfs(ns).type == "hartmann" ) {
      if (wfs(ns).disjointpup) {
        smes = sh_wfs(disjointpup(,,ns),phase,ns);
      } else {
        smes = sh_wfs(ipupil,phase,ns);
      }
    } else if (wfs(ns).type == "curvature") {
      smes = curv_wfs(pupil,phase,ns);
    } else if (wfs(ns).type == "pyramid") {
      smes = pyramid_wfs(pupil,phase,ns);
    } else if (wfs(ns).type == "zernike") {
      smes = zernike_wfs(ipupil,phase,ns);
    } else {
      // assign user_wfs to requested function/type:
      cmd = swrite(format="user_wfs = %s",wfs(ns).type);
      include,[cmd],1;
      smes = user_wfs(ipupil,phase,ns);
    }

    if (sim.svipc) {
      // hack, this has nothing to do here. FIXME
      shm_write,shmkey,swrite(format="wfs%d_image",ns),wfs(ns)._fimage;
    }

    // subtract the reference vector for this sensor:
    if (wfs(ns)._cyclecounter == 1) {
      smes = smes - *wfs(ns)._refmes;
    }
	  
    // compute the TT and subtract if required:
    wfs(ns)._tt(1) = sum( smes * (*wfs(ns)._tiprefvn) );
    wfs(ns)._tt(2) = sum( smes * (*wfs(ns)._tiltrefvn) );
    if (wfs(ns).filtertilt) {
      smes = smes - wfs(ns)._tt(1) * (*wfs(ns)._tiprefv) \
        - wfs(ns)._tt(2) * (*wfs(ns)._tiltrefv);
    }
    if (wfs(ns)._cyclecounter == 1) {
      wfs(ns)._lastvalidtt = wfs(ns)._tt;
    }

    grow,mes,smes;
  }
  return mes;
}
