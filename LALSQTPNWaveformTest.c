/**
 * @file LALSQTPNWaveformTest.c
 *		The user interface for the SpinQuadTaylor program.
 * This file is an example howto use the SpinQuadTaylor program.
 * The input parameters are:<br/>
 * <em>Waveform parameters:</em>(note:\f$i=1,2\f$ and \f$j=x,y,z\f$)
 *	<ul>
 *	<li>masses of the BHs (or NSs) \f$m_i\f$ in \f$M_\odot\f$</li>
 *	<li>the spin components \f$\chi_{ij}\f$, the values of \f$\sqrt{\sum_j\chi_{ij}}\f$, are between 0 and 1</li>
 *	<li>the quadrupole parameters \f$w_i\in(4,8)\f$ for NSs [1] and \f$w_i=1\f$ for BHs[2] are 1 (default 1)</li>
 *	<li>the inclination (angle between the line of sight and Newtonian orbital angular momentum) \f$\iota\f$ in \f$rad\f$</li>
 *	<li>the initial frequency \f$f_L\f$ in \f$Hz\f$</li>
 *	<li>the distance \f$d\f$ in \f$Mpc\f$</li>
 *	<li>the sampling time \f$t_s\f$ in \f$s\f$</li>
 *	</ul>
 *<em>Program parameters:</em>
 *	<ul>
 *	<li>the name of the output file (default out.dat)</li>
 *	<li>the PN order (newtonian, oneHalfPN, onePN, onePointFivePN, twoPN(default), twoPointFivePN, threePN, threePointFivePN)</li>
 *	<li>level of accuracy in including spin and quadrupole contributions
 *	(NoInter, SOInter, SSInter, SSselfInter, QMInter, AllInter(default))</li>
 *	</ul>
 *	The output file contains three coloums: elapsed time, \f$h_+\f$, \f$h_x\f$.
 *	\f{center}
 *	\begin{gather*}
 *		h_+=a_1\cos\left(2\alpha\right)\cos\left(2\Phi\right)-a_2\sin\left(2\alpha\right)\sin\left(2\Phi\right),\\
 *		h_x=a_1\sin\left(2\alpha\right)\cos\left(2\Phi\right)+a_2\cos\left(2\alpha\right)\sin\left(2\Phi\right)
 *		\end{gather*}
 *	\f}
 *	with \f$a_i\f$ amplitudes, \f$\alpha\f$ polarization shift, \f$\Phi\f$ phase (Eq. (4.28)-(4.30) of [3] up to leading order (The \f$\Phi\f$ is shifted by \f$\pi\f$ with respect to [3]). We note that \f$\Theta=0\f$ in leading order because we use radiation frame).<br />
 *	\f$a_1\f$, \f$a_2\f$ are defined in LALSQTPNGenerator() function, \f$\alpha\f$ and \f$\Phi\f$ phase is defined in LALSQTPNDerivator() function.<br />
 *	<b>References</b><br />
 *	[1] E. Poisson, Phys.Rev. D57 5287 (1998)<br />
 *	[2] K. S. Thorne, Rev.Mod.Phys. 52 299 (1980)<br />
 *	[3] L. E. Kidder, Phys.Rev.D 52, 821 (1995)<br /> 
 * @author László Veréb
 * @date 2010.06.27.
 */

#include <lal/GenerateInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/LALSQTPNWaveformInterface.h>

NRCSID(LALSQTPNWAVEFORMTESTC, "$Id: LALSQTPNWaveformTest.c,v 0.1 2010/05/21");

int lalDebugLevel = 7;	///< the debug level

/** The main program.
 */
int main(int argc, char *argv[]) {

	// variable declaration and initialization
	static LALStatus mystatus;
	CoherentGW thewaveform;
	SimInspiralTable injParams;
	PPNParamStruc ppnParams;
	//const char *filename;
	FILE *outputfile;
	INT4 i, length;
	REAL8 dt;
	char PNString[50];
	char filename[50];

	if (argc != 16) {
		printf(
				"                         1  2  3   4   5   6   7   8   9    10	     11	     12       13 14      15\n");
		printf(
				"Correct parameter order: m1 m2 S1x S1y S1z S2x S2y S2z incl f_lower f_final distance dt PNorder SpinInter\n");
		return (1);
	}
	
	sprintf(PNString, "SpinQuadTaylor%s%s", argv[14], argv[15]);
	sprintf(filename, "%s%s.out", argv[14], argv[15]);
	//filename = argv[15];
	memset(&mystatus, 0, sizeof(LALStatus));
	memset(&thewaveform, 0, sizeof(CoherentGW));
	memset(&injParams, 0, sizeof(SimInspiralTable));
	memset(&ppnParams, 0, sizeof(PPNParamStruc));

	//	setting the parameters
	injParams.mass1 = atof(argv[1]);
	injParams.mass2 = atof(argv[2]);
	injParams.spin1x = atof(argv[3]);
	injParams.spin1y = atof(argv[4]);
	injParams.spin1z = atof(argv[5]);
	injParams.spin2x = atof(argv[6]);
	injParams.spin2y = atof(argv[7]);
	injParams.spin2z = atof(argv[8]);
	injParams.qmParameter1 = 1.;//atof(argv[9]);
	injParams.qmParameter2 = 1.;//atof(argv[10]);
	injParams.inclination = atof(argv[9]);
	injParams.f_lower = atof(argv[10]);
	injParams.f_final = atof(argv[11]);
	injParams.distance = atof(argv[12]);
	ppnParams.deltaT = atof(argv[13]);
	injParams.polarization = 0;
	LALSnprintf(injParams.waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR),
			PNString);

	//interface(&mystatus, &thewaveform, &injParams, &ppnParams);
	LALGenerateInspiral(&mystatus, &thewaveform, &injParams, &ppnParams);
	if (mystatus.statusCode) {
		fprintf( stderr, "LALSQTPNWaveformTest: error generating waveform\n" );
		return mystatus.statusCode;
	}
	// and finally save in a file
	outputfile = fopen(filename, "w");

	length = thewaveform.f->data->length;
	dt      = thewaveform.phi->deltaT;
    REAL8       a1, a2, phi, shift;
    for(i = 0; i < length; i++) {
        a1  = thewaveform.a->data->data[2*i];
        a2  = thewaveform.a->data->data[2*i+1];
        phi     = thewaveform.phi->data->data[i] - thewaveform.phi->data->data[0];
        shift   = thewaveform.shift->data->data[i];

        fprintf(outputfile,"% 10.7e\t% 10.7e\t% 10.7e\n",
            i*dt,
            a1*cos(shift)*cos(phi) - a2*sin(shift)*sin(phi),
            a1*sin(shift)*cos(phi) + a2*cos(shift)*sin(phi));
    }
	fclose(outputfile);
	XLALSQTPNDestroyCoherentGW(&thewaveform);
	puts("Done.");
	LALCheckMemoryLeaks();
	return mystatus.statusCode;
}

