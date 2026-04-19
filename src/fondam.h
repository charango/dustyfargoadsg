/** \file fondam.h

Contains fondamental constants used thorough the code.
*/

#define		G	1.0
#define	      PI4	12.56637061435917295376
#define	       PI	3.14159265358979323844
#define   CPUOVERLAP    6   /* Zeus-like overlap kernel: value should
				be tested with an indiscernability
				test against varying number of
				CPUs. At least CPUOVERLAP should be 6
				when the energy equation, thermal
				diffusion, viscosity, and a second fluid are 
				included. Note however that inclusion of the gas 
				indirect term and/or gas self-gravity will disrupt 
				strict indiscernability (noise introduced at ~machine precision) */
#define      MU         1.0  /* Mean molecular weight */
#define      R          1.0  /* Universal Gas Constant in code units */
