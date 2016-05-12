#ifndef _DEFS
#define _DEFS

#define REFPROP_ACMODEL


#define _R410A
//#define _R407C
//#define _R22

#ifdef _R410A
#define _RefMix
#endif

#ifdef _R407C
#define _RefMix
#endif

// NIST chiller model expansion valve superheats.
// These should be transferred into component parameters.
#define TSH_LP   2.0
#define TSH_HP   13.0

#endif
