/*******************************************************************************
 *
 * NAME: DDAlgorithm.h
 *
 * DESCRIPTION: Base class for all algorithms.
 *
 *  Created on: October 19, 2020
 *      Author: Sam Anderson
 *
 *******************************************************************************/

#ifndef DDAlgorithm_H_
#define DDAlgorithm_H_

#include <boost/multi_array.hpp>

#include <DDataset.hpp>
#include <DDProcess.h>

using namespace std;

enum class flags {
	ATMFAIL=	1,
	LAND=		2,
	PRODWARN=	4,
	HIGLINT=	8,
	HILT=		16,
	HISATZEN= 	32,
	COASTZ=		64,
	SPARE=		128,
	STRAYLIGHT=	256,
	CLDICE=		512,
	COCCOLITH= 	1024,
	TURBIDW=	2048,
	HISOLZEN=	4096,
	SPARE1=		8192,
	LOWLW=		16384,
	CHLFAIL=	32768,
	NAVWARN=	65536,
	ABSAER=		131072,
	SPARE2=		262144,
	MAXAERITER= 524288,
	MODGLINT=	1048576,
	CHLWARN=	2097152,
	ATMWARN=	4194304,
	SPARE3=		8388608,
	SEAICE=		16777216,
	NAVFAIL=	33554432,
	FILTER=     67108864,
	SPARE5=     134217728,
	BOWTIEDEL=  268435456,
	HIPOL=      536870912,
	PRODFAIL=   1073741824
};

/**
 *  Algorithm -- virtual base class for algorithm processes
 */

class DDAlgorithm
{
public:

	DDAlgorithm();
	virtual ~DDAlgorithm ();
	virtual int initialize(map<string, ddata*> imap) {return 0;};
	virtual map<string, ddata*> process(vector<size_t> start, vector<size_t> count,
			map<string, ddata*> imap) {map<string, ddata*> map; return map;};
	virtual vector<string>  get_products() {vector<string> list; return list;};

protected:

	size_t  lines_;
	size_t  pixels_;
	float   threshsolz_;
	float   threshsenz_;
	float   threshglint_;
    float   rfl_[NTWL];
    float   rfla_[NTWL][3][3];
    float   gasc_[NTWL];
	float 	lat_;
	float 	lon_;
	float 	solz_;
	float 	senz_;
	float 	raa_;
	float 	height_;
	float 	ws_;
	float 	pwv_;
	float 	oz_;
	float 	ps_;
	int		month_;
	bool    bgascorrect_;
	bool    bmaskglint_;
	bool    bmaskcloud_;
	bool    bmasksolz_;
	bool    bmasksenz_;
	bool    btest_;

	unsigned char	cloud_mask_;
    short   		qual_flag_;
    short   		aerosol_type_;
    short   		error_flag_;
    unsigned int	l2_flags_;
    float   scatter_ang_;
    float   glint_ang_;
    float   sse_;
    float   fmf_;
    float   aot_550_;
    float   ae1_;
    float   ae2_;
    float   ndv_;
    float   chlor_;
    float   ssa_[NLWL+1];
    float   sr_[NLWL+1];
    float   aot_[NOWL+1];

	int get_inputs( vector<size_t> start, vector<size_t> count,
			map<string, ddata*> imap );

	map<string, ddata*> set_outputs();
	map<string, ddata*> set_fills();
};


#endif /* DDAlgorithm_H_ */
