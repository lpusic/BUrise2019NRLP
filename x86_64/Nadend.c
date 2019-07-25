/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__Nadend
#define _nrn_initial _nrn_initial__Nadend
#define nrn_cur _nrn_cur__Nadend
#define _nrn_current _nrn_current__Nadend
#define nrn_jacob _nrn_jacob__Nadend
#define nrn_state _nrn_state__Nadend
#define _net_receive _net_receive__Nadend 
#define _f_evaluate_fct _f_evaluate_fct__Nadend 
#define evaluate_fct evaluate_fct__Nadend 
#define states states__Nadend 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gnadend _p[0]
#define gl _p[1]
#define el _p[2]
#define ina _p[3]
#define il _p[4]
#define m _p[5]
#define h _p[6]
#define ena _p[7]
#define Dm _p[8]
#define Dh _p[9]
#define mexp _p[10]
#define _g _p[11]
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_Exp(void);
 static void _hoc_evaluate_fct(void);
 static void _hoc_states(void);
 static void _hoc_vtrap(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_Nadend", _hoc_setdata,
 "Exp_Nadend", _hoc_Exp,
 "evaluate_fct_Nadend", _hoc_evaluate_fct,
 "states_Nadend", _hoc_states,
 "vtrap_Nadend", _hoc_vtrap,
 0, 0
};
#define Exp Exp_Nadend
#define vtrap vtrap_Nadend
 extern double Exp( double );
 extern double vtrap( double , double );
 /* declare global and static user variables */
#define htau htau_Nadend
 double htau = 0;
#define hexp hexp_Nadend
 double hexp = 0;
#define hinf hinf_Nadend
 double hinf = 0;
#define mtau mtau_Nadend
 double mtau = 0;
#define minf minf_Nadend
 double minf = 0;
#define usetable usetable_Nadend
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_Nadend", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "mtau_Nadend", "ms",
 "htau_Nadend", "ms",
 "gnadend_Nadend", "mho/cm2",
 "gl_Nadend", "mho/cm2",
 "el_Nadend", "mV",
 "ina_Nadend", "mA/cm2",
 "il_Nadend", "mA/cm2",
 0,0
};
 static double delta_t = 1;
 static double h0 = 0;
 static double m0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "minf_Nadend", &minf_Nadend,
 "hinf_Nadend", &hinf_Nadend,
 "hexp_Nadend", &hexp_Nadend,
 "mtau_Nadend", &mtau_Nadend,
 "htau_Nadend", &htau_Nadend,
 "usetable_Nadend", &usetable_Nadend,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"Nadend",
 "gnadend_Nadend",
 "gl_Nadend",
 "el_Nadend",
 0,
 "ina_Nadend",
 "il_Nadend",
 0,
 "m_Nadend",
 "h_Nadend",
 0,
 0};
 static Symbol* _na_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 12, _prop);
 	/*initialize range parameters*/
 	gnadend = 0.0117;
 	gl = 5e-05;
 	el = -70;
 	_prop->param = _p;
 	_prop->param_size = 12;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _Nadend_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("na", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 12, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Nadend /Users/luka.pusic/Desktop/Turi_et_al_2018/x86_64/Nadend.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_minf;
 static double *_t_hinf;
 static double *_t_hexp;
 static double *_t_mtau;
 static double *_t_htau;
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_evaluate_fct(double);
static int evaluate_fct(double);
static int states();
 static void _n_evaluate_fct(double);
 
static int  states (  ) {
   evaluate_fct ( _threadargscomma_ v ) ;
   h = h + hexp * ( hinf - h ) ;
   
/*VERBATIM*/
	return 0;
  return 0; }
 
static void _hoc_states(void) {
  double _r;
   _r = 1.;
 states (  );
 hoc_retpushx(_r);
}
 static double _mfac_evaluate_fct, _tmin_evaluate_fct;
 static void _check_evaluate_fct();
 static void _check_evaluate_fct() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_dt;
  static double _sav_celsius;
  if (!usetable) {return;}
  if (_sav_dt != dt) { _maktable = 1;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_evaluate_fct =  - 200.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_evaluate_fct)/300.; _mfac_evaluate_fct = 1./_dx;
   for (_i=0, _x=_tmin_evaluate_fct; _i < 301; _x += _dx, _i++) {
    _f_evaluate_fct(_x);
    _t_minf[_i] = minf;
    _t_hinf[_i] = hinf;
    _t_hexp[_i] = hexp;
    _t_mtau[_i] = mtau;
    _t_htau[_i] = htau;
   }
   _sav_dt = dt;
   _sav_celsius = celsius;
  }
 }

 static int evaluate_fct(double _lv){ _check_evaluate_fct();
 _n_evaluate_fct(_lv);
 return 0;
 }

 static void _n_evaluate_fct(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_evaluate_fct(_lv); return; 
}
 _xi = _mfac_evaluate_fct * (_lv - _tmin_evaluate_fct);
 if (isnan(_xi)) {
  minf = _xi;
  hinf = _xi;
  hexp = _xi;
  mtau = _xi;
  htau = _xi;
  return;
 }
 if (_xi <= 0.) {
 minf = _t_minf[0];
 hinf = _t_hinf[0];
 hexp = _t_hexp[0];
 mtau = _t_mtau[0];
 htau = _t_htau[0];
 return; }
 if (_xi >= 300.) {
 minf = _t_minf[300];
 hinf = _t_hinf[300];
 hexp = _t_hexp[300];
 mtau = _t_mtau[300];
 htau = _t_htau[300];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 minf = _t_minf[_i] + _theta*(_t_minf[_i+1] - _t_minf[_i]);
 hinf = _t_hinf[_i] + _theta*(_t_hinf[_i+1] - _t_hinf[_i]);
 hexp = _t_hexp[_i] + _theta*(_t_hexp[_i+1] - _t_hexp[_i]);
 mtau = _t_mtau[_i] + _theta*(_t_mtau[_i+1] - _t_mtau[_i]);
 htau = _t_htau[_i] + _theta*(_t_htau[_i+1] - _t_htau[_i]);
 }

 
static int  _f_evaluate_fct (  double _lv ) {
   double _lq10 , _ltinc , _lalpha , _lbeta ;
 _lq10 = 1.0 ;
   _ltinc = - dt * _lq10 ;
   _lalpha = 0.1 * vtrap ( _threadargscomma_ - ( _lv + 45.0 ) , 10.0 ) ;
   _lbeta = 4.0 * exp ( - ( _lv + 70.0 ) / 18.0 ) ;
   mtau = 1.0 / ( _lalpha + _lbeta ) ;
   minf = _lalpha * mtau ;
   _lalpha = 0.07 * Exp ( _threadargscomma_ - ( _lv + 70.0 ) / 20.0 ) ;
   _lbeta = 1.0 / ( 1.0 + Exp ( _threadargscomma_ - ( _lv + 40.0 ) / 10.0 ) ) ;
   htau = 1.0 / ( _lalpha + _lbeta ) ;
   hinf = _lalpha * htau ;
   hexp = 1.0 - Exp ( _threadargscomma_ _ltinc / htau ) ;
    return 0; }
 
static void _hoc_evaluate_fct(void) {
  double _r;
    _r = 1.;
 evaluate_fct (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap (  double _lx , double _ly ) {
   double _lvtrap;
 if ( fabs ( _lx / _ly ) < 1e-6 ) {
     _lvtrap = _ly * ( 1.0 - _lx / _ly / 2.0 ) ;
     }
   else {
     _lvtrap = _lx / ( Exp ( _threadargscomma_ _lx / _ly ) - 1.0 ) ;
     }
   
return _lvtrap;
 }
 
static void _hoc_vtrap(void) {
  double _r;
   _r =  vtrap (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double Exp (  double _lx ) {
   double _lExp;
 if ( _lx < - 100.0 ) {
     _lExp = 0.0 ;
     }
   else {
     _lExp = exp ( _lx ) ;
     }
   
return _lExp;
 }
 
static void _hoc_Exp(void) {
  double _r;
   _r =  Exp (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ hoc_execerror("Nadend", "cannot be used with CVODE"); return 0;}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
 {
   m = minf ;
   h = hinf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ena = _ion_ena;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ina = gnadend * minf * minf * minf * h * ( v - ena ) ;
   il = gl * ( v - el ) ;
   }
 _current += ina;
 _current += il;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ena = _ion_ena;
 _g = _nrn_current(_v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ena = _ion_ena;
 { error =  states();
 if(error){fprintf(stderr,"at line 83 in file Nadend.mod:\n	SOLVE states\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
   _t_minf = makevector(301*sizeof(double));
   _t_hinf = makevector(301*sizeof(double));
   _t_hexp = makevector(301*sizeof(double));
   _t_mtau = makevector(301*sizeof(double));
   _t_htau = makevector(301*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/luka.pusic/Desktop/Turi_et_al_2018/mechanisms/Nadend.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "\n"
  "Sodium current for the dendrites\n"
  "\n"
  "References:\n"
  "\n"
  "1.	Martina, M., Vida, I., and Jonas, P.  Distal initiation and active\n"
  "	propagation of action potentials in interneuron dendrites,\n"
  "	Science, 287:295-300, 2000.\n"
  "\n"
  "			soma	axon-lacking dend	axon-bearing dend\n"
  "Na+	gmax	    107 ps/um2	   117 ps/um2		   107 ps/um2\n"
  "	slope 	    10.9 mV/e	   11.2 mV/e		   11.2 mV/e\n"
  "	V1/2        -37.8 mV       -45.6 mV                -45.6 mV\n"
  "\n"
  "\n"
  "\n"
  "2.	Marina, M. and Jonas, P.  Functional differences in Na+ channel\n"
  "	gating between fast-spiking interneurons and principal neurones of rat\n"
  "	hippocampus, J. Physiol., 505.3:593-603, 1997.\n"
  "\n"
  "*Note* The interneurons here are basket cells from the dentate gyrus.\n"
  "\n"
  "Na+	Activation V1/2				-25.1 mV\n"
  "	slope			 		11.5\n"
  "	Activation t (-20 mV)	 		0.16 ms\n"
  "	Deactivation t (-40 mV)	 		0.13 ms\n"
  " 	Inactivation V1/2			-58.3 mV\n"
  "	slope			 		6.7\n"
  "	onset of inactivation t (-20 mV)	1.34 ms\n"
  "	onset of inactivation t (-55 mV)	18.6 ms\n"
  "	recovery from inactivation t		2.0 ms\n"
  "	(30 ms conditioning pulse)\n"
  "	recovery from inactivation t		2.7 ms\n"
  "	(300 ms conditioning pulse)\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "UNITS {\n"
  "    (mA) = (milliamp)\n"
  "    (mV) = (millivolt)\n"
  "}\n"
  " \n"
  "NEURON {\n"
  "    SUFFIX Nadend\n"
  "    USEION na READ ena WRITE ina\n"
  "    NONSPECIFIC_CURRENT il\n"
  "    RANGE gnadend, gl, el, ina\n"
  "    GLOBAL minf, hinf, hexp, mtau, htau\n"
  "}\n"
  " \n"
  "INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }\n"
  " \n"
  "PARAMETER {\n"
  "    v                (mV)\n"
  "    celsius = 24     (degC)\n"
  "    dt               (ms)\n"
  "    gnadend = .0117  (mho/cm2)\n"
  "    ena     = 90     (mV)\n"
  "    gl      = .00005 (mho/cm2)\n"
  "    el      = -70    (mV)\n"
  "}\n"
  " \n"
  "STATE { m h }\n"
  " \n"
  "ASSIGNED {\n"
  "    ina (mA/cm2)\n"
  "    il  (mA/cm2)\n"
  "    minf \n"
  "	mexp \n"
  "	hinf \n"
  "	hexp\n"
  "	mtau (ms)\n"
  "	htau (ms)\n"
  "}\n"
  " \n"
  "INITIAL {\n"
  "	m = minf\n"
  "	h = hinf\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states\n"
  "	ina = gnadend*minf*minf*minf*h*(v - ena)    \n"
  "	il = gl*(v - el)\n"
  "}\n"
  "\n"
  "PROCEDURE states() {	:exact when v held constant\n"
  "	evaluate_fct(v)\n"
  "	h = h + hexp*(hinf - h)\n"
  "	VERBATIM\n"
  "	return 0;\n"
  "	ENDVERBATIM \n"
  "}\n"
  "UNITSOFF\n"
  "\n"
  "PROCEDURE evaluate_fct(v(mV)) {  :Computes rate and other constants at \n"
  "		:current v.\n"
  "        \n"
  "        :Call once from HOC to initialize inf at resting v.\n"
  "        LOCAL q10, tinc, alpha, beta\n"
  "        TABLE minf, hinf, hexp, mtau, htau DEPEND dt, celsius FROM -200 TO 100 WITH 300\n"
  "		\n"
  "		:q10   = 3^((celsius - 24)/10)\n"
  "		q10   = 1	: BPG\n"
  "		tinc  = -dt*q10\n"
  "		alpha = 0.1*vtrap(-(v+45),10)\n"
  "		beta  = 4*exp(-(v+70)/18)\n"
  "		mtau  = 1/(alpha + beta)\n"
  "		minf  = alpha*mtau\n"
  "		alpha = 0.07*Exp(-(v+70)/20)\n"
  "		beta  = 1/(1+Exp(-(v+40)/10))\n"
  "		htau  = 1/(alpha + beta)\n"
  "		hinf  = alpha*htau\n"
  "		hexp  = 1-Exp(tinc/htau)\n"
  "}\n"
  "\n"
  "FUNCTION vtrap(x,y) {\n"
  "	:Traps for 0 in denominator of rate eqns.\n"
  "	if (fabs(x/y) < 1e-6) {\n"
  "		vtrap = y*(1 - x/y/2)\n"
  "	}else{\n"
  "		vtrap = x/(Exp(x/y) - 1)\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION Exp(x) {\n"
  "	if (x < -100) {\n"
  "		Exp = 0\n"
  "	}else{\n"
  "		Exp = exp(x)\n"
  "	}\n"
  "}\n"
  "UNITSON\n"
  ;
#endif
