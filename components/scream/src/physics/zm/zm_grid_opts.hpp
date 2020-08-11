#ifndef ZM_GRID_OPTS_HPP
#define ZM_GRID_OPTS_HPP

#include "ekat/scream_assert.hpp"

#include <unordered_map>


#include <string>

const int SCALAR_3D_MID = 0;
const int SCALAR_3D_INT = 1;
const int VECTOR_3D_MID = 2;
const int TRACERS = 3;
const int LINEAR = 4;

using namespace std;
using namespace scream;
using namespace units;

struct GridOpts{
  string name;
  bool isOut;
  const scream::units::Units* unit;
  int field_idx;
};


unordered_map<string, GridOpts> opt_map;

void set_grid_opts_helper(GridOpts O, string n, bool out, const scream::units::Units* unit, int field_idx
                          ){
  O.name = n;
  O.isOut = out;
  O.field_idx = field_idx;
  O.unit = unit;
  opt_map.insert({O.name, O});
}

void set_grid_opts(){
  	
  //t(pcols,pver)  
  GridOpts t;
  set_grid_opts_helper(t, "t", true, &Pa, SCALAR_3D_MID); //temperature(K)

  GridOpts qh;
  set_grid_opts_helper(qh, "qh", true, NULL, SCALAR_3D_MID);

  //@Aaron: I created 'linear' layout for vars of the form x(y)
  GridOpts prec;
  set_grid_opts_helper(prec, "prec", true, NULL, LINEAR);

  GridOpts jctop;
  set_grid_opts_helper(jctop, "jctop", true, NULL, LINEAR);

  GridOpts jcbot;
  set_grid_opts_helper(jcbot, "jcbot", true, NULL, LINEAR);

  GridOpts pblh;
  set_grid_opts_helper(pblh, "pblh", true, NULL, LINEAR);

  GridOpts zm;
  set_grid_opts_helper(zm, "zm", true, NULL, SCALAR_3D_MID);

  GridOpts geos;
  set_grid_opts_helper(geos, "geos", true, NULL, LINEAR);

//zi(pcols,pver+1)
  GridOpts zi;
  set_grid_opts_helper(zi, "zi", true, NULL, SCALAR_3D_INT);

  GridOpts qtnd;
  set_grid_opts_helper(qtnd, "qtnd", true, NULL, SCALAR_3D_MID);

  GridOpts heat;
  set_grid_opts_helper(heat, "heat", true, NULL, SCALAR_3D_MID);

  GridOpts pap;
  set_grid_opts_helper(pap, "pap", true, NULL, SCALAR_3D_MID);

  GridOpts paph;
  set_grid_opts_helper(paph, "paph", true, NULL, SCALAR_3D_INT);

  GridOpts dpp;
  set_grid_opts_helper(dpp, "dpp", true, NULL, SCALAR_3D_MID);

  GridOpts mcon;
  set_grid_opts_helper(mcon, "mcon", true, NULL, SCALAR_3D_MID);

  GridOpts cme;
  set_grid_opts_helper(cme, "cme", true, NULL, SCALAR_3D_MID);

  GridOpts cape;
  set_grid_opts_helper(cape, "cape", true, NULL, LINEAR);

  GridOpts tpert;
  set_grid_opts_helper(tpert, "tpert", true, NULL, LINEAR);

  GridOpts dlf;
  set_grid_opts_helper(dlf, "dlf", true, NULL, SCALAR_3D_MID);

  GridOpts pflx;
  set_grid_opts_helper(pflx, "pflx", true, NULL, SCALAR_3D_MID);

  GridOpts zdu;
  set_grid_opts_helper(zdu, "zdu", true, NULL, SCALAR_3D_MID);

  GridOpts rprd;
  set_grid_opts_helper(rprd, "rprd", true, NULL, SCALAR_3D_MID);
  
  GridOpts mu;
  set_grid_opts_helper(mu, "mu", true, NULL, SCALAR_3D_MID);
  
  GridOpts md;
  set_grid_opts_helper(md, "md", true, NULL, SCALAR_3D_MID);
  
  GridOpts du;
  set_grid_opts_helper(du, "du", true, NULL, SCALAR_3D_MID);
  
  GridOpts eu;
  set_grid_opts_helper(eu, "eu", true, NULL, SCALAR_3D_MID);
  
  GridOpts ed;
  set_grid_opts_helper(ed, "ed", true, NULL, SCALAR_3D_MID);
  
  GridOpts dp;
  set_grid_opts_helper(dp, "dp", true, NULL, SCALAR_3D_MID);
  
  GridOpts dsubcld;
  set_grid_opts_helper(dsubcld, "dsubcld", true, NULL, LINEAR);
  
  GridOpts jt;
  set_grid_opts_helper(jt, "jt", true, NULL, LINEAR);
  
  GridOpts maxg;
  set_grid_opts_helper(maxg, "maxg", true, NULL, LINEAR);
  
  GridOpts ideep;
  set_grid_opts_helper(ideep, "ideep", true, NULL, LINEAR);
  
  GridOpts ql;
  set_grid_opts_helper(ql, "ql", true, NULL, SCALAR_3D_MID);
  
  GridOpts rliq;
  set_grid_opts_helper(rliq, "rliq", true, NULL, SCALAR_3D_MID);
  
  GridOpts landfrac;
  set_grid_opts_helper(landfrac, "landfrac", true, NULL, LINEAR);
  
  GridOpts hu_nm1;
  set_grid_opts_helper(hu_nm1, "hu_nm1", true, NULL, SCALAR_3D_MID);
  
  GridOpts cnv_nm1;
  set_grid_opts_helper(cnv_nm1, "cnv_nm1", true, NULL, SCALAR_3D_MID);
  
  GridOpts tm1;
  set_grid_opts_helper(tm1, "tm1", true, NULL, SCALAR_3D_MID);
  
  GridOpts qm1;
  set_grid_opts_helper(qm1, "qm1", true, NULL, SCALAR_3D_MID);
  
  GridOpts dcape;
  set_grid_opts_helper(dcape, "dcape", true, NULL, LINEAR);
  
  GridOpts q;
  set_grid_opts_helper(q, "q", true, NULL, SCALAR_3D_MID);
  
  GridOpts snow;
  set_grid_opts_helper(snow, "snow", true, NULL, LINEAR);
  
  GridOpts ntprprd;
  set_grid_opts_helper(ntprprd, "ntprprd", true, NULL, SCALAR_3D_MID);
  
  GridOpts ntsnprd;
  set_grid_opts_helper(ntsnprd, "ntsnprd", true, NULL, SCALAR_3D_MID);
  
  //@Aaron: Is vector 3d layout right for pguall(pcols,pver,2)?
  GridOpts pguall;
  set_grid_opts_helper(pguall, "pguall", true, NULL, VECTOR_3D_MID);
  
  GridOpts pgdall;
  set_grid_opts_helper(pgdall, "pgdall", true, NULL, VECTOR_3D_MID);
  
  GridOpts icwu;
  set_grid_opts_helper(icwu, "icwu", true, NULL, VECTOR_3D_MID);
  

} 
#endif
