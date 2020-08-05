#ifndef ZM_GRID_OPTS_HPP
#define ZM_GRID_OPTS_HPP


#include <unordered_map>


#include <string>


using namespace std;

struct GridOpts{
  string name;
  bool isOut;
//auto unit
//FieldLayout field
};

unordered_map<string, GridOpts> opt_map;

void set_grid_opts_helper(GridOpts O, string n, bool out //auto unit, FieldLayout field
                          ){
  O.name = n;
  O.isOut = out;
  opt_map.insert({O.name, O});
}

void set_grid_opts(){
  	
  GridOpts t;
  set_grid_opts_helper(t, "t", true);

  GridOpts qh;
  set_grid_opts_helper(qh, "qh", true);

  GridOpts prec;
  set_grid_opts_helper(prec, "prec", true);

  GridOpts jctop;
  set_grid_opts_helper(jctop, "jctop", true);

  GridOpts jcbot;
  set_grid_opts_helper(jcbot, "jcbot", true);

  GridOpts pblh;
  set_grid_opts_helper(pblh, "pblh", true);

  GridOpts zm;
  set_grid_opts_helper(zm, "zm", true);

  GridOpts geos;
  set_grid_opts_helper(geos, "geos", true);

  GridOpts zi;
  set_grid_opts_helper(zi, "zi", true);

  GridOpts qtnd;
  set_grid_opts_helper(qtnd, "qtnd", true);

  GridOpts heat;
  set_grid_opts_helper(heat, "heat", true);

  GridOpts pap;
  set_grid_opts_helper(pap, "pap", true);

  GridOpts paph;
  set_grid_opts_helper(paph, "paph", true);

  GridOpts dpp;
  set_grid_opts_helper(dpp, "dpp", true);

  GridOpts mcon;
  set_grid_opts_helper(mcon, "mcon", true);

  GridOpts cme;
  set_grid_opts_helper(cme, "cme", true);

  GridOpts cape;
  set_grid_opts_helper(cape, "cape", true);

  GridOpts tpert;
  set_grid_opts_helper(tpert, "tpert", true);

  GridOpts dlf;
  set_grid_opts_helper(dlf, "dlf", true);

  GridOpts pflx;
  set_grid_opts_helper(pflx, "pflx", true);

  GridOpts zdu;
  set_grid_opts_helper(zdu, "zdu", true);

  GridOpts rprd;
  set_grid_opts_helper(rprd, "rprd", true);
  
  GridOpts mu;
  set_grid_opts_helper(mu, "mu", true);
  
  GridOpts md;
  set_grid_opts_helper(md, "md", true);
  
  GridOpts du;
  set_grid_opts_helper(du, "du", true);
  
  GridOpts eu;
  set_grid_opts_helper(eu, "eu", true);
  
  GridOpts ed;
  set_grid_opts_helper(ed, "ed", true);
  
  GridOpts dp;
  set_grid_opts_helper(dp, "dp", true);
  
  GridOpts dsubcld;
  set_grid_opts_helper(dsubcld, "dsubcld", true);
  
  GridOpts jt;
  set_grid_opts_helper(jt, "jt", true);
  
  GridOpts maxg;
  set_grid_opts_helper(maxg, "maxg", true);
  
  GridOpts ideep;
  set_grid_opts_helper(ideep, "ideep", true);
  
  GridOpts ql;
  set_grid_opts_helper(ql, "ql", true);
  
  GridOpts rliq;
  set_grid_opts_helper(rliq, "rliq", true);
  
  GridOpts landfrac;
  set_grid_opts_helper(landfrac, "landfrac", true);
  
  GridOpts hu_nm1;
  set_grid_opts_helper(hu_nm1, "hu_nm1", true);
  
  GridOpts cnv_nm1;
  set_grid_opts_helper(cnv_nm1, "cnv_nm1", true);
  
  GridOpts tm1;
  set_grid_opts_helper(tm1, "tm1", true);
  
  GridOpts qm1;
  set_grid_opts_helper(qm1, "qm1", true);
  
  GridOpts dcape;
  set_grid_opts_helper(dcape, "dcape", true);
  
  GridOpts q;
  set_grid_opts_helper(q, "q", true);
  
  GridOpts snow;
  set_grid_opts_helper(snow, "snow", true);
  
  GridOpts ntprprd;
  set_grid_opts_helper(ntprprd, "ntprprd", true);
  
  GridOpts ntsnprd;
  set_grid_opts_helper(ntsnprd, "ntsnprd", true);
  
  GridOpts pguall;
  set_grid_opts_helper(pguall, "pguall", true);
  
  GridOpts pgdall;
  set_grid_opts_helper(pgdall, "pgdall", true);
  
  GridOpts icwu;
  set_grid_opts_helper(icwu, "icwu", true);
  

} 
#endif
