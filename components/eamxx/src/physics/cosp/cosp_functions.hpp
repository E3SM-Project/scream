#ifndef SCREAM_COSP_FUNCTIONS_HPP
#define SCREAM_COSP_FUNCTIONS_HPP
#include "share/scream_types.hpp"
using scream::Real;
extern "C" void cosp_c2f_init(int ncol, int nsubcol, int nlay);
extern "C" void cosp_c2f_final();
extern "C" void cosp_c2f_run(const int ncol, const int nsubcol, const int nlay, const int ntau, const int nctp, const int ncth, const int ndbze, const int nhgt,
    const Real emsfc_lw, const Real* sunlit, const Real* skt,
    const Real* T_mid, const Real* p_mid, const Real* p_int, const Real* z_mid, const Real* z_int, const Real* qv, const Real* qc, const Real* qi,
    const Real* cldfrac, const Real* reff_qc, const Real* reff_qi, const Real* dtau067, const Real* dtau105,
    Real* isccp_cldtot, Real* isccp_ctptau, Real* modis_ctptau, Real* misr_cthtau, Real* radar_dbze, Real* radar_cfad);

namespace scream {

    namespace CospFunc {
        using lview_host_1d = typename ekat::KokkosTypes<HostDevice>::template lview<Real*  >;
        using lview_host_2d = typename ekat::KokkosTypes<HostDevice>::template lview<Real** >;
        using lview_host_3d = typename ekat::KokkosTypes<HostDevice>::template lview<Real***>;
        template <typename S>
        using view_1d = typename ekat::KokkosTypes<HostDevice>::template view_1d<S>;
        template <typename S>
        using view_2d = typename ekat::KokkosTypes<HostDevice>::template view_2d<S>;
        template <typename S>
        using view_3d = typename ekat::KokkosTypes<HostDevice>::template view_3d<S>;

        inline void initialize(int ncol, int nsubcol, int nlay) {
            cosp_c2f_init(ncol, nsubcol, nlay);
        };
        inline void finalize() {
            cosp_c2f_final();
        };
        inline void main(
                const Int ncol, const Int nsubcol, const Int nlay, const Int ntau, const Int nctp, const Int ncth, const Int ndbze, const Int nhgt,
                const Real emsfc_lw,
                view_1d<const Real>& sunlit , view_1d<const Real>& skt,
                view_2d<const Real>& T_mid  , view_2d<const Real>& p_mid  , view_2d<const Real>& p_int,
                view_2d<const Real>& z_mid  , view_2d<const Real>& z_int  ,
                view_2d<const Real>& qv     , view_2d<const Real>& qc     , view_2d<const Real>& qi,
                view_2d<const Real>& cldfrac,
                view_2d<const Real>& reff_qc, view_2d<const Real>& reff_qi,
                view_2d<const Real>& dtau067, view_2d<const Real>& dtau105,
                view_1d<Real>& isccp_cldtot , view_3d<Real>& isccp_ctptau, view_3d<Real>& modis_ctptau, view_3d<Real>& misr_cthtau,
                view_3d<Real>& radar_dbze   , view_3d<Real>& radar_cfad) {

            // Make host copies and permute data as needed
            lview_host_2d
                  T_mid_h("T_mid_h", ncol, nlay), p_mid_h("p_mid_h", ncol, nlay), p_int_h("p_int_h", ncol, nlay+1),
                  z_mid_h("z_mid_h", ncol, nlay), z_int_h("z_int_h", ncol, nlay+1), 
                  qv_h("qv_h", ncol, nlay), qc_h("qc_h", ncol, nlay), qi_h("qi_h", ncol, nlay),
                  cldfrac_h("cldfrac_h", ncol, nlay),
                  reff_qc_h("reff_qc_h", ncol, nlay), reff_qi_h("reff_qi_h", ncol, nlay),
                  dtau067_h("dtau_067_h", ncol, nlay), dtau105_h("dtau105_h", ncol, nlay);
            lview_host_3d isccp_ctptau_h("isccp_ctptau_h", ncol, ntau, nctp);
            lview_host_3d modis_ctptau_h("modis_ctptau_h", ncol, ntau, nctp);
            lview_host_3d misr_cthtau_h("misr_cthtau_h", ncol, ntau, ncth);
            lview_host_3d radar_dbze_h("radar_dbze_h", ncol, nsubcol, nlay);
            lview_host_3d radar_cfad_h("radar_cfad_h", ncol, ndbze, nhgt);

            // Copy to layoutLeft host views
            for (int i = 0; i < ncol; i++) {
                for (int j = 0; j < nlay; j++) {
                    T_mid_h(i,j) = T_mid(i,j);
                    p_mid_h(i,j) = p_mid(i,j);
                    z_mid_h(i,j) = z_mid(i,j);
                    z_int_h(i,j) = z_int(i,j);
                    qv_h(i,j) = qv(i,j);
                    qc_h(i,j) = qc(i,j);
                    qi_h(i,j) = qi(i,j);
                    cldfrac_h(i,j) = cldfrac(i,j);
                    reff_qc_h(i,j) = reff_qc(i,j);
                    reff_qi_h(i,j) = reff_qi(i,j);
                    dtau067_h(i,j) = dtau067(i,j);
                    dtau105_h(i,j) = dtau105(i,j);
                }
            }
            for (int i = 0; i < ncol; i++) {
                for (int j = 0; j < nlay+1; j++) {
                    p_int_h(i,j) = p_int(i,j);
                }
            }

            // Subsample here?

            // Call COSP wrapper
            cosp_c2f_run(ncol, nsubcol, nlay, ntau, nctp, ncth, ndbze, nhgt,
                    emsfc_lw, sunlit.data(), skt.data(), T_mid_h.data(), p_mid_h.data(), p_int_h.data(),
                    z_mid_h.data(), z_int_h.data(), qv_h.data(), qc_h.data(), qi_h.data(),
                    cldfrac_h.data(), reff_qc_h.data(), reff_qi_h.data(), dtau067_h.data(), dtau105_h.data(),
                    isccp_cldtot.data(), isccp_ctptau_h.data(), modis_ctptau_h.data(), misr_cthtau_h.data(),
                    radar_dbze_h.data(), radar_cfad_h.data());

            // Copy outputs back to layoutRight views
            for (int i = 0; i < ncol; i++) {
                for (int j = 0; j < ntau; j++) {
                    for (int k = 0; k < nctp; k++) {
                        isccp_ctptau(i,j,k) = isccp_ctptau_h(i,j,k);
                        modis_ctptau(i,j,k) = modis_ctptau_h(i,j,k);
                    }
                    for (int k = 0; k < ncth; k++) {
                        misr_cthtau(i,j,k) = misr_cthtau_h(i,j,k);
                    }
                }
                for (int j = 0; j < nsubcol; j++) {
                    for (int k = 0; k < nlay; k++) {
                        radar_dbze(i,j,k) = radar_dbze_h(i,j,k);
                    }
                }
                for (int j = 0; j < ndbze; j++) {
                    for (int k = 0; k < nhgt; k++) {
                        radar_cfad(i,j,k) = radar_cfad_h(i,j,k);
                    }
                }
            }
        }
    }
}
#endif  /* SCREAM_COSP_FUNCTIONS_HPP */
