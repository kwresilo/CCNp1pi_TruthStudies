// Analysis macro for use in the CCNp0pi single transverse variable analysis
// Designed for use with the PeLEE group's "searchingfornues" ntuples
//
// Updated 22 April 2023
// Steven Gardiner <gardiner@fnal.gov>

// Standard library includes
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

// ROOT includes
#include "TChain.h"
#include "TFile.h"
#include "TParameter.h"
#include "TTree.h"
#include "TVector3.h"

// STV analysis includes
#include "EventCategory.hh"
#include "FiducialVolume.hh"
#include "TreeUtils.hh"

// Helper function that avoids NaNs when taking square roots of negative
// numbers
double real_sqrt( double x ) {
  if ( x < 0. ) return 0.;
  else return std::sqrt( x );
}

// Helper function that returns true if a given PDG code represents a meson or
// antimeson. Otherwise returns false. Based on points 10, 12, and 13 of the
// Particle Data Group's "Monte Carlo Particle Numbering Scheme"
// (2019 revision).
bool is_meson_or_antimeson( int pdg_code ) {

  // Ignore differences between mesons and antimesons for this test. Mesons
  // will have positive PDG codes, while antimesons will have negative ones.
  int abs_pdg = std::abs( pdg_code );

  // make exception for charged pions
  // we want to return false to satisfy defintion of mc_is_signal_
  if (abs_pdg == 211) return false;

  // Meson PDG codes have no more than seven digits. Seven-digit
  // codes beginning with "99" are reserved for generator-specific
  // particles
  if ( abs_pdg >= 9900000 ) return false;

  // Mesons have a value of zero for $n_{q1}$, the thousands digit
  int thousands_digit = ( abs_pdg / 1000 ) % 10;
  if ( thousands_digit != 0 ) return false;

  // They also have a nonzero value for $n_{q2}$, the hundreds digit
  int hundreds_digit = ( abs_pdg / 100 ) % 10;
  if ( hundreds_digit == 0 ) return false;

  // Reserved codes for Standard Model parton distribution functions
  if ( abs_pdg >= 901 && abs_pdg <= 930 ) return false;

  // Reggeon and pomeron
  if ( abs_pdg == 110 || abs_pdg == 990 ) return false;

  // Reserved codes for GEANT tracking purposes
  if ( abs_pdg == 998 || abs_pdg == 999 ) return false;

  // Reserved code for generator-specific pseudoparticles
  if ( abs_pdg == 100 ) return false;

  // If we've passed all of the tests above, then the particle is a meson
  return true;
}

// constants needed for multi-plane Bragg calculations
const float floatMinValue = -340282346638528859811704183484516925440.000000f;
const auto sin2AngleThreshold = 0.175f;
const auto piBy3 = std::acos(0.5f);
const auto lengthFraction = 0.333333333f;
const auto nHitsToSkip = 3;

// pion momentum calibration
constexpr float a = 1.6353;
constexpr float b = 0.0022385;
constexpr float c = 1.5668;
constexpr float d = 0.013327;


// A few helpful dummy constants
constexpr float BOGUS = 9999.;
constexpr int BOGUS_INT = 9999;
constexpr int BOGUS_INDEX = -1;
constexpr float LOW_FLOAT = -1e30;
constexpr float DEFAULT_WEIGHT = 1.;

// Integer representation of CC versus NC for the ccnc branch
constexpr int CHARGED_CURRENT = 0;
constexpr int NEUTRAL_CURRENT = 1;

// Useful PDG codes
constexpr int ELECTRON_NEUTRINO = 12;
constexpr int MUON = 13;
constexpr int MUON_NEUTRINO = 14;
constexpr int TAU_NEUTRINO = 16;
constexpr int PROTON = 2212;
constexpr int PI_ZERO = 111;
constexpr int PI_PLUS = 211;
constexpr int NEUTRON = 2112;
// Values of parameters to use in analysis cuts
constexpr float DEFAULT_PROTON_PID_CUT = 0.2;
constexpr float LEAD_P_MIN_MOM_CUT = 0.300; // GeV/c
constexpr float LEAD_P_MAX_MOM_CUT = 1.; // GeV/c
constexpr float MUON_P_MIN_MOM_CUT = 0.150; // GeV/c
constexpr float MUON_P_MAX_MOM_CUT = 1.200; // GeV/c
constexpr float CHARGED_PI_MOM_CUT = 0.05; // GeV/c
constexpr float MUON_MOM_QUALITY_CUT = 0.25; // fractional difference

constexpr float TOPO_SCORE_CUT = 0.1;
constexpr float COSMIC_IP_CUT = 10.; // cm

constexpr float MUON_TRACK_SCORE_CUT = 0.8;
constexpr float MUON_VTX_DISTANCE_CUT = 4.; // cm
constexpr float MUON_LENGTH_CUT = 10.; // cm
constexpr float MUON_PID_CUT = 0.2;

constexpr float TRACK_SCORE_CUT = 0.5;

constexpr float PROTON_BDT_CUT = 0.01;
constexpr float MUON_BDT_CUT = 0.0;
constexpr float TOPO_SCORE_CUT_2 = 0.67; 
constexpr float PFP_DISTANCE_CUT = 9.5; 
// Function that defines the track-length-dependent proton PID cut
double proton_pid_cut( double track_length ) {

  double cut = DEFAULT_PROTON_PID_CUT;

  // Piecewise cut removed 27 June 2021
  //// All track length values are in cm
  //if ( track_length >= 0. && track_length <= 10.5 ) {
  //  cut = -0.0034219305*std::pow( track_length, 2 );
  //  cut += 0.018436866*track_length + 0.062718401;
  //}
  //else if ( track_length > 10.5 && track_length <= 33.1776508 ) {
  //  cut = 0.014153245*( track_length - 10.5 ) - 0.12096235;
  //}

  return cut;

}

// Boundaries of the proton containment volume (used in reco only) in cm
double PCV_X_MIN =   10.;
double PCV_X_MAX =  246.35;

double PCV_Y_MIN = -106.5;
double PCV_Y_MAX =  106.5;

double PCV_Z_MIN =   10.;
double PCV_Z_MAX = 1026.8;

// Mass values from GENIE v3.0.6
constexpr double TARGET_MASS = 37.215526; // 40Ar, GeV
constexpr double NEUTRON_MASS = 0.93956541; // GeV
constexpr double PROTON_MASS = 0.93827208; // GeV
constexpr double MUON_MASS = 0.10565837; // GeV
constexpr double PI_PLUS_MASS = 0.13957000; // GeV

// This binding energy value is used in GENIE v3.0.6
//constexpr double BINDING_ENERGY = 0.0295; // 40Ar, GeV

// This value is the shell-occupancy-weighted mean of the $E_{\alpha}$ values
// listed for 40Ar in Table II of arXiv:1609.03530. MINERvA uses an identical
// procedure for 12C to obtain the binding energy value of 27.13 MeV, which is
// adopted in their STV analysis described in arXiv:1910.08658.
constexpr double BINDING_ENERGY = 0.34381; // 40Ar, GeV
constexpr double MEAN_EXCITATION_ENERGY = 0.0309; 

// Class used to hold information from the searchingfornues TTree branches and
// process it for our analysis
class AnalysisEvent {

  public:

    AnalysisEvent() {}

    ~AnalysisEvent() {}

    EventCategory categorize_event();
    void apply_selection();
    void get_multi_plane_Bragg_likelihood();
    void apply_numu_CC_selection();
    void find_pion_candidate();
    void find_muon_candidate(); 
    void find_lead_p_candidate();
    void compute_observables();
    void compute_mc_truth_observables();

    // Event scores needed for numu CC selection
    float topological_score_ = BOGUS;
    float cosmic_impact_parameter_ = BOGUS;

    // Backtracked purity and completeness of hits (MC only)
    float nu_completeness_from_pfp_ = BOGUS;
    float nu_purity_from_pfp_ = BOGUS;

      
    // Reco PDG code of the neutrino candidate
    int nu_pdg_ = BOGUS_INT;

    // Number of neutrino slices identified by the SliceID. Allowed values
    // are zero or one.
    int nslice_ = BOGUS_INT;

    // Reco neutrino vertex coordinates (cm). Space charge corrections have
    // been applied for these.
    float nu_vx_ = BOGUS;
    float nu_vy_ = BOGUS;
    float nu_vz_ = BOGUS;

    // Reconstructed object counts
    int num_pf_particles_ = BOGUS_INT;
    int num_tracks_ = BOGUS_INT;
    int num_showers_ = BOGUS_INT;

    // PFParticle properties
    MyPointer< std::vector<unsigned int> > pfp_generation_;
    MyPointer< std::vector<unsigned int> > pfp_trk_daughters_count_;
    MyPointer< std::vector<unsigned int> > pfp_shr_daughters_count_;

    MyPointer< std::vector<float> > pfp_track_score_;

    // Reco PDG code assigned by Pandora
    MyPointer< std::vector<int> > pfp_reco_pdg_;

    // Total number of wire plane hits associated with each PFParticle
    MyPointer< std::vector<int> > pfp_hits_;

    // Number of hits on the three individual planes
    // (Y is the collection plane)
    MyPointer< std::vector<int> > pfp_hitsU_;
    MyPointer< std::vector<int> > pfp_hitsV_;
    MyPointer< std::vector<int> > pfp_hitsY_;

    // True PDG code found using the backtracker
    MyPointer< std::vector<int> > pfp_true_pdg_;

    // True 4-momentum components found using the backtracker
    MyPointer< std::vector<float> > pfp_true_E_;
    MyPointer< std::vector<float> > pfp_true_px_;
    MyPointer< std::vector<float> > pfp_true_py_;
    MyPointer< std::vector<float> > pfp_true_pz_;

    // Shower properties
    MyPointer< std::vector<unsigned long> > shower_pfp_id_;
    MyPointer< std::vector<float> > shower_startx_;
    MyPointer< std::vector<float> > shower_starty_;
    MyPointer< std::vector<float> > shower_startz_;
    MyPointer< std::vector<float> > shower_start_distance_;

    // Track properties
    MyPointer< std::vector<unsigned long> > track_pfp_id_;
    MyPointer< std::vector<float> > track_length_;
    MyPointer< std::vector<float> > track_startx_;
    MyPointer< std::vector<float> > track_starty_;
    MyPointer< std::vector<float> > track_startz_;
    MyPointer< std::vector<float> > track_start_distance_;
    MyPointer< std::vector<float> > track_endx_;
    MyPointer< std::vector<float> > track_endy_;
    MyPointer< std::vector<float> > track_endz_;
    MyPointer< std::vector<float> > track_dirx_;
    MyPointer< std::vector<float> > track_diry_;
    MyPointer< std::vector<float> > track_dirz_;

    // Proton *kinetic* energy using range-based momentum reconstruction
    MyPointer< std::vector<float> > track_kinetic_energy_p_;

    MyPointer< std::vector<float> > track_range_mom_mu_;
    MyPointer< std::vector<float> > track_mcs_mom_mu_;
    MyPointer< std::vector<float> > track_chi2_proton_;

    // BDT scores
    MyPointer< std::vector<float> > proton_BDT_score_;
    MyPointer< std::vector<float> > muon_BDT_score_;
    // Log-likelihood ratio particle ID information

    // Product of muon/proton log-likelihood ratios from all wire three planes
    MyPointer< std::vector<float> > track_llr_pid_;

    // Individual wire plane muon/proton log-likelihood ratios
    MyPointer< std::vector<float> > track_llr_pid_U_;
    MyPointer< std::vector<float> > track_llr_pid_V_;
    MyPointer< std::vector<float> > track_llr_pid_Y_;

    // Rescaled overall PID score (all three planes) that lies
    // on the interval [-1, 1]
    MyPointer< std::vector<float> > track_llr_pid_score_;

    // Golden pion handles
    MyPointer< std::vector<float> > mc_end_p_;
    MyPointer< std::vector<int> > mc_n_inelastic_;
    MyPointer< std::vector<int> > mc_n_elastic_;
    
    MyPointer< std::vector<int> > pfp_n_descendents_;
    MyPointer< std::vector<int> > trk_end_spacepoints_;
    MyPointer< std::vector<float> > trk_avg_deflection_stdev_;
   
    // Calo features, 3-plane Bragg

    MyPointer< std::vector<float> > bragg_mip_uvw_;

    MyPointer< std::vector<float> > bragg_p_fwd_uvw_;
    MyPointer< std::vector<float> > bragg_p_bwd_uvw_;
    MyPointer< std::vector<float> > bragg_p_to_MIP_;
    MyPointer< std::vector<float> > bragg_p_fwd_2_bwd_;

    MyPointer< std::vector<float> > bragg_pion_fwd_uvw_;
    MyPointer< std::vector<float> > bragg_pion_bwd_uvw_;
    MyPointer< std::vector<float> > bragg_pion_to_MIP_;
    MyPointer< std::vector<float> > bragg_pion_fwd_2_bwd_;

    MyPointer< std::vector<float> > bragg_mu_fwd_uvw_;
    MyPointer< std::vector<float> > bragg_mu_bwd_uvw_;

    MyPointer< std::vector<float> > truncated_mean_dEdx_;

    // Calo feature single plane Bragg
    
    MyPointer< std::vector<float> > trk_trunk_dEdx_w_;
    MyPointer< std::vector<float> > trk_trunk_dEdx_u_;
    MyPointer< std::vector<float> > trk_trunk_dEdx_v_;

    MyPointer< std::vector<float> > bragg_p_fwd_w_;
    MyPointer< std::vector<float> > bragg_p_bwd_w_;
    
    MyPointer< std::vector<float> > bragg_mu_fwd_w_;
    MyPointer< std::vector<float> > bragg_mu_bwd_w_;

    MyPointer< std::vector<float> > bragg_pion_fwd_w_;
    MyPointer< std::vector<float> > bragg_pion_bwd_w_;

    MyPointer< std::vector<float> > bragg_p_fwd_v_;
    MyPointer< std::vector<float> > bragg_p_bwd_v_;

    MyPointer< std::vector<float> > bragg_mu_fwd_v_;
    MyPointer< std::vector<float> > bragg_mu_bwd_v_;

    MyPointer< std::vector<float> > bragg_pion_fwd_v_;
    MyPointer< std::vector<float> > bragg_pion_bwd_v_;

    MyPointer< std::vector<float> > bragg_p_fwd_u_;
    MyPointer< std::vector<float> > bragg_p_bwd_u_;

    MyPointer< std::vector<float> > bragg_mu_fwd_u_;
    MyPointer< std::vector<float> > bragg_mu_bwd_u_;

    MyPointer< std::vector<float> > bragg_pion_fwd_u_;
    MyPointer< std::vector<float> > bragg_pion_bwd_u_;

    MyPointer< std::vector<bool> > trk_bragg_p_fwd_preferred_w_;
    MyPointer< std::vector<float> > trk_bragg_p_w_;
    MyPointer< std::vector<float> > trk_bragg_p_alt_dir_w_;
    
    MyPointer< std::vector<bool> > trk_bragg_pion_fwd_preferred_w_;
    MyPointer< std::vector<float> > trk_bragg_pion_w_;
    MyPointer< std::vector<float> > trk_bragg_pion_alt_dir_w_;

    MyPointer< std::vector<bool> > trk_bragg_mu_fwd_preferred_w_;
    MyPointer< std::vector<float> > trk_bragg_mu_w_;
    MyPointer< std::vector<float> > trk_bragg_mu_alt_dir_w_;

    MyPointer< std::vector<bool> > trk_bragg_p_fwd_preferred_u_;
    MyPointer< std::vector<float> > trk_bragg_p_u_;
    MyPointer< std::vector<float> > trk_bragg_p_alt_dir_u_;

    MyPointer< std::vector<bool> > trk_bragg_pion_fwd_preferred_u_;
    MyPointer< std::vector<float> > trk_bragg_pion_u_;
    MyPointer< std::vector<float> > trk_bragg_pion_alt_dir_u_;

    MyPointer< std::vector<bool> > trk_bragg_mu_fwd_preferred_u_;
    MyPointer< std::vector<float> > trk_bragg_mu_u_;
    MyPointer< std::vector<float> > trk_bragg_mu_alt_dir_u_;
 
    MyPointer< std::vector<bool> > trk_bragg_p_fwd_preferred_v_;
    MyPointer< std::vector<float> > trk_bragg_p_v_;
    MyPointer< std::vector<float> > trk_bragg_p_alt_dir_v_;

    MyPointer< std::vector<bool> > trk_bragg_pion_fwd_preferred_v_;
    MyPointer< std::vector<float> > trk_bragg_pion_v_;
    MyPointer< std::vector<float> > trk_bragg_pion_alt_dir_v_;

    MyPointer< std::vector<bool> > trk_bragg_mu_fwd_preferred_v_;
    MyPointer< std::vector<float> > trk_bragg_mu_v_;
    MyPointer< std::vector<float> > trk_bragg_mu_alt_dir_v_;

    MyPointer< std::vector<float> > bragg_mip_w_;
    MyPointer< std::vector<float> > bragg_mip_v_;
    MyPointer< std::vector<float> > bragg_mip_u_;
    // True neutrino PDG code
    int mc_nu_pdg_ = BOGUS_INT;

    // True neutrino vertex coordinates (cm)
    float mc_nu_vx_ = BOGUS;
    float mc_nu_vy_ = BOGUS;
    float mc_nu_vz_ = BOGUS;

    // True neutrino 4-momentum
    float mc_nu_energy_ = BOGUS;

    // Whether the event is CC (0) or NC (1)
    int mc_nu_ccnc_ = false;

    // Interaction mode (QE, MEC, etc.)
    int mc_nu_interaction_type_ = BOGUS_INT;

    // Final-state particle PDG codes and energies (post-FSIs)
    MyPointer< std::vector<int> > mc_nu_daughter_pdg_;
    MyPointer< std::vector<float> > mc_nu_daughter_energy_;
    MyPointer< std::vector<float> > mc_nu_daughter_px_;
    MyPointer< std::vector<float> > mc_nu_daughter_py_;
    MyPointer< std::vector<float> > mc_nu_daughter_pz_;

    // General systematic weights
    MyPointer< std::map< std::string, std::vector<double> > > mc_weights_map_;
    // Map of pointers used to set output branch addresses for the elements
    // of the weights map. Hacky, but it works.
    // TODO: revisit this to make something more elegant
    std::map< std::string, std::vector<double>* > mc_weights_ptr_map_;

    // GENIE weights
    float spline_weight_ = DEFAULT_WEIGHT;
    float tuned_cv_weight_ = DEFAULT_WEIGHT;

    bool mc_golden_ = false;
    bool mc_pi_stopping_ = false;
    // Signal definition requirements
    bool is_mc_ = false;
    bool mc_neutrino_is_numu_ = false;
    bool mc_vertex_in_FV_ = false;
    bool mc_muon_in_mom_range_ = false;
    bool mc_lead_p_in_mom_range_ = false;
    bool mc_no_fs_mesons_ = false;
    bool mc_1cpi_= false;
    bool mc_no_neutrons_ = true;  
    // Intersection of all of these requirements
    bool mc_is_signal_ = false;

    // Extra flags for looking specifically at final-state pions
    bool mc_no_fs_pi0_ = false;
    bool mc_no_charged_pi_above_threshold_ = false;
	
    EventCategory category_ = kUnknown;

    // **** Reco selection requirements ****

    // Whether the event passed the numu CC selection (a subset of the cuts
    // used for the full analysis)
    bool sel_nu_mu_cc_ = false;

    // Whether the reconstructed neutrino vertex lies within the fiducial
    // volume
    bool sel_reco_vertex_in_FV_ = false;
    // Whether the event passed the topological score cut
    bool sel_topo_cut_passed_ = false;
    // Whether the event passed the cosmic impact parameter cut
    bool sel_cosmic_ip_cut_passed_ = false;
    // Whether the start points for all PFParticles lie within the
    // proton containment volume
    bool sel_pfp_starts_in_PCV_ = false;

    // Whether Pandora assigns PDG score that of a muon netrino
    bool sel_pandoraPDG_isnumu_ = false;

    // True if a generation == 2 muon candidate was identified
    bool sel_has_muon_candidate_ = false;

    // Whether all pfps are contained
    bool sel_all_pfp_contained_ = true;

    // Whether the end point of the muon candidate track is contained
    // in the "containment volume"
    bool sel_muon_contained_ = false;

    // Whether the muon candidate has MCS- and range-based reco momenta
    // that agree within a given tolerance
    bool sel_muon_quality_ok_ = false;

    // Whether the muon candidate has a reco momentum above threshold
    bool sel_muon_passed_mom_cuts_ = false;

    // False if at least one generation == 2 shower was reconstructed
    bool sel_no_reco_showers_ = false;

    // False if less than 3 tracks reconstructed
    bool sel_min_3_tracks_ = false;

    // False if number of non proton-like particles is not exactly two
    bool sel_2_non_proton_ = false;
   
    // All pfp must start within certain distance of reco vtx
    bool sel_all_pfp_in_vtx_proximity_ = true; 
   
    // Whether it has pion candidate
    bool sel_has_pion_candidate_ = false;
 
    // Whether at least one generation == 2 reco track exists that is not the
    // muon candidate
    bool sel_has_p_candidate_ = false;
    // Whether all proton candidates (i.e., all tracks which are not the muon
    // candidate) pass the proton PID cut or not
    bool sel_passed_proton_pid_cut_ = false;
    // Whether all proton candidates have track end coordinates that lie within
    // the "containment volume"
    bool sel_protons_contained_ = false;
    // Whether the leading proton candidate has a range-based reco momentum
    // above LEAD_P_MIN_MOM_CUT and below LEAD_P_MAX_MOM_CUT
    bool sel_lead_p_passed_mom_cuts_ = false;
    // Intersection of all of the above requirements
    bool sel_CCNp1pi_ = false;

    // Muon and leading proton candidate indices (BOGUS_INDEX if not present)
    // in the reco track arrays
    int muon_candidate_pid_idx_ = BOGUS_INDEX;
    int muon_candidate_length_idx_ = BOGUS_INDEX;
    int muon_candidate_bdt_idx_ = BOGUS_INDEX;
    int pion_candidate_idx_ = BOGUS_INDEX;
    int lead_p_candidate_idx_ = BOGUS_INDEX;

    // ** Reconstructed observables **

    // n non proton like
    int n_non_proton_like_ = BOGUS_INT;
    int n_reco_tracks_ = BOGUS_INT;

    // 3-momenta
    MyPointer< TVector3 > p3_mu_;
    MyPointer< TVector3 > p3_lead_p_;
    MyPointer< TVector3 > p3_cpi_;
    // Reconstructed 3-momenta for all proton candidates,
    // ordered from highest to lowest by magnitude
    MyPointer< std::vector<TVector3> > p3_p_vec_;

    // Reco STVs and other variables of interest
   
    // Longitudinal particle momenta
    float muE_ = BOGUS;
    float pE_ = BOGUS;
    float piE_ = BOGUS;
 
    float reco_Ecal_ = BOGUS;
    float pKE_ = BOGUS;

    MyPointer< TVector3 > q_ ; 

    float delta_pT_ = BOGUS;
    //float delta_phiT_ = BOGUS;
    //float delta_alphaT_ = BOGUS;
    float delta_pL_ = BOGUS;
    float delta_pL2_ = BOGUS;
    float pn_ = BOGUS;
    float pn2_ = BOGUS;
    //float delta_pTx_ = BOGUS;
    //float delta_pTy_ = BOGUS;
    float theta_mu_p_ = BOGUS;
    float theta_mu_cpi_ = BOGUS;
    
    float delta_alpha3D_ = BOGUS;
    float delta_phi3d_had_ = BOGUS;
    float delta_phi3d_mu_ = BOGUS;
    float delta_phi3d_p_ = BOGUS;
    float delta_phi3d_pi_ = BOGUS;
    // ** MC truth observables **
    // These are loaded for signal events whenever we have MC information
    // to use

    // 3-momenta
    MyPointer< TVector3 > mc_p3_mu_;
    MyPointer< TVector3 > mc_p3_lead_p_;
    MyPointer< TVector3 > mc_p3_cpi_;

    // True 3-momenta for all true MC protons and charged pions, ordered from highest to lowest
    // by magnitude
    MyPointer< std::vector<TVector3> > mc_p3_p_vec_;
    int mc_n_protons_ = 0;

    MyPointer< std::vector<TVector3> > mc_p3_cpi_vec_;
    int mc_n_cpi_ = 0; 
    // MC truth STVs and other variables of interest
    float mc_muE_ = BOGUS;
    float mc_pE_ = BOGUS;
    float mc_piE_ = BOGUS;

    float mc_Ecal_ = BOGUS;
    float mc_pKE_ = BOGUS;

    MyPointer< TVector3 > mc_q_;

    float mc_delta_pT_ = BOGUS;
    //float mc_delta_phiT_ = BOGUS;
    //float mc_delta_alphaT_ = BOGUS;
    float mc_delta_pL_ = BOGUS;
    float mc_delta_pL2_ = BOGUS;
    float mc_pn_ = BOGUS;
    float mc_pn2_ = BOGUS;
    //float mc_delta_pTx_ = BOGUS;
    //float mc_delta_pTy_ = BOGUS;
    float mc_theta_mu_p_ = BOGUS;
    float mc_theta_mu_cpi_ = BOGUS;

    float mc_delta_alpha3D_ = BOGUS;
    float mc_delta_phi3d_had_ = BOGUS;
    float mc_delta_phi3d_mu_ = BOGUS;
    float mc_delta_phi3d_p_ = BOGUS;
    float mc_delta_phi3d_pi_ = BOGUS;

    bool reco_vertex_inside_FV() {
      return point_inside_FV( nu_vx_, nu_vy_, nu_vz_ );
    }

    bool mc_vertex_inside_FV() {
      return point_inside_FV( mc_nu_vx_, mc_nu_vy_, mc_nu_vz_ );
    }

    bool in_proton_containment_vol( float x, float y, float z ) {
      bool x_inside_PCV = ( PCV_X_MIN < x ) && ( x < PCV_X_MAX );
      bool y_inside_PCV = ( PCV_Y_MIN < y ) && ( y < PCV_Y_MAX );
      bool z_inside_PCV = ( PCV_Z_MIN < z ) && ( z < PCV_Z_MAX );
      return ( x_inside_PCV && y_inside_PCV && z_inside_PCV );
    }

};

// Helper function to set branch addresses for reading information
// from the Event TTree
void set_event_branch_addresses(TTree& etree, AnalysisEvent& ev)
{
  // variables for identifying golden pion
  set_object_input_branch_address( etree, "mc_end_p", ev.mc_end_p_ );
  set_object_input_branch_address( etree, "mc_n_elastic", ev.mc_n_elastic_);
  set_object_input_branch_address( etree, "mc_n_inelastic", ev.mc_n_inelastic_);

  set_object_input_branch_address( etree, "pfp_n_descendents_v", ev.pfp_n_descendents_);
  set_object_input_branch_address( etree, "trk_end_spacepoints_v", ev.trk_end_spacepoints_);
  set_object_input_branch_address( etree, "trk_avg_deflection_stdev_v", ev.trk_avg_deflection_stdev_);

  set_object_input_branch_address( etree, "trk_trunk_rr_dEdx_y_v", ev.trk_trunk_dEdx_w_);
  set_object_input_branch_address( etree, "trk_trunk_rr_dEdx_u_v", ev.trk_trunk_dEdx_u_);
  set_object_input_branch_address( etree, "trk_trunk_rr_dEdx_v_v", ev.trk_trunk_dEdx_v_);

  set_object_input_branch_address( etree, "trk_bragg_p_fwd_preferred_v", ev.trk_bragg_p_fwd_preferred_w_);
  set_object_input_branch_address( etree, "trk_bragg_p_v", ev.trk_bragg_p_w_);
  set_object_input_branch_address( etree, "trk_bragg_p_alt_dir_v", ev.trk_bragg_p_alt_dir_w_);
  
  set_object_input_branch_address( etree, "trk_bragg_mu_fwd_preferred_v", ev.trk_bragg_mu_fwd_preferred_w_);
  set_object_input_branch_address( etree, "trk_bragg_mu_v", ev.trk_bragg_mu_w_);
  set_object_input_branch_address( etree, "trk_bragg_mu_alt_dir_v", ev.trk_bragg_mu_alt_dir_w_);

  set_object_input_branch_address( etree, "trk_bragg_pion_fwd_preferred_v", ev.trk_bragg_pion_fwd_preferred_w_);
  set_object_input_branch_address( etree, "trk_bragg_pion_v", ev.trk_bragg_pion_w_);
  set_object_input_branch_address( etree, "trk_bragg_pion_alt_dir_v", ev.trk_bragg_pion_alt_dir_w_);

  set_object_input_branch_address( etree, "trk_bragg_p_fwd_preferred_u_v", ev.trk_bragg_p_fwd_preferred_u_);
  set_object_input_branch_address( etree, "trk_bragg_p_u_v", ev.trk_bragg_p_u_);
  set_object_input_branch_address( etree, "trk_bragg_p_alt_dir_u_v", ev.trk_bragg_p_alt_dir_u_);

  set_object_input_branch_address( etree, "trk_bragg_mu_fwd_preferred_u_v", ev.trk_bragg_mu_fwd_preferred_u_);
  set_object_input_branch_address( etree, "trk_bragg_mu_u_v", ev.trk_bragg_mu_u_);
  set_object_input_branch_address( etree, "trk_bragg_mu_alt_dir_u_v", ev.trk_bragg_mu_alt_dir_u_);

  set_object_input_branch_address( etree, "trk_bragg_pion_fwd_preferred_u_v", ev.trk_bragg_pion_fwd_preferred_u_);
  set_object_input_branch_address( etree, "trk_bragg_pion_u_v", ev.trk_bragg_pion_u_);
  set_object_input_branch_address( etree, "trk_bragg_pion_alt_dir_u_v", ev.trk_bragg_pion_alt_dir_u_);

  set_object_input_branch_address( etree, "trk_bragg_p_fwd_preferred_v_v", ev.trk_bragg_p_fwd_preferred_v_);
  set_object_input_branch_address( etree, "trk_bragg_p_v_v", ev.trk_bragg_p_v_);
  set_object_input_branch_address( etree, "trk_bragg_p_alt_dir_v_v", ev.trk_bragg_p_alt_dir_v_);

  set_object_input_branch_address( etree, "trk_bragg_mu_fwd_preferred_v_v", ev.trk_bragg_mu_fwd_preferred_v_);
  set_object_input_branch_address( etree, "trk_bragg_mu_v_v", ev.trk_bragg_mu_v_);
  set_object_input_branch_address( etree, "trk_bragg_mu_alt_dir_v_v", ev.trk_bragg_mu_alt_dir_v_);

  set_object_input_branch_address( etree, "trk_bragg_pion_fwd_preferred_v_v", ev.trk_bragg_pion_fwd_preferred_v_);
  set_object_input_branch_address( etree, "trk_bragg_pion_v_v", ev.trk_bragg_pion_v_);
  set_object_input_branch_address( etree, "trk_bragg_pion_alt_dir_v_v", ev.trk_bragg_pion_alt_dir_v_);
  
  set_object_input_branch_address( etree, "trk_bragg_mip_v", ev.bragg_mip_w_);
  set_object_input_branch_address( etree, "trk_bragg_mip_v_v", ev.bragg_mip_v_);
  set_object_input_branch_address( etree, "trk_bragg_mip_u_v", ev.bragg_mip_u_); 

  // Reco PDG code of primary PFParticle in slice (i.e., the neutrino
  // candidate)
  etree.SetBranchAddress( "slpdg", &ev.nu_pdg_ );

  // Number of neutrino slices identified by the SliceID. Allowed values
  // are zero or one.
  etree.SetBranchAddress( "nslice", &ev.nslice_ );

  // Topological score
  etree.SetBranchAddress( "topological_score", &ev.topological_score_ );
  etree.SetBranchAddress( "CosmicIP", &ev.cosmic_impact_parameter_ );
  //etree.SetBranchAddress( "CosmicIPAll3D", &ev.CosmicIPAll3D_ );

  // Reconstructed neutrino vertex position (with corrections for
  // space charge applied)
  etree.SetBranchAddress( "reco_nu_vtx_sce_x", &ev.nu_vx_ );
  etree.SetBranchAddress( "reco_nu_vtx_sce_y", &ev.nu_vy_ );
  etree.SetBranchAddress( "reco_nu_vtx_sce_z", &ev.nu_vz_ );

  // Reconstructed object counts
  etree.SetBranchAddress( "n_pfps", &ev.num_pf_particles_ );
  etree.SetBranchAddress( "n_tracks", &ev.num_tracks_ );
  etree.SetBranchAddress( "n_showers", &ev.num_showers_ );

  // PFParticle properties
  set_object_input_branch_address( etree, "pfp_generation_v",
    ev.pfp_generation_ );

  set_object_input_branch_address( etree, "pfp_trk_daughters_v",
    ev.pfp_trk_daughters_count_ );

  set_object_input_branch_address( etree, "pfp_shr_daughters_v",
    ev.pfp_shr_daughters_count_ );

  set_object_input_branch_address( etree, "trk_score_v", ev.pfp_track_score_ );
  set_object_input_branch_address( etree, "pfpdg", ev.pfp_reco_pdg_ );
  set_object_input_branch_address( etree, "pfnhits", ev.pfp_hits_ );
  set_object_input_branch_address( etree, "pfnplanehits_U", ev.pfp_hitsU_ );
  set_object_input_branch_address( etree, "pfnplanehits_V", ev.pfp_hitsV_ );
  set_object_input_branch_address( etree, "pfnplanehits_Y", ev.pfp_hitsY_ );

  // Backtracked PFParticle properties
  set_object_input_branch_address( etree, "backtracked_pdg", ev.pfp_true_pdg_ );
  set_object_input_branch_address( etree, "backtracked_e", ev.pfp_true_E_ );
  set_object_input_branch_address( etree, "backtracked_px", ev.pfp_true_px_ );
  set_object_input_branch_address( etree, "backtracked_py", ev.pfp_true_py_ );
  set_object_input_branch_address( etree, "backtracked_pz", ev.pfp_true_pz_ );

  // Shower properties
  // These are excluded from some ntuples to ensure blindness for the LEE
  // analyses. We will skip them when not available.
  bool has_shower_branches = ( etree.GetBranch("shr_pfp_id_v") != nullptr );
  if ( has_shower_branches ) {
    set_object_input_branch_address( etree, "shr_pfp_id_v", ev.shower_pfp_id_ );
    set_object_input_branch_address( etree, "shr_start_x_v", ev.shower_startx_ );
    set_object_input_branch_address( etree, "shr_start_y_v", ev.shower_starty_ );
    set_object_input_branch_address( etree, "shr_start_z_v", ev.shower_startz_ );
    // Shower start distance from reco neutrino vertex (pre-calculated for
    // convenience)
    set_object_input_branch_address( etree, "shr_dist_v",
      ev.shower_start_distance_ );
  }
  else {
    // When the shower information is not available, delete the owned vectors
    // to signal that the associated branches should not be written to the
    // output TTree
    ev.shower_pfp_id_.reset( nullptr );
    ev.shower_startx_.reset( nullptr );
    ev.shower_starty_.reset( nullptr );
    ev.shower_startz_.reset( nullptr );
    ev.shower_start_distance_.reset( nullptr );
  }

  // Track properties
  set_object_input_branch_address( etree, "trk_pfp_id_v", ev.track_pfp_id_ );
  set_object_input_branch_address( etree, "trk_len_v", ev.track_length_ );
  set_object_input_branch_address( etree, "trk_sce_start_x_v", ev.track_startx_ );
  set_object_input_branch_address( etree, "trk_sce_start_y_v", ev.track_starty_ );
  set_object_input_branch_address( etree, "trk_sce_start_z_v", ev.track_startz_ );

  // Track start distance from reco neutrino vertex (pre-calculated for
  // convenience)
  set_object_input_branch_address( etree, "trk_distance_v",
    ev.track_start_distance_ );

  set_object_input_branch_address( etree, "trk_sce_end_x_v", ev.track_endx_ );
  set_object_input_branch_address( etree, "trk_sce_end_y_v", ev.track_endy_ );
  set_object_input_branch_address( etree, "trk_sce_end_z_v", ev.track_endz_ );

  set_object_input_branch_address( etree, "trk_dir_x_v", ev.track_dirx_ );
  set_object_input_branch_address( etree, "trk_dir_y_v", ev.track_diry_ );
  set_object_input_branch_address( etree, "trk_dir_z_v", ev.track_dirz_ );

  set_object_input_branch_address( etree, "trk_energy_proton_v",
    ev.track_kinetic_energy_p_ );

  set_object_input_branch_address( etree, "trk_range_muon_mom_v",
    ev.track_range_mom_mu_ );

  set_object_input_branch_address( etree, "trk_mcs_muon_mom_v",
    ev.track_mcs_mom_mu_ );

  // BDT scores
  set_object_input_branch_address( etree, "protonBDTResponses",
    ev.proton_BDT_score_);
  set_object_input_branch_address( etree, "muonBDTResponses",
    ev.muon_BDT_score_);
  // Some ntuples exclude the old proton chi^2 PID score. Only include it
  // in the output if this branch is available.
  bool has_chipr = ( etree.GetBranch("trk_pid_chipr_v") != nullptr );
  if ( has_chipr ) {
    set_object_input_branch_address( etree, "trk_pid_chipr_v",
      ev.track_chi2_proton_ );
  }
  else {
    ev.track_chi2_proton_.reset( nullptr );
  }

  // Log-likelihood-based particle ID information
  set_object_input_branch_address( etree, "trk_llr_pid_v", ev.track_llr_pid_ );

  set_object_input_branch_address( etree, "trk_llr_pid_u_v",
    ev.track_llr_pid_U_ );

  set_object_input_branch_address( etree, "trk_llr_pid_v_v",
    ev.track_llr_pid_V_ );

  set_object_input_branch_address( etree, "trk_llr_pid_y_v",
    ev.track_llr_pid_Y_ );

  set_object_input_branch_address( etree, "trk_llr_pid_score_v",
    ev.track_llr_pid_score_ );

  // MC truth information for the neutrino
  etree.SetBranchAddress( "nu_pdg", &ev.mc_nu_pdg_ );
  etree.SetBranchAddress( "true_nu_vtx_x", &ev.mc_nu_vx_ );
  etree.SetBranchAddress( "true_nu_vtx_y", &ev.mc_nu_vy_ );
  etree.SetBranchAddress( "true_nu_vtx_z", &ev.mc_nu_vz_ );
  etree.SetBranchAddress( "nu_e", &ev.mc_nu_energy_ );
  etree.SetBranchAddress( "ccnc", &ev.mc_nu_ccnc_ );
  etree.SetBranchAddress( "interaction", &ev.mc_nu_interaction_type_ );

  // MC truth information for the final-state primary particles
  set_object_input_branch_address( etree, "mc_pdg", ev.mc_nu_daughter_pdg_ );
  set_object_input_branch_address( etree, "mc_E", ev.mc_nu_daughter_energy_ );
  set_object_input_branch_address( etree, "mc_px", ev.mc_nu_daughter_px_ );
  set_object_input_branch_address( etree, "mc_py", ev.mc_nu_daughter_py_ );
  set_object_input_branch_address( etree, "mc_pz", ev.mc_nu_daughter_pz_ );

  // GENIE and other systematic variation weights
  bool has_genie_mc_weights = ( etree.GetBranch("weightSpline") != nullptr );
  if ( has_genie_mc_weights ) {
    etree.SetBranchAddress( "weightSpline", &ev.spline_weight_ );
    etree.SetBranchAddress( "weightTune", &ev.tuned_cv_weight_ );
  }

  bool has_weight_map = ( etree.GetBranch("weights") != nullptr );
  if ( has_weight_map ) {
    set_object_input_branch_address( etree, "weights", ev.mc_weights_map_ );
  }
  else {
    ev.mc_weights_map_.reset( nullptr );
  }

  // Purity and completeness of the backtracked hits in the neutrino slice
  bool has_pfp_backtracked_purity = ( etree.GetBranch("nu_purity_from_pfp")
    != nullptr );
  if ( has_pfp_backtracked_purity ) {

    etree.SetBranchAddress( "nu_completeness_from_pfp",
      &ev.nu_completeness_from_pfp_ );

    etree.SetBranchAddress( "nu_purity_from_pfp", &ev.nu_purity_from_pfp_ );

  }

}

// Helper function to set branch addresses for the output TTree
void set_event_output_branch_addresses(TTree& out_tree, AnalysisEvent& ev,
  bool create = false)
{
  set_output_branch_address( out_tree, "is_mc", &ev.mc_pi_stopping_, create, "mc_pi_stopping/O" );
  set_output_branch_address( out_tree, "is_mc", &ev.mc_golden_, create, "mc_golden/O" );
  // Signal definition flags
  set_output_branch_address( out_tree, "is_mc", &ev.is_mc_, create, "is_mc/O" );

  set_output_branch_address( out_tree, "mc_neutrino_is_numu",
    &ev.mc_neutrino_is_numu_, create, "mc_neutrino_is_numu/O" );

  set_output_branch_address( out_tree, "mc_vertex_in_FV",
    &ev.mc_vertex_in_FV_, create, "mc_vertex_in_FV/O" );

  set_output_branch_address( out_tree, "mc_muon_in_mom_range",
    &ev.mc_muon_in_mom_range_, create, "mc_muon_in_mom_range/O" );

  set_output_branch_address( out_tree, "mc_lead_p_in_mom_range",
    &ev.mc_lead_p_in_mom_range_, create, "mc_lead_p_in_mom_range/O" );

  set_output_branch_address( out_tree, "mc_no_fs_pi0",
    &ev.mc_no_fs_pi0_, create, "mc_no_fs_pi0/O" );

  set_output_branch_address( out_tree, "mc_no_charged_pi_above_threshold",
    &ev.mc_no_charged_pi_above_threshold_, create,
    "mc_no_charged_pi_above_threshold/O" );

  set_output_branch_address( out_tree, "mc_no_fs_mesons",
    &ev.mc_no_fs_mesons_, create, "mc_no_fs_mesons/O" );

  set_output_branch_address( out_tree, "mc_is_signal",
    &ev.mc_is_signal_, create, "mc_is_signal/O" );

  set_output_branch_address( out_tree, "mc_1cpi",
    &ev.mc_1cpi_, create, "mc_1cpi/O" );

  set_output_branch_address( out_tree, "mc_1cpi",
    &ev.mc_no_neutrons_, create, "mc_no_neutrons/O" );
  // MC event category
  set_output_branch_address( out_tree, "category",
    &ev.category_, create, "category/I" );

  // Event weights
  set_output_branch_address( out_tree, "spline_weight",
    &ev.spline_weight_, create, "spline_weight/F" );

  set_output_branch_address( out_tree, "tuned_cv_weight",
    &ev.tuned_cv_weight_, create, "tuned_cv_weight/F" );

  // If MC weights are available, prepare to store them in the output TTree
  if ( ev.mc_weights_map_ ) {

    // Make separate branches for the various sets of systematic variation
    // weights in the map
    for ( auto& pair : *ev.mc_weights_map_ ) {

      // Prepend "weight_" to the name of the vector of weights in the map
      std::string weight_branch_name = "weight_" + pair.first;

      // Store a pointer to the vector of weights (needed to set the branch
      // address properly) in the temporary map of pointers
      ev.mc_weights_ptr_map_[ weight_branch_name ] = &pair.second;

      // Set the branch address for this vector of weights
      set_object_output_branch_address< std::vector<double> >( out_tree,
        weight_branch_name, ev.mc_weights_ptr_map_.at(weight_branch_name),
        create );
    }
  }

  
  // Backtracked neutrino purity and completeness
  set_output_branch_address( out_tree, "nu_completeness_from_pfp",
    &ev.nu_completeness_from_pfp_, create, "nu_completeness_from_pfp/F" );

  set_output_branch_address( out_tree, "nu_purity_from_pfp",
    &ev.nu_purity_from_pfp_, create, "nu_purity_from_pfp/F" );

  // Number of neutrino slices identified by the SliceID
  set_output_branch_address( out_tree, "nslice", &ev.nslice_, create,
    "nslice/I" );

  // CCNp0pi selection criteria
  set_output_branch_address( out_tree, "sel_nu_mu_cc", &ev.sel_nu_mu_cc_,
    create, "sel_nu_mu_cc/O" );

  set_output_branch_address( out_tree, "sel_reco_vertex_in_FV",
    &ev.sel_reco_vertex_in_FV_, create, "sel_reco_vertex_in_FV/O" );

  set_output_branch_address( out_tree, "sel_topo_cut_passed",
    &ev.sel_topo_cut_passed_, create, "sel_topo_cut_passed/O" );

  set_output_branch_address( out_tree, "sel_cosmic_ip_cut_passed",
    &ev.sel_cosmic_ip_cut_passed_, create, "sel_cosmic_ip_cut_passed/O" );

  set_output_branch_address( out_tree, "sel_pfp_starts_in_PCV",
    &ev.sel_pfp_starts_in_PCV_, create, "sel_pfp_starts_in_PCV/O" );

  set_output_branch_address( out_tree, "sel_no_reco_showers",
    &ev.sel_no_reco_showers_, create, "sel_no_reco_showers/O" );

  set_output_branch_address( out_tree, "sel_min_3_tracks",
    &ev.sel_min_3_tracks_, create, "sel_min_3_tracks/O" );

  set_output_branch_address( out_tree, "sel_2_non_proton",
    &ev.sel_2_non_proton_, create, "sel_2_non_proton/O" );

  set_output_branch_address( out_tree, "sel_has_pion_candidate",
    &ev.sel_has_pion_candidate_, create, "sel_has_pion_candidate/O");
 
  set_output_branch_address( out_tree, "sel_all_pfp_contained",
    &ev.sel_all_pfp_contained_, create, "sel_all_pfp_contained/O" );
 
  set_output_branch_address( out_tree, "sel_all_pfp_in_vtx_proximity",
    &ev.sel_all_pfp_in_vtx_proximity_, create, "sel_all_pfp_in_vtx_proximity/O" );

  set_output_branch_address( out_tree, "sel_has_muon_candidate",
    &ev.sel_has_muon_candidate_, create,
    "sel_has_muon_candidate/O" );

  set_output_branch_address( out_tree, "sel_muon_contained",
    &ev.sel_muon_contained_, create, "sel_muon_contained/O" );

  set_output_branch_address( out_tree, "sel_muon_passed_mom_cuts",
    &ev.sel_muon_passed_mom_cuts_, create, "sel_muon_passed_mom_cuts/O" );

  set_output_branch_address( out_tree, "sel_muon_quality_ok",
    &ev.sel_muon_quality_ok_, create, "sel_muon_quality_ok/O" );

  set_output_branch_address( out_tree, "sel_has_p_candidate",
    &ev.sel_has_p_candidate_, create, "sel_has_p_candidate/O" );

  set_output_branch_address( out_tree, "sel_passed_proton_pid_cut",
    &ev.sel_passed_proton_pid_cut_, create, "sel_passed_proton_pid_cut/O" );

  set_output_branch_address( out_tree, "sel_protons_contained",
    &ev.sel_protons_contained_, create, "sel_protons_contained/O" );

  set_output_branch_address( out_tree, "sel_lead_p_passed_mom_cuts",
    &ev.sel_lead_p_passed_mom_cuts_, create, "sel_lead_p_passed_mom_cuts/O" );

  set_output_branch_address( out_tree, "sel_CCNp1pi",
    &ev.sel_CCNp1pi_, create, "sel_CCNp1pi/O" );

  // Index for the muon candidate in the vectors of PFParticles
  set_output_branch_address( out_tree, "muon_candidate_pid_idx",
    &ev.muon_candidate_pid_idx_, create, "muon_candidate_pid_idx/I" );

  set_output_branch_address( out_tree, "muon_candidate_length_idx",
    &ev.muon_candidate_length_idx_, create, "muon_candidate_length_idx/I" );

  set_output_branch_address( out_tree, "muon_candidate_bdt_idx",
    &ev.muon_candidate_bdt_idx_, create, "muon_candidate_bdt_idx/I" );
  // pion candidate index
  set_output_branch_address( out_tree, "pion_candidate_idx",
    &ev.pion_candidate_idx_, create, "pion_candidate_idx/I" );

  // Index for the leading proton candidate in the vectors of PFParticles
  set_output_branch_address( out_tree, "lead_p_candidate_idx",
    &ev.lead_p_candidate_idx_, create, "lead_p_candidate_idx/I" );

  // Reco n_non_proton_like
  set_output_branch_address( out_tree,
    "n_non_proton_like", &ev.n_non_proton_like_, create, "n_non_proton_like/I" );

  // N reco track-like
  set_output_branch_address( out_tree,
    "n_reco_tracks", &ev.n_reco_tracks_, create, "n_reco_tracks/I" );

  //Golden pion handles

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "mc_end_p", ev.mc_end_p_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "mc_n_elastic", ev.mc_n_elastic_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "mc_n_inelastic", ev.mc_n_inelastic_, create );

  set_object_output_branch_address< std::vector<int> > ( out_tree,
    "pfp_n_descendents", ev.pfp_n_descendents_, create );

  set_object_output_branch_address< std::vector<int> > ( out_tree,
    "pfp_n_spacepoints", ev.trk_end_spacepoints_, create );
   
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_avg_deflection_stdev", ev.trk_avg_deflection_stdev_, create );
 
  // mutiplane bragg likelihoods
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_mip_uvw", ev.bragg_mip_uvw_, create);

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_p_fwd_uvw", ev.bragg_p_fwd_uvw_, create);
  
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_p_bwd_uvw", ev.bragg_p_bwd_uvw_, create);

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_p_to_MIP", ev.bragg_p_to_MIP_, create);
 
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_p_fwd_2_bwd", ev.bragg_p_fwd_2_bwd_, create); 

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_pion_fwd_uvw", ev.bragg_pion_fwd_uvw_, create);

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_pion_bwd_uvw", ev.bragg_pion_bwd_uvw_, create);

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_pion_to_MIP", ev.bragg_pion_to_MIP_, create);

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_pion_fwd_2_bwd", ev.bragg_pion_fwd_2_bwd_, create);

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_mu_fwd_uvw", ev.bragg_mu_fwd_uvw_, create);

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_mu_bwd_uvw", ev.bragg_mu_bwd_uvw_, create);

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "truncated_mean_dEdx", ev.truncated_mean_dEdx_, create);

  // single plane bragg likelihoods

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_trunk_dEdx_w", ev.trk_trunk_dEdx_w_, create);

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_trunk_dEdx_u", ev.trk_trunk_dEdx_u_, create);

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_trunk_dEdx_v", ev.trk_trunk_dEdx_v_, create);

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_pion_fwd_w", ev.bragg_pion_fwd_w_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_pion_fwd_u", ev.bragg_pion_fwd_u_, create );
  
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_pion_fwd_v", ev.bragg_pion_fwd_v_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_mu_fwd_w", ev.bragg_mu_fwd_w_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_mu_fwd_u", ev.bragg_mu_fwd_u_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_mu_fwd_v", ev.bragg_mu_fwd_v_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_p_fwd_w", ev.bragg_p_fwd_w_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_p_fwd_u", ev.bragg_p_fwd_u_, create );
 
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_p_fwd_v", ev.bragg_p_fwd_v_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_pion_bwd_w", ev.bragg_pion_bwd_w_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_pion_bwd_u", ev.bragg_pion_bwd_u_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_pion_bwd_v", ev.bragg_pion_bwd_v_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_mu_bwd_w", ev.bragg_mu_bwd_w_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_mu_bwd_u", ev.bragg_mu_bwd_u_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_mu_bwd_v", ev.bragg_mu_bwd_v_, create );
 
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_p_bwd_w", ev.bragg_p_bwd_w_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_p_bwd_u", ev.bragg_p_bwd_u_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_p_bwd_v", ev.bragg_p_bwd_v_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_mip_w", ev.bragg_mip_w_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_mip_v", ev.bragg_mip_v_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "bragg_mip_u", ev.bragg_mip_u_, create );
  
  // BDT scores
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "muon_BDT_score", ev.muon_BDT_score_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "proton_BDT_score", ev.proton_BDT_score_, create );
  // Reco 3-momenta (muon, leading proton, pion)
  set_object_output_branch_address< TVector3 >( out_tree,
    "p3_mu", ev.p3_mu_, create );

  set_object_output_branch_address< TVector3 >( out_tree,
    "p3_lead_p", ev.p3_lead_p_, create );

  set_object_output_branch_address< TVector3 >( out_tree,
    "p3_cpi", ev.p3_cpi_, create );
  // Reco 3-momenta (all proton candidates, ordered from highest to lowest
  // magnitude)
  set_object_output_branch_address< std::vector<TVector3> >( out_tree,
    "p3_p_vec", ev.p3_p_vec_, create );

  // True 3-momenta (muon, leading proton, leading charged pion)
  set_object_output_branch_address< TVector3 >( out_tree,
    "mc_p3_mu", ev.mc_p3_mu_, create );

  set_object_output_branch_address< TVector3 >( out_tree,
    "mc_p3_lead_p", ev.mc_p3_lead_p_, create );

  set_object_output_branch_address< TVector3 >( out_tree,
    "mc_p3_cpi", ev.mc_p3_cpi_, create );

  // True 3-momenta (all protons, ordered from highest to lowest magnitude)
  set_object_output_branch_address< std::vector<TVector3> >( out_tree,
    "mc_p3_p_vec", ev.mc_p3_p_vec_, create );

  set_output_branch_address( out_tree, "mc_n_protons",
    &ev.mc_n_protons_, create, "mc_n_protons/I" );
  
  set_object_output_branch_address< std::vector<TVector3> >( out_tree,
    "mc_p3_cpi_vec", ev.mc_p3_cpi_vec_, create );
 
  set_output_branch_address( out_tree, "mc_n_cpi",
      &ev.mc_n_cpi_, create, "n_cpi/I" );
  // Reco STVs
  set_output_branch_address( out_tree, "muE",
    &ev.muE_, create, "muE/F" );

  set_output_branch_address( out_tree, "pE",
    &ev.pE_, create, "pE/F" );

  set_output_branch_address( out_tree, "piE",
    &ev.piE_, create, "piE/F" );

  set_output_branch_address( out_tree, "reco_Ecal",
    &ev.reco_Ecal_, create, "reco_Ecal/F" );

  set_output_branch_address( out_tree, "pKE",
    &ev.pKE_, create, "pKE/F" );
  
  set_object_output_branch_address< TVector3 >( out_tree, "q", ev.q_, create ); 

  set_output_branch_address( out_tree, "delta_pT",
    &ev.delta_pT_, create, "delta_pT/F" );
  
  set_output_branch_address( out_tree, "delta_pL",
    &ev.delta_pL_, create, "delta_pL/F" );

   set_output_branch_address( out_tree, "delta_pL2",
    &ev.delta_pL2_, create, "delta_pL2/F" );

  set_output_branch_address( out_tree, "pn",
    &ev.pn_, create, "pn/F" );

  set_output_branch_address( out_tree, "pn2",
    &ev.pn2_, create, "pn2/F" );

  set_output_branch_address( out_tree, "theta_mu_p",
    &ev.theta_mu_p_, create, "theta_mu_p/F" );

  set_output_branch_address( out_tree, "theta_mu_cpi",
    &ev.theta_mu_cpi_, create, "theta_mu_cpi/F" );

  set_output_branch_address( out_tree, "delta_alpha3D",
    &ev.delta_alpha3D_, create, "delta_alpha3D/F" );

  set_output_branch_address( out_tree, "delta_phi3d_had",
    &ev.delta_phi3d_had_, create, "delta_phi3d_had/F" );


  set_output_branch_address( out_tree, "delta_phi3d_mu",
    &ev.delta_phi3d_mu_, create, "delta_phi3d_mu/F" );
 
  set_output_branch_address( out_tree, "delta_phi3d_p",
    &ev.delta_phi3d_p_, create, "delta_phi3d_p/F" );

  set_output_branch_address( out_tree, "delta_phi3d_pi",
    &ev.delta_phi3d_pi_, create, "delta_phi3d_pi/F" ); 
  // MC STVs (only filled for signal events)
  set_output_branch_address( out_tree, "mc_muE",
    &ev.mc_muE_, create, "mc_muE/F" );

  set_output_branch_address( out_tree, "mc_pE",
    &ev.mc_pE_, create, "mc_pE/F" );

  set_output_branch_address( out_tree, "mc_piE",
    &ev.mc_piE_, create, "mc_piE/F" );

  set_output_branch_address( out_tree, "mc_Ecal",
    &ev.mc_Ecal_, create, "mc_Ecal/F" );

  set_output_branch_address( out_tree, "mc_pKE",
    &ev.mc_pKE_, create, "mc_pKE/F" );

  set_object_output_branch_address< TVector3 >( out_tree,
    "mc_q", ev.mc_q_, create );

  set_output_branch_address( out_tree, "mc_delta_pT",
    &ev.mc_delta_pT_, create, "mc_delta_pT/F" );

  //set_output_branch_address( out_tree, "mc_delta_phiT",
    //&ev.mc_delta_phiT_, create, "mc_delta_phiT/F" );

  //set_output_branch_address( out_tree, "mc_delta_alphaT",
    //&ev.mc_delta_alphaT_, create, "mc_delta_alphaT/F" );

  set_output_branch_address( out_tree, "mc_delta_pL",
    &ev.mc_delta_pL2_, create, "mc_delta_pL/F" );

  set_output_branch_address( out_tree, "mc_delta_pL2",
    &ev.mc_delta_pL_, create, "mc_delta_pL2/F" );

  set_output_branch_address( out_tree, "mc_pn",
    &ev.mc_pn_, create, "mc_pn/F" );

  set_output_branch_address( out_tree, "mc_pn2",
    &ev.mc_pn2_, create, "mc_pn2/F" );

  //set_output_branch_address( out_tree, "mc_delta_pTx",
    //&ev.mc_delta_pTx_, create, "mc_delta_pTx/F" );

  //set_output_branch_address( out_tree, "mc_delta_pTy",
    //&ev.mc_delta_pTy_, create, "mc_delta_pTy/F" );

  set_output_branch_address( out_tree, "mc_theta_mu_p",
    &ev.mc_theta_mu_p_, create, "mc_theta_mu_p/F" );

  set_output_branch_address( out_tree, "mc_theta_mu_cpi",
    &ev.mc_theta_mu_cpi_, create, "mc_theta_mu_cpi/F" );
  
  set_output_branch_address( out_tree, "mc_delta_alpha3D",
    &ev.mc_delta_alpha3D_, create, "mc_delta_alpha3D/F" );

  set_output_branch_address( out_tree, "mc_delta_phi3d_had",
    &ev.mc_delta_phi3d_had_, create, "mc_delta_phi3d_had/F" );


  set_output_branch_address( out_tree, "mc_delta_phi3d_mu",
    &ev.mc_delta_phi3d_mu_, create, "mc_delta_phi3d_mu/F" );

  set_output_branch_address( out_tree, "mc_delta_phi3d_p",
    &ev.mc_delta_phi3d_p_, create, "mc_delta_phi3d_p/F" );
  set_output_branch_address( out_tree, "mc_delta_phi3d_pi",
    &ev.mc_delta_phi3d_pi_, create, "mc_delta_phi3d_pi/F" );
  // *** Branches copied directly from the input ***

  // Cosmic rejection parameters for numu CC inclusive selection
  set_output_branch_address( out_tree, "topological_score",
    &ev.topological_score_, create, "topological_score/F" );

  set_output_branch_address( out_tree, "CosmicIP",
    &ev.cosmic_impact_parameter_, create, "CosmicIP/F" );

  // Reconstructed neutrino vertex position
  set_output_branch_address( out_tree, "reco_nu_vtx_sce_x",
    &ev.nu_vx_, create, "reco_nu_vtx_sce_x/F" );

  set_output_branch_address( out_tree, "reco_nu_vtx_sce_y",
    &ev.nu_vy_, create, "reco_nu_vtx_sce_y/F" );

  set_output_branch_address( out_tree, "reco_nu_vtx_sce_z",
    &ev.nu_vz_, create, "reco_nu_vtx_sce_z/F" );

  // MC truth information for the neutrino
  set_output_branch_address( out_tree, "mc_nu_pdg", &ev.mc_nu_pdg_,
    create, "mc_nu_pdg/I" );

  set_output_branch_address( out_tree, "mc_nu_vtx_x", &ev.mc_nu_vx_,
    create, "mc_nu_vtx_x/F" );

  set_output_branch_address( out_tree, "mc_nu_vtx_y", &ev.mc_nu_vy_,
    create, "mc_nu_vtx_y/F" );

  set_output_branch_address( out_tree, "mc_nu_vtx_z", &ev.mc_nu_vz_,
    create, "mc_nu_vtx_z/F" );

  set_output_branch_address( out_tree, "mc_nu_energy", &ev.mc_nu_energy_,
    create, "mc_nu_energy/F" );

  set_output_branch_address( out_tree, "mc_ccnc", &ev.mc_nu_ccnc_,
    create, "mc_ccnc/I" );

  set_output_branch_address( out_tree, "mc_interaction",
    &ev.mc_nu_interaction_type_, create, "mc_interaction/I" );

  // PFParticle properties
  set_object_output_branch_address< std::vector<unsigned int> >( out_tree,
    "pfp_generation_v", ev.pfp_generation_, create );

  set_object_output_branch_address< std::vector<unsigned int> >( out_tree,
    "pfp_trk_daughters_v", ev.pfp_trk_daughters_count_, create );

  set_object_output_branch_address< std::vector<unsigned int> >( out_tree,
    "pfp_shr_daughters_v", ev.pfp_shr_daughters_count_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_score_v", ev.pfp_track_score_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "pfpdg", ev.pfp_reco_pdg_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "pfnhits", ev.pfp_hits_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "pfnplanehits_U", ev.pfp_hitsU_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "pfnplanehits_V", ev.pfp_hitsV_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "pfnplanehits_Y", ev.pfp_hitsY_, create );

  // Backtracked PFParticle properties
  set_object_output_branch_address< std::vector<int> >( out_tree,
    "backtracked_pdg", ev.pfp_true_pdg_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "backtracked_e", ev.pfp_true_E_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "backtracked_px", ev.pfp_true_px_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "backtracked_py", ev.pfp_true_py_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "backtracked_pz", ev.pfp_true_pz_, create );

  // Shower properties
  // For some ntuples, reconstructed shower information is excluded.
  // In such cases, skip writing these branches to the output TTree.
  if ( ev.shower_startx_ ) {
    set_object_output_branch_address< std::vector<float> >( out_tree,
      "shr_start_x_v", ev.shower_startx_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "shr_start_y_v", ev.shower_starty_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "shr_start_z_v", ev.shower_startz_, create );

    // Shower start distance from reco neutrino vertex (pre-calculated for
    // convenience)
    set_object_output_branch_address< std::vector<float> >( out_tree,
      "shr_dist_v", ev.shower_start_distance_, create );
  }

  // Track properties
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_len_v", ev.track_length_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_sce_start_x_v", ev.track_startx_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_sce_start_y_v", ev.track_starty_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_sce_start_z_v", ev.track_startz_, create );

  // Track start distance from reco neutrino vertex (pre-calculated for
  // convenience)
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_distance_v", ev.track_start_distance_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_sce_end_x_v", ev.track_endx_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_sce_end_y_v", ev.track_endy_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_sce_end_z_v", ev.track_endz_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_dir_x_v", ev.track_dirx_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_dir_y_v", ev.track_diry_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_dir_z_v", ev.track_dirz_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_energy_proton_v", ev.track_kinetic_energy_p_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_range_muon_mom_v", ev.track_range_mom_mu_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_mcs_muon_mom_v", ev.track_mcs_mom_mu_, create );

  // Some ntuples exclude the old chi^2 proton PID score. Only include it in
  // the output if it is available.
  if ( ev.track_chi2_proton_ ) {
    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_pid_chipr_v", ev.track_chi2_proton_, create );
  }

  // Log-likelihood-based particle ID information
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_llr_pid_v", ev.track_llr_pid_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_llr_pid_u_v", ev.track_llr_pid_U_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_llr_pid_v_v", ev.track_llr_pid_V_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_llr_pid_y_v", ev.track_llr_pid_Y_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_llr_pid_score_v", ev.track_llr_pid_score_, create );

  // MC truth information for the final-state primary particles
  set_object_output_branch_address< std::vector<int> >( out_tree, "mc_pdg",
    ev.mc_nu_daughter_pdg_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree, "mc_E",
    ev.mc_nu_daughter_energy_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree, "mc_px",
    ev.mc_nu_daughter_px_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree, "mc_py",
    ev.mc_nu_daughter_py_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree, "mc_pz",
    ev.mc_nu_daughter_pz_, create );

}

void analyze(const std::vector<std::string>& in_file_names,
  const std::string& output_filename)
{
  // Get the TTrees containing the event ntuples and subrun POT information
  // Use TChain objects for simplicity in manipulating multiple files
  TChain events_ch( "NeutrinoSelectionFilter" );
  TChain subruns_ch( "SubRun" );

  for ( const auto& f_name : in_file_names ) {
    events_ch.Add( f_name.c_str() );
    subruns_ch.Add( f_name.c_str() );
  }

  // OUTPUT TTREE
  // Make an output TTree for plotting (one entry per event)
  TFile* out_file = new TFile( output_filename.c_str(), "recreate" );
  out_file->cd();
  TTree* out_tree = new TTree( "stv_tree", "STV analysis tree" );

  // Get the total POT from the subruns TTree. Save it in the output
  // TFile as a TParameter<float>. Real data doesn't have this TTree,
  // so check that it exists first.
  float pot;
  float summed_pot = 0.;
  bool has_pot_branch = ( subruns_ch.GetBranch("pot") != nullptr );
  if ( has_pot_branch ) {
    subruns_ch.SetBranchAddress( "pot", &pot );
    for ( int se = 0; se < subruns_ch.GetEntries(); ++se ) {
      subruns_ch.GetEntry( se );
      summed_pot += pot;
    }
  }

  TParameter<float>* summed_pot_param = new TParameter<float>( "summed_pot",
    summed_pot );

  summed_pot_param->Write();

  // EVENT LOOP
  // TChains can potentially be really big (and spread out over multiple
  // files). When that's the case, calling TChain::GetEntries() can be very
  // slow. I get around this by using a while loop instead of a for loop.
  bool created_output_branches = false;
  long events_entry = 0;
  while ( true ) {
    //if (events_entry == 100) break;

    if ( events_entry % 1000 == 0 ) {
      std::cout << "Processing event #" << events_entry << '\n';
    }

    // Create a new AnalysisEvent object. This will reset all analysis
    // variables for the current event.
    AnalysisEvent cur_event;

    // Set branch addresses for the member variables that will be read
    // directly from the Event TTree.
    set_event_branch_addresses( events_ch, cur_event );

    // TChain::LoadTree() returns the entry number that should be used with
    // the current TTree object, which (together with the TBranch objects
    // that it owns) doesn't know about the other TTrees in the TChain.
    // If the return value is negative, there was an I/O error, or we've
    // attempted to read past the end of the TChain.
    int local_entry = events_ch.LoadTree( events_entry );

    // If we've reached the end of the TChain (or encountered an I/O error),
    // then terminate the event loop
    if ( local_entry < 0 ) break;

    // Load all of the branches for which we've called
    // TChain::SetBranchAddress() above
    events_ch.GetEntry( events_entry );

    // Set the output TTree branch addresses, creating the branches if needed
    // (during the first event loop iteration)
    bool create_them = false;
    if ( !created_output_branches ) {
      create_them = true;
      created_output_branches = true;
    }

    set_event_output_branch_addresses( *out_tree, cur_event, create_them );
    
    //if (num_pf_particles_ == BOGUS_INT || num_pf_particles_ > 10000000)  
    // Apply the CCNp0pi selection criteria and categorize the event.
    cur_event.apply_selection();

    // Compute observables to save to the output TTree
    cur_event.compute_observables();

    // We're done. Save the results and move on to the next event.
    out_tree->Fill();
    ++events_entry;
  }

  out_tree->Write();
  out_file->Close();
  delete out_file;
}

// Sets the signal definition flags and returns an event category based on MC
// truth information
EventCategory AnalysisEvent::categorize_event() {

  // Real data has a bogus true neutrino PDG code that is not one of the
  // allowed values (±12, ±14, ±16)
  int abs_mc_nu_pdg = std::abs( mc_nu_pdg_ );
  is_mc_ = ( abs_mc_nu_pdg == ELECTRON_NEUTRINO
    || abs_mc_nu_pdg == MUON_NEUTRINO || abs_mc_nu_pdg == TAU_NEUTRINO );
  if ( !is_mc_ ) return kUnknown;

  mc_vertex_in_FV_ = mc_vertex_inside_FV();
  mc_neutrino_is_numu_ = ( mc_nu_pdg_ == MUON_NEUTRINO );

  if ( !mc_vertex_in_FV_ ) {
    mc_is_signal_ = false;
    return kOOFV;
  }
  else if ( mc_nu_ccnc_ == NEUTRAL_CURRENT ) {
    mc_is_signal_ = false;
    return kNC;
  }
  else if ( !mc_neutrino_is_numu_ ) {
    mc_is_signal_ = false;
    if ( mc_nu_pdg_ == ELECTRON_NEUTRINO
      && mc_nu_ccnc_ == CHARGED_CURRENT ) return kNuECC;
    else return kOther;
  }

  // Set flags to their default values here
  mc_muon_in_mom_range_ = false;
  mc_lead_p_in_mom_range_ = false;
  mc_no_fs_pi0_ = true;
  mc_no_charged_pi_above_threshold_ = true;
  mc_no_fs_mesons_ = true;
  mc_1cpi_ = false;
  mc_no_neutrons_ = true;
  mc_pi_stopping_ = false;
  mc_golden_ = false;

  double lead_p_mom = LOW_FLOAT;
  double lead_cpi_mom = LOW_FLOAT; //K

  for ( size_t p = 0u; p < mc_nu_daughter_pdg_->size(); ++p ) {
    int pdg = mc_nu_daughter_pdg_->at( p );
    float energy = mc_nu_daughter_energy_->at( p );

    // Do the general check for (anti)mesons first before considering
    // any individual PDG codes
    // make exception for charged pion
    if ( is_meson_or_antimeson(pdg) ) {
      mc_no_fs_mesons_ = false;
    }

    if (pdg == NEUTRON) {
       mc_no_neutrons_ = false;
    }

    // Check that the muon has a momentum within the allowed range
    if ( pdg == MUON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(MUON_MASS, 2) );
      if ( mom >= MUON_P_MIN_MOM_CUT && mom <= MUON_P_MAX_MOM_CUT ) {
        mc_muon_in_mom_range_ = true;
      }
    }
    else if ( pdg == PROTON ) {
      mc_n_protons_ += 1;
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PROTON_MASS, 2) );
      if ( mom > lead_p_mom ) lead_p_mom = mom;
    }
    else if ( pdg == PI_ZERO ) {
      mc_no_fs_pi0_ = false;
    }
    // check charged pion above threshold present
    else if ( std::abs(pdg) == PI_PLUS ) {
      // check if golden
      mc_pi_stopping_ = ( mc_end_p_->at(p) <= std::numeric_limits<float>::epsilon() );
      int n_scatters = mc_n_elastic_->at(p) + mc_n_inelastic_->at(p);

      mc_golden_ = (mc_pi_stopping_ && n_scatters == 0);

      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PI_PLUS_MASS, 2) );     
      if ( mom > CHARGED_PI_MOM_CUT ) {
        mc_no_charged_pi_above_threshold_ = false;
	mc_n_cpi_ += 1;
	lead_cpi_mom = mom; //K	
      }
    }
  }
  
  // Check number of charged pions
  if (mc_n_cpi_ == 1) {
    mc_1cpi_ = true;
  }

  // Check that the leading proton has a momentum within the allowed range
  if ( lead_p_mom >= LEAD_P_MIN_MOM_CUT && lead_p_mom <= LEAD_P_MAX_MOM_CUT  ) {
    mc_lead_p_in_mom_range_ = true;
  }

  mc_is_signal_ = mc_vertex_in_FV_ && mc_neutrino_is_numu_
    && mc_muon_in_mom_range_ && mc_lead_p_in_mom_range_
    && mc_no_fs_mesons_ && mc_no_fs_pi0_ && mc_1cpi_;
  


  // Sort signal by interaction mode
  if ( mc_is_signal_ ) {
    if ( mc_nu_interaction_type_ == 0 ) return kSignalCCQE; // QE
    else if ( mc_nu_interaction_type_ == 10 ) return kSignalCCMEC; // MEC
    else if ( mc_nu_interaction_type_ == 1 ) return kSignalCCRES; // RES
    else if ( mc_nu_interaction_type_ == 2 )  return kSignalCCDIS;// DIS
    else if ( mc_nu_interaction_type_ == 3 ) return kSignalCCCOH;// COH
    else return kSignalOther;
  }

  // any pion (charged or neutral) and proton inclusice (Xp, X=0,1,2...)
  /*else if ( !mc_no_fs_pi0_ &&  (mc_n_cpi_>=1) ) {
    C1pi0p
  }*/

  // pi0 proton inclusive, with pi0s any number of mesons allowed
  else if (!mc_no_fs_pi0_) {
    return kNuMuCCpi0;
  }

  // at least 1 proton above threshold, no charged pions 
  else if ( mc_muon_in_mom_range_ && mc_no_charged_pi_above_threshold_) {
    return kNuMuCC0piXp; 
  }
  
  // one charged pion above threhold, no protons, any number of mesons
  //else if (mc_1cpi_ && !mc_lead_p_in_mom_range_ ) {
  //  return kNuMuCC1pi0p;

  //}
  else return kNuMuCCOther;
}

float compute_multi_plane_Bragg_likelihood(float &likelihood_w, float &likelihood_u, float &likelihood_v, float yz_angle, int n_hits_u, int n_hits_v){

  bool trackAlongW = ( (std::pow(std::sin(yz_angle), 2) < sin2AngleThreshold ));
  
  bool hasW = (!(trackAlongW) && likelihood_w >= -1.f);

  if (hasW) return likelihood_w; 
  

  bool trackAlongU = ( (std::pow( std::sin(yz_angle - piBy3),2 ) < sin2AngleThreshold));
  bool trackAlongV = ( (std::pow( std::sin(yz_angle + piBy3),2 ) < sin2AngleThreshold));

  bool hasU = (!(trackAlongU) && likelihood_u >= -1.f);
  bool hasV = (!(trackAlongV) && likelihood_v >= -1.f);

  auto dofU = hasU ? static_cast<float>(n_hits_u) : 0.f;
  auto dofV = hasV ? static_cast<float>(n_hits_v) : 0.f;
  auto dofUV = dofU + dofV;

  auto wU = hasU ? (likelihood_u * dofU ) : 0.f;
  auto wV = hasV ? (likelihood_v *dofV ) :0.f;
  
  if ((!hasU && !hasV) || (n_hits_u == 0 && n_hits_v == 0) || dofUV <= std::numeric_limits<float>::epsilon()) return BOGUS;
  
  auto likelihood = (wU + wV) / dofUV;
  if (likelihood >= -1.f) return likelihood;
  else return BOGUS;
}

float compute_log_likelihood_ratio(float &likelihood1, float &likelihood2){
 
  if (likelihood1 == BOGUS || likelihood2 == BOGUS) return BOGUS;
  if (likelihood2 <= std::numeric_limits<float>::epsilon() ) return BOGUS;
  return std::log(likelihood1/likelihood2); 
  
}

float compute_softmax(float &fwd, float &bwd){

  if (fwd == BOGUS || bwd == BOGUS) return BOGUS;
  
  float exp_sum = std::exp(fwd) + std::exp(bwd);
  if (exp_sum <= std::numeric_limits<float>::epsilon()) return BOGUS;
  
  return std::exp(fwd) / exp_sum;
}

float compute_multi_plane_trunc_mean_dEdx(float yz_angle, float trunc_mean_w, float trunc_mean_u, float trunc_mean_v, int n_hits_u, int n_hits_v){

  bool trackAlongW = ( (std::pow(std::sin(yz_angle), 2) < sin2AngleThreshold ));
  if (trackAlongW && (trunc_mean_w != floatMinValue)) return trunc_mean_w;
  else{
    const auto uWeight = std::max(0.f, (n_hits_u * lengthFraction) - nHitsToSkip);
    const auto vWeight = std::max(0.f, (n_hits_v * lengthFraction) - nHitsToSkip);
    const auto hasU = (uWeight>0.f && trunc_mean_u != floatMinValue);
    const auto hasV = (vWeight>0.f && trunc_mean_v != floatMinValue);
    
    if (hasU || hasV){
      float truncatedMeandEdx = 0.f;
      if (hasU) truncatedMeandEdx += trunc_mean_u * uWeight;
      if (hasV) truncatedMeandEdx += trunc_mean_v * vWeight;
      truncatedMeandEdx /= (uWeight + vWeight);
      return truncatedMeandEdx;
    }
    else return BOGUS; 
  }
}

void AnalysisEvent::get_multi_plane_Bragg_likelihood(){

  for ( int p = 0; p < num_pf_particles_; ++p ) {
    unsigned int generation = pfp_generation_->at( p );
    if ( generation != 2u ) continue;
    //std::cout << trk_bragg_p_fwd_preferred_w_->at( p ) << std::endl;

    // particle YZ angle & hits in u,v planes
    float yz_angle = std::atan2( track_dirz_->at( p ), track_diry_->at( p ) );
    int n_hits_u = pfp_hitsU_->at( p );
    int n_hits_v = pfp_hitsV_->at( p );
  
    // mip likelikelihoods
    float mip_w = ( bragg_mip_w_->at( p ) >= -1.f) ? bragg_mip_w_->at( p ) : BOGUS;
    float mip_u = ( bragg_mip_u_->at( p ) >= -1.f) ? bragg_mip_u_->at( p ) : BOGUS;
    float mip_v = ( bragg_mip_v_->at( p ) >= -1.f) ? bragg_mip_v_->at( p ) : BOGUS;

    float mip_uvw = compute_multi_plane_Bragg_likelihood(mip_w, mip_u, mip_v, yz_angle, n_hits_u, n_hits_v);
   
    bragg_mip_uvw_->push_back(mip_uvw);
    // proton likelihoods
    float proton_fwd_w = trk_bragg_p_fwd_preferred_w_->at( p ) ? trk_bragg_p_w_->at( p ) : trk_bragg_p_alt_dir_w_->at( p );
    float proton_bwd_w = ( !(trk_bragg_p_fwd_preferred_w_->at( p )) ) ? trk_bragg_p_w_->at( p ) : trk_bragg_p_alt_dir_w_->at( p );
       
    float proton_fwd_u = trk_bragg_p_fwd_preferred_u_->at( p ) ? trk_bragg_p_u_->at( p ) : trk_bragg_p_alt_dir_u_->at( p );
    float proton_bwd_u = ( !(trk_bragg_p_fwd_preferred_u_->at( p )) ) ? trk_bragg_p_u_->at( p ) : trk_bragg_p_alt_dir_u_->at( p );
  
    float proton_fwd_v = trk_bragg_p_fwd_preferred_v_->at( p ) ? trk_bragg_p_v_->at( p ) : trk_bragg_p_alt_dir_v_->at( p );
    float proton_bwd_v = ( !(trk_bragg_p_fwd_preferred_v_->at( p )) ) ? trk_bragg_p_v_->at( p ) : trk_bragg_p_alt_dir_v_->at( p );

    bragg_p_fwd_w_->push_back(proton_fwd_w);
    bragg_p_bwd_w_->push_back(proton_bwd_w);
 
    bragg_p_fwd_u_->push_back(proton_fwd_u);
    bragg_p_bwd_v_->push_back(proton_bwd_u);

    bragg_p_fwd_w_->push_back(proton_fwd_v);
    bragg_p_bwd_w_->push_back(proton_bwd_v);


    float proton_fwd_uvw = compute_multi_plane_Bragg_likelihood(proton_fwd_w, proton_fwd_u, proton_fwd_v, yz_angle, n_hits_u, n_hits_v);
    float proton_bwd_uvw = compute_multi_plane_Bragg_likelihood(proton_bwd_w, proton_bwd_u, proton_bwd_v, yz_angle, n_hits_u, n_hits_v);

    bragg_p_fwd_uvw_->push_back( proton_fwd_uvw );
    bragg_p_bwd_uvw_->push_back( proton_bwd_uvw );

    bragg_p_to_MIP_->push_back( compute_log_likelihood_ratio(proton_fwd_uvw, mip_uvw) );

    bragg_p_fwd_2_bwd_->push_back( compute_softmax(proton_fwd_uvw, proton_bwd_uvw) );
       
    // muon likelihoods
    float mu_fwd_w = trk_bragg_mu_fwd_preferred_w_->at( p ) ? trk_bragg_mu_w_->at( p ) : trk_bragg_mu_alt_dir_w_->at( p );
    float mu_bwd_w = ( !(trk_bragg_mu_fwd_preferred_w_->at( p )) ) ? trk_bragg_mu_w_->at( p ) : trk_bragg_mu_alt_dir_w_->at( p );

    float mu_fwd_u = trk_bragg_mu_fwd_preferred_u_->at( p ) ? trk_bragg_mu_u_->at( p ) : trk_bragg_mu_alt_dir_u_->at( p );
    float mu_bwd_u = ( !(trk_bragg_mu_fwd_preferred_u_->at( p )) ) ? trk_bragg_mu_u_->at( p ) : trk_bragg_mu_alt_dir_u_->at( p );

    float mu_fwd_v = trk_bragg_mu_fwd_preferred_v_->at( p ) ? trk_bragg_mu_v_->at( p ) : trk_bragg_mu_alt_dir_v_->at( p );
    float mu_bwd_v = ( !(trk_bragg_mu_fwd_preferred_v_->at( p )) ) ? trk_bragg_mu_v_->at( p ) : trk_bragg_mu_alt_dir_v_->at( p );

    bragg_mu_fwd_w_->push_back(mu_fwd_w);
    bragg_mu_bwd_w_->push_back(mu_bwd_w);

    bragg_mu_fwd_u_->push_back(mu_fwd_u);
    bragg_mu_bwd_v_->push_back(mu_bwd_u);

    bragg_mu_fwd_w_->push_back(mu_fwd_v);
    bragg_mu_bwd_w_->push_back(mu_bwd_v);


    float mu_fwd_uvw = compute_multi_plane_Bragg_likelihood(mu_fwd_w, mu_fwd_u, mu_fwd_v, yz_angle, n_hits_u, n_hits_v);
    float mu_bwd_uvw = compute_multi_plane_Bragg_likelihood(mu_bwd_w, mu_bwd_u, mu_bwd_v, yz_angle, n_hits_u, n_hits_v);

    bragg_mu_fwd_uvw_->push_back( mu_fwd_uvw );
    bragg_mu_bwd_uvw_->push_back( mu_bwd_uvw );

    // pion likelihoods
    float pion_fwd_w = trk_bragg_pion_fwd_preferred_w_->at( p ) ? trk_bragg_pion_w_->at( p ) : trk_bragg_pion_alt_dir_w_->at( p );
    float pion_bwd_w = ( !(trk_bragg_pion_fwd_preferred_w_->at( p )) ) ? trk_bragg_pion_w_->at( p ) : trk_bragg_pion_alt_dir_w_->at( p );

    float pion_fwd_u = trk_bragg_pion_fwd_preferred_u_->at( p ) ? trk_bragg_pion_u_->at( p ) : trk_bragg_pion_alt_dir_u_->at( p );
    float pion_bwd_u = ( !(trk_bragg_pion_fwd_preferred_u_->at( p )) ) ? trk_bragg_pion_u_->at( p ) : trk_bragg_pion_alt_dir_u_->at( p );

    float pion_fwd_v = trk_bragg_pion_fwd_preferred_v_->at( p ) ? trk_bragg_pion_v_->at( p ) : trk_bragg_pion_alt_dir_v_->at( p );
    float pion_bwd_v = ( !(trk_bragg_pion_fwd_preferred_v_->at( p )) ) ? trk_bragg_pion_v_->at( p ) : trk_bragg_pion_alt_dir_v_->at( p );

    bragg_pion_fwd_w_->push_back(pion_fwd_w);
    bragg_pion_bwd_w_->push_back(pion_bwd_w);

    bragg_pion_fwd_u_->push_back(pion_fwd_u);
    bragg_pion_bwd_v_->push_back(pion_bwd_u);

    bragg_pion_fwd_w_->push_back(pion_fwd_v);
    bragg_pion_bwd_w_->push_back(pion_bwd_v);

    float pion_fwd_uvw = compute_multi_plane_Bragg_likelihood(pion_fwd_w, pion_fwd_u, pion_fwd_v, yz_angle, n_hits_u, n_hits_v);
    float pion_bwd_uvw = compute_multi_plane_Bragg_likelihood(pion_bwd_w, pion_bwd_u, pion_bwd_v, yz_angle, n_hits_u, n_hits_v);

    bragg_pion_fwd_uvw_->push_back( pion_fwd_uvw );
    bragg_pion_bwd_uvw_->push_back( pion_bwd_uvw );

    bragg_pion_to_MIP_->push_back( compute_log_likelihood_ratio(pion_fwd_uvw, mip_uvw) );

    bragg_pion_fwd_2_bwd_->push_back( compute_softmax(pion_fwd_uvw, pion_bwd_uvw) );

    truncated_mean_dEdx_->push_back( compute_multi_plane_trunc_mean_dEdx(yz_angle, trk_trunk_dEdx_w_->at(p), trk_trunk_dEdx_u_->at(p), trk_trunk_dEdx_v_->at(p), n_hits_u, n_hits_v ) );
  }
}



void AnalysisEvent::apply_numu_CC_selection() {
  //use nu_pdg_ == slpdg (from PeLEE ntuple)
  sel_pandoraPDG_isnumu_ = (nu_pdg_ == MUON_NEUTRINO);
  sel_reco_vertex_in_FV_ = this->reco_vertex_inside_FV();
  sel_topo_cut_passed_ = topological_score_ > TOPO_SCORE_CUT;
  sel_cosmic_ip_cut_passed_ = cosmic_impact_parameter_ > COSMIC_IP_CUT;

  // Apply the containment cut to the starting positions of all
  // reconstructed tracks and showers. Pass this cut by default.
  sel_pfp_starts_in_PCV_ = true;

  // Loop over each PFParticle in the event
  for ( int p = 0; p < num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    // Use the track reconstruction results to get the start point for
    // every PFParticle for the purpose of verifying containment. We could
    // in principle differentiate between tracks and showers here, but
    // (1) we cut out all showers later on in the selection anyway, and
    // (2) the blinded PeLEE data ntuples do not include shower information.
    // We therefore apply the track reconstruction here unconditionally.
    float x = track_startx_->at( p );
    float y = track_starty_->at( p );
    float z = track_startz_->at( p );

    // Verify that the start of the PFParticle lies within the containment
    // volume.
    // TODO: revisit which containment volume to use for PFParticle start
    // positions. See https://stackoverflow.com/a/2488507 for an explanation
    // of the use of &= here. Don't worry, it's type-safe since both operands
    // are bool.
    sel_pfp_starts_in_PCV_ &= in_proton_containment_vol( x, y, z );
  }

  // Sets the sel_has_muon_candidate_ flag as appropriate. The threshold check
  // is handled later.
  this->find_muon_candidate();

  sel_nu_mu_cc_ = sel_reco_vertex_in_FV_ && sel_pfp_starts_in_PCV_
    && sel_has_muon_candidate_ && sel_topo_cut_passed_; // && sel_pandoraPDG_isnumu_;
}

// Sets the index of the muon candidate in the track vectors, or BOGUS_INDEX if
// one could not be found. The sel_has_muon_candidate_ flag is also set by this
// function.
void AnalysisEvent::find_muon_candidate() {

  std::vector<int> muon_candidate_indices;
  std::vector<int> muon_pid_scores;
  std::vector<int> muon_lengths;
  std::vector<int> muon_bdt_scores;

  for ( int p = 0; p < num_pf_particles_; ++p ) {
    // Only direct neutrino daughters (generation == 2) will be considered as
    // possible muon candidates
    unsigned int generation = pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    float track_score = pfp_track_score_->at( p );
    float start_dist = track_start_distance_->at( p );
    float track_length = track_length_->at( p );
    float pid_score = track_llr_pid_score_->at( p );
    float bdt_score = muon_BDT_score_->at( p ); 

    if ( track_score > MUON_TRACK_SCORE_CUT
      && start_dist < MUON_VTX_DISTANCE_CUT
      && track_length > MUON_LENGTH_CUT
      && pid_score > MUON_PID_CUT )
    {
      muon_candidate_indices.push_back( p );
      muon_pid_scores.push_back( pid_score );
      muon_lengths.push_back( track_length );
      muon_bdt_scores.push_back( bdt_score );	
    }
  }

  size_t num_candidates = muon_candidate_indices.size();
  if ( num_candidates > 0u ) sel_has_muon_candidate_ = true;

  if ( num_candidates == 1u ) {
    muon_candidate_pid_idx_ = muon_candidate_indices.front();
    muon_candidate_length_idx_ = muon_candidate_indices.front();
    muon_candidate_bdt_idx_ = muon_candidate_indices.front();
  }
  else if ( num_candidates > 1u ) {
    // In the case of multiple muon candidates, choose the one with the highest
    // PID score (most muon-like) and the longest length (will compare which method does best downstream)
    float highest_score = LOW_FLOAT;
    int chosen_index_pid = BOGUS_INDEX;
    float longest = LOW_FLOAT;
    int chosen_index_length = BOGUS_INDEX;
    float most_muon_like = BOGUS + 1. ;
    int chosen_index_bdt = BOGUS_INDEX;

    for ( size_t c = 0; c < num_candidates; ++c ) {
      float score = muon_pid_scores.at( c );
      float length = muon_lengths.at( c );
      float muon_bdt_score = muon_bdt_scores.at( c );

      if ( highest_score < score ) {
        highest_score = score;
        chosen_index_pid = muon_candidate_indices.at( c );
      }

      if ( longest < length ) {
        longest = length;
        chosen_index_length = muon_candidate_indices.at( c );
      }
  
      if (most_muon_like > muon_bdt_score ) {
	most_muon_like = muon_bdt_score;
	chosen_index_bdt = muon_candidate_indices.at( c );
      }
    }
    muon_candidate_pid_idx_ = chosen_index_pid;
    muon_candidate_length_idx_ = chosen_index_length;
    muon_candidate_bdt_idx_ = chosen_index_bdt;
  }
  else {
    muon_candidate_pid_idx_ = BOGUS_INDEX;
    muon_candidate_length_idx_ = BOGUS_INDEX;
    muon_candidate_bdt_idx_ = BOGUS_INDEX;
  }
}

void AnalysisEvent::find_pion_candidate() {

  float current_proton_bdt_score = 1.1; 
  int current_pion_candidate_idx_ = BOGUS_INDEX;

  for ( int p = 0; p < num_pf_particles_; ++p ) {
    // Only direct neutrino daughters (generation == 2) will be considered
    unsigned int generation = pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    // Skip particles already identified as muons
    if ( p == muon_candidate_pid_idx_) continue;

    // Skip particles with bogus track score 
    float track_length = track_length_->at( p );
    if (track_length <= 0. ) continue;

    float proton_bdt_score = proton_BDT_score_->at( p );

    // Skip particles for which BDT score could not be calculated
    // Use BOGUS - 1 to avoid compring floating point numbers
    if ( proton_bdt_score > BOGUS - 1. ) continue; 
    
    if ( proton_bdt_score < current_proton_bdt_score)
    {
      sel_has_pion_candidate_ = true; 
      current_pion_candidate_idx_ = p;
      current_proton_bdt_score = proton_bdt_score;	
    }
  }
 
  if (current_pion_candidate_idx_ != BOGUS_INDEX) pion_candidate_idx_ = current_pion_candidate_idx_;
}  
// Sets the analysis cut flags and decides whether the MC truth information
/// matches our signal definition
void AnalysisEvent::apply_selection() {

  // If we're working with an MC event, then categorize the event and set the
  // MC signal flags before proceeding with the selection. This function
  // keeps the category as kUnknown for real data.
  category_ = this->categorize_event();

  // Set sel_nu_mu_cc_ by applying those criteria
  this->apply_numu_CC_selection();
  this->get_multi_plane_Bragg_likelihood();
  // Count number of track like objects reconstructed
  // Fail the shower cut if any showers were reconstructed
  // NOTE: We could do this quicker like this,
  //   sel_no_reco_showers_ = ( num_showers_ > 0 );
  // but it might be nice to be able to adjust the track score for this cut.
  // Thus, we do it the hard way.
  int reco_shower_count = 0;
  int reco_track_count = 0;
  int n_non_proton_like = 0;

  for ( int p = 0; p < num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    
    unsigned int generation = pfp_generation_->at( p );
    if ( generation != 2u ) continue;


    float proton_score = proton_BDT_score_->at ( p );


    if ( proton_score < PROTON_BDT_CUT) ++n_non_proton_like;

    float tscore = pfp_track_score_->at( p );


    if ( tscore <= TRACK_SCORE_CUT ) ++reco_shower_count;
    else ++reco_track_count;
  }
 
  // Check the shower cut
  sel_no_reco_showers_ = ( reco_shower_count == 0 );

  // Check the track count. Need at least 3 for CC1piNp
  sel_min_3_tracks_ = ( reco_track_count >=3 );
  n_reco_tracks_ = reco_track_count;
  
  // Check we have 2 non proton like daugthers
  sel_2_non_proton_ = (n_non_proton_like == 2);
  n_non_proton_like_ = n_non_proton_like;
  
  // Set flags that default to true here
  sel_all_pfp_contained_ = true;
  // Set flags that default to false here
  sel_muon_contained_ = false;

  // If we have at least 3 tracks, 2 being non-proton like, find pion candidate
  // Only bother to do it if two above conditions are left to save computational time
  if ( sel_no_reco_showers_ && sel_min_3_tracks_ && sel_2_non_proton_ ) this->find_pion_candidate(); 

  // If no pion candidate found, skip to next event
  if ( !sel_has_pion_candidate_)
  {
    sel_CCNp1pi_ = false;
    return;
  }  
  for ( int p = 0; p < num_pf_particles_; ++p ) {

    // Only worry about direct neutrino daughters (PFParticles considered
    // daughters of the reconstructed neutrino)
    unsigned int generation = pfp_generation_->at( p );
    if ( generation != 2u ) continue; 
 
    float start_dist = track_start_distance_->at( p );
    if ( start_dist > PFP_DISTANCE_CUT ) sel_all_pfp_in_vtx_proximity_ = false;

    // Check all reco PFPs are contained
    float endx = track_endx_->at( p );
    float endy = track_endy_->at( p );
    float endz = track_endz_->at( p );
    bool end_contained = this->in_proton_containment_vol( endx, endy, endz );

    if( !(end_contained) )  sel_all_pfp_contained_ = false;
    
    // Check that we can find a muon candidate in the event. If more than
    // one is found, also fail the cut.
    if ( p == muon_candidate_pid_idx_ ) {
      // Check whether the muon candidate is contained. Use the same
      // containment volume as the protons. TODO: revisit this as needed.
      float endx = track_endx_->at( p );
      float endy = track_endy_->at( p );
      float endz = track_endz_->at( p );
      bool end_contained = this->in_proton_containment_vol( endx, endy, endz );

      if ( end_contained ) sel_muon_contained_ = true;

      // Check that the muon candidate is above threshold. Use the best
      // momentum based on whether it was contained or not.
      float muon_mom = LOW_FLOAT;
      float range_muon_mom = track_range_mom_mu_->at( p );
      float mcs_muon_mom = track_mcs_mom_mu_->at( p );

      if ( sel_muon_contained_ ) muon_mom = range_muon_mom;
      else muon_mom = mcs_muon_mom;

      if ( muon_mom >= MUON_P_MIN_MOM_CUT && muon_mom <= MUON_P_MAX_MOM_CUT ) {
        sel_muon_passed_mom_cuts_ = true;
      }

      // Apply muon candidate quality cut by comparing MCS and range-based
      // momentum estimators. Default to failing the cut.
      sel_muon_quality_ok_ = false;

      double frac_diff_range_mcs = std::abs( range_muon_mom - mcs_muon_mom );
      if ( range_muon_mom > 0. ) {
        frac_diff_range_mcs /= range_muon_mom;
        if ( frac_diff_range_mcs < MUON_MOM_QUALITY_CUT ) {
          sel_muon_quality_ok_ = true;
        }
      }

    }
  } 
  // All that remains is to apply the leading proton candidate cuts. We could
  // search for it above, but doing it here makes the code more readable (with
  // likely negligible impact on performance)
  this->find_lead_p_candidate();

  // if no proton candidate present, fail the selection and skip to next event 
  if ( !(sel_has_p_candidate_) ) 
  {
    sel_CCNp1pi_ = false;
    return; 
  }

  // Check the range-based reco momentum for the leading proton candidate
  float lead_p_KE = track_kinetic_energy_p_->at( lead_p_candidate_idx_ );
  float range_mom_lead_p = real_sqrt( lead_p_KE*lead_p_KE
    + 2.*PROTON_MASS*lead_p_KE );
  if ( range_mom_lead_p >= LEAD_P_MIN_MOM_CUT
    && range_mom_lead_p <= LEAD_P_MAX_MOM_CUT )
  {
    sel_lead_p_passed_mom_cuts_ = true;
  }

  // All right, we've applied all selection cuts. Set the flag that indicates
  // whether all were passed (and thus the event is selected as a CCNp0pi
  // candidate)
  sel_CCNp1pi_ = sel_nu_mu_cc_ && sel_no_reco_showers_
    && sel_muon_passed_mom_cuts_ && sel_muon_contained_ && sel_muon_quality_ok_
    && sel_has_p_candidate_ && sel_min_3_tracks_ && sel_2_non_proton_ && sel_has_pion_candidate_ 
    && sel_all_pfp_contained_ && sel_lead_p_passed_mom_cuts_ && sel_all_pfp_in_vtx_proximity_;
}

void AnalysisEvent::find_lead_p_candidate() {
  float lead_p_track_length = LOW_FLOAT;
  size_t lead_p_index = 0u;
  for ( int p = 0; p < num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    // Skip the muon and pion candidate reco track (this function assumes that it has
    // already been found)
    if ( p == muon_candidate_pid_idx_  || p == pion_candidate_idx_) continue;
     
    // Skip PFParticles that are shower-like (track scores near 0)
    float track_score = pfp_track_score_->at( p );
    //if ( track_score <= TRACK_SCORE_CUT ) continue;
    if ( track_score <= TRACK_SCORE_CUT ) continue;
    // All non-muon/pion-candidate reco tracks are considered proton candidates
    float track_length = track_length_->at( p );
    if ( track_length <= 0. ) continue;
    
    // All non-muon/pion-candidate reco tracks are considered proton candidates
    sel_has_p_candidate_ = true;
    
    if ( track_length > lead_p_track_length ) {
      lead_p_track_length = track_length;
      lead_p_index = p;
    }
  }

  // If the leading proton track length changed from its initial
  // value, then we found one. Set the index appropriately.
  if ( lead_p_track_length != LOW_FLOAT ) lead_p_candidate_idx_ = lead_p_index;
  // Otherwise, set the index to BOGUS_INDEX
  else lead_p_candidate_idx_ = BOGUS_INDEX;
}

void compute_GKI(const TVector3& p3mu, std::vector<TVector3>& p3p_v,/* const TVector3& p3p,*/ const TVector3& p3pi, 
  float& muE, float& pE, float& piE, float& nuE, float& pKE, float& delta_pT, float& delta_pL, float& delta_pL2,
  TVector3& q, float& pn, float& pn2,  float& delta_alpha3D, float& delta_phi3d_had, float& delta_phi3d_mu, float& delta_phi3d_p, float& delta_phi3d_pi)
{
  // Energies of particles 
  muE = std::sqrt(std::pow(MUON_MASS, 2) + p3mu.Mag2());
  //pE = std::sqrt(std::pow(PROTON_MASS, 2) + p3p.Mag2());
  piE = std::sqrt(std::pow(PI_PLUS_MASS, 2) + p3pi.Mag2());
  pE = 0.;
  pKE = 0.;

  
  //std::cout << p3p_v.size() << std::endl; 
  // MyPointer< std::vector<TVector3> >
  TVector3 p3p_total;
  //std::cout << "pE " <<  pE << std::endl;
  //std::cout << "pKE " <<  pKE << std::endl;
  //std::cout << "p3p_tptal " << p3p_total.X() << " " << p3p_total.Y() << " " << p3p_total.Z() << std::endl;
  for (int i = 0; i < p3p_v.size(); i++){
    //std::cout << i << std::endl;
    TVector3 temp_p3 = p3p_v.at(i);
    p3p_total += temp_p3;
    float temp_pE = std::sqrt( std::pow(PROTON_MASS, 2) + temp_p3.Mag2() ) ; 
    pE += temp_pE;
    pKE += temp_pE - PROTON_MASS;
    /*std::cout << i << std::endl;
    std::cout << "temp_pE " <<  temp_pE << std::endl;
    std::cout << "pE " <<  pE << std::endl;
    std::cout << "temp_pKE " <<  temp_pE - PROTON_MASS << std::endl;
    std::cout << "pKE " <<  pKE << std::endl;
    std::cout << "temp_p3 " << temp_p3.X() << " " << temp_p3.Y() << " " << temp_p3.Z() << std::endl;
    std::cout << "p3p_tptal " << p3p_total.X() << " " << p3p_total.Y() << " " << p3p_total.Z() << std::endl;
  */}
  //pKE = pE - PROTON_MASS;
  // Transverse momenta
  TVector3 pTmu(p3mu.X(), p3mu.Y(), 0.);
  TVector3 pTp(p3p_total.X(), p3p_total.Y(), 0.);
  TVector3 pTpi(p3pi.X(), p3pi.Y(), 0.);
  
  // Transverse missing momentum magnitude
  delta_pT = (pTmu + pTp + pTpi).Mag();
  // Longitudinal momenta
  TVector3 pLmu(0., 0., p3mu.Z());
  TVector3 pLp(0., 0., p3p_total.Z());
  TVector3 pLpi(0., 0.,  p3pi.Z());
  
  // Define R to simplify formula
  // Define mass of remnant nucleus
  //float R = TARGET_MASS + pLmu.Mag() + pLp.Mag() + pLpi.Mag() - muE - pE - piE;
  float R = TARGET_MASS + (pLmu + pLp + pLpi).Mag() - muE - pE - piE;
  float mf = TARGET_MASS - PROTON_MASS + MEAN_EXCITATION_ENERGY;

  // Longitudinal missing momentum magnitude from T2K paper (PhysRevD.103.112009) 
  delta_pL = 0.5 * R - 0.5 * ( ( std::pow(mf, 2) + std::pow(delta_pT, 2) ) / R); 
  
  // Compute E_cal estimator and longitudinal missing momentum magnitude from UB GKI paper to compare later
  nuE = muE + piE + pKE + MEAN_EXCITATION_ENERGY; 
  //delta_pL2 = pLmu.Mag() + pLp.Mag() + pLpi.Mag() - nuE;
  delta_pL2 = (pLmu + pLp + pLpi).Mag() - nuE;
  // Compute q momentum tranfer Q
  TVector3 nuEz(0., 0., nuE); 
  q = nuEz - p3mu;

  // Compute GKI 
  pn = std::sqrt( std::pow(delta_pL, 2) + std::pow(delta_pT, 2) );
  pn2 =  std::sqrt( std::pow(delta_pL2, 2) + std::pow(delta_pT, 2) );
  
  TVector3 pn_vec = p3mu + p3p_total + p3pi - nuEz;
  TVector3 p3had = p3p_total + p3pi;
  delta_alpha3D = std::acos ( ( q.X()*pn_vec.X() + q.Y()*pn_vec.Y() + q.Z()*pn_vec.Z() ) / ( q.Mag() * pn_vec.Mag()) ) ;
  delta_phi3d_had = std::acos ( ( q.X()*p3had.X() + q.Y()*p3had.Y() + q.Z()*p3had.Z() ) / ( q.Mag() * p3had.Mag()) );
  delta_phi3d_p = std::acos ( ( q.X()*p3p_total.X() + q.Y()*p3p_total.Y() + q.Z()*p3p_total.Z() ) / ( q.Mag() * p3p_total.Mag()) );
  delta_phi3d_pi = std::acos ( ( q.X()*p3pi.X() + q.Y()*p3pi.Y() + q.Z()*p3pi.Z() ) / ( q.Mag() * p3pi.Mag()) );
  delta_phi3d_mu = std::acos ( ( q.X()*p3mu.X() + q.Y()*p3mu.Y() + q.Z()*p3mu.Z() ) / ( q.Mag() * p3mu.Mag()) );

}
// Helper function for computing STVs (either reco or true)
void compute_stvs( const TVector3& p3mu, const TVector3& p3p, float& delta_pT,
  float& delta_phiT, float& delta_alphaT, float& delta_pL, float& pn,
  float& delta_pTx, float& delta_pTy )
{
  delta_pT = (p3mu + p3p).Perp();

  delta_phiT = std::acos( (-p3mu.X()*p3p.X() - p3mu.Y()*p3p.Y())
    / (p3mu.XYvector().Mod() * p3p.XYvector().Mod()) );

  TVector2 delta_pT_vec = (p3mu + p3p).XYvector();
  delta_alphaT = std::acos( (-p3mu.X()*delta_pT_vec.X()
    - p3mu.Y()*delta_pT_vec.Y())
    / (p3mu.XYvector().Mod() * delta_pT_vec.Mod()) );

  float Emu = std::sqrt(std::pow(MUON_MASS, 2) + p3mu.Mag2());
  float Ep = std::sqrt(std::pow(PROTON_MASS, 2) + p3p.Mag2());
  float R = TARGET_MASS + p3mu.Z() + p3p.Z() - Emu - Ep;

  // Estimated mass of the final remnant nucleus (CCQE assumption)
  float mf = TARGET_MASS - NEUTRON_MASS + BINDING_ENERGY;
  delta_pL = 0.5*R - (std::pow(mf, 2) + std::pow(delta_pT, 2)) / (2.*R);

  pn = std::sqrt( std::pow(delta_pL, 2) + std::pow(delta_pT, 2) );

  // Components of the 2D delta_pT vector (see arXiv:1910.08658)

  // We assume that the neutrino travels along the +z direction (also done
  // in the other expressions above)
  TVector3 zUnit( 0., 0., 1. );

  // Defines the x direction for the components of the delta_pT vector
  TVector2 xTUnit = zUnit.Cross( p3mu ).XYvector().Unit();

  delta_pTx = xTUnit.X()*delta_pT_vec.X() + xTUnit.Y()*delta_pT_vec.Y();

  // Defines the y direction for the components of the delta_T vector
  TVector2 yTUnit = ( -p3mu ).XYvector().Unit();

  delta_pTy = yTUnit.X()*delta_pT_vec.X() + yTUnit.Y()*delta_pT_vec.Y();
}

void AnalysisEvent::compute_observables() {

  // First compute the MC truth observables (if this is a signal MC event)
  this->compute_mc_truth_observables();

  // In cases where we failed to find a muon candidate, check whether there are
  // at least two generation == 2 PFParticles. If there are, then compute the
  // usual observables using the longest track as the muon candidate and the
  // second-longest track as the leading proton candidate. This will enable
  // sideband studies of NC backgrounds in the STV phase space.
  /*if ( !sel_has_muon_candidate_ ) {

    float max_trk_len = LOW_FLOAT;
    int max_trk_idx = BOGUS_INDEX;

    float next_to_max_trk_len = LOW_FLOAT;
    int next_to_max_trk_idx = BOGUS_INDEX;

    for ( int p = 0; p < num_pf_particles_; ++p ) {

      // Only include direct neutrino daughters (generation == 2)
      unsigned int generation = pfp_generation_->at( p );
      if ( generation != 2u ) continue;

      float trk_len = track_length_->at( p );

      if ( trk_len > next_to_max_trk_len ) {

        next_to_max_trk_len = trk_len;
        next_to_max_trk_idx = p;

        if ( next_to_max_trk_len > max_trk_len ) {

          next_to_max_trk_len = max_trk_len;
          next_to_max_trk_idx = max_trk_idx;

          max_trk_len = trk_len;
          max_trk_idx = p;
        }
      }
    }

    // If we found at least two usable PFParticles, then assign the indices to
    // be used below
    if ( max_trk_idx != BOGUS_INDEX && next_to_max_trk_idx != BOGUS_INDEX ) {
      muon_candidate_pid_idx_ = max_trk_idx;
      lead_p_candidate_idx_ = next_to_max_trk_idx;
    }
  }*/

  // Abbreviate some of the calculations below by using these handy
  // references to the muon and leading proton 3-momenta
  auto& p3mu = *p3_mu_;
  auto& p3p = *p3_lead_p_;
  auto& p3pi = *p3_cpi_;

  // Set the reco 3-momentum of the muon candidate if we found one
  bool muon = muon_candidate_pid_idx_ != BOGUS_INDEX;
  if ( muon ) {
    float mu_dirx = track_dirx_->at( muon_candidate_pid_idx_ );
    float mu_diry = track_diry_->at( muon_candidate_pid_idx_ );
    float mu_dirz = track_dirz_->at( muon_candidate_pid_idx_ );

    // The selection flag indicating whether the muon candidate is contained
    // was already set when the selection was applied. Use it to choose the
    // best momentum estimator to use.
    float muon_mom = LOW_FLOAT;
    if ( sel_muon_contained_ ) {
      muon_mom = track_range_mom_mu_->at( muon_candidate_pid_idx_ );
    }
    else {
      muon_mom = track_mcs_mom_mu_->at( muon_candidate_pid_idx_ );
    }

    p3mu = TVector3( mu_dirx, mu_diry, mu_dirz );
    p3mu = p3mu.Unit() * muon_mom;
  }

  bool pion = pion_candidate_idx_ != BOGUS_INDEX;
  if ( pion ) {
    float pi_dirx = track_dirx_->at( pion_candidate_idx_ );
    float pi_diry = track_diry_->at( pion_candidate_idx_ );
    float pi_dirz = track_dirz_->at( pion_candidate_idx_ );

    // The selection flag indicating whether the muon candidate is contained
    // was already set when the selection was applied. Use it to choose the
    // best momentum estimator to use.
    float pion_mom = LOW_FLOAT;
    float trk_length = track_length_->at( pion_candidate_idx_ ); 
    pion_mom =  a + b*trk_length - c*std::pow(trk_length, -1.*d); 
    
    p3pi = TVector3( pi_dirx, pi_diry, pi_dirz );
    p3pi = p3pi.Unit() * pion_mom;
  }

  // Set the reco 3-momentum of the leading proton candidate if we found one
  bool lead_p = lead_p_candidate_idx_ != BOGUS_INDEX;
  if ( lead_p ) {

    float p_dirx = track_dirx_->at( lead_p_candidate_idx_ );
    float p_diry = track_diry_->at( lead_p_candidate_idx_ );
    float p_dirz = track_dirz_->at( lead_p_candidate_idx_ );
    float KEp = track_kinetic_energy_p_->at( lead_p_candidate_idx_ );
    float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );

    p3p = TVector3( p_dirx, p_diry, p_dirz );
    p3p = p3p.Unit() * p_mom;
  }

  // Reset the vector of reconstructed proton candidate 3-momenta
  p3_p_vec_->clear();

  // Set the reco 3-momenta of all proton candidates (i.e., all generation == 2
  // tracks except the muon candidate) assuming we found both a muon candidate
  // pion candidate and at least one proton candidate.
  if ( muon && pion && lead_p ) {
    for ( int p = 0; p < num_pf_particles_; ++p ) {
      // Skip the muon and pion candidate
      if ( p == muon_candidate_pid_idx_ || pion_candidate_idx_) continue;

      // Only include direct neutrino daughters (generation == 2)
      unsigned int generation = pfp_generation_->at( p );
      if ( generation != 2u ) continue;

      float p_dirx = track_dirx_->at( p );
      float p_diry = track_diry_->at( p );
      float p_dirz = track_dirz_->at( p );
      float KEp = track_kinetic_energy_p_->at( p );
      float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );

      TVector3 p3_temp( p_dirx, p_diry, p_dirz );
      p3_temp = p3_temp.Unit() * p_mom;

      p3_p_vec_->push_back( p3_temp );
    }

    // TODO: reduce code duplication by just getting the leading proton
    // 3-momentum from this sorted vector
    // Sort the reco proton 3-momenta in order from highest to lowest magnitude
    std::sort( p3_p_vec_->begin(), p3_p_vec_->end(), [](const TVector3& a,
      const TVector3& b) -> bool { return a.Mag() > b.Mag(); } );
  } 

  // Compute reco STVs if we have both a muon candidate
  // and a leading proton candidate in the event
  if ( muon && lead_p && pion ) {
    //compute_stvs( p3mu, p3p, delta_pT_, delta_phiT_,
      //delta_alphaT_, delta_pL_, pn_, delta_pTx_, delta_pTy_ );
    compute_GKI( p3mu, *p3_p_vec_, p3pi, muE_, pE_, piE_, reco_Ecal_, pKE_,
      delta_pT_, delta_pL_, delta_pL2_, *q_, pn_, pn2_, delta_alpha3D_, delta_phi3d_had_, delta_phi3d_mu_, delta_phi3d_p_, delta_phi3d_pi_ );
    theta_mu_p_ = std::acos( p3mu.Dot(p3p) / p3mu.Mag() / p3p.Mag() );
    theta_mu_cpi_ = std::acos( p3mu.Dot(p3pi) / p3mu.Mag() / p3pi.Mag() );
  }
}

void AnalysisEvent::compute_mc_truth_observables() {

  // If this is not an MC event, then just return without doing anything
  if ( !is_mc_ ) return;

  size_t num_mc_daughters = mc_nu_daughter_pdg_->size();

  // Set the true 3-momentum of the final-state muon if there is one
  bool true_muon = ( mc_neutrino_is_numu_ && mc_nu_ccnc_ == CHARGED_CURRENT );
  if ( true_muon ) {
    // Loop over the MC neutrino daughters, find the muon, and get its
    // true 3-momentum. Note that we assume there is only one muon in
    // this loop.
    bool found_muon = false;
    for ( size_t d = 0u; d < num_mc_daughters; ++d ) {
      int pdg = mc_nu_daughter_pdg_->at( d );
      if ( pdg == MUON ) {
        found_muon = true;
        float px = mc_nu_daughter_px_->at( d );
        float py = mc_nu_daughter_py_->at( d );
        float pz = mc_nu_daughter_pz_->at( d );
        *mc_p3_mu_ = TVector3( px, py, pz );
        break;
      }
    }

    if ( !found_muon ) {
      std::cout << "WARNING: Missing muon in MC signal event!\n";
      return;
    }
  }

  // Reset the vector of true MC proton 3-momenta
  mc_p3_p_vec_->clear();
  mc_p3_cpi_vec_->clear(); 

  int mc_n_protons_temp = 0;
  // Set the true 3-momentum of the leading proton (if there is one)
  float max_mom_p = LOW_FLOAT;
  for ( int p = 0; p < num_mc_daughters; ++p ) {
    int pdg = mc_nu_daughter_pdg_->at( p );
    if ( pdg == PROTON )
    {

      mc_n_protons_temp += 1;
      float px = mc_nu_daughter_px_->at( p );
      float py = mc_nu_daughter_py_->at( p );
      float pz = mc_nu_daughter_pz_->at( p );
      TVector3 temp_p3 = TVector3( px, py, pz );

      mc_p3_p_vec_->push_back( temp_p3 );

      float mom = temp_p3.Mag();
      if ( mom > max_mom_p ) {
        max_mom_p = mom;
        *mc_p3_lead_p_ = temp_p3;
      }
    }
  }

  mc_n_protons_ = mc_n_protons_temp;
  // Three mom of leading pion i
  float max_mom_cpi = LOW_FLOAT;
  for ( int p = 0; p < num_mc_daughters; ++p){
    int pdg = mc_nu_daughter_pdg_->at(p);
    if ( std::abs(pdg) == PI_PLUS ){
      float px = mc_nu_daughter_px_->at( p );
      float py = mc_nu_daughter_py_->at( p );
      float pz = mc_nu_daughter_pz_->at( p );
      TVector3 temp_p3 = TVector3( px, py, pz );

      mc_p3_cpi_vec_->push_back( temp_p3 );
      float mom = temp_p3.Mag();
      if ( mom > max_mom_cpi ) {
      	max_mom_cpi = mom;
	*mc_p3_cpi_ = temp_p3;
      }
    }
  }

  // TODO: reduce code duplication by just getting the leading proton
  // 3-momentum from this sorted vector
  // Sort the true proton 3-momenta in order from highest to lowest magnitude
  std::sort( mc_p3_p_vec_->begin(), mc_p3_p_vec_->end(), [](const TVector3& a,
    const TVector3& b) -> bool { return a.Mag() > b.Mag(); } );

  // If the event contains a leading proton, then set the 3-momentum
  // accordingly
  bool true_lead_p = ( max_mom_p != LOW_FLOAT );
  if ( !true_lead_p && mc_is_signal_ ) {
    // If it doesn't for a signal event, then something is wrong.
    std::cout << "WARNING: Missing leading proton in MC signal event!\n";
    return;
  }

  // Compute angle between muon and pion if event contains a muon pion and a leading proton
  // Compute GKI 
  bool true_lead_cpi = ( max_mom_cpi != LOW_FLOAT ); 
  if (true_muon && true_lead_cpi && true_lead_p) {

    mc_theta_mu_cpi_ = std::acos( mc_p3_mu_->Dot(*mc_p3_cpi_)
      / mc_p3_mu_->Mag() / mc_p3_cpi_->Mag() );
  
    mc_theta_mu_p_ = std::acos( mc_p3_mu_->Dot(*mc_p3_lead_p_)
      / mc_p3_mu_->Mag() / mc_p3_lead_p_->Mag() );

    compute_GKI( *mc_p3_mu_, *mc_p3_p_vec_, *mc_p3_cpi_, mc_muE_, mc_pE_, mc_piE_, mc_Ecal_, mc_pKE_,
      mc_delta_pT_, mc_delta_pL_, mc_delta_pL2_, *mc_q_, mc_pn_, mc_pn2_, mc_delta_alpha3D_, mc_delta_phi3d_had_, mc_delta_phi3d_mu_, mc_delta_phi3d_p_, mc_delta_phi3d_pi_);

  }

}

void analyzer(const std::string& in_file_name,
 const std::string& output_filename)
{
  std::vector<std::string> in_files = { in_file_name };
  analyze( in_files, output_filename );
}

int main( int argc, char* argv[] ) {

  if ( argc != 3 ) {
    std::cout << "Usage: analyzer INPUT_PELEE_NTUPLE_FILE OUTPUT_FILE\n";
    return 1;
  }

  std::string input_file_name( argv[1] );
  std::string output_file_name( argv[2] );

  analyzer( input_file_name, output_file_name );

  return 0;
}
