#pragma once

// Standard library includes
#include <cmath>
#include <limits>
#include <memory>

// ROOT includes
#include "TH1.h"
#include "Math/SpecFuncMathCore.h" // Needed for ROOT::Math::inc_gamma_c()

// STV analysis includes
#include "MatrixUtils.hh"
#include "SliceBinning.hh"
#include "SystematicsCalculator.hh"

class SliceHistogram {

  public:

    SliceHistogram() {}

    static SliceHistogram* make_slice_histogram( TH1D& reco_bin_histogram,
      const Slice& slice, const CovMatrix* input_cov_mat = nullptr );

    // TODO: revisit this implementation
    static SliceHistogram* make_slice_efficiency_histogram(
      const TH1D& true_bin_histogram, const TH2D& hist_2d, const Slice& slice );

    struct Chi2Result {
      Chi2Result( double chi2, int nbins, int dof, double pval )
        : chi2_( chi2 ), num_bins_( nbins ), dof_( dof ), p_value_( pval ) {}
      double chi2_;
      int num_bins_;
      int dof_;
      double p_value_;
    };

    Chi2Result get_chi2( const SliceHistogram& other );

    std::unique_ptr< TH1 > hist_;
    CovMatrix cmat_;
};

// Creates a new event histogram and an associated covariance matrix for a
// particular slice of phase space. The histogram is filled from the
// appropriate bin(s) of a 1D histogram of reco bin event counts. The mapping
// from reco bin number to the slice histogram bins is described by the input
// Slice object. Bin errors are set according to the reco-bin-space CovMatrix
// object pointed to by the input_cov_mat argument. If it is null, the bin
// errors are set to a default value of zero, and the output CovMatrix object
// owns a nullptr.
SliceHistogram* SliceHistogram::make_slice_histogram( TH1D& reco_bin_histogram,
  const Slice& slice, const CovMatrix* input_cov_mat )
{
  // Get the binning and axis labels for the current slice by cloning the
  // (empty) histogram owned by the Slice object
  TH1* slice_hist = dynamic_cast< TH1* >(
    slice.hist_->Clone("slice_hist") );

  slice_hist->SetDirectory( nullptr );

  // Fill the slice bins based on the input reco bins
  for ( const auto& pair : slice.bin_map_ ) {

    // One-based index for the global TH1 bin number in the slice
    int slice_bin_idx = pair.first;

    const auto& reco_bin_set = pair.second;

    double slice_bin_content = 0.;
    for ( const auto& rb_idx : reco_bin_set ) {
      // The ResponseMatrixMaker reco bin indices are zero-based, so I correct
      // for this here when pulling values from the one-based input ROOT
      // histogram
      slice_bin_content += reco_bin_histogram.GetBinContent( rb_idx + 1 );
    }

    slice_hist->SetBinContent( slice_bin_idx, slice_bin_content );

  } // slice bins

  // If we've been handed a non-null pointer to a CovMatrix object, then
  // we will use it to propagate uncertainties.
  TH2D* covmat_hist = nullptr;
  if ( input_cov_mat ) {

    // Create a new TH2D to hold the covariance matrix elements associated with
    // the slice histogram.
    // NOTE: I assume here that every slice bin is represented in the bin_map.
    // If this isn't the case, the bin counting will be off.
    // TODO: revisit this assumption and perhaps do something better
    int num_slice_bins = slice.bin_map_.size();

    covmat_hist = new TH2D( "covmat_hist", "covariance; slice bin;"
      " slice bin; covariance", num_slice_bins, 0., num_slice_bins,
      num_slice_bins, 0., num_slice_bins );
    covmat_hist->SetDirectory( nullptr );
    covmat_hist->SetStats( false );

    // We're ready. Populate the new covariance matrix using the elements
    // of the one for the reco bin space
    for ( const auto& pair_a : slice.bin_map_ ) {
      // Global slice bin index
      int sb_a = pair_a.first;
      // Set of reco bins that correspond to slice bin sb_a
      const auto& rb_set_a = pair_a.second;
      for ( const auto& pair_b : slice.bin_map_ ) {
        int sb_b = pair_b.first;
        const auto& rb_set_b = pair_b.second;

        double cov = 0.;
        const TH2D* cmat = input_cov_mat->cov_matrix_.get();
        for ( const auto& rb_m : rb_set_a ) {
          for ( const auto& rb_n : rb_set_b ) {
            // The covariance matrix TH2D uses one-based indices even though
            // the ResponseMatrixMaker numbering scheme is zero-based. I
            // correct for this here.
            cov += cmat->GetBinContent( rb_m + 1, rb_n + 1 );
          } // reco bin index m
        } // reco bin index n
        covmat_hist->SetBinContent( sb_a, sb_b, cov );
      } // slice bin index b
    } // slice bin index a


    // We have a finished covariance matrix for the slice. Use it to set
    // the bin errors on the slice histogram.
    for ( const auto& pair : slice.bin_map_ ) {

      int slice_bin_idx = pair.first;
      double bin_variance = covmat_hist->GetBinContent( slice_bin_idx,
        slice_bin_idx );
      double bin_error = std::sqrt( std::max(0., bin_variance) );

      // This works for a multidimensional slice because a global bin index
      // (as returned by TH1::GetBin) is used for slice_bin_idx.
      slice_hist->SetBinError( slice_bin_idx, bin_error );

    } // slice bins

  } // non-null input_cov_mat

  // We're done. Prepare the SliceHistogram object and return it.
  auto* result = new SliceHistogram;
  result->hist_.reset( slice_hist );
  result->cmat_.cov_matrix_.reset( covmat_hist );

  return result;
}

// TODO: revisit this rough draft. Right now, an assumption is made that the
// true and reco bins are defined in the same way with the same indices. This
// isn't enforced by the ResponseMatrixMaker configuration itself, although
// it is currently consistent with what you've done so far.
SliceHistogram* SliceHistogram::make_slice_efficiency_histogram(
  const TH1D& true_bin_histogram, const TH2D& hist_2d, const Slice& slice )
{
  // Get the binning and axis labels for the current slice by cloning the
  // (empty) histogram owned by the Slice object
  TH1* slice_hist = dynamic_cast< TH1* >(
    slice.hist_->Clone("slice_hist") );

  slice_hist->SetDirectory( nullptr );
  slice_hist->GetYaxis()->SetTitle( "efficiency" );
  slice_hist->GetYaxis()->SetRangeUser( 0., 1. );

  // Fill the slice bins based on the input reco bins
  for ( const auto& pair : slice.bin_map_ ) {

    // One-based index for the global TH1 bin number in the slice
    int slice_bin_idx = pair.first;

    const auto& reco_bin_set = pair.second;

    double selected_signal_evts = 0.;
    double all_signal_evts = 0.;
    for ( const auto& rb_idx : reco_bin_set ) {
      // The ResponseMatrixMaker reco bin indices are zero-based, so I correct
      // for this here when pulling values from the one-based input ROOT
      // histogram.
      all_signal_evts += true_bin_histogram.GetBinContent( rb_idx + 1 );

      // Include selected signal events in the current true bin that fall into
      // any of the reco bins
      selected_signal_evts += hist_2d.Integral( rb_idx + 1, rb_idx + 1,
        1, hist_2d.GetNbinsY() );
    }

    double bin_efficiency = selected_signal_evts / all_signal_evts;
    // See DocDB #32401, Eq. (5.2)
    double bin_stat_err = std::sqrt( std::max(0., bin_efficiency
      * (1. - bin_efficiency) / all_signal_evts) );
    slice_hist->SetBinContent( slice_bin_idx, bin_efficiency );
    slice_hist->SetBinError( slice_bin_idx, bin_stat_err );

  } // slice bins

  TH2D* covmat_hist = nullptr;

  // We're done. Prepare the SliceHistogram object and return it.
  auto* result = new SliceHistogram;
  result->hist_.reset( slice_hist );
  result->cmat_.cov_matrix_.reset( nullptr );

  return result;
}

SliceHistogram::Chi2Result SliceHistogram::get_chi2(
  const SliceHistogram& other )
{
  int num_bins = hist_->GetNbinsX();
  if ( other.hist_->GetNbinsX() != num_bins ) {
    throw std::runtime_error( "Incompatible vector sizes in chi^2"
      " calculation" );
  }

  int my_cov_mat_x_bins = cmat_.cov_matrix_->GetNbinsX();
  int my_cov_mat_y_bins = cmat_.cov_matrix_->GetNbinsY();
  int other_cov_mat_x_bins = other.cmat_.cov_matrix_->GetNbinsY();
  int other_cov_mat_y_bins = other.cmat_.cov_matrix_->GetNbinsY();
  if ( my_cov_mat_x_bins != num_bins
    || my_cov_mat_y_bins != num_bins
    || other_cov_mat_x_bins != num_bins
    || other_cov_mat_y_bins != num_bins )
  {
    throw std::runtime_error( "Invalid covariance matrix dimensions"
      " encountered in chi^2 calculation" );
  }

  // The total covariance matrix on the difference between the
  // two histograms is just the sum of each individual SliceHistogram's
  // owned covariance matrix.
  CovMatrix cov_mat;
  cov_mat += cmat_;
  cov_mat += other.cmat_;

  // Get access to a TMatrixD object representing the covariance matrix.
  auto cov_matrix = cov_mat.get_matrix();

  // Invert the covariance matrix
  auto inverse_cov_matrix = invert_matrix( *cov_matrix );

  // Create a 1D vector containing the difference between the two slice
  // histograms in each bin
  TMatrixD diff_vec( 1, num_bins );
  for ( int a = 0; a < num_bins; ++a ) {
    // Note the one-based bin indices used for ROOT histograms
    double counts = hist_->GetBinContent( a + 1 );
    double other_counts = other.hist_->GetBinContent( a + 1 );
    diff_vec( 0, a ) = counts - other_counts;
  }

  // Multiply diff * covMat^{-1} * diff^{T} to get chi-squared
  TMatrixD temp1( diff_vec, TMatrixD::kMult, *inverse_cov_matrix );
  TMatrixD temp2( temp1, TMatrixD::kMult, diff_vec.T() );

  // We'll now have a 1x1 matrix containing the chi-squared value
  double chi2 = temp2( 0, 0 );

  // Assume that parameter fitting is not done, so that the relevant degrees
  // of freedom for the chi^2 test is just the number of bins minus one
  int dof = num_bins - 1;

  // Calculate a p-value for observing a chi^2 value at least as large as the
  // one actually obtained
  double p_value = ROOT::Math::inc_gamma_c( dof / 2., chi2 / 2. );

  Chi2Result result( chi2, num_bins, dof, p_value );
  return result;

}