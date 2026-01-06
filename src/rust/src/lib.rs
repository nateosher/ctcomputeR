use ctcompute::{
    duration::types::EnrollmentRate, spending::types::SpendingFcn,
    trial::compute_trial::compute_trial, trial_characteristics::compute_ss_range::compute_ss_range,
};
use extendr_api::prelude::*;

/// Computes characteristics of a clinical trial with specified
/// power, sample size, etc.
/// @export
#[extendr]
fn ctcompute(
    n_patients: usize,
    alpha: f64,
    power: f64,
    maybe_lower_spending_fcn: Option<String>,
    maybe_upper_spending_fcn: Option<String>,
    maybe_look_fractions: Option<Vec<f64>>,
    prop_treated: f64,
    lambda_event_trt: f64,
    lambda_event_ctrl: f64,
    maybe_lambda_dropout: Option<f64>,
    enrollment_rates: Vec<f64>,
    enrollment_times: Vec<f64>,
    maybe_custom_alpha_spend: Option<Vec<f64>>,
    r: usize,
    tol: f64,
) -> List {
    let maybe_lower_spending_fcn = match (
        maybe_lower_spending_fcn.as_deref(),
        maybe_custom_alpha_spend.as_deref(),
    ) {
        (Some("LDOF"), _) => Some(SpendingFcn::LDOF),
        (Some("custom"), Some(custom_alpha_spend)) => Some(SpendingFcn::Custom {
            cumulative_spend: custom_alpha_spend.into(),
        }),
        (None, _) => None,
        (Some(unknown_spend), _) => {
            extendr_api::throw_r_error(format!("invalid spending function: `{}`", unknown_spend))
        }
    };
    let maybe_lower_spending_fcn_ref = maybe_lower_spending_fcn.as_ref();

    let maybe_upper_spending_fcn = match (
        maybe_upper_spending_fcn.as_deref(),
        maybe_custom_alpha_spend.as_deref(),
    ) {
        (Some("LDOF"), _) => Some(SpendingFcn::LDOF),
        (Some("custom"), Some(custom_alpha_spend)) => Some(SpendingFcn::Custom {
            cumulative_spend: custom_alpha_spend.into(),
        }),
        (None, _) => None,
        (Some(unknown_spend), _) => {
            extendr_api::throw_r_error(format!("invalid spending function: `{}`", unknown_spend))
        }
    };
    let maybe_upper_spending_fcn_ref = maybe_upper_spending_fcn.as_ref();

    let maybe_look_fractions_ref = maybe_look_fractions.as_ref();

    let enrollment_rate = match EnrollmentRate::new(enrollment_times, enrollment_rates) {
        Ok(er) => er,
        Err(e) => {
            rprintln!("");
            extendr_api::throw_r_error(&e.to_string());
        }
    };

    let ct = match compute_trial(
        n_patients,
        alpha,
        power,
        maybe_lower_spending_fcn_ref,
        maybe_upper_spending_fcn_ref,
        maybe_look_fractions_ref,
        prop_treated,
        lambda_event_trt,
        lambda_event_ctrl,
        maybe_lambda_dropout,
        &enrollment_rate,
        r,
        tol,
    ) {
        Ok(ct) => ct,
        Err(e) => {
            rprintln!("");
            extendr_api::throw_r_error(&e.to_string())
        }
    };

    list!(
        accrual_duration = ct.accrual_duration,
        trial_duration = ct.trial_duration,
        n_events = ct.n_events,
        n_patients = ct.n_patients,
        h0_expected_accrual_duration = ct.h0_expected_accrual_duration,
        h1_expected_accrual_duration = ct.h1_expected_accrual_duration,
        h0_expected_trial_duration = ct.h0_expected_trial_duration,
        h1_expected_trial_duration = ct.h1_expected_trial_duration,
        h0_expected_sample_size = ct.h0_expected_sample_size,
        h1_expected_sample_size = ct.h1_expected_sample_size,
    )
}

/// Computes range of sample sizes appropriate for given trial based on
/// heuristic of diminishing returns
/// @export
#[extendr]
pub fn ss_range(
    alpha: f64,
    power: f64,
    maybe_lower_spending_fcn: Option<String>,
    maybe_upper_spending_fcn: Option<String>,
    maybe_look_fractions: Option<Vec<f64>>,
    prop_treated: f64,
    lambda_event_trt: f64,
    lambda_event_ctrl: f64,
    maybe_lambda_dropout: Option<f64>,
    enrollment_rates: Vec<f64>,
    enrollment_times: Vec<f64>,
    maybe_custom_alpha_spend: Option<Vec<f64>>,
    tol: f64,
    delta: f64,
    min_perc_change: f64,
) -> extendr_api::Result<Vec<usize>> {
    let maybe_lower_spending_fcn = match (
        maybe_lower_spending_fcn.as_deref(),
        maybe_custom_alpha_spend.as_deref(),
    ) {
        (Some("LDOF"), _) => extendr_api::Result::Ok(Some(SpendingFcn::LDOF)),
        (Some("custom"), Some(custom_alpha_spend)) => {
            extendr_api::Result::Ok(Some(SpendingFcn::Custom {
                cumulative_spend: custom_alpha_spend.into(),
            }))
        }
        (None, _) => extendr_api::Result::Ok(None),
        _ => extendr_api::Result::Err("invalid spending function".into()),
    }?;
    let maybe_lower_spending_fcn_ref = maybe_lower_spending_fcn.as_ref();

    let maybe_upper_spending_fcn = match (
        maybe_upper_spending_fcn.as_deref(),
        maybe_custom_alpha_spend.as_deref(),
    ) {
        (Some("LDOF"), _) => extendr_api::Result::Ok(Some(SpendingFcn::LDOF)),
        (Some("custom"), Some(custom_alpha_spend)) => {
            extendr_api::Result::Ok(Some(SpendingFcn::Custom {
                cumulative_spend: custom_alpha_spend.into(),
            }))
        }
        (None, _) => extendr_api::Result::Ok(None),
        _ => extendr_api::Result::Err("invalid spending function".into()),
    }?;

    let maybe_upper_spending_fcn_ref = maybe_upper_spending_fcn.as_ref();

    let maybe_look_fractions_ref = maybe_look_fractions.as_ref();

    let enrollment_rate = match EnrollmentRate::new(enrollment_times, enrollment_rates) {
        Ok(er) => er,
        Err(e) => {
            rprintln!("");
            extendr_api::throw_r_error(&e.to_string());
        }
    };

    let ss_range_tup = match compute_ss_range(
        alpha,
        power,
        maybe_lower_spending_fcn_ref,
        maybe_upper_spending_fcn_ref,
        maybe_look_fractions_ref,
        prop_treated,
        lambda_event_trt,
        lambda_event_ctrl,
        maybe_lambda_dropout,
        &enrollment_rate,
        tol,
        delta,
        min_perc_change,
    ) {
        Ok(ss_tup) => ss_tup,
        Err(e) => {
            rprintln!("");
            extendr_api::throw_r_error(&e.to_string());
        }
    };

    Ok(vec![ss_range_tup.0, ss_range_tup.1])
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod ctcomputeR;
    fn ctcompute;
    fn ss_range;
}
