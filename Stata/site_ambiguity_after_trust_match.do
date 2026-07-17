* =============================================================================
* site ambiguity after trust-matching
* -----------------------------------------------------------------------------
* Amanda's question, from her email of 17 July: once cosd rows whose trust does
* not match the registry's diagnosis_trust are thrown out, are the patients left
* with more than one cosd row mostly disagreeing on trust anyway - which no
* amount of snomed work can fix - or are they the same trust submitting more
* than one 5-character site code, which is the case snomed (or morphology, or
* just picking one) is actually needed for.
*
* This does not build anything or save a file. It only answers the question, so
* it can be run on its own against the raw extracts.
*
* There is an R version of this same check at
* R/reference/12_site_ambiguity_after_trust_match.R. The two should agree on
* every headline number; if they do not, one of them is reading a field
* differently to the other and that is worth chasing down before trusting
* either.
*
* Edit the two paths below, then run the whole file.
* =============================================================================

clear all
set more off

local path_rapid "W:\_DATA\IainTimmins\2026 OG SOTN data\20260212_Rapidtumour_linked_2026SOTN_clean_OG_postPT.dta"
local path_cosd  "W:\_DATA\IainTimmins\2026 OG SOTN data\20260212_all_cosddiagnosis_rapid_202601_OG.dta"

* -----------------------------------------------------------------------------
* the registry's trust for each patient
* -----------------------------------------------------------------------------
use patient_pseudo_id diagnosis_trust using "`path_rapid'", clear
gen reg_trust3 = upper(substr(trim(diagnosis_trust), 1, 3))
keep patient_pseudo_id reg_trust3
duplicates drop patient_pseudo_id, force
tempfile reg
save "`reg'"

* -----------------------------------------------------------------------------
* every cosd row offering a usable site code, joined to the registry's trust
* -----------------------------------------------------------------------------
use patient_pseudo_id site_code_of_diagnosis using "`path_cosd'", clear
gen site = upper(trim(site_code_of_diagnosis))
drop if site == "" | length(site) != 5
gen site_trust3 = substr(site, 1, 3)

merge m:1 patient_pseudo_id using "`reg'"
* a cosd row for a patient the registry extract does not have says nothing
* about this question, so it is set aside rather than counted either way
drop if _merge == 1
drop _merge

count
display as text "cosd rows with a five-character site code: " r(N)

* -----------------------------------------------------------------------------
* before amanda's filter: how many patients offer more than one code at all
* -----------------------------------------------------------------------------
preserve
    duplicates drop patient_pseudo_id site, force
    bysort patient_pseudo_id: gen n_codes_before = _N
    by patient_pseudo_id: keep if _n == 1
    gen has_choice = n_codes_before > 1

    display as text _newline "before the trust filter:"
    tab has_choice
restore

* -----------------------------------------------------------------------------
* amanda's filter: keep only rows whose trust matches the registry
* -----------------------------------------------------------------------------
count if site_trust3 == reg_trust3
display as text _newline "rows left after dropping trust mismatches: " r(N)

keep if site_trust3 == reg_trust3
duplicates drop patient_pseudo_id site, force
tempfile matched_pairs
save "`matched_pairs'"

bysort patient_pseudo_id: gen n_codes = _N
by patient_pseudo_id: gen first = (_n == 1)

count if first & n_codes > 1
display as text _newline r(N) ///
    " patients still offer more than one site code from the same trust -" ///
    _newline "different hospital, same trust, or a coding slip within one hospital."

preserve
    keep if first
    gen has_choice = n_codes > 1

    display as text _newline "of the patients with a trust-matched cosd row:"
    tab has_choice

    display as text _newline "how many different codes these patients choose between:"
    tab n_codes if has_choice == 1
restore

* -----------------------------------------------------------------------------
* a handful of ambiguous patients, to look at directly
* -----------------------------------------------------------------------------
preserve
    use "`matched_pairs'", clear
    bysort patient_pseudo_id: gen n_codes = _N
    keep if n_codes > 1
    egen pid_group = group(patient_pseudo_id)
    keep if pid_group <= 5
    sort patient_pseudo_id site

    display as text _newline "a sample, to look at directly:"
    list patient_pseudo_id site, noobs sepby(patient_pseudo_id)
restore

display as text _newline ///
    "the gap between the two 'has_choice' tables above is what trust-matching" ///
    _newline "alone settles. what's left in the second one is the true test of" ///
    _newline "whether snomed, morphology, or picking the more common code, is" ///
    _newline "pulling its weight."
