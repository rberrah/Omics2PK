"""Synthetic-only fail-closed biomarker -> observed PK chronology tests."""
import copy
import unittest
from omics2pk.pk_link_contract import PKLinkError, qualify_clinical_link


def fixture():
    omics = [
        dict(sample_id="OM-RNA-01", subject_id="SYNTH-01",
             modality="RNA", collected_at="2026-01-01T10:00:00+01:00"),
        dict(sample_id="OM-RNA-02", subject_id="SYNTH-02",
             modality="RNA", collected_at="2026-01-02T07:00:00+01:00"),
    ]
    doses = [
        dict(dose_id="DO-01", subject_id="SYNTH-01", drug="tacrolimus",
             administered_at="2026-01-02T10:00:00+01:00", amount=2,
             dose_unit="mg", route="oral"),
        dict(dose_id="DO-02", subject_id="SYNTH-02", drug="tacrolimus",
             administered_at="2026-01-02T08:00:00+01:00", amount=3,
             dose_unit="mg", route="oral"),
    ]
    concentrations = [
        dict(observation_id="PK-01", subject_id="SYNTH-01", drug="tacrolimus",
             dose_id="DO-01", sampled_at="2026-01-02T12:00:00+01:00",
             concentration=8.0, concentration_unit="ng/mL"),
        dict(observation_id="PK-02", subject_id="SYNTH-02", drug="tacrolimus",
             dose_id="DO-02", sampled_at="2026-01-02T09:00:00+01:00",
             concentration=7.0, concentration_unit="ng/mL"),
    ]
    return omics, doses, concentrations


def qc(parts=None, modality="RNA"):
    return qualify_clinical_link(*(parts or fixture()),
                                 drug="tacrolimus", max_omics_lag_hours=48,
                                 selected_modality=modality)


class TemporalPKLinkTests(unittest.TestCase):
    def test_two_subjects_with_no_exported_identifiers(self):
        result = qc()
        self.assertEqual(result["n_dose_subjects_target"], 2)
        self.assertEqual(result["n_target_concentrations_linked"], 2)
        self.assertEqual(result["n_matched_assay_links"], 2)
        self.assertNotIn("SYNTH-01", repr(result))
        self.assertNotIn("OM-RNA-01", repr(result))
        self.assertEqual(result["status"], "TEMPORAL_IDENTITY_GATE_ONLY_NOT_POPPK_VALIDATION")

    def test_rejects_post_dose_transcriptome_as_predictor(self):
        parts = list(fixture())
        parts[0][0]["collected_at"] = "2026-01-02T11:00:00+01:00"
        with self.assertRaisesRegex(PKLinkError, "No pre-dose"):
            qc(parts)

    def test_rejects_mismatched_subject_dose_concentration(self):
        parts = list(fixture())
        parts[2][0]["subject_id"] = "SYNTH-02"
        with self.assertRaisesRegex(PKLinkError, "IDs disagree"):
            qc(parts)

    def test_rejects_ambiguous_latest_biological_sample(self):
        parts = list(fixture())
        extra = copy.deepcopy(parts[0][0])
        extra["sample_id"] = "NEW_TECH_REPLICATE"
        parts[0].append(extra)
        with self.assertRaisesRegex(PKLinkError, "Ambiguous"):
            qc(parts)

    def test_multimodal_requires_all_modalities_at_same_individual_not_necessarily_same_time(self):
        parts = list(fixture())
        with self.assertRaisesRegex(PKLinkError, "No pre-dose"):
            qc(parts, "MULTIOMICS")
        for subj, date in (("SYNTH-01", "2026-01-01T10:00:00+01:00"),
                           ("SYNTH-02", "2026-01-02T07:00:00+01:00")):
            for modality in ("PROTEIN", "METABOLITE"):
                parts[0].append(dict(sample_id=f"{subj}-{modality}", subject_id=subj,
                                     collected_at=date, modality=modality))
        result = qc(parts, "MULTIOMICS")
        self.assertEqual(result["n_matched_assay_links"], 6)

    def test_too_old_molecular_profile_does_not_link(self):
        with self.assertRaisesRegex(PKLinkError, "No pre-dose"):
            qualify_clinical_link(*fixture(), drug="tacrolimus",
                                  max_omics_lag_hours=2, selected_modality="RNA")

    def test_timezone_awareness_not_assumed(self):
        parts = list(fixture())
        parts[0][0]["collected_at"] = "2026-01-01T10:00:00"
        with self.assertRaisesRegex(PKLinkError, "time zone"):
            qc(parts)

    def test_does_not_accept_unmatched_dose_id(self):
        parts = list(fixture())
        parts[2][1]["dose_id"] = "UNKNOWN"
        with self.assertRaisesRegex(PKLinkError, "unknown"):
            qc(parts)

    def test_duplicate_concentration_id_rejected(self):
        parts = list(fixture())
        parts[2].append(copy.deepcopy(parts[2][0]))
        with self.assertRaisesRegex(PKLinkError, "Duplicate concentration"):
            qc(parts)

    def test_imputed_concentration_cannot_be_negative_or_nan(self):
        for value in (-1.0, 0, float("nan")):
            parts = list(fixture())
            parts[2][0]["concentration"] = value
            with self.assertRaisesRegex(PKLinkError, "strictly positive"):
                qc(parts)


if __name__ == "__main__":
    unittest.main()
