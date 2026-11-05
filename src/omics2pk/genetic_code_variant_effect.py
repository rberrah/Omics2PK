"""Reference-checked DNA -> transcript -> genetic code -> assay interpretation.

This is a small, deterministic sequence CONSEQUENCE validator, not a genome
annotation engine, splice prediction tool, protein-function predictor, PGx
haplotype caller, clinical interpretation, or patient-specific PK model.

For genome-scale variant-to-transcript mapping use a pinned Ensembl VEP +
MANE/RefSeq/PharmVar source and preserve all transcript consequences. The
reference transcript and coordinates supplied here are CALLER-ATTESTED, not
independently retrieved. Splice effects and protein function require external
RNA/protein/activity evidence; variants need not affect every transcript.

Strictly 1-based genomic and 0-based spliced cDNA coordinates (forward
genomic alleles even for a negative-strand transcript); no silent HGVS ↔
genomic conversions. Standard NUCLEAR genetic code only (not mitochondria).
"""
from __future__ import annotations

from collections.abc import Mapping
from urllib.parse import urlparse

from .assay_aware_physio_triage import _classical_nmd_rule

class GeneticCodeEffectError(ValueError):
    pass

_ORDER="TCAG"
_AMINO="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TO_AA={a+b+c:aa for (a,b,c),aa in zip(
    ((a,b,c) for a in _ORDER for b in _ORDER for c in _ORDER),_AMINO)}
assert len(CODON_TO_AA)==64 and CODON_TO_AA["ATG"]=="M" and CODON_TO_AA["TGA"]=="*"
_COMPLEMENT=str.maketrans("ACGT","TGCA")
_ALLOWED_BUILDS={"GRCh37","GRCh38"}
_NMD_CLASSES={"CLASSICAL_EJC_NMD_PLAUSIBLE_NOT_CONFIRMED",
              "NMD_BOUNDARY_51_TO_54NT_INDETERMINATE",
              "CLASSICAL_EJC_NMD_NOT_PREDICTED_WITHIN_LAST_50NT",
              "CLASSICAL_EJC_NMD_NOT_PREDICTED_LAST_EXON_OR_AFTER_FINAL_JUNCTION"}

def _reqstr(v,key):
    if not isinstance(v,str) or not v.strip() or len(v)>2048:
        raise GeneticCodeEffectError("Explicit nonempty "+key+" required")
    return v.strip()

def _reference(rec):
    if not isinstance(rec,Mapping):
        raise GeneticCodeEffectError("Reference transcript metadata must be a mapping")
    transcript=_reqstr(rec.get("transcript_id"),"versioned transcript_id")
    if (not transcript.startswith(("NM_","ENST")) or
            (transcript.startswith("NM_") and "." not in transcript)
            or (transcript.startswith("ENST") and "." not in transcript)):
        raise GeneticCodeEffectError("Versioned NM_ or ENST transcript reference required")
    gene=_reqstr(rec.get("gene"),"HGNC gene").upper()
    assembly=rec.get("assembly")
    if assembly not in _ALLOWED_BUILDS:
        raise GeneticCodeEffectError("Explicit GRCh37/GRCh38 reference assembly required")
    chrom=_reqstr(rec.get("chromosome"),"chromosome")
    strand=rec.get("strand")
    if strand not in {"+","-"}:
        raise GeneticCodeEffectError("Original annotated transcript strand required")
    source=_reqstr(rec.get("reference_source_url"),"reference_source_url")
    if urlparse(source).scheme!="https":
        raise GeneticCodeEffectError("HTTPS reference source required")
    cdna=_reqstr(rec.get("spliced_cdna"),"original oriented transcript sequence").upper()
    if any(c not in "ACGT" for c in cdna):
        raise GeneticCodeEffectError("Only fully observed ACGT DNA reference accepted")
    exons=rec.get("exons")
    if not isinstance(exons,(list,tuple)) or len(exons)==0:
        raise GeneticCodeEffectError("Source transcript exon genomic intervals required")
    checked=[]
    for e in exons:
        if (not isinstance(e,(list,tuple)) or len(e)!=2 or
            any(not isinstance(z,int) or isinstance(z,bool) or z<1 for z in e) or e[1]<e[0]):
            raise GeneticCodeEffectError("Exon intervals require inclusive positive genomic 1-based start,end")
        checked.append((e[0],e[1]))
    if strand=="+" and any(checked[i][1]>=checked[i+1][0] for i in range(len(checked)-1)):
        raise GeneticCodeEffectError("Positive strand exons must be ordered, separated, non-overlapping")
    if strand=="-" and any(checked[i][0]<=checked[i+1][1] for i in range(len(checked)-1)):
        raise GeneticCodeEffectError("Negative strand exons must be ordered, separated, non-overlapping")
    if sum(e-s+1 for s,e in checked)!=len(cdna):
        raise GeneticCodeEffectError("Exon genomic sizes disagree with original spliced cDNA length")
    start=rec.get("cds_start_offset")
    end=rec.get("cds_end_offset")
    if any(not isinstance(v,int) or isinstance(v,bool) for v in (start,end)):
        raise GeneticCodeEffectError("Coding sequence requires 0-based spliced transcript CDS offsets")
    if not (0<=start<end<=len(cdna) and (end-start)%3==0):
        raise GeneticCodeEffectError("CDS interval and frame do not match original spliced cDNA")
    cds=cdna[start:end]
    if cds[:3]!="ATG" or CODON_TO_AA[cds[-3:]]!="*":
        raise GeneticCodeEffectError("Complete original nuclear CDS must have ATG and native terminal stop")
    if any(CODON_TO_AA[cds[i:i+3]]=="*" for i in range(0,len(cds)-3,3)):
        raise GeneticCodeEffectError("Reference CDS contains internal stop; review actual transcript")
    return {**rec,"gene":gene,"assembly":assembly,"chromosome":chrom,
            "strand":strand,"spliced_cdna":cdna,"exons":checked,
            "cds_start_offset":start,"cds_end_offset":end,
            "cds":cds}

def _genomic_coordinate_to_transcript(ref,pos):
    """Map single 1-based forward-genomic position to spliced oriented offset."""
    shift=0
    for s,e in ref["exons"]:
        if s<=pos<=e:
            return shift+(pos-s if ref["strand"]=="+" else e-pos)
        shift+=e-s+1
    return None

def _intron_location(ref,pos):
    exons=ref["exons"]
    for i in range(len(exons)-1):
        left,right=exons[i],exons[i+1]
        if ref["strand"]=="+":
            low,high=left[1],right[0]
        else:
            low,high=right[1],left[0]
        if not low<pos<high:
            continue
        donor=abs(pos-(left[1] if ref["strand"]=="+" else left[0]))
        acceptor=abs(pos-(right[0] if ref["strand"]=="+" else right[1]))
        if donor<=2:
            kind="splice_donor_variant"
        elif acceptor<=2:
            kind="splice_acceptor_variant"
        elif donor<=8 or acceptor<=20:
            kind="splice_region_variant"
        else:
            kind="intron_variant"
        return {"region":"intron","sequence_consequence":kind,
                "distance_to_donor_nt":donor,
                "distance_to_acceptor_nt":acceptor,
                "splicing_observed_in_RNA":False}
    return None

def _background(genotype):
    if not isinstance(genotype,Mapping):
        raise GeneticCodeEffectError("Individual DNA assay coverage must be declared")
    status=genotype.get("status")
    if status not in {"not_assessed","partial_assay","complete_target_assay"}:
        raise GeneticCodeEffectError("No inferred normal/wildtype genotype from absent calls")
    if status=="complete_target_assay":
        targets=genotype.get("covered_transcript_ids")
        if not isinstance(targets,(list,tuple)) or not targets or any(
            not isinstance(x,str) or not x for x in targets):
            raise GeneticCodeEffectError("Reference target coverage must be source-attested")
        url=_reqstr(genotype.get("assay_provenance_url"),"DNA assay provenance")
        if urlparse(url).scheme!="https":
            raise GeneticCodeEffectError("Complete-coverage source must be HTTPS")
    return {"original_genetic_assay_status":status,
            "no_variant_detected_is_not_proof_of_unmodified_function":True,
            "all_intronic_regulatory_and_structural_variants_excluded":False,
            "star_allele_diplotype_verified":False,
            "absence_of_called_variant_is_not_universal_reference_genome":True}

def annotate_coding_edit(*, ref, transcript_offset, ref_allele, alt_allele):
    """Consequence of transcript-oriented replacement, including coding indels.

    Exact REF matching mandatory; no splice consequences inferred from an
    exonic synonymous change. For genomic indels use externally normalized
    transcript HGVS, not this transcript-edit helper by row position alone.
    """
    r=_reference(ref)
    p=transcript_offset
    if not isinstance(p,int) or isinstance(p,bool) or p<0 or p>=len(r["spliced_cdna"]):
        raise GeneticCodeEffectError("0-based oriented original spliced transcript edit offset required")
    before=_reqstr(ref_allele,"transcript reference allele").upper()
    after=alt_allele.upper() if isinstance(alt_allele,str) else None
    if (after is None or len(before)>1000 or len(after)>1000 or
            any(x not in "ACGT" for x in before+after) or before==after):
        raise GeneticCodeEffectError("Distinct ACGT transcript REF/ALT with short explicit normalized edit required")
    if r["spliced_cdna"][p:p+len(before)]!=before:
        raise GeneticCodeEffectError("Original transcript reference allele does not match the supplied sequence")
    if p<r["cds_start_offset"] or p+len(before)>r["cds_end_offset"]:
        return {
            "sequence_consequence":("5_prime_UTR_variant" if p<r["cds_start_offset"]
                                    else "3_prime_UTR_variant"),
            "predicted_protein_amino_acid_change":None,
            "splice_effect_verified":False,
            "clinical_protein_function_estimated":False,
        }
    start=r["cds_start_offset"]
    cdspos=p-start
    original=r["cds"]
    mutated=original[:cdspos]+after+original[cdspos+len(before):]
    oldcodon=original[(cdspos//3)*3:(cdspos//3+1)*3]
    newcodon=mutated[(cdspos//3)*3:(cdspos//3+1)*3]
    if len(after)!=len(before):
        delta=len(after)-len(before)
        effect="frameshift_variant" if delta%3 else "inframe_insertion_or_deletion"
    else:
        if len(before)!=1:
            effect="multi_nucleotide_variant_needs_haplotype_review"
        else:
            oldaa=CODON_TO_AA[oldcodon]
            newaa=CODON_TO_AA[newcodon]
            if cdspos//3==0 and oldaa=="M" and newaa!="M":
                effect="start_lost"
            elif oldaa!="*" and newaa=="*":
                effect="stop_gained"
            elif oldaa=="*" and newaa!="*":
                effect="stop_lost"
            elif oldaa==newaa:
                effect="synonymous_variant"
            else:
                effect="missense_variant"
    premature=None
    if effect in {"stop_gained","frameshift_variant","inframe_insertion_or_deletion",
                  "multi_nucleotide_variant_needs_haplotype_review"}:
        for i in range((cdspos//3)*3,len(mutated)-2,3):
            if CODON_TO_AA[mutated[i:i+3]]=="*":
                # If this is the unchanged natural stop, not a new PTC.
                if i < len(original)-3 or len(mutated)!=len(original):
                    premature=start+i
                break
    last_junction=None
    if len(r["exons"])>=2:
        last_junction=sum(e-s+1 for s,e in r["exons"][:-1])
    # NMD classification is only meaningful for a genuinely new stop;
    # direct translation of exonic cDNA DOES NOT establish a splice outcome.
    if premature is not None and last_junction is not None:
        nmd=_classical_nmd_rule({
            "consequence":"premature_stop",
            "premature_stop_nt":premature,
            "last_exon_junction_nt":last_junction,
            "transcript_id":r["transcript_id"],
        })
    elif premature is not None:
        nmd="NO_DOWNSTREAM_EXON_JUNCTION_SINGLE_EXON_NMD_RULE_UNAVAILABLE"
    else:
        nmd="NO_SEQUENCE_CONFIRMED_NEW_PREMATURE_STOP"
    aa0=CODON_TO_AA[oldcodon]
    aa1=CODON_TO_AA[newcodon] if len(newcodon)==3 else None
    return {
        "sequence_consequence":effect,
        "reference_codon":oldcodon,
        "observed_alternate_codon":newcodon,
        "reference_aa":aa0,
        "alternate_aa":aa1,
        "aa_position_1based":cdspos//3+1,
        "reference_codon_translation_verified":True,
        "candidate_NMD_not_observed":nmd,
        "new_stop_cDNA_offset":premature,
        "last_exon_junction_cDNA_offset":last_junction,
        "splice_effect_verified":False,
        "clinical_protein_function_estimated":False,
        "predicted_protein_amino_acid_change":(
            f"{aa0}{cdspos//3+1}{aa1 or '?'}"
            if effect in {"missense_variant","synonymous_variant","stop_gained",
                           "stop_lost","start_lost"} else None),
        "protein_misfolding_trafficking_and_catalytic_activity_unknown":True,
    }

def annotate_genomic_snv(*, reference_transcript, variant, genome_assay):
    """Transcript-oriented SNV effect from exact forward-genome input.

    Genome assembly, chromosome, reference allele, strand and *original
    transcript sequence* verified for consistency. A one-gene/one-transcript
    result cannot replace VEP all-transcript consequences or PGx diplotypes.
    """
    r=_reference(reference_transcript)
    background=_background(genome_assay)
    if not isinstance(variant,Mapping):
        raise GeneticCodeEffectError("Original genomic variation record required")
    if variant.get("assembly")!=r["assembly"] or variant.get("chromosome")!=r["chromosome"]:
        raise GeneticCodeEffectError("Genome assembly and chromosome must match transcript")
    p=variant.get("genomic_position_1based")
    if not isinstance(p,int) or isinstance(p,bool) or p<=0:
        raise GeneticCodeEffectError("Exact 1-based genomic SNV position required")
    first=variant.get("genomic_ref_allele")
    alt=variant.get("genomic_alt_allele")
    if not isinstance(first,str) or not isinstance(alt,str) or first not in "ACGT" or alt not in "ACGT" or len(first)!=1 or len(alt)!=1 or first==alt:
        raise GeneticCodeEffectError("Original forward genomic single-base REF and ALT required")
    if variant.get("gene")!=r["gene"]:
        raise GeneticCodeEffectError("Gene cannot be silently assigned from genome position")
    original_url=_reqstr(variant.get("variant_source_url"),"source genomic variant")
    if urlparse(original_url).scheme!="https":
        raise GeneticCodeEffectError("Source genomic variant requires HTTPS provenance")
    idx=_genomic_coordinate_to_transcript(r,p)
    base={
        "gene":r["gene"],"transcript_id":r["transcript_id"],
        "assembly":r["assembly"],"chromosome":r["chromosome"],"strand":r["strand"],
        "original_source_attested_not_live_reference_verified":True,
        "real_patient_germline_or_somatic_status_not_inferred":True,
        "complete_star_allele_haplotype_not_inferred":True,
        **background,
    }
    if idx is None:
        intron=_intron_location(r,p)
        if intron is None:
            return {**base,"genomic_region":"outside_declared_transcript_intervals",
                    "sequence_consequence":"intergenic_or_uncovered_unknown",
                    "coding_aa_change":None,"functional_protein_measured":False}
        return {**base,"genomic_region":"intron",**intron,
                "coding_aa_change":None,
                "genomic_reference_base_checked_against_original_genome":False,
                "intronic_variant_alone_does_not_prove_splicing_disruption":True,
                "functional_protein_measured":False}
    t_ref=first if r["strand"]=="+" else first.translate(_COMPLEMENT)
    t_alt=alt if r["strand"]=="+" else alt.translate(_COMPLEMENT)
    if r["spliced_cdna"][idx]!=t_ref:
        raise GeneticCodeEffectError("Genomic reference allele conflicts with original transcript: strand/build/isoform?")
    seq=annotate_coding_edit(ref=r,transcript_offset=idx,ref_allele=t_ref,alt_allele=t_alt)
    return {**base,"genomic_region":"exon","transcript_offset_0based":idx,
            "transcript_reference_allele":t_ref,"transcript_alt_allele":t_alt,
            "genomic_reference_base_checked_against_original_genome":False,
            **seq}

def reconcile_curated_cyp3a5_allele(*, sequence_effect, source_variant_id,
                                   hgvs_c, star_allele, source_url,
                                   variant_call_complete=False):
    """Reconcile coding consequences with independently published CYP3A5 mechanisms.

    CYP3A5 *6 is synonymous by translation yet can disrupt splicing. Only
    defining variants are crosschecked; not a complete haplotype/diplotype.
    """
    from .functional_knowledge import CYP3A5_ALLELES
    variants={
        "*3":("rs776746","NM_000777.5:c.219-237A>G",
              {"intron_variant","splice_region_variant"},
              "CURATED_CRYPTIC_INTRON3_ACCEPTOR_ABERRANT_SPLICING"),
        "*6":("rs10264272","NM_000777.5:c.624G>A",
              {"synonymous_variant"},
              "CURATED_SYNONYMOUS_EXON7_CHANGE_ABERRANT_SPLICING"),
        "*7":("rs41303343","NM_000777.5:c.1035dup",
              {"frameshift_variant"},
              "CURATED_CODING_DUPLICATION_FRAMESHIFT_PREMATURE_STOP"),
    }
    if star_allele not in variants or star_allele not in CYP3A5_ALLELES:
        raise GeneticCodeEffectError("No qualified curated functional CYP3A5 variant for this identity")
    expected_id,expected_hgvs,expected_effect,mechanism=variants[star_allele]
    if (source_variant_id!=expected_id or hgvs_c!=expected_hgvs or
            sequence_effect.get("gene")!="CYP3A5" or
            sequence_effect.get("transcript_id")!="NM_000777.5" or
            sequence_effect.get("sequence_consequence") not in expected_effect or
            variant_call_complete is not True):
        raise GeneticCodeEffectError("CYP3A5 defining variant transcript/code effect conflict or incomplete call")
    if urlparse(_reqstr(source_url,"curated variant source")).scheme!="https":
        raise GeneticCodeEffectError("Curated CYP3A5 allele needs HTTPS evidence")
    return {
        "gene":"CYP3A5",
        "source_variant_id":source_variant_id,
        "source_transcript_hgvs":hgvs_c,
        "curated_star_allele":star_allele,
        "genetic_code_result":sequence_effect["sequence_consequence"],
        "independent_literature_mechanism":mechanism,
        "published_allele_function":"no_function",
        "in_this_subject_splicing_measured":False,
        "in_this_subject_protein_activity_measured":False,
        "full_diplotype_or_variant_phasing_established":False,
        "patient_specific_activity_loss_percent":None,
        "curated_source_url":source_url,
    }


def omic_implications(effect, *, genotype_status="partial_assay"):
    """Variant-consequence-specific, *nonquantitative* interpretation limits.

    'No variants seen' ≠ no variation: small panels cannot establish absence
    outside the assayed sites, and missense effects require functional assays.
    """
    if genotype_status not in {"not_assessed","partial_assay","complete_target_assay"}:
        raise GeneticCodeEffectError("Explicit genome assay coverage status required")
    consequence=effect.get("sequence_consequence")
    source={
        "missense_variant": (
            "ARN abundance can be normal despite altered folding, stability, localization or catalytic function",
            ["protein_variant_specific","protein_structure_or_localization","validated_functional_activity"]),
        "synonymous_variant": (
            "Same amino acid by genetic code; splice regulation, RNA stability or translation may still change",
            ["splice_junction_RNA","protein_abundance","validated_functional_activity"]),
        "stop_gained":(
            "A new stop is sequence-confirmed on this transcript; NMD and active protein loss remain unmeasured",
            ["isoform_RNA","variant_position_specific_peptides","intact_proteoform","validated_functional_activity"]),
        "stop_lost":(
            "Reference stop abolished; downstream alternative stop, readthrough and protein function require assay",
            ["full_length_transcript_RNA","intact_proteoform","validated_functional_activity"]),
        "frameshift_variant":(
            "Reading frame shifts; resulting protein/truncation and NMD require verified transcript consequences",
            ["isoform_RNA","peptide_mapping","intact_proteoform","validated_functional_activity"]),
        "start_lost":(
            "Canonical ATG changed; downstream alternative translation initiation cannot be ruled out",
            ["isoform_RNA","N_terminal_protein","validated_functional_activity"]),
        "inframe_insertion_or_deletion":(
            "Protein sequence length changes; function requires structure or activity evidence",
            ["protein_variant_specific","intact_proteoform","validated_functional_activity"]),
        "splice_donor_variant":(
            "Intron adjacent to splice donor: site effect hypothesis only; translated codon sequence cannot prove splice outcome",
            ["splice_junction_RNA","long_read_isoform_RNA","validated_functional_activity"]),
        "splice_acceptor_variant":(
            "Intron adjacent to splice acceptor: site effect hypothesis only; translated codon sequence cannot prove splice outcome",
            ["splice_junction_RNA","long_read_isoform_RNA","validated_functional_activity"]),
        "splice_region_variant":(
            "Near splice boundary but not canonical donor/acceptor: no automatic damaging call",
            ["splice_junction_RNA","long_read_isoform_RNA","validated_functional_activity"]),
        "intron_variant":(
            "Intronic region without direct protein codon change; may affect regulation/splicing",
            ["splice_junction_RNA","allele_specific_RNA","validated_functional_activity"]),
        "5_prime_UTR_variant":(
            "UTR sequence may alter regulatory and translation features; coding amino acid sequence unchanged",
            ["quantified_RNA","ribosome_or_translation_assay","protein_abundance"]),
        "3_prime_UTR_variant":(
            "UTR sequence may alter stability/localization; no direct coding amino acid change",
            ["quantified_RNA","RNA_stability_assay","protein_abundance"]),
    }
    message,assays=source.get(consequence,(
        "Unresolved genetic event: no protein function prediction",
        ["variant_annotation_and_reference_validation"]))
    return {
        "sequence_consequence":consequence,
        "genotype_status":genotype_status,
        "interpretation_hypothesis":message,
        "experimentally_informative_assay_classes":assays,
        "RNA_abundance_is_not_functional_protein":True,
        "mutation_causes_loss_of_function_proven":False,
        "mutation_absent_implies_wildtype_function":False,
        "no_predictive_omic_superiority_inferred":True,
        "functional_activity_and_drug_specificity_must_be_measured":True,
        "variant_annotation_requires_external_crosscheck":True,
    }
