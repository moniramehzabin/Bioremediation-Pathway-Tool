"""Conservative annotation-seeded DyP family validation, independent of DIAMOND.

Annotations only nominate proteins for checking. Conventional, sequence-matched
family and both terminal-domain signatures are required for family support.
This is not a substrate-specific activity assignment.
"""
from collections import defaultdict
import hashlib
from pathlib import Path
import re


def normalized(value):
    return re.sub(r'\s+', ' ', value.strip().lower())


def nominated_dyp(proteins, active_reactions, resolved_pairs=()):
    if 'BRXN_DYP' not in active_reactions:
        return []
    products = {'dyp-type peroxidase', 'dye-decolorizing peroxidase',
                'dye decolorizing peroxidase', 'multifunctional dye peroxidase'}
    return sorted(tag for tag, p in proteins.items()
                  if (normalized(p.get('product', '')) in products
                      or normalized(p.get('gene', '')) in {'dyp', 'dyp2'})
                  and ('BRXN_DYP', tag) not in resolved_pairs)


def evaluate_dyp(tsv, sequences, nominated):
    rows = defaultdict(list)
    if tsv and Path(tsv).is_file():
        for line in Path(tsv).read_text(encoding='utf-8', errors='replace').splitlines():
            c = line.split('\t')
            if len(c) >= 12 and c[0] in nominated:
                rows[c[0]].append(c)
    output = []
    for tag in nominated:
        sequence = sequences[tag]
        digest = hashlib.md5(sequence.upper().encode()).hexdigest()
        valid = []
        mismatches = 0
        for c in rows[tag]:
            # InterProScan 6 uses column 10 to distinguish InterPro-N predictions.
            if c[9].strip().lower() == 'interpro-n':
                continue
            if c[1].lower() != digest or c[2] != str(len(sequence)):
                mismatches += 1
                continue
            try:
                start, end, score = int(c[6]), int(c[7]), float(c[8])
            except ValueError:
                continue
            if not (1 <= start <= end <= len(sequence)) or end-start+1 < 30:
                continue
            # Only conventional HMM E-values, not PROSITE scores or AI confidence.
            if c[3].lower() not in {'pfam', 'ncbifam', 'panther'} or not (0 <= score <= 1e-5):
                continue
            valid.append(c)
        nterm = [c for c in valid if c[3].lower() == 'pfam' and c[4] == 'PF04261']
        cterm = [c for c in valid if c[3].lower() == 'pfam' and c[4] == 'PF20628']
        family = [c for c in valid if c[11] == 'IPR006314' and
                  ((c[3].lower() == 'ncbifam' and c[4] == 'TIGR01413') or
                   (c[3].lower() == 'panther' and c[4].split(':')[0] == 'PTHR30521'))]
        ordered = any(int(n[7]) < int(c[6]) for n in nterm for c in cterm)
        supported = bool(family and ordered)
        matches = family + nterm + cterm
        detail = '; '.join(f'{c[3]}:{c[4]}:{c[6]}-{c[7]}:E={c[8]}' for c in matches)
        if supported:
            state = 'SUPPORTED_DYP_FAMILY'
        elif not rows[tag]:
            state = 'REVIEW_DOMAIN_RESULTS_UNAVAILABLE'
        elif mismatches and not valid:
            state = 'REVIEW_SEQUENCE_MISMATCH'
        else:
            state = 'REVIEW_INCOMPLETE_DYP_SIGNATURES'
        output.append({'locus_tag': tag, 'family_state': state,
                       'countable_family': supported, 'domain_evidence': detail,
                       'sequence_mismatch_rows': mismatches,
                       'interpretation': 'DyP family prediction only; activity against a tested dye is unproven.'})
    return output


def candidate_row(result, proteins, reactions):
    tag = result['locus_tag']; p = proteins[tag]; supported = result['countable_family']
    return dict(reaction_id='BRXN_DYP', gene_family=reactions['BRXN_DYP'].get('gene_family', ''),
        enzyme_name=reactions['BRXN_DYP'].get('enzyme_name', ''), locus_tag=tag,
        gene=p.get('gene', ''), product=p.get('product', ''), contig=p.get('contig', ''),
        reference_id='', reference_name='', reference_accession='', component_id='',
        component_gene_family='', identity_pct='', query_coverage_pct='',
        reference_coverage_pct='', evalue='', bitscore='', diamond_confidence='Not assessed',
        policy_interpro_status='SPECIFIC' if supported else 'NO_MATCH',
        policy_interpro_detail=result['domain_evidence'] or result['family_state'],
        genomic_context_support=False, genomic_context_detail='',
        pathway_context_support=False, pathway_context_fraction=0.0, pathway_context_detail='',
        adjudication='REVIEW', count_toward_pathway=supported, candidate_detected=True,
        exact_locus_supported=supported, confidence='Moderate' if supported else 'Review',
        evidence_warning='FAMILY_ONLY_SUBSTRATE_UNVALIDATED',
        evidence_source='GenBank nomination+InterPro family/domains' if supported else 'GenBank nomination only')
