version 1.0

# Purpose: build reference PCA axes and project the 05a cohort onto them.
# Phenotype: none. Unit: one cohort/sample_set; all study samples are retained.
# Inputs: study/1KG Phase 3 PGEN trios and official population panel.
# Both datasets must already use GRCh37 and the same forward strand.
# Outputs: PC1_SUM-PC20_SUM on the same projection-score scale, labels, model and QC.
# Scores are not rescaled to .eigenvec coordinates.
# Parameters: autosomal non-palindromic biallelic SNPs; geno <= 0.01 in each
# dataset; official second-degree relative exclusions; reference MAF >= 0.05; 20 PCs.
# Reference exclusions are applied before missingness, MAF, LD and PCA.
# Downloads PLINK2 deg2_phase3.king.cutoff.out.id (0.0884; not hg38).
# LD pruning: 500 kb windows, step 1, r2 <= 0.2, after long-LD exclusion.
# Require at least 1000 SNPs at matching, joint QC and final PCA stages.
# No HWE filtering, study-relative removal, strand guessing or ancestry assignment.

workflow AncestryProjection {
    input {
        File study_pgen
        File study_pvar
        File study_psam
        File reference_pgen
        File reference_pvar
        File reference_psam
        File population_panel
    }

    call ProjectOnto1kg {
        input:
            study_pgen = study_pgen, study_pvar = study_pvar, study_psam = study_psam,
            reference_pgen = reference_pgen, reference_pvar = reference_pvar,
            reference_psam = reference_psam, population_panel = population_panel
    }

    output {
        File study_pcs = ProjectOnto1kg.study_pcs
        File reference_pcs = ProjectOnto1kg.reference_pcs
        File reference_exclusions = ProjectOnto1kg.reference_exclusions
        File reference_samples = ProjectOnto1kg.reference_samples
        File pca_variants = ProjectOnto1kg.pca_variants
        File allele_weights = ProjectOnto1kg.allele_weights
        File allele_counts = ProjectOnto1kg.allele_counts
        File eigenvalues = ProjectOnto1kg.eigenvalues
        File matching_summary = ProjectOnto1kg.matching_summary
        File matching_sites = ProjectOnto1kg.matching_sites
        File study_duplicate_sites = ProjectOnto1kg.study_duplicate_sites
        File qc_summary = ProjectOnto1kg.qc_summary
        Array[File] plink_logs = ProjectOnto1kg.plink_logs
    }
}

task ProjectOnto1kg {
    input {
        File study_pgen
        File study_pvar
        File study_psam
        File reference_pgen
        File reference_pvar
        File reference_psam
        File population_panel
    }

    command <<<
        set -euo pipefail

        # Part 1. Preflight, download and small Python helpers
        echo "GRCh37 autosomal SNP PCA; geno <= 0.01 in each dataset"
        echo "Reference: official deg2 exclusions, MAF >= 0.05, long-LD exclusion, pruning 500kb 1 0.2"
        echo "20 PCs; all study samples retained; reference-mean imputation; no ancestry labels"
        echo "Resources: 2 CPUs, 16 GB RAM, 80 GB disk; minimum SNP count: 1000"
        command -v plink2 >/dev/null
        command -v python3 >/dev/null
        for input_file in "~{study_pgen}" "~{study_pvar}" "~{study_psam}" \
            "~{reference_pgen}" "~{reference_pvar}" "~{reference_psam}" "~{population_panel}"; do
            [[ -s "$input_file" ]] || { echo "ERROR: missing/empty input" >&2; exit 1; }
        done

        cat > qc_helpers.py <<'PY'
        import csv

        def rows(path):
            with open(path) as handle:
                for line in handle:
                    if line.strip() and not line.startswith('##'):
                        header = line.lstrip('#').split()
                        break
                else:
                    raise ValueError('Missing table header')
                for line in handle:
                    if line.strip():
                        fields = line.split()
                        if len(fields) != len(header):
                            raise ValueError('Table column count mismatch')
                        yield dict(zip(header, fields))

        def sample_id(row):
            return row.get('FID', '0'), row['IID']

        def write_table(path, header, records):
            with open(path, 'w') as handle:
                writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
                writer.writerow(header)
                writer.writerows(records)

        def variant_ids(path):
            return {row['ID'] for row in rows(path)}

        def position(row):
            chrom = row['CHROM'].removeprefix('chr')
            return chrom, int(row['POS'])

        def snp_id(row):
            chrom, pos = position(row)
            alleles = sorted([row['REF'], row['ALT']])
            return f'{chrom}:{pos}:{alleles[0]}:{alleles[1]}'

        def sample_table(path):
            records = list(rows(path))
            table = {sample_id(row): row for row in records}
            if len(table) != len(records):
                raise ValueError('Duplicate FID/IID in sample report')
            return table

        PY

        python3 -u - <<'PY'
        import hashlib
        import json
        import time
        import urllib.request
        from pathlib import Path
        from qc_helpers import rows, sample_id, write_table

        # Download only the small exclusion list; genotype inputs are already in Terra.
        remove_url = 'https://www.dropbox.com/s/zj8d14vv9mp6x3c/deg2_phase3.king.cutoff.out.id?dl=1'
        remove_file = 'deg2_phase3.king.cutoff.out.id'
        for attempt in range(3):
            try:
                with urllib.request.urlopen(remove_url, timeout=60) as response:
                    content = response.read()
                    expected = response.headers.get('Content-Length')
                if not content.strip() or (expected and len(content) != int(expected)):
                    raise ValueError('Empty or incomplete exclusion-list download')
                if content.lstrip().startswith(b'<'):
                    raise ValueError('Received HTML instead of the exclusion list')
                Path(remove_file).write_bytes(content)
                break
            except (OSError, ValueError):
                if attempt == 2:
                    raise
                print('Retrying exclusion-list download', flush=True)
                time.sleep(5)
        with open('qc_summary.json', 'w') as handle:
            json.dump({'reference_exclusion_url': remove_url,
                       'reference_exclusion_sha256': hashlib.sha256(content).hexdigest()}, handle)

        # Strip pedigree fields in temporary PSAMs so all retained references
        # contribute to frequency estimation and LD pruning. Preserve FID/IID.
        for prefix, path in [('study', '~{study_psam}'), ('reference', '~{reference_psam}')]:
            samples = list(rows(path))
            ids = [sample_id(row) for row in samples]
            if len(set(ids)) != len(ids) or any(row.get('SID', '0') != '0' for row in samples):
                raise ValueError('Expected unique FID/IID and no nonzero SID')
            write_table(prefix + '.ids.psam', ['#FID', 'IID'], ids)
        PY

        # Part 2. Select reference samples before variant QC
        # Official Phase 3 list: KING cutoff 0.0884 (14 exclusions).
        # https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3
        plink2 --psam reference.ids.psam --remove deg2_phase3.king.cutoff.out.id \
            --write-samples --threads 2 --memory 12000 --out reference.pruned

        python3 -u - <<'PY'
        from qc_helpers import rows, sample_id, write_table

        ids = [sample_id(row) for row in rows('reference.ids.psam')]
        kept = {sample_id(row) for row in rows('reference.pruned.id')}
        if not kept < set(ids):
            raise ValueError('Official exclusion list did not remove reference samples; check IDs/version')
        if len(kept) <= 20:
            raise ValueError('Too few unrelated reference samples for 20 PCs')
        panel = list(rows('~{population_panel}'))
        labels = {row['sample']: (row['pop'], row['super_pop']) for row in panel}
        if len(labels) != len(panel) or any(iid not in labels for _, iid in ids):
            raise ValueError('Reference IDs are not uniquely covered by the official panel')
        write_table('reference_samples.tsv', ['FID', 'IID', 'POP', 'SUPER_POP', 'KEPT'],
                    [(fid, iid, *labels[iid], int((fid, iid) in kept)) for fid, iid in ids])
        PY

        # Part 3. Match and filter SNPs
        ## Step 1. PLINK selects candidate SNPs; this PVAR is for matching only.
        # Never pair this filtered PVAR directly with the original study PGEN.
        plink2 --pvar "~{study_pvar}" --allow-extra-chr --autosome \
            --snps-only just-acgt --min-alleles 2 --max-alleles 2 --exclude-palindromic-snps \
            --make-just-pvar --threads 2 --memory 12000 --out study.candidates

        python3 -u - <<'PY'
        import json
        from collections import Counter
        from qc_helpers import rows, write_table, position, snp_id

        study = list(rows('~{study_pvar}'))
        counts = Counter(position(row) for row in study)
        # PLINK supplies eligible SNPs; reject unknown/same alleles still allowed by --snps-only.
        candidates = {position(row): snp_id(row) for row in rows('study.candidates.pvar')
                      if counts[position(row)] == 1 and row['REF'] in {'A', 'C', 'G', 'T'}
                      and row['ALT'] in {'A', 'C', 'G', 'T'} and row['REF'] != row['ALT']}
        reference_alleles = {}
        # Rewrite every PVAR row in its original order. Only IDs/chromosome labels
        # change: REF, ALT and record order must stay aligned with the input PGEN.
        with open('reference.named.pvar', 'w') as handle:
            handle.write('#CHROM\tPOS\tID\tREF\tALT\n')
            for i, row in enumerate(rows('~{reference_pvar}')):
                pos = position(row)
                name = f'reference_row_{i}'
                if pos in candidates:
                    reference_alleles.setdefault(pos, []).append((row['REF'], row['ALT']))
                    if snp_id(row) == candidates[pos]:
                        name = candidates[pos]
                handle.write(f"{pos[0]}\t{pos[1]}\t{name}\t{row['REF']}\t{row['ALT']}\n")
        # Classify unique study candidates; duplicate reference positions stay ambiguous.
        complement = str.maketrans('ACGT', 'TGCA')
        categories = ['direct', 'complement', 'incompatible', 'no_position', 'duplicate_reference']
        category_counts = Counter()
        matched = set()
        with open('matching_sites.tsv', 'w') as handle:
            handle.write('CHROM\tPOS\tID\tSTUDY_ALLELES\tREFERENCE_ALLELES\tCATEGORY\n')
            for pos, name in candidates.items():
                alleles = tuple(name.rsplit(':', 2)[1:])
                reference = reference_alleles.get(pos, [])
                if not reference:
                    category = 'no_position'
                elif len(reference) > 1:
                    category = 'duplicate_reference'
                elif tuple(sorted(reference[0])) == alleles:
                    category = 'direct'
                    matched.add(name)
                elif tuple(sorted(reference[0])) == tuple(sorted(a.translate(complement) for a in alleles)):
                    category = 'complement'
                else:
                    category = 'incompatible'
                category_counts[category] += 1
                ref_text = ';'.join('/'.join(pair) for pair in reference) or '.'
                handle.write(f"{pos[0]}\t{pos[1]}\t{name}\t{'/'.join(alleles)}\t{ref_text}\t{category}\n")
        total = len(candidates)
        write_table('matching_summary.tsv', ['CATEGORY', 'COUNT', 'DENOMINATOR', 'PROPORTION'],
                    [(category, category_counts[category], total,
                      category_counts[category] / total if total else 'NA') for category in categories])
        # All original records at duplicated study positions; excluded from the denominator.
        write_table('study_duplicate_sites.tsv', ['CHROM', 'POS', 'ID', 'REF', 'ALT'],
                    [(position(row)[0], row['POS'], row['ID'], row['REF'], row['ALT'])
                     for row in study if counts[position(row)] > 1])
        print('Matching denominator (unique eligible study SNPs):', total)
        for category in categories:
            print(category, category_counts[category])
        write_table('study.named.pvar', ['#CHROM', 'POS', 'ID', 'REF', 'ALT'],
                    [(position(row)[0], row['POS'], candidates.get(position(row), f'study_row_{i}'),
                      row['REF'], row['ALT']) for i, row in enumerate(study)])
        with open('matched.snps', 'w') as handle:
            handle.write('\n'.join(sorted(matched)) + '\n')

        with open('qc_summary.json') as handle:
            summary = json.load(handle)
        summary.update({'matched_snps': len(matched),
                   'study_duplicate_positions': sum(n > 1 for n in counts.values()),
                   'reference_duplicate_candidate_positions': category_counts['duplicate_reference'],
                   'matching_denominator': total,
                   'matching_counts': {category: category_counts[category] for category in categories},
                   'unmatched_unique_study_candidates': len(candidates) - len(matched)})
        with open('qc_summary.json', 'w') as handle:
            json.dump(summary, handle)
        if len(matched) < 1000:
            raise ValueError('Fewer than 1000 direct-matching SNPs; review matching reports')
        PY

        ## Step 2. Missingness on all study samples and retained reference samples.
        plink2 --pgen "~{study_pgen}" --pvar study.named.pvar --psam study.ids.psam \
            --allow-extra-chr --autosome --snps-only just-acgt --max-alleles 2 --exclude-palindromic-snps \
            --extract matched.snps --geno 0.01 --make-pgen \
            --threads 2 --memory 12000 --out study.qc
        plink2 --pgen "~{reference_pgen}" --pvar reference.named.pvar --psam reference.ids.psam \
            --allow-extra-chr --autosome --snps-only just-acgt --max-alleles 2 --exclude-palindromic-snps \
            --keep reference.pruned.id --extract matched.snps --geno 0.01 --make-pgen \
            --threads 2 --memory 12000 --out reference.qc

        python3 -u - <<'PY'
        from qc_helpers import variant_ids

        shared = variant_ids('study.qc.pvar') & variant_ids('reference.qc.pvar')
        if len(shared) < 1000:
            raise ValueError('Fewer than 1000 shared SNPs after missingness filtering')
        with open('shared.snps', 'w') as handle:
            handle.write('\n'.join(sorted(shared)) + '\n')
        PY

        ## Step 3. Reference MAF filtering, long-LD exclusion and LD pruning.
        # Frozen GRCh37 ranges, 1-based inclusive:
        # https://github.com/meyer-lab-cshl/plinkQC/blob/master/inst/extdata/high-LD-regions-hg19-GRCh37.txt
        cat > long_ld.txt <<'REGIONS'
        1 48227412 52227412 1
        2 86146488 101133568 2
        2 89563673 89564169 3
        2 98080299 98080324 4
        2 183291754 190291755 5
        3 47524995 50024996 6
        3 83417309 86917310 7
        5 44464242 50464243 8
        5 128972100 131972101 9
        6 25392020 33392022 10
        6 56892040 63942041 11
        6 139958306 142458307 12
        7 55032505 66362565 13
        8 7962589 11962591 14
        8 42880842 49837447 15
        8 111930823 114930824 16
        10 36959993 43679994 17
        11 87860351 90860352 18
        12 33108732 41713733 19
        20 32536338 35066586 20
        REGIONS
        plink2 --pfile reference.qc --extract shared.snps --maf 0.05 \
            --exclude range long_ld.txt --make-pgen --sort-vars \
            --threads 2 --memory 12000 --out reference.common
        plink2 --pfile reference.common --indep-pairwise 500kb 1 0.2 --indep-order 2 \
            --threads 2 --memory 12000 --out reference.ld
        plink2 --pfile reference.common --extract reference.ld.prune.in --make-pgen \
            --threads 2 --memory 12000 --out reference.pca_input

        python3 -u - <<'PY'
        from qc_helpers import variant_ids

        if len(variant_ids('reference.pca_input.pvar')) < 1000:
            raise ValueError('Fewer than 1000 PCA SNPs remain; review matching and QC')
        PY

        # Part 4. Reference PCA and projection with identical frequencies/weights
        plink2 --pfile reference.pca_input --freq counts \
            --pca 20 allele-wts meanimpute vcols=chrom,ref,alt \
            --threads 2 --memory 12000 --out reference.pca
        # Score sums avoid sample-specific denominators. Missing calls use reference means.
        for dataset in reference study; do
            input_prefix=study.qc
            if [[ "$dataset" == reference ]]; then input_prefix=reference.pca_input; fi
            plink2 --pfile "$input_prefix" --extract reference.ld.prune.in \
                --read-freq reference.pca.acount \
                --score reference.pca.eigenvec.allele 2 5 header-read variance-standardize list-variants cols=+scoresums \
                --score-col-nums 6-25 --missing sample-only \
                --threads 2 --memory 12000 --out "$dataset.projection"
        done

        # Part 5. Check projection and write unscaled score sums
        python3 -u - <<'PY'
        import json
        import math
        from collections import Counter
        from qc_helpers import rows, sample_id, sample_table, write_table, variant_ids

        with open('qc_summary.json') as handle:
            summary = json.load(handle)
        expected = variant_ids('reference.pca_input.pvar')
        pc_columns = [f'PC{k}_SUM' for k in range(1, 21)]
        reference_samples = list(rows('reference_samples.tsv'))
        labels = {sample_id(row): (row['POP'], row['SUPER_POP']) for row in reference_samples}
        sample_counts = {}
        for dataset in ['reference', 'study']:
            with open(dataset + '.projection.sscore.vars') as handle:
                scored = {line.strip() for line in handle if line.strip()}
            if scored != expected:
                raise ValueError('Projection skipped or added SNPs')
            scores = sample_table(dataset + '.projection.sscore')
            missing = sample_table(dataset + '.projection.smiss')
            source = 'reference.pca_input.psam' if dataset == 'reference' else 'study.ids.psam'
            if set(scores) != {sample_id(row) for row in rows(source)}:
                raise ValueError('Projection sample set changed')
            if set(missing) != set(scores):
                raise ValueError('Missingness/sample set mismatch')
            if dataset == 'reference' and set(sample_table('reference.pca.eigenvec')) != set(scores):
                raise ValueError('PCA/reference sample set mismatch')
            output = []
            for sid, row in scores.items():
                if not all(column in row for column in pc_columns):
                    raise ValueError('Missing one or more PC1_SUM-PC20_SUM columns')
                pcs = [float(row[column]) for column in pc_columns]
                if not all(math.isfinite(pc) for pc in pcs):
                    raise ValueError('Non-finite projected coordinate')
                f_miss = float(missing[sid]['F_MISS'])
                if not math.isfinite(f_miss) or not 0 <= f_miss <= 1:
                    raise ValueError('Invalid PCA SNP missingness rate')
                label = list(labels[sid]) if dataset == 'reference' else []
                output.append([*sid, *label, f_miss, *pcs])
            header = ['FID', 'IID'] + (['POP', 'SUPER_POP'] if dataset == 'reference' else [])
            write_table(dataset + '_pcs.tsv', header + ['PCA_SNP_F_MISS'] + pc_columns, output)
            sample_counts[dataset] = len(scores)

        study_snps = variant_ids('study.qc.pvar')
        reference_snps = variant_ids('reference.qc.pvar')
        summary.update({
            'study_samples_projected': sample_counts['study'],
            'reference_samples_before': len(reference_samples),
            'reference_samples_kept': sample_counts['reference'],
            'reference_populations_before': dict(Counter(row['POP'] for row in reference_samples)),
            'reference_populations_kept': dict(Counter(row['POP'] for row in reference_samples if row['KEPT'] == '1')),
            'study_snps_after_geno': len(study_snps),
            'reference_snps_after_geno': len(reference_snps),
            'shared_snps_after_geno': len(study_snps & reference_snps),
            'reference_snps_after_maf_and_long_ld': len(variant_ids('reference.common.pvar')),
            'pca_snps': len(expected),
            'pc_score_scale': 'PLINK2 PC1_SUM-PC20_SUM; no eigenvector rescaling',
            'parameters': {'build': 'GRCh37', 'geno': 0.01, 'king_cutoff': 0.0884,
                           'reference_maf': 0.05, 'ld': '500kb 1 0.2', 'pcs': 20,
                           'reference_exclusion_source': 'PLINK2 deg2_phase3.king.cutoff.out.id'},
            'ancestry_labels_assigned_to_study': False,
            'projection_shrinkage_corrected': False})
        with open('qc_summary.json', 'w') as handle:
            json.dump(summary, handle, indent=2)
            handle.write('\n')
        print('Reference PCA and study projection: complete')
        print('Study samples:', sample_counts['study'], '; PCA SNPs:', len(expected))
        PY
    >>>

    output {
        File study_pcs = "study_pcs.tsv"
        File reference_pcs = "reference_pcs.tsv"
        File reference_exclusions = "deg2_phase3.king.cutoff.out.id"
        File reference_samples = "reference_samples.tsv"
        File pca_variants = "reference.ld.prune.in"
        File allele_weights = "reference.pca.eigenvec.allele"
        File allele_counts = "reference.pca.acount"
        File eigenvalues = "reference.pca.eigenval"
        File matching_summary = "matching_summary.tsv"
        File matching_sites = "matching_sites.tsv"
        File study_duplicate_sites = "study_duplicate_sites.tsv"
        File qc_summary = "qc_summary.json"
        Array[File] plink_logs = glob("*.log")
    }

    runtime {
        docker: "ytian820/chip_cad_docker@sha256:fdee1d62eb43ee1e78a2ecf9df55e44a41af6b2b7eb48904ff1db19fa6c57963"
        cpu: 2
        memory: "16 GB"
        disks: "local-disk 80 SSD"
    }
}
