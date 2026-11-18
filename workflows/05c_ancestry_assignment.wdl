version 1.0

# Purpose: assign reference-based genetic ancestry labels to the 05b cohort.
# Phenotype: none. Unit: participant; no study samples are removed.
# Inputs: 05b reference/study PC tables on the same PC_SUM scale.
# Outputs: assignments, validation predictions/report, PCA plots and fitted model.
# Parameters: PC1_SUM-PC10_SUM, 500 trees, seed 42, stratified 80/20 split,
# and maximum classification probability >= 0.90; otherwise unassigned.
# Execution only: input integrity and scientific review are performed manually.
# Validation assesses the classifier in the existing reference PCA space;
# PCA itself was fitted to all retained reference samples in 05b.
# Probabilities are classifier scores, not ancestry proportions.

workflow AncestryAssignment {
    input {
        File reference_pcs
        File study_pcs
    }

    call AssignAncestry {
        input: reference_pcs = reference_pcs, study_pcs = study_pcs
    }

    output {
        File ancestry_assignments = AssignAncestry.ancestry_assignments
        File validation_predictions = AssignAncestry.validation_predictions
        File validation_report = AssignAncestry.validation_report
        File ancestry_pca = AssignAncestry.ancestry_pca
        File ancestry_model = AssignAncestry.ancestry_model
    }
}

task AssignAncestry {
    input {
        File reference_pcs
        File study_pcs
    }

    command <<<
        set -euo pipefail
        python -m pip install --quiet --no-cache-dir \
            numpy==1.26.4 scipy==1.13.1 pandas==2.2.3 scikit-learn==1.5.2 \
            matplotlib==3.9.2 joblib==1.4.2 threadpoolctl==3.5.0

        python - <<'PY'
        # Part 1. Configuration and input tables
        import json
        import platform
        import joblib
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        import sklearn
        from matplotlib.backends.backend_pdf import PdfPages
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.metrics import classification_report, confusion_matrix
        from sklearn.model_selection import train_test_split

        pc_columns = [f'PC{k}_SUM' for k in range(1, 11)]
        groups = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
        min_probability = 0.90
        seed = 42
        test_fraction = 0.20
        model_parameters = dict(n_estimators=500, random_state=seed, n_jobs=2)
        print('Ancestry assignment: PC1_SUM-PC10_SUM; 500 trees; seed 42')
        print('Reference split: stratified 80/20; assignment threshold: 0.90')
        print('All study samples retained; no automatic QC exclusions')

        reference = pd.read_csv('~{reference_pcs}', sep='\t', dtype={'FID': str, 'IID': str})
        study = pd.read_csv('~{study_pcs}', sep='\t', dtype={'FID': str, 'IID': str})

        def predict_table(model, samples):
            probabilities = model.predict_proba(samples[pc_columns])
            best = probabilities.argmax(axis=1)
            result = samples[['FID', 'IID', 'PCA_SNP_F_MISS']].copy()
            result['predicted_group'] = model.classes_[best]
            result['max_probability'] = probabilities.max(axis=1)
            result['assigned_group'] = np.where(
                result['max_probability'] >= min_probability,
                result['predicted_group'], 'unassigned')
            for index, group in enumerate(model.classes_):
                result['prob_' + group] = probabilities[:, index]
            return result

        # Part 2. Held-out reference classification
        train, validation = train_test_split(
            reference, test_size=test_fraction, random_state=seed,
            stratify=reference['SUPER_POP'])
        model = RandomForestClassifier(**model_parameters)
        model.fit(train[pc_columns], train['SUPER_POP'])
        predictions = predict_table(model, validation)
        predictions['true_group'] = validation['SUPER_POP']
        predictions.to_csv('validation_predictions.tsv', sep='\t', index=False)

        accepted = predictions['assigned_group'] != 'unassigned'
        matrix_labels = groups + ['unassigned']
        report = {
            'pc_columns': pc_columns,
            'model_parameters': model.get_params(),
            'test_fraction': test_fraction,
            'split_seed': seed,
            'min_probability': min_probability,
            'reference_training_samples': len(train),
            'reference_validation_samples': len(validation),
            'validation_scope': 'Classifier holdout; reference PCA was fitted in 05b using all reference samples',
            'classification_before_threshold': classification_report(
                predictions['true_group'], predictions['predicted_group'],
                labels=groups, output_dict=True, zero_division=0),
            'confusion_matrix_labels': matrix_labels,
            'confusion_matrix_rows_true_columns_assigned': confusion_matrix(
                predictions['true_group'], predictions['assigned_group'],
                labels=matrix_labels).tolist(),
            'validation_unassigned_fraction': float((~accepted).mean()),
            'validation_accepted_accuracy': float((
                predictions.loc[accepted, 'assigned_group'] ==
                predictions.loc[accepted, 'true_group']).mean()) if accepted.any() else None,
            'software': {'python': platform.python_version(), 'scikit_learn': sklearn.__version__,
                         'numpy': np.__version__, 'pandas': pd.__version__,
                         'matplotlib': matplotlib.__version__, 'joblib': joblib.__version__}
        }

        # Part 3. Refit on all reference samples and assign study labels
        model.fit(reference[pc_columns], reference['SUPER_POP'])
        assignments = predict_table(model, study)
        assignments.to_csv('ancestry_assignments.tsv', sep='\t', index=False)
        report['reference_samples_final_fit'] = len(reference)
        report['study_samples'] = len(study)
        report['study_assignment_counts'] = assignments['assigned_group'].value_counts().reindex(
            matrix_labels, fill_value=0).to_dict()
        with open('validation_report.json', 'w') as handle:
            json.dump(report, handle, indent=2, allow_nan=False)
            handle.write('\n')
        joblib.dump({'model': model, 'pc_columns': pc_columns,
                     'min_probability': min_probability, 'software': report['software']},
                    'ancestry_model.joblib')

        # Part 4. Reference/study overlays on the unchanged PC_SUM scale
        colors = dict(zip(matrix_labels,
                          ['#0072B2', '#E69F00', '#009E73', '#CC79A7', '#D55E00', '#000000']))
        with PdfPages('ancestry_pca.pdf') as pdf:
            for x_pc, y_pc in [(1, 2), (3, 4), (5, 6)]:
                x, y = f'PC{x_pc}_SUM', f'PC{y_pc}_SUM'
                fig, ax = plt.subplots(figsize=(10, 7))
                for group in groups:
                    subset = reference.loc[reference['SUPER_POP'] == group]
                    ax.scatter(subset[x], subset[y], c=colors[group], s=12, alpha=0.30,
                               marker='o', label='Reference ' + group)
                for group in matrix_labels:
                    subset = study.loc[assignments['assigned_group'] == group]
                    if not subset.empty:
                        ax.scatter(subset[x], subset[y], c=colors[group], s=32, alpha=0.85,
                                   marker='x' if group == 'unassigned' else '^',
                                   label='Study ' + group)
                ax.set(xlabel=x, ylabel=y, title='1KG reference and study ancestry assignments')
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)
        print('Ancestry assignment complete; study samples:', len(study))
        print('Study assignment counts:', report['study_assignment_counts'])
        PY
    >>>

    output {
        File ancestry_assignments = "ancestry_assignments.tsv"
        File validation_predictions = "validation_predictions.tsv"
        File validation_report = "validation_report.json"
        File ancestry_pca = "ancestry_pca.pdf"
        File ancestry_model = "ancestry_model.joblib"
    }

    runtime {
        docker: "python:3.11.11-slim-bookworm"
        cpu: 2
        memory: "4 GB"
        disks: "local-disk 10 SSD"
    }
}
