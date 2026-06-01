from pathlib import Path
import pandas as pd
from data_preprocess import load_molecules, load_adverse_effects
from config import BASE_DIR, DATABASE_PATH, ENGINE, MOLECULES_FILE, ADVERSE_EFFECTS_FILE
from sqlalchemy import text

def build_database(mol_tsv: Path, adverse_tsv: Path, database_path: Path) -> None:
    # Load the data from TSV files
    molecules = load_molecules(mol_tsv, write_processed=False)
    adverse = load_adverse_effects(adverse_tsv, write_processed=False)

    # Prepare the drugs table
    drugs = molecules[
        [
            "chembl_id",
            "name",
            "drugType",
            "blackBoxWarning",
            "yearOfFirstApproval",
            "maximumClinicalTrialPhase",
            "hasBeenWithdrawn",
            "isApproved",
            "description",
            "target_count",
        ]
    ].copy()

    # Follow the snake_case convention for column names in the database, for easier to query and maintain consistency across the database schema.
    drugs = drugs.rename(
        columns={
            "drugType": "drug_type",
            "blackBoxWarning": "black_box_warning",
            "yearOfFirstApproval": "year_of_first_approval",
            "maximumClinicalTrialPhase": "maximum_clinical_trial_phase",
            "hasBeenWithdrawn": "has_been_withdrawn",
            "isApproved": "is_approved",
        }
    )
    
    print(f"drugs.head(3): {drugs.head(3).to_string()}")
    
    # Collect one row for each drug-target relationship
    target_rows = []
    for _, row in molecules.iterrows():
        target_ids = row["linked_target_ids"] # list of Ensembl IDs for the targets associated with the drug.
        target_symbols = row["linked_target_symbols"] # list of gene symbols corresponding to the Ensembl IDs, if available.
        
        for index, ensembl_id in enumerate(target_ids):
            if index < len(target_symbols):
                gene_symbol = target_symbols[index]
            else:
                gene_symbol = None
                
            target_rows.append(
                {
                    "chembl_id": row["chembl_id"],
                    "ensembl_id": ensembl_id,
                    "gene_symbol": gene_symbol,
                }
            )

    drug_targets = pd.DataFrame(target_rows)
    print(f"drug_targets.head(3): {drug_targets.head(3).to_string()}")
    
    adverse_summary = (
        adverse.groupby("chembl_id")
        .agg(
            adverse_event_count=("event", "count"),
            median_llr=("llr", "median"),
            total_reports=("count", "sum"),
        )
        .reset_index()
    )

    with ENGINE.begin() as connection:
        drugs.to_sql(
            "drugs",
            connection,
            if_exists="replace",
            index=False,
        )

        drug_targets.to_sql(
            "drug_targets",
            connection,
            if_exists="replace",
            index=False,
        )

        adverse.to_sql(
            "adverse_effects",
            connection,
            if_exists="replace",
            index=False,
        )

        adverse_summary.to_sql(
            "adverse_effect_summary",
            connection,
            if_exists="replace",
            index=False,
        )

        connection.execute(
            text("CREATE INDEX IF NOT EXISTS idx_drugs_chembl ON drugs(chembl_id)")
        )

        connection.execute(
            text("CREATE INDEX IF NOT EXISTS idx_targets_chembl ON drug_targets(chembl_id)")
        )

        connection.execute(
            text("CREATE INDEX IF NOT EXISTS idx_targets_gene ON drug_targets(gene_symbol)")
        )

        connection.execute(
            text("CREATE INDEX IF NOT EXISTS idx_adverse_chembl ON adverse_effects(chembl_id)")
        ) 
        
        connection.execute(
            text(
                "CREATE INDEX IF NOT EXISTS idx_adverse_summary_chembl "
                "ON adverse_effect_summary(chembl_id)"
            )
        )

    print(f"Database built successfully: path: {database_path}")
     
if __name__ == "__main__":

    mol_tsv = Path(MOLECULES_FILE)
    adverse_tsv = Path(ADVERSE_EFFECTS_FILE)
    database_path = Path(DATABASE_PATH)
    
    build_database(mol_tsv, adverse_tsv, database_path)
