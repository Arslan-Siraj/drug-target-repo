import re
import requests
from pathlib import Path
from biothings_client import get_client
import mygene
import pandas as pd

def load_tsv_file(file_path: Path) -> pd.DataFrame:
    """
    Read a data file in the data directory (.tsv format).
    Example: molecules.tsv
    """
    if not file_path.exists():
        raise FileNotFoundError(
            f"Could not find '{file_path}'. Please check that the file exists."
        )

    if file_path.suffix != ".tsv":
        raise ValueError(
            f"Expected a '.tsv' file, but received '{file_path.suffix or 'no extension'}'."
        )

    return pd.read_csv(file_path, sep="\t")
    
def refine_chembl_id(value: str | int | float) -> str:
    """
    Refine ChEMBL identifiers.
    adverseEffects.tsv define the CHEMBL identifier as 1200632
    molecules.tsv define the CHEMBL identifier as CHEMBL1200632
    This function makes both formats consistent.
    Args:
      value: The input value for ChEMBL identifier.
    Returns:
      A refined ChEMBL identifier.
    """
    chembl_value = str(value).strip()

    if chembl_value.startswith("CHEMBL"):
        return chembl_value

    return f"CHEMBL{chembl_value}"
    
def refine_ensembl_targets(value) -> list[str]:
    """
    Refine linkedTargets values.
    Handles formats such as:
        [ENSG00000133019,ENSG00000181072]
        []
        NaN
        ENSG00000133019
    Args:
      value: The input value for linkedTargets.
    Returns:
      Returns a list of Ensembl gene IDs.
    """
    if value is None or pd.isna(value):
        return []

    if isinstance(value, list):
        return [str(v).strip() for v in value if str(v).strip()]

    text_value = str(value).strip()

    if text_value in {"", "[]", "nan", "None"}:
        return []

    return re.findall(r"ENSG\d+", text_value)
    
def map_ensembl_to_symbol(ensembl_ids: list[str], batch_size: int = 100) -> dict[str, str]:
    """
    Map Ensembl gene IDs to HUGO gene symbols.
    The function first tries the BioThings/MyGene Python client with caching
    disabled. If that fails, it falls back to the MyGene HTTP API directly.
    """
    unique_ids = sorted(set(ensembl_ids))

    # check if there are any unique ids to process
    if not unique_ids:
        return {}
    # start with an empty mapping dictionary
    mapping: dict[str, str] = {}

    # First attempt: BioThings Python client, with cache disabled.
    try:
        mg = get_client("gene")

        if hasattr(mg, "stop_caching"):
            mg.stop_caching()

        for start in range(0, len(unique_ids), batch_size):
            batch = unique_ids[start : start + batch_size]

            try:
                results = mg.querymany(
                    batch,
                    scopes="ensembl.gene",
                    fields="symbol", # we only need the symbol field for mapping
                    species="human", # we are only interested in human genes
                    as_dataframe=False,
                    verbose=False,
                )
            except Exception as exc:
                print(
                    f"Warning: BioThings client mapping failed for batch "
                    f"{start // batch_size + 1}: {exc}"
                )
                break

            for row in results:
                query = row.get("query")
                symbol = row.get("symbol")

                if query and symbol and not row.get("Not Found", False):
                    mapping[query] = symbol

        if mapping:
            print(
                f"Successfully mapped {len(mapping)} of {len(unique_ids)} "
                "Ensembl IDs to HUGO gene symbols using BioThings client."
            )
            return mapping

    except Exception as exc:
        print(f"Warning: BioThings client initialization failed: {exc}")

    # Second attempt: direct HTTP request to MyGene.info.
    print("Trying direct MyGene.info HTTP mapping fallback...")

    url = "https://mygene.info/v3/query"

    for start in range(0, len(unique_ids), batch_size):
        batch = unique_ids[start : start + batch_size]

        payload = {
            "q": ",".join(batch),
            "scopes": "ensembl.gene",
            "fields": "symbol",
            "species": "human",
        }

        try:
            response = requests.post(url, data=payload, timeout=120)
            response.raise_for_status()
            results = response.json()
        except Exception as exc:
            print(
                f"Warning: MyGene HTTP mapping failed for batch "
                f"{start // batch_size + 1}: {exc}"
            )
            continue

        for row in results:
            query = row.get("query")
            symbol = row.get("symbol")

            if query and symbol and not row.get("Not Found", False):
                mapping[query] = symbol

    if mapping:
        print(
            f"Successfully mapped {len(mapping)} of {len(unique_ids)} Ensembl IDs to HUGO gene symbols using HTTP fallback."
        )
    else:
        print(
            "Warning: no Ensembl IDs were mapped to HUGO gene symbols. Gene-symbol search may not work."
        )

    return mapping

def load_molecules(file_path: str, write_processed: bool = True) -> pd.DataFrame:
    """
    Load and process the molecules data (.tsv).
    Args:
      file_path: Path to the molecules.tsv file.
      write_processed: Whether to write the processed DataFrame to a .tsv file.
    Returns:
      A DataFrame with processed molecule data, including refined ChEMBL IDs and linked target information.
     """

    file_path = Path(file_path)
    print(f"Loading molecules from: {file_path}")

    molecules_df = load_tsv_file(file_path)

    # check for required columns
    required_columns = {"id", "linkedTargets"}
    missing_columns = required_columns - set(molecules_df.columns)

    if missing_columns:
        raise ValueError(
            f"molecules file is missing required columns: {missing_columns}"
        )
    
    # take the refine chembl_id add them to the molecules dataframe
    molecules_df["chembl_id"] = molecules_df["id"].map(refine_chembl_id)

    # take the target_ids from refine linkedTargets column to extract Ensembl IDs (e-g: some are empty lists)
    molecules_df["linked_target_ids"] = molecules_df["linkedTargets"].apply(
        refine_ensembl_targets
    )

    # count the number of linked targets for each molecule
    molecules_df["target_count"] = molecules_df["linked_target_ids"].apply(len)

    # extract all unique Ensembl target IDs across all molecules
    all_targets = [
        target
        for targets in molecules_df["linked_target_ids"]
        for target in targets
    ]

    print(f"Found {len(set(all_targets))} unique Ensembl target IDs.")
    
    # Map Ensembl IDs to HUGO gene symbols.
    target_symbol_map = map_ensembl_to_symbol(all_targets)
    
    if target_symbol_map:
        print(
            f"Successfully mapped {len(target_symbol_map)} "
            "Ensembl IDs to HUGO gene symbols."
        )
    else:
        print("Warning: no Ensembl IDs were mapped to HUGO gene symbols.")
    
    # Keep the same order as linked_target_ids.
    # Unmapped Ensembl IDs remain as None.
    molecules_df["linked_target_symbols"] = molecules_df["linked_target_ids"].apply(
        lambda ids: [target_symbol_map.get(ensembl_id) for ensembl_id in ids]
    )

    if write_processed:
       molecules_df.to_csv("processed_molecules.tsv", sep="\t", index=False)

    return molecules_df
    
def load_adverse_effects(file_path: str, write_processed: bool = True) -> pd.DataFrame:
    """
    Load and process the adverse effects data (.tsv).
    Args:
      file_path: Path to the adverseEffects.tsv file.
      write_processed: Whether to write the processed DataFrame to a CSV file.
    Returns:
      A DataFrame with processed adverse effect data, including refined ChEMBL IDs and numeric LLR and count values.
    """
    print(f"Loading adverse effects from: {file_path}")

    adverse = load_tsv_file(Path(file_path))
    
    # makwe sure the required columns are present
    required_columns = {"chembl_id", "event", "count", "llr"}
    missing_columns = required_columns - set(adverse.columns)

    if missing_columns:
        raise ValueError(
            f"adverseEffects file is missing required columns: {missing_columns}"
        )

    # refine chembl_id and convert llr and count to numeric values 
    adverse["chembl_id"] = adverse["chembl_id"].map(refine_chembl_id)
    adverse["llr"] = pd.to_numeric(adverse["llr"], errors="coerce") # coercing errors to NaN
    adverse["count"] = pd.to_numeric(adverse["count"], errors="coerce") # coercing errors to NaN

    if write_processed:
        adverse.to_csv("processed_adverse_effects.tsv", sep="\t", index=False)

    return adverse

# For testing individual test
#Molecules_file = load_molecules("data/molecules.tsv", write_processed=True)
#AdverseEffects_file = load_adverse_effects("data/adverseEffects.tsv", write_processed=True)
#print(f"Loaded {Molecules_file.shape[0]} molecules and {AdverseEffects_file.shape[0]} adverse effects.")
