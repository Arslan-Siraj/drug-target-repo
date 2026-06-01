import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import requests
from config import API_URL

def fetch_drugs_for_target(target: str, limit: int = 20) -> pd.DataFrame:
    """
    Query the API for drugs associated with a target gene symbol or Ensembl ID.
    """
    url = f"{API_URL}/targets/{target}/drugs"

    response = requests.get(
        url,
        params={"limit": limit},
        timeout=30,
    )

    if response.status_code == 404:
        raise ValueError(f"No drugs found for target '{target}'.")

    response.raise_for_status()

    data = response.json()
    return pd.DataFrame(data["results"])


def plot_and_csv_ranked_drugs(df: pd.DataFrame, target: str, output_dir: Path) -> Path:
    """
    Save a bar plot of the top-ranked drugs.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    output_base = output_dir / f"{target.upper()}_drug_ranking"
    csv_path = output_base.with_suffix(".csv")
    png_path = output_base.with_suffix(".png")
    
    df.to_csv(csv_path, index=False)
    
    plot_df = df.sort_values("ranking_score").head(20).copy()
    plot_df["label"] = plot_df["name"].fillna(plot_df["chembl_id"])
    
    plt.figure(figsize=(10, 6))
    plt.barh(plot_df["label"], plot_df["ranking_score"])
    plt.xlabel("Ranking score")
    plt.ylabel("Drug")
    plt.title(f"Top candidate drugs for {target.upper()}")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(png_path, dpi=300)
    plt.close()

    return output_base


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Retrieve and visualize ranked drugs for a target gene symbol "
            "or Ensembl ID."
        )
    )

    parser.add_argument(
        "--target",
        required=True,
        help="Target gene symbol or Ensembl ID, for example BRAF or ENSG00000149295.",
    )

    parser.add_argument(
        "--limit",
        type=int,
        default=20,
        help="Maximum number of drugs to retrieve.",
    )

    args = parser.parse_args()

    df = fetch_drugs_for_target(args.target, args.limit)

    columns_to_show = [
        "chembl_id",
        "name",
        "drug_type",
        "is_approved",
        "maximum_clinical_trial_phase",
        "target_count",
        "ensembl_id",
        "gene_symbol",
        "median_llr",
        "adverse_event_count",
        "ranking_score",
    ]

    available_columns = [col for col in columns_to_show if col in df.columns]

    print(f"\nTop drugs for target: {args.target.upper()}\n")
    print(df[available_columns].to_string(index=False))

    output_path = plot_and_csv_ranked_drugs(
        df=df,
        target=args.target,
        output_dir=Path("outputs"),
    )

    print(f"\nSaved results to: {output_path}")
    

if __name__ == "__main__":
    main()