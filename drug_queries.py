import math
import pandas as pd
from sqlalchemy import text
from config import ENGINE

def dataframe_to_json_safe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Replace NaN, Infinity, and -Infinity values with None-safe values.
    """
    df = df.copy()
    df = df.astype(object).where(pd.notnull(df), None)

    for column in df.columns:
        df[column] = df[column].apply(
            lambda value: None
            if isinstance(value, float) and not math.isfinite(value)
            else value
        )

    return df

def dataframe_to_json_records(df: pd.DataFrame) -> list[dict]:
    """
    Convert a pandas DataFrame into JSON-safe records.
    """
    return dataframe_to_json_safe(df).to_dict(orient="records")
    
def get_ranked_drugs_for_target(target: str, limit: int = 20) -> pd.DataFrame:
    """
    Retrieve and rank drugs associated with a HUGO gene symbol or Ensembl ID.

    Ranking score:
        normalized_target_count
      + normalized_median_llr
      + missing_adverse_data_penalty

    Lower score is better.
    """
    # Column issing_adverse_data_penalty: If median_llr is missing, give penalty = 1. If median_llr exists, give penalty = 0.
    query = text(
        """
        WITH candidate_drugs AS (
            SELECT
                d.chembl_id,
                d.name,
                d.drug_type,
                d.is_approved,
                d.maximum_clinical_trial_phase,
                d.target_count,
                t.ensembl_id,
                t.gene_symbol,
                s.adverse_event_count,
                s.median_llr,
                s.total_reports,

                CASE
                    WHEN s.median_llr IS NULL THEN 1
                    ELSE 0
                END AS missing_adverse_data_penalty

            FROM drug_targets AS t
            JOIN drugs AS d
                ON t.chembl_id = d.chembl_id
            LEFT JOIN adverse_effect_summary AS s
                ON d.chembl_id = s.chembl_id
            WHERE UPPER(t.gene_symbol) = UPPER(:target)
               OR UPPER(t.ensembl_id) = UPPER(:target)
        ),
        filled AS (
            SELECT
                *,
                COALESCE(
                    median_llr,
                    MAX(median_llr) OVER (),
                    0
                ) AS median_llr_for_ranking
            FROM candidate_drugs
        ),
        normalized AS (
            SELECT
                *,
                CASE
                    WHEN MAX(target_count) OVER () = MIN(target_count) OVER ()
                    THEN 0
                    ELSE
                        1.0 * (target_count - MIN(target_count) OVER ())
                        / (MAX(target_count) OVER () - MIN(target_count) OVER ())
                END AS normalized_target_count,

                CASE
                    WHEN MAX(median_llr_for_ranking) OVER () = MIN(median_llr_for_ranking) OVER ()
                    THEN 0
                    ELSE
                        1.0 * (median_llr_for_ranking - MIN(median_llr_for_ranking) OVER ())
                        / (MAX(median_llr_for_ranking) OVER () - MIN(median_llr_for_ranking) OVER ())
                END AS normalized_median_llr

            FROM filled
        )
        SELECT
            chembl_id,
            name,
            drug_type,
            is_approved,
            maximum_clinical_trial_phase,
            target_count,
            ensembl_id,
            gene_symbol,
            adverse_event_count,
            median_llr,
            total_reports,
            missing_adverse_data_penalty,
            normalized_target_count,
            normalized_median_llr,
            normalized_target_count
            + normalized_median_llr
            + missing_adverse_data_penalty AS ranking_score
        FROM normalized
        ORDER BY
            ranking_score ASC,
            missing_adverse_data_penalty ASC,
            normalized_target_count ASC,
            normalized_median_llr ASC
        LIMIT :limit
        """
    )

    with ENGINE.begin() as connection:
        result = pd.read_sql(
            query,
            connection,
            params={
                "target": target,
                "limit": limit,
            },
        )

    return dataframe_to_json_safe(result)
    
def list_available_targets(limit: int = 20) -> pd.DataFrame:
    """
    List targets available in the database.
    """
    query = text(
        """
        SELECT
            ensembl_id,
            gene_symbol,
            COUNT(DISTINCT chembl_id) AS drug_count
        FROM drug_targets
        GROUP BY ensembl_id, gene_symbol
        ORDER BY drug_count DESC
        LIMIT :limit
        """
    )

    with ENGINE.begin() as connection:
        result = pd.read_sql(
            query,
            connection,
            params={"limit": limit},
        )

    return dataframe_to_json_safe(result)
