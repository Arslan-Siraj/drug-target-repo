import pandas as pd
from fastapi import FastAPI, HTTPException, Query
from sqlalchemy import text
from config import DATABASE_PATH, ENGINE
from drug_queries import dataframe_to_json_records, get_ranked_drugs_for_target, list_available_targets

app = FastAPI(
    title="Drug-Target Example API",
    description=(
        "Retrieve compounds for a target gene symbol or Ensembl ID, with ranking score."
    ),
)

@app.get("/")
def root():
    return {
        "message": "Drug-Target API is running.",
        "examples": [
            "/targets/BRAF/drugs",
            "/targets/ENSG00000149295/drugs",
            "/targets/ENSG00000149295/drugs?limit=5",
            "/targets?limit=10",
            "/health",
        ],
    }
    
@app.get("/health")
def health_check():
    """
    Check whether the API can connect to the SQLite database.
    """
    if not DATABASE_PATH.exists():
        raise HTTPException(
            status_code=500,
            detail=f"Database file not found: {DATABASE_PATH}",
        )

    try:
        with ENGINE.begin() as connection:
            # Get all tables from database
            tables = pd.read_sql(
                text(
                    """
                    SELECT name
                    FROM sqlite_master
                    WHERE type = 'table'
                    ORDER BY name
                    """
                ),
                connection,
            )
            
    except Exception as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Database connection failed: {exc}",
        )

    return {
    "status": "ok",
    "database": str(DATABASE_PATH),
    "tables": tables["name"].tolist(),
    }

# Example: /targets?limit=10
@app.get("/targets")
def available_targets(
    limit: int = Query(default=20, ge=1, le=100), # ge: lower limit, le: upper limit (max allowed)
):
    """
    List available targets in the database.
    """
    try:
        result = list_available_targets(limit=limit)
        
    except Exception as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Database query failed: {exc}",
        )
        
    # Transform to json format, FASTAPI retreive in json
    records = dataframe_to_json_records(result)

    return {
        "n_results": len(records),
        "results": records, 
    }

# Example: /targets/EGFR/drugs?limit=10
# /targets/ENSG00000149295/drugs?limit=5
@app.get("/targets/{target}/drugs")
def get_drugs_for_target(
    target: str,
    limit: int = Query(default=20, ge=1, le=100),  
):
    """
    Retrieve ranked drugs for a HUGO gene symbol or Ensembl gene ID.
    """
    try:
        result = get_ranked_drugs_for_target(
            target=target,
            limit=limit,
        )
        
    except Exception as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Database query failed: {exc}",
        )

    if result.empty:
        raise HTTPException(
            status_code=404,
            detail=f"No drugs found for this target '{target}'.",
        )

    records = dataframe_to_json_records(result)

    return {
        "target": target.upper(),  # search in upper case
        "n_results": len(records),
        "ranking_method": (
            "ranking_score = normalized_target_count + normalized_median_llr + missing_adverse_data_penalty. Lower scores are better. "
        ),
        "results": records,
    }
    
