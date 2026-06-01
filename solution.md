# Drug-Target Investigation

This work provides a simple pipeline to query drug-target data and rank candidate drugs for a target of interest.
The solution include retrieve via command-line API client or in Streamlit web application

## Environment

Create and activate the environment:

```bash
conda create -n drugs-target python=3.11
conda activate drugs-target
```

Install dependencies:

```bash
pip install -r requirements.txt
```

## Build the database

Run:

```bash
python build_database.py
```

This script loads the data from the `data/` directory, preprocesses the drug and adverse-effect files, and creates:

```text
drug_database.db
```

The database contains:

```text
drugs
drug_targets
adverse_effects
adverse_effect_summary
```

## Run the API

Start the FastAPI server:

```bash
uvicorn api:app --reload
```

The API runs at:

```text
http://127.0.0.1:8000
```

Check that the API is working:

```bash
curl "http://127.0.0.1:8000/health"
```

List available targets:

```bash
curl "http://127.0.0.1:8000/targets?limit=10"
```

Query ranked drugs for a target:

```bash
curl "http://127.0.0.1:8000/targets/ENSG00000149295/drugs?limit=20"
```

If Ensembl-to-HUGO mapping worked, HUGO gene symbols can also be used:

```bash
curl "http://127.0.0.1:8000/targets/BRAF/drugs?limit=20"
```

## Run the command-line client

Keep the API running in one terminal.

In another terminal, run:

```bash
python run_API.py --target ENSG00000149295 --limit 20
```

This retrieves ranked drug candidates from the API.
Also, could see the results in output folder.

## Run the Streamlit web app

Run:

```bash
streamlit run webapp.py
```

The web app lets the user enter a target gene symbol or Ensembl ID, fetch ranked drugs, view a ranking plot, and download the results as a CSV file.

## Ranking score

The same ranking method is used in both the API and the Streamlit app.

```text
ranking_score = normalized_target_count
              + normalized_median_llr
              + missing_adverse_data_penalty
```

Lower scores are better.

#### columns

`target_count` is the number of targets associated with a drug.

`median_llr` is the median log-likelihood ratio across the drug's adverse-effect records.

`missing_adverse_data_penalty` is added when a drug has no adverse-effect summary.

---

#### Normalization

`target_count` and `median_llr` are normalized between 0 and 1 among the candidate drugs for the queried target:

```text
normalized_value = (value - minimum_value) / (maximum_value - minimum_value)
```

If all candidates have the same value, the normalized value is set to `0`.

This makes target specificity and adverse-effect burden comparable.

---

## Ensembl and HUGO identifiers

The preprocessing step attempts to map Ensembl gene IDs to HUGO gene symbols using BioThings/MyGene.

If mapping succeeds, users can query with gene symbols such as:

```text
BRAF
EGFR
JAK2
```

If mapping fails, users can still query by Ensembl ID, for example:

```text
ENSG00000149295
```

---

## Example workflow

```bash
python build_database.py
uvicorn api:app --reload 
python run_API.py --target ENSG00000149295 --limit 20 # From API
streamlit run webapp.py # Via WebApp
```
## Solution
Please checkout: [solution.md](https://github.com/Arslan-Siraj/drug-target-repo/blob/main/solution.md)