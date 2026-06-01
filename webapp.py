from pathlib import Path

import streamlit as st

from config import DATABASE_NAME, DATABASE_PATH
from drug_queries import get_ranked_drugs_for_target, list_available_targets 


def check_database() -> bool:
    """
    Check whether the SQLite database exists.
    """
    return DATABASE_PATH.exists()


def main() -> None:
    st.set_page_config(
        page_title="Drug-Target Example Task",
        layout="wide",
    )

    st.title("Drug-Target Example Task")
    st.write(
        "Find clinical or marketed compounds for a target and rank them by "
        "target specificity and adverse-effect association."
    )

    with st.sidebar:
        st.header("Database")
        st.write(f"Database name: `{DATABASE_NAME}`")
        st.write(f"Database path: `{DATABASE_PATH}`")

        if check_database():
            st.success("Database connected")
        else:
            st.error(f"Database file not found: {DATABASE_PATH}")
            st.stop()

        st.header("Search")

        target = st.text_input(
            "Target gene symbol or Ensembl ID",
            value="ENSG00000149295",
            help=(
                "Use a HUGO gene symbol such as BRAF if mapping worked, "
                "or an Ensembl ID such as ENSG00000149295."
            ),
        )

        limit = st.slider(
            "Number of drugs to return",
            min_value=5,
            max_value=100,
            value=20,
            step=5,
        )

        search_button = st.button("Fetch ranked drugs")

    st.subheader(
        "Available targets",
        help="If gene_symbol is empty, query using the Ensembl ID.",
    )

    st.write("These are the most common targets in the database.")

    try:
        available_targets = list_available_targets(limit=10)
        st.dataframe(available_targets, width="stretch")
    except Exception as exc:
        st.error(f"Could not load available targets: {exc}")
        st.stop()

    st.divider()

    if search_button:
        if not target.strip():
            st.warning("Please enter a target gene symbol or Ensembl ID.")
            st.stop()

        try:
            result = get_ranked_drugs_for_target(
                target=target.strip(),
                limit=limit,
            )
        except Exception as exc:
            st.error(f"Could not retrieve ranked drugs: {exc}")
            st.stop()

        if result.empty:
            st.warning(
                f"No drugs found for target '{target}'. "
                "Try an Ensembl ID from the available targets table."
            )
            st.stop()

        st.subheader(
            f"Ranked drug candidates for {target.upper()}",
            help=(
                "Ranking score = normalized_target_count + normalized_median_llr "
                "+ missing_adverse_data_penalty. Lower scores are better. "
                "The score prioritizes drugs with fewer targets and weaker "
                "typical adverse-effect associations. Drugs without adverse-effect "
                "summary data receive a penalty so that missing safety data is not "
                "treated as evidence of low adverse-effect risk."
            ),
        )

        display_columns = [
            "chembl_id",
            "name",
            "drug_type",
            "is_approved",
            "maximum_clinical_trial_phase",
            "target_count",
            "ensembl_id",
            "gene_symbol",
            "median_llr",
            "total_reports",
            "missing_adverse_data_penalty",
            "normalized_target_count",
            "normalized_median_llr",
            "ranking_score",
        ]

        existing_columns = [
            column for column in display_columns if column in result.columns
        ]

        st.dataframe(result[existing_columns], width="stretch")

        st.subheader("Ranking score plot")

        plot_df = result.copy()
        plot_df["drug_label"] = plot_df["name"].fillna(plot_df["chembl_id"])
        plot_df = plot_df.sort_values("ranking_score").head(50)

        st.bar_chart(
            plot_df.set_index("drug_label")["ranking_score"]
        )

        csv = result.to_csv(index=False).encode("utf-8")

        st.download_button(
            label="Download results as CSV",
            data=csv,
            file_name=f"{target.upper()}_ranked_drugs.csv",
            mime="text/csv",
        )


if __name__ == "__main__":
    main()