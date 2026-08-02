#!/usr/bin/env python3
"""Module 3: Drug-Gene Interaction Mining via DGIdb GraphQL API.

Queries DGIdb v5 (GraphQL) for all significant CRISPR screen hits and compiles
a drug-gene interaction table with approval status and references.
"""

import csv
import os
import sys
import time
from pathlib import Path

import httpx

# --- Configuration ---
ALL_GENES = [
    "ASL", "GBA3", "ATP6V0A4", "TPCN1",          # depleted
    "PLA2G4E", "SLC1A1", "RRM2", "ARSD", "GRIN1", # enriched
    "NAT2", "SLC25A20", "ALOX15"                   # other
]

GENE_CLASS = {
    "ASL": "Depleted", "GBA3": "Depleted", "ATP6V0A4": "Depleted",
    "TPCN1": "Depleted", "PLA2G4E": "Enriched", "SLC1A1": "Enriched",
    "RRM2": "Enriched", "ARSD": "Enriched", "GRIN1": "Enriched",
    "NAT2": "Other", "SLC25A20": "Other", "ALOX15": "Other"
}

RESULTS_DIR = Path("../results/03_drug_interactions")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DGIDB_GRAPHQL_URL = "https://dgidb.org/api/graphql"

# GraphQL query for gene interactions (reverse-engineered from dgidb.org JS bundle)
GENE_INTERACTIONS_QUERY = """
query drugGeneInteractions($names: [String!]!) {
    genes(names: $names) {
        nodes {
            name
            conceptId
            interactions {
                id
                interactionScore
                interactionTypes {
                    type
                    directionality
                }
                drug {
                    name
                    conceptId
                    approved
                    drugApprovalRatings {
                        rating
                    }
                }
                publications {
                    pmid
                }
                sources {
                    fullName
                }
            }
        }
    }
}
"""


def query_dgidb_graphql(gene: str, client: httpx.Client) -> dict:
    """Query DGIdb GraphQL API for a single gene's interactions."""
    variables = {"names": [gene]}
    payload = {"query": GENE_INTERACTIONS_QUERY, "variables": variables}
    resp = client.post(DGIDB_GRAPHQL_URL, json=payload, timeout=30)
    resp.raise_for_status()
    return resp.json()


def parse_interactions(data: dict, gene: str) -> list[dict]:
    """Parse DGIdb GraphQL response into flat interaction records."""
    results = []
    try:
        nodes = data["data"]["genes"]["nodes"]
    except (KeyError, TypeError):
        return results

    for node in nodes:
        if node is None or node.get("name") != gene:
            continue
        interactions = node.get("interactions") or []

        for ix in interactions:
            drug_info = ix.get("drug") or {}
            drug_name = drug_info.get("name", "N/A")
            approved = drug_info.get("approved", False)
            approval_ratings = ", ".join(
                r["rating"] for r in (drug_info.get("drugApprovalRatings") or [])
            )

            interaction_types = ", ".join(
                it["type"] for it in (ix.get("interactionTypes") or [])
            )

            sources = ", ".join(
                s["fullName"] for s in (ix.get("sources") or [])
            )

            pmids = ", ".join(
                str(p["pmid"]) for p in (ix.get("publications") or []) if p.get("pmid")
            )

            score = ix.get("interactionScore", "N/A")

            # Determine approval status label
            if approved:
                approval_status = "Approved"
            elif approval_ratings:
                approval_status = approval_ratings
            else:
                approval_status = "Not Approved"

            results.append({
                "gene": gene,
                "gene_class": GENE_CLASS.get(gene, "N/A"),
                "drug_name": drug_name,
                "interaction_types": interaction_types,
                "sources": sources,
                "score": score,
                "approval_status": approval_status,
                "pmids": pmids,
            })

    return results


def main():
    print(f"[1/3] Querying DGIdb v5 GraphQL API for {len(ALL_GENES)} genes...")
    all_interactions = []

    with httpx.Client(verify=True) as client:
        for i, gene in enumerate(ALL_GENES):
            print(f"  [{i+1}/{len(ALL_GENES)}] {gene}...", end=" ")
            try:
                data = query_dgidb_graphql(gene, client)
                gene_interactions = parse_interactions(data, gene)
                all_interactions.extend(gene_interactions)
                print(f"{len(gene_interactions)} interactions")
            except Exception as e:
                print(f"ERROR: {e}")
            time.sleep(0.3)  # rate limit

    # --- Save full table ---
    print(f"\n[2/3] Saving {len(all_interactions)} interactions...")
    output_file = RESULTS_DIR / "drug_gene_interactions.csv"
    if all_interactions:
        fieldnames = [
            "gene", "gene_class", "drug_name", "interaction_types",
            "sources", "score", "approval_status", "pmids"
        ]
        with open(output_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_interactions)
        print(f"  -> {output_file}")
    else:
        print("  WARNING: No interactions found for any gene")
        with open(output_file, "w") as f:
            f.write("gene,gene_class,drug_name,interaction_types,sources,score,approval_status,pmids\n")
            f.write("No interactions found,,,,,,\n")

    # --- Summary table (gene-level) ---
    print("[3/3] Generating gene-level summary...")
    gene_summary = []
    for gene in ALL_GENES:
        gene_ix = [ix for ix in all_interactions if ix["gene"] == gene]
        approved = [ix for ix in gene_ix if ix["approval_status"] == "Approved"]
        gene_summary.append({
            "gene": gene,
            "class": GENE_CLASS.get(gene, "N/A"),
            "total_interactions": len(gene_ix),
            "approved_drugs": len(approved),
            "unique_drugs": len(set(ix["drug_name"] for ix in gene_ix)) if gene_ix else 0,
            "sources": ", ".join(sorted(set(
                s for ix in gene_ix for s in ix["sources"].split(", ") if s
            ))) if gene_ix else "",
        })

    summary_file = RESULTS_DIR / "drug_summary.csv"
    if gene_summary:
        with open(summary_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=gene_summary[0].keys())
            writer.writeheader()
            writer.writerows(gene_summary)
        print(f"  -> {summary_file}")

    # --- Print highlights ---
    genes_with_drugs = [gs for gs in gene_summary if gs["total_interactions"] > 0]
    print(f"\n[done] {len(genes_with_drugs)}/{len(ALL_GENES)} genes have known drug interactions")
    total_drugs = sum(gs["unique_drugs"] for gs in gene_summary)
    total_approved = sum(gs["approved_drugs"] for gs in gene_summary)
    print(f"       {len(all_interactions)} total interactions, {total_drugs} unique drugs, {total_approved} approved")


if __name__ == "__main__":
    main()
