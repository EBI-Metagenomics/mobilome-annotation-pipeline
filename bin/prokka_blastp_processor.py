#!/usr/bin/env python3
"""
Process BLASTP tabular output exactly like Prokka does.
This script replicates Prokka's cleanup_product() function and UniProt parsing
for tabular format with stitle field.
"""

import sys
import argparse
import re
from pathlib import Path
import pandas as pd

# Prokka constants
HYPO = "hypothetical protein"
GOOD_PROD = {
    "rep",
    "Conserved virulence factor B",
    "Cypemycin N-terminal methyltransferase",
}


def cleanup_product(product):
    """
    Exact replica of Prokka's cleanup_product() function.
    """
    if not product or product.strip() == "":
        return HYPO

    p = product.strip()

    # Check the whitelist first
    if p in GOOD_PROD:
        return p

    # Return hypothetical if matches certain patterns
    if re.search(
        r"DUF\d|UPF\d|conserved|domain of unknown|\b[CN]\.term|paralog",
        p,
        re.IGNORECASE,
    ):
        return HYPO

    # Return hypothetical if no lowercase letters
    if not re.search(r"[a-z]", p):
        return HYPO

    # Clean up the product name
    # eg. Leu/Ile/Val-binding protein homolog 3
    p = re.sub(r"\bhomolog( \d)?\b", "", p)
    p = re.sub(r"^arCOG\d+\s+", "", p)
    p = re.sub(r"\((EC|COG).*?\)", "", p)

    # Remove possible locus tags unless it's an IS element
    if not re.search(r"\bIS\d+\b", p):
        p = re.sub(r"\s+\w+\d{4,}c?", "", p)

    p = re.sub(r" and (inactivated|related) \w+", "", p)
    p = re.sub(r",\s*family$", "", p)

    # Replace various prefixes with "putative"
    p = re.sub(
        r"^(potential|possible|probable|predicted|uncharacteri.ed)",
        "putative",
        p,
        flags=re.IGNORECASE,
    )

    # Add "protein" suffix if needed
    if (
        re.search(
            r"(domain|family|binding|fold|like|motif|repeat)\s*$", p, re.IGNORECASE
        )
        and "," not in p
    ):
        p += " protein"

    # Collapse multiple spaces
    p = re.sub(r"\s+", " ", p)

    return p.strip()


def parse_prokka_uniprot_description(description):
    """
    Parse Prokka's UniProt format with ~~~ separators.
    Format: "EC_number~~~gene_name~~~product_description~~~COG_id"
    """
    gene = ""
    product = description
    ec_number = ""
    cog = ""

    # Check for Prokka's ~~~ format
    if "~~~" in description:
        parts = description.split("~~~")
        if len(parts) >= 3:
            # Format: EC~~~gene~~~product~~~COG (COG is optional)
            ec_number = parts[0].strip() if parts[0].strip() else ""
            gene = parts[1].strip() if parts[1].strip() else ""
            product = parts[2].strip() if parts[2].strip() else ""

            # COG is optional (4th field)
            if len(parts) >= 4:
                cog = parts[3].strip() if parts[3].strip() else ""

            # Collapse transitionary EC numbers (n1, n2, etc. -> -)
            if ec_number:
                ec_number = re.sub(r"n\d+", "-", ec_number)
                # Clean up EC numbers that are just dashes or incomplete
                if ec_number in ["-", "-.-.-.", ""] or ec_number.endswith(".-"):
                    ec_number = ""
        else:
            # Fallback if format is unexpected
            product = description
    else:
        # No ~~~ separators, treat entire description as product
        product = description

    return gene, product, ec_number, cog


def process_prokka_blastp_tabular(input_file, output_file, raw_product=False):
    """
    Process BLASTP tabular results exactly like Prokka does.

    Args:
        input_file: Path to BLAST tabular output file with stitle field
        output_file: Path to output annotation file
        raw_product: If True, don't clean product names
    """

    # Column names for the specific BLAST output format
    columns = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qlen",
        "slen",
        "stitle",
    ]

    try:
        # Read BLAST results
        df = pd.read_csv(input_file, sep="\t", names=columns, comment="#")
        print(f"Loaded {len(df)} BLAST hits (already filtered by BLAST)")

        # Since BLAST already applied -evalue 1E-9 and -qcov_hsp_perc 80,
        # we don't need to filter again

        # Group by query and take the best hit (should be only one due to -num_descriptions 1)
        # But just in case, take the one with lowest e-value, highest bitscore
        best_hits = df.loc[df.groupby("qseqid")["evalue"].idxmin()]

        annotations = []
        num_cleaned = 0

        for _, hit in best_hits.iterrows():
            qseqid = hit["qseqid"]
            sseqid = hit["sseqid"]
            stitle = hit["stitle"]  # This contains the full description with ~~~ format

            # Parse the description using Prokka's format
            gene, product, ec_number, cog = parse_prokka_uniprot_description(stitle)

            # Clean product name unless raw_product is True
            note = ""
            if not raw_product:
                original_product = product
                clean_product = cleanup_product(product)
                if clean_product != original_product:
                    print(f"Modify product: {original_product} => {clean_product}")

                    # If product becomes hypothetical, clear other annotations and add note
                    if clean_product == HYPO:
                        note = original_product
                        gene = ""
                        ec_number = ""
                        cog = ""

                    num_cleaned += 1
                product = clean_product

            # Calculate actual query coverage (like Prokka does)
            query_coverage = ((hit["qend"] - hit["qstart"] + 1) / hit["qlen"]) * 100

            # Create annotation record matching Prokka's output
            annotation = {
                "query_id": qseqid,
                "subject_id": sseqid,
                "product": product,
                "gene": gene,
                "ec_number": ec_number,
                "cog": cog,
                "note": note,
                "evalue": hit["evalue"],
                "identity": hit["pident"],
                "coverage": round(query_coverage, 1),
                "inference": f"similar to AA sequence:UniProtKB:{sseqid}",
            }
            annotations.append(annotation)

        # Create output DataFrame
        result_df = pd.DataFrame(annotations)

        # Remove empty columns for cleaner output
        for col in ["gene", "ec_number", "cog", "note"]:
            if result_df[col].eq("").all():
                result_df = result_df.drop(columns=[col])

        # Write results
        result_df.to_csv(output_file, sep="\t", index=False)
        print(f"Wrote {len(result_df)} annotations to {output_file}")
        print(f"Cleaned {num_cleaned} product names")

        # Print summary statistics
        hypothetical_count = len(result_df[result_df["product"] == HYPO])
        annotated_count = len(result_df) - hypothetical_count

        print(f"\nSummary:")
        print(f"  Total annotations: {len(result_df)}")
        print(f"  Functional annotations: {annotated_count}")
        print(f"  Hypothetical proteins: {hypothetical_count}")

        if "ec_number" in result_df.columns:
            print(f"  With EC numbers: {len(result_df[result_df['ec_number'] != ''])}")
        if "gene" in result_df.columns:
            print(f"  With gene names: {len(result_df[result_df['gene'] != ''])}")
        if "cog" in result_df.columns:
            print(f"  With COG assignments: {len(result_df[result_df['cog'] != ''])}")

        return result_df

    except Exception as e:
        print(f"Error processing BLASTP results: {e}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Process BLASTP tabular output exactly like Prokka does",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python prokka_blastp_tabular.py -i blastp_results.tsv -o annotations.tsv
  python prokka_blastp_tabular.py -i results.tsv -o anno.tsv --raw-product
  
Note: Input must be BLAST tabular format with stitle field:
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle"
        """,
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input BLAST tabular file with stitle field",
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output annotation file (TSV format)"
    )
    parser.add_argument(
        "--raw-product",
        action="store_true",
        help="Do not clean product names (like Prokka --rawproduct)",
    )

    args = parser.parse_args()

    # Check input file exists
    if not Path(args.input).exists():
        print(f"Error: Input file {args.input} does not exist")
        sys.exit(1)

    # Process the results
    process_prokka_blastp_tabular(
        input_file=args.input, output_file=args.output, raw_product=args.raw_product
    )


if __name__ == "__main__":
    main()
