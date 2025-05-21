#!/usr/bin/env python3

from Bio import Entrez, SeqIO
import time
import os
import csv
import matplotlib.pyplot as plt

class NCBIRetriever:
    def __init__(self, email, api_key):
        self.email = email
        self.api_key = api_key
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = "BioScriptEx10"

    def search_taxid(self, taxid):
        print(f"Searching for records with taxID: {taxid}")
        try:
            # Get taxonomy information
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            organism_name = records[0]["ScientificName"]
            print(f"Organism: {organism_name} (TaxID: {taxid})")

            # Search for records
            search_term = f"txid{taxid}[Organism]"
            handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
            search_results = Entrez.read(handle)
            count = int(search_results["Count"])

            if count == 0:
                print(f"No records found for {organism_name}")
                return None

            print(f"Found {count} records")
            self.webenv = search_results["WebEnv"]
            self.query_key = search_results["QueryKey"]
            self.count = count

            return count

        except Exception as e:
            print(f"Error searching TaxID {taxid}: {e}")
            return None

    def fetch_records(self, start=0, max_records=100):
        if not hasattr(self, 'webenv') or not hasattr(self, 'query_key'):
            print("No search results to fetch. Run search_taxid() first.")
            return []

        try:
            batch_size = min(max_records, 500)
            handle = Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=self.webenv,
                query_key=self.query_key
            )
            records_text = handle.read()
            return records_text
        except Exception as e:
            print(f"Error fetching records: {e}")
            return ""

def filter_and_parse_sequences(gb_text, min_len, max_len):
    from io import StringIO
    handle = StringIO(gb_text)
    filtered_data = []

    for record in SeqIO.parse(handle, "genbank"):
        seq_len = len(record.seq)
        if min_len <= seq_len <= max_len:
            filtered_data.append({
                "accession": record.id,
                "length": seq_len,
                "description": record.description
            })

    return filtered_data

def save_csv(data, filename):
    """Save filtered sequence data to CSV."""
    with open(filename, mode="w", newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=["accession", "length", "description"])
        writer.writeheader()
        for row in data:
            writer.writerow(row)
    print(f"CSV saved to {filename}")

def plot_data(data, filename):
    """Plot a line chart of sequence lengths and save as PNG."""
    data_sorted = sorted(data, key=lambda x: x["length"], reverse=True)
    accessions = [d["accession"] for d in data_sorted]
    lengths = [d["length"] for d in data_sorted]

    plt.figure(figsize=(12, 6))
    plt.plot(accessions, lengths, marker='o')
    plt.xticks(rotation=90)
    plt.xlabel("GenBank Accession Number")
    plt.ylabel("Sequence Length")
    plt.title("Sequence Lengths by GenBank Accession")
    plt.tight_layout()
    plt.savefig(filename)
    print(f"Chart saved to {filename}")
    plt.close()

def main():
    email = input("Enter your email address for NCBI: ")
    api_key = input("Enter your NCBI API key: ")
    taxid = input("Enter taxonomic ID (taxid) of the organism: ")
    min_len = int(input("Enter minimum sequence length: "))
    max_len = int(input("Enter maximum sequence length: "))

    retriever = NCBIRetriever(email, api_key)

    count = retriever.search_taxid(taxid)
    if not count:
        print("No records found. Exiting.")
        return

    print("\nFetching records from GenBank...")
    all_records_text = retriever.fetch_records(start=0, max_records=100)

    print("Filtering and parsing sequences...")
    filtered_data = filter_and_parse_sequences(all_records_text, min_len, max_len)

    if not filtered_data:
        print("No sequences matched the length filter.")
        return

    csv_file = f"taxid_{taxid}_filtered.csv"
    chart_file = f"taxid_{taxid}_chart.png"

    save_csv(filtered_data, csv_file)
    plot_data(filtered_data, chart_file)

    print("\nTask complete!")

if __name__ == "__main__":
    main()
